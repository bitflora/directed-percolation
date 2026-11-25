//! Low-density series for directed percolation on the directed square lattice.
//!
//! Non-nodal path-union transfer matrix (after Jensen, J. Phys. A 32, 5233
//! (1999)), with these changes:
//!
//! * Translation reduction. States are cluster *shapes* (occupied columns
//!   relative to the leftmost one at the start of each row). Absolute position
//!   is carried as three moment polynomials M_k = sum_o o^k P(shape, o), which
//!   is all the x^2 moment needs. Translates of a shape share one entry.
//! * Mirror reduction. A shape and its mirror image are stored once (a
//!   symmetric shape with weight 1/2), canonicalized at the end of every row.
//! * Per-state coefficient windows [lo, hi]: lo bounds the bonds (sites)
//!   already placed, n - hi the bonds still needed to finish a non-nodal graph.
//!   Both come from an interval-cover argument (see `bound`).
//! * Sharded, multithreaded steps: all targets of one source share a group key,
//!   so groups are partitioned across threads without locking.
//!
//! Coordinates: row s has columns c; site (s, c) has bonds to (s+1, c) and
//! (s+1, c+1); x = 2c - s.
//!
//! Key layout while processing column j of row s (kink at j):
//!   bits 0..j-1 : row s+1 columns 0..j-1, final
//!   bit  j      : row s+1 column j, partial (received a right bond from j-1)
//!   bit  c+1    : row s column c, for c >= j, not yet processed
//! Processing column j rewrites bits (j, j+1) from (q, a) to (b_j, q_{j+1}).

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::io::Write;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

mod modp;
mod prof;
#[cfg(feature = "gpu")]
mod gpu;
mod spill;
mod table;

use modp::{Mod, C};
use table::{shard_of, Table};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Model {
    Bond,
    Site,
}

pub struct Params {
    pub n: usize,
    pub model: Model,
    pub tmax: usize,
    pub moments: usize, // 1 (S, mu01, mu02 only) or 3 (also mu20)
    pub verbose: bool,
    pub profile: bool, // DP_PROFILE=1: detailed counters and per-bound timing
    pub loose: bool, // testing only: disable order pruning (keep only the row limit)
    pub mirror: bool,
    pub threads: usize,
    pub shards: usize,
    pub slack: bool,
    pub nosign: bool, // diagnostics only: count graphs without signs (no cancellation)
    pub cheap_lo: bool, // experiment: lower end 2s + fixed for hashed targets, rely on trimming
    pub block: usize,   // parallel mode: kink steps per synchronization
    pub park: bool,     // single-thread mode: parked pass (--park): ~11% faster, ~46% more memory
    pub spill: Option<std::path::PathBuf>, // --spill DIR: out-of-core sweep (`run_spill`)
    pub spill_above: u64, // bytes: generations larger than this go to disk
    pub batch: bool, // two-pass steps with batched bounds (`step_batch`), the default; --no-batch for `step_one`
    pub gpu: bool,   // --gpu: batched bounds on the GPU (implies --batch)
    pub gpu_min: usize, // batches smaller than this stay on the CPU (DP_GPU_MIN)
    pub key_lo: bool, // DP_KEY_LO=1: batched steps take hashed targets' lower ends from the interval cover
}

/// Heap accounting (current and peak bytes), so reported peaks are measured
/// rather than estimated from table sizes.
pub static HEAP_CUR: AtomicUsize = AtomicUsize::new(0);
pub static HEAP_PEAK: AtomicUsize = AtomicUsize::new(0);
struct CountingAlloc;
#[inline(always)]
fn heap_add(n: usize) {
    let c = HEAP_CUR.fetch_add(n, Ordering::Relaxed) + n;
    HEAP_PEAK.fetch_max(c, Ordering::Relaxed);
}
unsafe impl GlobalAlloc for CountingAlloc {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let q = System.alloc(l);
        if !q.is_null() {
            heap_add(l.size());
        }
        q
    }
    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        let q = System.alloc_zeroed(l);
        if !q.is_null() {
            heap_add(l.size());
        }
        q
    }
    unsafe fn dealloc(&self, q: *mut u8, l: Layout) {
        System.dealloc(q, l);
        HEAP_CUR.fetch_sub(l.size(), Ordering::Relaxed);
    }
    unsafe fn realloc(&self, q: *mut u8, l: Layout, new: usize) -> *mut u8 {
        let r = System.realloc(q, l, new);
        if !r.is_null() {
            if new >= l.size() {
                heap_add(new - l.size());
            } else {
                HEAP_CUR.fetch_sub(l.size() - new, Ordering::Relaxed);
            }
        }
        r
    }
}
#[global_allocator]
static GLOBAL: CountingAlloc = CountingAlloc;
/// Return freed heap memory to the OS. The Windows process heap keeps freed
/// small and medium blocks committed, so without this the working set runs
/// well above the live heap after the many short-lived pieces of a parallel
/// step. Off with DP_NO_TRIM=1 (diagnostics).
#[cfg(windows)]
fn trim_heap() {
    use std::ffi::c_void;
    #[repr(C)]
    struct OptimizeInfo {
        version: u32,
        flags: u32,
    }
    extern "system" {
        fn GetProcessHeap() -> *mut c_void;
        fn HeapCompact(h: *mut c_void, flags: u32) -> usize;
        fn HeapSetInformation(h: *mut c_void, class: i32, info: *mut c_void, len: usize) -> i32;
    }
    static OFF: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    if *OFF.get_or_init(|| std::env::var("DP_NO_TRIM").is_ok()) {
        return;
    }
    unsafe {
        let h = GetProcessHeap();
        HeapCompact(h, 0);
        let mut info = OptimizeInfo { version: 1, flags: 0 };
        // HeapOptimizeResources = 3 (Windows 8.1 and later): decommit free memory
        HeapSetInformation(h, 3, &mut info as *mut OptimizeInfo as *mut c_void, std::mem::size_of::<OptimizeInfo>());
    }
}
#[cfg(not(windows))]
fn trim_heap() {}

fn heap_peak_mib() -> f64 {
    HEAP_PEAK.load(Ordering::Relaxed) as f64 / 1048576.0
}

pub static BOUND_CALLS: AtomicU64 = AtomicU64::new(0);
pub static STATE_STEPS: AtomicU64 = AtomicU64::new(0);
/// Diagnostics (DP_SLACK=1): window slack at row ends.
pub static SLACK: [AtomicU64; 8] = [const { AtomicU64::new(0) }; 8];
/// Diagnostics: nanoseconds in step loops, in finish(), and identity copies.
pub static PROF: [AtomicU64; 4] = [const { AtomicU64::new(0) }; 4];
thread_local! {
    // Four profile counters plus the bound-call count. Merge once per worker,
    // rather than contending on shared cache lines for each state.
    static LOCAL_PROF: Cell<[u64; 5]> = const { Cell::new([0; 5]) };
}

#[inline]
fn profile_add(k: usize, value: u64) {
    LOCAL_PROF.with(|local| {
        let mut counters = local.get();
        counters[k] += value;
        local.set(counters);
    });
}

fn flush_profile() {
    let counters = LOCAL_PROF.with(|local| local.replace([0; 5]));
    for k in 0..4 {
        if counters[k] != 0 {
            PROF[k].fetch_add(counters[k], Ordering::Relaxed);
        }
    }
    if counters[4] != 0 {
        BOUND_CALLS.fetch_add(counters[4], Ordering::Relaxed);
    }
}

/// Minimum number of windows of length d covering intervals [l, r] (r - l <= 1),
/// a lower bound on the number of distinct descendants (ancestors) after d
/// rows. `pts` must be sorted by r with l also nondecreasing (true for both
/// point lists built in `bound`). Greedy by right end is optimal; each window
/// jumps straight to the first point it does not cover, found by binary
/// search, so the cost is O(G log k) instead of O(k).
#[inline]
fn cover(pts: &[(i32, i32)], d: i32) -> u32 {
    let k = pts.len();
    let mut count = 0u32;
    let mut i = 0usize;
    while i < k {
        count += 1;
        let end = pts[i].1 + d;
        // first index > i with l > end
        let rest = &pts[i + 1..];
        i = i + 1 + rest.partition_point(|&(l, _)| l <= end);
    }
    count
}

#[cfg(test)]
fn cover_linear(pts: &[(i32, i32)], d: i32) -> u32 {
    let mut count = 0u32;
    let mut end = i32::MIN;
    for &(l, r) in pts {
        if l > end {
            count += 1;
            end = r + d;
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;
    fn params(n: usize, model: Model) -> Params {
        let tmax = if model == Model::Bond { n / 2 } else { (n + 1) / 2 };
        Params { n, model, tmax, moments: 3, verbose: false, profile: false, loose: false, mirror: true,
                 threads: 1, shards: 1, slack: false, nosign: false, cheap_lo: false, block: 8, park: true,
                 spill: None, spill_above: 0, batch: false, gpu: false, gpu_min: 0, key_lo: false }
    }

    #[test]
    fn mask_bound_matches_reference() {
        let mut seed = 987654321u64;
        let mut rnd = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            seed >> 33
        };
        for model in [Model::Bond, Model::Site] {
            for _ in 0..200000 {
                let n = 20 + (rnd() % 200) as usize;
                let p = params(n, model);
                let width = 1 + (rnd() % 60) as usize;
                let mut key: u128 = 0;
                for b in 0..width {
                    if rnd() % 3 == 0 {
                        key |= 1u128 << b;
                    }
                }
                if key == 0 {
                    continue;
                }
                let j = (rnd() % (width as u64 + 2)) as usize;
                let s = (rnd() % (p.tmax as u64)) as usize;
                assert_eq!(bound_inner(key, j, s, &p, None), bound_ref(key, j, s, &p), "key {key:b} j {j} s {s} n {n}");
            }
        }
    }

    #[test]
    fn cover_matches_linear() {
        let mut seed = 12345u64;
        for _ in 0..20000 {
            let mut pts = Vec::new();
            let mut pos = 0i32;
            let k = 1 + (seed % 20) as usize;
            for _ in 0..k {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                pos += 1 + ((seed >> 33) % 3) as i32;
                let fuzzy = (seed >> 40) & 1 == 1;
                pts.push(if fuzzy { (pos - 1, pos) } else { (pos, pos) });
            }
            for d in 0..30 {
                assert_eq!(cover(&pts, d), cover_linear(&pts, d), "{:?} d={}", pts, d);
            }
        }
    }
}

/// Bounds for a state at kink j of row s (transition s -> s+1).
/// Returns None if the state cannot contribute at order <= n or row <= tmax.
pub fn bound(key: u128, j: usize, s: usize, p: &Params) -> Option<(u32, u32)> {
    bound_with(key, j, s, p, None)
}

/// As `bound`, but with the lower end already known exactly (a target that
/// only one source produces inherits the source's trimmed lower end plus the
/// degree shift), which skips the past-cover loop.
pub fn bound_with(key: u128, j: usize, s: usize, p: &Params, lo_known: Option<u32>) -> Option<(u32, u32)> {
    if p.profile {
        let t0 = Instant::now();
        let r = bound_fast(key, j, s, p, lo_known);
        profile_add(3, t0.elapsed().as_nanos() as u64);
        profile_add(4, 1);
        return r;
    }
    prof::add_fine(prof::N_BOUND, 1);
    prof::time_fine(prof::BOUND, || bound_fast(key, j, s, p, lo_known))
}

/// Greedy window cover on bit masks. Points are intervals [l, l] (bit l of
/// `fixed`) or [l, l+1] (bit l of `fuzzy`). Equivalent to `cover`, but each
/// window costs a few u128 operations instead of a scan over the points.
#[inline(always)]
fn cover_mask(mut fixed: u128, mut fuzzy: u128, d: u32) -> u32 {
    let mut count = 0u32;
    loop {
        let r = fixed | (fuzzy << 1); // right ends
        if r == 0 {
            return count;
        }
        count += 1;
        let e = r.trailing_zeros() + d; // window [r0, r0 + d]
        if e >= 127 {
            return count;
        }
        // drop every point whose left end is <= e
        let keep = !((2u128 << e) - 1);
        fixed &= keep;
        fuzzy &= keep;
    }
}

/// Mask-based bound (see `bound_ref` for the point-list version it replaces).
/// Positions are shifted by +1 so that the fuzzy parent interval of column 0,
/// [-1, 0], starts at bit 0. With that shift:
///   past:   fixed = key bits >= j, fuzzy = key bits < j
///   future: fixed = (key bits <= j) << 1, fuzzy = key bits > j
#[inline]
fn bound_inner(key: u128, j: usize, s: usize, p: &Params, lo_known: Option<u32>) -> Option<(u32, u32)> {
    let below_j: u128 = if j >= 128 { !0 } else { (1u128 << j) - 1 };
    let upto_j: u128 = if j >= 127 { !0 } else { (2u128 << j) - 1 };
    let past_fixed = key & !below_j;
    let past_fuzzy = key & below_j;
    let fut_fixed = (key & upto_j) << 1;
    let fut_fuzzy = key & !upto_j;
    let nfixed = (key & upto_j).count_ones();
    let nunproc = (key & !upto_j).count_ones();

    let mut lo: u32 = nfixed;
    let mut d = 0u32;
    if let Some(l) = lo_known {
        lo = l;
        d = s as u32;
    } else if p.cheap_lo {
        lo += 2 * s as u32;
        d = s as u32;
    }
    while (d as usize) < s {
        let g = cover_mask(past_fixed, past_fuzzy, d);
        if g <= 2 {
            lo += 2 * (s as u32 - d);
            break;
        }
        lo += g;
        d += 1;
    }
    if lo as usize > p.n {
        return None;
    }
    let g0 = cover_mask(fut_fixed, fut_fuzzy, 0);
    let mut extra: u32 = match p.model {
        Model::Bond => nunproc,
        Model::Site => g0 - nfixed.min(g0),
    };
    let mut rmin = 0u32;
    loop {
        let g = if rmin == 0 { g0 } else { cover_mask(fut_fixed, fut_fuzzy, rmin) };
        if g <= 1 {
            break;
        }
        if p.model == Model::Bond || rmin >= 1 {
            extra += g;
        }
        rmin += 1;
    }
    if p.model == Model::Site && rmin >= 1 {
        extra += 1;
    }
    if s + 1 + rmin as usize > p.tmax {
        return None;
    }
    if p.loose {
        return Some((0, p.n as u32));
    }
    if (lo + extra) as usize > p.n {
        return None;
    }
    Some((lo, p.n as u32 - extra))
}

/// `cover_mask` on 64-bit masks; requires `fixed | fuzzy << 1` to fit in 64
/// bits. A window ending at e >= 63 covers every remaining point, where the
/// u128 version would see them dropped and return on the next pass.
#[inline(always)]
fn cover_mask64(mut fixed: u64, mut fuzzy: u64, d: u32) -> u32 {
    let mut count = 0u32;
    loop {
        let r = fixed | (fuzzy << 1);
        if r == 0 {
            return count;
        }
        count += 1;
        let e = r.trailing_zeros() + d;
        if e >= 63 {
            return count;
        }
        let keep = !((2u64 << e) - 1);
        fixed &= keep;
        fuzzy &= keep;
    }
}

/// Same result as `bound_inner`; keys below 2^63 take a 64-bit path.
#[inline]
fn bound_fast(key: u128, j: usize, s: usize, p: &Params, lo_known: Option<u32>) -> Option<(u32, u32)> {
    if key >> 63 != 0 {
        return bound_inner(key, j, s, p, lo_known);
    }
    let key = key as u64;
    let below_j: u64 = if j >= 64 { !0 } else { (1u64 << j) - 1 };
    let upto_j: u64 = if j >= 63 { !0 } else { (2u64 << j) - 1 };
    let past_fixed = key & !below_j;
    let past_fuzzy = key & below_j;
    let fut_fixed = (key & upto_j) << 1;
    let fut_fuzzy = key & !upto_j;
    let nfixed = (key & upto_j).count_ones();
    let nunproc = (key & !upto_j).count_ones();

    let mut lo: u32 = nfixed;
    let mut d = 0u32;
    if let Some(l) = lo_known {
        lo = l;
        d = s as u32;
    } else if p.cheap_lo {
        lo += 2 * s as u32;
        d = s as u32;
    }
    while (d as usize) < s {
        let g = cover_mask64(past_fixed, past_fuzzy, d);
        if g <= 2 {
            lo += 2 * (s as u32 - d);
            break;
        }
        lo += g;
        d += 1;
    }
    if lo as usize > p.n {
        return None;
    }
    let g0 = cover_mask64(fut_fixed, fut_fuzzy, 0);
    let mut extra: u32 = match p.model {
        Model::Bond => nunproc,
        Model::Site => g0 - nfixed.min(g0),
    };
    let mut rmin = 0u32;
    loop {
        let g = if rmin == 0 { g0 } else { cover_mask64(fut_fixed, fut_fuzzy, rmin) };
        if g <= 1 {
            break;
        }
        if p.model == Model::Bond || rmin >= 1 {
            extra += g;
        }
        rmin += 1;
    }
    if p.model == Model::Site && rmin >= 1 {
        extra += 1;
    }
    if s + 1 + rmin as usize > p.tmax {
        return None;
    }
    if p.loose {
        return Some((0, p.n as u32));
    }
    if (lo + extra) as usize > p.n {
        return None;
    }
    Some((lo, p.n as u32 - extra))
}

/// Reference implementation of the bound on explicit point lists (kept for tests).
#[allow(dead_code)]
fn bound_ref(key: u128, j: usize, s: usize, p: &Params) -> Option<(u32, u32)> {
    // future points: row s+1 positions (fixed) or fuzzy images of unprocessed vertices
    // past points: row s positions (fixed unprocessed, or fuzzy parents of fixed b's)
    let mut fut: [(i32, i32); 130] = [(0, 0); 130];
    let mut past: [(i32, i32); 130] = [(0, 0); 130];
    let mut nf = 0usize;
    let mut np = 0usize;
    let mut nfixed = 0u32;
    let mut nunproc = 0u32;
    let mut k = key;
    let mut pos = 0usize;
    while k != 0 {
        let tz = k.trailing_zeros() as usize;
        pos += tz;
        k >>= tz;
        let c = pos as i32;
        if pos < j {
            fut[nf] = (c, c);
            past[np] = (c - 1, c);
            nfixed += 1;
        } else if pos == j {
            fut[nf] = (c, c);
            past[np] = (c - 1, c - 1);
            nfixed += 1;
        } else {
            let c = c - 1; // row s column
            fut[nf] = (c, c + 1);
            past[np] = (c, c);
            nunproc += 1;
        }
        nf += 1;
        np += 1;
        k >>= 1;
        pos += 1;
    }
    let fut = &mut fut[..nf];
    let past = &mut past[..np];
    insertion_sort(fut);
    insertion_sort(past);

    // lower bound on degree already placed
    let mut lo: u32 = nfixed; // bonds into row s+1 so far (bond) / sites of row s+1 (site)
    let mut d = 0i32;
    while (d as usize) < s {
        let g = cover(past, d);
        if g <= 2 {
            lo += 2 * (s as u32 - d as u32);
            break;
        }
        lo += g;
        d += 1;
    }
    if lo as usize > p.n {
        return None;
    }

    // lower bound on degree still needed, and minimal extra rows
    let g0 = cover(fut, 0);
    let mut extra: u32 = match p.model {
        Model::Bond => nunproc,
        Model::Site => g0 - nfixed.min(g0),
    };
    let mut rmin = 0i32;
    loop {
        let g = if rmin == 0 { g0 } else { cover(fut, rmin) };
        if g <= 1 {
            break;
        }
        match p.model {
            Model::Bond => extra += g,
            Model::Site => {
                if rmin >= 1 {
                    extra += g
                }
            }
        }
        rmin += 1;
    }
    if p.model == Model::Site && rmin >= 1 {
        extra += 1; // the sink row s+1+rmin
    }
    if s + 1 + rmin as usize > p.tmax {
        return None;
    }
    if p.loose {
        return Some((0, p.n as u32));
    }
    if (lo + extra) as usize > p.n {
        return None;
    }
    Some((lo, p.n as u32 - extra))
}

#[inline]
fn insertion_sort(v: &mut [(i32, i32)]) {
    for i in 1..v.len() {
        let x = v[i];
        let mut k = i;
        while k > 0 && v[k - 1].1 > x.1 {
            v[k] = v[k - 1];
            k -= 1;
        }
        v[k] = x;
    }
}

/// Result of one run modulo one prime: non-nodal sums per row.
pub struct Raw {
    pub sn: Vec<Vec<C>>, // sn[t][deg]
    pub xn: Vec<Vec<C>>, // xn[t][deg]
}

type Layer = Vec<Table>;

fn layer_len(l: &Layer) -> usize {
    l.iter().map(|t| t.len()).sum()
}
fn layer_bytes(l: &Layer) -> usize {
    l.iter().map(|t| t.bytes()).sum()
}

/// Run `f(i)` for i in 0..count on `threads` threads; results in order.
fn par_map<T: Send, F: Fn(usize) -> T + Sync>(count: usize, threads: usize, f: F) -> Vec<T> {
    use std::sync::atomic::AtomicUsize;
    use std::sync::Mutex;
    if threads <= 1 || count <= 1 {
        return (0..count).map(f).collect();
    }
    let next = AtomicUsize::new(0);
    let out: Mutex<Vec<Option<T>>> = Mutex::new((0..count).map(|_| None).collect());
    std::thread::scope(|sc| {
        for _ in 0..threads.min(count) {
            sc.spawn(|| {
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= count {
                        break;
                    }
                    let r = f(i);
                    out.lock().unwrap()[i] = Some(r);
                }
                flush_profile();
                prof::flush();
            });
        }
    });
    out.into_inner().unwrap().into_iter().map(|x| x.unwrap()).collect()
}

/// Move the states of `layer` into pieces by target shard of their group key,
/// consuming the layer: each source shard is freed as soon as it has been
/// split, and each piece as soon as its thread has stepped it, so a parallel
/// step holds about one layer plus the pieces in flight instead of the whole
/// old layer next to the whole new one. Pieces are sized exactly (two passes).
fn split_layer<F: Fn(u128) -> u128 + Sync>(
    layer: Layer,
    nshards: usize,
    threads: usize,
    nm: usize,
    group: F,
) -> Vec<std::sync::Mutex<Vec<Table>>> {
    use std::sync::Mutex;
    let srcs: Vec<Mutex<Option<Table>>> = layer.into_iter().map(|t| Mutex::new(Some(t))).collect();
    let parts: Vec<Vec<Table>> = par_map(srcs.len(), threads, |si| {
        let _ts = prof::tsc();
        let t = srcs[si].lock().unwrap().take().unwrap();
        let mut dest: Vec<u32> = Vec::with_capacity(t.len());
        let mut nk = vec![0usize; nshards];
        let mut nc = vec![0usize; nshards];
        for i in 0..t.len() {
            let d = shard_of(group(t.key(i)), nshards);
            dest.push(d as u32);
            nk[d] += 1;
            nc[d] += t.block(i).len();
        }
        let mut v: Vec<Table> = (0..nshards).map(|d| Table::with_exact(nm, nk[d], nc[d])).collect();
        for i in 0..t.len() {
            let (lo, hi) = t.window(i);
            v[dest[i] as usize].push_copy(t.key(i), lo, hi, t.block(i));
        }
        prof::add(prof::SPLIT, prof::tsc() - _ts);
        v
    });
    let mut out: Vec<Vec<Table>> = (0..nshards).map(|_| Vec::with_capacity(parts.len())).collect();
    for part in parts {
        for (k, t) in part.into_iter().enumerate() {
            if t.len() > 0 {
                out[k].push(t);
            }
        }
    }
    out.into_iter().map(Mutex::new).collect()
}

/// Canonical member of {shape, mirror(shape)}.
/// Returns (canonical shape, input was mirrored, input is symmetric).
#[inline]
fn canonical(shape: u128) -> (u128, bool, bool) {
    let w = 127 - shape.leading_zeros();
    let rev = shape.reverse_bits() >> (127 - w);
    if rev == shape {
        (shape, false, true)
    } else if rev < shape {
        (rev, true, false)
    } else {
        (shape, false, false)
    }
}

/// a += mirror(v): offsets o -> L - o with L = t - w, so
/// M0 -> M0, M1 -> L M0 - M1, M2 -> L^2 M0 - 2 L M1 + M2.
#[allow(clippy::too_many_arguments)]
#[inline]
fn add_mirror(a: &mut [Vec<C>; 3], v: &[Vec<C>; 3], len: usize, l: i64, nm: usize, m: Mod) {
    let lm = m.from_i64(l);
    let l2 = m.from_i64(l * l);
    let two_l = m.from_i64(2 * l);
    for k in 0..len {
        a[0][k] = m.add(a[0][k], v[0][k]);
        if nm == 3 {
            a[1][k] = m.add(a[1][k], m.sub(m.mul(lm, v[0][k]), v[1][k]));
            a[2][k] = m.add(a[2][k], m.add(m.sub(m.mul(l2, v[0][k]), m.mul(two_l, v[1][k])), v[2][k]));
        }
    }
}

pub fn run(p: &Params, m: Mod) -> Raw {
    // Each prime reports its own diagnostics, including single-thread runs.
    for c in &PROF { c.store(0, Ordering::Relaxed); }
    BOUND_CALLS.store(0, Ordering::Relaxed);
    STATE_STEPS.store(0, Ordering::Relaxed);
    let n = p.n;
    let nm = p.moments;
    let nshards = p.shards;
    let threads = p.threads;
    let inv2 = (m.p + 1) / 2;
    let mirror = p.mirror;
    let mut sn = vec![vec![0 as C; n + 1]; p.tmax + 1];
    let mut xn = vec![vec![0 as C; n + 1]; p.tmax + 1];

    // row 0: origin, key bit 1 (column 0 of row 0, unprocessed), q = 0.
    // The origin is mirror-symmetric, so it is stored with weight 1/2.
    let mut layer: Layer = (0..nshards).map(|_| Table::new(nm, 16)).collect();
    {
        let key: u128 = 0b10;
        let (lo, hi) = bound(key, 0, 0, p).expect("origin must survive");
        let t = &mut layer[0];
        let idx = t.push_unique(key, lo, hi);
        t.coef_mut(idx, 0)[0] = if p.mirror { inv2 } else { 1 };
        t.finish();
    }
    HEAP_PEAK.store(HEAP_CUR.load(Ordering::Relaxed), Ordering::Relaxed);
    let t0 = Instant::now();
    let mut peak_states = 0usize;
    let mut peak_bytes = 0usize;
    for s in 0..p.tmax {
        let wmax = layer.iter().map(|t| t.max_width()).max().unwrap_or(0);
        let mut j = 0usize;
        if threads <= 1 && p.park {
            // Parked single-thread pass: every state waits in the table of the
            // next column that changes it, so states a column does not touch
            // are never copied (about half of all state-steps otherwise).
            let done = wmax + 2;
            let mut pend: Vec<Table> = (0..=done).map(|_| Table::new(nm, 16)).collect();
            for si in 0..layer.len() {
                let src = std::mem::replace(&mut layer[si], Table::new(nm, 0));
                for i in 0..src.len() {
                    let key = src.key(i);
                    let na = next_active(key, 0, done);
                    let (lo, hi) = src.window(i);
                    if let Some((idx, nlo, nhi)) = pend[na].get_or_insert(key, |_| Some((lo, hi))) {
                        add_block(&src, i, &mut pend[na], idx, nlo, nhi, 0, false, nm, m);
                    }
                }
            }
            for jj in 0..done {
                let mut t = std::mem::replace(&mut pend[jj], Table::new(nm, 0));
                t.finish();
                if p.profile { STATE_STEPS.fetch_add(t.len() as u64, Ordering::Relaxed); }
                for i in 0..t.len() {
                    step_parked(&t, i, &mut pend, jj, s, p, m, done);
                }
                let b: usize = pend.iter().map(|x| x.bytes()).sum::<usize>() + t.bytes();
                peak_bytes = peak_bytes.max(b);
                peak_states = peak_states.max(pend.iter().map(|x| x.len()).sum::<usize>());
            }
            let mut fin = std::mem::replace(&mut pend[done], Table::new(nm, 0));
            fin.finish();
            layer = vec![fin];
            j = wmax + 2; // skip the per-column loop below
        }
        while j <= wmax + 1 {
            let bj: u128 = 1 << j;
            let bj1: u128 = 1 << (j + 1);
            let mask = !(bj | bj1);
            if p.profile { STATE_STEPS.fetch_add(layer_len(&layer) as u64, Ordering::Relaxed); }
            if threads <= 1 {
                // streaming: free each source shard as soon as it is processed,
                // so the peak is about max(old, new) instead of old + new
                let est = layer_len(&layer) / nshards + 16;
                let ts = p.profile.then(Instant::now);
                let mut new: Layer = (0..nshards).map(|_| Table::new(nm, est)).collect();
                for si in 0..layer.len() {
                    let src = std::mem::replace(&mut layer[si], Table::new(nm, 0));
                    for i in 0..src.len() {
                        let sh = shard_of(src.key(i) & mask, nshards);
                        step_one(&src, i, &mut new[sh], j, s, p, m);
                    }
                    drop(src);
                    if si % 8 == 7 {
                        peak_bytes = peak_bytes.max(layer_bytes(&layer) + layer_bytes(&new));
                    }
                }
                if let Some(ts) = ts { profile_add(0, ts.elapsed().as_nanos() as u64); }
                let tf = p.profile.then(Instant::now);
                for t in new.iter_mut() {
                    t.finish();
                }
                if let Some(tf) = tf { profile_add(1, tf.elapsed().as_nanos() as u64); }
                peak_states = peak_states.max(layer_len(&new));
                layer = new;
                j += 1;
                continue;
            }
            // Parallel: a block of steps j..jend-1 rewrites only bits j..=jend,
            // so states that agree outside those bits never interact within the
            // block. Shard by the untouched bits and let each thread run the
            // whole block on its shard without synchronizing.
            let _ = mask;
            let jend = (j + p.block).min(wmax + 2);
            let touched: u128 = (if jend >= 127 { !0u128 } else { (2u128 << jend) - 1 }) & !((1u128 << j) - 1);
            let pieces = split_layer(std::mem::take(&mut layer), nshards, threads, nm, |k| k & !touched);
            let new: Layer = par_map(nshards, threads, |ti| {
                let _tb = prof::tsc();
                let list = std::mem::take(&mut *pieces[ti].lock().unwrap());
                let cnt: usize = list.iter().map(|t| t.len()).sum();
                let nc: usize = list.iter().map(|t| t.coefs()).sum();
                let mut cur = Table::new(nm, cnt + cnt / 2 + 16);
                cur.reserve_coefs(nc + nc / 8);
                {
                    let refs: Vec<&Table> = list.iter().collect();
                    step_all(&refs, &mut cur, j, s, p, m);
                }
                drop(list);
                cur.finish();
                for jj in (j + 1)..jend {
                    let mut nt = Table::new(nm, cur.len() + cur.len() / 2 + 16);
                    nt.reserve_coefs(cur.coefs() + cur.coefs() / 8);
                    step_all(&[&cur], &mut nt, jj, s, p, m);
                    nt.finish();
                    cur = nt;
                }
                prof::add(prof::BUSY, prof::tsc() - _tb);
                cur
            });
            peak_states = peak_states.max(layer_len(&new));
            peak_bytes = peak_bytes.max(layer_bytes(&new));
            layer = new;
            trim_heap();
            j = jend;
        }
        // row end: harvest singletons, normalize and canonicalize into row s+1
        let t = s + 1;
        let group = |k: u128| {
            let sh = k >> k.trailing_zeros();
            if mirror { canonical(sh).0 } else { sh }
        };
        let mut new: Layer = Vec::with_capacity(nshards);
        if threads <= 1 {
            let est = layer_len(&layer) / nshards + 16;
            let mut nts: Layer = (0..nshards).map(|_| Table::new(nm, est)).collect();
            let mut sc = Scratch::new(n);
            for si in 0..layer.len() {
                let src = std::mem::replace(&mut layer[si], Table::new(nm, 0));
                for i in 0..src.len() {
                    let sh = shard_of(group(src.key(i)), nshards);
                    row_end_one(&src, i, &mut nts[sh], &mut sc, t, p, m);
                }
            }
            for nt in nts.iter_mut() {
                nt.finish();
            }
            for d in 0..=n {
                sn[t][d] = m.add(sn[t][d], sc.sn[d]);
                xn[t][d] = m.add(xn[t][d], sc.xn[d]);
            }
            new = nts;
        } else {
            let pieces = split_layer(std::mem::take(&mut layer), nshards, threads, nm, group);
            let results: Vec<(Table, Scratch)> = par_map(nshards, threads, |ti| {
                let _tb = prof::tsc();
                let list = std::mem::take(&mut *pieces[ti].lock().unwrap());
                let cnt: usize = list.iter().map(|t| t.len()).sum();
                let nc: usize = list.iter().map(|t| t.coefs()).sum();
                let mut nt = Table::new(nm, cnt + 16);
                nt.reserve_coefs(nc + nc / 8);
                let mut sc = Scratch::new(n);
                prof::time(prof::ROWEND, || for src in list {
                    for i in 0..src.len() {
                        row_end_one(&src, i, &mut nt, &mut sc, t, p, m);
                    }
                });
                nt.finish();
                prof::add(prof::BUSY, prof::tsc() - _tb);
                (nt, sc)
            });
            for (nt, sc) in results {
                for d in 0..=n {
                    sn[t][d] = m.add(sn[t][d], sc.sn[d]);
                    xn[t][d] = m.add(xn[t][d], sc.xn[d]);
                }
                new.push(nt);
            }
        }
        layer = new;
        trim_heap();
        if p.verbose {
            eprintln!(
                "row {:3} -> states {:>11}  peak {:>11}  peak MiB {:>8.1}  heap peak MiB {:>8.1}  {:9.2}s",
                t,
                layer_len(&layer),
                peak_states,
                peak_bytes as f64 / 1048576.0,
                heap_peak_mib(),
                t0.elapsed().as_secs_f64()
            );
        }
    }
    flush_profile();
    if p.verbose {
        eprintln!(
            "done: peak states {} peak MiB {:.1} heap peak MiB {:.1} {:.2}s",
            peak_states,
            peak_bytes as f64 / 1048576.0,
            heap_peak_mib(),
            t0.elapsed().as_secs_f64()
        );
    }
    if p.profile {
        let g = |k: usize| PROF[k].load(Ordering::Relaxed);
        eprintln!(
            "profile: step loops {:.2}s  finish {:.2}s  identity copies {}  in bound {:.2}s  bound calls {}  state-steps {}",
            g(0) as f64 * 1e-9, g(1) as f64 * 1e-9, g(2), g(3) as f64 * 1e-9,
            BOUND_CALLS.load(Ordering::Relaxed), STATE_STEPS.load(Ordering::Relaxed)
        );
    }
    if p.slack {
        let g = |k: usize| SLACK[k].load(Ordering::Relaxed);
        eprintln!(
            "slack: states {} | lo slack 0:{} 1:{} 2:{} 3+:{} | all-zero {} | mean hi slack {:.2} of mean len {:.2}",
            g(0), g(1), g(2), g(3), g(4), g(7),
            g(5) as f64 / (g(0) - g(7)).max(1) as f64,
            g(6) as f64 / (g(0) - g(7)).max(1) as f64
        );
    }
    Raw { sn, xn }
}

/// One phase of the spilled sweep over row s: a block of column steps j..jend,
/// or the row end.
#[derive(Clone, Copy, PartialEq, Debug)]
enum Phase {
    Block(usize, usize),
    RowEnd,
}

/// Phases of row s whose layer has maximum width wmax (as `run`).
fn row_phases(wmax: usize, block: usize) -> Vec<Phase> {
    let mut v = Vec::new();
    let mut j = 0;
    while j <= wmax + 1 {
        let jend = (j + block).min(wmax + 2);
        v.push(Phase::Block(j, jend));
        j = jend;
    }
    v.push(Phase::RowEnd);
    v
}

fn touched_bits(j: usize, jend: usize) -> u128 {
    (if jend >= 127 { !0u128 } else { (2u128 << jend) - 1 }) & !((1u128 << j) - 1)
}

/// Shard grouping for the input of a phase: states in one group interact only
/// with each other. Before row s + 1's width is known, its first block is
/// grouped by the bits of the widest possible first block, 0..=block: a
/// coarser grouping, so still exact.
fn phase_group(ph: Phase, mirror: bool) -> impl Fn(u128) -> u128 {
    move |k: u128| match ph {
        Phase::Block(j, jend) => k & !touched_bits(j, jend),
        Phase::RowEnd => {
            let sh = k >> k.trailing_zeros();
            if mirror {
                canonical(sh).0
            } else {
                sh
            }
        }
    }
}

fn enum_ckpt_tag(p: &Params, prime: u64) -> String {
    format!(
        "dpenum model={:?} n={} prime={} moments={} mirror={} loose={} block={} shards={} coef_bits={}",
        p.model,
        p.n,
        prime,
        p.moments,
        p.mirror,
        p.loose,
        p.block,
        p.shards,
        8 * std::mem::size_of::<C>()
    )
}

fn enum_ckpt_path(p: &Params, prime: u64) -> std::path::PathBuf {
    p.spill.as_deref().unwrap().join(format!("dpenum-{:?}-n{}-p{}.ckpt", p.model, p.n, prime))
}

/// As `run` in parallel mode, but out of core: each phase (a block of steps,
/// or a row end) reads its input one target shard at a time and splits its
/// output by the next phase's grouping straight into the next generation,
/// which goes to disk once a generation exceeds `spill_above` and stays there
/// (see src/spill.rs). Memory holds only the shards in flight.
///
/// Pausing, as in dpwall: after every phase whose output is on disk a
/// checkpoint is written before the input generation is deleted. A file
/// `PAUSE` in the spill directory stops the run at the next phase boundary
/// (exit code 3); a kill or crash loses at most the phase in progress.
/// Rerunning the same command resumes.
pub fn run_spill(p: &Params, m: Mod, prime: u64) -> Raw {
    let n = p.n;
    let nm = p.moments;
    let ns = p.shards;
    let dir = p.spill.as_deref().unwrap();
    let inv2 = (m.p + 1) / 2;
    let t0 = Instant::now();
    let tag = enum_ckpt_tag(p, prime);
    let ckpt = enum_ckpt_path(p, prime);
    let gen_dir = |s: usize, pi: usize| dir.join(format!("dpenum-{:?}-n{}-p{}-r{}-{}", p.model, n, prime, s, pi));
    let mut sn = vec![vec![0 as C; n + 1]; p.tmax + 1];
    let mut xn = vec![vec![0 as C; n + 1]; p.tmax + 1];
    let mut stats = [0u64; 4]; // peak states, peak layer bytes, peak shard bytes, disk I/O bytes
    let (mut gen, s0, mut pi0, mut wmax, secs0) = match spill::load_text(&ckpt, &tag) {
        Some(kv) => {
            let num = |k: &str| -> u64 { kv[k].parse().unwrap() };
            let us = |k: &str| spill::split_nums(&kv[k]).into_iter().map(|x| x as usize).collect::<Vec<_>>();
            let (s, pi) = (num("row") as usize, num("phase") as usize);
            let d: std::path::PathBuf = kv["dir"].clone().into();
            let gen = spill::Gen::<u128>::reopen(d, pi, &us("nkeys"), &us("ncoefs"), num("bytes"));
            for t in 0..=p.tmax {
                for (k, v) in spill::split_nums(&kv[&format!("sn{t}")]).into_iter().enumerate() {
                    sn[t][k] = v as C;
                }
                for (k, v) in spill::split_nums(&kv[&format!("xn{t}")]).into_iter().enumerate() {
                    xn[t][k] = v as C;
                }
            }
            let st = spill::split_nums(&kv["stats"]);
            stats.copy_from_slice(&st);
            let secs: f64 = kv["secs"].parse().unwrap();
            eprintln!("resuming prime {} at row {} phase {} ({:.1}s done)", prime, s, pi, secs);
            (gen, s, pi, num("wmax") as usize, secs)
        }
        None => {
            // row 0: the origin, mirror-symmetric, stored with weight 1/2
            let mut first = Table::new(nm, 16);
            let key: u128 = 0b10;
            let (lo, hi) = bound(key, 0, 0, p).expect("origin must survive");
            let idx = first.push_unique(key, lo, hi);
            first.coef_mut(idx, 0)[0] = if p.mirror { inv2 } else { 1 };
            first.finish();
            let wmax = first.max_width();
            let gen = spill::Gen::<u128>::new(None, 0, ns);
            spill::write_split(&first, phase_group(row_phases(wmax, p.block)[0], p.mirror), &gen);
            (gen, 0, 0, wmax, 0.0)
        }
    };
    let [mut peak_states, mut peak_gen, mut peak_shard, mut io] = stats;
    let mut on_disk = gen.on_disk;
    HEAP_PEAK.store(HEAP_CUR.load(Ordering::Relaxed), Ordering::Relaxed);
    for s in s0..p.tmax {
        let phases = row_phases(wmax, p.block);
        let mut next_wmax = 0usize;
        for pi in pi0..phases.len() {
            let ph = phases[pi];
            let last_phase = s + 1 == p.tmax && ph == Phase::RowEnd;
            let next = if ph == Phase::RowEnd {
                // the first block of row s + 1, before its width is known
                Phase::Block(0, p.block)
            } else {
                phases[pi + 1]
            };
            let group = phase_group(next, p.mirror);
            on_disk |= gen.bytes() > p.spill_above;
            let (ns_next, pi_next) = if ph == Phase::RowEnd { (s + 1, 0) } else { (s, pi + 1) };
            let out: spill::Gen<u128> = if on_disk && !last_phase {
                spill::Gen::at(gen_dir(ns_next, pi_next), pi_next, ns)
            } else {
                spill::Gen::new(None, pi_next, ns)
            };
            let t = s + 1;
            prof::add(prof::N_PHASE, 1);
            let _tp = prof::tsc();
            let res: Vec<(usize, usize, usize, Option<Scratch>)> = par_map(ns, p.threads, |d| {
                let _tb = prof::tsc();
                let Some(src) = prof::time(prof::TAKE, || gen.take(d, nm).map(|t| t.into_indexed())) else {
                    return (0, 0, 0, None);
                };
                let (cnt, nc) = (src.len(), src.coefs());
                let (cur, sc, mem) = match ph {
                    Phase::Block(j, jend) => {
                        let mut cur = Table::new(nm, cnt + cnt / 2 + 16);
                        cur.reserve_coefs(nc + nc / 8);
                        step_all(&[&src], &mut cur, j, s, p, m);
                        cur.finish();
                        let mem = src.bytes() + cur.bytes();
                        drop(src);
                        for jj in (j + 1)..jend {
                            let mut nt = Table::new(nm, cur.len() + cur.len() / 2 + 16);
                            nt.reserve_coefs(cur.coefs() + cur.coefs() / 8);
                            step_all(&[&cur], &mut nt, jj, s, p, m);
                            nt.finish();
                            cur = nt;
                        }
                        (cur, None, mem)
                    }
                    Phase::RowEnd => {
                        let mut nt = Table::new(nm, cnt + 16);
                        nt.reserve_coefs(nc + nc / 8);
                        let mut sc = Scratch::new(n);
                        prof::time(prof::ROWEND, || for i in 0..src.len() {
                            row_end_one(&src, i, &mut nt, &mut sc, t, p, m);
                        });
                        nt.finish();
                        let mem = src.bytes() + nt.bytes();
                        (nt, Some(sc), mem)
                    }
                };
                if !last_phase {
                    prof::time(prof::WRITE, || spill::write_split(&cur, &group, &out));
                }
                prof::add(prof::BUSY, prof::tsc() - _tb);
                (cur.len(), mem.max(cur.bytes()), cur.max_width(), sc)
            });
            prof::add(prof::M_PAR, prof::tsc() - _tp);
            prof::time(prof::M_CLOSE, || out.close());
            if gen.on_disk {
                io += gen.bytes();
            }
            if out.on_disk {
                io += out.bytes();
            }
            let mut states = 0u64;
            for (st, mem, w, sc) in res {
                states += st as u64;
                peak_shard = peak_shard.max(mem as u64);
                if ph == Phase::RowEnd {
                    next_wmax = next_wmax.max(w); // the width of row s + 1
                }
                if let Some(sc) = sc {
                    for d in 0..=n {
                        sn[t][d] = m.add(sn[t][d], sc.sn[d]);
                        xn[t][d] = m.add(xn[t][d], sc.xn[d]);
                    }
                }
            }
            peak_states = peak_states.max(states);
            peak_gen = peak_gen.max(out.bytes());
            let secs = secs0 + t0.elapsed().as_secs_f64();
            if p.verbose && ph == Phase::RowEnd {
                eprintln!(
                    "row {:3} -> states {:>11} layer MiB {:8.1} {} peak states {:>11} peak layer MiB {:8.1} peak shard MiB {:6.1} heap peak MiB {:8.1} disk I/O GiB {:7.2} {:9.1}s",
                    t,
                    states,
                    out.bytes() as f64 / 1048576.0,
                    if out.on_disk { "disk" } else { "mem " },
                    peak_states,
                    peak_gen as f64 / 1048576.0,
                    peak_shard as f64 / 1048576.0,
                    heap_peak_mib(),
                    io as f64 / 1073741824.0,
                    secs
                );
            }
            let saved = out.on_disk;
            if saved {
                let (nkeys, ncoefs) = out.counts();
                let w = if ph == Phase::RowEnd { next_wmax } else { wmax };
                let mut text = format!("{tag}\nrow={ns_next}\nphase={pi_next}\nwmax={w}\n");
                text += &format!("dir={}\nbytes={}\n", out.directory().unwrap().display(), out.bytes());
                text += &format!("nkeys={}\nncoefs={}\n", spill::join(&nkeys), spill::join(&ncoefs));
                text += &format!("stats={}\nsecs={}\n", spill::join(&[peak_states, peak_gen, peak_shard, io]), secs);
                for tt in 0..=p.tmax {
                    text += &format!("sn{tt}={}\nxn{tt}={}\n", spill::join(&sn[tt]), spill::join(&xn[tt]));
                }
                prof::time(prof::M_CKPT, || spill::save_text(&ckpt, &text));
            }
            prof::time(prof::M_DROP, || gen = out); // drops the spent generation and deletes its files
            prof::time(prof::M_TRIM, trim_heap);
            if !last_phase && dir.join("PAUSE").exists() {
                if saved {
                    eprintln!("paused at row {} phase {}, {:.1}s: rerun the same command to resume", ns_next, pi_next, secs);
                } else {
                    eprintln!("paused before the layer reached disk: nothing saved, a rerun starts this prime over");
                }
                prof::report(p.threads);
                std::process::exit(3); // skips Drop, so the checkpointed generation stays
            }
        }
        wmax = next_wmax;
        pi0 = 0;
    }
    let _ = std::fs::remove_file(&ckpt);
    if p.verbose {
        eprintln!(
            "done: peak states {} peak layer MiB {:.1} peak shard MiB {:.1} heap peak MiB {:.1} disk I/O GiB {:.2} {:.2}s",
            peak_states,
            peak_gen as f64 / 1048576.0,
            peak_shard as f64 / 1048576.0,
            heap_peak_mib(),
            io as f64 / 1073741824.0,
            secs0 + t0.elapsed().as_secs_f64()
        );
    }
    Raw { sn, xn }
}

/// Per-thread scratch space and harvested non-nodal sums for one row end.
pub struct Scratch {
    sn: Vec<C>,
    xn: Vec<C>,
    v: [Vec<C>; 3],
    a: [Vec<C>; 3],
}

impl Scratch {
    fn new(n: usize) -> Self {
        let z = || vec![0 as C; n + 1];
        Scratch { sn: z(), xn: z(), v: [z(), z(), z()], a: [z(), z(), z()] }
    }
}

/// Row end for one state: shift to the shape frame, canonicalize under the
/// mirror, harvest singletons (non-nodal pieces ending at row t) and insert
/// the rest into the table for row t.
#[inline]
fn row_end_one(tab: &Table, i: usize, nt: &mut Table, sc: &mut Scratch, t: usize, p: &Params, m: Mod) {
    let nm = p.moments;
    let tt = t as i64;
    let inv2 = (m.p + 1) / 2;
    let key = tab.key(i);
    let (lo, hi) = tab.window(i);
    let len = (hi - lo + 1) as usize;
    let tz = key.trailing_zeros();
    let shape = key >> tz;
    let is_single = shape == 1;
    if p.slack {
        let c = tab.coef(i, 0);
        let first = c.iter().position(|&x| x != 0);
        let last = c.iter().rposition(|&x| x != 0);
        SLACK[0].fetch_add(1, Ordering::Relaxed);
        // slack of the trimmed window against the bound the state was created with
        let blo = bound(key, 127, t - 1, p).map(|b| b.0).unwrap_or(lo);
        match (first, last) {
            (Some(f), Some(l)) => {
                let f = f + (lo - blo.min(lo)) as usize;
                SLACK[1 + f.min(3)].fetch_add(1, Ordering::Relaxed);
                SLACK[5].fetch_add((len - 1 - l) as u64, Ordering::Relaxed);
                SLACK[6].fetch_add(len as u64, Ordering::Relaxed);
            }
            _ => {
                SLACK[7].fetch_add(1, Ordering::Relaxed);
            }
        }
    }
    if !is_single && t >= p.tmax {
        return;
    }
    let w = 127 - shape.leading_zeros() as i64;
    let v = &mut sc.v;
    let a = &mut sc.a;
    for mm in 0..nm {
        v[mm][..len].copy_from_slice(tab.coef(i, mm));
    }
    if nm == 3 && tz == 1 {
        for k in 0..len {
            // offset o -> o + 1
            v[2][k] = m.add(m.add(v[2][k], m.add(v[1][k], v[1][k])), v[0][k]);
            v[1][k] = m.add(v[1][k], v[0][k]);
        }
    }
    let (canon, mirrored, symmetric) = if p.mirror { canonical(shape) } else { (shape, false, false) };
    if mirrored {
        for mm in 0..nm {
            a[mm][..len].iter_mut().for_each(|x| *x = 0);
        }
        add_mirror(a, v, len, tt - w, nm, m);
    } else {
        for mm in 0..nm {
            a[mm][..len].copy_from_slice(&v[mm][..len]);
        }
        if symmetric {
            add_mirror(a, v, len, tt - w, nm, m);
            if !is_single {
                for mm in 0..nm {
                    for x in a[mm][..len].iter_mut() {
                        *x = m.mul(*x, inv2);
                    }
                }
            }
        }
    }
    if is_single {
        // full weight T of the singleton at column o: x = 2o - t
        let c1 = m.from_i64(-4 * tt);
        let c0 = m.from_i64(tt * tt);
        for k in 0..len {
            let deg = lo as usize + k;
            sc.sn[deg] = m.add(sc.sn[deg], a[0][k]);
            if nm == 3 {
                let x2 = m.add(m.add(m.mul(4, a[2][k]), m.mul(c1, a[1][k])), m.mul(c0, a[0][k]));
                sc.xn[deg] = m.add(sc.xn[deg], x2);
            }
        }
        return;
    }
    let nkey = canon << 1;
    let Some((idx, nlo, nhi)) = nt.get_or_insert(nkey, |k| bound(k, 0, t, p)) else {
        return;
    };
    let a_lo = lo.max(nlo);
    let a_hi = hi.min(nhi);
    if a_lo > a_hi {
        return;
    }
    let so = (a_lo - lo) as usize;
    let dof = (a_lo - nlo) as usize;
    let l = (a_hi - a_lo + 1) as usize;
    for mm in 0..nm {
        nt.accumulate(idx, mm, dof, &a[mm][so..so + l], false, m);
    }
}

/// First column >= k at which `key` changes (bit k or k+1 set), or `done`.
#[inline(always)]
fn next_active(key: u128, k: usize, done: usize) -> usize {
    if k >= 127 {
        return done;
    }
    let masked = (key >> k) << k;
    if masked == 0 {
        return done;
    }
    let m = masked.trailing_zeros() as usize;
    let ja = if m == 0 { 0 } else { (m - 1).max(k) };
    ja.min(done)
}

/// Add state i of `src`, shifted by e degrees (negated if `neg`), into entry
/// idx of `nt` whose window is [nlo, nhi].
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn add_block(src: &Table, i: usize, nt: &mut Table, idx: usize, nlo: u32, nhi: u32, e: u32, neg: bool, nm: usize, m: Mod) {
    let (lo, hi) = src.window(i);
    let a0 = (lo + e).max(nlo);
    let a1 = (hi + e).min(nhi);
    if a0 > a1 {
        return;
    }
    let so = (a0 - e - lo) as usize;
    let dof = (a0 - nlo) as usize;
    let len = (a1 - a0 + 1) as usize;
    for mm in 0..nm {
        nt.accumulate(idx, mm, dof, &src.coef(i, mm)[so..so + len], neg, m);
    }
}

/// Column j of row s for one state in the parked pass: targets go to the table
/// of their next active column. States here always have bit j or j+1 set.
#[allow(clippy::too_many_arguments)]
#[inline]
fn step_parked(tab: &Table, i: usize, pend: &mut [Table], j: usize, s: usize, p: &Params, m: Mod, done: usize) {
    let nm = p.moments;
    let bj: u128 = 1 << j;
    let bj1: u128 = 1 << (j + 1);
    let key = tab.key(i);
    let (lo, _hi) = tab.window(i);
    let q = key & bj != 0;
    let a = key & bj1 != 0;
    let base = key & !(bj | bj1);
    let mut tr: [(u128, u32, bool, bool); 3] = [(0, 0, false, false); 3];
    let ntr = match (q, a, p.model) {
        (true, false, _) => {
            tr[0] = (bj, 0, false, false);
            1
        }
        (false, true, _) => {
            tr[0] = (bj, 1, false, false);
            tr[1] = (bj1, 1, false, true);
            tr[2] = (bj | bj1, 2, false, false);
            3
        }
        (true, true, Model::Bond) => {
            tr[0] = (bj, 1, true, false);
            tr[1] = (bj | bj1, 1, false, false);
            tr[2] = (bj | bj1, 2, true, false);
            3
        }
        (true, true, Model::Site) if p.nosign => {
            tr[0] = (bj, 0, false, false);
            tr[1] = (bj | bj1, 1, false, false);
            tr[2] = (bj | bj1, 1, false, false);
            3
        }
        (true, true, Model::Site) => {
            tr[0] = (bj, 0, true, false);
            1
        }
        (false, false, _) => unreachable!("parked state inactive at its column"),
    };
    for &(bits, e, neg, unique) in &tr[..ntr] {
        let neg = neg && !p.nosign;
        let nkey = base | bits;
        let na = next_active(nkey, j + 1, done);
        let nt = &mut pend[na];
        let (idx, nlo, nhi) = if unique {
            match prof::time_fine(prof::BOUND_U, || bound_with(nkey, j + 1, s, p, Some(lo + e))) {
                Some((nlo, nhi)) => (nt.push_unique(nkey, nlo, nhi), nlo, nhi),
                None => continue,
            }
        } else {
            match nt.get_or_insert(nkey, |k| bound(k, j + 1, s, p)) {
                Some(x) => x,
                None => continue,
            }
        };
        add_block(tab, i, nt, idx, nlo, nhi, e, neg, nm, m);
    }
}

/// Targets of a state at column j: (target bits, degree shift, negate, unique),
/// or None for the identity (vertex absent, nothing pending). Same moves as
/// `step_one`.
#[inline(always)]
fn moves(key: u128, j: usize, p: &Params) -> Option<([(u128, u32, bool, bool); 3], usize)> {
    let bj: u128 = 1 << j;
    let bj1: u128 = 1 << (j + 1);
    let q = key & bj != 0;
    let a = key & bj1 != 0;
    let mut tr: [(u128, u32, bool, bool); 3] = [(0, 0, false, false); 3];
    let ntr = match (q, a, p.model) {
        (false, false, _) => return None,
        (true, false, _) => {
            tr[0] = (bj, 0, false, false);
            1
        }
        (false, true, _) => {
            tr[0] = (bj, 1, false, false);
            tr[1] = (bj1, 1, false, true);
            tr[2] = (bj | bj1, 2, false, false);
            3
        }
        (true, true, Model::Bond) => {
            tr[0] = (bj, 1, true, false);
            tr[1] = (bj | bj1, 1, false, false);
            tr[2] = (bj | bj1, 2, true, false);
            3
        }
        (true, true, Model::Site) if p.nosign => {
            tr[0] = (bj, 0, false, false);
            tr[1] = (bj | bj1, 1, false, false);
            tr[2] = (bj | bj1, 1, false, false);
            3
        }
        (true, true, Model::Site) => {
            tr[0] = (bj, 0, true, false);
            1
        }
    };
    let base = key & !(bj | bj1);
    for t in tr[..ntr].iter_mut() {
        t.0 |= base;
        t.2 = t.2 && !p.nosign;
    }
    Some((tr, ntr))
}

/// Bounds for a batch of new states at column j1 of row s: (lo, hi), or
/// (1, 0) for a dead state. `lo_known[k] >= 0` passes a known lower end
/// (`bound_with`). On the GPU with `--gpu`, unless the batch is small.
fn bound_batch(keys: &[u128], lo_known: &[i32], j1: usize, s: usize, p: &Params, out: &mut Vec<(u32, u32)>) {
    out.clear();
    #[cfg(feature = "gpu")]
    if p.gpu && keys.len() >= p.gpu_min {
        prof::add(prof::N_GPU, keys.len() as u64);
        prof::add(prof::N_GPU_CALLS, 1);
        prof::time(prof::GPU, || gpu::bound(keys, lo_known, j1, s, p, out));
        return;
    }
    dump_bounds(keys, lo_known, j1, s);
    let t = prof::tsc();
    out.extend(keys.iter().zip(lo_known).map(|(&k, &l)| {
        let lk = if l >= 0 { Some(l as u32) } else { None };
        bound_with(k, j1, s, p, lk).unwrap_or((1, 0))
    }));
    if !prof::FINE {
        prof::add(prof::BOUND, prof::tsc() - t);
        prof::add(prof::N_BOUND, keys.len() as u64);
    }
}

/// Development aid (DP_DUMP_BOUND=FILE, with --batch): append every
/// DP_DUMP_EVERY-th batch (default 64) of bound inputs as 24-byte records
/// (key u128, lo_known i32, j u16, s u16) for `--bench-bound`.
fn dump_bounds(keys: &[u128], lo_known: &[i32], j1: usize, s: usize) {
    use std::sync::{Mutex, OnceLock};
    static DUMP: OnceLock<Option<(Mutex<std::fs::File>, u64)>> = OnceLock::new();
    static SEEN: AtomicU64 = AtomicU64::new(0);
    let d = DUMP.get_or_init(|| {
        let path = std::env::var("DP_DUMP_BOUND").ok()?;
        let every = std::env::var("DP_DUMP_EVERY").ok().and_then(|v| v.parse().ok()).unwrap_or(64);
        Some((Mutex::new(std::fs::File::create(path).unwrap()), every))
    });
    let Some((f, every)) = d else { return };
    if SEEN.fetch_add(1, Ordering::Relaxed) % every != 0 {
        return;
    }
    let mut buf = Vec::with_capacity(keys.len() * 24);
    for (&k, &l) in keys.iter().zip(lo_known) {
        buf.extend_from_slice(&k.to_le_bytes());
        buf.extend_from_slice(&l.to_le_bytes());
        buf.extend_from_slice(&(j1 as u16).to_le_bytes());
        buf.extend_from_slice(&(s as u16).to_le_bytes());
    }
    f.lock().unwrap().write_all(&buf).unwrap();
}

/// `--bench-bound FILE`: replay dumped bound inputs through `bound_inner` and
/// `bound_fast`, check that they agree exactly, and time both.
fn bench_bound(path: &str, p: &Params) {
    let raw = std::fs::read(path).unwrap();
    let recs: Vec<(u128, Option<u32>, usize, usize)> = raw
        .chunks_exact(24)
        .map(|c| {
            let k = u128::from_le_bytes(c[0..16].try_into().unwrap());
            let l = i32::from_le_bytes(c[16..20].try_into().unwrap());
            let j = u16::from_le_bytes(c[20..22].try_into().unwrap()) as usize;
            let s = u16::from_le_bytes(c[22..24].try_into().unwrap()) as usize;
            (k, if l >= 0 { Some(l as u32) } else { None }, j, s)
        })
        .collect();
    let wide = recs.iter().filter(|r| r.0 >> 63 != 0).count();
    let known = recs.iter().filter(|r| r.1.is_some()).count();
    eprintln!("{} records, {} keys over 63 bits, {} with known lo", recs.len(), wide, known);
    let mut bad = 0usize;
    for &(k, l, j, s) in &recs {
        if bound_inner(k, j, s, p, l) != bound_fast(k, j, s, p, l) {
            bad += 1;
            if bad <= 5 {
                eprintln!("MISMATCH key {k:#x} j {j} s {s} lo_known {l:?}: {:?} vs {:?}", bound_inner(k, j, s, p, l), bound_fast(k, j, s, p, l));
            }
        }
    }
    eprintln!("mismatches: {bad}");
    for round in 0..3 {
        for (name, f) in [("inner", bound_inner as fn(u128, usize, usize, &Params, Option<u32>) -> Option<(u32, u32)>), ("fast", bound_fast)] {
            let t = Instant::now();
            let mut acc = 0u64;
            for &(k, l, j, s) in &recs {
                if let Some((a, b)) = f(std::hint::black_box(k), j, s, p, l) {
                    acc = acc.wrapping_add((a ^ b) as u64);
                }
            }
            let ns = t.elapsed().as_nanos() as f64 / recs.len() as f64;
            eprintln!("round {round} {name:6} {ns:7.1} ns/call  (check {acc})");
        }
        let t = Instant::now();
        let mut acc = 0u64;
        for &(k, _l, j, s) in &recs {
            if let Some((a, b)) = bound_fast(std::hint::black_box(k), j, s, p, Some(0)) {
                acc = acc.wrapping_add((a ^ b) as u64);
            }
        }
        eprintln!("round {round} future-only (fast) {:7.1} ns/call", t.elapsed().as_nanos() as f64 / recs.len() as f64);
    }
}

/// Column j of row s for every state of `srcs`, into `nt` (empty), in two
/// passes so that bounds are computed in one batch: pass 1 inserts every
/// target as a pending entry and collects the new ones; the batch gives their
/// windows and the arena is allocated; pass 2 adds the coefficients. Produces
/// the same table as `step_one` over the same sources (up to dead entries,
/// which `finish` drops).
fn step_batch(srcs: &[&Table], nt: &mut Table, j: usize, s: usize, p: &Params, m: Mod) {
    assert_eq!(nt.len(), 0);
    let nm = p.moments;
    let total: usize = srcs.iter().map(|t| t.len()).sum();
    let mut tgt: Vec<u32> = Vec::with_capacity(total * 2);
    let mut windows: Vec<(u32, u32)> = Vec::with_capacity(total * 2);
    let mut bkeys: Vec<u128> = Vec::with_capacity(total);
    let mut blo: Vec<i32> = Vec::with_capacity(total);
    let mut bpos: Vec<u32> = Vec::with_capacity(total);
    // Lower end of each new state: the least (source lower end + shift) over
    // its sources. Source windows are trimmed to their nonzero range, so each
    // source has a configuration at its lower end, and the target one at that
    // plus the shift: this is a lower bound on the target's orders, and never
    // below the interval-cover bound (which bounds every configuration), so it
    // replaces the cover's past loop. DP_KEY_LO=1 uses the cover instead.
    let key_lo = p.key_lo;
    let mut lomin: Vec<u32> = Vec::with_capacity(total * 2);
    // pass 1: targets
    for src in srcs {
        for i in 0..src.len() {
            let key = src.key(i);
            let (lo, hi) = src.window(i);
            let Some((tr, ntr)) = moves(key, j, p) else {
                if p.profile { profile_add(2, 1); }
                tgt.push(nt.pend_push(key) as u32);
                windows.push((lo, hi));
                lomin.push(lo);
                continue;
            };
            for &(nkey, e, _neg, unique) in &tr[..ntr] {
                let (idx, new) = if unique { (nt.pend_push(nkey), true) } else { nt.pend_insert(nkey) };
                if new {
                    bpos.push(windows.len() as u32);
                    windows.push((1, 0));
                    lomin.push(lo + e);
                    bkeys.push(nkey);
                    blo.push(if unique || !key_lo { -2 } else { -1 }); // -2: filled from lomin below
                } else {
                    lomin[idx] = lomin[idx].min(lo + e);
                }
                tgt.push(idx as u32);
            }
        }
    }
    for (l, &pos) in blo.iter_mut().zip(&bpos) {
        if *l == -2 {
            *l = lomin[pos as usize] as i32;
        }
    }
    drop(lomin);
    // windows of the new states, then the arena
    let mut bw = Vec::new();
    bound_batch(&bkeys, &blo, j + 1, s, p, &mut bw);
    drop(bkeys);
    drop(blo);
    for (&pos, &w) in bpos.iter().zip(&bw) {
        windows[pos as usize] = w;
    }
    drop(bpos);
    drop(bw);
    nt.alloc_pending(0, &windows);
    drop(windows);
    // pass 2: coefficients
    let mut c = 0usize;
    for src in srcs {
        for i in 0..src.len() {
            let key = src.key(i);
            let (lo, hi) = src.window(i);
            let Some((tr, ntr)) = moves(key, j, p) else {
                nt.put_block(tgt[c] as usize, src.block(i));
                c += 1;
                continue;
            };
            for &(_nkey, e, neg, _unique) in &tr[..ntr] {
                let idx = tgt[c] as usize;
                c += 1;
                let (nlo, nhi) = nt.window(idx);
                let a0 = (lo + e).max(nlo);
                let a1 = (hi + e).min(nhi);
                if a0 > a1 {
                    continue;
                }
                let src_off = (a0 - e - lo) as usize;
                let dst_off = (a0 - nlo) as usize;
                let len = (a1 - a0 + 1) as usize;
                for mm in 0..nm {
                    nt.accumulate(idx, mm, dst_off, &src.coef(i, mm)[src_off..src_off + len], neg, m);
                }
            }
        }
    }
}

/// Column j of row s for every state of `srcs` into `nt`: batched (`--batch`,
/// `--gpu`) or one state at a time.
fn step_all(srcs: &[&Table], nt: &mut Table, j: usize, s: usize, p: &Params, m: Mod) {
    prof::time(prof::STEP, || {
        if p.batch {
            step_batch(srcs, nt, j, s, p, m);
        } else {
            for src in srcs {
                for i in 0..src.len() {
                    step_one(src, i, nt, j, s, p, m);
                }
            }
        }
    })
}

/// Process column j of row s for one source state into the target table.
#[inline]
fn step_one(tab: &Table, i: usize, nt: &mut Table, j: usize, s: usize, p: &Params, m: Mod) {
    let nm = p.moments;
    let bj: u128 = 1 << j;
    let bj1: u128 = 1 << (j + 1);
    // (all-zero states were removed by Table::finish, so no zero check here)
    let key = tab.key(i);
    let (lo, hi) = tab.window(i);
    let q = key & bj != 0;
    let a = key & bj1 != 0;
    let base = key & !(bj | bj1);
    if !q && !a {
        // vertex absent and nothing pending: identity, same bound, unique target
        if p.profile { profile_add(2, 1); }
        nt.push_copy(key, lo, hi, tab.block(i));
        return;
    }
    // (target bits, degree shift, negate, unique)
    let mut tr: [(u128, u32, bool, bool); 3] = [(0, 0, false, false); 3];
    let ntr = match (q, a, p.model) {
        (true, false, _) => {
            tr[0] = (bj, 0, false, false);
            1
        }
        (false, true, _) => {
            // L: (1,0) p ; R: (0,1) p (unique) ; LR: (1,1) p^2
            tr[0] = (bj, 1, false, false);
            tr[1] = (bj1, 1, false, true);
            tr[2] = (bj | bj1, 2, false, false);
            3
        }
        (true, true, Model::Bond) => {
            // L: (1,0) -p ; R: (1,1) +p ; LR: (1,1) -p^2
            tr[0] = (bj, 1, true, false);
            tr[1] = (bj | bj1, 1, false, false);
            tr[2] = (bj | bj1, 2, true, false);
            3
        }
        (true, true, Model::Site) if p.nosign => {
            tr[0] = (bj, 0, false, false);
            tr[1] = (bj | bj1, 1, false, false);
            tr[2] = (bj | bj1, 1, false, false);
            3
        }
        (true, true, Model::Site) => {
            // L: (1,0) -1 ; R: (1,1) +p ; LR: (1,1) -p : the last two cancel
            tr[0] = (bj, 0, true, false);
            1
        }
        (false, false, _) => unreachable!(),
    };
    for &(bits, e, neg, unique) in &tr[..ntr] {
        let neg = neg && !p.nosign;
        let nkey = base | bits;
        let (idx, nlo, nhi) = if unique {
            match prof::time_fine(prof::BOUND_U, || bound_with(nkey, j + 1, s, p, Some(lo + e))) {
                Some((nlo, nhi)) => (nt.push_unique(nkey, nlo, nhi), nlo, nhi),
                None => continue,
            }
        } else {
            match nt.get_or_insert(nkey, |k| bound(k, j + 1, s, p)) {
                Some(x) => x,
                None => continue,
            }
        };
        let a0 = (lo + e).max(nlo);
        let a1 = (hi + e).min(nhi);
        if a0 > a1 {
            continue;
        }
        let src_off = (a0 - e - lo) as usize;
        let dst_off = (a0 - nlo) as usize;
        let len = (a1 - a0 + 1) as usize;
        for mm in 0..nm {
            let src = &tab.coef(i, mm)[src_off..src_off + len];
            nt.accumulate(idx, mm, dst_off, src, neg, m);
        }
    }
}

/// Assemble the full moment series modulo the prime (Jensen 1999, eqs. 6-9).
pub fn assemble(raw: &Raw, p: &Params, m: Mod) -> [Vec<C>; 4] {
    let n = p.n;
    let mut s_n = vec![0 as C; n + 1];
    let mut n01 = vec![0 as C; n + 1];
    let mut n02 = vec![0 as C; n + 1];
    let mut n20 = vec![0 as C; n + 1];
    for t in 1..=p.tmax {
        let tt = m.from_i64(t as i64);
        let t2 = m.from_i64((t * t) as i64);
        for d in 0..=n {
            let v = raw.sn[t][d];
            s_n[d] = m.add(s_n[d], v);
            n01[d] = m.add(n01[d], m.mul(tt, v));
            n02[d] = m.add(n02[d], m.mul(t2, v));
            n20[d] = m.add(n20[d], raw.xn[t][d]);
        }
    }
    let mut one_minus = vec![0 as C; n + 1];
    one_minus[0] = 1;
    for d in 1..=n {
        one_minus[d] = m.neg(s_n[d]);
    }
    let s = m.series_inv(&one_minus);
    let s2 = m.series_mul(&s, &s);
    let mu01 = m.series_mul(&n01, &s2);
    let tmp = m.series_mul(&m.series_mul(&n01, &n01), &s);
    let mut inner = n02.clone();
    for d in 0..=n {
        inner[d] = m.add(inner[d], m.add(tmp[d], tmp[d]));
    }
    let mu02 = m.series_mul(&inner, &s2);
    let mu20 = m.series_mul(&n20, &s2);
    [s, mu01, mu02, mu20]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut n = 20usize;
    let mut model = Model::Bond;
    let mut nprimes = 1usize;
    let mut first_prime = 0usize;
    let mut moments = 3usize;
    let mut out: Option<String> = None;
    let mut bench: Option<String> = None;
    let mut verbose = false;
    let mut loose = false;
    let mut mirror = true;
    let mut threads = 4usize;
    let mut block = 8usize;
    let mut park = false;
    let mut batch = true; // two-pass steps (multithreaded and spilled runs); --no-batch for the one-pass step
    let mut gpu = false;
    let gpu_min: usize = std::env::var("DP_GPU_MIN").ok().and_then(|v| v.parse().ok()).unwrap_or(4096);
    let mut shards: Option<usize> = None;
    let mut spill: Option<std::path::PathBuf> = None;
    let mut spill_above_mib = 256u64;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-n" => {
                n = args[i + 1].parse().unwrap();
                i += 1;
            }
            "--site" => model = Model::Site,
            "--bond" => model = Model::Bond,
            "--primes" => {
                nprimes = args[i + 1].parse().unwrap();
                i += 1;
            }
            "--first-prime" => {
                first_prime = args[i + 1].parse().unwrap();
                i += 1;
            }
            "--threads" => {
                threads = args[i + 1].parse().unwrap();
                i += 1;
            }
            "--no-x" => moments = 1,
            "--loose" => loose = true,
            "--no-mirror" => mirror = false,
            "--park" => park = true,
            "--batch" => batch = true,
            "--no-batch" => batch = false,
            "--gpu" => {
                if !cfg!(feature = "gpu") {
                    eprintln!("--gpu needs a build with --features gpu");
                    std::process::exit(2);
                }
                gpu = true;
                batch = true;
            }
            "--block" => {
                block = args[i + 1].parse().unwrap();
                i += 1;
            }
            "--shards" => {
                shards = Some(args[i + 1].parse().unwrap());
                i += 1;
            }
            "--spill" => {
                spill = Some(args[i + 1].clone().into());
                i += 1;
            }
            "--spill-above" => {
                spill_above_mib = args[i + 1].parse().unwrap();
                i += 1;
            }
            "-o" => {
                out = Some(args[i + 1].clone());
                i += 1;
            }
            "-v" => verbose = true,
            "--bench-bound" => {
                bench = Some(args[i + 1].clone());
                i += 1;
            }
            a => panic!("unknown argument {a}"),
        }
        i += 1;
    }
    assert!(n <= 250, "keys are u128: n <= 250");
    let tmax = match model {
        Model::Bond => n / 2,
        Model::Site => (n + 1) / 2,
    }
    .max(1);
    let shards = shards.unwrap_or(if spill.is_some() || threads <= 1 { 64 } else { threads * 8 });
    if let Some(d) = &spill {
        assert!(!park, "--park does not combine with --spill");
        std::fs::create_dir_all(d).expect("create spill directory");
        // starting the program is the request to run: clear an old pause
        if std::fs::remove_file(d.join("PAUSE")).is_ok() {
            eprintln!("removed {:?}", d.join("PAUSE"));
        }
    }
    let slack = std::env::var("DP_SLACK").is_ok();
    let nosign = std::env::var("DP_NOSIGN").is_ok();
    let cheap_lo = std::env::var("DP_CHEAP_LO").is_ok();
    let profile = std::env::var("DP_PROFILE").is_ok();
    let p = Params { n, model, tmax, moments, verbose, profile, loose, mirror, threads, shards, slack, nosign, cheap_lo, block, park, spill, spill_above: spill_above_mib << 20, batch, gpu, gpu_min, key_lo: std::env::var("DP_KEY_LO").is_ok() };
    if let Some(b) = bench {
        bench_bound(&b, &p);
        return;
    }
    let primes = modp::primes(first_prime + nprimes);
    // spill mode: primes already in the output file are done, so rerunning
    // the same command resumes a queue
    let done = if p.spill.is_some() { out.as_ref().and_then(|o| std::fs::read_to_string(o).ok()).unwrap_or_default() } else { String::new() };
    let mname = if model == Model::Bond { "bond" } else { "site" };
    for &pr in &primes[first_prime..] {
        let tagline = format!("\"model\":\"{}\",\"n\":{},\"prime\":{},\"moments\":{},", mname, n, pr, moments);
        if done.lines().any(|l| l.contains(&tagline)) {
            eprintln!("prime {} already in the output: skipped", pr);
            continue;
        }
        let m = Mod::new(pr);
        prof::start();
        let t0 = Instant::now();
        let raw = if p.spill.is_some() { run_spill(&p, m, pr) } else { run(&p, m) };
        let series = assemble(&raw, &p, m);
        let names = ["S", "mu01", "mu02", "mu20"];
        let mut line = format!(
            "{{\"model\":\"{}\",\"n\":{},\"prime\":{},\"moments\":{},\"seconds\":{:.3}",
            if model == Model::Bond { "bond" } else { "site" },
            n,
            pr,
            moments,
            t0.elapsed().as_secs_f64()
        );
        for (k, name) in names.iter().enumerate() {
            if moments == 1 && *name == "mu20" {
                continue;
            }
            let v: Vec<String> = series[k].iter().map(|x| x.to_string()).collect();
            line.push_str(&format!(",\"{}\":[{}]", name, v.join(",")));
        }
        line.push('}');
        eprintln!("prime {} done in {:.2}s", pr, t0.elapsed().as_secs_f64());
        prof::report(threads);
        #[cfg(feature = "gpu")]
        if p.gpu && p.verbose {
            let (h2d, k, d2h) = gpu::stats();
            eprintln!("gpu: stream time h2d {:.2}s kernel {:.2}s d2h {:.2}s (summed over streams)", h2d, k, d2h);
        }
        // write each prime as soon as it is done so long runs can be resumed
        match &out {
            Some(path) => {
                let mut f = std::fs::OpenOptions::new().create(true).append(true).open(path).unwrap();
                writeln!(f, "{}", line).unwrap();
            }
            None => println!("{}", line),
        }
    }
}
