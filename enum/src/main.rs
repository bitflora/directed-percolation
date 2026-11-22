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

use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

mod modp;
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
    pub loose: bool, // testing only: disable order pruning (keep only the row limit)
    pub mirror: bool,
    pub threads: usize,
    pub shards: usize,
    pub slack: bool,
    pub nosign: bool, // diagnostics only: count graphs without signs (no cancellation)
}

pub static BOUND_CALLS: AtomicU64 = AtomicU64::new(0);
pub static STATE_STEPS: AtomicU64 = AtomicU64::new(0);
/// Diagnostics (DP_SLACK=1): window slack at row ends.
pub static SLACK: [AtomicU64; 8] = [const { AtomicU64::new(0) }; 8];
/// Diagnostics: nanoseconds in step loops, in finish(), and identity copies.
pub static PROF: [AtomicU64; 4] = [const { AtomicU64::new(0) }; 4];

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
        Params { n, model, tmax, moments: 3, verbose: false, loose: false, mirror: true,
                 threads: 1, shards: 1, slack: false, nosign: false }
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
                assert_eq!(bound_inner(key, j, s, &p), bound_ref(key, j, s, &p), "key {key:b} j {j} s {s} n {n}");
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
    if p.verbose {
        let t0 = Instant::now();
        let r = bound_inner(key, j, s, p);
        PROF[3].fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
        return r;
    }
    bound_inner(key, j, s, p)
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
fn bound_inner(key: u128, j: usize, s: usize, p: &Params) -> Option<(u32, u32)> {
    BOUND_CALLS.fetch_add(1, Ordering::Relaxed);
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
            sc.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= count {
                    break;
                }
                let r = f(i);
                out.lock().unwrap()[i] = Some(r);
            });
        }
    });
    out.into_inner().unwrap().into_iter().map(|x| x.unwrap()).collect()
}

/// For each target shard, the list of (source shard, index) whose group key
/// hashes to it.
fn bucket<F: Fn(u128) -> u128 + Sync>(
    layer: &Layer,
    nshards: usize,
    threads: usize,
    group: F,
) -> Vec<Vec<(u32, u32)>> {
    if nshards == 1 {
        let mut v = Vec::with_capacity(layer_len(layer));
        for (si, t) in layer.iter().enumerate() {
            for i in 0..t.len() {
                v.push((si as u32, i as u32));
            }
        }
        return vec![v];
    }
    let parts: Vec<Vec<Vec<(u32, u32)>>> = par_map(layer.len(), threads, |si| {
        let t = &layer[si];
        let mut v: Vec<Vec<(u32, u32)>> = (0..nshards).map(|_| Vec::new()).collect();
        for i in 0..t.len() {
            let g = group(t.key(i));
            v[shard_of(g, nshards)].push((si as u32, i as u32));
        }
        v
    });
    let mut out: Vec<Vec<(u32, u32)>> = (0..nshards).map(|_| Vec::new()).collect();
    for part in parts {
        for (k, v) in part.into_iter().enumerate() {
            if out[k].is_empty() {
                out[k] = v;
            } else {
                out[k].extend_from_slice(&v);
            }
        }
    }
    out
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
    let t0 = Instant::now();
    let mut peak_states = 0usize;
    let mut peak_bytes = 0usize;
    for s in 0..p.tmax {
        let wmax = layer.iter().map(|t| t.max_width()).max().unwrap_or(0);
        for j in 0..=(wmax + 1) {
            let bj: u128 = 1 << j;
            let bj1: u128 = 1 << (j + 1);
            let mask = !(bj | bj1);
            STATE_STEPS.fetch_add(layer_len(&layer) as u64, Ordering::Relaxed);
            if threads <= 1 {
                // streaming: free each source shard as soon as it is processed,
                // so the peak is about max(old, new) instead of old + new
                let est = layer_len(&layer) / nshards + 16;
                let ts = Instant::now();
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
                PROF[0].fetch_add(ts.elapsed().as_nanos() as u64, Ordering::Relaxed);
                let tf = Instant::now();
                for t in new.iter_mut() {
                    t.finish();
                }
                PROF[1].fetch_add(tf.elapsed().as_nanos() as u64, Ordering::Relaxed);
                peak_states = peak_states.max(layer_len(&new));
                layer = new;
                continue;
            }
            let buckets = bucket(&layer, nshards, threads, |k| k & mask);
            let src = &layer;
            let new: Layer = par_map(nshards, threads, |ti| {
                let list = &buckets[ti];
                let mut nt = Table::new(nm, list.len() + list.len() / 2 + 16);
                for &(si, i) in list {
                    step_one(&src[si as usize], i as usize, &mut nt, j, s, p, m);
                }
                nt.finish();
                nt
            });
            peak_states = peak_states.max(layer_len(&layer) + layer_len(&new));
            peak_bytes = peak_bytes.max(layer_bytes(&layer) + layer_bytes(&new));
            drop(buckets);
            layer = new;
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
            let buckets = bucket(&layer, nshards, threads, group);
            let src = &layer;
            let results: Vec<(Table, Scratch)> = par_map(nshards, threads, |ti| {
                let list = &buckets[ti];
                let mut nt = Table::new(nm, list.len() + 16);
                let mut sc = Scratch::new(n);
                for &(si, i) in list {
                    row_end_one(&src[si as usize], i as usize, &mut nt, &mut sc, t, p, m);
                }
                nt.finish();
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
        if p.verbose {
            eprintln!(
                "row {:3} -> states {:>11}  peak {:>11}  peak MiB {:>8.1}  {:9.2}s",
                t,
                layer_len(&layer),
                peak_states,
                peak_bytes as f64 / 1048576.0,
                t0.elapsed().as_secs_f64()
            );
        }
    }
    if p.verbose {
        eprintln!(
            "done: peak states {} peak MiB {:.1} {:.2}s  bound calls {}  state-steps {}",
            peak_states,
            peak_bytes as f64 / 1048576.0,
            t0.elapsed().as_secs_f64(),
            BOUND_CALLS.load(Ordering::Relaxed),
            STATE_STEPS.load(Ordering::Relaxed)
        );
    }
    if p.verbose {
        let g = |k: usize| PROF[k].load(Ordering::Relaxed);
        eprintln!(
            "profile: step loops {:.2}s  finish {:.2}s  identity copies {}  in bound {:.2}s",
            g(0) as f64 * 1e-9, g(1) as f64 * 1e-9, g(2), g(3) as f64 * 1e-9
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
    if tab.is_zero(i) {
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

/// Process column j of row s for one source state into the target table.
#[inline]
fn step_one(tab: &Table, i: usize, nt: &mut Table, j: usize, s: usize, p: &Params, m: Mod) {
    let nm = p.moments;
    let bj: u128 = 1 << j;
    let bj1: u128 = 1 << (j + 1);
    if tab.is_zero(i) {
        return; // exact cancellation (or zero mod p): contributes nothing
    }
    let key = tab.key(i);
    let (lo, hi) = tab.window(i);
    let q = key & bj != 0;
    let a = key & bj1 != 0;
    let base = key & !(bj | bj1);
    if !q && !a {
        // vertex absent and nothing pending: identity, same bound, unique target
        PROF[2].fetch_add(1, Ordering::Relaxed);
        let idx = nt.push_unique(key, lo, hi);
        for mm in 0..nm {
            nt.coef_mut(idx, mm).copy_from_slice(tab.coef(i, mm));
        }
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
            match bound(nkey, j + 1, s, p) {
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
    let mut verbose = false;
    let mut loose = false;
    let mut mirror = true;
    let mut threads = 4usize;
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
            "-o" => {
                out = Some(args[i + 1].clone());
                i += 1;
            }
            "-v" => verbose = true,
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
    let shards = if threads <= 1 { 64 } else { threads * 8 };
    let slack = std::env::var("DP_SLACK").is_ok();
    let nosign = std::env::var("DP_NOSIGN").is_ok();
    let p = Params { n, model, tmax, moments, verbose, loose, mirror, threads, shards, slack, nosign };
    let primes = modp::primes(first_prime + nprimes);
    for &pr in &primes[first_prime..] {
        let m = Mod::new(pr);
        let t0 = Instant::now();
        let raw = run(&p, m);
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
