//! Cycle-counter profile (`--features prof`): per-thread rdtsc totals by stage,
//! merged once per worker. Compiled to nothing without the feature.

use std::cell::Cell;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;
use std::time::Instant;

pub const ON: bool = cfg!(feature = "prof");
/// Per-call timers (hash, accumulate, push, single bounds): `prof-fine`. They
/// cost about 25% of run time, so stage totals are distorted with them on.
pub const FINE: bool = cfg!(feature = "prof-fine");

pub const TAKE: usize = 0; // read a shard (disk or memory) and index it
pub const STEP: usize = 1; // step loops, everything inside
pub const FINISH: usize = 2; // Table::finish / finish_seq
pub const WRITE: usize = 3; // write_split (disk or memory)
pub const ROWEND: usize = 4; // row-end loop (harvest, canonicalize)
pub const BOUND: usize = 5; // bound / window computations, all callers
pub const GETINS: usize = 6; // get_or_insert, including the bound for new keys
pub const ACC: usize = 7; // accumulate (shifted mod-p adds)
pub const PUSH: usize = 8; // push_unique / push_copy
pub const BUSY: usize = 9; // whole worker task
pub const SPLIT: usize = 10; // in-memory split_layer
pub const N_GETINS: usize = 11; // counts
pub const N_ACC: usize = 12;
pub const N_ACC_ELEMS: usize = 13;
pub const N_BOUND: usize = 14;
pub const N_PUSH: usize = 15;
pub const BOUND_U: usize = 16; // bound for unique (unhashed) targets
pub const M_PAR: usize = 17; // main thread: wall time inside par_map
pub const M_CLOSE: usize = 18; // main thread: flush/close output files
pub const M_CKPT: usize = 19; // main thread: checkpoint text
pub const M_DROP: usize = 20; // main thread: drop spent generation (delete files)
pub const M_TRIM: usize = 21;
pub const N_PHASE: usize = 22;
pub const GPU: usize = 23; // waiting on GPU bound batches
pub const N_GPU: usize = 24; // states bounded on the GPU
pub const N_GPU_CALLS: usize = 25;
pub const W_DEST: usize = 26; // write_split: destination shards and permutation
pub const W_GATHER: usize = 27; // write_split: gather each destination's chunk
pub const W_OUT: usize = 28; // write_split: Gen::write (lock, buffer, syscalls)
pub const T_READ: usize = 29; // take: open and read a shard file
pub const T_INDEX: usize = 30; // take: rebuild offsets (into_indexed)
const NB: usize = 31;
const NAMES: [&str; 11] = ["take", "step", "finish", "write", "rowend", "bound", "get_or_insert", "accumulate", "push", "busy", "split"];

static TOTAL: [AtomicU64; NB] = [const { AtomicU64::new(0) }; NB];
thread_local! {
    static LOCAL: Cell<[u64; NB]> = const { Cell::new([0; NB]) };
}
static START: OnceLock<(Instant, u64)> = OnceLock::new();

#[inline(always)]
pub fn tsc() -> u64 {
    if ON {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            return core::arch::x86_64::_rdtsc();
        }
    }
    #[allow(unreachable_code)]
    0
}

#[inline(always)]
pub fn add(k: usize, v: u64) {
    if ON {
        LOCAL.with(|l| {
            let mut c = l.get();
            c[k] += v;
            l.set(c);
        });
    }
}

/// Per-call counterpart of `add` (only with `prof-fine`).
#[inline(always)]
pub fn add_fine(k: usize, v: u64) {
    if FINE {
        add(k, v);
    }
}

/// Per-call counterpart of `time` (only with `prof-fine`).
#[inline(always)]
pub fn time_fine<T>(k: usize, f: impl FnOnce() -> T) -> T {
    if FINE { time(k, f) } else { f() }
}

/// Time `f` into bucket k.
#[inline(always)]
pub fn time<T>(k: usize, f: impl FnOnce() -> T) -> T {
    if !ON {
        return f();
    }
    let t = tsc();
    let r = f();
    add(k, tsc() - t);
    r
}

pub fn start() {
    if ON {
        START.get_or_init(|| (Instant::now(), tsc()));
    }
}

pub fn flush() {
    if ON {
        let c = LOCAL.with(|l| l.replace([0; NB]));
        for k in 0..NB {
            if c[k] != 0 {
                TOTAL[k].fetch_add(c[k], Ordering::Relaxed);
            }
        }
    }
}

/// Print totals as thread-seconds and as a share of threads x wall time.
pub fn report(threads: usize) {
    if !ON {
        return;
    }
    flush();
    let Some(&(t0, c0)) = START.get() else { return };
    let wall = t0.elapsed().as_secs_f64();
    let hz = (tsc() - c0) as f64 / wall;
    let cap = wall * threads as f64;
    let g = |k: usize| TOTAL[k].load(Ordering::Relaxed);
    eprintln!("prof: wall {:.2}s, {} threads, tsc {:.3} GHz", wall, threads, hz / 1e9);
    for (k, name) in NAMES.iter().enumerate() {
        let s = TOTAL[k].load(Ordering::Relaxed) as f64 / hz;
        eprintln!("prof: {:14} {:10.2} thread-s  {:5.1}% of capacity", name, s, 100.0 * s / cap);
    }
    eprintln!("prof: {:14} {:10.2} thread-s  (subset of bound)", "bound(unique)", g(BOUND_U) as f64 / hz);
    let hash = (g(GETINS) as f64 - (g(BOUND) - g(BOUND_U)) as f64) / hz;
    eprintln!("prof: {:14} {:10.2} thread-s  {:5.1}% (get_or_insert minus its bounds)", "hash probe", hash, 100.0 * hash / cap);
    for (k, name) in [(M_PAR, "par_map"), (M_CLOSE, "close"), (M_CKPT, "checkpoint"), (M_DROP, "drop gen"), (M_TRIM, "trim heap")] {
        let v = g(k) as f64 / hz;
        eprintln!("prof: main {:10} {:8.2} s  {:5.1}% of wall", name, v, 100.0 * v / wall);
    }
    eprintln!("prof: phases {}", g(N_PHASE));
    for (k, name) in [(W_DEST, "write: dest"), (W_GATHER, "write: gather"), (W_OUT, "write: out"), (T_READ, "take: read"), (T_INDEX, "take: index")] {
        let v = g(k) as f64 / hz;
        eprintln!("prof: {:14} {:10.2} thread-s  {:5.1}% of capacity", name, v, 100.0 * v / cap);
    }
    let gs = g(GPU) as f64 / hz;
    eprintln!("prof: gpu wait {:.2} thread-s, {} batches, {} states ({:.1} M/thread-s, {:.0} per batch)", gs, g(N_GPU_CALLS), g(N_GPU), g(N_GPU) as f64 / gs.max(1e-9) / 1e6, g(N_GPU) as f64 / g(N_GPU_CALLS).max(1) as f64);
    eprintln!(
        "prof: calls get_or_insert {}  accumulate {} ({} elems, {:.1}/call)  bound {}  push {}",
        g(N_GETINS),
        g(N_ACC),
        g(N_ACC_ELEMS),
        g(N_ACC_ELEMS) as f64 / g(N_ACC).max(1) as f64,
        g(N_BOUND),
        g(N_PUSH)
    );
}
