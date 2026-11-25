//! dpwall: low-density series for bond directed percolation next to an
//! impenetrable wall (Essam et al. 1996; Jensen 1999 II), modulo primes.
//!
//! Geometry: sites (t, x), t + x even, bonds (t, x) -> (t+1, x -+ 1); the source
//! is (0, 0) and bonds into x < 0 are absent. Series, coefficients of p^0..p^n:
//!   T   = sum_{t>=1} P(t), P(t) = Prob(the cluster reaches row t)  (tau_1)
//!   S   = sum C(t,x)                                                (gamma_1)
//!   mx1 = sum x C, mx2 = sum x^2 C                (gamma_1 + k nu_perp)
//!   mt1 = sum t C, mt2 = sum t^2 C                (gamma_1 + k nu_par)
//!
//! Method: the path-union expansion. A state is the set of boundary sites that
//! paths pass through; its weight is the signed sum over unions of directed
//! paths, p^(bonds) times -1 for every site entered by two bonds. Every path
//! continues to the current row, so C(t, x) is the weight of the state {x} on
//! row t, and P(t) is the sum over row-t states B of (-1)^(|B|-1) weight(B)
//! (inclusion-exclusion over the wet sites of row t). Checked against the
//! exact Domany-Kinzel transfer matrix in analysis/wall_ground_truth.py.
//!
//! Unlike C(t, x), the survival probability has no requirement that paths
//! merge again, so the collapse pruning of the bulk engine does not apply and
//! the wall breaks the translation reduction; only the order bound prunes.
//!
//! Sweep: row t -> t+1 one parent at a time. Row t has parents j = 0..=t/2 at
//! x = 2j + (t mod 2); with e = 1 for even t and 0 for odd t, parent j has
//! children j - e (x - 1) and j - e + 1 (x + 1) on row t+1, and child -1 is
//! the wall. Key layout at the start of step j: finished children c < j - e at
//! bit c + 1, the pending child j - e at bit P - 1, parent q >= j at bit
//! q - e + 2, where P = j - e + 2. Step j only reads and writes bits P - 1 and
//! P (after it, bit P holds the pending child j - e + 1), so states that agree
//! outside those bits form independent groups: blocks of steps run in
//! parallel on shards keyed by the untouched bits.

#![allow(dead_code)]

use std::io::Write;
use std::sync::atomic::Ordering;
use std::time::Instant;

mod modp;
mod prof;
mod spill;
mod table;

use modp::{Mod, C};
use table::{shard_of, Table};

const NAMES: [&str; 6] = ["T", "S", "mx1", "mx2", "mt1", "mt2"];

struct Params {
    n: usize,
    threads: usize,
    shards: usize,
    block: usize,
    verbose: bool,
    spill: Option<std::path::PathBuf>,
    spill_above: u64, // bytes: generations larger than this go to disk
}

/// Row t uses key bits up to t/2 + 3: u64 keys allow n <= 121, and the
/// `wide` feature (u128) n <= 249 (checked in main).
#[cfg(not(feature = "wide"))]
type Key = u64;
#[cfg(feature = "wide")]
type Key = u128;
type Layer = Vec<Table<Key>>;

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
                prof::flush();
            });
        }
    });
    out.into_inner().unwrap().into_iter().map(|x| x.unwrap()).collect()
}

/// Move the states into per-target-shard pieces (grouped by `group(key)`),
/// freeing each source shard once it is split.
fn split_layer<F: Fn(Key) -> Key + Sync>(layer: Layer, nshards: usize, threads: usize, group: F) -> Vec<std::sync::Mutex<Layer>> {
    use std::sync::Mutex;
    let srcs: Vec<Mutex<Option<Table<Key>>>> = layer.into_iter().map(|t| Mutex::new(Some(t))).collect();
    let parts: Vec<Layer> = par_map(srcs.len(), threads, |si| {
        let t = srcs[si].lock().unwrap().take().unwrap();
        let mut dest: Vec<u32> = Vec::with_capacity(t.len());
        let mut nk = vec![0usize; nshards];
        let mut nc = vec![0usize; nshards];
        for (key, _, _, b) in t.iter() {
            let d = shard_of(group(key), nshards);
            dest.push(d as u32);
            nk[d] += 1;
            nc[d] += b.len();
        }
        let mut v: Layer = (0..nshards).map(|d| Table::with_exact_seq(1, nk[d], nc[d])).collect();
        for (i, (key, lo, hi, b)) in t.iter().enumerate() {
            v[dest[i] as usize].push_copy(key, lo, hi, b);
        }
        v
    });
    let mut out: Vec<Layer> = (0..nshards).map(|_| Vec::with_capacity(parts.len())).collect();
    for part in parts {
        for (k, t) in part.into_iter().enumerate() {
            if t.len() > 0 {
                out[k].push(t);
            }
        }
    }
    out.into_iter().map(Mutex::new).collect()
}

/// Order window of a state on the boundary between rows t and t+1 whose
/// unprocessed parents are the bits at or above `above`.
/// lo: the union of paths to k distinct boundary sites has >= t + k - 1 bonds.
/// hi: every unprocessed occupied parent still needs at least one bond.
#[inline]
fn window(key: Key, t: usize, above: u32, n: usize) -> Option<(u32, u32)> {
    let k = key.count_ones() as usize;
    if k == 0 {
        return None;
    }
    let lo = t + k - 1;
    let pend = (key >> above).count_ones() as usize;
    if lo + pend > n {
        return None;
    }
    Some((lo as u32, (n - pend) as u32))
}

/// Add src (coefficients from order slo) times p^shift, negated if `neg`,
/// into the state with key `key` of `nt`.
#[inline]
fn add_to(nt: &mut Table<Key>, key: Key, src: &[C], slo: u32, shift: u32, neg: bool, t: usize, above: u32, n: usize, m: Mod) {
    let Some((idx, lo, hi)) = nt.get_or_insert(key, |k| prof::time(prof::BOUND, || window(k, t, above, n))) else {
        return;
    };
    let a = slo + shift; // order of src[0] after the shift
    let from = lo.max(a);
    let to = hi.min(a + src.len() as u32 - 1);
    if from > to {
        return;
    }
    let s = &src[(from - a) as usize..=(to - a) as usize];
    nt.accumulate(idx, 0, (from - lo) as usize, s, neg, m);
}

/// One step (parent j of row t) on every state of `srcs`, into a new table.
fn step(srcs: &[Table<Key>], t: usize, j: usize, p: &Params, m: Mod) -> Table<Key> {
    let e = if t % 2 == 0 { 1 } else { 0 };
    let pb = (j + 2 - e) as u32; // bit P: parent j
    let wall = j < e; // the left child would be child -1
    let above = pb + 1; // unprocessed parents after this step
    let n = p.n;
    let total: usize = srcs.iter().map(|s| s.len()).sum();
    let mut nt = Table::new(1, total + total / 4);
    let nc: usize = srcs.iter().map(|s| s.coefs()).sum();
    nt.reserve_coefs(nc + nc / 8);
    let bp: Key = 1 << pb;
    let bl: Key = 1 << (pb - 1);
    let _ts = prof::tsc();
    for s in srcs {
        for (key, slo, shi, w) in s.iter() {
            // Targets whose bit P-1 is clear have exactly one source (the
            // identity for bits 00, the R move for bits 01): append them
            // without hashing. The windows need no clipping: P was not a
            // parent (00), or one fewer parent is pending afterwards (01).
            let unique = key & bl == 0;
            if key & bp == 0 {
                if unique {
                    nt.push_copy(key, slo, shi, w);
                } else {
                    add_to(&mut nt, key, w, slo, 0, false, t, above, n, m);
                }
                continue;
            }
            // R: bond to child j-e+1 only; the key is unchanged
            if unique {
                nt.push_copy(key, slo + 1, shi + 1, w);
            } else {
                add_to(&mut nt, key, w, slo, 1, false, t, above, n, m);
            }
            if wall {
                continue;
            }
            let neg = key & bl != 0; // the pending child already has a bond
            // L: bond to child j-e only
            add_to(&mut nt, (key | bl) & !bp, w, slo, 1, neg, t, above, n, m);
            // LR: both bonds
            add_to(&mut nt, key | bl, w, slo, 2, neg, t, above, n, m);
        }
    }
    prof::add(prof::STEP, prof::tsc() - _ts);
    nt.finish_seq();
    nt
}

/// Bits touched by steps j..jend of row t: P-1 of the first to P of the last.
fn touched(t: usize, j: usize, jend: usize) -> Key {
    let e = if t % 2 == 0 { 1 } else { 0 };
    let lo_bit = j + 1 - e;
    let hi_bit = jend - 1 + 2 - e;
    ((1 << (hi_bit + 1)) - 1) & !((1 << lo_bit) - 1)
}

/// Blocks of steps (t, j, jend) in sweep order.
fn blocks(n: usize, block: usize) -> Vec<(usize, usize, usize)> {
    let mut v = Vec::new();
    for t in 0..n {
        let last = t / 2;
        let mut j = 0;
        while j <= last {
            let jend = (j + block).min(last + 1);
            v.push((t, j, jend));
            j = jend;
        }
    }
    v
}

/// Contributions of the completed row t1 (children c at bit c + 1) to the
/// six series.
fn accumulate(s: &Table<Key>, t1: usize, n: usize, m: Mod) -> [Vec<C>; 6] {
    let mut acc: [Vec<C>; 6] = std::array::from_fn(|_| vec![0 as C; n + 1]);
    for (key, lo, _, w) in s.iter() {
        let k = key.count_ones();
        let neg = k % 2 == 0;
        for (d, &c) in w.iter().enumerate() {
            let o = lo as usize + d;
            acc[0][o] = if neg { m.sub(acc[0][o], c) } else { m.add(acc[0][o], c) };
        }
        if k == 1 {
            let c = key.trailing_zeros() as u64 - 1;
            let x = 2 * c + (t1 as u64 % 2);
            let f = [1, x, x * x, t1 as u64, (t1 * t1) as u64];
            for (mm, &fac) in f.iter().enumerate() {
                let fm = m.from_i64(fac as i64);
                for (d, &c) in w.iter().enumerate() {
                    let o = lo as usize + d;
                    acc[mm + 1][o] = m.add(acc[mm + 1][o], m.mul(c, fm));
                }
            }
        }
    }
    acc
}

fn add_series(ser: &mut [Vec<C>; 6], a: &[Vec<C>; 6], m: Mod) {
    for (s, v) in ser.iter_mut().zip(a.iter()) {
        for (x, &y) in s.iter_mut().zip(v) {
            *x = m.add(*x, y);
        }
    }
}

fn first_table() -> Table<Key> {
    // row 0: parent 0 at bit 0 - 1 + 2 = 1
    let mut first = Table::new(1, 1);
    let (idx, _, _) = first.get_or_insert(1 << 1, |_| Some((0, 0))).unwrap();
    first.coef_mut(idx, 0)[0] = 1;
    first.finish_seq();
    first
}

/// Checkpoint of a spilled run, written after every block whose output is on
/// disk: the series so far, the next block, and the output generation's
/// directory and shard counts. Text, one `key=value` per line.
struct Ckpt {
    next: usize,
    dir: std::path::PathBuf,
    bytes: u64,
    nkeys: Vec<usize>,
    ncoefs: Vec<usize>,
    stats: [u64; 4], // peak states, peak layer bytes, peak shard bytes, disk I/O bytes
    secs: f64,
    ser: [Vec<C>; 6],
}

/// Everything that must match for a checkpoint to be resumed.
fn ckpt_tag(p: &Params, prime: u64) -> String {
    format!(
        "dpwall n={} prime={} block={} shards={} key_bits={} coef_bits={}",
        p.n,
        prime,
        p.block,
        p.shards,
        8 * std::mem::size_of::<Key>(),
        8 * std::mem::size_of::<C>()
    )
}

fn ckpt_path(p: &Params, prime: u64) -> std::path::PathBuf {
    p.spill.as_deref().unwrap().join(format!("dpwall-n{}-p{}.ckpt", p.n, prime))
}

fn join<T: ToString>(v: &[T]) -> String {
    v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

fn ckpt_write(p: &Params, prime: u64, c: &Ckpt) {
    let mut t = format!("{}\nnext={}\ndir={}\nbytes={}\n", ckpt_tag(p, prime), c.next, c.dir.display(), c.bytes);
    t += &format!("nkeys={}\nncoefs={}\nstats={}\nsecs={}\n", join(&c.nkeys), join(&c.ncoefs), join(&c.stats), c.secs);
    for (k, v) in c.ser.iter().enumerate() {
        t += &format!("ser{}={}\n", k, join(v));
    }
    let path = ckpt_path(p, prime);
    let tmp = path.with_extension("tmp");
    {
        let mut f = std::fs::File::create(&tmp).expect("write checkpoint");
        f.write_all(t.as_bytes()).expect("write checkpoint");
        f.sync_all().expect("write checkpoint");
    }
    std::fs::rename(&tmp, &path).expect("replace checkpoint"); // atomic on one volume
}

fn ckpt_read(p: &Params, prime: u64) -> Option<Ckpt> {
    let path = ckpt_path(p, prime);
    let text = std::fs::read_to_string(&path).ok()?;
    let mut lines = text.lines();
    let tag = lines.next()?;
    assert_eq!(tag, ckpt_tag(p, prime), "checkpoint {path:?} is from a different run: delete it to start over");
    let kv: std::collections::HashMap<&str, &str> = lines.filter_map(|l| l.split_once('=')).collect();
    let nums = |k: &str| -> Vec<u64> {
        let v = kv[k];
        if v.is_empty() {
            Vec::new()
        } else {
            v.split(',').map(|x| x.parse().unwrap()).collect()
        }
    };
    let us = |k: &str| nums(k).into_iter().map(|x| x as usize).collect::<Vec<_>>();
    let st = nums("stats");
    Some(Ckpt {
        next: kv["next"].parse().unwrap(),
        dir: kv["dir"].into(),
        bytes: kv["bytes"].parse().unwrap(),
        nkeys: us("nkeys"),
        ncoefs: us("ncoefs"),
        stats: [st[0], st[1], st[2], st[3]],
        secs: kv["secs"].parse().unwrap(),
        ser: std::array::from_fn(|k| nums(&format!("ser{k}")).into_iter().map(|x| x as C).collect()),
    })
}

/// As `run`, but each block reads its input one shard at a time and splits
/// its output straight into the next generation, which goes to disk once a
/// generation exceeds `spill_above` (and stays on disk from then on, so every
/// later block is checkpointed). Memory holds only the shards in flight
/// (about threads x 2 x layer / shards) plus the write buffers.
///
/// Pausing: after each block with output on disk a checkpoint is written, and
/// only then is the input generation deleted. A file `PAUSE` in the spill
/// directory stops the run at the next block boundary (exit code 3); a kill or
/// crash loses at most the block in progress. Rerunning the same command
/// resumes.
fn run_spill(p: &Params, m: Mod, prime: u64) -> [Vec<C>; 6] {
    let n = p.n;
    let dir = p.spill.as_deref().unwrap();
    let bl = blocks(n, p.block);
    let ns = p.shards;
    let t0 = Instant::now();
    let gen_dir = |bi: usize| dir.join(format!("dpwall-n{}-p{}-g{}", n, prime, bi));
    let (mut ser, mut gen, start, stats, secs0) = match ckpt_read(p, prime) {
        Some(c) => {
            let gen = spill::Gen::<Key>::reopen(c.dir.clone(), c.next, &c.nkeys, &c.ncoefs, c.bytes);
            eprintln!("resuming prime {} at block {} of {} ({:.1}s done)", prime, c.next, bl.len(), c.secs);
            (c.ser, gen, c.next, c.stats, c.secs)
        }
        None => {
            let mut ser: [Vec<C>; 6] = std::array::from_fn(|_| vec![0 as C; n + 1]);
            ser[1][0] = 1; // S: the origin
            let gen: spill::Gen<Key> = spill::Gen::new(None, 0, ns);
            let (t, j, jend) = bl[0];
            let tc = touched(t, j, jend);
            spill::write_split(&first_table(), |k| k & !tc, &gen);
            (ser, gen, 0, [0u64; 4], 0.0)
        }
    };
    let [mut peak_states, mut peak_gen, mut peak_shard, mut io] = stats;
    let mut on_disk = gen.on_disk;
    for bi in start..bl.len() {
        let (t, j, jend) = bl[bi];
        let e = if t % 2 == 0 { 1 } else { 0 };
        let row_end = jend == t / 2 + 1;
        let next = bl.get(bi + 1).map(|&(t2, j2, je2)| touched(t2, j2, je2));
        on_disk |= gen.bytes() > p.spill_above;
        let out: spill::Gen<Key> =
            if on_disk && next.is_some() { spill::Gen::at(gen_dir(bi + 1), bi + 1, ns) } else { spill::Gen::new(None, bi + 1, ns) };
        prof::add(prof::N_PHASE, 1);
        let _tp = prof::tsc();
        let res: Vec<(usize, usize, Option<[Vec<C>; 6]>)> = par_map(ns, p.threads, |d| {
            let _tb = prof::tsc();
            let Some(src) = prof::time(prof::TAKE, || gen.take(d, 1)) else {
                return (0, 0, None);
            };
            let mut cur = step(std::slice::from_ref(&src), t, j, p, m);
            let mem = src.bytes() + cur.bytes();
            drop(src);
            for jj in j + 1..jend {
                cur = step(std::slice::from_ref(&cur), t, jj, p, m);
            }
            let acc = if row_end { Some(prof::time(prof::ROWEND, || accumulate(&cur, t + 1, n, m))) } else { None };
            if row_end && e == 1 {
                cur.map_keys(|k| k << 1); // relabel for row t+1
            }
            if let Some(tc) = next {
                prof::time(prof::WRITE, || spill::write_split(&cur, |k| k & !tc, &out));
            }
            prof::add(prof::BUSY, prof::tsc() - _tb);
            (cur.len(), mem.max(cur.bytes()), acc)
        });
        prof::add(prof::M_PAR, prof::tsc() - _tp);
        prof::time(prof::M_CLOSE, || out.close());
        if gen.on_disk {
            io += gen.bytes();
        }
        if out.on_disk {
            io += out.bytes();
        }
        let mut states = 0;
        for (st, mem, acc) in res {
            states += st;
            peak_shard = peak_shard.max(mem as u64);
            if let Some(a) = acc {
                add_series(&mut ser, &a, m);
            }
        }
        peak_states = peak_states.max(states as u64);
        peak_gen = peak_gen.max(out.bytes());
        let secs = secs0 + t0.elapsed().as_secs_f64();
        if p.verbose && row_end {
            eprintln!(
                "row {:3} states {:10} layer MiB {:8.1} {} peak states {:10} peak layer MiB {:8.1} peak shard MiB {:6.1} disk I/O GiB {:7.2} {:8.1}s",
                t + 1,
                states,
                out.bytes() as f64 / 1048576.0,
                if out.on_disk { "disk" } else { "mem " },
                peak_states,
                peak_gen as f64 / 1048576.0,
                peak_shard as f64 / 1048576.0,
                io as f64 / 1073741824.0,
                secs
            );
        }
        let saved = out.on_disk;
        if saved {
            let (nkeys, ncoefs) = out.counts();
            let c = Ckpt {
                next: bi + 1,
                dir: out.directory().unwrap().to_path_buf(),
                bytes: out.bytes(),
                nkeys,
                ncoefs,
                stats: [peak_states, peak_gen, peak_shard, io],
                secs,
                ser: ser.clone(),
            };
            prof::time(prof::M_CKPT, || ckpt_write(p, prime, &c));
        }
        prof::time(prof::M_DROP, || gen = out); // drops the spent generation and deletes its files
        if bi + 1 < bl.len() && dir.join("PAUSE").exists() {
            if saved {
                eprintln!("paused before block {} of {} (row {}), {:.1}s: rerun the same command to resume", bi + 1, bl.len(), t + 1, secs);
            } else {
                eprintln!("paused before the layer reached disk: nothing saved, a rerun starts this prime over");
            }
            prof::report(p.threads);
            std::process::exit(3); // skips Drop, so the checkpointed generation stays
        }
    }
    let _ = std::fs::remove_file(ckpt_path(p, prime));
    if p.verbose {
        eprintln!(
            "done: peak states {} peak layer MiB {:.1} peak shard MiB {:.1} disk I/O GiB {:.2} {:.2}s",
            peak_states,
            peak_gen as f64 / 1048576.0,
            peak_shard as f64 / 1048576.0,
            io as f64 / 1073741824.0,
            secs0 + t0.elapsed().as_secs_f64()
        );
    }
    ser
}

fn run(p: &Params, m: Mod, prime: u64) -> [Vec<C>; 6] {
    if p.spill.is_some() {
        return run_spill(p, m, prime);
    }
    let n = p.n;
    let mut ser: [Vec<C>; 6] = std::array::from_fn(|_| vec![0 as C; n + 1]);
    ser[1][0] = 1; // S: the origin
    let mut layer: Layer = vec![first_table()];
    let t0 = Instant::now();
    let mut peak_states = 0usize;
    let mut peak_bytes = 0usize;
    for t in 0..n {
        let e = if t % 2 == 0 { 1 } else { 0 };
        let last = t / 2; // last parent index of row t
        let mut j = 0;
        while j <= last {
            let jend = (j + p.block).min(last + 1);
            let tc = touched(t, j, jend);
            let pieces = split_layer(std::mem::take(&mut layer), p.shards, p.threads, |k| k & !tc);
            layer = par_map(pieces.len(), p.threads, |si| {
                let mut srcs = std::mem::take(&mut *pieces[si].lock().unwrap());
                let mut cur = step(&srcs, t, j, p, m);
                srcs.clear();
                for jj in j + 1..jend {
                    cur = step(std::slice::from_ref(&cur), t, jj, p, m);
                }
                cur
            });
            let states: usize = layer.iter().map(|s| s.len()).sum();
            let bytes: usize = layer.iter().map(|s| s.bytes()).sum();
            peak_states = peak_states.max(states);
            peak_bytes = peak_bytes.max(bytes);
            j = jend;
        }
        // row t+1 is complete: children c at bit c + 1
        let t1 = t + 1;
        let parts: Vec<[Vec<C>; 6]> = par_map(layer.len(), p.threads, |si| accumulate(&layer[si], t1, n, m));
        for a in parts {
            add_series(&mut ser, &a, m);
        }
        // relabel for row t+1: parent q at bit q - e' + 2, e' = 1 - e
        if e == 1 {
            for s in layer.iter_mut() {
                s.map_keys(|k| k << 1);
            }
        }
        if p.verbose {
            let states: usize = layer.iter().map(|s| s.len()).sum();
            let coefs: usize = layer.iter().map(|s| s.coefs()).sum();
            eprintln!(
                "row {:3} states {:10} coefs {:12} peak states {:10} peak MiB {:8.1} {:8.1}s",
                t1,
                states,
                coefs,
                peak_states,
                peak_bytes as f64 / 1048576.0,
                t0.elapsed().as_secs_f64()
            );
        }
    }
    if p.verbose {
        eprintln!("done: peak states {} peak MiB {:.1} {:.2}s", peak_states, peak_bytes as f64 / 1048576.0, t0.elapsed().as_secs_f64());
    }
    ser
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut n = 20usize;
    let mut nprimes = 1usize;
    let mut first_prime = 0usize;
    let mut out: Option<String> = None;
    let mut verbose = false;
    let mut threads = 1usize;
    let mut block = 6usize;
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
            a => panic!("unknown argument {a}"),
        }
        i += 1;
    }
    // row t uses key bits up to t/2 + 3 (relabelled), and rows run to n - 1
    let key_bits = 8 * std::mem::size_of::<Key>();
    assert!((n - 1) / 2 + 3 < key_bits, "keys are u{key_bits}: n too large (build with --features wide)");
    assert!(n <= 249, "windows are u8");
    // spill mode: many shards keep the shards in flight small
    let shards = shards.unwrap_or(if spill.is_some() { 64 } else if threads <= 1 { 1 } else { threads * 8 });
    if let Some(d) = &spill {
        std::fs::create_dir_all(d).expect("create spill directory");
        // starting the program is the request to run: clear an old pause
        if std::fs::remove_file(d.join("PAUSE")).is_ok() {
            eprintln!("removed {:?}", d.join("PAUSE"));
        }
    }
    let p = Params { n, threads, shards, block, verbose, spill, spill_above: spill_above_mib << 20 };
    let primes = modp::primes(first_prime + nprimes);
    // primes already in the output file are done (so a rerun resumes a queue)
    let done = out.as_ref().and_then(|o| std::fs::read_to_string(o).ok()).unwrap_or_default();
    for &pr in &primes[first_prime..] {
        if done.lines().any(|l| l.contains(&format!("\"n\":{},\"prime\":{},", n, pr))) {
            eprintln!("prime {} already in the output: skipped", pr);
            continue;
        }
        let m = Mod::new(pr);
        prof::start();
        let t0 = Instant::now();
        let ser = run(&p, m, pr);
        let mut line = format!("{{\"model\":\"bond_wall\",\"n\":{},\"prime\":{},\"seconds\":{:.3}", n, pr, t0.elapsed().as_secs_f64());
        for (k, name) in NAMES.iter().enumerate() {
            let v: Vec<String> = ser[k].iter().map(|x| x.to_string()).collect();
            line.push_str(&format!(",\"{}\":[{}]", name, v.join(",")));
        }
        line.push('}');
        eprintln!("prime {} done in {:.2}s", pr, t0.elapsed().as_secs_f64());
        prof::report(threads);
        match &out {
            Some(path) => {
                let mut f = std::fs::OpenOptions::new().create(true).append(true).open(path).unwrap();
                writeln!(f, "{}", line).unwrap();
            }
            None => println!("{}", line),
        }
    }
}
