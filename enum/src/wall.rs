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
}

/// Keys are u64: row t uses bits up to t/2 + 3, so n <= 120 (checked in main).
type Key = u64;
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
    let Some((idx, lo, hi)) = nt.get_or_insert(key, |k| window(k, t, above, n)) else {
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
    let bp: Key = 1 << pb;
    let bl: Key = 1 << (pb - 1);
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
    nt.finish();
    nt.drop_offsets();
    nt
}

fn run(p: &Params, m: Mod) -> [Vec<C>; 6] {
    let n = p.n;
    let mut ser: [Vec<C>; 6] = std::array::from_fn(|_| vec![0 as C; n + 1]);
    ser[1][0] = 1; // S: the origin
    // row 0: parent 0 at bit 0 - 1 + 2 = 1
    let mut first = Table::new(1, 1);
    let (idx, _, _) = first.get_or_insert(1 << 1, |_| Some((0, 0))).unwrap();
    first.coef_mut(idx, 0)[0] = 1;
    first.finish();
    first.drop_offsets();
    let mut layer: Layer = vec![first];
    let t0 = Instant::now();
    let mut peak_states = 0usize;
    let mut peak_bytes = 0usize;
    for t in 0..n {
        let e = if t % 2 == 0 { 1 } else { 0 };
        let last = t / 2; // last parent index of row t
        let mut j = 0;
        while j <= last {
            let jend = (j + p.block).min(last + 1);
            // bits touched by steps j..jend: P-1 of the first to P of the last
            let lo_bit = j + 1 - e;
            let hi_bit = jend - 1 + 2 - e;
            let touched: Key = ((1 << (hi_bit + 1)) - 1) & !((1 << lo_bit) - 1);
            let pieces = split_layer(std::mem::take(&mut layer), p.shards, p.threads, |k| k & !touched);
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
        let parts: Vec<[Vec<C>; 6]> = par_map(layer.len(), p.threads, |si| {
            let s = &layer[si];
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
        });
        for a in parts {
            for (s, v) in ser.iter_mut().zip(a.iter()) {
                for (x, &y) in s.iter_mut().zip(v) {
                    *x = m.add(*x, y);
                }
            }
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
    assert!((n - 1) / 2 + 3 < 64, "keys are u64: n <= 121");
    let shards = if threads <= 1 { 1 } else { threads * 8 };
    let p = Params { n, threads, shards, block, verbose };
    let primes = modp::primes(first_prime + nprimes);
    for &pr in &primes[first_prime..] {
        let m = Mod::new(pr);
        let t0 = Instant::now();
        let ser = run(&p, m);
        let mut line = format!("{{\"model\":\"bond_wall\",\"n\":{},\"prime\":{},\"seconds\":{:.3}", n, pr, t0.elapsed().as_secs_f64());
        for (k, name) in NAMES.iter().enumerate() {
            let v: Vec<String> = ser[k].iter().map(|x| x.to_string()).collect();
            line.push_str(&format!(",\"{}\":[{}]", name, v.join(",")));
        }
        line.push('}');
        eprintln!("prime {} done in {:.2}s", pr, t0.elapsed().as_secs_f64());
        match &out {
            Some(path) => {
                let mut f = std::fs::OpenOptions::new().create(true).append(true).open(path).unwrap();
                writeln!(f, "{}", line).unwrap();
            }
            None => println!("{}", line),
        }
    }
}
