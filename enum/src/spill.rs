//! Layer generations for dpwall and dpenum, held in memory or spilled to disk.
//!
//! A block of steps reads its input layer one target shard at a time and
//! splits its output by the next block's shard grouping straight into the
//! next generation. On disk a generation is one file per shard, appended in
//! chunks (one per writing shard): a header of two u64 (keys, coefficients),
//! then the keys, lo bytes, hi bytes and coefficients, in native byte order.
//! A generation's files are deleted when it is dropped, not as its shards are
//! read, so the input stays intact until the block that reads it has
//! finished and checkpointed (see `Gen::at` and `Gen::reopen`). So a block is
//! one sequential read and one sequential write of the layer, and memory
//! holds only the shards in flight.

use std::fs::{self, File};
use std::io::{BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

use crate::modp::C;
use crate::prof;
use crate::table::{shard_of, Key, Table};

fn bytes_of<T: Copy>(v: &[T]) -> &[u8] {
    // only plain integers are passed: every byte pattern is valid
    unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) }
}

fn read_into<T: Copy, R: Read>(r: &mut R, v: &mut Vec<T>, n: usize) {
    v.reserve_exact(n);
    unsafe {
        let dst = std::slice::from_raw_parts_mut(v.as_mut_ptr().add(v.len()) as *mut u8, n * std::mem::size_of::<T>());
        r.read_exact(dst).expect("spill read");
        v.set_len(v.len() + n);
    }
}

struct Parts<K> {
    keys: Vec<K>,
    lo: Vec<u8>,
    hi: Vec<u8>,
    coefs: Vec<C>,
}

impl<K> Parts<K> {
    fn new() -> Self {
        Parts { keys: Vec::new(), lo: Vec::new(), hi: Vec::new(), coefs: Vec::new() }
    }
}

/// Pending bytes of a disk shard are handed to the writer threads once they
/// reach this size (so up to nshards x this is buffered, plus the queue).
const FLUSH_BYTES: usize = 1 << 18;
/// Buffers queued for the writer threads before senders block.
const QUEUE: usize = 32;
const WRITERS: usize = 2;

enum Sink<K> {
    Memory(Parts<K>),
    /// serialized chunks not yet handed to a writer thread
    Disk(Vec<u8>),
}

/// Background writer threads for a generation on disk: workers hand them
/// full shard buffers, so file writes overlap computation.
struct Writer {
    tx: std::sync::mpsc::SyncSender<(usize, Vec<u8>)>,
    handles: Vec<std::thread::JoinHandle<()>>,
}

/// One generation of the layer, split into target shards.
pub struct Gen<K: Key> {
    directory: Option<PathBuf>,
    paths: Vec<PathBuf>,
    sinks: Vec<Mutex<Sink<K>>>,
    files: std::sync::Arc<Vec<Mutex<Option<File>>>>,
    writer: Mutex<Option<Writer>>,
    nkeys: Vec<AtomicUsize>,
    ncoefs: Vec<AtomicUsize>,
    bytes: AtomicU64,
    pub on_disk: bool,
}

impl<K: Key> Gen<K> {
    fn blank(nshards: usize, on_disk: bool) -> Self {
        Gen {
            directory: None,
            paths: Vec::new(),
            sinks: (0..nshards).map(|_| Mutex::new(if on_disk { Sink::Disk(Vec::new()) } else { Sink::Memory(Parts::new()) })).collect(),
            files: std::sync::Arc::new((0..nshards).map(|_| Mutex::new(None)).collect()),
            writer: Mutex::new(None),
            nkeys: (0..nshards).map(|_| AtomicUsize::new(0)).collect(),
            ncoefs: (0..nshards).map(|_| AtomicUsize::new(0)).collect(),
            bytes: AtomicU64::new(0),
            on_disk,
        }
    }
    /// `dir`: spill to files there; None: keep the generation in memory.
    pub fn new(dir: Option<&Path>, id: usize, nshards: usize) -> Self {
        let on_disk = dir.is_some();
        // Atomically reserve a private namespace, even when processes share
        // the spill root or a previous process left files behind.
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let directory = dir.map(|root| loop {
            let serial = NEXT.fetch_add(1, Ordering::Relaxed);
            let path = root.join(format!("dpwall-{}-{serial}-g{id}", std::process::id()));
            match fs::create_dir(&path) {
                Ok(()) => break path,
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => panic!("create spill namespace {path:?}: {e}"),
            }
        });
        let mut g = Gen::blank(nshards, on_disk);
        if let Some(d) = &directory {
            g.paths = (0..nshards).map(|s| d.join(format!("g{id}_{s}.bin"))).collect();
        }
        g.directory = directory;
        g
    }
    /// A generation on disk at a fixed directory (for checkpointed runs);
    /// anything already there, such as a block interrupted by a crash, is
    /// removed first.
    pub fn at(dir: PathBuf, id: usize, nshards: usize) -> Self {
        if dir.exists() {
            fs::remove_dir_all(&dir).unwrap_or_else(|e| panic!("clear spill directory {dir:?}: {e}"));
        }
        fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("create spill directory {dir:?}: {e}"));
        let mut g = Gen::blank(nshards, true);
        g.paths = (0..nshards).map(|s| dir.join(format!("g{id}_{s}.bin"))).collect();
        g.directory = Some(dir);
        g
    }
    /// Reopen a closed generation written by `at` (resuming a checkpoint).
    pub fn reopen(dir: PathBuf, id: usize, nkeys: &[usize], ncoefs: &[usize], bytes: u64) -> Self {
        let mut g = Gen::blank(nkeys.len(), true);
        g.paths = (0..nkeys.len()).map(|s| dir.join(format!("g{id}_{s}.bin"))).collect();
        for (d, (&k, &c)) in nkeys.iter().zip(ncoefs).enumerate() {
            g.nkeys[d].store(k, Ordering::Relaxed);
            g.ncoefs[d].store(c, Ordering::Relaxed);
            if k > 0 {
                assert!(g.paths[d].exists(), "checkpoint shard missing: {:?}", g.paths[d]);
            }
        }
        g.bytes.store(bytes, Ordering::Relaxed);
        g.directory = Some(dir);
        g
    }
    pub fn directory(&self) -> Option<&Path> {
        self.directory.as_deref()
    }
    /// Keys and coefficients per shard.
    pub fn counts(&self) -> (Vec<usize>, Vec<usize>) {
        (self.nkeys.iter().map(|a| a.load(Ordering::Relaxed)).collect(), self.ncoefs.iter().map(|a| a.load(Ordering::Relaxed)).collect())
    }
    pub fn nshards(&self) -> usize {
        self.sinks.len()
    }
    /// Size of the generation in serialized form.
    pub fn bytes(&self) -> u64 {
        self.bytes.load(Ordering::Relaxed)
    }
    fn account(&self, d: usize, nk: usize, nc: usize) {
        self.nkeys[d].fetch_add(nk, Ordering::Relaxed);
        self.ncoefs[d].fetch_add(nc, Ordering::Relaxed);
        self.bytes.fetch_add(chunk_bytes::<K>(nk, nc) as u64, Ordering::Relaxed);
    }
    /// Hand a full shard buffer to the writer threads (starting them).
    fn send(&self, d: usize, buf: Vec<u8>) {
        let mut guard = self.writer.lock().unwrap();
        let w = guard.get_or_insert_with(|| {
            let (tx, rx) = std::sync::mpsc::sync_channel::<(usize, Vec<u8>)>(QUEUE);
            let rx = std::sync::Arc::new(Mutex::new(rx));
            let handles = (0..WRITERS)
                .map(|_| {
                    let rx = rx.clone();
                    let files = self.files.clone();
                    let paths = self.paths.clone();
                    std::thread::spawn(move || loop {
                        let msg = rx.lock().unwrap().recv();
                        let Ok((d, buf)) = msg else { return };
                        let mut f = files[d].lock().unwrap();
                        if f.is_none() {
                            *f = Some(File::create(&paths[d]).unwrap_or_else(|e| panic!("spill create {:?}: {e}", paths[d])));
                        }
                        f.as_mut().unwrap().write_all(&buf).expect("spill write (disk full?)");
                    })
                })
                .collect();
            Writer { tx, handles }
        });
        let tx = w.tx.clone();
        drop(guard); // do not hold the lock while the queue applies backpressure
        tx.send((d, buf)).expect("spill writer thread died");
    }
    /// Append one chunk to shard d.
    fn write(&self, d: usize, keys: &[K], lo: &[u8], hi: &[u8], coefs: &[C]) {
        let mut g = self.sinks[d].lock().unwrap();
        match &mut *g {
            Sink::Memory(p) => {
                p.keys.extend_from_slice(keys);
                p.lo.extend_from_slice(lo);
                p.hi.extend_from_slice(hi);
                p.coefs.extend_from_slice(coefs);
            }
            Sink::Disk(buf) => {
                let hdr = [keys.len() as u64, coefs.len() as u64];
                for b in [bytes_of(&hdr), bytes_of(keys), lo, hi, bytes_of(coefs)] {
                    buf.extend_from_slice(b);
                }
                if buf.len() >= FLUSH_BYTES {
                    let full = std::mem::take(buf);
                    drop(g);
                    self.send(d, full);
                }
            }
        }
        self.account(d, keys.len(), coefs.len());
    }
    /// Flush and close every file (call once all writers are done).
    pub fn close(&self) {
        for d in 0..self.sinks.len() {
            let mut g = self.sinks[d].lock().unwrap();
            if let Sink::Disk(buf) = &mut *g {
                if !buf.is_empty() {
                    let full = std::mem::take(buf);
                    drop(g);
                    self.send(d, full);
                }
            }
        }
        if let Some(w) = self.writer.lock().unwrap().take() {
            drop(w.tx);
            for h in w.handles {
                h.join().expect("spill writer thread panicked");
            }
        }
        for f in self.files.iter() {
            f.lock().unwrap().take();
        }
    }
    /// Take shard d as a sequential table with `nm` moments, or None if empty.
    /// Its file stays until the generation is dropped.
    pub fn take(&self, d: usize, nm: usize) -> Option<Table<K>> {
        let nk = self.nkeys[d].load(Ordering::Relaxed);
        if nk == 0 {
            return None;
        }
        let nc = self.ncoefs[d].load(Ordering::Relaxed);
        let mut g = self.sinks[d].lock().unwrap();
        match &mut *g {
            Sink::Memory(p) => {
                let p = std::mem::replace(p, Parts::new());
                Some(Table::from_seq_parts(nm, p.keys, p.lo, p.hi, p.coefs))
            }
            Sink::Disk(buf) => {
                assert!(buf.is_empty(), "spill shard {d} read before close");
                let _t = prof::tsc();
                let path = &self.paths[d];
                let file = File::open(path).unwrap_or_else(|e| panic!("spill open {path:?}: {e}"));
                let mut r = BufReader::with_capacity(1 << 20, file);
                let mut keys: Vec<K> = Vec::with_capacity(nk);
                let mut lo: Vec<u8> = Vec::with_capacity(nk);
                let mut hi: Vec<u8> = Vec::with_capacity(nk);
                let mut coefs: Vec<C> = Vec::with_capacity(nc);
                let mut hdr: Vec<u64> = Vec::with_capacity(2);
                while keys.len() < nk {
                    hdr.clear();
                    read_into(&mut r, &mut hdr, 2);
                    let (a, b) = (hdr[0] as usize, hdr[1] as usize);
                    read_into(&mut r, &mut keys, a);
                    read_into(&mut r, &mut lo, a);
                    read_into(&mut r, &mut hi, a);
                    read_into(&mut r, &mut coefs, b);
                }
                assert_eq!(coefs.len(), nc, "spill shard {d}: coefficient count");
                prof::add(prof::T_READ, prof::tsc() - _t);
                Some(Table::from_seq_parts(nm, keys, lo, hi, coefs))
            }
        }
    }
}

fn chunk_bytes<K: Key>(nk: usize, nc: usize) -> usize {
    16 + nk * (std::mem::size_of::<K>() + 2) + nc * std::mem::size_of::<C>()
}

impl<K: Key> Drop for Gen<K> {
    fn drop(&mut self) {
        // Stop the writers and close the files before deleting them (Windows).
        if let Some(w) = self.writer.get_mut().unwrap().take() {
            drop(w.tx);
            for h in w.handles {
                let _ = h.join();
            }
        }
        for f in self.files.iter() {
            f.lock().unwrap().take();
        }
        for p in &self.paths {
            let _ = fs::remove_file(p);
        }
        if let Some(d) = &self.directory {
            let _ = fs::remove_dir(d);
        }
    }
}

/// Reads the shards of a generation ahead of the workers, in shard order
/// (the order `par_map` hands them out), at most `depth` shards ahead of
/// those taken. Run `reader` on a few threads and `get` from the workers.
pub struct Prefetch<'a, K: Key> {
    gen: &'a Gen<K>,
    nm: usize,
    depth: usize,
    state: Mutex<PrefetchState<K>>,
    cv: std::sync::Condvar,
}

struct PrefetchState<K: Key> {
    next: usize,
    taken: usize,
    ready: Vec<Option<Option<Table<K>>>>,
}

impl<'a, K: Key> Prefetch<'a, K> {
    pub fn new(gen: &'a Gen<K>, nm: usize, depth: usize) -> Self {
        let ns = gen.nshards();
        Prefetch {
            gen,
            nm,
            depth: depth.max(1),
            state: Mutex::new(PrefetchState { next: 0, taken: 0, ready: (0..ns).map(|_| None).collect() }),
            cv: std::sync::Condvar::new(),
        }
    }
    /// Reader thread body: read shards in order until all are read.
    pub fn reader(&self) {
        loop {
            let r = {
                let mut st = self.state.lock().unwrap();
                while st.next < st.ready.len() && st.next >= st.taken + self.depth {
                    st = self.cv.wait(st).unwrap();
                }
                if st.next >= st.ready.len() {
                    return;
                }
                st.next += 1;
                st.next - 1
            };
            let t = self.gen.take(r, self.nm);
            self.state.lock().unwrap().ready[r] = Some(t);
            self.cv.notify_all();
        }
    }
    /// Shard d (None if empty), waiting for its reader if needed.
    pub fn get(&self, d: usize) -> Option<Table<K>> {
        let mut st = self.state.lock().unwrap();
        loop {
            if let Some(t) = st.ready[d].take() {
                st.taken += 1;
                drop(st);
                self.cv.notify_all();
                return t;
            }
            st = self.cv.wait(st).unwrap();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generations_sharing_a_root_do_not_overwrite_or_delete_each_other() {
        let root = std::env::temp_dir();
        let a = Gen::<u64>::new(Some(&root), 1, 1);
        let b = Gen::<u64>::new(Some(&root), 1, 1);
        for (out, value) in [(&a, 11 as C), (&b, 22 as C)] {
            out.write(0, &[2], &[0], &[0], &[value]);
            out.close();
        }
        assert_eq!(a.take(0, 1).unwrap().arena, vec![11 as C]);
        let a_dir = a.directory.clone().unwrap();
        drop(a);
        assert!(!a_dir.exists());
        assert_eq!(b.take(0, 1).unwrap().arena, vec![22 as C]);
        let b_dir = b.directory.clone().unwrap();
        drop(b);
        assert!(!b_dir.exists());
    }
}

/// Split a finished (or append-only) table by `group` into the shards of
/// `out`: one chunk per target shard, gathered through a counting sort. For a
/// generation on disk, each chunk is serialized straight into the pending
/// buffer of its shard.
pub fn write_split<K: Key, F: Fn(K) -> K>(t: &Table<K>, group: F, out: &Gen<K>) {
    let _t0 = prof::tsc();
    let ns = out.nshards();
    let n = t.len();
    let mut dest: Vec<u32> = Vec::with_capacity(n);
    let mut offs: Vec<u32> = Vec::with_capacity(n);
    let mut lens: Vec<u32> = Vec::with_capacity(n);
    let mut cnt = vec![0usize; ns];
    let mut ncs = vec![0usize; ns];
    let mut o = 0usize;
    for (key, _, _, b) in t.iter() {
        let d = shard_of(group(key), ns);
        dest.push(d as u32);
        offs.push(o as u32);
        lens.push(b.len() as u32);
        o += b.len();
        cnt[d] += 1;
        ncs[d] += b.len();
    }
    assert!(o <= u32::MAX as usize, "shard arena exceeds 2^32 coefficients: use more shards");
    let mut start = vec![0usize; ns + 1];
    for d in 0..ns {
        start[d + 1] = start[d] + cnt[d];
    }
    let mut fill = start.clone();
    let mut perm = vec![0u32; n];
    for (i, &d) in dest.iter().enumerate() {
        perm[fill[d as usize]] = i as u32;
        fill[d as usize] += 1;
    }
    drop(dest);
    prof::add(prof::W_DEST, prof::tsc() - _t0);
    if out.on_disk {
        let ks = std::mem::size_of::<K>();
        let cs = std::mem::size_of::<C>();
        for d in 0..ns {
            if cnt[d] == 0 {
                continue;
            }
            let _tg = prof::tsc();
            let idx = &perm[start[d]..start[d + 1]];
            let (nk, nc) = (cnt[d], ncs[d]);
            let size = chunk_bytes::<K>(nk, nc);
            let mut g = out.sinks[d].lock().unwrap();
            let Sink::Disk(buf) = &mut *g else { unreachable!() };
            let base = buf.len();
            buf.resize(base + size, 0);
            let (hdr, rest) = buf[base..].split_at_mut(16);
            hdr[..8].copy_from_slice(&(nk as u64).to_ne_bytes());
            hdr[8..].copy_from_slice(&(nc as u64).to_ne_bytes());
            let (kb, rest) = rest.split_at_mut(nk * ks);
            let (lb, rest) = rest.split_at_mut(nk);
            let (hb, cb) = rest.split_at_mut(nk);
            let mut co = 0usize;
            for (k, &i) in idx.iter().enumerate() {
                let i = i as usize;
                let (l, h) = t.window(i);
                kb[k * ks..(k + 1) * ks].copy_from_slice(bytes_of(&[t.key(i)]));
                lb[k] = l as u8;
                hb[k] = h as u8;
                let s = offs[i] as usize;
                let len = lens[i] as usize;
                cb[co * cs..(co + len) * cs].copy_from_slice(bytes_of(&t.arena[s..s + len]));
                co += len;
            }
            let full = if buf.len() >= FLUSH_BYTES { Some(std::mem::take(buf)) } else { None };
            drop(g);
            out.account(d, nk, nc);
            prof::add(prof::W_GATHER, prof::tsc() - _tg);
            if let Some(full) = full {
                let _tw = prof::tsc();
                out.send(d, full);
                prof::add(prof::W_OUT, prof::tsc() - _tw);
            }
        }
        return;
    }
    let mut p: Parts<K> = Parts::new();
    for d in 0..ns {
        if cnt[d] == 0 {
            continue;
        }
        p.keys.clear();
        p.lo.clear();
        p.hi.clear();
        p.coefs.clear();
        for &i in &perm[start[d]..start[d + 1]] {
            let i = i as usize;
            let (l, h) = t.window(i);
            p.keys.push(t.key(i));
            p.lo.push(l as u8);
            p.hi.push(h as u8);
            let s = offs[i] as usize;
            p.coefs.extend_from_slice(&t.arena[s..s + lens[i] as usize]);
        }
        out.write(d, &p.keys, &p.lo, &p.hi, &p.coefs);
    }
}

/// Write a small text checkpoint atomically: to a temporary file, synced, then
/// renamed over the old one (atomic on one volume).
pub fn save_text(path: &Path, text: &str) {
    let tmp = path.with_extension("tmp");
    {
        let mut f = File::create(&tmp).unwrap_or_else(|e| panic!("write checkpoint {tmp:?}: {e}"));
        f.write_all(text.as_bytes()).expect("write checkpoint");
        f.sync_all().expect("write checkpoint");
    }
    fs::rename(&tmp, path).unwrap_or_else(|e| panic!("replace checkpoint {path:?}: {e}"));
}

/// Read a checkpoint written by `save_text` as `key=value` lines after a tag
/// line, which must equal `tag` (else the checkpoint is from another run).
pub fn load_text(path: &Path, tag: &str) -> Option<std::collections::HashMap<String, String>> {
    let text = fs::read_to_string(path).ok()?;
    let mut lines = text.lines();
    let first = lines.next()?;
    assert_eq!(first, tag, "checkpoint {path:?} is from a different run: delete it to start over");
    Some(lines.filter_map(|l| l.split_once('=')).map(|(k, v)| (k.to_string(), v.to_string())).collect())
}

/// Comma-separated numbers (empty for an empty list).
pub fn join<T: ToString>(v: &[T]) -> String {
    v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

pub fn split_nums(v: &str) -> Vec<u64> {
    if v.is_empty() {
        Vec::new()
    } else {
        v.split(',').map(|x| x.parse().expect("checkpoint number")).collect()
    }
}
