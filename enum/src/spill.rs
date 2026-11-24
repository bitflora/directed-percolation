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
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

use crate::modp::C;
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

enum Sink<K> {
    Memory(Parts<K>),
    Disk(Option<BufWriter<File>>),
}

/// One generation of the layer, split into target shards.
pub struct Gen<K: Key> {
    directory: Option<PathBuf>,
    paths: Vec<PathBuf>,
    sinks: Vec<Mutex<Sink<K>>>,
    nkeys: Vec<AtomicUsize>,
    ncoefs: Vec<AtomicUsize>,
    bytes: AtomicU64,
    pub on_disk: bool,
}

impl<K: Key> Gen<K> {
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
        Gen {
            paths: match &directory {
                Some(d) => (0..nshards).map(|s| d.join(format!("g{id}_{s}.bin"))).collect(),
                None => Vec::new(),
            },
            directory,
            sinks: (0..nshards).map(|_| Mutex::new(if on_disk { Sink::Disk(None) } else { Sink::Memory(Parts::new()) })).collect(),
            nkeys: (0..nshards).map(|_| AtomicUsize::new(0)).collect(),
            ncoefs: (0..nshards).map(|_| AtomicUsize::new(0)).collect(),
            bytes: AtomicU64::new(0),
            on_disk,
        }
    }
    /// A generation on disk at a fixed directory (for checkpointed runs);
    /// anything already there, such as a block interrupted by a crash, is
    /// removed first.
    pub fn at(dir: PathBuf, id: usize, nshards: usize) -> Self {
        if dir.exists() {
            fs::remove_dir_all(&dir).unwrap_or_else(|e| panic!("clear spill directory {dir:?}: {e}"));
        }
        fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("create spill directory {dir:?}: {e}"));
        let mut g = Gen::new(None, id, nshards);
        g.paths = (0..nshards).map(|s| dir.join(format!("g{id}_{s}.bin"))).collect();
        g.sinks = (0..nshards).map(|_| Mutex::new(Sink::Disk(None))).collect();
        g.directory = Some(dir);
        g.on_disk = true;
        g
    }
    /// Reopen a closed generation written by `at` (resuming a checkpoint).
    pub fn reopen(dir: PathBuf, id: usize, nkeys: &[usize], ncoefs: &[usize], bytes: u64) -> Self {
        let mut g = Gen::new(None, id, nkeys.len());
        g.paths = (0..nkeys.len()).map(|s| dir.join(format!("g{id}_{s}.bin"))).collect();
        g.sinks = (0..nkeys.len()).map(|_| Mutex::new(Sink::Disk(None))).collect();
        for (d, (&k, &c)) in nkeys.iter().zip(ncoefs).enumerate() {
            g.nkeys[d].store(k, Ordering::Relaxed);
            g.ncoefs[d].store(c, Ordering::Relaxed);
            if k > 0 {
                assert!(g.paths[d].exists(), "checkpoint shard missing: {:?}", g.paths[d]);
            }
        }
        g.bytes.store(bytes, Ordering::Relaxed);
        g.directory = Some(dir);
        g.on_disk = true;
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
    fn write(&self, d: usize, keys: &[K], lo: &[u8], hi: &[u8], coefs: &[C]) {
        let mut g = self.sinks[d].lock().unwrap();
        match &mut *g {
            Sink::Memory(p) => {
                p.keys.extend_from_slice(keys);
                p.lo.extend_from_slice(lo);
                p.hi.extend_from_slice(hi);
                p.coefs.extend_from_slice(coefs);
            }
            Sink::Disk(f) => {
                if f.is_none() {
                    let file = File::create(&self.paths[d]).unwrap_or_else(|e| panic!("spill create {:?}: {e}", self.paths[d]));
                    *f = Some(BufWriter::with_capacity(1 << 16, file));
                }
                let w = f.as_mut().unwrap();
                let hdr = [keys.len() as u64, coefs.len() as u64];
                for b in [bytes_of(&hdr), bytes_of(keys), lo, hi, bytes_of(coefs)] {
                    w.write_all(b).expect("spill write (disk full?)");
                }
            }
        }
        self.nkeys[d].fetch_add(keys.len(), Ordering::Relaxed);
        self.ncoefs[d].fetch_add(coefs.len(), Ordering::Relaxed);
        let b = 16 + std::mem::size_of_val(keys) + lo.len() + hi.len() + std::mem::size_of_val(coefs);
        self.bytes.fetch_add(b as u64, Ordering::Relaxed);
    }
    /// Flush and close every file (call once all writers are done).
    pub fn close(&self) {
        for s in &self.sinks {
            if let Sink::Disk(f) = &mut *s.lock().unwrap() {
                if let Some(mut w) = f.take() {
                    w.flush().expect("spill flush (disk full?)");
                }
            }
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
            Sink::Disk(_) => {
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
                Some(Table::from_seq_parts(nm, keys, lo, hi, coefs))
            }
        }
    }
}

impl<K: Key> Drop for Gen<K> {
    fn drop(&mut self) {
        // Close any remaining writers before deleting files on Windows.
        self.sinks.clear();
        for p in &self.paths {
            let _ = fs::remove_file(p);
        }
        if let Some(d) = &self.directory {
            let _ = fs::remove_dir(d);
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
/// `out`: one chunk per target shard, gathered through a counting sort.
pub fn write_split<K: Key, F: Fn(K) -> K>(t: &Table<K>, group: F, out: &Gen<K>) {
    let ns = out.nshards();
    let n = t.len();
    let mut dest: Vec<u32> = Vec::with_capacity(n);
    let mut offs: Vec<u32> = Vec::with_capacity(n);
    let mut lens: Vec<u32> = Vec::with_capacity(n);
    let mut cnt = vec![0usize; ns];
    let mut o = 0usize;
    for (key, _, _, b) in t.iter() {
        let d = shard_of(group(key), ns);
        dest.push(d as u32);
        offs.push(o as u32);
        lens.push(b.len() as u32);
        o += b.len();
        cnt[d] += 1;
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
