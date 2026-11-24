//! State table: open-addressing hash from u128 keys to per-state coefficient
//! windows stored contiguously in one arena.
//!
//! Keys that only one source can produce are appended without hashing
//! (`push_unique`). A state whose bound says it cannot contribute is remembered
//! as dead (hi < lo) so its bound is computed only once per key per step.

use crate::modp::{Mod, C};

/// A state key: u128 for the bulk engine, u64 where the boundary fits.
pub trait Key: Copy + Eq {
    const BYTES: usize;
    fn hash64(self) -> u64;
}

impl Key for u128 {
    const BYTES: usize = 16;
    #[inline(always)]
    fn hash64(self) -> u64 {
        let x = (self as u64) ^ ((self >> 64) as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        let x = x.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        x ^ (x >> 29)
    }
}

impl Key for u64 {
    const BYTES: usize = 8;
    #[inline(always)]
    fn hash64(self) -> u64 {
        let x = self.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        x ^ (x >> 29)
    }
}

pub struct Table<K: Key = u128> {
    nm: usize,
    keys: Vec<K>,
    lo: Vec<u8>,  // window ends: n <= 250, so a byte suffices
    hi: Vec<u8>,
    off: Vec<u32>, // per-table (per-shard) arena offset in coefficients; empty when `seq`
    seq: bool,     // blocks are read in order only (`iter`), so no offsets are kept
    pub arena: Vec<C>,
    slots: Vec<u32>, // index+1, 0 = empty
    mask: usize,
    hashed: usize,
    is_hashed: Vec<bool>,
}

#[inline(always)]
fn hash<K: Key>(k: K) -> u64 {
    k.hash64()
}

#[inline(always)]
pub fn shard_of<K: Key>(base: K, nshards: usize) -> usize {
    let h = hash(base).wrapping_mul(0x94D0_49BB_1331_11EB);
    ((h >> 32) as usize * nshards) >> 32
}

impl<K: Key> Table<K> {
    pub fn new(nm: usize, expect: usize) -> Self {
        let cap = (expect * 2).next_power_of_two().max(16);
        Table {
            nm,
            keys: Vec::with_capacity(expect),
            lo: Vec::with_capacity(expect),
            hi: Vec::with_capacity(expect),
            off: Vec::with_capacity(expect),
            arena: Vec::new(),
            slots: vec![0; cap],
            mask: cap - 1,
            hashed: 0,
            is_hashed: Vec::with_capacity(expect),
            seq: false,
        }
    }
    /// An append-only table (filled by `push_copy` only) with exact capacity.
    pub fn with_exact(nm: usize, keys: usize, coefs: usize) -> Self {
        Table {
            nm,
            keys: Vec::with_capacity(keys),
            lo: Vec::with_capacity(keys),
            hi: Vec::with_capacity(keys),
            off: Vec::with_capacity(keys),
            arena: Vec::with_capacity(coefs),
            slots: Vec::new(),
            mask: 0,
            hashed: 0,
            is_hashed: Vec::with_capacity(keys),
            seq: false,
        }
    }
    /// As `with_exact`, for a table that is only ever read in order (`iter`):
    /// no per-state offsets are stored.
    #[allow(dead_code)] // used by dpwall only
    pub fn with_exact_seq(nm: usize, keys: usize, coefs: usize) -> Self {
        let mut t = Self::with_exact(nm, keys, coefs);
        t.off = Vec::new();
        t.seq = true;
        t
    }
    /// Drop the per-state offsets of a finished table; afterwards its blocks
    /// can only be read in order, with `iter`.
    /// A sequential table from its parts (read back from a spill file).
    #[allow(dead_code)] // used by dpwall only
    pub fn from_seq_parts(nm: usize, keys: Vec<K>, lo: Vec<u8>, hi: Vec<u8>, arena: Vec<C>) -> Self {
        Table { nm, keys, lo, hi, off: Vec::new(), arena, slots: Vec::new(), mask: 0, hashed: 0, is_hashed: Vec::new(), seq: true }
    }
    #[allow(dead_code)] // used by dpwall only
    pub fn drop_offsets(&mut self) {
        assert!(self.slots.is_empty(), "call after finish");
        self.off = Vec::new();
        self.seq = true;
    }
    /// (key, lo, hi, block) for every state in order. Blocks are contiguous
    /// in index order in finished and append-only tables.
    #[allow(dead_code)] // used by dpwall only
    pub fn iter(&self) -> impl Iterator<Item = (K, u32, u32, &[C])> + '_ {
        assert!(self.slots.is_empty(), "iter needs a finished or append-only table");
        let mut o = 0usize;
        (0..self.keys.len()).map(move |i| {
            let (lo, hi) = (self.lo[i] as u32, self.hi[i] as u32);
            let len = (hi - lo + 1) as usize * self.nm;
            let b = &self.arena[o..o + len];
            o += len;
            (self.keys[i], lo, hi, b)
        })
    }
    /// Coefficients stored (all moments).
    pub fn coefs(&self) -> usize {
        self.arena.len()
    }
    /// Reserve the arena up front (estimate from the sources), so it does not
    /// grow by repeated reallocation, each briefly holding two copies.
    pub fn reserve_coefs(&mut self, n: usize) {
        self.arena.reserve_exact(n.saturating_sub(self.arena.len()));
    }
    pub fn len(&self) -> usize {
        self.keys.len()
    }
    pub fn bytes(&self) -> usize {
        self.arena.capacity() * std::mem::size_of::<C>() + self.keys.capacity() * (K::BYTES + 2) + self.off.capacity() * 4 + self.is_hashed.capacity() + self.slots.capacity() * 4
    }
    #[inline(always)]
    pub fn key(&self, i: usize) -> K {
        self.keys[i]
    }
    #[inline(always)]
    pub fn window(&self, i: usize) -> (u32, u32) {
        (self.lo[i] as u32, self.hi[i] as u32)
    }
    fn grow(&mut self) {
        let cap = self.slots.len() * 2;
        self.slots = vec![0; cap];
        self.mask = cap - 1;
        for i in 0..self.keys.len() {
            if !self.is_hashed[i] {
                continue; // appended by push_unique
            }
            let k = self.keys[i];
            let mut h = hash(k) as usize & self.mask;
            while self.slots[h] != 0 {
                h = (h + 1) & self.mask;
            }
            self.slots[h] = i as u32 + 1;
        }
    }
    fn push(&mut self, key: K, lo: u32, hi: u32) -> usize {
        let idx = self.keys.len();
        self.is_hashed.push(false);
        self.keys.push(key);
        self.lo.push(lo as u8);
        self.hi.push(hi as u8);
        assert!(self.arena.len() < u32::MAX as usize, "shard arena exceeds 2^32 coefficients: use more shards");
        self.off.push(self.arena.len() as u32);
        if hi >= lo {
            let len = (hi - lo + 1) as usize * self.nm;
            let need = self.arena.len() + len;
            if need > self.arena.capacity() {
                // grow by 25% instead of doubling to limit slack at the peak
                let extra = (self.arena.capacity() / 4).max(len).max(1024);
                self.arena.reserve_exact(extra);
            }
            self.arena.resize(need, 0);
        }
        idx
    }
    /// Append a key that no other source can produce in this step (no hashing).
    #[inline]
    pub fn push_unique(&mut self, key: K, lo: u32, hi: u32) -> usize {
        self.push(key, lo, hi)
    }
    /// Append a unique key whose coefficient block (all moments) is copied
    /// verbatim, without zero-filling the arena first.
    #[inline]
    pub fn push_copy(&mut self, key: K, lo: u32, hi: u32, block: &[C]) {
        self.is_hashed.push(false);
        self.keys.push(key);
        self.lo.push(lo as u8);
        self.hi.push(hi as u8);
        if !self.seq {
            assert!(self.arena.len() < u32::MAX as usize, "shard arena exceeds 2^32 coefficients: use more shards");
            self.off.push(self.arena.len() as u32);
        }
        let need = self.arena.len() + block.len();
        if need > self.arena.capacity() {
            let extra = (self.arena.capacity() / 4).max(block.len()).max(1024);
            self.arena.reserve_exact(extra);
        }
        self.arena.extend_from_slice(block);
    }
    /// The whole coefficient block of state i (all moments, contiguous).
    #[inline(always)]
    pub fn block(&self, i: usize) -> &[C] {
        let len = (self.hi[i] - self.lo[i] + 1) as usize * self.nm;
        let o = self.off[i] as usize;
        &self.arena[o..o + len]
    }
    /// Index and window of an existing hashed key, or compute its bound and
    /// insert it. Returns None if the state is dead. Probes the table once.
    #[inline]
    pub fn get_or_insert<F: FnOnce(K) -> Option<(u32, u32)>>(
        &mut self,
        key: K,
        f: F,
    ) -> Option<(usize, u32, u32)> {
        if (self.hashed + 1) * 2 > self.slots.len() {
            self.grow();
        }
        let mut h = hash(key) as usize & self.mask;
        loop {
            let s = self.slots[h];
            if s == 0 {
                break;
            }
            let i = s as usize - 1;
            if self.keys[i] == key {
                let (lo, hi) = self.window(i);
                return if hi >= lo { Some((i, lo, hi)) } else { None };
            }
            h = (h + 1) & self.mask;
        }
        let (lo, hi, live) = match f(key) {
            Some((lo, hi)) => (lo, hi, true),
            None => (1, 0, false),
        };
        let idx = self.push(key, lo, hi);
        self.is_hashed[idx] = true;
        self.hashed += 1;
        self.slots[h] = idx as u32 + 1;
        if live {
            Some((idx, lo, hi))
        } else {
            None
        }
    }
    /// Drop dead entries and the hash index (call when construction is done).
    pub fn finish(&mut self) {
        self.slots = Vec::new();
        self.mask = 0;
        self.is_hashed = Vec::new();
        // Compact: drop dead and all-zero states, and trim every window to the
        // range where some moment is nonzero. Blocks were appended in index
        // order, so moving them down never overwrites unread data.
        let nm = self.nm;
        let mut w = 0;
        let mut dst = 0usize;
        for i in 0..self.keys.len() {
            if self.hi[i] < self.lo[i] {
                continue;
            }
            let len = (self.hi[i] - self.lo[i] + 1) as usize;
            let o = self.off[i] as usize;
            let mut first = len;
            let mut last = 0usize;
            for mm in 0..nm {
                let c = &self.arena[o + mm * len..o + (mm + 1) * len];
                if let Some(f) = c.iter().position(|&x| x != 0) {
                    first = first.min(f);
                    last = last.max(c.iter().rposition(|&x| x != 0).unwrap());
                }
            }
            if first == len {
                continue; // all zero
            }
            let nlen = last - first + 1;
            for mm in 0..nm {
                let src = o + mm * len + first;
                self.arena.copy_within(src..src + nlen, dst + mm * nlen);
            }
            let lo_i = self.lo[i];
            self.keys[w] = self.keys[i];
            self.lo[w] = lo_i + first as u8;
            self.hi[w] = lo_i + last as u8;
            self.off[w] = dst as u32;
            dst += nm * nlen;
            w += 1;
        }
        self.arena.truncate(dst);
        self.arena.shrink_to_fit();
        self.keys.truncate(w);
        self.lo.truncate(w);
        self.hi.truncate(w);
        self.off.truncate(w);
        self.keys.shrink_to_fit();
        self.lo.shrink_to_fit();
        self.hi.shrink_to_fit();
        self.off.shrink_to_fit();
    }
    /// Rewrite every key. Only for tables without a hash index (after
    /// `finish`, or built by `push_copy`).
    #[allow(dead_code)]
    pub fn map_keys<F: Fn(K) -> K>(&mut self, f: F) {
        assert!(self.slots.is_empty());
        for k in self.keys.iter_mut() {
            *k = f(*k);
        }
    }
    #[inline(always)]
    pub fn coef(&self, i: usize, mm: usize) -> &[C] {
        let len = (self.hi[i] - self.lo[i] + 1) as usize;
        let o = self.off[i] as usize + mm * len;
        &self.arena[o..o + len]
    }
    #[inline(always)]
    pub fn coef_mut(&mut self, i: usize, mm: usize) -> &mut [C] {
        let len = (self.hi[i] - self.lo[i] + 1) as usize;
        let o = self.off[i] as usize + mm * len;
        &mut self.arena[o..o + len]
    }
    #[inline(always)]
    pub fn accumulate(&mut self, i: usize, mm: usize, off: usize, src: &[C], neg: bool, m: Mod) {
        let dst = &mut self.coef_mut(i, mm)[off..off + src.len()];
        if neg {
            for (d, &s) in dst.iter_mut().zip(src) {
                *d = m.sub(*d, s);
            }
        } else {
            for (d, &s) in dst.iter_mut().zip(src) {
                *d = m.add(*d, s);
            }
        }
    }
}

impl Table<u128> {
    pub fn max_width(&self) -> usize {
        let mut w = 0;
        for &k in &self.keys {
            w = w.max(128 - k.leading_zeros() as usize);
        }
        w
    }
}
