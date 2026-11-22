//! Arithmetic modulo a prime below 2^62, values stored as u64.
//!
//! One 62-bit prime carries as many bits as four 15-bit primes (Jensen 1999
//! used 16-bit residues), so a series needs a quarter as many runs for the
//! same total memory-time.

/// Coefficient type. The default stores residues modulo 62-bit primes in u64.
/// `--features u32` stores residues modulo 31-bit primes in u32: half the memory
/// per coefficient, twice as many primes. Prefer it when memory, not time, binds.
/// `--features u16` stores residues modulo 15-bit primes in u16 (as Jensen
/// 1999 did): a quarter of the u64 coefficient memory, four times the primes.
#[cfg(all(not(feature = "u32"), not(feature = "u16")))]
pub type C = u64;
#[cfg(feature = "u32")]
pub type C = u32;
#[cfg(feature = "u16")]
pub type C = u16;

#[cfg(all(feature = "u32", feature = "u16"))]
compile_error!("features u32 and u16 are mutually exclusive");

#[cfg(all(not(feature = "u32"), not(feature = "u16")))]
const PRIME_BITS: u32 = 62;
#[cfg(feature = "u32")]
const PRIME_BITS: u32 = 31;
#[cfg(feature = "u16")]
const PRIME_BITS: u32 = 15;

#[derive(Clone, Copy, Debug)]
pub struct Mod {
    pub p: C,
}

impl Mod {
    pub fn new(p: u64) -> Self {
        assert!(p < (1u64 << PRIME_BITS));
        Mod { p: p as C }
    }
    #[inline(always)]
    pub fn add(&self, a: C, b: C) -> C {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    #[inline(always)]
    pub fn sub(&self, a: C, b: C) -> C {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }
    #[inline(always)]
    pub fn neg(&self, a: C) -> C {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }
    #[inline(always)]
    pub fn mul(&self, a: C, b: C) -> C {
        ((a as u128 * b as u128) % self.p as u128) as C
    }
    pub fn from_i64(&self, x: i64) -> C {
        (x as i128).rem_euclid(self.p as i128) as C
    }
    pub fn pow(&self, mut a: C, mut e: u64) -> C {
        let mut r = 1;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }
    pub fn inv(&self, a: C) -> C {
        assert!(a != 0);
        self.pow(a, self.p as u64 - 2)
    }
    pub fn series_mul(&self, a: &[C], b: &[C]) -> Vec<C> {
        let n = a.len().min(b.len());
        let mut out = vec![0; n];
        for i in 0..n {
            if a[i] == 0 {
                continue;
            }
            for k in 0..(n - i) {
                out[i + k] = self.add(out[i + k], self.mul(a[i], b[k]));
            }
        }
        out
    }
    /// Inverse of a power series with nonzero constant term.
    pub fn series_inv(&self, a: &[C]) -> Vec<C> {
        let n = a.len();
        let mut inv = vec![0; n];
        let c = self.inv(a[0]);
        inv[0] = c;
        for k in 1..n {
            let mut s = 0;
            for i in 1..=k {
                s = self.add(s, self.mul(a[i], inv[k - i]));
            }
            inv[k] = self.mul(self.neg(s), c);
        }
        inv
    }
}

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for &sp in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % sp == 0 {
            return n == sp;
        }
    }
    let mut d = n - 1;
    let mut r = 0;
    while d % 2 == 0 {
        d /= 2;
        r += 1;
    }
    let mulmod = |a: u64, b: u64| ((a as u128 * b as u128) % n as u128) as u64;
    let powmod = |mut a: u64, mut e: u64| {
        let mut res = 1u64;
        a %= n;
        while e > 0 {
            if e & 1 == 1 {
                res = mulmod(res, a);
            }
            a = mulmod(a, a);
            e >>= 1;
        }
        res
    };
    // these bases are deterministic for all n < 3.3e24
    'outer: for &a in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 0..r - 1 {
            x = mulmod(x, x);
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

/// The k largest primes below 2^PRIME_BITS, in decreasing order.
pub fn primes(k: usize) -> Vec<u64> {
    let mut out = Vec::new();
    let mut c: u64 = (1u64 << PRIME_BITS) - 1;
    while out.len() < k {
        if is_prime(c) {
            out.push(c);
        }
        c -= 1;
    }
    out
}
