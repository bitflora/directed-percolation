//! Batched bounds on the GPU (`--features gpu`, `--gpu`): the kernel in
//! `gpu_bound.cu` is a port of `bound_inner`, linked statically by build.rs.

use crate::{Model, Params};

extern "C" {
    fn dp_gpu_bound(
        keys: *const u128,
        lo_known: *const i32,
        out: *mut u32,
        count: usize,
        j: u32,
        s: u32,
        n: u32,
        tmax: u32,
        bond: i32,
        cheap_lo: i32,
        loose: i32,
    ) -> i32;
    fn dp_gpu_error(code: i32) -> *const std::os::raw::c_char;
    fn dp_gpu_stats(out3: *mut u64);
}

/// Summed GPU stream time in seconds: (host-to-device, kernel, device-to-host).
pub fn stats() -> (f64, f64, f64) {
    let mut v = [0u64; 3];
    unsafe { dp_gpu_stats(v.as_mut_ptr()) };
    (v[0] as f64 / 1e6, v[1] as f64 / 1e6, v[2] as f64 / 1e6)
}

/// Same results as `bound_with` for each key: (lo, hi), or (1, 0) if dead.
pub fn bound(keys: &[u128], lo_known: &[i32], j: usize, s: usize, p: &Params, out: &mut Vec<(u32, u32)>) {
    assert_eq!(keys.len(), lo_known.len());
    let mut packed = vec![0u32; keys.len()];
    let err = unsafe {
        dp_gpu_bound(
            keys.as_ptr(),
            lo_known.as_ptr(),
            packed.as_mut_ptr(),
            keys.len(),
            j as u32,
            s as u32,
            p.n as u32,
            p.tmax as u32,
            (p.model == Model::Bond) as i32,
            p.cheap_lo as i32,
            p.loose as i32,
        )
    };
    if err != 0 {
        let msg = unsafe { std::ffi::CStr::from_ptr(dp_gpu_error(err)) };
        panic!("GPU bound batch failed: {} ({err})", msg.to_string_lossy());
    }
    out.extend(packed.iter().map(|&w| (w & 0xffff, w >> 16)));
}
