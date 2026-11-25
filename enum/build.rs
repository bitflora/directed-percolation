// With `--features gpu`, compile src/gpu_bound.cu with nvcc into a static
// library and link it with the static CUDA runtime. Nothing to do otherwise.
//
// Environment: CUDA_PATH (set by the CUDA installer); NVCC_CCBIN, the
// directory of a cl.exe that nvcc supports (defaults to the newest Visual
// Studio 2022 MSVC found), and CUDA_ARCH (default sm_86).

use std::path::PathBuf;
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    if std::env::var_os("CARGO_FEATURE_GPU").is_none() {
        return;
    }
    println!("cargo:rerun-if-changed=src/gpu_bound.cu");
    println!("cargo:rerun-if-env-changed=NVCC_CCBIN");
    println!("cargo:rerun-if-env-changed=CUDA_ARCH");
    let cuda = PathBuf::from(std::env::var("CUDA_PATH").expect("CUDA_PATH not set: install the CUDA toolkit"));
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let arch = std::env::var("CUDA_ARCH").unwrap_or_else(|_| "sm_86".into());
    let lib = out.join("dpgpu.lib");
    let mut cmd = Command::new(cuda.join("bin").join("nvcc"));
    cmd.args(["-O3", "-lib", "-Xcompiler", "/MD,/O2", &format!("-arch={arch}"), "-o"]).arg(&lib).arg("src/gpu_bound.cu");
    if let Some(ccbin) = std::env::var_os("NVCC_CCBIN").map(PathBuf::from).or_else(find_msvc_2022) {
        cmd.arg("-ccbin").arg(ccbin);
    }
    let status = cmd.status().expect("failed to run nvcc");
    assert!(status.success(), "nvcc failed");
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=dpgpu");
    println!("cargo:rustc-link-search=native={}", cuda.join("lib").join("x64").display());
    println!("cargo:rustc-link-lib=static=cudart_static");
}

/// The newest MSVC host-x64 compiler directory under Visual Studio 2022.
fn find_msvc_2022() -> Option<PathBuf> {
    for edition in ["Community", "Professional", "Enterprise", "BuildTools"] {
        let root = PathBuf::from(format!(r"C:\Program Files\Microsoft Visual Studio\2022\{edition}\VC\Tools\MSVC"));
        let Ok(rd) = std::fs::read_dir(&root) else { continue };
        let mut vers: Vec<PathBuf> = rd.filter_map(|e| e.ok()).map(|e| e.path()).collect();
        vers.sort();
        if let Some(v) = vers.pop() {
            let bin = v.join(r"bin\Hostx64\x64");
            if bin.join("cl.exe").exists() {
                return Some(bin);
            }
        }
    }
    None
}
