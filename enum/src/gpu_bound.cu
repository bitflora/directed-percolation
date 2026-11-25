// Batched interval-cover bounds on the GPU: a line-for-line port of
// `bound_inner` and `cover_mask` in main.rs, one thread per state.
// 128-bit masks are two u64 words (MSVC hosts have no __int128).
//
// Built by build.rs with nvcc into a static library when the `gpu` feature
// is on. Each calling CPU thread borrows a context (stream and device
// buffers) from a pool, so worker threads run batches concurrently.

#include <cuda_runtime.h>
#include <cstdint>
#include <mutex>
#include <vector>

typedef unsigned long long u64;

struct U128 {
    u64 lo, hi;
};

__device__ __forceinline__ U128 u_and(U128 a, U128 b) { return {a.lo & b.lo, a.hi & b.hi}; }
__device__ __forceinline__ U128 u_or(U128 a, U128 b) { return {a.lo | b.lo, a.hi | b.hi}; }
__device__ __forceinline__ U128 u_not(U128 a) { return {~a.lo, ~a.hi}; }
__device__ __forceinline__ bool u_zero(U128 a) { return (a.lo | a.hi) == 0; }
__device__ __forceinline__ U128 u_shl1(U128 a) { return {a.lo << 1, (a.hi << 1) | (a.lo >> 63)}; }
__device__ __forceinline__ unsigned u_popc(U128 a) { return __popcll(a.lo) + __popcll(a.hi); }
// trailing zeros of a nonzero value
__device__ __forceinline__ unsigned u_tz(U128 a) {
    return a.lo ? (unsigned)(__ffsll((long long)a.lo) - 1) : 64u + (unsigned)(__ffsll((long long)a.hi) - 1);
}
// bits with index >= k (k <= 128; k == 128 gives 0)
__device__ __forceinline__ U128 ge_mask(unsigned k) {
    if (k >= 128) return {0, 0};
    if (k >= 64) return {0, ~0ull << (k - 64)};
    return {~0ull << k, ~0ull};
}

__device__ unsigned cover_mask(U128 fixed, U128 fuzzy, unsigned d) {
    unsigned count = 0;
    for (;;) {
        U128 r = u_or(fixed, u_shl1(fuzzy));
        if (u_zero(r)) return count;
        count += 1;
        unsigned e = u_tz(r) + d;
        if (e >= 127) return count;
        U128 keep = ge_mask(e + 1);  // drop every point whose left end is <= e
        fixed = u_and(fixed, keep);
        fuzzy = u_and(fuzzy, keep);
    }
}

struct Par {
    unsigned j, s, n, tmax;
    int bond, cheap_lo, loose;
};

// (lo, hi) packed as lo | hi << 16; a dead state is (1, 0).
__device__ uint32_t bound_one(U128 key, int lo_known, Par p) {
    const uint32_t DEAD = 1;
    U128 below_j = u_not(ge_mask(p.j));     // bits < j (all if j >= 128)
    U128 upto_j = u_not(ge_mask(p.j + 1));  // bits <= j (all if j >= 127)
    U128 past_fixed = u_and(key, u_not(below_j));
    U128 past_fuzzy = u_and(key, below_j);
    U128 fut_fixed = u_shl1(u_and(key, upto_j));
    U128 fut_fuzzy = u_and(key, u_not(upto_j));
    unsigned nfixed = u_popc(u_and(key, upto_j));
    unsigned nunproc = u_popc(u_and(key, u_not(upto_j)));

    unsigned lo = nfixed;
    unsigned d = 0;
    if (lo_known >= 0) {
        lo = (unsigned)lo_known;
        d = p.s;
    } else if (p.cheap_lo) {
        lo += 2 * p.s;
        d = p.s;
    }
    while (d < p.s) {
        unsigned g = cover_mask(past_fixed, past_fuzzy, d);
        if (g <= 2) {
            lo += 2 * (p.s - d);
            break;
        }
        lo += g;
        d += 1;
    }
    if (lo > p.n) return DEAD;
    unsigned g0 = cover_mask(fut_fixed, fut_fuzzy, 0);
    unsigned extra = p.bond ? nunproc : g0 - (nfixed < g0 ? nfixed : g0);
    unsigned rmin = 0;
    for (;;) {
        unsigned g = rmin == 0 ? g0 : cover_mask(fut_fixed, fut_fuzzy, rmin);
        if (g <= 1) break;
        if (p.bond || rmin >= 1) extra += g;
        rmin += 1;
    }
    if (!p.bond && rmin >= 1) extra += 1;
    if (p.s + 1 + rmin > p.tmax) return DEAD;
    if (p.loose) return 0u | (p.n << 16);
    if (lo + extra > p.n) return DEAD;
    return lo | ((p.n - extra) << 16);
}

__global__ void bound_kernel(const U128* keys, const int* lo_known, uint32_t* out, size_t count, Par p) {
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < count) out[i] = bound_one(keys[i], lo_known[i], p);
}

#include <atomic>
static std::atomic<unsigned long long> stat_kernel_us{0}, stat_h2d_us{0}, stat_d2h_us{0};

struct Ctx {
    cudaStream_t stream = nullptr;
    cudaEvent_t ev[4] = {nullptr, nullptr, nullptr, nullptr};
    U128* d_keys = nullptr;
    int* d_lo = nullptr;
    uint32_t* d_out = nullptr;
    // pinned staging, so copies are asynchronous DMA
    U128* h_keys = nullptr;
    int* h_lo = nullptr;
    uint32_t* h_out = nullptr;
    size_t cap = 0;
};

static std::mutex pool_mu;
static std::vector<Ctx*> pool;

static int grow(Ctx* c, size_t n) {
    if (n <= c->cap) return 0;
    size_t cap = n + n / 4;
    cudaFree(c->d_keys);
    cudaFree(c->d_lo);
    cudaFree(c->d_out);
    cudaFreeHost(c->h_keys);
    cudaFreeHost(c->h_lo);
    cudaFreeHost(c->h_out);
    c->cap = 0;
    cudaError_t e;
    if ((e = cudaMalloc(&c->d_keys, cap * sizeof(U128))) != cudaSuccess) return e;
    if ((e = cudaMalloc(&c->d_lo, cap * sizeof(int))) != cudaSuccess) return e;
    if ((e = cudaMalloc(&c->d_out, cap * sizeof(uint32_t))) != cudaSuccess) return e;
    if ((e = cudaMallocHost(&c->h_keys, cap * sizeof(U128))) != cudaSuccess) return e;
    if ((e = cudaMallocHost(&c->h_lo, cap * sizeof(int))) != cudaSuccess) return e;
    if ((e = cudaMallocHost(&c->h_out, cap * sizeof(uint32_t))) != cudaSuccess) return e;
    c->cap = cap;
    return 0;
}

extern "C" int dp_gpu_bound(const void* keys, const int* lo_known, uint32_t* out, size_t count, unsigned j,
                            unsigned s, unsigned n, unsigned tmax, int bond, int cheap_lo, int loose) {
    if (count == 0) return 0;
    Ctx* c;
    {
        std::lock_guard<std::mutex> g(pool_mu);
        if (pool.empty()) {
            c = new Ctx();
        } else {
            c = pool.back();
            pool.pop_back();
        }
    }
    int err = 0;
    if (!c->stream) {
        err = cudaStreamCreateWithFlags(&c->stream, cudaStreamNonBlocking);
        for (int k = 0; k < 4 && !err; k++) err = cudaEventCreate(&c->ev[k]);
    }
    if (!err) err = grow(c, count);
    if (!err) {
        memcpy(c->h_keys, keys, count * sizeof(U128));
        memcpy(c->h_lo, lo_known, count * sizeof(int));
        cudaEventRecord(c->ev[0], c->stream);
        cudaMemcpyAsync(c->d_keys, c->h_keys, count * sizeof(U128), cudaMemcpyHostToDevice, c->stream);
        cudaMemcpyAsync(c->d_lo, c->h_lo, count * sizeof(int), cudaMemcpyHostToDevice, c->stream);
        cudaEventRecord(c->ev[1], c->stream);
        Par p{j, s, n, tmax, bond, cheap_lo, loose};
        unsigned threads = 128;
        unsigned blocks = (unsigned)((count + threads - 1) / threads);
        bound_kernel<<<blocks, threads, 0, c->stream>>>(c->d_keys, c->d_lo, c->d_out, count, p);
        cudaEventRecord(c->ev[2], c->stream);
        cudaMemcpyAsync(c->h_out, c->d_out, count * sizeof(uint32_t), cudaMemcpyDeviceToHost, c->stream);
        cudaEventRecord(c->ev[3], c->stream);
        err = cudaStreamSynchronize(c->stream);
        float ms[3];
        for (int k = 0; k < 3; k++) cudaEventElapsedTime(&ms[k], c->ev[k], c->ev[k + 1]);
        stat_h2d_us += (unsigned long long)(ms[0] * 1000);
        stat_kernel_us += (unsigned long long)(ms[1] * 1000);
        stat_d2h_us += (unsigned long long)(ms[2] * 1000);
        if (!err) err = cudaGetLastError();
        if (!err) memcpy(out, c->h_out, count * sizeof(uint32_t));
    }
    {
        std::lock_guard<std::mutex> g(pool_mu);
        pool.push_back(c);
    }
    return err;
}

extern "C" const char* dp_gpu_error(int code) { return cudaGetErrorString((cudaError_t)code); }

// Summed stream time (microseconds): host-to-device copies, kernels, device-to-host copies.
extern "C" void dp_gpu_stats(unsigned long long* out3) {
    out3[0] = stat_h2d_us;
    out3[1] = stat_kernel_us;
    out3[2] = stat_d2h_us;
}
