// GPU operations on CSR pieces in device memory. Indices are 1-based, as in
// Fortran. Functions return 0 on success.
//
// Values are stored as T (double, complex, float or complex float) and
// combined in the arithmetic type A (double for real, complex for complex):
// single precision values are widened, and results are rounded once when
// stored.
#include <algorithm>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <cub/block/block_reduce.cuh>
#include <cub/block/block_scan.cuh>
#include <cuda_runtime.h>
#include <thrust/complex.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/scan.h>
#include <thrust/transform_reduce.h>

namespace {

using complex = thrust::complex<double>;
using complexf = thrust::complex<float>;
constexpr int block = 256;

// At least one block: an empty local matrix would otherwise launch a grid of
// zero blocks, an invalid configuration that Thrust later reports as a throw.
int blocks(long n) { return int(std::max<long>(1, (n + block - 1) / block)); }
int status() { return cudaDeviceSynchronize() == cudaSuccess ? 0 : 1; }


__device__ double magnitude(double x) { return fabs(x); }
__device__ double magnitude(complex x) { return thrust::abs(x); }
__device__ double magnitude(float x) { return fabs(double(x)); }
__device__ double magnitude(complexf x) { return thrust::abs(complex(x)); }
__device__ double conjugate(double x) { return x; }
__device__ complex conjugate(complex x) { return thrust::conj(x); }
__device__ float conjugate(float x) { return x; }
__device__ complexf conjugate(complexf x) { return thrust::conj(x); }

template <class T>
constexpr bool is_complex = std::is_same<T, complex>::value ||
                            std::is_same<T, complexf>::value;

// The coefficient re + i im in the arithmetic type A.
template <class A> A coefficient(double re, double im);
template <> double coefficient<double>(double re, double) { return re; }
template <> complex coefficient<complex>(double re, double im) {
  return complex(re, im);
}

// Row r of y + alpha * x, keeping nonzero values with |value| >= threshold.
// The counting pass stores the kept count in out_p[r + 1]; the filling pass
// writes entries from out_p[r].
template <class T, class A>
__global__ void add_rows(int rows, const int* yp, const int* yc, const T* yv,
                         const int* xp, const int* xc, const T* xv, A alpha,
                         double threshold, int* out_p, int* out_c, T* out_v,
                         bool fill) {
  int r = blockIdx.x * blockDim.x + threadIdx.x;
  if (r >= rows) return;
  int p = yp[r] - 1, pe = yp[r + 1] - 1, q = xp[r] - 1, qe = xp[r + 1] - 1;
  int w = fill ? out_p[r] - 1 : 0, kept = 0;
  while (p < pe || q < qe) {
    int c;
    A v;
    if (p >= pe || (q < qe && xc[q] < yc[p])) {
      c = xc[q];
      v = alpha * A(xv[q++]);
    } else if (q >= qe || yc[p] < xc[q]) {
      c = yc[p];
      v = A(yv[p++]);
    } else {
      c = yc[p];
      v = A(yv[p++]) + alpha * A(xv[q++]);
    }
    if (magnitude(v) == 0 || magnitude(v) < threshold) continue;
    if (fill) {
      out_c[w] = c;
      out_v[w++] = T(v);
    }
    ++kept;
  }
  if (!fill) out_p[r + 1] = kept;
}

template <class T> struct row_dot {
  const int *ap, *ac, *bp, *bc;
  const T *av, *bv;
  __device__ complex operator()(int r) const {
    complex sum = 0;
    int p = ap[r] - 1, q = bp[r] - 1;
    while (p < ap[r + 1] - 1 && q < bp[r + 1] - 1) {
      if (ac[p] < bc[q]) ++p;
      else if (bc[q] < ac[p]) ++q;
      else sum += complex(conjugate(av[p++])) * complex(bv[q++]);
    }
    return sum;
  }
};

template <class T> struct row_trace {
  const int *rp, *col;
  const T* val;
  int row0, col0;
  __device__ complex operator()(int r) const {
    complex sum = 0;
    for (int p = rp[r] - 1; p < rp[r + 1] - 1; ++p)
      if (col0 + col[p] == row0 + r + 1) sum += complex(val[p]);
    return sum;
  }
};

template <class T> struct squared {
  __device__ double operator()(T x) const { return magnitude(x) * magnitude(x); }
};

template <class T, class A> __global__ void scale_values(long n, T* v, A alpha) {
  long i = blockIdx.x * long(blockDim.x) + threadIdx.x;
  if (i < n) v[i] = T(alpha * A(v[i]));
}

__global__ void set_ints(int* p, int n, int value) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) p[i] = value;
}

// Band row r of `parts` pieces side by side: piece p's row r, columns
// shifted by shifts[p]. rp holds every piece's row offsets (rows + 1 each),
// starts[p] where its entries begin in col and val.
template <int bytes> struct value_bytes { unsigned char b[bytes]; };

template <class V>
__global__ void concat_offsets(int rows, int parts, const int* rp, int* out_p) {
  int r = blockIdx.x * blockDim.x + threadIdx.x;
  if (r > rows) return;
  int n = 0;
  for (int p = 0; p < parts; ++p) n += rp[p * (rows + 1) + r] - 1;
  out_p[r] = n + 1;
}

template <class V>
__global__ void concat_entries(int rows, int parts, const int* rp,
                               const int* starts, const int* shifts,
                               const int* col, const V* val, const int* out_p,
                               int* out_c, V* out_v) {
  int r = blockIdx.x * blockDim.x + threadIdx.x;
  if (r >= rows) return;
  int w = out_p[r] - 1;
  for (int p = 0; p < parts; ++p) {
    const int* prp = rp + p * (rows + 1);
    for (int q = starts[p] + prp[r] - 1; q < starts[p] + prp[r + 1] - 1; ++q) {
      out_c[w] = col[q] + shifts[p];
      out_v[w++] = val[q];
    }
  }
}

template <class V>
int concat(int rows, int parts, const int* rp, const int* host_starts,
           const int* host_shifts, const int* col, const V* val, int* out_p,
           int* out_c, V* out_v) {
  int* small;
  if (cudaMallocAsync(&small, 2 * parts * sizeof(int), 0) != 0) return 1;
  cudaMemcpy(small, host_starts, parts * sizeof(int), cudaMemcpyDefault);
  cudaMemcpy(small + parts, host_shifts, parts * sizeof(int), cudaMemcpyDefault);
  concat_offsets<V><<<blocks(rows + 1), block>>>(rows, parts, rp, out_p);
  concat_entries<V><<<blocks(rows), block>>>(rows, parts, rp, small,
                                             small + parts, col, val, out_p,
                                             out_c, out_v);
  cudaFreeAsync(small, 0);
  return status();
}

// Band row offsets of pieces stacked in order: piece p's offsets, shifted
// by the entries of the pieces before it. rp holds every piece's offsets
// (heights[p] + 1 each).
__global__ void stack_offsets(int total_rows, int parts, const int* heights,
                              const int* starts, const int* rp, int* out_p) {
  int r = blockIdx.x * blockDim.x + threadIdx.x;
  if (r > total_rows) return;
  int first_row = 0, first_offset = 0;
  for (int p = 0; p < parts; ++p) {
    if (r < first_row + heights[p] || p == parts - 1) {
      out_p[r] = rp[first_offset + r - first_row] + starts[p];
      return;
    }
    first_row += heights[p];
    first_offset += heights[p] + 1;
  }
}

template <class From, class To>
__global__ void convert_values(long n, const From* x, To* y) {
  long i = blockIdx.x * long(blockDim.x) + threadIdx.x;
  if (i < n) y[i] = To(x[i]);
}

// Values of type From converted to the value type numbered `to` (0 double,
// 1 complex, 2 float, 3 complex float). Complex values cannot become real.
template <class From> int convert_from(long n, const From* x, int to, void* y) {
  auto run = [&](auto* out) {
    using To = std::remove_pointer_t<decltype(out)>;
    if constexpr (is_complex<From> && !is_complex<To>) {
      return 1;
    } else {
      if (n > 0) convert_values<<<blocks(n), block>>>(n, x, out);
      return status();
    }
  };
  switch (to) {
    case 0: return run(static_cast<double*>(y));
    case 1: return run(static_cast<complex*>(y));
    case 2: return run(static_cast<float*>(y));
    default: return run(static_cast<complexf*>(y));
  }
}

// out[i] = in[i] + shift.
__global__ void shifted(int n, const int* in, int shift, int* out) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) out[i] = in[i] + shift;
}

__device__ void atomic_add(double* p, double v) { atomicAdd(p, v); }

__device__ void atomic_add(complex* p, complex v) {
  atomicAdd(reinterpret_cast<double*>(p), v.real());
  atomicAdd(reinterpret_cast<double*>(p) + 1, v.imag());
}

__device__ bool keep(double x, double threshold) { return x > 0 && x >= threshold; }

constexpr int row_threads = 128;

// Block r accumulates alpha * a(r0 + r, :) * b into the dense row acc(r, :)
// and marks each touched column in the bitmap touched(r, :). Each warp takes
// its own entries of a, so the warps of a block work on different rows of b
// at once; the lanes of a warp read along one row of b.
template <class T, class A>
__global__ void accumulate(int r0, const int* ap, const int* ac, const T* av,
                           const int* bp, const int* bc, const T* bv, A alpha,
                           int n, int words, A* acc, unsigned* touched) {
  A* row = acc + size_t(blockIdx.x) * n;
  unsigned* mark = touched + size_t(blockIdx.x) * words;
  const int r = r0 + blockIdx.x;
  for (int p = ap[r] - 1 + threadIdx.x / 32; p < ap[r + 1] - 1;
       p += blockDim.x / 32) {
    const int k = ac[p];
    const A a = alpha * A(av[p]);
    for (int q = bp[k - 1] - 1 + threadIdx.x % 32; q < bp[k] - 1; q += 32) {
      const int c = bc[q] - 1;
      atomic_add(row + c, a * A(bv[q]));
      atomicOr(mark + c / 32, 1u << (c % 32));
    }
  }
}

// Block r starts the dense row acc(r, :) as beta * c(r0 + r, :) and marks
// its columns. Columns of a row of c are unique, so plain stores suffice.
template <class T, class A>
__global__ void seed(int r0, const int* cp, const int* cc, const T* cv, A beta,
                     int n, int words, A* acc, unsigned* touched) {
  A* row = acc + size_t(blockIdx.x) * n;
  unsigned* mark = touched + size_t(blockIdx.x) * words;
  const int r = r0 + blockIdx.x;
  for (int q = cp[r] - 1 + threadIdx.x; q < cp[r + 1] - 1; q += blockDim.x) {
    const int c = cc[q] - 1;
    row[c] = beta * A(cv[q]);
    atomicOr(mark + c / 32, 1u << (c % 32));
  }
}

// Number of kept entries of each accumulated row.
template <class A>
__global__ void count_kept(int n, int words, const A* acc,
                           const unsigned* touched, double threshold,
                           int* counts) {
  const A* row = acc + size_t(blockIdx.x) * n;
  const unsigned* mark = touched + size_t(blockIdx.x) * words;
  int kept = 0;
  for (int w = threadIdx.x; w < words; w += blockDim.x)
    for (unsigned bits = mark[w]; bits; bits &= bits - 1)
      kept += keep(magnitude(row[w * 32 + __ffs(bits) - 1]), threshold);
  using reduce = cub::BlockReduce<int, row_threads>;
  __shared__ typename reduce::TempStorage storage;
  const int total = reduce(storage).Sum(kept);
  if (threadIdx.x == 0) counts[blockIdx.x] = total;
}

// Write the kept entries of each accumulated row in column order, starting
// at offsets(r) (0-based), and reset the accumulator and bitmap to zero.
// Each thread owns one bitmap word of a tile; a block scan of the per-word
// counts gives every word its output position.
template <class T, class A>
__global__ void compact(int n, int words, A* acc, unsigned* touched,
                        double threshold, const int* offsets, int* out_c,
                        T* out_v) {
  A* row = acc + size_t(blockIdx.x) * n;
  unsigned* mark = touched + size_t(blockIdx.x) * words;
  using scan = cub::BlockScan<int, row_threads>;
  __shared__ typename scan::TempStorage storage;
  int base = offsets[blockIdx.x];
  for (int tile = 0; tile < words; tile += row_threads) {
    const int w = tile + threadIdx.x;
    const unsigned bits = w < words ? mark[w] : 0;
    int kept = 0;
    for (unsigned b = bits; b; b &= b - 1)
      kept += keep(magnitude(row[w * 32 + __ffs(b) - 1]), threshold);
    int at, total;
    scan(storage).ExclusiveSum(kept, at, total);
    at += base;
    for (unsigned b = bits; b; b &= b - 1) {
      const int c = w * 32 + __ffs(b) - 1;
      if (keep(magnitude(row[c]), threshold)) {
        out_c[at] = c + 1;
        out_v[at++] = T(row[c]);
      }
      row[c] = 0;
    }
    if (w < words) mark[w] = 0;
    base += total;
    __syncthreads();
  }
}

// Rows r0 .. r1 - 1 of a product, held in pool memory until finish.
template <class T> struct chunk {
  int r0, r1;
  long nnz = 0;
  int *rp = nullptr, *col = nullptr;
  T* val = nullptr;
};

// A product computed by start and copied into the caller's arrays by finish.
// Rows are accumulated in blocks whose dense accumulator fits `budget`
// bytes; the bitmap of touched columns keeps counting and compaction
// proportional to the output, not to the row width.
template <class T> struct product {
  std::vector<chunk<T>> chunks;
};

constexpr size_t budget = size_t(1) << 30;

// alpha * a * b + beta * c for the m rows of a, n columns of b, dropping
// exact zeros and |value| < threshold. With beta zero, c is not read. Both
// terms meet in the accumulator, so each result is rounded once.
template <class T, class A>
int product_start(int m, int n, int* ap, int* ac, T* av, int* bp, int* bc,
                  T* bv, A alpha, int* cp, int* cc, T* cv, A beta,
                  double threshold, void** state, long* nnz) {
  auto* s = new product<T>;
  *state = s;
  *nnz = 0;
  const int words = (n + 31) / 32;
  const size_t row_bytes = size_t(n) * sizeof(A) + size_t(words) * 4;
  const int block_rows = int(std::min<size_t>(m, std::max<size_t>(1, budget / row_bytes)));
  A* acc;
  unsigned* touched;
  int* offsets;
  if (cudaMallocAsync(&acc, size_t(block_rows) * n * sizeof(A), 0) != 0 ||
      cudaMallocAsync(&touched, size_t(block_rows) * words * 4, 0) != 0 ||
      cudaMallocAsync(&offsets, (block_rows + 1) * sizeof(int), 0) != 0)
    return 1;
  cudaMemsetAsync(acc, 0, size_t(block_rows) * n * sizeof(A), 0);
  cudaMemsetAsync(touched, 0, size_t(block_rows) * words * 4, 0);
  int code = 0;
  for (int r0 = 0; r0 < m && code == 0; r0 += block_rows) {
    chunk<T> c{r0, std::min(m, r0 + block_rows)};
    const int rows = c.r1 - c.r0;
    int total;
    if (beta != A(0))
      seed<<<rows, row_threads>>>(r0, cp, cc, cv, beta, n, words, acc,
                                  touched);
    accumulate<<<rows, row_threads>>>(r0, ap, ac, av, bp, bc, bv, alpha, n,
                                      words, acc, touched);
    cudaMemsetAsync(offsets, 0, sizeof(int), 0);
    count_kept<<<rows, row_threads>>>(n, words, acc, touched, threshold,
                                      offsets + 1);
    thrust::inclusive_scan(thrust::device, offsets + 1, offsets + rows + 1,
                           offsets + 1);
    cudaMemcpy(&total, offsets + rows, sizeof(int), cudaMemcpyDefault);
    c.nnz = total;
    if (cudaMallocAsync(&c.rp, (rows + 1) * sizeof(int), 0) != 0 ||
        cudaMallocAsync(&c.col, std::max(total, 1) * sizeof(int), 0) != 0 ||
        cudaMallocAsync(&c.val, std::max(total, 1) * sizeof(T), 0) != 0)
      code = 2;
    shifted<<<blocks(rows + 1), block>>>(rows + 1, offsets, 1, c.rp);
    compact<<<rows, row_threads>>>(n, words, acc, touched, threshold, offsets,
                                   c.col, c.val);
    s->chunks.push_back(c);
    *nnz += total;
  }
  cudaFreeAsync(acc, 0);
  cudaFreeAsync(touched, 0);
  cudaFreeAsync(offsets, 0);
  return code ? code : status();
}

template <class T>
int product_finish(void* state, int* cp, int* cc, T* cv) {
  auto* s = static_cast<product<T>*>(state);
  long offset = 0;
  bool ok = true;
  for (auto& c : s->chunks) {
    int rows = c.r1 - c.r0;
    shifted<<<blocks(rows + 1), block>>>(rows + 1, c.rp, int(offset),
                                         cp + c.r0);
    ok = ok &&
         cudaMemcpy(cc + offset, c.col, c.nnz * sizeof(int),
                    cudaMemcpyDefault) == 0 &&
         cudaMemcpy(cv + offset, c.val, c.nnz * sizeof(T),
                    cudaMemcpyDefault) == 0;
    offset += c.nnz;
    cudaFreeAsync(c.rp, 0);
    cudaFreeAsync(c.col, 0);
    cudaFreeAsync(c.val, 0);
  }
  delete s;
  return ok ? status() : 1;
}

template <class T, class A>
int add(int rows, int* yp, int* yc, T* yv, int* xp, int* xc, T* xv, A alpha,
        double threshold, int* out_p, int* out_c, T* out_v, int fill,
        int* total) {
  if (!fill) set_ints<<<1, 1>>>(out_p, 1, 1);
  add_rows<<<blocks(rows), block>>>(rows, yp, yc, yv, xp, xc, xv, alpha,
                                    threshold, out_p, out_c, out_v, fill);
  if (!fill)
    thrust::inclusive_scan(thrust::device, out_p, out_p + rows + 1, out_p);
  cudaMemcpy(total, out_p + rows, sizeof(int), cudaMemcpyDefault);
  *total -= 1;
  return status();
}

template <class T>
void dot(int rows, int* ap, int* ac, T* av, int* bp, int* bc, T* bv,
         double* out) {
  complex d = thrust::transform_reduce(
      thrust::device, thrust::counting_iterator<int>(0),
      thrust::counting_iterator<int>(rows), row_dot<T>{ap, ac, bp, bc, av, bv},
      complex(0), thrust::plus<complex>());
  out[0] = d.real();
  out[1] = d.imag();
}

template <class T>
void trace(int rows, int row0, int col0, int* rp, int* col, T* val,
           double* out) {
  complex t = thrust::transform_reduce(
      thrust::device, thrust::counting_iterator<int>(0),
      thrust::counting_iterator<int>(rows),
      row_trace<T>{rp, col, val, row0, col0}, complex(0),
      thrust::plus<complex>());
  out[0] = t.real();
  out[1] = t.imag();
}

template <class T> double norm2(long n, T* v) {
  return thrust::transform_reduce(thrust::device, v, v + n, squared<T>(), 0.0,
                                  thrust::plus<double>());
}

// Arrays of pieces and bands come from cudaMalloc, which CUDA-aware MPI can
// register (stream-ordered pool memory it cannot). Freed blocks are cached
// by power-of-two size and reused, so steady-state allocation is free.
std::unordered_map<size_t, std::vector<void*>> cached;
std::unordered_map<void*, size_t> block_size;
long live = 0;

size_t bucket(size_t bytes) {
  size_t b = 256;
  while (b < bytes) b *= 2;
  return b;
}

void release_cache() {
  for (auto& [size, list] : cached)
    for (void* p : list) {
      cudaFree(p);
      block_size.erase(p);
    }
  cached.clear();
}

}  // namespace

extern "C" {

void* msp_gpu_alloc(size_t bytes) {
  const size_t b = bucket(bytes);
  void* p = nullptr;
  if (!cached[b].empty()) {
    p = cached[b].back();
    cached[b].pop_back();
  } else if (cudaMalloc(&p, b) != cudaSuccess) {
    release_cache();
    if (cudaMalloc(&p, b) != cudaSuccess) return nullptr;
  }
  block_size[p] = b;
  ++live;
  return p;
}

void msp_gpu_free(void* p) {
  if (!p) return;
  cached[block_size[p]].push_back(p);
  --live;
}

long msp_gpu_allocations() { return live; }

int msp_gpu_set_ints(int* p, int n, int value) {
  set_ints<<<blocks(n), block>>>(p, n, value);
  return status();
}

int msp_gpu_concat_rows(int rows, int parts, int* rp, int* starts, int* shifts,
                        int bytes, int* col, void* val, int* out_p, int* out_c,
                        void* out_v) {
  if (bytes == 4)
    return concat(rows, parts, rp, starts, shifts, col,
                  static_cast<value_bytes<4>*>(val), out_p, out_c,
                  static_cast<value_bytes<4>*>(out_v));
  if (bytes == 8)
    return concat(rows, parts, rp, starts, shifts, col,
                  static_cast<value_bytes<8>*>(val), out_p, out_c,
                  static_cast<value_bytes<8>*>(out_v));
  return concat(rows, parts, rp, starts, shifts, col,
                static_cast<value_bytes<16>*>(val), out_p, out_c,
                static_cast<value_bytes<16>*>(out_v));
}

int msp_gpu_stack_rows(int parts, int* heights, int* starts, int* rp,
                       int* out_p) {
  int total_rows = 0;
  for (int p = 0; p < parts; ++p) total_rows += heights[p];
  int* small;
  if (cudaMallocAsync(&small, 2 * parts * sizeof(int), 0) != 0) return 1;
  cudaMemcpy(small, heights, parts * sizeof(int), cudaMemcpyDefault);
  cudaMemcpy(small + parts, starts, parts * sizeof(int), cudaMemcpyDefault);
  stack_offsets<<<blocks(total_rows + 1), block>>>(total_rows, parts, small,
                                                   small + parts, rp, out_p);
  cudaFreeAsync(small, 0);
  return status();
}

int msp_gpu_init(int local_rank) {
  int count = 0, device;
  if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0) return 1;
  if (cudaSetDevice(local_rank % count) != cudaSuccess) return 1;
  // Work buffers come from the stream-ordered pool, which keeps freed memory
  // for reuse instead of returning it to the driver.
  cudaMemPool_t pool;
  uint64_t keep_all = UINT64_MAX;
  cudaGetDevice(&device);
  cudaDeviceGetDefaultMemPool(&pool, device);
  cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &keep_all);
  return 0;
}

int msp_gpu_prefetch(const void* p, size_t bytes) {
  if (bytes == 0) return 0;
  int device;
  cudaGetDevice(&device);
  cudaMemLocation location{cudaMemLocationTypeDevice, device};
  return cudaMemPrefetchAsync(p, bytes, location, 0, 0) == cudaSuccess ? 0 : 1;
}

int msp_gpu_copy(void* dst, const void* src, size_t bytes) {
  if (bytes == 0) return 0;
  return cudaMemcpy(dst, src, bytes, cudaMemcpyDefault) == cudaSuccess ? 0 : 1;
}

// Values of the type numbered `from` converted to the type numbered `to`
// (0 double, 1 complex, 2 float, 3 complex float).
int msp_gpu_convert(long n, int from, const void* x, int to, void* y) {
  switch (from) {
    case 0: return convert_from(n, static_cast<const double*>(x), to, y);
    case 1: return convert_from(n, static_cast<const complex*>(x), to, y);
    case 2: return convert_from(n, static_cast<const float*>(x), to, y);
    default: return convert_from(n, static_cast<const complexf*>(x), to, y);
  }
}

// The operations of one value type K: values T, arithmetic A.
#define MSP_GPU_KIND(K, T, A)                                                \
  int msp_gpu_scale_##K(long n, T* v, double re, double im) {                \
    if (n > 0)                                                               \
      scale_values<<<blocks(n), block>>>(n, v, coefficient<A>(re, im));      \
    return status();                                                         \
  }                                                                          \
  int msp_gpu_add_##K(int rows, int* yp, int* yc, T* yv, int* xp, int* xc,   \
                      T* xv, double re, double im, double threshold,         \
                      int* out_p, int* out_c, T* out_v, int fill,            \
                      int* total) {                                          \
    return add(rows, yp, yc, yv, xp, xc, xv, coefficient<A>(re, im),         \
               threshold, out_p, out_c, out_v, fill, total);                 \
  }                                                                          \
  void msp_gpu_dot_##K(int rows, int* ap, int* ac, T* av, int* bp, int* bc,  \
                       T* bv, double* out) {                                 \
    dot(rows, ap, ac, av, bp, bc, bv, out);                                  \
  }                                                                          \
  void msp_gpu_trace_##K(int rows, int row0, int col0, int* rp, int* col,    \
                         T* val, double* out) {                              \
    trace(rows, row0, col0, rp, col, val, out);                              \
  }                                                                          \
  double msp_gpu_norm2_##K(long n, T* v) { return norm2(n, v); }             \
  int msp_gpu_product_start_##K(int m, int n, int* ap, int* ac, T* av,       \
                                int* bp, int* bc, T* bv, double re,          \
                                double im, int* cp, int* cc, T* cv,          \
                                double beta_re, double beta_im,              \
                                double threshold, void** state, long* nnz) { \
    return product_start(m, n, ap, ac, av, bp, bc, bv,                       \
                         coefficient<A>(re, im), cp, cc, cv,                 \
                         coefficient<A>(beta_re, beta_im), threshold, state, \
                         nnz);                                               \
  }                                                                          \
  int msp_gpu_product_finish_##K(void* state, int* cp, int* cc, T* cv) {     \
    return product_finish(state, cp, cc, cv);                               \
  }

MSP_GPU_KIND(r, double, double)
MSP_GPU_KIND(c, complex, complex)
MSP_GPU_KIND(rs, float, double)
MSP_GPU_KIND(cs, complexf, complex)

}  // extern "C"
