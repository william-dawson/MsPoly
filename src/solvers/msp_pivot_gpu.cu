// Sparse input columns and a dense thin factor on one GPU. Fortran owns
// the sequential pivot loop; each step uses GEMV for the prior columns.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <memory>
#include <vector>
#include <cub/device/device_reduce.cuh>
#include <cuComplex.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <math_constants.h>

namespace {
struct Maximum {
  double value;
  int64_t index;
};
constexpr int threads = 256;
int blocks(int n) { return (n + threads - 1) / threads; }

int check(int code, const char* what) {
  if (code) std::fprintf(stderr, "msp_pivot_gpu: %s failed (%d)\n", what,
                         code);
  return code;
}
#define TRY(call) do { int rc = check(int(call), #call); \
  if (rc) return rc; } while (0)

struct Context {
  int n;
  bool complex;
  cublasHandle_t blas = nullptr;
  int *ptr = nullptr, *index = nullptr;
  cuDoubleComplex* value = nullptr;
  double* diagonal = nullptr;
  Maximum* maximum = nullptr;
  void *lower = nullptr, *coeff = nullptr, *w = nullptr, *temp = nullptr;
  size_t temp_bytes = 0;
  std::vector<void*> allocations;

  template <class T> int allocate(T*& p, size_t bytes) {
    void* q = nullptr;
    int code = int(cudaMallocAsync(&q, std::max<size_t>(bytes, 1), 0));
    if (code) return code;
    allocations.push_back(q);
    p = static_cast<T*>(q);
    return 0;
  }
  ~Context() {
    for (void* p : allocations) cudaFreeAsync(p, 0);
    cudaStreamSynchronize(0);
    if (blas) cublasDestroy(blas);
  }
};

__global__ void init_diagonal(int n, const int* ptr, const int* index,
                              const cuDoubleComplex* value, double* d) {
  int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= n) return;
  double v = 0;
  for (int q = ptr[row]; q < ptr[row + 1]; ++q)
    if (index[q] == row) v = value[q].x;
  d[row] = v;
}

template <class T> __device__ T cast(cuDoubleComplex x);
template <> __device__ double cast<double>(cuDoubleComplex x) { return x.x; }
template <> __device__ cuDoubleComplex cast<cuDoubleComplex>(
    cuDoubleComplex x) { return x; }
__device__ double conjugate(double x) { return x; }
__device__ cuDoubleComplex conjugate(cuDoubleComplex x) { return cuConj(x); }
__device__ double divide(double x, double d) { return x / d; }
__device__ cuDoubleComplex divide(cuDoubleComplex x, double d) {
  return make_cuDoubleComplex(x.x / d, x.y / d);
}
__device__ double norm_squared(double x) { return x * x; }
__device__ double norm_squared(cuDoubleComplex x) {
  return x.x * x.x + x.y * x.y;
}

template <class T>
__global__ void prepare(int n, int j, int pivot, const T* lower, T* coeff,
                        T* w) {
  int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row < n) w[row] = cast<T>(make_cuDoubleComplex(0, 0));
  if (row < j) coeff[row] = conjugate(lower[pivot + size_t(row) * n]);
}

template <class T>
__global__ void load_column(int pivot, const int* ptr,
                            const int* index, const cuDoubleComplex* value,
                            T* w) {
  int q = ptr[pivot] + blockIdx.x * blockDim.x + threadIdx.x;
  if (q < ptr[pivot + 1]) w[index[q]] = cast<T>(value[q]);
}

template <class T>
__global__ void finish_column(int n, int j, int pivot, double value,
                              const T* w, T* lower, double* d) {
  int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= n) return;
  T v = cast<T>(make_cuDoubleComplex(0, 0));
  if (d[row] != -CUDART_INF) {
    v = divide(w[row], sqrt(value));
    d[row] -= norm_squared(w[row]) / value;
  }
  if (row == pivot) d[row] = -CUDART_INF;
  lower[row + size_t(j) * n] = v;
}

int gemv(cublasHandle_t blas, int n, int j, const double* lower,
          const double* coeff, double* w) {
  const double alpha = -1, beta = 1;
  return int(cublasDgemv(blas, CUBLAS_OP_N, n, j, &alpha, lower, n,
                         coeff, 1, &beta, w, 1));
}
int gemv(cublasHandle_t blas, int n, int j, const cuDoubleComplex* lower,
          const cuDoubleComplex* coeff, cuDoubleComplex* w) {
  const cuDoubleComplex alpha = make_cuDoubleComplex(-1, 0);
  const cuDoubleComplex beta = make_cuDoubleComplex(1, 0);
  return int(cublasZgemv(blas, CUBLAS_OP_N, n, j, &alpha, lower, n,
                         coeff, 1, &beta, w, 1));
}

template <class T> int step(Context& c, int j, int pivot, double value) {
  T* lower = static_cast<T*>(c.lower);
  T* coeff = static_cast<T*>(c.coeff);
  T* w = static_cast<T*>(c.w);
  prepare<<<blocks(c.n), threads>>>(c.n, j, pivot, lower, coeff, w);
  load_column<<<blocks(c.n), threads>>>(pivot, c.ptr, c.index,
                                       c.value, w);
  if (j) TRY(gemv(c.blas, c.n, j, lower, coeff, w));
  finish_column<<<blocks(c.n), threads>>>(c.n, j, pivot, value, w, lower,
                                         c.diagonal);
  return int(cudaGetLastError());
}
}  // namespace

extern "C" int msp_pivot_gpu_open(int n, int rank, int complex,
                                   const int* ptr, const int* index,
                                   const cuDoubleComplex* value,
                                   void** handle) {
  *handle = nullptr;
  std::unique_ptr<Context> c(new Context);
  c->n = n;
  c->complex = complex;
  int count = ptr[n];
  size_t width = complex ? sizeof(cuDoubleComplex) : sizeof(double);
  TRY(cublasCreate(&c->blas));
  TRY(c->allocate(c->ptr, size_t(n + 1) * sizeof(int)));
  TRY(c->allocate(c->index, size_t(count) * sizeof(int)));
  TRY(c->allocate(c->value, size_t(count) * sizeof(cuDoubleComplex)));
  TRY(c->allocate(c->diagonal, size_t(n) * sizeof(double)));
  TRY(c->allocate(c->maximum, sizeof(Maximum)));
  TRY(c->allocate(c->lower, size_t(n) * rank * width));
  TRY(c->allocate(c->coeff, size_t(rank) * width));
  TRY(c->allocate(c->w, size_t(n) * width));
  TRY(cudaMemcpyAsync(c->ptr, ptr, size_t(n + 1) * sizeof(int),
                       cudaMemcpyHostToDevice, 0));
  if (count) {
    TRY(cudaMemcpyAsync(c->index, index, size_t(count) * sizeof(int),
                         cudaMemcpyHostToDevice, 0));
    TRY(cudaMemcpyAsync(c->value, value,
                         size_t(count) * sizeof(cuDoubleComplex),
                         cudaMemcpyHostToDevice, 0));
  }
  init_diagonal<<<blocks(n), threads>>>(n, c->ptr, c->index, c->value,
                                       c->diagonal);
  TRY(cudaGetLastError());
  TRY(cub::DeviceReduce::ArgMax(nullptr, c->temp_bytes, c->diagonal,
                                 &c->maximum->value, &c->maximum->index, n));
  TRY(c->allocate(c->temp, c->temp_bytes));
  *handle = c.release();
  return 0;
}

extern "C" int msp_pivot_gpu_step(void* handle, int column,
                                   double tolerance, int* pivot) {
  Context& c = *static_cast<Context*>(handle);
  *pivot = 0;
  TRY(cub::DeviceReduce::ArgMax(c.temp, c.temp_bytes, c.diagonal,
                                 &c.maximum->value, &c.maximum->index, c.n));
  Maximum maximum;
  TRY(cudaMemcpyAsync(&maximum, c.maximum, sizeof(Maximum),
                       cudaMemcpyDeviceToHost, 0));
  TRY(cudaStreamSynchronize(0));
  if (maximum.value <= tolerance) return 0;
  *pivot = int(maximum.index) + 1;
  if (c.complex) return step<cuDoubleComplex>(c, column, int(maximum.index),
                                             maximum.value);
  return step<double>(c, column, int(maximum.index), maximum.value);
}

extern "C" int msp_pivot_gpu_close(void* handle, int columns,
                                    void* output) {
  std::unique_ptr<Context> c(static_cast<Context*>(handle));
  size_t width = c->complex ? sizeof(cuDoubleComplex) : sizeof(double);
  if (columns) TRY(cudaMemcpyAsync(output, c->lower,
                                    size_t(c->n) * columns * width,
                                    cudaMemcpyDeviceToHost, 0));
  TRY(cudaStreamSynchronize(0));
  return 0;
}
