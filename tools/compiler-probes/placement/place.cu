// Placement of fresh system-memory arrays on a GB200 node: read bandwidth
// after (a) GPU first touch, (b) CPU fill + prefetch, (c) copy from device.
#include <cuda_runtime.h>
#include <thrust/execution_policy.h>
#include <thrust/reduce.h>
extern "C" {
static float elapsed(cudaEvent_t a, cudaEvent_t b) { float ms; cudaEventElapsedTime(&ms, a, b); return ms; }
__global__ void fill(double* x, long n) { long i = blockIdx.x * long(blockDim.x) + threadIdx.x; if (i < n) x[i] = 1.0; }
int place_init(int rank) { int c; cudaGetDeviceCount(&c); cudaSetDevice(rank % c); int d; cudaGetDevice(&d); return d; }
// Seconds per read of n doubles (average of 10).
double place_read(const double* x, long n) {
  cudaEvent_t a, b; cudaEventCreate(&a); cudaEventCreate(&b);
  thrust::reduce(thrust::device, x, x + n, 0.0);
  cudaEventRecord(a);
  for (int r = 0; r < 10; ++r) thrust::reduce(thrust::device, x, x + n, 0.0);
  cudaEventRecord(b); cudaEventSynchronize(b);
  return elapsed(a, b) / 1000.0 / 10;
}
double place_gpu_write(double* x, long n) {
  cudaEvent_t a, b; cudaEventCreate(&a); cudaEventCreate(&b);
  cudaEventRecord(a); fill<<<(n + 255) / 256, 256>>>(x, n); cudaEventRecord(b); cudaEventSynchronize(b);
  return elapsed(a, b) / 1000.0;
}
double place_prefetch(double* x, long n) {
  int d; cudaGetDevice(&d); cudaMemLocation l{cudaMemLocationTypeDevice, d};
  cudaEvent_t a, b; cudaEventCreate(&a); cudaEventCreate(&b);
  cudaEventRecord(a); cudaMemPrefetchAsync(x, n * 8, l, 0, 0); cudaEventRecord(b); cudaEventSynchronize(b);
  return elapsed(a, b) / 1000.0;
}
double place_copy_from_device(double* x, long n) {
  double* d; cudaMalloc(&d, n * 8); fill<<<(n + 255) / 256, 256>>>(d, n); cudaDeviceSynchronize();
  cudaEvent_t a, b; cudaEventCreate(&a); cudaEventCreate(&b);
  cudaEventRecord(a); cudaMemcpy(x, d, n * 8, cudaMemcpyDefault); cudaEventRecord(b); cudaEventSynchronize(b);
  cudaFree(d);
  return elapsed(a, b) / 1000.0;
}
}
