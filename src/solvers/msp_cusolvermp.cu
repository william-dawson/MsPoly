// Hermitian eigendecomposition with cuSOLVERMp. The caller holds the matrix
// in ScaLAPACK's block-cyclic layout (square blocks, source process 0, a
// row-major process grid) on host memory; this file moves it to the GPU,
// solves, and copies the eigenvalues and eigenvectors back.
//
// A context (stream, NCCL communicator, solver handle and process grid) is
// costly to create, so the caller opens one per process grid, built from an
// NCCL unique id that it shares among the ranks, and reuses it.
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cusolverMp.h>
#include <nccl.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <type_traits>
#include <vector>

namespace {

struct Context {
  cudaStream_t stream = nullptr;
  ncclComm_t nccl = nullptr;
  cusolverMpHandle_t handle = nullptr;
  cusolverMpGrid_t grid = nullptr;
  // Two ints on the device for agreeing on a status across the ranks.
  int* flag = nullptr;
};

// Contexts by key; a closed context stays as an empty slot.
std::vector<Context> contexts;

// Device memory of one solve, released on return.
struct Buffers {
  cudaStream_t stream;
  cusolverMpMatrixDescriptor_t desc = nullptr;
  void *a = nullptr, *z = nullptr, *work = nullptr;
  double* w = nullptr;
  int* info = nullptr;

  explicit Buffers(cudaStream_t s) : stream(s) {}
  ~Buffers() {
    cudaStreamSynchronize(stream);
    for (void* p : {a, z, work, static_cast<void*>(w),
                    static_cast<void*>(info)})
      if (p) cudaFree(p);
    if (desc) cusolverMpDestroyMatrixDesc(desc);
  }
};

// Nonzero when a call failed, after printing which one.
int failed(int status, const char* what) {
  if (status) std::fprintf(stderr, "msp_cusolvermp: %s failed (%d)\n", what,
                           status);
  return status;
}
#define TRY(call) \
  if (failed(int(call), #call)) return 1

void release(Context& c) {
  if (c.stream) cudaStreamSynchronize(c.stream);
  if (c.grid) cusolverMpDestroyGrid(c.grid);
  if (c.handle) cusolverMpDestroy(c.handle);
  if (c.nccl) ncclCommDestroy(c.nccl);
  if (c.flag) cudaFree(c.flag);
  if (c.stream) cudaStreamDestroy(c.stream);
  c = Context{};
}

// The largest status over the ranks of the context, so that a failure on
// one rank stops them all before the next collective call.
int agree(const Context& c, int status, int* worst) {
  TRY(cudaMemcpyAsync(c.flag, &status, sizeof(int), cudaMemcpyHostToDevice,
                      c.stream));
  TRY(ncclAllReduce(c.flag, c.flag + 1, 1, ncclInt, ncclMax, c.nccl,
                    c.stream));
  TRY(cudaMemcpyAsync(worst, c.flag + 1, sizeof(int), cudaMemcpyDeviceToHost,
                      c.stream));
  TRY(cudaStreamSynchronize(c.stream));
  return 0;
}

// Descriptor, device buffers and workspace sizes for one solve, with the
// matrix on its way to the device. Nothing here is collective.
template <class T>
int prepare(const Context& c, Buffers& b, cudaDataType type, char jobz,
            int64_t n, int nb, int64_t lld, size_t bytes, const T* a,
            size_t* device_bytes, size_t* host_bytes) {
  TRY(cusolverMpCreateMatrixDesc(&b.desc, c.grid, type, n, n, nb, nb, 0, 0,
                                 lld));
  TRY(cudaMalloc(&b.a, bytes));
  TRY(cudaMalloc(&b.z, bytes));
  TRY(cudaMalloc(&b.w, size_t(n) * sizeof(double)));
  TRY(cudaMalloc(&b.info, sizeof(int)));
  TRY(cudaMemcpyAsync(b.a, a, bytes, cudaMemcpyHostToDevice, c.stream));
  TRY(cusolverMpSyevd_bufferSize(c.handle, &jobz, CUBLAS_FILL_MODE_LOWER, n,
                                 b.a, 1, 1, b.desc, b.w, b.z, 1, 1, b.desc,
                                 type, device_bytes, host_bytes));
  TRY(cudaMalloc(&b.work, std::max<size_t>(*device_bytes, 1)));
  return 0;
}

// Without vectors, z is not touched.
template <class T>
int eigh(int key, bool vectors, int64_t n, int nb, int64_t lld,
         int64_t local_cols, const T* a, double* w, T* z) {
  constexpr cudaDataType type =
      std::is_same<T, double>::value ? CUDA_R_64F : CUDA_C_64F;
  const size_t bytes = size_t(lld) * size_t(local_cols) * sizeof(T);
  char jobz = vectors ? 'V' : 'N';
  const Context& c = contexts.at(key);
  Buffers b(c.stream);
  size_t device_bytes = 0, host_bytes = 0;
  int worst = 0;
  TRY(agree(c, prepare(c, b, type, jobz, n, nb, lld, bytes, a, &device_bytes,
                       &host_bytes),
            &worst));
  if (worst) return 1;
  std::vector<char> host_work(std::max<size_t>(host_bytes, 1));
  TRY(cusolverMpSyevd(c.handle, &jobz, CUBLAS_FILL_MODE_LOWER, n, b.a, 1, 1,
                      b.desc, b.w, b.z, 1, 1, b.desc, type, b.work,
                      device_bytes, host_work.data(), host_bytes, b.info));
  int info = 0;
  TRY(cudaMemcpyAsync(&info, b.info, sizeof(int), cudaMemcpyDeviceToHost,
                      c.stream));
  TRY(cudaMemcpyAsync(w, b.w, size_t(n) * sizeof(double),
                      cudaMemcpyDeviceToHost, c.stream));
  if (vectors)
    TRY(cudaMemcpyAsync(z, b.z, bytes, cudaMemcpyDeviceToHost, c.stream));
  TRY(cudaStreamSynchronize(c.stream));
  return failed(info, "cusolverMpSyevd (info)");
}

}  // namespace

extern "C" {

// Fill id with a new NCCL unique id (NCCL_UNIQUE_ID_BYTES bytes).
int msp_cusolvermp_id(char* id) {
  ncclUniqueId unique;
  TRY(ncclGetUniqueId(&unique));
  std::copy(unique.internal, unique.internal + NCCL_UNIQUE_ID_BYTES, id);
  return 0;
}

// Open a context for rank of size ranks on a rows x cols process grid on
// the current device, sharing the NCCL unique id; key identifies it.
// Collective over the ranks.
int msp_cusolvermp_open(const char* id, int size, int rank, int rows,
                        int cols, int* key) {
  ncclUniqueId unique;
  int device;
  std::copy(id, id + NCCL_UNIQUE_ID_BYTES, unique.internal);
  auto slot = std::find_if(contexts.begin(), contexts.end(),
                           [](const Context& c) {
                             return !c.stream && !c.nccl && !c.handle &&
                                    !c.grid;
                           });
  *key = int(slot - contexts.begin());
  if (slot == contexts.end()) contexts.push_back(Context{});
  Context& s = contexts.at(*key);
  struct Guard {
    Context& context;
    bool complete = false;
    ~Guard() { if (!complete) release(context); }
  } guard{s};
  TRY(cudaGetDevice(&device));
  TRY(cudaStreamCreate(&s.stream));
  TRY(ncclCommInitRank(&s.nccl, size, unique, rank));
  TRY(cusolverMpCreate(&s.handle, device, s.stream));
  TRY(cusolverMpCreateDeviceGrid(s.handle, &s.grid, s.nccl, rows, cols,
                                 CUSOLVERMP_GRID_MAPPING_ROW_MAJOR));
  TRY(cudaMalloc(&s.flag, 2 * sizeof(int)));
  guard.complete = true;
  return 0;
}

// Release the context key. Collective over its ranks.
void msp_cusolvermp_close(int key) { release(contexts.at(key)); }

int msp_cusolvermp_eigh_r(int key, int vectors, int64_t n, int nb, int64_t lld,
                          int64_t local_cols, const double* a, double* w,
                          double* z) {
  return eigh(key, vectors != 0, n, nb, lld, local_cols, a, w, z);
}

int msp_cusolvermp_eigh_c(int key, int vectors, int64_t n, int nb, int64_t lld,
                          int64_t local_cols, const cuDoubleComplex* a,
                          double* w, cuDoubleComplex* z) {
  return eigh(key, vectors != 0, n, nb, lld, local_cols, a, w, z);
}

}  // extern "C"
