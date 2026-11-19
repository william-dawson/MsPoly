// Gathering of device buffers over a communicator with NCCL. The sparse
// products of the GPU matrices collect the pieces of their operands from
// every rank of a row or column of the process grid, and MPI collectives on
// GPU memory are far slower than NCCL's, increasingly so with the number of
// ranks. A communicator is opened once per MPI communicator, from an NCCL
// unique id that the ranks share, and reused.
#include <cuda_runtime.h>
#include <nccl.h>
#include <sched.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

struct Comm {
  cudaStream_t stream = nullptr;
  ncclComm_t nccl = nullptr;
  int rank = 0, size = 0;
};

// Communicators by key.
std::vector<Comm> comms;

int failed(int status, const char* what) {
  if (status) std::fprintf(stderr, "msp_gpu_nccl: %s failed (%d)\n", what,
                           status);
  return status;
}
#define TRY(call) \
  if (failed(int(call), #call)) return 1

// The cores that the process was allowed at start, before an OpenMP runtime
// binds its threads.
cpu_set_t startup_mask;
bool have_startup_mask = false;
__attribute__((constructor)) void remember_startup_mask() {
  have_startup_mask = !sched_getaffinity(0, sizeof(startup_mask),
                                         &startup_mask);
}

// NCCL starts helper threads that poll for the network, and they inherit the
// affinity of the thread that creates the communicator. When the OpenMP
// runtime has bound that thread to a single core (OMP_PROC_BIND), they end up
// on the core that the same thread spins on while it waits for the GPU, and
// every collective crawls. While the communicator is created, the thread
// gets back the cores that the process was given at start (with
// MSP_GPU_NCCL_DEBUG set, the change is reported).
struct WideAffinity {
  cpu_set_t saved;
  bool active = false;

  WideAffinity() {
    CPU_ZERO(&saved);
    if (!have_startup_mask || sched_getaffinity(0, sizeof(saved), &saved))
      return;
    if (CPU_COUNT(&startup_mask) > CPU_COUNT(&saved))
      active = !sched_setaffinity(0, sizeof(startup_mask), &startup_mask);
    if (std::getenv("MSP_GPU_NCCL_DEBUG"))
      std::fprintf(stderr, "msp_gpu_nccl: thread on %d cores, process had %d%s\n",
                   CPU_COUNT(&saved), CPU_COUNT(&startup_mask),
                   active ? ", widened" : "");
  }
  ~WideAffinity() {
    if (active) sched_setaffinity(0, sizeof(saved), &saved);
  }
};

}  // namespace

extern "C" {

// Fill id with a new NCCL unique id (NCCL_UNIQUE_ID_BYTES bytes).
int msp_gpu_nccl_id(char* id) {
  ncclUniqueId unique;
  TRY(ncclGetUniqueId(&unique));
  std::copy(unique.internal, unique.internal + NCCL_UNIQUE_ID_BYTES, id);
  return 0;
}

// Open a communicator for rank of size ranks on the current device, sharing
// the NCCL unique id. Collective over the ranks. A nonzero return (NCCL
// cannot make the communicator, for instance with several ranks on one GPU
// and no multi-rank support) leaves nothing open.
int msp_gpu_nccl_open(const char* id, int size, int rank, int* key) {
  ncclUniqueId unique;
  std::copy(id, id + NCCL_UNIQUE_ID_BYTES, unique.internal);
  Comm c;
  c.rank = rank;
  c.size = size;
  WideAffinity wide;
  TRY(cudaStreamCreate(&c.stream));
  if (failed(int(ncclCommInitRank(&c.nccl, size, unique, rank)),
             "ncclCommInitRank")) {
    cudaStreamDestroy(c.stream);
    return 1;
  }
  comms.push_back(c);
  *key = int(comms.size()) - 1;
  return 0;
}

// Release the communicator key. Collective over its ranks.
void msp_gpu_nccl_close(int key) {
  Comm& c = comms.at(key);
  if (c.stream) cudaStreamSynchronize(c.stream);
  if (c.nccl) ncclCommDestroy(c.nccl);
  if (c.stream) cudaStreamDestroy(c.stream);
  c = Comm{};
}

// Allgatherv of bytes: every rank contributes the counts[rank] bytes at send
// and receives the pieces of all ranks in recv, at displs (in bytes, rank
// order). Both are device pointers. Collective; returns when the data has
// arrived.
int msp_gpu_nccl_gatherv(int key, const void* send, void* recv,
                         const int64_t* counts, const int64_t* displs) {
  Comm& c = comms.at(key);
  char* out = static_cast<char*>(recv);
  // The work already queued on the device has to be finished first.
  TRY(cudaDeviceSynchronize());
  TRY(ncclGroupStart());
  for (int r = 0; r < c.size; ++r) {
    if (counts[r] == 0) continue;
    const void* source = r == c.rank ? send : out + displs[r];
    TRY(ncclBroadcast(source, out + displs[r], size_t(counts[r]), ncclUint8,
                      r, c.nccl, c.stream));
  }
  TRY(ncclGroupEnd());
  TRY(cudaStreamSynchronize(c.stream));
  return 0;
}

}  // extern "C"
