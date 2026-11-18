---
name: fortran-style
description: Language level, naming, comments, error handling and file layout for MsPoly Fortran code. Use when writing or editing any .f90/.fypp source or test in this repository.
---

# MsPoly Fortran style

## Language

- Fortran 2008. gfortran builds with `-std=f2008`. Do not use Fortran 2018
  features (variable stop codes, assumed rank, `select rank`, teams).
- `implicit none` everywhere; every module is `private` by default with an
  explicit `public` list.
- Kinds come from `msp_kinds_mod`: `dp` (real/complex storage), `i32` (local
  indices), `i64` (sizes, counts, offsets).
- MPI through `mpi_f08` only. Matrices are distributed over a 2D
  `msp_grid`, with one local piece per rank.
- Types that own resources (`msp_grid`) or contain polymorphic values
  (`msp_dict`) cannot be assigned: copy explicitly (`copy(dst, src)`), and
  release grids with `free`. Do not add final procedures to types with
  default initialization: under `-std=f2008`, gfortran warns that such
  locals are finalized on entry.
- Derived types with polymorphic components follow the
  `polymorphic-fortran` skill.
- Pass `real(z, dp)`, not the part designator `z%re`, to generic
  procedures: nvfortran resolves `z%re` to the complex specific (NV14).
- Select entries with a compacting loop, not `pack` with a mask assigned
  from an array expression: nvfortran -O2 miscomputes the result size
  (NV15).

## Names

- Modules: `msp_<name>_mod`, one per file, file named after the module.
- Public types and procedures: `msp_` prefix (`msp_dict`, `msp_matrix`).
  Public constants: `MSP_` prefix. Do not reuse a type's name in another
  case for a constant: Fortran is case-insensitive.
- Generated sources are `.fypp`. Keep fypp to type-generic repetition;
  everything else is plain Fortran.

## Comments

- Say what a thing is and what it does, only where the code is not already
  clear. Stay minimal.
- When an algorithm or data-structure choice is made, explain the choice
  and why (for example: "entries are sorted by key, so output and hashing
  are deterministic").
- Never write history in code: no "this used to", "X did not work", or
  compiler-bug workaround notes. Compiler knowledge goes in
  `tools/compiler-probes/BUGS.md` (see the `compiler-bugs` skill), and the
  rules it implies go in the `polymorphic-fortran` skill.
- `!>` documents a procedure, type or component. `!!` is an in-body note.

## Errors

- A procedure that can fail takes optional `stat` and `errmsg` arguments,
  like `ALLOCATE`. Call `msp_clear(stat)` on entry and
  `msp_raise(msg, stat, errmsg)` on failure, then return.
- Without `stat`, `msp_raise` prints and stops. Library code never calls
  `stop`, `error stop` or `MPI_Abort` directly.
- Error messages name the full key or object path involved.
- When only some ranks can detect a failure, each rank fills an error
  record (an `msp_dict` with at least `message`) and all ranks call
  `grid%agree(err, stat, errmsg)` before continuing. Every rank then sees
  the same outcome and the same record, so no rank continues alone into a
  collective.

## Layout

- `src/<layer>/` holds one library per layer, with CMake target
  `mspoly_<layer>` and alias `MsPoly::<layer>`. Dependencies only point
  downward.
- Apply compiler flags with `msp_target_flags(<target>)`. Add new fypp
  sources with `msp_fypp(...)`.

## GPU code

- GPU pieces (`msp_gpu_csr_r/c`) hold device pointers (`type(c_ptr)`) to
  memory from `msp_gpu_alloc`: `cudaMalloc` blocks, cached by power-of-two
  size, which CUDA-aware MPI can send directly. Kernel scratch that MPI never
  sees uses the stream-ordered pool (`cudaMallocAsync`).
- Fresh system memory is not an option for GPU data: on Grace-Blackwell it
  costs 5-10 GB/s to create and reads at a third of device bandwidth
  (`tools/compiler-probes/placement`).
- Device memory is released explicitly. Drop every piece with `msp_free`
  (local temporaries, bands, converted copies), and free matrices with
  `free`. Storage arguments that receive a new piece are
  `allocatable, intent(inout)` and are freed first; never `intent(out)`, which
  would deallocate without releasing. Tests check `msp_gpu_allocations() ==
  0` at the end.
- Kernels live in `src/<layer>/*.cu` as `extern "C"` functions over raw,
  1-based device arrays, templated over `double` and `thrust::complex<double>`.
  They return 0 on success; call them through `check(code, what)`.
- When the output size is unknown, use a count pass, allocate, then a fill
  pass.
- Exchange device data with CUDA-aware MPI. NCCL is not used until
  measurements show it is faster.
- Programs call `msp_select_gpu()` before `MPI_Init`. UCX sets up its GPU
  transports at `MPI_Init` only if a CUDA context exists; without one,
  sends from device memory between nodes fail ("cannot find remote
  protocol ... from cuda"). Within one node CUDA IPC hides this.
- Run with `UCC_TL_CUDA_TUNE=0` (HPC-X). Once a CUDA context exists at
  `MPI_Init`, UCC's CUDA collectives take over device `Allgatherv`: about
  2x slower within a node and 1.4x slower across nodes than UCC over UCX.
  Disabling UCC altogether is worse: Open MPI's own device `Allgatherv`
  across nodes is 15x slower.
- Code that only a GPU build compiles sits inside
  `#:if defined('MSPOLY_GPU')`.
