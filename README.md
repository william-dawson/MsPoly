MsPoly
================================================================================
![Build Status](https://github.com/william-dawson/MsPoly/workflows/CI/badge.svg)

MsPoly is a massively parallel library for computing functions of sparse,
Hermitian matrices. It is the second version of
[NTPoly](https://github.com/william-dawson/NTPoly), rewritten from scratch in
pure Fortran 2008 with GPUs in mind. For sufficiently sparse matrices, most of
the matrix functions in MsPoly can be computed in linear time.


Set Up Guide
--------------------------------------------------------------------------------
MsPoly requires:

* A Fortran 2008 compiler (tested with gfortran 13 and newer, and nvfortran).
* An MPI installation (MPI-3, with the `mpi_f08` module).
* CMake (3.24+).
* [fypp](https://github.com/aradi/fypp): a preprocessor, `pip install fypp`.
* ScaLAPACK (with BLACS, LAPACK and BLAS): for the dense eigensolvers.

For GPU support you also need CUDA (nvcc, cuBLAS), cuSOLVERMp and NCCL, and a
CUDA-aware MPI when running on more than one GPU.

To build and test:

> cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFYPP=$(command -v fypp)  
> cmake --build build  
> ctest --test-dir build --output-on-failure

The options are:

* `-DMSPOLY_SCALAPACK_LIBRARIES=` the libraries to link for ScaLAPACK (default
  `scalapack`; separate libraries go in a `;` separated list).
* `-DMSPOLY_GPU=ON` to compute on NVIDIA GPUs. Point `CMAKE_PREFIX_PATH` at
  cuSOLVERMp and NCCL if they are not in a default location.

The libraries are `mspoly_solvers`, `mspoly_matrix` and `mspoly_utils`; link
all three, and add `build/include` for the module files.

Calling Convention
--------------------------------------------------------------------------------
Every public routine looks the same. Matrices are typed arguments. Everything
else is a keyword in an `args` dictionary, and everything the routine has to
say comes back in an `info` dictionary. A parameter you do not set takes its
default, so a call only names what it changes. The same holds for every
function in the library, and for the routines they call: a nested routine's
report appears as a section of the caller's report.

```fortran
program example
  use mpi_f08, only: MPI_Init, MPI_Finalize, MPI_COMM_WORLD
  use msp_kinds_mod, only: dp
  use msp_dict_mod, only: msp_dict
  use msp_grid_mod, only: msp_grid
  use msp_matrix_mod, only: msp_matrix
  use msp_density_mod, only: msp_density
  use msp_exponential_mod, only: msp_exp
  implicit none

  integer, parameter :: n = 200
  type(msp_grid), target :: grid
  type(msp_matrix) :: h, density, e
  type(msp_dict) :: args, info
  integer, allocatable :: rows(:), cols(:)
  real(dp), allocatable :: values(:)
  integer :: i

  call MPI_Init()
  call grid%init(MPI_COMM_WORLD)

  ! A tight-binding chain. Every rank may contribute any entries; here
  ! rank 0 contributes all of them.
  call h%init(grid, n)
  if (grid%rank == 0) then
     rows = [(i, i=1, n), (i, i=1, n-1), (i+1, i=1, n-1)]
     cols = [(i, i=1, n), (i+1, i=1, n-1), (i, i=1, n-1)]
     values = [(2.0_dp, i=1, n), (-1.0_dp, i=1, 2*(n-1))]
  else
     rows = [integer ::]
     cols = [integer ::]
     values = [real(dp) ::]
  end if
  call h%fill(rows, cols, values)

  ! Density matrix with 100 electrons: inputs go in, a report comes out.
  call args%set("trace", 100.0_dp)
  call args%set("method", "trs4")
  call args%set("threshold", 1.0e-8_dp)
  call msp_density(h, density, args, info)
  call grid%log(info)

  ! The same dictionaries carry arguments to any function.
  call args%clear()
  call info%clear()
  call args%set("threshold", 1.0e-8_dp)
  call msp_exp(h, e, args, info)
  call grid%log(info)

  call h%free()
  call density%free()
  call e%free()
  call grid%free()
  call MPI_Finalize()
end program example
```

The report is YAML, printed once by `grid%log`. It holds the result values,
the iteration histories, whether the run converged, and the references to cite
for the methods that were used (the keys match `Citations.bib`):

```yaml
bounds:
  max: 4.00000000000000000E+000
  min: 0.00000000000000000E+000
chemical_potential: 1.99999850988388062E+000
citations:
  kruchinina2016parameterless: true
  niklasson2002: true
converged: true
energy: 7.30381231829969977E+001
iterations: 16
stop: "tolerance"
```

Inside a program, read the same report with `info%get_logical("converged")`,
`info%get_real("energy")` and so on. To move the calculation to the GPU, move
the input with `msp_to_device(h)`; nothing else changes.

Feature Outline
--------------------------------------------------------------------------------
The following functions of a matrix are implemented, for real and complex
matrices, on the CPU and the GPU:

* Density matrix by purification (`msp_density`: `trs2`, `trs4`, `hpcp`, `pm`),
  with a generalized overlap matrix.
* Square root and inverse square root (`msp_sqrt`, `msp_inverse_sqrt`),
  general roots and inverse roots (`msp_root`, `msp_inverse_root`).
* Inverse and pseudo-inverse (`msp_inverse`, `msp_pseudo_inverse`).
* Sign function and polar decomposition (`msp_sign`, `msp_polar`).
* Exponential and logarithm (`msp_exp`, `msp_log`).
* Sine and cosine (`msp_sin`, `msp_cos`).
* Fermi operator (`msp_fermi`).
* General polynomials, in the monomial, Chebyshev or Hermite basis
  (`msp_polynomial`).

All of them take `method = "dense"`, which computes the function from the
eigendecomposition (ScaLAPACK, or cuSOLVERMp on the GPU) instead of by
polynomial expansion or iteration, for matrices that are too full for the
sparse methods to pay off.

The following tools are also available:

* Eigenvalue bounds (`msp_bounds`), and HOMO-LUMO gap estimates
  (`msp_estimate_gap`).
* Linear solvers: conjugate gradient, and Cholesky factorization with pivoting
  (`msp_cg`, `msp_cholesky`, `msp_pivoted_cholesky`).
* Eigendecomposition and singular value decomposition (`msp_eigen`,
  `msp_svd`), and any scalar function of a matrix through the
  eigendecomposition (`msp_dense_function`).
