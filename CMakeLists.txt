cmake_minimum_required(VERSION 3.24)
project(MsPoly VERSION 2.0.0 LANGUAGES Fortran)

option(MSPOLY_BUILD_TESTS "Build and register the test suite" ON)
option(MSPOLY_BUILD_EXAMPLES "Build the examples" OFF)
option(MSPOLY_GPU "Compute on NVIDIA GPUs" OFF)
option(MSPOLY_BUILD_BENCH "Build the benchmarks" OFF)
set(MSPOLY_SCALAPACK_LIBRARIES "scalapack" CACHE STRING
  "Libraries for ScaLAPACK, with BLACS, LAPACK and BLAS if separate")
set(MSPOLY_SCALAPACK "AUTO" CACHE STRING
  "Dense solvers on the CPU through ScaLAPACK: AUTO (when it links), ON or OFF")
set_property(CACHE MSPOLY_SCALAPACK PROPERTY STRINGS AUTO ON OFF)

# The libraries are static; they can go into shared libraries.
if(NOT DEFINED CMAKE_POSITION_INDEPENDENT_CODE)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(MsPolyCompilers)
include(MsPolyFypp)

find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(OpenMP REQUIRED COMPONENTS Fortran)
set(MSP_FYPP_FLAGS)
if(MSPOLY_GPU)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 100)
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  find_path(CUSOLVERMP_INCLUDE_DIR cusolverMp.h REQUIRED)
  find_library(CUSOLVERMP_LIBRARY cusolverMp PATH_SUFFIXES lib64 REQUIRED)
  find_path(NCCL_INCLUDE_DIR nccl.h REQUIRED)
  find_library(NCCL_LIBRARY nccl PATH_SUFFIXES lib64 REQUIRED)
  list(APPEND MSP_FYPP_FLAGS -DMSPOLY_GPU)
endif()
# ScaLAPACK serves the dense solvers (eigenvalues, and the matrix functions
# through them). Without it MsPoly still builds, with the polynomial solvers,
# and a call to a dense solver stops with a message.
set(MSPOLY_HAVE_SCALAPACK OFF)
if(NOT MSPOLY_SCALAPACK STREQUAL "OFF")
  include(CheckFortranSourceCompiles)
  set(CMAKE_REQUIRED_LIBRARIES ${MSPOLY_SCALAPACK_LIBRARIES} MPI::MPI_Fortran)
  check_fortran_source_compiles("
program probe
  implicit none
  external blacs_gridinit, descinit, pdsyevd
  integer :: context, desc(9), info
  double precision :: a(1), w(1), z(1), work(1)
  call blacs_gridinit(context, 'R', 1, 1)
  call descinit(desc, 1, 1, 1, 1, 0, 0, context, 1, info)
  call pdsyevd('V', 'L', 1, a, 1, 1, desc, w, z, 1, 1, desc, work, 1, &
       & info, 1, info)
end program probe" MSPOLY_SCALAPACK_LINKS SRC_EXT f90)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(MSPOLY_SCALAPACK_LINKS)
    set(MSPOLY_HAVE_SCALAPACK ON)
    list(APPEND MSP_FYPP_FLAGS -DMSPOLY_SCALAPACK)
  elseif(MSPOLY_SCALAPACK STREQUAL "ON")
    message(FATAL_ERROR "ScaLAPACK does not link with MSPOLY_SCALAPACK_LIBRARIES="
      "${MSPOLY_SCALAPACK_LIBRARIES}")
  else()
    message(WARNING "ScaLAPACK does not link with MSPOLY_SCALAPACK_LIBRARIES="
      "${MSPOLY_SCALAPACK_LIBRARIES}: the CPU dense solvers are not built "
      "(set MSPOLY_SCALAPACK=OFF to silence this, or ON to make it an error)")
  endif()
endif()
# The tests of the dense solvers need a dense backend.
if(MSPOLY_HAVE_SCALAPACK OR MSPOLY_GPU)
  set(MSPOLY_HAVE_DENSE ON)
else()
  set(MSPOLY_HAVE_DENSE OFF)
endif()
if(NOT MPI_Fortran_HAVE_F08_MODULE)
  message(FATAL_ERROR "MsPoly needs an MPI library with the mpi_f08 module")
endif()
set(MSPOLY_MPIEXEC_FLAGS "" CACHE STRING
  "Extra flags for mpiexec in tests (for example --oversubscribe)")

# All module files go to one include directory.
set(MSP_MODULE_DIR ${PROJECT_BINARY_DIR}/include)

add_subdirectory(src)
include(MsPolyInstall)
if(MSPOLY_BUILD_BENCH)
  add_subdirectory(bench)
endif()
if(MSPOLY_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(MSPOLY_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# make package_source: the source tarball that the BigDFT installer asks for.
set(CPACK_PACKAGE_NAME "mspoly")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "mspoly-${PROJECT_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES "/\\\\.git/;/\\\\.claude/;/build[^/]*/;\\\\.swp$;~$")
include(CPack)
