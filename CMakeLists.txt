cmake_minimum_required(VERSION 3.24)
project(MsPoly VERSION 2.0.0 LANGUAGES Fortran)

option(MSPOLY_BUILD_TESTS "Build and register the test suite" ON)
option(MSPOLY_GPU "Compute on NVIDIA GPUs" OFF)
option(MSPOLY_BUILD_BENCH "Build the benchmarks" OFF)
set(MSPOLY_SCALAPACK_LIBRARIES "scalapack" CACHE STRING
  "Libraries for ScaLAPACK, with BLACS, LAPACK and BLAS if separate")

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
if(NOT MPI_Fortran_HAVE_F08_MODULE)
  message(FATAL_ERROR "MsPoly needs an MPI library with the mpi_f08 module")
endif()
set(MSPOLY_MPIEXEC_FLAGS "" CACHE STRING
  "Extra flags for mpiexec in tests (for example --oversubscribe)")

# All module files go to one include directory.
set(MSP_MODULE_DIR ${PROJECT_BINARY_DIR}/include)

add_subdirectory(src)
if(MSPOLY_BUILD_BENCH)
  add_subdirectory(bench)
endif()

if(MSPOLY_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
