# mspoly_matrix: distributed square matrices and their algebra.
msp_fypp(matrix_generated msp_csr_mod.fypp msp_matrix_mod.fypp)

add_library(mspoly_matrix
  msp_storage_mod.f90
  msp_io_mod.f90
  msp_maps_mod.f90
  ${matrix_generated})
add_library(MsPoly::matrix ALIAS mspoly_matrix)
msp_target_flags(mspoly_matrix)
set_target_properties(mspoly_matrix PROPERTIES
  Fortran_MODULE_DIRECTORY ${MSP_MODULE_DIR})
target_link_libraries(mspoly_matrix PUBLIC MsPoly::utils OpenMP::OpenMP_Fortran)
if(MSPOLY_GPU)
  msp_fypp(gpu_generated msp_gpu_csr_mod.fypp)
  target_sources(mspoly_matrix PRIVATE ${gpu_generated} msp_gpu_kernels.cu)
  target_link_libraries(mspoly_matrix PUBLIC CUDA::cudart stdc++)
endif()
