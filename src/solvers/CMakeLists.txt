# mspoly_solvers: matrix functions built on the matrix layer.
msp_fypp(solvers_generated msp_eigen_mod.fypp msp_linear_mod.fypp)

add_library(mspoly_solvers
  msp_bounds_mod.f90
  msp_convergence_mod.f90
  msp_polynomial_mod.f90
  msp_sqrt_mod.f90
  msp_inverse_mod.f90
  msp_sign_mod.f90
  msp_exponential_mod.f90
  msp_trigonometry_mod.f90
  msp_density_mod.f90
  msp_gap_mod.f90
  msp_fermi_mod.f90
  msp_extrapolate_mod.f90
  msp_analysis_mod.f90
  msp_svd_mod.f90
  ${solvers_generated})
add_library(MsPoly::solvers ALIAS mspoly_solvers)
msp_target_flags(mspoly_solvers)
set_target_properties(mspoly_solvers PROPERTIES
  Fortran_MODULE_DIRECTORY ${MSP_MODULE_DIR})
target_link_libraries(mspoly_solvers PUBLIC MsPoly::matrix)
if(MSPOLY_HAVE_SCALAPACK)
  target_link_libraries(mspoly_solvers PRIVATE ${MSPOLY_SCALAPACK_LIBRARIES})
endif()
if(MSPOLY_GPU)
  target_sources(mspoly_solvers PRIVATE msp_pivot_gpu.cu msp_cusolvermp.cu)
  target_link_libraries(mspoly_solvers PRIVATE CUDA::cublas)
  target_include_directories(mspoly_solvers PRIVATE ${CUSOLVERMP_INCLUDE_DIR}
    ${NCCL_INCLUDE_DIR})
  target_link_libraries(mspoly_solvers PRIVATE ${CUSOLVERMP_LIBRARY}
    ${NCCL_LIBRARY} CUDA::cudart CUDA::cuda_driver)
endif()
