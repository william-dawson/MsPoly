# mspoly_utils: kinds, errors, keyword dictionaries, process grid, logging.
msp_fypp(utils_generated msp_dict_mod.fypp)

add_library(mspoly_utils
  msp_kinds_mod.f90
  msp_error_mod.f90
  msp_grid_mod.f90
  msp_citations_mod.f90
  ${utils_generated})
add_library(MsPoly::utils ALIAS mspoly_utils)
msp_target_flags(mspoly_utils)
set_target_properties(mspoly_utils PROPERTIES
  Fortran_MODULE_DIRECTORY ${MSP_MODULE_DIR})
target_include_directories(mspoly_utils PUBLIC
  $<BUILD_INTERFACE:${MSP_MODULE_DIR}>)
target_link_libraries(mspoly_utils PUBLIC MPI::MPI_Fortran)
