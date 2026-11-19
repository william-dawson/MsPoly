# Install the libraries, the module files and a pkg-config file, whose link
# line carries what the static libraries need: ScaLAPACK, OpenMP and, with the
# GPU, the CUDA libraries.
include(GNUInstallDirs)

install(TARGETS mspoly_utils mspoly_matrix mspoly_solvers
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ${MSP_MODULE_DIR}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mspoly
  FILES_MATCHING PATTERN "*.mod")

set(MSP_PC_LIBS ${OpenMP_Fortran_FLAGS})
if(MSPOLY_HAVE_SCALAPACK)
  foreach(lib IN LISTS MSPOLY_SCALAPACK_LIBRARIES)
    if(lib MATCHES "^[-/]")
      string(APPEND MSP_PC_LIBS " ${lib}")
    else()
      string(APPEND MSP_PC_LIBS " -l${lib}")
    endif()
  endforeach()
endif()
if(MSPOLY_GPU)
  foreach(lib IN ITEMS ${CUSOLVERMP_LIBRARY} ${NCCL_LIBRARY})
    get_filename_component(dir ${lib} DIRECTORY)
    get_filename_component(name ${lib} NAME_WE)
    string(REGEX REPLACE "^lib" "" name ${name})
    string(APPEND MSP_PC_LIBS " -L${dir} -l${name}")
  endforeach()
  # The CUDA libraries, and the stub of the driver library to link against.
  string(APPEND MSP_PC_LIBS " -L${CUDAToolkit_LIBRARY_DIR}"
    " -L${CUDAToolkit_LIBRARY_DIR}/stubs -lcublas -lcudart -lcuda -lstdc++")
endif()

configure_file(${PROJECT_SOURCE_DIR}/cmake/mspoly.pc.in
  ${PROJECT_BINARY_DIR}/mspoly.pc @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/mspoly.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
