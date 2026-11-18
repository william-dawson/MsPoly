add_library(msp_testing support/testing_mod.f90)
add_library(msp_spectral support/spectral_mod.f90)
target_link_libraries(msp_spectral PUBLIC MsPoly::matrix)
set_target_properties(msp_spectral PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/tests/include)
msp_target_flags(msp_spectral)
set_target_properties(msp_testing PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/tests/include)
target_include_directories(msp_testing PUBLIC ${PROJECT_BINARY_DIR}/tests/include)

# msp_add_test(<layer> <name> [RANKS <n>...] LIBRARIES <library>...)
# Builds <layer>/<name>.f90, or <layer>/<name>.fypp through fypp. With
# RANKS, registers one MPI run per rank count (<layer>.<name>.np<n>).
function(msp_add_test layer name)
  cmake_parse_arguments(T "" "" "RANKS;LIBRARIES" ${ARGN})
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${layer}/${name}.fypp)
    msp_fypp(source ${layer}/${name}.fypp)
  else()
    set(source ${layer}/${name}.f90)
  endif()
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE msp_testing ${T_LIBRARIES})
  msp_target_flags(${name})
  target_compile_options(${name} PRIVATE ${MSP_TEST_FLAGS})
  if(NOT T_RANKS)
    add_test(NAME ${layer}.${name} COMMAND ${name})
    set_tests_properties(${layer}.${name} PROPERTIES LABELS ${layer})
    return()
  endif()
  separate_arguments(flags NATIVE_COMMAND "${MSPOLY_MPIEXEC_FLAGS}")
  foreach(n IN LISTS T_RANKS)
    add_test(NAME ${layer}.${name}.np${n}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${n} ${flags}
              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${name}> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(${layer}.${name}.np${n} PROPERTIES
      LABELS "${layer};mpi" PROCESSORS ${n} ENVIRONMENT OMP_NUM_THREADS=2)
  endforeach()
endfunction()

foreach(name test_values test_structure test_views test_nested test_citations)
  msp_add_test(utils ${name} LIBRARIES MsPoly::utils)
endforeach()
add_test(NAME utils.test_citations.conflict COMMAND test_citations conflict)
set_tests_properties(utils.test_citations.conflict PROPERTIES LABELS utils
  PASS_REGULAR_EXPRESSION "citations must be a dictionary")
msp_add_test(utils test_grid RANKS 1 2 4 LIBRARIES MsPoly::utils)
msp_add_test(matrix test_matrix RANKS 1 2 4 LIBRARIES MsPoly::matrix)
msp_add_test(matrix test_io RANKS 1 2 4 LIBRARIES MsPoly::matrix msp_spectral)
msp_add_test(solvers test_solvers RANKS 1 2 4 LIBRARIES MsPoly::solvers
  msp_spectral)
msp_add_test(solvers test_functions RANKS 1 2 4 LIBRARIES MsPoly::solvers
  msp_spectral)
msp_add_test(solvers test_eigen RANKS 1 2 4 LIBRARIES MsPoly::solvers
  msp_spectral)
msp_add_test(solvers test_gap RANKS 1 2 4 LIBRARIES MsPoly::solvers
  msp_spectral)
foreach(spec "limit;max_iterations must be positive"
             "tolerance;tolerance must be positive"
             "trace;density.trace must be an interior integer"
             "sampling;sampling_tolerance must be positive"
             "shape;density differs in grid or dimension")
  list(GET spec 0 mode)
  list(GET spec 1 expected)
  separate_arguments(flags NATIVE_COMMAND "${MSPOLY_MPIEXEC_FLAGS}")
  add_test(NAME solvers.test_gap.${mode}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${flags}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_gap> ${mode})
  set_tests_properties(solvers.test_gap.${mode} PROPERTIES
    LABELS solvers PASS_REGULAR_EXPRESSION "${expected}")
endforeach()
# Calls that must stop: (argument, message the run must print).
foreach(spec "basis;unknown basis 'legendre'"
             "method;unknown method 'estrin'"
             "exp_method;unknown method 'krylov'"
             "bounds_method;unknown method 'lanczos'"
             "log_roots;max_square_roots"
             "root;root must be at least 1"
             "cholesky;not positive definite"
             "cholesky_backend;backend 'missing' is not built in")
  list(GET spec 0 mode)
  list(GET spec 1 expected)
  separate_arguments(flags NATIVE_COMMAND "${MSPOLY_MPIEXEC_FLAGS}")
  add_test(NAME solvers.test_functions.${mode}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${flags}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_functions> ${mode})
  set_tests_properties(solvers.test_functions.${mode} PROPERTIES
    LABELS solvers PASS_REGULAR_EXPRESSION "${expected}")
endforeach()
separate_arguments(flags NATIVE_COMMAND "${MSPOLY_MPIEXEC_FLAGS}")
add_test(NAME solvers.test_solvers.method
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${flags}
          ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_solvers> method)
set_tests_properties(solvers.test_solvers.method PROPERTIES LABELS solvers
  PASS_REGULAR_EXPRESSION "unknown method 'trs3'")
foreach(spec "backend;backend 'magma' is not built in"
             "domain;is not positive"
             "trace;trace must lie between"
             "temperature;inverse_temperature must be positive")
  list(GET spec 0 mode)
  list(GET spec 1 expected)
  if(mode STREQUAL "backend" OR MSPOLY_SCALAPACK)
    add_test(NAME solvers.test_eigen.${mode}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${flags}
              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_eigen> ${mode})
    set_tests_properties(solvers.test_eigen.${mode} PROPERTIES
      LABELS solvers PASS_REGULAR_EXPRESSION "${expected}")
  endif()
endforeach()

# Expected stops: (mode, ranks, message the run must print).
add_executable(test_stops utils/test_stops.f90)
target_link_libraries(test_stops PRIVATE MsPoly::utils)
msp_target_flags(test_stops)
foreach(spec "dict;0;dictionary cannot be assigned"
             "raise;0;missing key 'missing'"
             "grid;1;grid cannot be copied"
             "agree;2;message: \"boom\"")
  list(GET spec 0 mode)
  list(GET spec 1 ranks)
  list(GET spec 2 expected)
  if(ranks EQUAL 0)
    add_test(NAME utils.test_stops.${mode} COMMAND test_stops ${mode})
  else()
    separate_arguments(flags NATIVE_COMMAND "${MSPOLY_MPIEXEC_FLAGS}")
    add_test(NAME utils.test_stops.${mode}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${flags}
              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_stops> ${mode})
  endif()
  set_tests_properties(utils.test_stops.${mode} PROPERTIES LABELS utils
    PASS_REGULAR_EXPRESSION "${expected}")
endforeach()

msp_add_test(matrix test_storage LIBRARIES MsPoly::matrix)
set(storage_kinds r c rs cs)
if(MSPOLY_GPU)
  list(APPEND storage_kinds gr gc grs gcs)
endif()
foreach(mode add dot multiply_a multiply_b)
  foreach(kind ${storage_kinds})
    add_test(NAME matrix.test_storage.${mode}_${kind}
      COMMAND test_storage ${mode} ${kind})
    set_tests_properties(matrix.test_storage.${mode}_${kind} PROPERTIES
      LABELS matrix PASS_REGULAR_EXPRESSION "different storage types")
  endforeach()
endforeach()
foreach(kind c cs)
  add_test(NAME matrix.test_storage.real_${kind} COMMAND test_storage real ${kind})
  set_tests_properties(matrix.test_storage.real_${kind} PROPERTIES LABELS matrix
    PASS_REGULAR_EXPRESSION "cannot become real")
endforeach()
add_test(NAME matrix.test_storage.matrix COMMAND test_storage matrix)
set_tests_properties(matrix.test_storage.matrix PROPERTIES LABELS matrix
  PASS_REGULAR_EXPRESSION "matrix cannot be assigned")

if(MSPOLY_GPU)
  # The matrix and solver tests again on the GPU (label gpu).
  separate_arguments(flags NATIVE_COMMAND "${MSPOLY_MPIEXEC_FLAGS}")
  foreach(test matrix/test_matrix solvers/test_solvers solvers/test_functions
                solvers/test_eigen solvers/test_gap)
    get_filename_component(layer ${test} DIRECTORY)
    get_filename_component(name ${test} NAME)
    foreach(n 1 2 4)
      add_test(NAME ${layer}.${name}.gpu.np${n}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${n} ${flags}
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${name}> device)
      set_tests_properties(${layer}.${name}.gpu.np${n} PROPERTIES
        LABELS "${layer};mpi;gpu" PROCESSORS ${n} ENVIRONMENT OMP_NUM_THREADS=2)
    endforeach()
  endforeach()
endif()
# Without GPU support, or with no visible GPU, moving to the device stops.
separate_arguments(flags NATIVE_COMMAND "${MSPOLY_MPIEXEC_FLAGS}")
add_test(NAME matrix.test_storage.device
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${flags}
          ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_storage> device)
set_tests_properties(matrix.test_storage.device PROPERTIES LABELS matrix
  ENVIRONMENT CUDA_VISIBLE_DEVICES=
  PASS_REGULAR_EXPRESSION "built without GPU support|no usable GPU")
