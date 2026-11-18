# Compiler flags shared by every MsPoly target.
#   MSP_FLAGS        always
#   MSP_DEBUG_FLAGS  Debug builds
#   MSP_TEST_FLAGS   test programs only
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Exact comparisons with zero are intended (a zero coefficient means the
  # term is skipped), so -Wcompare-reals is off.
  set(MSP_FLAGS -std=f2008 -fimplicit-none -Wall -Wextra
    -Wno-maybe-uninitialized -Wno-compare-reals)
  set(MSP_DEBUG_FLAGS -fcheck=all -ffpe-trap=zero,overflow -g)
  set(MSP_TEST_FLAGS)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
  set(MSP_FLAGS -Mstandard)
  set(MSP_DEBUG_FLAGS -Mbounds -Mchkptr -g)
  set(MSP_TEST_FLAGS)
endif()

# msp_target_flags(<target>) applies the flags above to Fortran sources.
function(msp_target_flags target)
  target_compile_options(${target} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${MSP_FLAGS}>
    $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:${MSP_DEBUG_FLAGS}>)
endfunction()
