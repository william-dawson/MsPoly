# msp_fypp(<out_var> <source.fypp>...)
# Preprocess fypp sources into the current binary directory and return the
# generated .f90 paths in <out_var>.
find_program(FYPP fypp REQUIRED)

function(msp_fypp out_var)
  set(generated)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fypp)
  foreach(src IN LISTS ARGN)
    get_filename_component(stem ${src} NAME_WE)
    set(out ${CMAKE_CURRENT_BINARY_DIR}/fypp/${stem}.f90)
    add_custom_command(OUTPUT ${out}
      COMMAND ${FYPP} -n ${MSP_FYPP_FLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/${src} ${out}
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${src}
      COMMENT "fypp ${src}")
    list(APPEND generated ${out})
  endforeach()
  set(${out_var} ${generated} PARENT_SCOPE)
endfunction()
