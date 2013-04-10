set(GRIDDIM  3)

set(TESTS bndsegcheck decomposition test-spgrid)

add_executable(bndsegcheck bndsegcheck.cc)
add_executable(decomposition decomposition.cc)
add_executable(test-spgrid test-spgrid.cc)
set_property(TARGET bndsegcheck decomposition test-spgrid
             APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}")

add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS}) 
