set(HEADERS
  checkbndsegiterator.hh
  checkidcommunication.hh
  checkseiterator.hh
  checktree.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/test)

set(GRIDFILES
  1dcube.dgf
  2dcube.dgf
  3dcube.dgf
  4dcube.dgf
  5dcube.dgf
  6dcube.dgf
)

foreach(gridfile ${GRIDFILES})
  configure_file(${gridfile} ${CMAKE_CURRENT_BINARY_DIR}/${gridfile} COPYONLY)
endforeach()

foreach(test bndsegcheck decomposition)
  foreach(dimgrid RANGE 1 6)
    set(testtarget ${test}-${dimgrid})
    add_executable(${testtarget} EXCLUDE_FROM_ALL ${test}.cc)
    target_compile_definitions(${testtarget} PRIVATE "DIMGRID=${dimgrid}")
    add_dune_mpi_flags(${testtarget})
    target_link_dune_default_libraries(${testtarget})
  endforeach()
endforeach()

foreach(test test-spgrid)
  foreach(dimgrid RANGE 1 6)
    dune_add_test(
        NAME ${test}-${dimgrid}
        SOURCES ${test}.cc
        COMPILE_DEFINITIONS "DIMGRID=${dimgrid}"
        TIMEOUT 500
      )
  endforeach()
endforeach()

dune_add_test(SOURCES test-jacobians.cc)
