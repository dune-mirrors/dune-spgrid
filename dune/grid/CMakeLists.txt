add_subdirectory(extensions)
add_subdirectory(spgrid)
add_subdirectory(test)

install(FILES spgrid.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid)
