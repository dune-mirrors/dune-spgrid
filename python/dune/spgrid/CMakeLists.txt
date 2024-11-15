add_python_targets(spgrid
  __init__
)

dune_add_pybind11_module(NAME _spgrid SOURCES _spgrid.cc)
if( DUNE_ENABLE_PYTHONMODULE_PRECOMPILE )
  set(DIMS 1 2 3 4)
  foreach( DIM ${DIMS} )
    add_library(spgrid_${DIM} OBJECT registergrid.cc)
    target_compile_definitions(spgrid_${DIM} PRIVATE DIM=${DIM})
    dune_target_enable_all_packages( spgrid_${DIM} )

    target_link_libraries(_spgrid PRIVATE spgrid_${DIM})
  endforeach()
endif()

if(SKBUILD)
  install(TARGETS _spgrid LIBRARY DESTINATION python/dune/spgrid)
endif()
