AC_DEFUN([DUNE_SPGRID_CHECKS],[
  AC_REQUIRE([DUNE_RESOLVED_ABS_TOPSRCDIR])
  DUNE_DEFINE_GRIDTYPE([SPGRID],[GRIDDIM == WORLDDIM],[Dune::SPGrid< double, dimgrid >],[dune/grid/spgrid.hh],[dune/grid/spgrid/dgfparser.hh])
  DUNE_DEFINE_GRIDTYPE([SPGRID_SERIAL],[GRIDDIM == WORLDDIM],[Dune::SPGrid< double, dimgrid, SPIsotropicRefinement, No_Comm >],[dune/grid/spgrid.hh],[dune/grid/spgrid/dgfparser.hh])

  DUNE_DEFINE_GRIDTYPE([SPGRID_ISOTROPIC],[GRIDDIM == WORLDDIM],[Dune::SPGrid< double, dimgrid, SPIsotropicRefinement >],[dune/grid/spgrid.hh],[dune/grid/spgrid/dgfparser.hh])
  DUNE_DEFINE_GRIDTYPE([SPGRID_ANISOTROPIC],[GRIDDIM == WORLDDIM],[Dune::SPGrid< double, dimgrid, SPAnisotropicRefinement >],[dune/grid/spgrid.hh],[dune/grid/spgrid/dgfparser.hh])
  DUNE_DEFINE_GRIDTYPE([SPGRID_BISECTION],[GRIDDIM == WORLDDIM],[Dune::SPGrid< double, dimgrid, SPBisectionRefinement >],[dune/grid/spgrid.hh],[dune/grid/spgrid/dgfparser.hh])
  
  DUNE_DEFINE_GRIDTYPE([SPGRID_ISOTROPIC_SERIAL],[GRIDDIM == WORLDDIM],[Dune::SPGrid< double, dimgrid, SPIsotropicRefinement, No_Comm >],[dune/grid/spgrid.hh],[dune/grid/spgrid/dgfparser.hh])
  DUNE_DEFINE_GRIDTYPE([SPGRID_ANISOTROPIC_SERIAL],[GRIDDIM == WORLDDIM],[Dune::SPGrid< double, dimgrid, SPAnisotropicRefinement, No_Comm >],[dune/grid/spgrid.hh],[dune/grid/spgrid/dgfparser.hh])
  DUNE_DEFINE_GRIDTYPE([SPGRID_BISECTION_SERIAL],[GRIDDIM == WORLDDIM],[Dune::SPGrid< double, dimgrid, SPBisectionRefinement, No_Comm >],[dune/grid/spgrid.hh],[dune/grid/spgrid/dgfparser.hh])
])

AC_DEFUN([DUNE_SPGRID_CHECK_MODULE],[
  AC_MSG_NOTICE([Searching for dune-spgrid...])
  DUNE_CHECK_MODULES([dune-spgrid],[grid/spgrid.hh])
])
