AC_DEFUN([DUNE_SPGRID_CHECKS],[
  DUNE_DEFINE_GRIDTYPE([SPGRID],[GRIDDIM == WORLDDIM],[Dune::SPGrid< double, dimgrid >],[dune/grid/spgrid.hh],[dune/grid/spgrid/dgfparser.hh])
  DUNE_DEFINE_GRIDTYPE([SPGRID_ISOTROPIC],[GRIDDIM == WORLDDIM],[Dune::SPGrid< double, dimgrid, SPIsotropicRefinement >],[dune/grid/spgrid.hh],[dune/grid/spgrid/dgfparser.hh])
  DUNE_DEFINE_GRIDTYPE([SPGRID_ANISOTROPIC],[GRIDDIM == WORLDDIM],[Dune::SPGrid< double, dimgrid, SPAnisotropicRefinement >],[dune/grid/spgrid.hh],[dune/grid/spgrid/dgfparser.hh])
  DUNE_DEFINE_GRIDTYPE([SPGRID_BISECTION],[GRIDDIM == WORLDDIM],[Dune::SPGrid< double, dimgrid, SPBisectionRefinement >],[dune/grid/spgrid.hh],[dune/grid/spgrid/dgfparser.hh])
])

AC_DEFUN([DUNE_SPGRID_CHECK_MODULE],[
  AC_MSG_NOTICE([Searching for dune-spgrid...])
  DUNE_CHECK_MODULES([dune-spgrid],[grid/spgrid.hh])
])
