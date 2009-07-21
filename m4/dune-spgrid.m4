AC_DEFUN([DUNE_SPGRID_CHECKS])

AC_DEFUN([DUNE_SPGRID_CHECK_MODULE],[
  AC_MSG_NOTICE([Searching for dune-spgrid...])
  DUNE_CHECK_MODULES([dune-spgrid], [grid/spgrid.hh])
])
