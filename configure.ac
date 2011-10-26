AC_PREREQ(2.50)
DUNE_AC_INIT
AM_INIT_AUTOMAKE([foreign 1.5])
AC_CONFIG_SRCDIR([dune/grid/spgrid.hh])
AM_CONFIG_HEADER([config.h])


DUNE_CHECK_ALL

AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  dune/Makefile
  dune/common/Makefile
  dune/grid/Makefile
  dune/grid/cartesiangrid/Makefile
  dune/grid/cartesiangrid/hostgridinfo/Makefile
  dune/grid/extensions/Makefile
  dune/grid/freiburg/Makefile
  dune/grid/spgrid/Makefile
  dune/grid/test/Makefile
  doc/Makefile
  doc/progs/Makefile
  m4/Makefile
  test/Makefile
  dune/grid/spgrid/version.hh
  doc/Doxyfile
  doc/footer.html
  dune-spgrid.pc
])

AC_OUTPUT
DUNE_SUMMARY_ALL
