AC_PREREQ(2.50)
DUNE_AC_INIT
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([dune/grid/spgrid.hh])
AM_CONFIG_HEADER([config.h])


DUNE_CHECK_ALL

AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  dune/Makefile
  dune/grid/Makefile
  dune/grid/spgrid/Makefile
  doc/Makefile
  m4/Makefile
  dune-spgrid.pc
])

AC_OUTPUT
DUNE_SUMMARY_ALL
