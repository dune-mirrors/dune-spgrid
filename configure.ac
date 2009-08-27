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
  dune/grid/Makefile
  dune/grid/extensions/Makefile
  dune/grid/io/Makefile
  dune/grid/io/file/Makefile
  dune/grid/io/file/dgfparser/Makefile
  dune/grid/spgrid/Makefile
  dune/grid/utility/Makefile
  doc/Makefile
  m4/Makefile
  test/Makefile
  dune-spgrid.pc
])

AC_OUTPUT
DUNE_SUMMARY_ALL
