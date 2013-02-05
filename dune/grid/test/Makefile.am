GRIDDIM = 3

AM_CXXFLAGS = -g -Wall -Wfatal-errors
AM_CPPFLAGS = $(DUNE_CPPFLAGS) $(DUNEMPICPPFLAGS) -DDIMGRID=$(GRIDDIM)
LDADD = $(DUNE_LIBS) $(DUNEMPILIBS)

check_PROGRAMS = decomposition test-spgrid

decomposition_SOURCES = decomposition.cc
test_spgrid_SOURCES = test-spgrid.cc

DGFFILES = 2dcube.dgf  2dtorus.dgf  3dcube.dgf  4dcube.dgf

EXTRA_DIST = checkidcommunication.hh  checkseiterator.hh $(DGFFILES)

CLEANFILES = dgfparser.log *.spgrid *.vtu

include $(top_srcdir)/am/global-rules
