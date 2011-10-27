GRIDTYPE = ALUGRID_CUBE
#GRIDTYPE = UGGRID

GRIDDIM = 3

AM_CXXFLAGS = -g -Wall -Wfatal-errors
AM_CPPFLAGS = $(DUNE_CPPFLAGS) $(DUNEMPICPPFLAGS) -DDIMGRID=$(GRIDDIM)
LDADD = $(DUNE_LIBS) $(DUNEMPILIBS)

NORMALTESTS = test-cartesiangrid

TESTS = $(NORMALTESTS)

check_PROGRAMS = $(NORMALTESTS) decomposition seitcheck test-spgrid

test_cartesiangrid_SOURCES = test-cartesiangrid.cc
test_cartesiangrid_CPPFLAGS = $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
test_cartesiangrid_LDFLAGS = $(AM_LDFLAGS) $(ALL_PKG_LDFLAGS)
test_cartesiangrid_LDADD = $(ALL_PKG_LIBS) $(LDADD)

decomposition_SOURCES = decomposition.cc
test_spgrid_SOURCES = test-spgrid.cc
seitcheck_SOURCES = seitcheck.cc

EXTRA_DIST = checkidcommunication.hh checkseiterator.hh

CLEANFILES = dgfparser.log *.spgrid *.vtu

include $(top_srcdir)/am/global-rules
