GRIDTYPE = ALUGRID_CUBE
#GRIDTYPE = UGGRID

GRIDDIM = 3

AM_CXXFLAGS = -g -Wall -Wfatal-errors

NORMALTESTS = test-cartesiangrid

TESTS = $(NORMALTESTS)

check_PROGRAMS = $(NORMALTESTS)

test_cartesiangrid_SOURCES = test-cartesiangrid.cc
test_cartesiangrid_CPPFLAGS = $(AM_CPPFLAGS) $(ALL_PKG_CPPFLAGS)
test_cartesiangrid_LDFLAGS = $(AM_LDFLAGS) $(ALL_PKG_LDFLAGS)
test_cartesiangrid_LDADD = $(ALL_PKG_LIBS) $(LDADD)

EXTRA_DIST = $(SOURCES)

CLEANFILES = dgfparser.log

include $(top_srcdir)/am/global-rules
