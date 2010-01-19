DIMGRID = 2

AM_CPPFLAGS = $(DUNE_CPPFLAGS) $(DUNEMPICPPFLAGS) -DDIMGRID=$(DIMGRID)
LDADD = $(DUNE_LIBS) $(DUNEMPILIBS)

check_PROGRAMS = decomposition gridcheck seitcheck

decomposition_SOURCES = decomposition.cc
gridcheck_SOURCES = gridcheck.cc
seitcheck_SOURCES = seitcheck.cc

include $(top_srcdir)/am/global-rules
