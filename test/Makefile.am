DIMGRID = 2

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(MPI_CPPFLAGS) -DDIMGRID=$(DIMGRID)
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

check_PROGRAMS = gridcheck

gridcheck_SOURCES = gridcheck.cc

include $(top_srcdir)/am/global-rules

