DIMGRID = 2

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(MPI_CPPFLAGS) -DDIMGRID=$(DIMGRID)
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS) $(MPI_LIBS)

check_PROGRAMS = gridcheck seitcheck

gridcheck_SOURCES = gridcheck.cc
seitcheck_SOURCES = seitcheck.cc

include $(top_srcdir)/am/global-rules

