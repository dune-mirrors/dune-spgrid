AM_CPPFLAGS = $(DUNE_CPPFLAGS) $(DUNEMPICPPFLAGS)
LDADD = $(DUNE_LIBS) $(DUNEMPILIBS)

noinst_PROGRAMS = iterators

iterators_SOURCES = iterators.cc

CLEANFILES = vtu/*.vtu vtu/*.pvtu

include $(top_srcdir)/am/global-rules
