EXTRA_DIST = dune.module

SUBDIRS = doc m4 dune test

DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

