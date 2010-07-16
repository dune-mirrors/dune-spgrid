EXTRA_DIST = dune.module

SUBDIRS = doc m4 dune test

DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT)  CXX="$(CXX)" CC="$(CC)"

SOURCEFILE_MASK = -name \*.h -o -name \*.hh -o -name \*.c -o -name \*.cc

# add library name, authors list, and license to source files
dist-hook:
	@for SOURCEFILE in `find $(distdir) \( $(SOURCEFILE_MASK) \) -print` ; do \
	  TMPFILE=`mktemp authors.XXXXXX` ; \
	  echo "/**************************************************************************" > $$TMPFILE ; \
	  echo "  " >> $$TMPFILE ; \
	  cat $(top_srcdir)/AUTHORS >> $$TMPFILE ; \
	  echo "  " >> $$TMPFILE ; \
	  echo "**************************************************************************/" >> $$TMPFILE ; \
	  cat $$SOURCEFILE >> $$TMPFILE ; \
	  cp $$TMPFILE $$SOURCEFILE ; \
	  rm -f $$TMPFILE ; \
	done

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules
