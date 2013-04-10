SUBDIRS = progs

SRCDEPEND = -name \*.hh -o -name \*.png
DOCDEPEND = -name \*.page

doxygendir = $(datadir)/doc/dune-spgrid/

PAGES = main.page
EXTRA_DIST = CMakeLists.txt Doxyfile.in header.html footer.html.in $(PAGES)
CLEANFILES = stamp doxygen.log doxygen.out

if BUILD_DOCS
all: html/index.html
else
all:
endif

# pseudo target to allow writing 'make doc' instead of 'make html/index.html'
doc: html/index.html

FORCE:

stamp: FORCE
	@set -e; \
	if test ! -e stamp \
	|| test x"`find . -maxdepth 1 \( $(DOCDEPEND) \) -a -newer stamp`" != x \
	|| test x"`find $(top_srcdir)/dune \( $(SRCDEPEND) \) -a -newer stamp`" != x; \
	then \
	  touch stamp; \
	fi

html/index.html: Doxyfile stamp
	@echo "Creating documentation (running doxygen) ..."
	@$(DOXYGEN) Doxyfile > doxygen.log 2>&1
	@cat doxygen.out

if BUILD_DOCS
install-data-local: html
	@echo "Installing documentation..."
	@set -e; \
	$(mkinstalldirs) $(DESTDIR)/$(doxygendir); \
	list="$(srcdir)/html/*.html $(srcdir)/html/*.css $(srcdir)/html/*.png $(srcdir)/html/*.gif"; \
	for p in $$list; do \
	  $(install_sh_DATA) $$p $(DESTDIR)/$(doxygendir); \
	done

uninstall-local:
	@echo "Uninstalling documentation..."
	@rm -f $(DESTDIR)/$(doxygendir)/*.html
	@rm -f $(DESTDIR)/$(doxygendir)/*.css
	@rm -f $(DESTDIR)/$(doxygendir)/*.png
	@rm -f $(DESTDIR)/$(doxygendir)/*.gif
endif

clean-local:
	rm -rf html

include $(top_srcdir)/am/global-rules
