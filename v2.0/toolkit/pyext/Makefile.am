#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 08.11.2013
#
# fastNLO: C++ version
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

if ENABLE_PYEXT

EXTRA_DIST = fastnlo.i
fastnlo_wrap.cc fastnlo.py: fastnlo.i
	swig -c++ -python -I${srcdir}/../fastnlotoolkit/include/fastnlotk -I${top_builddir}/fastnlotoolkit/include/fastnlotk -I${top_srcdir}/fastnlotoolkit/include/fastnlotk -o fastnlo_wrap.cc $<

AM_CXXFLAGS = -Wno-long-long

all-local: fastnlo_wrap.cc
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

uninstall-local:
	rm -rf $(DESTDIR)$(libdir)/python*/site-packages/fastnlo*

clean-local:
	$(PYTHON) setup.py clean --all
	@rm -f $(top_builddir)/*.pyc
	@rm -f fastnlo_wrap.cc fastnlo.py

endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in fastnlo.py fastnlo_wrap.cc
