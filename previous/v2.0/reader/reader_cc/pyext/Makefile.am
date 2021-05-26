#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 08.11.2013
#
# fastNLO_reader: C++ version
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

if ENABLE_PYEXT

EXTRA_DIST = fastnloreader.i
fastnloreader_wrap.cpp fastnloreader.py: fastnloreader.i
	swig -c++ -python -I$(top_srcdir)/include -o fastnloreader_wrap.cpp $<

AM_CXXFLAGS = -Wno-long-long

all-local: fastnloreader_wrap.cpp
	$(PYTHON) setup.py build

install-exec-local:
	$(PYTHON) setup.py install --prefix=$(DESTDIR)$(prefix)

uninstall-local:
	rm -rf $(DESTDIR)$(libdir)/python*/site-packages/lhapdf*

clean-local:
	$(PYTHON) setup.py clean --all
	@rm -f $(top_builddir)/*.pyc

endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in fastnloreader.py fastnloreader_wrap.cpp
