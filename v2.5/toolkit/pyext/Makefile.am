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

EXTRA_DIST = fastnlo.i

if ENABLE_PYEXT

python_PYTHON = fastnlo.py
pyexec_LTLIBRARIES = _fastnlo.la

_fastnlo_la_SOURCES = fastnlo_wrap.cc
_fastnlo_la_CPPFLAGS = $(PYTHON_CPPFLAGS) -I$(top_srcdir)/fastnlotoolkit/include/fastnlotk -I$(LHAPDFINCLUDEPATH)
_fastnlo_la_LDFLAGS = -module -shared -avoid-version
_fastnlo_la_LIBADD = $(top_srcdir)/fastnlotoolkit/libfastnlotoolkit.la

include ./$(DEPDIR)/fastnlo_wrap.Plo

$(DEPDIR)/fastnlo_wrap.Plo: fastnlo_wrap.cc

fastnlo_wrap.cc: fastnlo.i
	$(SWIG) -c++ -python -I$(top_srcdir)/fastnlotoolkit/include/fastnlotk -o $@ -MD -MF $(DEPDIR)/fastnlo_wrap.Plo $<

clean-local:
	rm -f fastnlo_wrap.cc fastnlo.py
#$(DEPDIR)/fastnlo_wrap.Plo

endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in fastnlo.py fastnlo_wrap.cc
