#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# fastNLO toolsfor NLOJet++
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Targets
#
# fastNLO interface library with the jet algos for NLOJet++
pkglib_LTLIBRARIES = libnlojetpptools.la

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
# Common sources
AM_SOURCES =
# Individual sources
libnlojetpptools_la_SOURCES = pdf-hhc-dummy.cc fnlo_int_hhc_nlojet.cc

#
# Preprocessor flags
#
# All targets
# AM_CPPFLAGS =

#
# C++ compiler flags
#
# top_srcdir for direct compilation before install; pkgincludedir for after install
# All targets
AM_CXXFLAGS = -g -O2 -Wall -I${NLOJETINCLUDEPATH} `${FNLOTOOLKITBINPATH}/fnlo-tk-config --cppflags` -I$(srcdir)/../include -I$(top_builddir)/../include -I$(top_srcdir)/../include
# Debug flags
# AM_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall -I${NLOJETINCLUDEPATH} `${FNLOTOOLKITBINPATH}/fnlo-tk-config --cppflags` -I$(srcdir)/../include -I$(top_builddir)/../include -I$(top_srcdir)/../include

#
# Linker flags
# top_builddir for distcheck; top_srcdir for direct linking before install; pkglibdir for after install
#
# All targets
AM_LDFLAGS =

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = .deps .obj .libs *.la \
	Makefile.in Makefile.gen Makefile
