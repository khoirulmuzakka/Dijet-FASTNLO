# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main targets, here the libtool plugin libraries to build
lib_LTLIBRARIES = libfnh1001.la libfnh5001.la libfnh5002.la libfnhd5102.la libfnhd0001.la
# libfnh2101d.la

#
# All source files which go into the main targets (and have to be distributed)
#
# All targets
AM_SOURCES = pdf-cteq6.h pdf-dis-dummy.h
# Individual sources
libfnh1001_la_SOURCES = fnh1001.cc alew.f $(AM_SOURCES)
libfnh5001_la_SOURCES = fnh5001.cc alew.f $(AM_SOURCES)
libfnh5002_la_SOURCES = fnh5002.cc alew.f $(AM_SOURCES)
libfnhd5102_la_SOURCES = fnhd5102.cc alew.f $(AM_SOURCES)
libfnhd0001_la_SOURCES = fnhd0001.cc alew.f $(AM_SOURCES)
#libfnh2101d_la_SOURCES = fnh2101d.cc alew.f $(AM_SOURCES)

#
# C++ compiler flags to be used for the main target
# top_srcdir for direct compilation before install; pkgincludedir for after install
#
# All targets
AM_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
# Example of individual specification of CXXFLAGS
#libfnh1001_la_CXXFLAGS = -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
#libfnh2101d_la_CXXFLAGS = -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}

#
# Fortran compiler flags to be used for the main target
#
# All targets
AM_FFLAGS = -O0 -Wall
# Example of individual specification of FFLAGS
#libfnh1001_la_FFLAGS = -O0 -Wall
#libfnh2101d_la_FFLAGS = -O0 -Wall

#
# Linker flags to be used for the main target
# top_builddir for distcheck; top_srcdir for direct linking before install; pkglibdir for after install
#
# All targets
#AM_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPAT#H} -llhpdf_cteq6
AM_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6
# Example of individual specification of LDFLAGS
libfnh1001_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6
libfnh5001_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6
libfnh5002_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6
#libfnh2101d_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution.
MAINTAINERCLEANFILES = .deps .obj .libs *.la \
	Makefile.in Makefile.gen Makefile
