# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main target, here the libtool plugin library to build
lib_LTLIBRARIES = libfnh1001.la
# libfnh2101d.la

# All source files which go into the main target (and have to be distributed)
libfnh1001_la_SOURCES = fnh1001.cc alew.f pdf-cteq6.h pdf-dis-dummy.h
#libfnh2101d_la_SOURCES = fnh2101d.cc alew.f pdf-cteq6.h pdf-photo-dummy.h

# C++ compiler flags to be used for the main target
# top_srcdir for direct compilation before install; pkgincludedir for after install
libfnh1001_la_CXXFLAGS = -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
#libfnh2101d_la_CXXFLAGS = -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}

# Fortran compiler flags to be used for the main target
libfnh1001_la_FFLAGS = -O0 -Wall
#libfnh2101d_la_FFLAGS = -O0 -Wall

# Linker flags to be used for the main target
# top_builddir for distcheck; top_srcdir for direct linking before install; pkglibdir for after install
libfnh1001_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6
#libfnh2101d_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution.
MAINTAINERCLEANFILES = .deps .obj .libs *.la \
	Makefile.in Makefile.gen Makefile
