AUTOMAKE_OPTIONS = gnu

# the main target, here the libtool plugin library to build
lib_LTLIBRARIES = libfnh1001.la

# all source files which go into the main target (and have to be distributed)
libfnh1001_la_SOURCES = fnh1001.cc alew.f pdf-cteq6.h pdf-dis-dummy.h

# C++ compiler flags to be used for the main target
#libfnh1001_la_CXXFLAGS = -O0 -Wall -I../../common -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
libfnh1001_la_CXXFLAGS = -O0 -Wall -I${srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}

# Fortran compiler flags to be used for the main target
libfnh1001_la_FFLAGS = -O0 -Wall

# Linker flags to be used for the main target
# libfastnloNlojet in pkglibdir only after make install ...
libfnh1001_la_LDFLAGS = -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6
#libfnh1001_la_LDFLAGS = ../../common/libfastnloNlojet.la -L${LHPDFLIBPATH} -llhpdf_cteq6

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = .deps .obj .libs *.la \
	Makefile.in Makefile.gen Makefile
