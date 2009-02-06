AUTOMAKE_OPTIONS = foreign

# the main target, here the libtool plugin library to build
lib_LTLIBRARIES = libfnt1001midp.la

# all source files which go into the main target (and have to be distributed)
libfnt1001midp_la_SOURCES = fnt1001midp.cc cone-et-07.cc

# C++ compiler flags to be used for the main target
libfnt1001midp_la_CXXFLAGS = -O0 -Wall -I../../common -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}

# Linker flags to be used for the main target
#libfnt1001midp_la_LDFLAGS = -L${HOME}/local/lib/fastnlo -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6
libfnt1001midp_la_LDFLAGS = -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = .deps .obj .libs \
	Makefile.in Makefile.gen Makefile
