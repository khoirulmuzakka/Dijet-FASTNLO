AUTOMAKE_OPTIONS = gnu

# the main target, here the libtool plugin library to build
#lib_LTLIBRARIES = libfnt1001midp.la libfnl1310.la libfnl1308.la
lib_LTLIBRARIES = libfnt1001midp.la

# all source files which go into the main target (and have to be distributed)
libfnt1001midp_la_SOURCES = fnt1001midp.cc cone-et-07.h cone-et-07.cc pdf-cteq6.h pdf-hhc-dummy.h cteq6.h
#libfnl1310_la_SOURCES = fnl1310.cc fj-kt-06.h fj-kt-06.cc pdf-cteq6.h pdf-hhc-dummy.h cteq6.h
#libfnl1308_la_SOURCES = fnl1308.cc fj-sc-07.h fj-sc-07.cc pdf-cteq6.h pdf-hhc-dummy.h cteq6.h

# C++ compiler flags to be used for the main target
libfnt1001midp_la_CXXFLAGS = -g -O0 -Wall -I../../common -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
#libfnl1310_la_CXXFLAGS = -g -O0 -Wall -I../../common -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH} -I${FASTJETINCLUDEPATH}
#libfnl1308_la_CXXFLAGS = -g -O0 -Wall -I../../common -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH} -I${FASTJETINCLUDEPATH}

# Linker flags to be used for the main target
# libfastnloNlojet in pkglibdir only after make install ...
#libfnt1001midp_la_LDFLAGS = -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6
libfnt1001midp_la_LDFLAGS = ../../common/libfastnloNlojet.la -L${LHPDFLIBPATH} -llhpdf_cteq6
#libfnl1310_la_LDFLAGS = ../../common/libfastnloNlojet.la -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet
#libfnl1308_la_LDFLAGS = ../../common/libfastnloNlojet.la -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = .deps .obj .libs *.la \
	Makefile.in Makefile.gen Makefile
