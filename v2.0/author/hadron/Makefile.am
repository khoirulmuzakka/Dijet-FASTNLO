# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main target, here the libtool plugin library to build
lib_LTLIBRARIES = libfnt1001midp.la libfnt2009midp.la libfnt2012test.la libfnt2012.la libfnt2013norm.la libfnt2014norm.la libfnl1001.la
# libfnl1310.la libfnl1308.la libfnt2015norm.la

# All source files which go into the main target (and have to be distributed)
libfnt1001midp_la_SOURCES = fnt1001midp.cc cone-et-07.h cone-et-07.cc pdf-cteq6.h pdf-hhc-dummy.h
libfnt2009midp_la_SOURCES = fnt2009midp.cc cone-e-07.cc pdf-cteq6.h pdf-hhc-dummy.h
libfnt2012test_la_SOURCES = fnt2012test.cc cone-e-07.cc pdf-cteq6.h pdf-hhc-dummy.h
#libfnt2012_la_SOURCES = fnt2012.cc fj-sc-07.h fj-sc-07.cc pdf-cteq6.h pdf-hhc-dummy.h cteq6.h
libfnt2012_la_SOURCES = fnt2012.cc fj-sc-07.h fj-sc-07.cc pdf-cteq6.h pdf-hhc-dummy.h
libfnt2013norm_la_SOURCES = fnt2013norm.cc cone-e-07.cc pdf-cteq6.h pdf-hhc-dummy.h
libfnt2014norm_la_SOURCES = fnt2014norm.cc cone-e-07.cc pdf-cteq6.h pdf-hhc-dummy.h
#libfnt2015norm_la_SOURCES = fnt2015norm.cc cone-e-07.h cone-e-07.cc pdf-cteq6.h pdf-hhc-dummy.h
libfnl1001_la_SOURCES = fnl1001.cc cone-e-07.h cone-e-07.cc pdf-cteq6.h pdf-hhc-dummy.h

#libfnl1310_la_SOURCES = fnl1310.cc fj-kt-06.h fj-kt-06.cc pdf-cteq6.h pdf-hhc-dummy.h cteq6.h
#libfnl1308_la_SOURCES = fnl1308.cc fj-sc-07.h fj-sc-07.cc pdf-cteq6.h pdf-hhc-dummy.h cteq6.h

# C++ compiler flags to be used for the main target
# top_srcdir for direct compilation before install; pkgincludedir for after install
libfnt1001midp_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
libfnt2009midp_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
libfnt2012test_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
libfnt2012_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH} -I${FASTJETINCLUDEPATH}
libfnt2013norm_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
libfnt2014norm_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
# libfnt2015norm_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
libfnl1001_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}

#libfnl1310_la_CXXFLAGS = -g -O0 -Wall -I../../common -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH} -I${FASTJETINCLUDEPATH}
#libfnl1308_la_CXXFLAGS = -g -O0 -Wall -I../../common -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH} -I${FASTJETINCLUDEPATH}

# Linker flags to be used for the main target
# top_builddir for distcheck; top_srcdir for direct linking before install; pkglibdir for after install
libfnt1001midp_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6
libfnt2009midp_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6
libfnt2012test_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6
libfnt2012_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone
libfnt2013norm_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6
libfnt2014norm_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6
# libfnt2015norm_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6
libfnl1001_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6
#libfnl1310_la_LDFLAGS = ../../common/libfastnloNlojet.la -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet
#libfnl1308_la_LDFLAGS = ../../common/libfastnloNlojet.la -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = .deps .obj .libs *.la \
	Makefile.in Makefile.gen Makefile
