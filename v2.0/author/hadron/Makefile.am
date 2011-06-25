# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main targets, here the libtool plugin libraries to build
lib_LTLIBRARIES = libfnt1001midp.la libfnt2009midp.la libfnt2012test.la libfnt2013norm.la libfnt2014norm.la libfnl1001.la
#
# Require fastjet
#
if HAVE_FASTJET
   lib_LTLIBRARIES += libfnl2222b.la libfnl2342b.la libfnl2412c.la libfnl2822num.la libfnl2822den.la libfnl1001ak07.la libfnl1001sc07.la libfnl1001d0r207.la libfnt2012.la
endif

#
# All source files which go into the main targets (and have to be distributed)
# Note: All standard jet algo interfaces like fj-ak.cc and fj-ak.h reside in the
#       "common" directory and are collected into the common library libfastnloJetAlgos.
#       An example how to add one locally here can be found below. 
#
# All targets
AM_SOURCES = pdf-cteq6.h pdf-hhc-dummy.h
# Individual sources
libfnt1001midp_la_SOURCES   = fnt1001midp.cc   $(AM_SOURCES)
libfnt2009midp_la_SOURCES   = fnt2009midp.cc   $(AM_SOURCES)
libfnt2012test_la_SOURCES   = fnt2012test.cc   $(AM_SOURCES)
libfnt2013norm_la_SOURCES   = fnt2013norm.cc   $(AM_SOURCES)
libfnt2014norm_la_SOURCES   = fnt2014norm.cc   $(AM_SOURCES)
#libfnt2015norm_la_SOURCES   = fnt2015norm.cc   $(AM_SOURCES)
libfnl1001_la_SOURCES       = fnl1001.cc       $(AM_SOURCES)
# Require fastjet
libfnl1001ak07_la_SOURCES   = fnl1001ak07.cc   $(AM_SOURCES)
libfnl1001sc07_la_SOURCES   = fnl1001sc07.cc   $(AM_SOURCES)
libfnl1001d0r207_la_SOURCES = fnl1001d0r207.cc $(AM_SOURCES)
libfnt2012_la_SOURCES       = fnt2012.cc       $(AM_SOURCES)
libfnl2222b_la_SOURCES      = fnl2222b.cc      $(AM_SOURCES)
libfnl2342b_la_SOURCES      = fnl2342b.cc      $(AM_SOURCES)
libfnl2412c_la_SOURCES      = fnl2412c.cc      $(AM_SOURCES)
libfnl2822den_la_SOURCES    = fnl2822den.cc    $(AM_SOURCES)
libfnl2822num_la_SOURCES    = fnl2822num.cc    $(AM_SOURCES)
# Example how to add a jet algo interface locally in this directory
# Should only be done temporarily
#libfnl1001ak07_la_SOURCES   = fnl1001ak07.cc   fj-ak.h      fj-ak.cc      $(AM_SOURCES)

#
# C++ compiler flags to be used for the main target
# top_srcdir for direct compilation before install; pkgincludedir for after install
#
# All targets
AM_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
if HAVE_FASTJETCONFIG 
   AM_CXXFLAGS += `${FASTJETBINPATH}/fastjet-config --plugins --cxxflags`
else
   if HAVE_FASTJET
      AM_CXXFLAGS += -I${FASTJETINCLUDEPATH}
   endif
endif

# Example of individual specification of CXXFLAGS
#libfnl1001_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
# Require fastjet, use fastjet-config if possible
#if HAVE_FASTJETCONFIG 
#   libfnl1001sc07_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH} `${FASTJETBINPATH}/fastjet-config --plugins --cxxflags`
#else
#   libfnl1001sc07_la_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH} -I${FASTJETINCLUDEPATH} 
#endif

#
# Linker flags to be used for the main target
# top_builddir for distcheck; top_srcdir for direct linking before install; pkglibdir for after install
#
# All targets
AM_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -lfastnloJetAlgos -L${LHPDFLIBPATH} -llhpdf_cteq6
if HAVE_FASTJETCONFIG 
   AM_LDFLAGS += `${FASTJETBINPATH}/fastjet-config --libs --plugins`
else
   if HAVE_FASTJET
      AM_LDFLAGS += -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone -lsiscone_spherical -lCDFConesPlugin -lPxConePlugin -lD0RunIIConePlugin -lNestedDefsPlugin -lTrackJetPlugin -lATLASConePlugin -lCMSIterativeConePlugin -lEECambridgePlugin -lJadePlugin
   endif
endif

# Example of individual specification of LDFLAGS
#libfnl1001_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6
# Require fastjet, use fastjet-config if possible
#if HAVE_FASTJETCONFIG 
#   libfnl1001sc07_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6 `${FASTJETBINPATH}/fastjet-config --libs --plugins`
#else
#   libfnl1001sc07_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone -lsiscone_spherical -lCDFConesPlugin -lPxConePlugin -lD0RunIIConePlugin -lNestedDefsPlugin -lTrackJetPlugin -lATLASConePlugin -lCMSIterativeConePlugin -lEECambridgePlugin -lJadePlugin
#endif

#libfnl1310_la_LDFLAGS = ../../common/libfastnloNlojet.la -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet
#libfnl1308_la_LDFLAGS = ../../common/libfastnloNlojet.la -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = .deps .obj .libs *.la \
	Makefile.in Makefile.gen Makefile
