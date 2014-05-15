# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main targets, here the libtool plugin libraries to build
# TBD: Better select the "main" targets ...
# TBD: Uncomment here or below to install additional scenarios
#
# No fastjet, have at least for starter
#
lib_LTLIBRARIES =
#
# Require fastjet
#
if HAVE_FASTJET
#   lib_LTLIBRARIES += libfnl2352v22.la libfnl5350eta0v22.la
   lib_LTLIBRARIES += libfnl5350eta0v22.la
endif

#
# All source files which go into the main targets (and have to be distributed)
# Note: All standard jet algo interfaces like fj-ak.cc and fj-ak.h reside in the
#       "common" directory and are collected into the common library libfastnlojetalgos.
#       An example how to add one locally here can be found below.
#
# All targets
AM_SOURCES = pdf-cteq6.h pdf-hhc-dummy.h
# cteq6.h is included from lhpdf installation which is different from original
# Individual sources
libfnl2352v22_la_SOURCES       = fnl2352v22.cc       $(AM_SOURCES)
libfnl5350eta0v22_la_SOURCES   = fnl5350eta0v22.cc   $(AM_SOURCES)

# Example how to add a jet algo interface locally in this directory
# Should only be done temporarily
#libfnl1001ak07_la_SOURCES   = fnl1001ak07.cc   fj-ak.h      fj-ak.cc      $(AM_SOURCES)

#
# Preprocessor flags
#
my_include_dirs = -I${srcdir}/include -I${top_builddir}/include -I${top_srcdir}/include
AM_CPPFLAGS = $(my_include_dirs) -I${includedir} -I${NLOJETINCLUDEPATH}

#
# C++ compiler flags
#
# top_srcdir for direct compilation before install; pkgincludedir for after install
# All targets
#AM_CXXFLAGS = -g -O0 -Wall -I${top_srcdir}/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
AM_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall -I${top_srcdir}/fastnlojetalgos/include/fastnloja -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH} -I${FNLOTOOLKITINCLUDEPATH}
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
#AM_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnloNlojet -lfastnloJetAlgo
AM_LDFLAGS = -L${top_builddir}/interface/hadron/.libs -L${top_srcdir}/interface/hadron
# -L$(pkglibdir) -lfastnlo -lfastnlojetalgos
# -L${LHPDFLIBPATH} -llhpdf_cteq6
if HAVE_FASTJETCONFIG
   AM_LDFLAGS += `${FASTJETBINPATH}/fastjet-config --libs --plugins`
else
if HAVE_FASTJET
   AM_LDFLAGS += -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone -lsiscone_spherical -lCDFConesPlugin -lPxConePlugin -lD0RunIIConePlugin -lNestedDefsPlugin -lTrackJetPlugin -lATLASConePlugin -lCMSIterativeConePlugin -lEECambridgePlugin -lJadePlugin
endif
endif

#v22 force
libfnl2352v22_la_LDFLAGS     = ${AM_LDFLAGS} -lfastnlotoolkit
libfnl5350eta0v22_la_LDFLAGS = ${AM_LDFLAGS} -lfastnlotoolkit

# Example of individual specification of LDFLAGS
#libfnl1001_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6
# Require fastjet, use fastjet-config if possible
#if HAVE_FASTJETCONFIG
#   libfnl1001sc07_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6 `${FASTJETBINPATH}/fastjet-config --libs --plugins`
#else
#   libfnl1001sc07_la_LDFLAGS = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo -L${LHPDFLIBPATH} -llhpdf_cteq6 -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone -lsiscone_spherical -lCDFConesPlugin -lPxConePlugin -lD0RunIIConePlugin -lNestedDefsPlugin -lTrackJetPlugin -lATLASConePlugin -lCMSIterativeConePlugin -lEECambridgePlugin -lJadePlugin
#endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = .deps .obj .libs *.la \
	Makefile.in Makefile.gen Makefile
