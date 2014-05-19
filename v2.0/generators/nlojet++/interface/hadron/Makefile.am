#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 19.05.2014
#
# fastNLO interface library with scenarios for NLOJet++
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Targets
#
# Without fastjet
#
lib_LTLIBRARIES = libfnl5350eta0v22conee.la
#
# Require fastjet
#
if HAVE_FASTJET
   lib_LTLIBRARIES += libfnl5350eta0v22new.la
endif

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
# Common sources
AM_SOURCES =
# Individual sources
libfnl5350eta0v22conee_la_SOURCES = fnl5350eta0v22conee.cc cone-e.cc $(AM_SOURCES)
libfnl5350eta0v22new_la_SOURCES   = fnl5350eta0v22.cc      fj-ak.cc  $(AM_SOURCES)

#
# Preprocessor flags
#
#my_include_dirs = -I${srcdir}/fastnlojetalgos/include/fastnloja -I${top_builddir}/fastnlojetalgos/include/fastnloja -I${top#_srcdir}/fastnlojetalgos/include/fastnloja
#AM_CPPFLAGS = $(my_include_dirs) -I${includedir} -I${NLOJETINCLUDEPATH}

#
# C++ compiler flags
#
# top_srcdir for direct compilation before install; pkgincludedir for after install
# All targets
AM_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall -I${NLOJETINCLUDEPATH} -I${FNLOTOOLKITINCLUDEPATH} -I../include

# -I${top_srcdir}/../../common -I${pkgincludedir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH}
#AM_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall -I${top_srcdir}/fastnlojetalgos/include/fastnloja
# -I${pkginclud#edir} -I${NLOJETINCLUDEPATH} -I${LHPDFINCLUDEPATH} -I${FNLOTOOLKITINCLUDEPATH} -I${FNLOTOOLKITINCLUDEPATH}
if HAVE_FASTJETCONFIG
   AM_CXXFLAGS += `${FASTJETBINPATH}/fastjet-config --plugins --cxxflags`
else
if HAVE_FASTJET
   AM_CXXFLAGS += -I${FASTJETINCLUDEPATH}
endif
endif

#
# Linker flags to be used for the main target
# top_builddir for distcheck; top_srcdir for direct linking before install; pkglibdir for after install
#
# All targets

AM_LDFLAGS =
if HAVE_FASTJETCONFIG
   AM_LDFLAGS += `${FASTJETBINPATH}/fastjet-config --libs --plugins`
else
if HAVE_FASTJET
   AM_LDFLAGS += -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone -lsiscone_spherical -lCDFConesPlugin -lPxConePlugin -lD0RunIIConePlugin -lNestedDefsPlugin -lTrackJetPlugin -lATLASConePlugin -lCMSIterativeConePlugin -lEECambridgePlugin -lJadePlugin
endif
endif

#AM_LDFLAGS += -L$(pkglibdir) -lfastnlo -lfastnloJetAlgos
# AM_LDFLAGS += -L${top_builddir}/../../common/.libs -L${top_srcdir}/../../common -L$(pkglibdir) -lfastnlo -lfastnloJetAlgos
# -L${LHPDFLIBPATH} -llhpdf_cteq6
#AM_LDFLAGS =

#AM_LDFLAGS += -L${top_builddir}/interface/hadron/.libs -L${top_builddir}/fastnlojetalgos/.libs -L${top_srcdir}/interface/hadron -L${top_srcdir}/fastnlojetalgos -L$(pkglibdir)
# -L$(pkglibdir) -lfastnlo -lfastnlojetalgos
# -L${LHPDFLIBPATH} -llhpdf_cteq6
#v22 force
libfnl5350eta0v22new_la_LDFLAGS = ${AM_LDFLAGS} -L${FNLOTOOLKITLIBPATH} -lfastnlotoolkit
libfnl5350eta0v22conee_la_LDFLAGS = ${AM_LDFLAGS} -L${FNLOTOOLKITLIBPATH} -lfastnlotoolkit


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
