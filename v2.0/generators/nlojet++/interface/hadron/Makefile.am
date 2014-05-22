#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 19.05.2014
#
# fastNLO interface library with example scenarios for NLOJet++
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Targets
#
# fastNLO libtool libraries for use with NLOJet++
#
# Without fastjet
#
lib_LTLIBRARIES = libfnl5350eta0v22conee.la
#
# Require fastjet
#
if HAVE_FASTJET
   lib_LTLIBRARIES += libfnl5350eta0v22.la libInclusiveJets.la
# libfnl2342bv22.la
endif

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
# Common sources
AM_SOURCES =
# Individual sources
libfnl5350eta0v22conee_la_SOURCES = fnl5350eta0v22conee.cc $(AM_SOURCES)
libfnl5350eta0v22_la_SOURCES      = fnl5350eta0v22.cc      $(AM_SOURCES)
libInclusiveJets_la_SOURCES       = InclusiveJets.cc       $(AM_SOURCES)
# libfnl2342bv22_la_SOURCES         = fnl2342bv22.cc      $(AM_SOURCES)

#
# Preprocessor flags
#
# All targets
# AM_CPPFLAGS =

#
# C++ compiler flags
#
# top_srcdir for direct compilation before install; pkgincludedir for after install
#
# All targets
AM_CXXFLAGS = -g -O2 -Wall -I${NLOJETINCLUDEPATH} -I${FNLOTOOLKITINCLUDEPATH} -I$(srcdir)/../include -I$(top_builddir)/../include -I$(top_srcdir)/../include
# Debug flags
# AM_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall -I${NLOJETINCLUDEPATH} -I${FNLOTOOLKITINCLUDEPATH} -I../include
# With FastJet
if HAVE_FASTJETCONFIG
   AM_CXXFLAGS += `${FASTJETBINPATH}/fastjet-config --plugins --cxxflags`
else
if HAVE_FASTJET
   AM_CXXFLAGS += -I${FASTJETINCLUDEPATH}
endif
endif

#
# Linker flags
# top_builddir for distcheck; top_srcdir for direct linking before install; pkglibdir for after install
#
# All targets
AM_LDFLAGS = -L${top_builddir}/interface/jetalgos/.libs -L${top_srcdir}/interface/jetalgos -L$(pkglibdir) -ljetalgos -L${FNLOTOOLKITLIBPATH} -lfastnlotoolkit
if HAVE_FASTJETCONFIG
   AM_LDFLAGS += `${FASTJETBINPATH}/fastjet-config --libs --plugins`
else
if HAVE_FASTJET
   AM_LDFLAGS += -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone -lsiscone_spherical -lCDFConesPlugin -lPxConePlugin -lD0RunIIConePlugin -lNestedDefsPlugin -lTrackJetPlugin -lATLASConePlugin -lCMSIterativeConePlugin -lEECambridgePlugin -lJadePlugin
endif
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = .deps .obj .libs *.la \
	Makefile.in Makefile.gen Makefile
