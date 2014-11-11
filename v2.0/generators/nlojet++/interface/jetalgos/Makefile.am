#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 19.05.2014
#
# fastNLO interface library with jet algos for NLOJet++
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Targets
#
# fastNLO interface library with the jet algos for NLOJet++
pkglib_LTLIBRARIES = libjetalgos.la

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
# Common sources
AM_SOURCES =
# Individual sources
libjetalgos_la_SOURCES = cone-e.cc cone-et.cc kt-e.cc $(AM_SOURCES)
# With FastJet
if HAVE_FASTJET
   libjetalgos_la_SOURCES += fj-jets.cc
endif

#
# Preprocessor flags
#
# All targets
# AM_CPPFLAGS =

#
# C++ compiler flags
#
# top_srcdir for direct compilation before install; pkgincludedir for after install
# All targets
AM_CXXFLAGS = -g -O2 -Wall -I${NLOJETINCLUDEPATH} `${FNLOTOOLKITBINPATH}/fnlo-tk-config --cppflags` -I$(srcdir)/../include -I$(top_builddir)/../include -I$(top_srcdir)/../include
# Debug flags
# AM_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall -I${NLOJETINCLUDEPATH} `${FNLOTOOLKITBINPATH}/fnlo-tk-config --cppflags` -I$(srcdir)/../include -I$(top_builddir)/../include -I$(top_srcdir)/../include
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
AM_LDFLAGS =
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
