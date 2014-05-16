#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 13.12.2011
#
# Common library
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# Subdirectories
SUBDIRS = include

#
# Targets
#
# The fastNLO library with the jet algos and interfaces inside
pkglib_LTLIBRARIES = libfastnlojetalgos.la

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
libfastnlojetalgos_la_SOURCES = cone-e.cc cone-et.cc kt-e.cc

#
# Preprocessor flags
#
my_include_dirs = -I${srcdir}/include -I${top_builddir}/include -I${top_srcdir}/include
libfastnlojetalgos_la_CPPFLAGS = $(my_include_dirs) -I${includedir} -I${NLOJETINCLUDEPATH}
libfastnlojetalgos_la_LDFLAGS =

#
# C++ compiler flags
#
libfastnlojetalgos_la_CXXFLAGS =
#-g -O2
# -Wall
# Debug flags
#libfastnlojetalgos_la_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall

#Compile only if FastJet is available
if HAVE_FASTJET
  libfastnlojetalgos_la_CPPFLAGS += -I${FASTJETINCLUDEPATH}
  libfastnlojetalgos_la_SOURCES  += fj-ak.cc fj-kt.cc fj-sc.cc fj-d0r2.cc
endif
if HAVE_FASTJETCONFIG
   libfastnlojetalgos_la_LDFLAGS += `${FASTJETBINPATH}/fastjet-config --libs`
else
   libfastnlojetalgos_la_LDFLAGS += -L${FASTJETLIBPATH} -lfastjet -lSISConePlugin -lsiscone -lsiscone_spherical -lCDFConesPlugin -lPxConePlugin -lD0RunIIConePlugin -lNestedDefsPlugin -lTrackJetPlugin -lATLASConePlugin -lCMSIterativeConePlugin -lEECambridgePlugin -lJadePlugin
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = configure Makefile.in aclocal.m4 depcomp missing install-sh
