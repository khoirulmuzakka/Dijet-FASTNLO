#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 22.11.2012
#
# fastNLO_reader: C++ version 
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Targets
#
# The C++ table reader example program
#
if HAVE_LHAPDFCONFIG
   bin_PROGRAMS = fnlo-cppread
endif

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
fnlo_cppread_SOURCES = main.cc fnlo-cppread.cc

#
# Preprocessor flags
#
fnlo_cppread_CPPFLAGS = -I${LHAPDFINCLUDEPATH}

#
# C++ compiler flags
#
fnlo_cppread_CXXFLAGS = -g -O0 -Wall -I${includedir} -I${builddir}/../include
# Debug flags
#fnlo_cppread_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall -I${includedir} -I${builddir}/../include

#
# Linker flags
#
# Require lhapdf, use lhapdf-config if possible
if HAVE_LHAPDFCONFIG 
   fnlo_cppread_LDFLAGS = `${LHAPDFBINPATH}/lhapdf-config --ldflags` -L${libdir} -lfastnloreader -L${builddir}/../lib -lfastnloreader
else
   fnlo_cppread_LDFLAGS = -L$(LHAPDFLIBPATH) -lLHAPDF -L${libdir} -lfastnloreader -L${builddir}/../lib -lfastnloreader
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in
