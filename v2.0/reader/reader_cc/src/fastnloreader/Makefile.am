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

# Subdirectories
SUBDIRS = include

#
# Targets
#
# The C++ table reader library
if HAVE_LHAPDFCONFIG
   lib_LTLIBRARIES = libfastnloreader.la
endif

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
libfastnloreader_la_SOURCES = FastNLOReader.cc FastNLODiffReader.cc FastNLOBlockB.cc Alphas.cc speaker.cc CRunDec.cc

#
# Preprocessor flags
#
my_include_dirs = -I${srcdir}/include -I${top_builddir}/include -I${top_srcdir}/include
libfastnloreader_la_CPPFLAGS = -I${includedir} $(my_include_dirs)

#
# C++ compiler flags
#
libfastnloreader_la_CXXFLAGS = -g -O0 -Wall
# Debug flags
#libfastnloreader_la_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in
