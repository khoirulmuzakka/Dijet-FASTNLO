#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 06.10.2013
#
# fastNLO_toolkit
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# Subdirectories
SUBDIRS = include

#
# Targets
#
# The C++ fastnlo toolkit library
if HAVE_LHAPDFCONFIG
   lib_LTLIBRARIES = libfastnlotoolkit.la
endif

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
libfastnlotoolkit_la_SOURCES = speaker.cc fastNLOTools.cc fastNLOInterpolBase.cc fastNLOInterpolCatmullRom.cc fastNLOInterpolLagrange.cc fastNLOInterpolLinear.cc fastNLOInterpolOneNode.cc fastNLOBase.cc fastNLOTable.cc fastNLOCoefficients.cc fastNLOCreate.cc read_steer.cc fastNLOCoeffBase.cc fastNLOCoeffData.cc fastNLOCoeffMult.cc fastNLOCoeffAddBase.cc fastNLOCoeffAddFix.cc fastNLOCoeffAddFlex.cc fastNLOReader.cc fastNLOPDFLinearCombinations.cc Alphas.cc CRunDec.cc

#
# Preprocessor flags
#
my_include_dirs = -I${srcdir}/include -I${top_builddir}/include -I${top_srcdir}/include
libfastnlotoolkit_la_CPPFLAGS = $(my_include_dirs) -I${includedir}

#
# C++ compiler flags
#
libfastnlotoolkit_la_CXXFLAGS =
#-g -O2
# -Wall
# Debug flags
#libfastnlotoolkit_la_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
