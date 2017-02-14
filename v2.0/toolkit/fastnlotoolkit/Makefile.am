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
#
if HAVE_LHAPDFCONFIG
   lib_LTLIBRARIES = libfastnlotoolkit.la
endif

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
# Common sources
AM_SOURCES =
# Individual sources
libfastnlotoolkit_la_SOURCES = speaker.cc fastNLOTools.cc fastNLOInterpolBase.cc fastNLOInterpolCatmullRom.cc fastNLOInterpolLagrange.cc fastNLOInterpolLinear.cc fastNLOInterpolOneNode.cc fastNLOTable.cc fastNLOCoefficients.cc fastNLOCreate.cc read_steer.cc fastNLOCoeffBase.cc fastNLOCoeffData.cc fastNLOCoeffMult.cc fastNLOCoeffAddBase.cc fastNLOCoeffAddFix.cc fastNLOCoeffAddFlex.cc fastNLOReader.cc fastNLOPDFLinearCombinations.cc fastNLODiffReader.cc fastNLODiffUser.cc Alphas.cc CRunDec.cc $(AM_SOURCES)

#
# Preprocessor flags
#
# All targets, -fext-numeric-literals needed for conflict with Boost
AM_CPPFLAGS = -std=c++11 -lstdc++ -fext-numeric-literals
# Individual targets
libfastnlotoolkit_la_CPPFLAGS = $(AM_CPPFLAGS)
#
# Use precompiler-defined char constant FNLO_SOMETHING in fastNLOConstants.h.in
# instead of precompiler idfdef's via -DHAVESOMETHING in the code
# With HOPPET
#if HAVE_HOPPET
#   libfastnlotoolkit_la_CPPFLAGS += -DHAVEHOPPET
#endif
# With QCDNUM
#if HAVE_QCDNUM
#   libfastnlotoolkit_la_CPPFLAGS += -DHAVEQCDNUM
#endif

#
# C++ compiler flags
#
# top_srcdir for direct compilation before install; pkgincludedir for after install
#
# All targets
AM_CXXFLAGS = -g -O2 -Wall  -Wdeprecated-declarations -I$(top_srcdir)/include -I$(srcdir)/include -I$(top_builddir)/include
# Debug flags
# AM_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall  -I$(top_srcdir)/include -I$(srcdir)/include -I$(top_builddir)/include
# Individual targets
libfastnlotoolkit_la_CXXFLAGS = $(AM_CXXFLAGS)
# Require LHAPDF
if HAVE_LHAPDF
   libfastnlotoolkit_la_SOURCES  += fastNLOLHAPDF.cc fastNLOAlphas.cc fastNLOCRunDec.cc
   libfastnlotoolkit_la_CXXFLAGS += -I${LHAPDFINCLUDEPATH}
endif
# With HOPPET
if HAVE_HOPPET
   libfastnlotoolkit_la_SOURCES  += fastNLOHoppet.cc fastNLOHoppetAs.cc HoppetInterface.cc
   libfastnlotoolkit_la_CXXFLAGS += ${HOPPETCXXFLAGS}
endif
# With QCDNUM
if HAVE_QCDNUM
   libfastnlotoolkit_la_SOURCES  += fastNLOQCDNUMAS.cc
endif

#
# Linker flags
# top_builddir for distcheck; top_srcdir for direct linking before install; pkglibdir for after install
#
# All targets
AM_LDFLAGS =
# Individual targets
libfastnlotoolkit_la_LDFLAGS = $(AM_LDFLAGS)
# Require LHAPDF
if HAVE_LHAPDF
if HAVE_LHAPDFCONFIG
   libfastnlotoolkit_la_LDFLAGS  += `${LHAPDFBINPATH}/lhapdf-config --ldflags`
else
   libfastnlotoolkit_la_LDFLAGS  += -L$(LHAPDFLIBPATH) -lLHAPDF
endif
endif
# With HOPPET
if HAVE_HOPPET
   libfastnlotoolkit_la_LDFLAGS  += -lgfortran ${HOPPETLDFLAGS}
endif
# With QCDNUM
if HAVE_QCDNUM
   libfastnlotoolkit_la_LDFLAGS  += -lgfortran ${QCDNUMLDFLAGS}
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
