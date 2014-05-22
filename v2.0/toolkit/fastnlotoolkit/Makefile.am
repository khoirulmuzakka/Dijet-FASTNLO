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
libfastnlotoolkit_la_SOURCES = speaker.cc fastNLOTools.cc fastNLOInterpolBase.cc fastNLOInterpolCatmullRom.cc fastNLOInterpolLagrange.cc fastNLOInterpolLinear.cc fastNLOInterpolOneNode.cc fastNLOBase.cc fastNLOTable.cc fastNLOCoefficients.cc fastNLOCreate.cc read_steer.cc fastNLOCoeffBase.cc fastNLOCoeffData.cc fastNLOCoeffMult.cc fastNLOCoeffAddBase.cc fastNLOCoeffAddFix.cc fastNLOCoeffAddFlex.cc fastNLOReader.cc fastNLOPDFLinearCombinations.cc fastNLODiffReader.cc fastNLODiffUser.cc Alphas.cc CRunDec.cc $(AM_SOURCES)

#
# Preprocessor flags
#
# All targets
# AM_CPPFLAGS =

#
# C++ compiler and linker flags
#
# top_srcdir for direct compilation before install; pkgincludedir for after install
#
# All targets
# AM_CXXFLAGS = -g -O2 -Wall
# Debug flags
# AM_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall
# Individual targets
libfastnlotoolkit_la_CXXFLAGS = -g -O2 -Wall
# Debug flags
#libfastnlotoolkit_la_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall
# Require LHAPDF
if HAVE_LHAPDF
   libfastnlotoolkit_la_SOURCES  += fastNLOLHAPDF.cc fastNLOAlphas.cc fastNLOCRunDec.cc
   libfastnlotoolkit_la_CXXFLAGS += -I${LHAPDFINCLUDEPATH}
if HAVE_LHAPDFCONFIG
   libfastnlotoolkit_la_LDFLAGS += `${LHAPDFBINPATH}/lhapdf-config --ldflags`
else
   libfastnlotoolkit_la_LDFLAGS += -L$(LHAPDFLIBPATH) -lLHAPDF
endif

# With QCDNUM
if HAVE_QCDNUM
   libfastnlotoolkit_la_CPPFLAGS += -DHAVEQCDNUM
   libfastnlotoolkit_la_SOURCES  += fastNLOQCDNUMAS.cc
   libfastnlotoolkit_la_LDFLAGS  += -lgfortran -L${QCDNUMLIBPATH} ${QCDNUMLIBS}
endif

# With HOPPET
if HAVE_HOPPET
   libfastnlotoolkit_la_CPPFLAGS += -DHAVEHOPPET
   libfastnlotoolkit_la_SOURCES  += fastNLOHoppet.cc fastNLOHoppetAs.cc
   libfastnlotoolkit_la_LDFLAGS  += -lgfortran ${HOPPETLDFLAGS}
   libfastnlotoolkit_la_CPPFLAGS += ${HOPPETCPPFLAGS}
endif
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
