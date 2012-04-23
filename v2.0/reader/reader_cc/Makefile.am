#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 10.11.2011
#
# fastNLO_reader: C++ version 
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main target: The C++ table reader
if HAVE_LHAPDFCONFIG
   bin_PROGRAMS = fnlo-cppread
endif

# All source files which go into the main target (and have to be distributed)
fnlo_cppread_SOURCES = main.cc fnlo-cppread.cc FastNLOReader.cc FastNLOReader.h FastNLOUser.h \
	FastNLOBlockB.cc FastNLOBlockB.h Alphas.cc Alphas.h \
	FastNLODiffReader.h FastNLODiffReader.cc FastNLODiffUser.h
#	diffpdf.f 
#	i_2006_fitb.f i_2006_fita.f h12006flux.f qcd_2006.f pion_stand_alone.f 

# Preprocessor flags
AM_CPPFLAGS = -I${LHAPDFINCLUDEPATH}

# Libraries to be linked to the main target
# Require lhapdf, use lhapdf-config if possible
if HAVE_LHAPDFCONFIG 
   fnlo_cppread_LDFLAGS = `${LHAPDFBINPATH}/lhapdf-config --ldflags`
else
   fnlo_cppread_LDFLAGS = -L$(LHAPDFLIBPATH) -lLHAPDF
endif


# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in
