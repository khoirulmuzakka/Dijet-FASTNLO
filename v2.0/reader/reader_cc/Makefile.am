# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main target: The C++ table reader
if HAVE_LHAPDFCONFIG
   bin_PROGRAMS = ReadFNLOTable
endif

# All source files which go into the main target (and have to be distributed)
ReadFNLOTable_SOURCES = ReadFNLOTable.cc FastNLOReader.cc FastNLOReader.h \
	FastNLOBlockB.cc FastNLOBlockB.h qcdnumdummy.f Alphas.cc Alphas.h

# Preprocessor flags
AM_CPPFLAGS = -I${LHAPDFINCLUDEPATH}

# Libraries to be linked to the main target
# Require lhapdf, use lhapdf-config if possible
if HAVE_LHAPDFCONFIG 
   ReadFNLOTable_LDFLAGS = `${LHAPDFBINPATH}/lhapdf-config --ldflags`
else
   ReadFNLOTable_LDFLAGS = -L$(LHAPDFLIBPATH) -lLHAPDF
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in
