#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 10.11.2011
#
# fastNLO_reader: Fortran version 
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main target: The Fortran table reader
if HAVE_LHAPDFCONFIG
   bin_PROGRAMS = fnlo-fread
endif

# All source files which go into the main target (and have to be distributed)
fnlo_fread_SOURCES = FReadFNLOTable.f fnx9999.f fnx9999.inc strings.inc \
	fx9999rw.f fnio.f fnset.f fn-interface.f alphas-grv.f alphas-grv-fnlo14.f 

# FORTRAN compiler flags to be used for the main target
# top_srcdir for direct compilation before install; pkgincludedir for after install
#F77 = gfortran
AM_FFLAGS= -O -fno-automatic -ffixed-line-length-132
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic -Wunused

# Libraries to be linked to the main target
# Require lhapdf, use lhapdf-config if possible
if HAVE_LHAPDFCONFIG 
   fnlo_fread_LDFLAGS = `${LHAPDFBINPATH}/lhapdf-config --ldflags`
else
   fnlo_fread_LDFLAGS = -L$(LHAPDFLIBPATH) -lLHAPDF
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in
