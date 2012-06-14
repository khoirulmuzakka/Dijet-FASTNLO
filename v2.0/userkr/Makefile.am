#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 08.06.2012
#
# My table reader to derive all uncertainties: Fortran version 
#
# -----------------------------------------------------------------------------
# Provide local m4 macros in dir m4, copy also system-wide macros
ACLOCAL_AMFLAGS = -I m4 --install
#
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main target: My Fortran table reader
if HAVE_LHAPDFCONFIG
if HAVE_CERNLIB
   bin_PROGRAMS = allunc
endif
endif

# All source files which go into the main target (and have to be distributed)
allunc_SOURCES = main.f allunc.f fnx9999.f fn-interface.f alphas-grv.f alphas-grv-fnlo14.f alphst.f \
	pyalps.f fx9999rw.f fnio.f uncert.f fnset.f fnx9999.inc strings.inc uncert.inc

# FORTRAN compiler flags to be used for the main target
# top_srcdir for direct compilation before install; pkgincludedir for after install
#F77 = gfortran
AM_FFLAGS= -O -fno-automatic -ffixed-line-length-132
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic -Wunused

# Libraries to be linked to the main target
# Require cernlib and lhapdf, use lhapdf-config if possible
allunc_LDFLAGS = -L$(CERNLIBLIBPATH) -lmathlib -lpacklib -lkernlib
if HAVE_LHAPDFCONFIG 
   allunc_LDFLAGS += `${LHAPDFBINPATH}/lhapdf-config --ldflags`
else
   allunc_LDFLAGS += -L$(LHAPDFLIBPATH) -lLHAPDF
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = INSTALL config.h.in config.guess config.sub configure Makefile.in \
	depcomp install-sh ltmain.sh missing aclocal.m4 SVNREV \
	m4/libtool.m4 m4/ltoptions.m4 m4/ltsugar.m4 m4/ltversion.m4 m4/lt~obsolete.m4

# Add dependencies on include files by hand
main.o:                   strings.inc
allunc.o:     fnx9999.inc strings.inc uncert.inc
uncert.o:     fnx9999.inc             uncert.inc
fx9999rw.o:   fnx9999.inc strings.inc
fnx9999.o:    fnx9999.inc strings.inc
fnset.o:      fnx9999.inc
alphas-grv.o:             strings.inc
