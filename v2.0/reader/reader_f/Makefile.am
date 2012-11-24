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

#
# Targets
#
# The Fortran table reader example program
#
if HAVE_LHAPDFCONFIG
   bin_PROGRAMS = fnlo-fread
endif

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
fnlo_fread_SOURCES = main.f fnlo-fread.f fnx9999.f fnx9999.inc strings.inc \
	fx9999rw.f fnio.f fnset.f fn-interface.f alphas-grv.f alphas-grv-fnlo14.f 

#
# FORTRAN compiler flags
#
my_include_dirs = -I${srcdir} -I${top_builddir} -I${top_srcdir}
fnlo_fread_FFLAGS= -O -fno-automatic -ffixed-line-length-132 $(my_include_dirs)
# Debug flags
#fnlo_fread_FFLAGS= -O -fno-automatic -ffixed-line-length-72 $(my_include_dirs)
#fnlo_fread_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic $(my_include_dirs)
#fnlo_fread_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic -Wunused $(my_include_dirs)

#
# Linker flags
#
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

# Add dependencies on include files by hand
main.o:                   strings.inc
fnlo-fread.o: fnx9999.inc strings.inc
fx9999rw.o:   fnx9999.inc strings.inc
fnx9999.o:    fnx9999.inc strings.inc
fnset.o:      fnx9999.inc
alphas-grv.o:             strings.inc
