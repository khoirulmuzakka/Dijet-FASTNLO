# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main target: My default user code
if HAVE_LHAPDFCONFIG
   bin_PROGRAMS = myexample myexample-ref
endif

# All source files which go into the main target (and have to be distributed)
myexample_SOURCES = myexample.f fnx9999.f fnx9999.inc strings.inc \
	fnset.f fn-alphas-demo.f fn-interface.f fx9999rw.f fnio.f
myexample_ref_SOURCES = myexample-ref.f fnx9999.f fnx9999.inc strings.inc \
	fnset.f fn-alphas-demo.f fn-interface.f fx9999rw.f fnio.f

#
# FORTRAN compiler flags to be used for the main target
# top_srcdir for direct compilation before install; pkgincludedir for after install
#
#F77 = gfortran
AM_FFLAGS = -o myexample -O -fno-automatic -ffixed-line-length-132
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic -Wunused

# C++ compiler flags to be used for the main target

# Libraries to be linked to the main target
# Require lhapdf, use lhapdf-config if possible
if HAVE_LHAPDFCONFIG 
   myexample_LDFLAGS = `${LHAPDFBINPATH}/lhapdf-config --ldflags`
   myexample_ref_LDFLAGS = `${LHAPDFBINPATH}/lhapdf-config --ldflags`
else
   myexample_LDFLAGS = -L$(LHAPDFLIBPATH) -lLHAPDF
   myexample_ref_LDFLAGS = -L$(LHAPDFLIBPATH) -lLHAPDF
endif

# cpp trick for path to LHAPDF pdfsets
#%.f: %.F
#	cpp -ansi -DLHAPDF_PATH=$(LHAPDFSETPATH) $*.F -o $*.f
# Trial to use in fx9999rw.F 
#	cpp -ansi -traditional-ccp -DLHAPDF_PATH=$(LHAPDFSETPATH) $*.F -o $*.f

# Clean files that are automatically produced by the preprocessor cpp
#CLEANFILES = example.f fnio.f

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure missing install-sh
