AUTOMAKE_OPTIONS = gnu

# The main target: My allunc user code
if HAVE_LHAPDFCONFIG
#   if HAVE_CERNLIB
      bin_PROGRAMS = allunc
#   endif
endif

# All source files which go into the main target (and have to be distributed)
allunc_SOURCES = allunc.f statcode.f fnx9999.f fx9999rw.f fnio.f fnset.f \
	fn-interface.f \
	fn-alphas-demo.f alphst.f pyalps.f  \
	fnx9999.inc strings.inc

#
# FORTRAN compiler flags to be used for the main target
# top_srcdir for direct compilation before install; pkgincludedir for after install
#
#F77 = gfortran
AM_FFLAGS = -o allunc -O -g -ff2c -fno-automatic -ffixed-line-length-132
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic -Wunused

# C++ compiler flags to be used for the main target

# Libraries to be linked to the main target
# Require lhapdf, use lhapdf-config if possible
if HAVE_LHAPDFCONFIG 
   allunc_LDFLAGS = `lhapdf-config --ldflags`
else
   allunc_LDFLAGS = -L$(LHAPDFLIBPATH) -lLHAPDF
endif
if HAVE_CERNLIB
   allunc_LDADD = $(CERNLIBS) 
endif
# allunc_LDFLAGS = -lg2c

# Clean files that are automatically produced by the preprocessor cpp
#CLEANFILES = example.f fnio.f

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure missing install-sh
