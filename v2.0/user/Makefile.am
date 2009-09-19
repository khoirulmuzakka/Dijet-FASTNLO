AUTOMAKE_OPTIONS = foreign

# the main target, i.e. the executable to build
bin_PROGRAMS = example allunc

# all source files which go into the main target (and have to be distributed)
example_SOURCES = example.f fnx9999.f fnx9999.inc fx9999rw.f fnio.f fnset.f fn-alphas-demo.f fn-interface.f alphst.f pyalps.f
allunc_SOURCES = allunc.f statcode.f fnx9999.f fx9999rw.f fnio.f fnset.f \
	fn-interface.f \
	fn-alphas-demo.f alphst.f pyalps.f  \
	fnx9999.inc strings.inc

# FORTRAN compiler flags to be used for the main target
AM_FFLAGS = -O -g -fno-automatic -ffixed-line-length-132
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic -Wunused

# libraries to be linked to the main target
example_LDADD = -L$(LHAPDFLIBPATH) -lLHAPDF $(CERNLIBS)
allunc_LDADD = -L$(LHAPDFLIBPATH) -lLHAPDF $(CERNLIBS)
allunc_LDFLAGS = -lg2c

# path to LHAPDF pdfsets
%.f: %.F
	cpp -ansi -DLHAPDF_PATH=$(LHAPDFSETPATH) $*.F -o $*.f
# Trial to use in fx9999rw.F 
#	cpp -ansi -traditional-ccp -DLHAPDF_PATH=$(LHAPDFSETPATH) $*.F -o $*.f

# Clean files that are automatically produced by the preprocessor cpp
CLEANFILES = example.f fnio.f

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure missing install-sh
