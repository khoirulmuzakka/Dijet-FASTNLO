AUTOMAKE_OPTIONS = foreign

# the main target, i.e. the executable to build
bin_PROGRAMS = example

# all source files which go into the main target (and have to be distributed)
example_SOURCES = example.f fnx9999.f fnx9999.inc fx9999rw.f fnio.f fnset.f fn-alphas-demo.f fn-interface.f  

# FORTRAN compiler flags to be used for the main target
AM_FFLAGS = -O -fno-automatic -ffixed-line-length-132

# libraries to be linked to the main target
LDADD = -L$(LHAPDFLIBPATH) -lLHAPDF $(CERNLIBS)

# path to LHAPDF pdfsets
%.f: %.F
	cpp -ansi -DLHAPDF_PATH=$(LHAPDFSETPATH) $*.F -o $*.f

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = example.f \
	Makefile.in aclocal.m4 configure missing install-sh
