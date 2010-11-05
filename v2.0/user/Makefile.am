# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main target: Markus' default user code
bin_PROGRAMS = example example-ref

# All source files which go into the main target (and have to be distributed)
example_SOURCES = example.f fnx9999.f fnset.f fn-alphas-demo.f fn-interface.f fx9999rw.f fnio.f \
	fnx9999.inc strings.inc
example_ref_SOURCES = example-ref.f fnx9999.f fnset.f fn-alphas-demo.f fn-interface.f fx9999rw.f fnio.f \
	fnx9999.inc strings.inc

# FORTRAN compiler flags to be used for the main target
# top_srcdir for direct compilation before install; pkgincludedir for after install
AM_FFLAGS = -o example -O -fno-automatic -ffixed-line-length-132
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic -Wunused

# C++ compiler flags to be used for the main target

# Libraries to be linked to the main target
example_LDADD = -L$(LHAPDFLIBPATH) -lLHAPDF $(CERNLIBS)
example_ref_LDADD = -L$(LHAPDFLIBPATH) -lLHAPDF $(CERNLIBS)
#allunc_LDFLAGS = -lg2c

# path to LHAPDF pdfsets
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
