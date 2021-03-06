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
if HAVE_LHAPDF
if HAVE_CERNLIB
bin_PROGRAMS = allunc
endif
endif

# All source files which go into the main target (and have to be distributed)
allunc_SOURCES = main.f allunc.f fnx9999.f fn-interface.f alphas-grv.f alphas-grv-fnlo14.f alphst.f \
	pyalps.f fx9999rw.f fnio.f uncert.f uncert2.f fnset.f fnx9999.inc strings.inc uncert.inc v14unc.inc \
	v22tools.f

# FORTRAN compiler flags to be used for the main target
# top_srcdir for direct compilation before install; pkgincludedir for after install
#F77 = gfortran
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-132
AM_FFLAGS= -O -fno-automatic -ffixed-line-length-132 -mcmodel=medium
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic
#AM_FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic -Wunused

# Libraries to be linked to the main target
# Require lhapdf and cernlib, use lhapdf-config and cernlib command if possible
allunc_LDFLAGS =
if HAVE_LHAPDFCONFIG
allunc_LDFLAGS += -Wl,--static `${LHAPDFBINPATH}/lhapdf-config --ldflags`
else
#allunc_LDFLAGS  += -L$(LHAPDFLIBPATH) -Wl,-static -lLHAPDF
allunc_LDFLAGS  += -L/home/rabbertz/tmp/2005/lib -Wl,-static -lmathlib -lkernlib -lpawlib -lpacklib -Wl,-Bdynamic
#allunc_LDFLAGS  += -L/usr/lib/x86_64-linux-gnu -Wl,-static -lmathlib -lkernlib -lpawlib -lpacklib -Wl,-Bdynamic
allunc_LDFLAGS  += -L$(LHAPDFLIBPATH)  -lLHAPDF
#allunc_LDFLAGS  += -L$(CERNLIBLIBPATH) -Wl,-static -lmathlib -lkernlib -lpawlib -lpacklib -Wl,-Bdynamic
endif
# To avoid 32-64-bit compatibility problems 64-bit systems require static linking of at least kernlib
# (Tried to use "--as-needed" option instead of "--static" but didn't work as intended.)
#if HAVE_64BIT
#allunc_LDFLAGS += $(CERNLIBLIBPATH)/libmathlib.a $(CERNLIBLIBPATH)/libpacklib.a $(CERNLIBLIBPATH)/libkernlib.a
#else
if HAVE_CERNLIBCONFIG
allunc_LDFLAGS += -L$(CERNLIBLIBPATH) -Wl,--static `${CERNLIBCMD} mathlib packlib kernlib`
else
allunc_LDFLAGS += -L$(CERNLIBLIBPATH) -Wl,-static -lmathlib -lpacklib -lkernlib -Wl,-Bdynamic
#allunc_LDFLAGS += -L$(CERNLIBLIBPATH) -Wl,-static libmathlib.a libpacklib.a libkernlib.a -Wl,-Bdynamic
endif
#endif


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
allunc.o:     fnx9999.inc strings.inc uncert.inc v14unc.inc
uncert.o:     fnx9999.inc             uncert.inc v14unc.inc
fx9999rw.o:   fnx9999.inc strings.inc
fnx9999.o:    fnx9999.inc strings.inc
fnset.o:      fnx9999.inc
alphas-grv.o:             strings.inc
