#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 06.10.2013
#
# fastNLO_toolkit
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Targets
#
# The C++ fastnlo toolkit example program
#
if HAVE_LHAPDFCONFIG
   bin_PROGRAMS = fnlo
endif

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
fnlo_SOURCES = main.cc

#
# Preprocessor flags
#
my_include_dirs = -I${srcdir}/../fastnlotoolkit/include -I${top_builddir}/fastnlotoolkit/include -I${top_srcdir}/fastnlotoolkit/include
fnlo_CPPFLAGS = $(my_include_dirs) -I${LHAPDFINCLUDEPATH} -I${includedir}

#
# C++ compiler flags
#
fnlo_CXXFLAGS = -g -O0 -Wall
# Debug flags
#fnlo_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall

#
# Linker flags
#
# Require lhapdf, use lhapdf-config if possible
my_lib_dirs = -L${srcdir}/../fastnlotoolkit -lfastnlotoolkit -L${top_builddir}/fastnlotoolkit -lfastnlotoolkit -L${top_srcdir}/fastnlotoolkit -lfastnlotoolkit
if HAVE_LHAPDFCONFIG
   fnlo_LDFLAGS = $(my_lib_dirs) `${LHAPDFBINPATH}/lhapdf-config --ldflags` -L${libdir} -lfastnlotoolkit
else
   fnlo_LDFLAGS = $(my_lib_dirs) -L$(LHAPDFLIBPATH) -lLHAPDF -L${libdir} -lfastnlotoolkit
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
