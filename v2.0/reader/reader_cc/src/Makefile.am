#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 22.11.2012
#
# fastNLO_reader: C++ version
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Targets
#
# The C++ table reader example program
#
if HAVE_LHAPDFCONFIG
   bin_PROGRAMS = fnlo-cppread
endif

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
fnlo_cppread_SOURCES = main.cc

#
# Preprocessor flags
#
my_include_dirs = -I${srcdir}/../fastnloreader/include -I${top_builddir}/fastnloreader/include -I${top_srcdir}/fastnloreader/include
fnlo_cppread_CPPFLAGS = $(my_include_dirs) -I${LHAPDFINCLUDEPATH} -I${includedir}

#
# C++ compiler flags
#
fnlo_cppread_CXXFLAGS = -g -O0 -Wall
# Debug flags
#fnlo_cppread_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall

#
# Linker flags
#
# Require lhapdf, use lhapdf-config if possible
my_lib_dirs = -L${srcdir}/../fastnloreader -lfastnloreader -L${top_builddir}/reader_cc/fastnloreader -lfastnloreader -L${top_srcdir}/reader_cc/fastnloreader -lfastnloreader
if HAVE_LHAPDFCONFIG
fnlo_cppread_LDFLAGS = $(my_lib_dirs) `${LHAPDFBINPATH}/lhapdf-config --ldflags` -L${libdir} -lfastnloreader
else
   fnlo_cppread_LDFLAGS = $(my_lib_dirs) -L$(LHAPDFLIBPATH) -lLHAPDF -L${libdir} -lfastnloreader
endif
if HAVE_QCDNUM
   fnlo_cppread_LDFLAGS += -L${QCDNUMLIBPATH} ${QCDNUMLIBS}
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
