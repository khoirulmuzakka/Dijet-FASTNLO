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
   bin_PROGRAMS = fnlo-tk-cppread fnlo-tk-merge fastNLO
endif

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
fnlo_tk_cppread_SOURCES = fnlo-tk-cppread.cc
fnlo_tk_merge_SOURCES = fnlo-tk-merge.cc
fastNLO_SOURCES = fastNLO.cc

#
# Preprocessor flags
#
my_include_dirs = -I${srcdir}/../fastnlotoolkit/include -I${top_builddir}/fastnlotoolkit/include -I${top_srcdir}/fastnlotoolkit/include
fnlo_tk_cppread_CPPFLAGS = $(my_include_dirs) -I${LHAPDFINCLUDEPATH} -I${includedir}
fnlo_tk_merge_CPPFLAGS = $(my_include_dirs) -I${includedir}
fastNLO_CPPFLAGS = $(my_include_dirs) -I${LHAPDFINCLUDEPATH} -I${includedir}

#
# C++ compiler flags
#
fnlo_tk_cppread_CXXFLAGS = -g -O2 -Wall
fnlo_tk_merge_CXXFLAGS = -g -O2 -Wall
fastNLO_CXXFLAGS = -g -O2 -Wall
# Debug flags
#fnlo_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall

#
# Linker flags
#
# Require lhapdf, use lhapdf-config if possible
my_lib_dirs = -L${srcdir}/../fastnlotoolkit -lfastnlotoolkit -L${top_builddir}/fastnlotoolkit -lfastnlotoolkit -L${top_srcdir}/fastnlotoolkit -lfastnlotoolkit
fnlo_tk_merge_LDFLAGS = $(my_lib_dirs) -L${libdir} -lfastnlotoolkit
if HAVE_LHAPDFCONFIG
   fnlo_tk_cppread_LDFLAGS = $(my_lib_dirs) `${LHAPDFBINPATH}/lhapdf-config --ldflags` -L${libdir} -lfastnlotoolkit
   fastNLO_LDFLAGS = $(my_lib_dirs) `${LHAPDFBINPATH}/lhapdf-config --ldflags` -L${libdir} -lfastnlotoolkit
else
   fnlo_tk_cppread_LDFLAGS = $(my_lib_dirs) -L$(LHAPDFLIBPATH) -lLHAPDF -L${libdir} -lfastnlotoolkit
   fastNLO_LDFLAGS = $(my_lib_dirs) -L$(LHAPDFLIBPATH) -lLHAPDF -L${libdir} -lfastnlotoolkit
endif
if HAVE_QCDNUM
   fnlo_tk_cppread_LDFLAGS  += -lgfortran -L${QCDNUMLIBPATH} ${QCDNUMLIBS}
endif
if HAVE_HOPPET
   fnlo_tk_cppread_LDFLAGS  += -lgfortran ${HOPPETLDFLAGS}
   fnlo_tk_cppread_CPPFLAGS += ${HOPPETCPPFLAGS}
endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
