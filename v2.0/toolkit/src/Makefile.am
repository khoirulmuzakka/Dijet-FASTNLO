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
# The C++ fastnlo toolkit example programs
#
if HAVE_LHAPDF
   bin_PROGRAMS = fnlo-tk-cppread fnlo-tk-example fnlo-tk-merge fnlo-tk-append
   if HAVE_YODA
      bin_PROGRAMS += fnlo-tk-yodaout
   endif
endif

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
fnlo_tk_cppread_SOURCES = fnlo-tk-cppread.cc
fnlo_tk_example_SOURCES = fnlo-tk-example.cc
fnlo_tk_yodaout_SOURCES = fnlo-tk-yodaout.cc
fnlo_tk_merge_SOURCES   = fnlo-tk-merge.cc
fnlo_tk_append_SOURCES  = fnlo-tk-append.cc

#
# Preprocessor flags
#
my_include_dirs = -I${srcdir}/../fastnlotoolkit/include -I${top_builddir}/fastnlotoolkit/include -I${top_srcdir}/fastnlotoolkit/include
fnlo_tk_cppread_CPPFLAGS = $(my_include_dirs) -I${includedir}
fnlo_tk_example_CPPFLAGS = $(my_include_dirs) -I${includedir}
fnlo_tk_yodaout_CPPFLAGS = $(my_include_dirs) -I${includedir}
fnlo_tk_merge_CPPFLAGS   = $(my_include_dirs) -I${includedir}
fnlo_tk_append_CPPFLAGS  = $(my_include_dirs) -I${includedir}

#
# C++ compiler flags
#
fnlo_tk_cppread_CXXFLAGS = -g -O2 -Wall
fnlo_tk_merge_CXXFLAGS   = -g -O2 -Wall
fnlo_tk_append_CXXFLAGS  = -g -O2 -Wall
fnlo_tk_example_CXXFLAGS = -g -O2 -Wall
fnlo_tk_yodaout_CXXFLAGS = -g -O2 -Wall
# Debug flags
#fnlo_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall

#
# Linker flags
#
# Require lhapdf, use lhapdf-config if possible
my_lib_dirs = -L${srcdir}/../fastnlotoolkit -lfastnlotoolkit -L${top_builddir}/fastnlotoolkit -lfastnlotoolkit -L${top_srcdir}/fastnlotoolkit -lfastnlotoolkit
fnlo_tk_cppread_LDFLAGS = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit
fnlo_tk_example_LDFLAGS = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit
fnlo_tk_yodaout_LDFLAGS = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit
fnlo_tk_merge_LDFLAGS   = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit
fnlo_tk_append_LDFLAGS  = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
