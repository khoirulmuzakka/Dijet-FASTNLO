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
# The C++ fastnlo toolkit utility and example programs
#
if HAVE_LHAPDF
   bin_PROGRAMS = fnlo-tk-cppread fnlo-tk-example fnlo-tk-merge fnlo-tk-merge2 fnlo-tk-append fnlo-tk-modify fnlo-tk-cat
   bin_SCRIPTS  = fnlo-tk-config
# YODA features in fnlo-tk-statunc & fnlo-tk-yodaout are now conditionally compiled
#if HAVE_YODA
   bin_PROGRAMS += fnlo-tk-statunc
   bin_PROGRAMS += fnlo-tk-yodaout
#endif
if HAVE_ROOT
   bin_PROGRAMS += fnlo-tk-rootout
endif
endif

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
# Scripts
EXTRA_DIST = fnlo-tk-config.in
# Programs to be compiled
fnlo_tk_cppread_SOURCES = fnlo-tk-cppread.cc
fnlo_tk_example_SOURCES = fnlo-tk-example.cc
fnlo_tk_yodaout_SOURCES = fnlo-tk-yodaout.cc
fnlo_tk_statunc_SOURCES = fnlo-tk-statunc.cc
fnlo_tk_rootout_SOURCES = fnlo-tk-rootout.cc
fnlo_tk_merge_SOURCES   = fnlo-tk-merge.cc
fnlo_tk_merge2_SOURCES  = fnlo-tk-merge2.cc
fnlo_tk_append_SOURCES  = fnlo-tk-append.cc
fnlo_tk_modify_SOURCES  = fnlo-tk-modify.cc
fnlo_tk_cat_SOURCES   = fnlo-tk-cat.cc

#
# Preprocessor flags
#
my_include_dirs = -I${srcdir}/../fastnlotoolkit/include -I${top_builddir}/fastnlotoolkit/include -I${top_srcdir}/fastnlotoolkit/include

if HAVE_LHAPDF
   my_include_dirs += -I${LHAPDFINCLUDEPATH}
endif

fnlo_tk_cppread_CPPFLAGS = $(my_include_dirs) -I${includedir} $(LHAPDFCPPFLAGS)
fnlo_tk_example_CPPFLAGS = $(my_include_dirs) -I${includedir} $(ROOTCPPFLAGS)
if HAVE_ROOT
   fnlo_tk_example_CPPFLAGS += $(ROOTCPPFLAGS)
endif
fnlo_tk_yodaout_CPPFLAGS = $(my_include_dirs) -I${includedir} $(YODACPPFLAGS)
fnlo_tk_statunc_CPPFLAGS = $(my_include_dirs) -I${includedir} $(YODACPPFLAGS)
fnlo_tk_rootout_CPPFLAGS = $(my_include_dirs) -I${includedir} $(ROOTCPPFLAGS)
fnlo_tk_merge_CPPFLAGS   = $(my_include_dirs) -I${includedir}
fnlo_tk_merge2_CPPFLAGS  = $(my_include_dirs) -I${includedir}
fnlo_tk_append_CPPFLAGS  = $(my_include_dirs) -I${includedir}
fnlo_tk_modify_CPPFLAGS  = $(my_include_dirs) -I${includedir}
fnlo_tk_cat_CPPFLAGS     = $(my_include_dirs) -I${includedir}

#
# C++ compiler flags
#
# Recall: -std=c++11 is requested by default through configure.ac macro AX_CXX_COMPILE_STDCXX_11
#         and works with g++, clang etc.!
#         -std=gnu++11 is restricted to g++ only and not the best option
fnlo_tk_cppread_CXXFLAGS = -g -O2 -Wall
fnlo_tk_example_CXXFLAGS = -g -O2 -Wall -Wno-deprecated-declarations
if HAVE_ROOT
   fnlo_tk_example_CXXFLAGS += $(ROOTCXXFLAGS)
endif
fnlo_tk_yodaout_CXXFLAGS = -g -O2 -Wall
fnlo_tk_statunc_CXXFLAGS = -g -O2 -Wall
fnlo_tk_rootout_CXXFLAGS = -g -O2 -Wall $(ROOTCXXFLAGS) -Wdeprecated-declarations
fnlo_tk_merge_CXXFLAGS   = -g -O2 -Wall
fnlo_tk_merge2_CXXFLAGS  = -g -O2 -Wall
fnlo_tk_append_CXXFLAGS  = -g -O2 -Wall
fnlo_tk_modify_CXXFLAGS  = -g -O2 -Wall
fnlo_tk_cat_CXXFLAGS     = -g -O2 -Wall
# Debug flags
#fnlo_CXXFLAGS = -g -ftrapping-math -fsignaling-nans -O0 -Wall

#
# Linker flags
#
# Require lhapdf, use lhapdf-config if possible
my_lib_dirs = -L${srcdir}/../fastnlotoolkit -lfastnlotoolkit -L${top_builddir}/fastnlotoolkit -lfastnlotoolkit -L${top_srcdir}/fastnlotoolkit -lfastnlotoolkit
fnlo_tk_cppread_LDFLAGS = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit $(LHAPDFLDFLAGS)
fnlo_tk_example_LDFLAGS = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit
if HAVE_ROOT
   fnlo_tk_example_LDFLAGS += $(ROOTLDFLAGS)
endif
fnlo_tk_yodaout_LDFLAGS = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit $(YODALDFLAGS)
fnlo_tk_statunc_LDFLAGS = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit $(YODALDFLAGS)
fnlo_tk_rootout_LDFLAGS = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit $(ROOTLDFLAGS)
fnlo_tk_merge_LDFLAGS   = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit
fnlo_tk_merge2_LDFLAGS  = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit
fnlo_tk_append_LDFLAGS  = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit
fnlo_tk_modify_LDFLAGS  = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit
fnlo_tk_cat_LDFLAGS     = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
