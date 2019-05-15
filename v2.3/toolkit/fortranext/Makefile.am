#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 08.11.2016
#
# fastNLO_toolkit
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# Optional extension
if ENABLE_FORTRANEXT

#
# Targets
#
# The fastnlo toolkit example for use with Fortran
#
   bin_PROGRAMS = fnlo-tk-fortex

#
# Sources
#
# All source files which go into the targets (and have to be distributed)
# Programs to be compiled
   fnlo_tk_fortex_SOURCES = fastNLOFI_cc.cc fastNLOFI_f.f

#
# Preprocessor flags
#
# All targets
   AM_CPPFLAGS =
   my_include_dirs = -I${srcdir}/../fastnlotoolkit/include -I${top_builddir}/fastnlotoolkit/include -I${top_srcdir}/fastnlotoolkit/include

# Individual targets
   fnlo_tk_fortex_CPPFLAGS = $(my_include_dirs) -I${includedir} $(AM_CPPFLAGS)

#
# Fortran compiler flags
#

#
# Linker flags
#
# All targets
# FRTLDFLAGS is determined and set in configure.ac
# This is required to make things work on Mac
   AM_LDFLAGS = $(FRTLDFLAGS)
   my_lib_dirs = -L${srcdir}/../fastnlotoolkit -lfastnlotoolkit -L${top_builddir}/fastnlotoolkit -lfastnlotoolkit -L${top_srcdir}/fastnlotoolkit -lfastnlotoolkit
# Individual targets
   fnlo_tk_fortex_LDFLAGS = $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit $(AM_LDFLAGS)

endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
