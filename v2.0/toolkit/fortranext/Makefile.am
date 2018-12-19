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
my_include_dirs = -I${srcdir}/../fastnlotoolkit/include -I${top_builddir}/fastnlotoolkit/include -I${top_srcdir}/fastnlotoolkit/include
fnlo_tk_fortex_CPPFLAGS = $(my_include_dirs) -I${includedir}

#
# Linker flags
#
# Require lhapdf, use lhapdf-config if possible
# The following construct is required to make things work on Mac
FRTLIB   = $(shell gfortran -print-file-name=libgfortran.a)
FRTLLIB  = -L$(subst /libgfortran.a, ,$(FRTLIB) ) -lgfortran
my_lib_dirs = -L${srcdir}/../fastnlotoolkit -lfastnlotoolkit -L${top_builddir}/fastnlotoolkit -lfastnlotoolkit -L${top_srcdir}/fastnlotoolkit -lfastnlotoolkit
fnlo_tk_fortex_LDFLAGS = $(FRTLLIB) $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit
#fnlo_tk_fortex_LDFLAGS =  -lgfortran $(my_lib_dirs) -L${pkglibdir} -lfastnlotoolkit

endif

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
