# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# The main target, i.e. the executable to build
bin_PROGRAMS = fnlo-modify

# All source files which go into the main target
fnlo_modify_SOURCES = fnlo-modify.cc read_steer.cc 

# Libraries to be linked to the main target
# top_builddir for distcheck; top_srcdir for direct linking before install; pkglibdir for after install
LDADD = -L${top_builddir}/common/.libs -L${top_srcdir}/common -L$(pkglibdir) -lfastnlo

# C++ compiler flags to be used for the main target
AM_CXXFLAGS = 

# Preprocessor flags
# top_srcdir for direct compilation before install; pkgincludedir for after install
AM_CPPFLAGS = -I$(top_srcdir)/common -I${pkgincludedir} -I${NLOJETINCLUDEPATH}

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure depcomp missing install-sh
