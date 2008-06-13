AUTOMAKE_OPTIONS = foreign

# the main target, i.e. the executable to build
bin_PROGRAMS = fnlo-merge

# all source files which go into the main target
fnlo_merge_SOURCES = main.cc entry.h

# libraries to be linked to the main target
LDADD = ../common/.libs/libfastnlo.a

# C++ compiler flags to be used for the main target
AM_CXXFLAGS = 

# Preprocessor flags
AM_CPPFLAGS = -I$(srcdir)/../common -I${HOME}/include/nlo++


# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure depcomp missing install-sh



