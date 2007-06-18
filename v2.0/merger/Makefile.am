AUTOMAKE_OPTIONS = foreign

# the main target, i.e. the executable to build
bin_PROGRAMS = fnlo-merge

# all source files which go into the main target
fnlo_merge_SOURCES = main.cc fnloTable.cc fnloBlockA1.cc fnloBlockA2.cc fnloBlockB.cc\
  fnloconstants.h fnloTable.h fnloBlockA1.h  fnloBlockA2.h  fnloBlockB.h

# C++ compiler flags to be used for the main target
AM_CXXFLAGS = 

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure depcomp missing install-sh



