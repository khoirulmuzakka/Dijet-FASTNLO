AUTOMAKE_OPTIONS = foreign

# the main library with the fastNLO classes inside
pkglib_LTLIBRARIES = libfastnloAuthor.la

# all source files which go into the main target
libfastnloAuthor_la_SOURCES = fnloTable.cc fnloBlockA1.cc fnloBlockA2.cc fnloBlockB.cc fnloBlockBNlojet.cc \
  fnloconstants.h fnloTable.h fnloBlockA1.h  fnloBlockA2.h  fnloBlockB.h fnloBlockBNlojet.h fnloTableDict.cc

# C++ compiler flags to be used for the main target
AM_CXXFLAGS = 

# Preprocessor flags
AM_CPPFLAGS = -I${HOME}/include/nlo++ -I${ROOTSYS}/include

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure depcomp missing install-sh



