AUTOMAKE_OPTIONS = foreign

# the main target, i.e. the executable to build
bin_PROGRAMS = conv20

# all source files which go into the main target (and have to be distributed)
conv20_SOURCES = conv20.inc conv20.f  

# FORTRAN compiler flags to be used for the main target
AM_FFLAGS =  -O -fno-automatic -ffixed-line-length-132

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution 
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure missing install-sh
