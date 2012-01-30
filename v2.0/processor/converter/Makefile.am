#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 26.11.2011
#
# Create separate distribution and installation package for fastNLO_converter
#
# -----------------------------------------------------------------------------
# Provide local m4 macros in dir m4, copy also system-wide macros
ACLOCAL_AMFLAGS = -I m4 --install
#
# Use GNU standards
AUTOMAKE_OPTIONS = gnu
#
EXTRA_DIST =
if HAVE_SUBPROJECT
# Add subversion revision number, recipe from Autotools by J. Calcote 
# Add source code of NLOJet++ 4.1.3 to be changed for use with fastNLO
   EXTRA_DIST += m4 SVNREV
# This gives e.g. fastnlo-2.0.0-643M for svn revision 643 locally modified
   distdir = $(PACKAGE)-$(VERSION)-$(SVNREV)
# Don't put .svn subdirs into distribution package
dist-hook:
	rm -rf `find $(distdir) -name .svn`
endif
#
# The main target: The Fortran table converter
bin_PROGRAMS = fnlo-convert

# All source files which go into the main target (and have to be distributed)
fnlo_convert_SOURCES = fnlo-convert.f fnlo-convert.inc strings.inc

# FORTRAN compiler flags to be used for the main target
AM_FFLAGS = -O -fno-automatic -ffixed-line-length-132

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"! 
MAINTAINERCLEANFILES = INSTALL config.h.in config.guess config.sub configure Makefile.in \
	depcomp install-sh ltmain.sh missing SVNREV aclocal.m4 \
	m4/libtool.m4 m4/ltoptions.m4 m4/ltsugar.m4 m4/ltversion.m4 m4/lt~obsolete.m4

# Additional scripts
bin_SCRIPTS = tablefix.pl
EXTRA_DIST  += tablefix.pl
