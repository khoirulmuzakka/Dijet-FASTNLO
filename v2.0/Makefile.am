#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 06.02.2009
#
# Create distribution and installation package for fastNLO v2 
#
# -----------------------------------------------------------------------------
# Provide local m4 macros in dir m4, copy also system-wide macros
#ACLOCAL_AMFLAGS = -I m4 --install
# Provide local m4 macros in dir m4, don't copy system-wide macros
ACLOCAL_AMFLAGS = -I m4
#
# Use GNU standards
AUTOMAKE_OPTIONS = gnu
#
# Add subversion revision number, recipe from Autotools by J. Calcote 
# Add source code of NLOJet++ 4.1.3 to be changed for use with fastNLO
EXTRA_DIST = m4 SVNREV nlojet++-4.1.3
# This gives e.g. fastnlo-2.0.0-643M for svn revision 643 locally modified
distdir = $(PACKAGE)-$(VERSION)-$(SVNREV)
# Don't put .svn subdirs into distribution package
dist-hook:
	rm -rf `find $(distdir) -name .svn`
#
SUBDIRS = common creator processor reader
#userkr
#
# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"! 
MAINTAINERCLEANFILES = INSTALL config.h.in config.guess config.sub configure Makefile.in \
	depcomp install-sh ltmain.sh missing \
	aclocal.m4 m4/libtool.m4 m4/ltoptions.m4 m4/ltsugar.m4 m4/ltversion.m4 m4/lt~obsolete.m4 \
	SVNREV
