# Provide local m4 macros in dir m4, copy also system-wide macros
#ACLOCAL_AMFLAGS = -I m4 --install
# Provide local m4 macros in dir m4, don't copy system-wide macros
ACLOCAL_AMFLAGS = -I m4
#
AUTOMAKE_OPTIONS = gnu
#
# Add subversion revision number, recipe from Autotools by J. Calcote 
EXTRA_DIST = m4 SVNREV
#distdir = $(PACKAGE)-$(VERSION).$(SVNREV)
#
SUBDIRS = common merger tableconvert author/ep author/hadron user
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
