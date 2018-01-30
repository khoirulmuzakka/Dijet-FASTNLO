#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 15.05.2014
#
# NLOJet++ interface to fastNLO toolkit
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# Subdirectories
SUBDIRS = interface

# Provide local m4 macros in dir m4, copy also system-wide macros
ACLOCAL_AMFLAGS = -I m4 --install

# Add git equivalent of subversion revision number (original svn recipe from Autotools by J. Calcote).
EXTRA_DIST = m4 GITREV
# This gives e.g. fastnlo-2.0.0-643M for svn revision 643 locally modified
distdir = $(PACKAGE)-$(VERSION)-$(GITREV)
# Don't put .git subdirs into distribution package
dist-hook:
	rm -rf `find $(distdir) -name .git`

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before git commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = INSTALL config.h.in config.guess config.sub configure Makefile.in \
	depcomp install-sh ltmain.sh missing compile \
	aclocal.m4 m4/libtool.m4 m4/ltoptions.m4 m4/ltsugar.m4 m4/ltversion.m4 m4/lt~obsolete.m4 \
	GITREV test-driver
