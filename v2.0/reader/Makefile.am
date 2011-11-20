#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 10.11.2011
#
# Create separate distribution and installation package for fastNLO_reader 
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
EXTRA_DIST = m4 SVNREV
# This gives e.g. fastnlo-2.0.0-643M for svn revision 643 locally modified
distdir = $(PACKAGE)-$(VERSION)-$(SVNREV)
# Don't put .svn subdirs into distribution package
dist-hook:
	rm -rf `find $(distdir) -name .svn`
#
SUBDIRS = reader_cc reader_f
