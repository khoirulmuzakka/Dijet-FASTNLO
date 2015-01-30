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
# Nothing to do here, just go to subdirectories
SUBDIRS = data include tools jetalgos hadron
#
# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = Makefile.in
