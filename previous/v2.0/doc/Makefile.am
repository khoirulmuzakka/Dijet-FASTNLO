#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 28.11.2011
#
# Include documentation into distribution tarball
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu
# Only include arXiv_1109.1310 subdirectory for now
SUBDIRS = arXiv_1109.1310
#
# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = Makefile.in
