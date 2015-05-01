#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 18.12.2014
#
# Include documentation into distribution tarball
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# Subdirectories
#if HAVE_PDFLATEX
#SUBDIRS = paper
#endif
# What about HAVE_LATEXDVIPDFM ? TBD

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = Makefile.in
