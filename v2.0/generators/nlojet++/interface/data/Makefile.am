#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 21.06.2014
#
# NLOJet++ interface to fastNLO toolkit
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Data files to be distributed
#
fnlo_steering_files = fnl2342b_v22_fix.str fnl2342b_v22_flex.str
pkgdata_DATA = $(fnlo_steering_files)
EXTRA_DIST =   $(fnlo_steering_files)

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = Makefile.in
