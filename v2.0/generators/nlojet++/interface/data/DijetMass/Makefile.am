#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 17.01.2015
#
# NLOJet++ interface to fastNLO toolkit
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Data files to be distributed
#
fnlo_steering_files  = fnl2412c_I895742_v23_fix.str fnl2412c_I895742_v23_flex.str
fnlo_steering_files += fnl2452_I1082936_v23_fix_ptmax.str fnl2452_I1082936_v23_fix_ptmaxexpystar.str fnl2452_I1082936_v23_fix_pt12ave.str fnl2452_I1082936_v23_flex.str
pkgdata_DATA = $(fnlo_steering_files)
EXTRA_DIST =   $(fnlo_steering_files)

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = Makefile.in
