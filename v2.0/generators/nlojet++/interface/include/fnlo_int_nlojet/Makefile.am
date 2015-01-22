#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 15.05.2014
#
# Header files for fastNLO interface to NLOJet++
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Include files to be distributed
#
# pkginclude_HEADERS would install to $(includedir)/$(PACKAGE_NAME) instead
fastnloincludedir = $(includedir)/fnlo_int_nlojet
fastnloinclude_HEADERS = fnlo_int_hhc_nlojet.h pdf-hhc-dummy.h cone-e.h cone-et.h kt-e.h fj-jets.h fastjet-jets.h

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
