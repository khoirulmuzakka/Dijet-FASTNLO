#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 15.05.2014
#
# fastNLO_toolkit
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Include files to be distributed
#
fastnloincludedir = $(includedir)
fastnloinclude_HEADERS = cone-e.h pdf-hhc-dummy.h fastNLOInterfaceToNLOJET.cc

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
