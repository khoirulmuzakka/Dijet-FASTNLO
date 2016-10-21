#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 20.10.2016
#
# fastNLO_toolkit
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in



check_SCRIPTS = fnlo-tk-cpptest.pl fnlo-tk-cattest.pl
TESTS = $(check_SCRIPTS)
fnlo-tk-cpptest.pl:
	./fnlo-tk-cpptest.pl
fnlo-tk-cattest.pl:
	./fnlo-tk-cattest.pl

CLEANFILES = *.log *.tab
