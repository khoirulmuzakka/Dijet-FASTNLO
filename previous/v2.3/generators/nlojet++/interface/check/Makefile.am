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

#
# Targets
#
# perl script targets to check the package (distributed, but not installed)
#
dist_check_SCRIPTS = fnlo-nj-wrmtest.pl fnlo-nj-lotest.pl fnlo-nj-nlotest.pl
TESTS = $(dist_check_SCRIPTS)

# Additional clean up for check targets
CLEANFILES = *.log *.tab *.str *.txt *.wrm
clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf output

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but not for distribution
MAINTAINERCLEANFILES = Makefile.in
