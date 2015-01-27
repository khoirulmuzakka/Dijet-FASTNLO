#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 27.01.2015
#
# Additions to provide fastNLO vs. data comparisons even though
# Rivet analysis not available
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

#
# Data files to be distributed
#
fnlorivetdir = $(datadir)/fastnlo_rivet_additions
rivet_analysis_files = ATLAS_2013_I1228693.info ATLAS_2013_I1228693.plot ATLAS_2013_I1228693.yoda \
	CMS_2011_I895742.info CMS_2011_I895742.plot CMS_2011_I895742.yoda
fnlorivet_DATA = $(rivet_analysis_files)
EXTRA_DIST     = $(rivet_analysis_files)

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = Makefile.in
