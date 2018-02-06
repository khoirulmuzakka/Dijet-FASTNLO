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

# Distribute for checks, but do not install
# Link filenames to avoid cmdline lengths > 99 for tar --> error in distribution
#   steer.str   --> fnr0001midpHT_I723509_v23_fix.str
#   warmup.txt  --> fnr0001midpHT_I723509_v23_fix_InclusiveNJets_warmup.txt
#   born.tab.gz --> fnr0001midpHT_I723509_v23_fix-hhc-born-2jet.tab.gz
#   nlo.tab.gz  --> fnr0001midpHT_I723509_v23_fix-hhc-nlo-2jet.tab.gz
dist_check_DATA = steer.str warmup.txt born.tab.gz nlo.tab.gz

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = Makefile.in
