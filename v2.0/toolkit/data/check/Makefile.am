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
dist_check_DATA = NJetEvents_0-0-2.tab.gz NJetEvents_0-0-2_cpp.log NJetEvents_2-2-2.tab.gz NJetEvents_2-2-2_cpp.log \
NJetEvents_2-2.tab.gz NJetEvents_2-2_cpp.log NJetEvents_norm_1-12.tab.gz NJetEvents_norm_1-12_cppnorm.log \
NJets.tab.gz NJets_born-2jet_0000.tab.gz NJets_born-2jet_0001.tab.gz NJets_born-2jet_0002.tab.gz \
NJets_born-2jet_stat.log NJets_born-2jet_stat_0000.tab.gz NJets_born-2jet_stat_0001.tab.gz NJets_born-2jet_stat_0002.tab.gz \
NJets_nlo-2jet_0100.tab.gz NJets_nlo-2jet_0101.tab.gz NJets_nlo-2jet_0102.tab.gz \
NJets_nlo-2jet_stat.log NJets_nlo-2jet_stat_0100.tab.gz NJets_nlo-2jet_stat_0101.tab.gz NJets_nlo-2jet_stat_0102.tab.gz \
SteerModify_EraseBins.str SteerModify_EraseBins2.str fnl1014_I902309_2-2.tab.gz fnl1014_I902309_2-2_cpp.log \
fnl1014_I902309_2-2_cppscales.log fnl1014_I902309_2.tab.gz zipappdiff.log

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = Makefile.in
