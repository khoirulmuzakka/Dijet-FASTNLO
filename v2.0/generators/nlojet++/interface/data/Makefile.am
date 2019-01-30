#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 12.02.2015
#
# NLOJet++ interface to fastNLO toolkit
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# Subdirectories
SUBDIRS = AzimuthalDecorrelation DijetAngular DijetMass Examples InclusiveJets InclusiveNJetDiffs InclusiveNJetEvents RivetAdditions check
#SUBDIRS = Examples InclusiveJets DijetMass RivetAdditions

# Install all data files from gzipped tar archive to avoid limits in command line length!
FIRST   = RivetAdditions
MYDIRS  = Examples InclusiveJets DijetMass DijetAngular AzimuthalDecorrelation InclusiveNJetEvents InclusiveNJetDiffs
ARCHIVE = Steering
TAR     = $(ARCHIVE).tar
TARGZ   = $(TAR).gz
# and distribute ...
pkgdata_DATA = $(TARGZ)
EXTRA_DIST = $(TARGZ)

# Create archive
$(TARGZ):
	tar cf $(TAR) $(FIRST)/*.info $(FIRST)/*.plot $(FIRST)/*.yoda; \
	for i in $(MYDIRS); do \
		echo $$i; \
		tar rf $(TAR) $$i/*.str; \
	done; \
	gzip $(TAR)

# Install archive hook
install-data-hook:
	mv $(TARGZ) $(DESTDIR)$(pkgdatadir); \
	cd $(DESTDIR)$(pkgdatadir); \
	tar xfz $(TARGZ); \
	rm -f $(TARGZ);

# Uninstall archive hook
uninstall-local:
	cd $(DESTDIR)$(pkgdatadir); \
	for i in $(FIRST) $(MYDIRS); do \
		rm -f $$i/*.info $$i/*.plot $$i/*.yoda $$i/*.str; \
		rmdir $$i; \
	done

# Built by make --> delete with make clean
MOSTLYCLEANFILES = $(TARGZ)

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = Makefile.in
