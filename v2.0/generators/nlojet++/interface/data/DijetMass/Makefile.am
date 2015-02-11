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
fnlo_steering_files = fnl2412c_I895742_v23_fix.str fnl2412c_I895742_v23_flex.str fnl2452_I1082936_v23_fix.str fnl2452_I1082936_v23_flex.str fnl2452_I1082936_v23_fix_pt12ave.str fnl2452_I1082936_v23_fix_ptmaxexpystar.str
dist_pkgdata_DATA = $(fnlo_steering_files)
EXTRA_DIST = $(fnlo_steering_files)

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = Makefile.in

# Preps to avoid line length limit ...!
#
#STEERARCHIVE = steering.tar.gz
#
#pkgdata_DATA = $(STEERARCHIVE)
#EXTRA_DIST = $(STEERARCHIVE)#
#
#install-data-hook:
#	cd $(DESTDIR)$(pkgdatadir); \
#	cat $(GIANTARCHIVE) | gunzip | tar --list > uninstall_manifest.txt; \
#	cat $(GIANTARCHIVE) | gunzip | tar --no-same-owner --extract; \
#	rm --force $(GIANTARCHIVE); \
#	cat uninstall_manifest.txt | sed --expression='s/^\|$$/"/g' | xargs chmod a=rX,u+w
#
#uninstall-local:
#	cd $(DESTDIR)$(pkgdatadir); \
#	cat uninstall_manifest.txt | sed --expression='s/ /\\ /g' | xargs rm --force; \
#	rm --force uninstall_manifest.txt
