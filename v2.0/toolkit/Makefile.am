#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 06.10.2013
#
# fastNLO_toolkit
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

# Subdirectories
SUBDIRS = data doc fastnlotoolkit scripts src pyext

# Provide local m4 macros in dir m4, copy also system-wide macros
ACLOCAL_AMFLAGS = -I m4 --install

# Add subversion revision number, recipe from Autotools by J. Calcote
EXTRA_DIST = m4 SVNREV
# This gives e.g. fastnlo-2.0.0-643M for svn revision 643 locally modified
distdir = $(PACKAGE)-$(VERSION)-$(SVNREV)
# Don't put .svn subdirs into distribution package
dist-hook:
	rm -rf `find $(distdir) -name .svn`

# Installation of DOXYGEN documentation
# TBD: Docu is always remade! Check tag file mechanism.
# TBD: Remove doxygen directory @ maintainer-clean.
if HAVE_DOXYGEN
if DX_COND_doc
include aminclude.am
EXTRA_DIST += aminclude.am doxygen.cfg
all-local-doxygen: doxygen-doc
	touch doc/doxygen/fastnlo_toolkit.tag
.PHONY: all-local-doxygen
all-local: all-local-doxygen
MOSTLYCLEANFILES = $(DX_CLEANFILES) all-local-doxygen

## This rule requires autoconf 2.62 or newer
install-data-local:
	@$(NORMAL_INSTALL)
	$(MKDIR_P) $(DESTDIR)$(datadir)/doc/$(PACKAGE)/html/search
	$(MKDIR_P) $(DESTDIR)$(datadir)/doc/$(PACKAGE)/latex
	$(INSTALL_DATA) doc/doxygen/html/search/* $(DESTDIR)$(datadir)/doc/$(PACKAGE)/html/search
	$(INSTALL_DATA) doc/doxygen/html/*.* $(DESTDIR)$(datadir)/doc/$(PACKAGE)/html
	$(INSTALL_DATA) doc/doxygen/latex/* $(DESTDIR)$(datadir)/doc/$(PACKAGE)/latex
	$(INSTALL_DATA) doc/doxygen/*.* $(DESTDIR)$(datadir)/doc/$(PACKAGE)

uninstall-local:
	@$(NORMAL_UNINSTALL)
	rm -f $(DESTDIR)$(datadir)/doc/$(PACKAGE)/html/search/*
	rmdir $(DESTDIR)$(datadir)/doc/$(PACKAGE)/html/search
	rm -f $(DESTDIR)$(datadir)/doc/$(PACKAGE)/html/*
	rmdir $(DESTDIR)$(datadir)/doc/$(PACKAGE)/html
	rm -f $(DESTDIR)$(datadir)/doc/$(PACKAGE)/latex/*
	rmdir $(DESTDIR)$(datadir)/doc/$(PACKAGE)/latex
	rm -f $(DESTDIR)$(datadir)/doc/$(PACKAGE)/*
	rmdir $(DESTDIR)$(datadir)/doc/$(PACKAGE)

endif DX_COND_doc
endif HAVE_DOXYGEN

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = INSTALL config.h.in config.guess config.sub configure Makefile.in \
	depcomp install-sh ltmain.sh missing compile \
	aclocal.m4 m4/libtool.m4 m4/ltoptions.m4 m4/ltsugar.m4 m4/ltversion.m4 m4/lt~obsolete.m4 \
	SVNREV py-compile
