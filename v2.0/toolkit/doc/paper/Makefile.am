#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 06.04.2014
#
# Include documentation into distribution tarball
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu

EXTRA_DIST =

# Normally, it's not necessary to define docdir! Accommodate feature
# lacking in autoconf-2.59, automake-1.9.6
#docdir = ${datadir}/doc/${PACKAGE}
docfiles = fnlo-paper.pdf

DOCNAM = fnlo-paper
DOCSRC = $(DOCNAM).tex
DOCAUX = $(DOCNAM).aux
DOCPDF = $(DOCNAM).pdf

EXTRA_DIST += $(DOCSRC) figs fnlo-paper.pdf
CLEANFILES = $(DOCAUX) $(DOCPDF) $(DOCNAM).log

#
# Requires pdflatex or latex & dvipdfm
#
if HAVE_PDFLATEX
$(DOCPDF): $(srcdir)/$(DOCSRC)
	$(PDFLATEX) $<
	$(PDFLATEX) $<
	$(PDFLATEX) $<
else
if HAVE_LATEXDVIPDFM
$(DOCPDF): $(srcdir)/$(DOCSRC)
	$(LATEX) $<
	$(LATEX) $<
	$(LATEX) $<
	$(DVIPDFM) $<
endif
endif

dist_doc_DATA = $(docfiles)

# Files that are required to build the project but that can be
# recreated IF the necessary autotools etc. are available by:
# "autoreconf -i; ./configure"
# Use "make maintainer-clean" (on a build machine) before svn commits,
# but NOT before distribution with "make dist"!
MAINTAINERCLEANFILES = Makefile.in
