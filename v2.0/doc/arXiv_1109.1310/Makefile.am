#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 28.11.2011
#
# Include documentation into distribution tarball
#
# -----------------------------------------------------------------------------
# Use GNU standards
AUTOMAKE_OPTIONS = gnu
#
EXTRA_DIST =
#
# Normally, it's not necessary to define docdir! Accommodate feature
# lacking in autoconf-2.59, automake-1.9.6
docdir = ${datadir}/doc/${PACKAGE}
docfiles = alljet.pdf

DOCNAM = alljet
DOCSRC = $(DOCNAM).tex
DOCAUX = $(DOCNAM).aux
DOCPDF = $(DOCNAM).pdf

EXTRA_DIST += $(DOCSRC) psfig.sty figs 1109.1310v1.pdf
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
