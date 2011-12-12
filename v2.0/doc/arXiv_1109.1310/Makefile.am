#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 28.11.2011
#
# Include documentation into distribution tarball
#
# -----------------------------------------------------------------------------
docfiles = alljet.pdf

DOCNAM = alljet
DOCSRC = $(DOCNAM).tex
DOCAUX = $(DOCNAM).aux
DOCPDF = $(DOCNAM).pdf

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
