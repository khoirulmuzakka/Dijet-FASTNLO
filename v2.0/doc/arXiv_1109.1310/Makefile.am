#                                               -*- Automake -*-
# -----------------------------------------------------------------------------
#
# K. Rabbertz, 28.11.2011
#
# Include documentation into distribution tarball
#
# -----------------------------------------------------------------------------
docfiles = alljet.tex alljet.pdf

#
# Require pdflatex
#
if HAVE_PDFLATEX
DOCNAM = alljet
DOCSRC = $(DOCNAM).tex
DOCAUX = $(DOCNAM).aux
DOCPDF = $(DOCNAM).pdf

CLEANFILES = $(DOCAUX) $(DOCPDF) $(DOCNAM).log

$(DOCPDF): $(srcdir)/$(DOCSRC)
	$(PDFLATEX) $<
endif

dist_doc_DATA = $(docfiles)

#%pdflatexopts = -file-line-error -halt-on-error

#ALL:	$(DOC).pdf
#	make pdf

#$(DOC).pdf:	$(DOC).tex
#	latex $(DOC)
#	latex $(DOC)
#	latex $(DOC)
#	dvipdfm $(DOC)

#	bibtex   $(DOC).aux
#	pdflatex $(pdflatexopts) $(DOC).tex

#clean:
#	rm -f $(DOC).dvi $(DOC).ps $(DOC).pdf $(DOC).bbl $(DOC).blg \
#	$(DOC).log $(DOC).out *.toc *.aux
