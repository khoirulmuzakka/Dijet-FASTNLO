#
# Makefile for threshold correction code for fastNLO
#  M. Wobisch Dec 12, 2005
#

NAME = fnthresh

FFILES= fnthresh.f Cteq61Pdf.f alpha_cwz.o 
FFLAGS= -O -fno-automatic -ffixed-line-length-132
FC=g77

OBJFILES = $(FFILES:.f=.o)


# -> new standalone code 
$(NAME): $(OBJFILES)
	$(FC) -o $@ $(FFLAGS) $(OBJFILES) \
        -L $(CERNLIB)  -lmathlib -lpacklib -lkernlib

