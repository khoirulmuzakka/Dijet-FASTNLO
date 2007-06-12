#
# Makefile for threshold correction code for fastNLO v1c
#  M. Wobisch Dec 12, 2005
#

NAME = fnthresh-v14


FFILES= fnthresh-v14.f Cteq61Pdf.f alpha_cwz.o 
FFLAGS= -O -fno-automatic -ffixed-line-length-132
FC=g77

OBJFILES = $(FFILES:.f=.o)

FFILES2= fnthresh-v14-fnt1004.f Cteq61Pdf.f alpha_cwz.o 
OBJFILES2 = $(FFILES2:.f=.o)


# -> new standalone code 
$(NAME): $(OBJFILES)
	$(FC) -o $@ $(FFLAGS) $(OBJFILES) \
        -L $(CERNLIB)  -lmathlib -lpacklib -lkernlib



# only for scenario fnt1004  where 2nd rap bin
#     contains weighted x-sect
fnthresh-v14-fnt1004: $(OBJFILES2)
	$(FC) -o $@ $(FFLAGS) $(OBJFILES2) \
        -L $(CERNLIB)  -lmathlib -lpacklib -lkernlib


