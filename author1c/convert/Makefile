#
# Makefile for table conversion tools
#  M. Wobisch Feb 22, 2006
#


FFLAGS= -O -fno-automatic -ffixed-line-length-132
FC=g77

OBJFILES = $(FFILES:.f=.o)

FFILES= convert2003.f
OBJFILES2 = $(FFILES2:.f=.o)


# -> new standalone code 
convert2003: $(OBJFILES)
	$(FC) -o $@ $(FFLAGS) $(OBJFILES) 


