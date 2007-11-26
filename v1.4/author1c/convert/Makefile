#
# Makefile for table conversion/fix tools
#  M. Wobisch Feb 22, 2006
#
# modified by: K. Rabbertz, 05.09.2007
#
FFLAGS= -O -fno-automatic -ffixed-line-length-132
FC=g77

OBJFILES = $(FFILES:.f=.o)

#FFILES= convert1004.f
#FFILES= convert1004b.f
#FFILES= convert2003.f
#FFILES= convert2004.f
# Remove rap-pt bins without statistic (phase space)
FFILES= tablefix0002.f

# -> new standalone code 
convert: $(OBJFILES)
	$(FC) -o $@ $(FFLAGS) $(OBJFILES) 
