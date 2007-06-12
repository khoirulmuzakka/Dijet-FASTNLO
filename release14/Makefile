#
# Makefile for fastNLO example program
# 
# the user needs to edit the variables LHAPDF and CERNLIB
#
#
#

LHAPDF=/work/shootingstar-clued0/wobisch/lhapdf500/lib
CERNLIB=/D0/ups/cern/Linux-2-4/2002/lib

NAME=fastnlo

FF0= fn-example.f fn-alphas-demo.f fn-interface.f
FF1= fnt1001.f fnt1002.f fnt1003.f fnt1004.f fnt1005.f fnt1007.f fnt1008.f
FF2= fnt2001.f fnt2002.f fnt2003.f  fnr0001.f fnl0004.f 
FF3= fnh1001.f fnh1002.f fnh1003.f fnh1004.f

FFLAGS= -O -fno-automatic -ffixed-line-length-132

FC=f77

OBJFILES = $(FF0:.f=.o) $(FF1:.F=.o) $(FF2:.F=.o) $(FF3:.F=.o)

$(NAME): $(OBJFILES)
	 $(FC) -o $@ $(FFLAGS) $(OBJFILES) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz

cln:
	'rm' -f *.o
	'rm' -f *.so

clean:
	'rm' -f *.o
	'rm' -f *.so
