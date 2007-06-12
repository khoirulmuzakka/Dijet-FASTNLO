# Makefile for fastNLO beta-release


NAME=example1


# set the path to LHAPDF and CERNLIB
LHAPDF=/disk2/work/wobisch/LHAPDFv4.1/lib
CERNLIB=/D0/ups/cern/Linux-2-4/2002/lib



FFILES= example1.f fnb0001.f fn-pdf.f fn-alphas-demo.f  fn-interface.f
FFILES2= gz.F 
CFILES= ftn_gzio.c

FFLAGS= -O -fno-automatic -ffixed-line-length-132
FC=f77


SRCFILES = $(FFILES) $(FFILES2) $(CFILES)
OBJFILES = $(FFILES:.f=.o) $(FFILES2:.F=.o) $(CFILES:.c=.o)
OBJFILESD = $(FFILESD:.f=.o) $(FFILES2:.F=.o) $(CFILES:.c=.o)

$(NAME): $(OBJFILES)
	 $(FC) -o $@ $(FFLAGS) $(OBJFILES) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB)  -lmathlib -lpacklib -lkernlib -lz
cln:
	'rm' -f *.o
	'rm' -f *.so

clean:
	'rm' -f *.o
	'rm' -f *.so
