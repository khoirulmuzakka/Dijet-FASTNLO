#######################################################
#
#     Makefile for fastNLO usercode  v2.0
#
#######################################################


NAME=example



FFILES= example.f fnx9999.f fnset.f fn-alphas-demo.f fn-interface.f
FFLAGS= -O -fno-automatic -ffixed-line-length-132

# --- use -pg   for "gprof"
#FFLAGS= -O -fno-automatic -ffixed-line-length-132 -pg

# FC=g77
# FC=gfortran
FC=f95 

SRCFILES = $(FFILES) 
OBJFILES = $(FFILES:.f=.o)

$(NAME): $(OBJFILES)
	 $(FC) -o $@ $(FFLAGS) $(OBJFILES) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib

cln:
	'rm' -f *.o
	'rm' -f *.so

clean:
	'rm' -f *.o
	'rm' -f *.so
