# Makefile for FASTNLO

NAME=example1


#
# to use the CTEQ PDF code instead of LHAPDF:
#  - replace  fn-interface.f  with   fn-interface-cteq6.f
#  - include Cteq61Pdf.f  (from authorcode directory)
#  - remove LHAPDF from the Makefile below
#

FFILES= example01.f fnx9999.f fastnlo00-hbook.f fn-pdf.f fn-alphas-demo.f  fn-interface.f tsh0001c.f
FFILESD= exampled01.f fastnlod00.f fn-pdf.f fn-alphas-demo.f  fn-interface.f
FFILES2= gz.F 
CFILES= ftn_gzio.c

FFILESERR= staterr.f fnx9999.f fastnlo00-hbook.f fn-pdf.f fn-alphas-demo.f  fn-interface.f


FFLAGS= -O -fno-automatic -ffixed-line-length-132

FC=f77

SRCFILES = $(FFILES) $(FFILES2) $(CFILES)
OBJFILES = $(FFILES:.f=.o) $(FFILES2:.F=.o) $(CFILES:.c=.o)
OBJFILESD = $(FFILESD:.f=.o) $(FFILES2:.F=.o) $(CFILES:.c=.o)
OBJFILESERR = $(FFILESERR:.f=.o) $(FFILES2:.F=.o) $(CFILES:.c=.o)

# --- hadron-hadron
$(NAME): $(OBJFILES)
	 $(FC) -o $@ $(FFLAGS) $(OBJFILES) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz


# --- for statistical error analysis
staterr: $(OBJFILESERR)
	 $(FC) -o $@ $(FFLAGS) $(OBJFILESERR) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz

# --- DIS
exampled1: $(OBJFILESD)
	 $(FC) -o $@ $(FFLAGS) $(OBJFILESD) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz

cln:
	'rm' -f *.o
	'rm' -f *.so

clean:
	'rm' -f *.o
	'rm' -f *.so
