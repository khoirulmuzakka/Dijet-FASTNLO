# Makefile for FASTNLO

NAME=example-v14a


#
# to use the CTEQ PDF code instead of LHAPDF:
#  - replace  fn-interface.f  with   fn-interface-cteq6.f
#  - include Cteq61Pdf.f  (from authorcode directory)
#  - remove LHAPDF from the Makefile below
#
# ---> maybe no longer needed with LHAPDF v5?? -> check
#

FFILES= example.f fnx9999.f fn-hbook.f fn-alphas-demo.f fn-interface.f
FFILESD= exampled01.f fastnlod00.f fn-pdf.f fn-alphas-demo.f  fn-interface.f
FFILES2= gz.F 
CFILES= ftn_gzio.c

FFILESERR= staterr.f fnx9999.f fn-hbook.f fn-alphas-demo.f  fn-interface.f

# --- executable at CEDAR web page - so far w/o LHAPDF and w/o alphas choice
FFILESCEDAR=example-cedar.f fnx9999.f fn-alphas-demo.f fn-interface-cedar.f

# --- for reference -> use original CTEQ code - works now with v1.4
REFFILES=example-ref.f fnx9999.f fn-hbook.f fn-alphas-demo.f  fn-interface-cteq6.f Cteq61Pdf.f

PDFFILES= pdfuncertainty.f fnx9999.f fn-hbook.f fn-alphas-demo.f fn-interface.f

FFLAGS= -O -fno-automatic -ffixed-line-length-132

#FC=f77
FC=g77
#FC=gfortran

SRCFILES = $(FFILES) $(FFILES2) $(CFILES)
OBJFILES = $(FFILES:.f=.o) $(FFILES2:.F=.o) $(CFILES:.c=.o)
OBJFILESD = $(FFILESD:.f=.o) $(FFILES2:.F=.o) $(CFILES:.c=.o)
OBJFILESERR = $(FFILESERR:.f=.o) $(FFILES2:.F=.o) $(CFILES:.c=.o)
OBJFILESREF = $(REFFILES:.f=.o) $(FFILES2:.F=.o) $(CFILES:.c=.o)
OBJFILESPDF = $(PDFFILES:.f=.o) $(CFILES:.c=.o)
OBJFILESCEDAR = $(FFILESCEDAR:.f=.o) 

# --- hadron-hadron
# - standard job
$(NAME): $(OBJFILES)
	 $(FC) -o $@ $(FFLAGS) $(OBJFILES) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz


# --- for statistical error analysis
staterr-v14a: $(OBJFILESERR)
	 $(FC) -o $@ $(FFLAGS) $(OBJFILESERR) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz

# --- compute PDF uncertainties
pdfunc:$(OBJFILESPDF)
	 $(FC) -o $@ $(FFLAGS) $(OBJFILESPDF) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz

#
# ----- stuff below not yet in v14a
#

# - for reference jobs one needs to use the original CTEQ PDF code
example-ref: $(OBJFILESREF)
	 $(FC) -o $@ $(FFLAGS) $(OBJFILESREF) \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz

example-cedar: $(OBJFILESCEDAR)
	 $(FC) -o $@ $(FFLAGS) $(OBJFILESCEDAR) \
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
