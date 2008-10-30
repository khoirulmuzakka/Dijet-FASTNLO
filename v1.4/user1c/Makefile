#
# My (KR) Makefile for algorithmic, statistical and PDF uncertainties
# Requirements:
#   - Sum table of scenario from reference calculations, e.g. fnl1007ref.tab,
#     ==> algorithmic uncertainty in histo nnnnnn4
#   - Collection of sum tables for each NLO job combined with all LO jobs  
#     ==> statistical uncertainty of NLO calculation in histo nnnnnn3
#   - Normal sum table of all LO, NLO, ... jobs and PDF with unc. sets
#     ==> PDF uncertainty, summed up upper, lower in histos nnnnnn1, nnnnnn2
#
NAME= allunc

CFILES= ftn_gzio.c
REFFILES= example-ref.f fnx9999.f fn-hbook.f fn-alphas-demo.f \
	fn-interface-cteq6.f Cteq61Pdf.f
STATFILES= staterr_kr.f fnx9999.f fn-interface.f fn-alphas-demo.f alphst.f
PDFFILES= pdfuncertainty.f fnx9999.f fn-interface.f fn-alphas-demo.f alphst.f
ALLFILES= allunc.f fnx9999.f fn-interface.f fn-alphas-demo.f alphst.f \
	fn-hbook.f statcode.f

FFLAGS= -O -fno-automatic -ffixed-line-length-72
#FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic
#FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic -Wunused
#FFLAGS= -O -fno-automatic -ffixed-line-length-132

FC=g77
#FC=f77
#FC=gfortran

#OBJFILESCOM  = $(COMFILES:.f=.o)
OBJFILESREF  = $(REFFILES:.f=.o) $(CFILES:.c=.o)
OBJFILESSTAT = $(STATFILES:.f=.o) $(CFILES:.c=.o)
OBJFILESPDF  = $(PDFFILES:.f=.o) $(CFILES:.c=.o)
OBJFILESALL  = $(ALLFILES:.f=.o) $(CFILES:.c=.o)

# Derive all uncertainties in one go
allunc:	$(OBJFILESALL) fnx9999.inc strings.inc
	$(FC) -o $@ $(FFLAGS) $(OBJFILESALL) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz

# Analyze algorithmic uncertainty from reference jobs with
# original CTEQ PDF code
example-ref: $(OBJFILESREF) fnx9999.inc
	 $(FC) -o $@ $(FFLAGS) $(OBJFILESREF) \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz

# My statistical error analysis (kr)
staterr: $(OBJFILESSTAT)
	 $(FC) -o $@ $(FFLAGS) $(OBJFILESSTAT) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz

# Compute PDF uncertainties
pdfunc: $(OBJFILESPDF) fnx9999.inc strings.inc
	$(FC) -o $@ $(FFLAGS) $(OBJFILESPDF) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz



# Todo, clean up common files.o
clean:
	'rm' -f *.o
	'rm' -f *.so
