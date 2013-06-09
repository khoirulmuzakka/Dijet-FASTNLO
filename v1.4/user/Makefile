#
# My (KR) Makefile for all uncertainties with allunc
# Requirements:
#   - Normal sum table, e.g. fnl????.tab, of all LO, NLO, ... jobs and
#     PDF with uncertainty sets from LHAPDF
#     ==> cross sections, PDF uncertainty, scale uncertainties
# Extra:
#   - Sum table of scenario from reference calculations, e.g. fnl1007ref.tab,
#     ==> algorithmic uncertainty
#   - Collection of sum tables for each NLO job combined with all LO jobs
#     ==> statistical uncertainty of NLO calculation
#
# 10.06.2011 kr: Replace CERNLIB function LENOCC by f95 Standard LEN_TRIM
#                CERNLIB only needed for HBOOK histograms
#
NAME= allunc.old

ALLFILES= allunc.f fnx9999.f fn-interface.f fn-alphas-demo.f alphst.f \
	pyalps.f fn-hbook.f uncert.f
EXAFILES= example.f fnx9999.f fn-interface.f fn-alphas-demo.f alphst.f \
	pyalps.f fn-hbook.f
REFFILES= example-ref.f fnx9999.f fn-interface-cteq6.f fn-alphas-demo.f alphst.f \
	pyalps.f fn-hbook.f Cteq61Pdf.f

#FFLAGS= -O -fno-automatic
FFLAGS= -O -fno-automatic -ffixed-line-length-72
#FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic
#FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic -Wunused
#FFLAGS= -O -fno-automatic -ffixed-line-length-132

# Scientific Linux 4
#FC=g77
# Scientific Linux 5
#FC=???
# Ubuntu 10.04
FC=gfortran

OBJFILESALL  = $(ALLFILES:.f=.o)
OBJFILESEXA  = $(EXAFILES:.f=.o)
OBJFILESREF  = $(REFFILES:.f=.o)

# Derive all uncertainties in one go
allunc.old:	$(OBJFILESALL) fnx9999.inc strings.inc uncert.inc
	$(FC) -o $@ $(FFLAGS) $(OBJFILESALL) \
	-L $(LHAPDF)/lib -lLHAPDF \
	-L $(CERNLIB)/lib -lmathlib -lpacklib -lkernlib

# Example job
example: $(OBJFILESEXA) fnx9999.inc strings.inc
	 $(FC) -o $@ $(FFLAGS) $(OBJFILESEXA) \
	-L $(LHAPDF)/lib -lLHAPDF \
	-L $(CERNLIB)/lib -lmathlib -lpacklib -lkernlib

# Analyze algorithmic uncertainty from reference jobs with
# original CTEQ PDF code
example-ref: $(OBJFILESREF) fnx9999.inc
	 $(FC) -o $@ $(FFLAGS) $(OBJFILESREF) \
	-L $(CERNLIB)/lib -lmathlib -lpacklib -lkernlib

# Build everyhing
EXEFILES= allunc.old example example-ref
all: $(EXEFILES)



clean:
	'rm' -f *.o
	'rm' -f *.so

maintainer-clean: clean
	'rm' -f $(EXEFILES)
