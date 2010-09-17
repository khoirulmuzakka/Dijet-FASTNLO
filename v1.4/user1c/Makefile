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
NAME= allunc

ALLFILES= allunc.f fnx9999.f fn-interface.f fn-alphas-demo.f alphst.f \
	pyalps.f fn-hbook.f uncert.f

FFLAGS= -O -fno-automatic -ffixed-line-length-72
#FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic
#FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic -Wunused
#FFLAGS= -O -fno-automatic -ffixed-line-length-132

FC=g77
#FC=f77
#FC=gfortran

OBJFILESALL  = $(ALLFILES:.f=.o)

# Derive all uncertainties in one go
allunc:	$(OBJFILESALL) fnx9999.inc strings.inc uncert.inc
	$(FC) -o $@ $(FFLAGS) $(OBJFILESALL) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib

clean:
	'rm' -f *.o
	'rm' -f *.so
