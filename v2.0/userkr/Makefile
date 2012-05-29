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
NAME= allunc

ALLFILES= allunc.f fnx9999.f fn-interface.f alphas-grv.f alphas-grv-fnlo14.f alphst.f \
	pyalps.f fx9999rw.f fnio.f uncert.f fnset.f

FFLAGS= -O -fno-automatic
#FFLAGS= -O -fno-automatic -ffixed-line-length-132 -pedantic -std=f95
#FFLAGS= -O -fno-automatic -ffixed-line-length-132 -pedantic -Wunused
#FFLAGS= -O -fno-automatic -ffixed-line-length-72 -pedantic -Wunused

# Ubuntu 10.04 cms with Ubuntu CERNLIB
FC=gfortran
# Scientific Linux 4
#FC=g77 ?
# Scientific Linux 5
#FC=???

OBJFILESALL  = $(ALLFILES:.f=.o)

# Derive all uncertainties in one go
allunc:	$(OBJFILESALL) fnx9999.inc strings.inc uncert.inc
	$(FC) -o $@ $(FFLAGS) $(OBJFILESALL) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib

clean:
	'rm' -f *.o
	'rm' -f *.so
