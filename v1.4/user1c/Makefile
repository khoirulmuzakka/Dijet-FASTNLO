#
# Makefile for PDF uncertainties
#
NAME=pdfunc

PDFFILES= pdfuncertainty.f fnx9999.f fn-hbook.f fn-interface.f \
	fn-alphas-demo.f alphst.f 
CFILES= ftn_gzio.c

FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic
#FFLAGS= -O -fno-automatic -ffixed-line-length-72 -fpedantic -Wunused
#FFLAGS= -O -fno-automatic -ffixed-line-length-132

FC=g77
#FC=f77
#FC=gfortran

OBJFILESPDF = $(PDFFILES:.f=.o) $(CFILES:.c=.o)

# Compute PDF uncertainties
pdfunc: $(OBJFILESPDF) fnx9999.inc strings.inc
	$(FC) -o $@ $(FFLAGS) $(OBJFILESPDF) \
	-L $(LHAPDF) -lLHAPDF \
	-L $(CERNLIB) -lmathlib -lpacklib -lkernlib -lz

clean:
	'rm' -f *.o
	'rm' -f *.so
