#
# Makefile for fastNLO interactive web code at HepForge
#

NAME=fastnloweb
FC=gfortran
#FC=f77
LIBTOOL=libtool
FFLAGS= -O -fno-automatic -ffixed-line-length-132

# At HepForge
LHAPDF=$(HOME)/fastnlo/local/lib/libLHAPDF.la
# Local test
#LHAPDF=$(HOME)/local/lhapdf-5.8.4/lib/libLHAPDF.la

FFILES= main-cedar.f fnx9999.f fn-pdf.f fn-alphas-demo.f fn-interface-cedar.f
SRCFILES = $(FFILES) $(CFILES)
OBJFILES = $(FFILES:.f=.o) $(CFILES:.c=.o)

$(NAME): $(OBJFILES)
	 $(LIBTOOL) --mode=link $(FC) -o $@ $(FFLAGS) $(OBJFILES) \
	 $(LHAPDF)

cln:
	'rm' -f *.o
	'rm' -f *.so

clean:
	'rm' -f *.o
	'rm' -f *.so
