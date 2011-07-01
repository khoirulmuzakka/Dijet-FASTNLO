#
# Makefile for fastNLO web code at HepForge
#

NAME=fastnloweb

FFILES= main-cedar.f fnx9999.f fn-pdf.f fn-alphas-demo.f fn-interface-cedar.f
FFLAGS= -O -fno-automatic -ffixed-line-length-132

# At HepForge
LHAPDF=$(HOME)/fastnlo/local/lib/libLHAPDF.la
# Local test
#LHAPDF=$(HOME)/local/lhapdf-5.8.4/lib/libLHAPDF.la

FC=gfortran
LIBTOOL=libtool

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
