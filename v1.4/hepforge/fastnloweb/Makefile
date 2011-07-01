#
# Makefile for fastNLO usercode
#

NAME=fastnloweb

#
# to use the CTEQ PDF code instead of LHAPDF:
#  - replace  fn-interface.f  with   fn-interface-cteq6.f
#  - include Cteq61Pdf.f  (from authorcode directory)
#  - remove LHAPDF from the Makefile below
#

FFILES= main-cedar.f fnx9999.f fn-pdf.f fn-alphas-demo.f fn-interface-cedar.f
FFLAGS= -O -fno-automatic -ffixed-line-length-132
# --- use -gp   for "gprof"
#FFLAGS= -O -fno-automatic -ffixed-line-length-132 -pg

#LHAPDF=/hepforge/home/fastnlo/public_html/lhapdf-5.7.0/src/libLHAPDFCore.la
LHAPDF=/hepforge/home/fastnlo/local/lib/libLHAPDF.la

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
