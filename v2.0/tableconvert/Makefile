# ###############################################
#
# Makefile for fastNLO convert program
#                   v1.4 to v2.0 table 
#
# ###############################################



NAME=conv20

FF= conv20.f 
FFLAGS= -O -fno-automatic -ffixed-line-length-132

FC=f95
#FC=gfortran


OBJFILES = $(FF:.f=.o) 

$(NAME): $(OBJFILES)
#	 $(FC) -o $@ $(FFLAGS) $(OBJFILES) 
	 $(FC) -o $@ $(FFLAGS) $(OBJFILES) 

cln:
	'rm' -f *.o
	'rm' -f *.so

clean:
	'rm' -f *.o
	'rm' -f *.so
