###########################################################################################
#
#  simple Makefile for FastNLOReader
#
###########################################################################################

OBJS 	= main.cc speaker.cc fastNLOInterpolBase.cc fastNLOInterpolCatmulRom.cc fastNLOBase.cc fastNLOTable.cc fastNLOCoefficients.cc fastNLOCreate.cc read_steer.cc
OBJS	+= fastNLOCoeffBase.cc fastNLOCoeffData.cc fastNLOCoeffMult.cc fastNLOCoeffAddBase.cc fastNLOCoeffAddPert.cc fastNLOCoeffAddFlex.cc
OBJS	+= fastNLOReader.cc fastNLOCoeffAddEvalBase.cc fastNLOCoeffAddFlexEval.cc fastNLOPDFLinearCombinations.cc
OBJS_LIB = speaker.o fastNLOInterpolBase.o fastNLOInterpolCatmulRom.o fastNLOBase.o fastNLOTable.o fastNLOCoefficients.o fastNLOCreate.o read_steer.o
OBJS_LIB +=  fastNLOCoeffBase.o fastNLOCoeffData.o fastNLOCoeffMult.o fastNLOCoeffAddBase.o fastNLOCoeffAddPert.o fastNLOCoeffAddFlex.o
OBJS_LIB += fastNLOReader.o fastNLOCoeffAddEvalBase.o fastNLOCoeffAddFlexEval.o fastNLOPDFLinearCombinations.o
#PlotContributions.o FastNLOBlockB.o FastNLOReader.o Alphas.o qcdnumdummy.o
#OBJS 	= PlotContributions.o FastNLOBlockB.cc FastNLOReader.cc Alphas.cc qcdnumdummy.f
CC 	= gcc
LDFLAGS	= -lm -lpthread `gtk-config --cflags` `gtk-config --libs` -lgthread

CFLAGS 	= -Wall -g -I$(LHAPDFINCLUDEPATH) -I$(ROOTSYS)/include -lstdc++ -O1

### reader sachen
LHAPDFVERSION	= 5.8.6
LHAPDFDIR	= /afs/desy.de/group/alliance/mcg/public/MCGenerators/lhapdf
LHAPDFLIBPATH	= $(LHAPDFDIR)/$(LHAPDFVERSION)/$(SYSNAME)/lib/
LHAPDFLIB	= LHAPDF
LHAPDFINCLUDEPATH	= $(LHAPDFDIR)/$(LHAPDFVERSION)/$(SYSNAME)/include


fnlo: $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) \
	-I$(LHAPDFINCLUDEPATH) \
	-L$(LHAPDFLIBPATH) -l$(LHAPDFLIB) \
	-Wl,-rpath -Wl, \
	-o fnlo

libfnlo.so: $(OBJS_LIB)
	$(CXX_COMP) $(OBJS_LIB) -shared -o ./$@ 
	cp libfnlo.so ../hadron/.

#	-L$(LHAPDFLIBPATH) -l$(LHAPDFLIB) \
#	-L$(ROOTLIB) -lCore -lGpad -lHist \
#$(LHAPDFLIBPATH) \

clean:
	rm -f *.o *.a *.so *Dict* 2>/dev/null

%.o: %.cc
	$(CC) -c $(CFLAGS) -o $@ $^


all: fnlo libfnlo.so

