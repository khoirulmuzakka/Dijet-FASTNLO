ObjSuf        = o
SrcSuf        = C
DepSuf        = D
ExeSuf        =
DllSuf        = so
LibDir        = $(HOME)/lib
BinDir        = $(HOME)/bin
HtmlDir       = $(HOME)/www/root

OutPutOpt     = -o

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --new --libs)
ROOTGLIBS     = $(shell root-config --glibs)

GSLINC        = -I$(shell gsl-config --prefix)/include

CXX           = g++
CXXFLAGS      = -O2 -fPIC $(GSLINC)
LD            = g++
SOFLAGS       = -shared
F2CFLAG       = -lg2c

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS)
GSLLIBS       = $(shell gsl-config --libs) -lgslblas

LHAPDFLIB     = /h1/h1gen/lhapdf/lib/liblhapdf.a
#------------------------------------------------------------------------------
 
FASTNLOO  = TIncljets.$(ObjSuf) \
	LHAPDFWrap.$(ObjSuf) \
	LHAPDFfw.$(ObjSuf) \
	fastnloDict.$(ObjSuf)

FASTNLOS  = TIncljets.$(SrcSuf) \
	fastnloDict.$(SrcSuf)

FASTNLODICTH  = TIncljets.h

FASTNLOSO  = $(LibDir)/libfastnlo.$(DllSuf)

.PHONY: lib  all html clean clobber

lib:            $(FASTNLOSO)

all:            lib  html

$(FASTNLOSO):       $(FASTNLOO)
	        $(LD) $(SOFLAGS) $(FASTNLOO) $(LHAPDFLIB) $(F2CFLAG) $(OutPutOpt) $(FASTNLOSO)  
		@echo "$(FASTNLOSO) done"

html:
	@echo "Generating web pages..."
	@root -l -b -q makehtml.C >/dev/null && cp -r html $(HtmlDir)  && rm -rf html && echo "done"

clean:
	@rm -f *.o *Dict.*

clobber: clean
	@rm -f $(FASTNLOSO)

#This replaces make depend, see also %.$(DepSuf): %.$(SrcSuf) rule below
include $(FASTNLOS:.$(SrcSuf)=.$(DepSuf))


.SUFFIXES: .$(SrcSuf)

###
fastnloDict.$(SrcSuf): $(FASTNLODICTH) fastnloLinkdef.h
	@echo "Generating dictionary fastnloDict..."
	@$(ROOTSYS)/bin/rootcint  -f fastnloDict.$(SrcSuf) -c $(FASTNLODICTH) 

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<

%.$(DepSuf): %.$(SrcSuf)
	@set -e; $(CXX) -MM $(CXXFLAGS) $< \
	| sed 's/\($*\)\.o[:]*/\1.o $@ : /g' > $@; \
	[ -s $@ ] || rm -f $@

.f.c:
	f2c -r8 -C++ $<  
