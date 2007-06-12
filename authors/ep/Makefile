#
# ***  create fastNLO user job
#

CommonPath        = ../common

cteq5pdf = $(CommonPath)/Ctq5Par.f \
	   $(CommonPath)/cteq5m.cc

cteq6pdf = $(CommonPath)/Cteq61Pdf.f \
	   $(CommonPath)/cteq6.cc

export CPATH = $(CommonPath)

lib: fnh2001-tst01

CommonDir: ./common

./common: 
	ln -s $(CommonPath) ./common

# H1 incl. jets 820 GeV
fnh1001-tst01: CommonDir  fnh1001-tst01.cc kt-et-10.cc
	$(NLOJET)/bin/create-nlojet-user fnh1001-tst01.cc \
	kt-et-10.cc $(CommonPath)/alew.f -o fnh1001-tst01 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# H1 incl. jets 920 GeV
fnh2001-tst01: CommonDir  fnh2001-tst01.cc kt-et-10.cc
	$(NLOJET)/bin/create-nlojet-user fnh2001-tst01.cc \
	kt-et-10.cc $(CommonPath)/alew.f -o fnh2001-tst01 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# H1 dijets cross section 920 GeV
fnh2002-tst01: CommonDir  fnh2002-tst01.cc kt-et-10.cc
	$(NLOJET)/bin/create-nlojet-user fnh2002-tst01.cc \
	kt-et-10.cc $(CommonPath)/alew.f -o fnh2002-tst01 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
# H1 trijets cross section 920 GeV
fnh2003-tst01: CommonDir  fnh2003-tst01.cc kt-et-10.cc
	$(NLOJET)/bin/create-nlojet-user fnh2003-tst01.cc \
	kt-et-10.cc $(CommonPath)/alew.f -o fnh2003-tst01 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

nlofast-add-dis: nlofast-add-dis.o $(CommonPath)/gzstream.o
	g++ -o $(NLOJET)/bin/nlofast-add-dis nlofast-add-dis.o $(CommonPath)/gzstream.o -lz

clean:
	rm -rf .libs .obj *~ *.la *.o
