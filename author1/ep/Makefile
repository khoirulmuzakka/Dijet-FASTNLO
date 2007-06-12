#
# ***  create fastNLO user job
#

CommonPath        = ../common

cteq5pdf = $(CommonPath)/Ctq5Par.f \
	   $(CommonPath)/cteq5m.cc

cteq6pdf = $(CommonPath)/Cteq61Pdf.f \
	   $(CommonPath)/cteq6.cc

export CPATH = $(CommonPath)

lib: fnh2001

CommonDir: ./common

./common: 
	ln -s $(CommonPath) ./common

# H1 incl. jets 920 GeV
fnh2001: CommonDir  fnh2001.cc kt-et-10.cc
	$(NLOJET)/bin/create-nlojet-user fnh2001.cc \
	kt-et-10.cc $(CommonPath)/alew.f -o fnh2001 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

nlofast-add-dis: nlofast-add-dis.o $(CommonPath)/gzstream.o
	g++ -o $(NLOJET)/bin/nlofast-add-dis nlofast-add-dis.o $(CommonPath)/gzstream.o -lz

clean:
	rm -rf .libs .obj *~ *.la *.o
