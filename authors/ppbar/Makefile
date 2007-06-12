#
# ***  create fastNLO author or standalone job 
#

CommonPath        = ../common

cteq5pdf = $(CommonPath)/Ctq5Par.f \
	   $(CommonPath)/cteq5m.cc

cteq6pdf = $(CommonPath)/Cteq61Pdf.f \
	   $(CommonPath)/cteq6.cc

export CPATH = $(CommonPath)

links: ./common ctq61.00.tbl 

./common:
	ln -s $(CommonPath) ./common
ctq61.00.tbl:
	ln -s $(CommonPath)/ctq61.00.tbl ctq61.00.tbl

# ======================================================================
# ===== Tevatron - Scenario FNT1001 ====================================
# ======================================================================
# -------------------------------- with Rsep=1.3
# Tevatron - Scenario FNT1001 (with Rsep=1.3) - standalone version
fnt1001rsep-stdl: links fnt1001rsep-stdl.cc $(CommonPath)/rsep-et-07.cc \
	$(cteq5pdf) $(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1001rsep-stdl.cc \
	$(CommonPath)/rsep-et-07.cc $(cteq5pdf) $(cteq6pdf) \
	-o fnt1001rsep-stdl -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# -------------------------------- with Rsep=1.3 
# Tevatron - Scenario FNT1001 (with Rsep=1.3) - fast table version
fnt1001rsep: links fnt1001rsep.cc $(CommonPath)/rsep-et-07.cc 
	$(NLOJET)/bin/create-nlojet-user fnt1001rsep.cc \
	$(CommonPath)/rsep-et-07.cc \
	-o fnt1001rsep -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# -------------------------------- with Midpoint Algo
# Tevatron - Scenario FNT1001 (no Rsep) - standalone version
fnt1001midp-stdl: links fnt1001midp-stdl.cc $(CommonPath)/cone-et-07.cc \
	$(cteq5pdf) $(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1001midp-stdl.cc \
	$(CommonPath)/cone-et-07.cc  $(cteq5pdf) $(cteq6pdf) \
	-o fnt1001midp-stdl -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# -------------------------------- with Midpoint Algo
# Tevatron - Scenario FNT1001 (no Rsep) - fast table version
fnt1001midp: links fnt1001midp.cc $(CommonPath)/cone-et-07.cc 
	$(NLOJET)/bin/create-nlojet-user fnt1001midp.cc \
	$(CommonPath)/cone-et-07.cc \
	-o fnt1001midp -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# nlofast add for table addition
nlofast-add-v1: nlofast-add-v1.o $(CommonPath)/gzstream.o
	g++ -o $(NLOJET)/bin/nlofast-add-v1 nlofast-add-v1.o $(CommonPath)/gzstream.o -lz

nlofast-add: nlofast-add.o $(CommonPath)/gzstream.o
	g++ -o $(NLOJET)/bin/nlofast-add nlofast-add.o $(CommonPath)/gzstream.o -lz

clean:
	rm -rf .libs .obj *~ *.la *.o $(CommonPath)/*.o ctq61.00.tbl ./common
