#
# ***  create fastNLO author or standalone job 
#

CommonPath        = ../common

cteq5pdf = $(CommonPath)/Ctq5Par.f \
	   $(CommonPath)/cteq5m.cc

cteq6pdf = $(CommonPath)/Cteq61Pdf.f \
	   $(CommonPath)/cteq6.cc

search_07 = $(CommonPath)/search-e-07.cc \
	    $(CommonPath)/search-int.f \
	    $(CommonPath)/pxcone_mod_new.f

searchcone = $(CommonPath)/MidPointAlgorithm.cc \
	$(CommonPath)/PhysicsTower.cc \
	$(CommonPath)/LorentzVector.cc \
	$(CommonPath)/ClusterComparisons.cc\
	$(CommonPath)/Cluster.cc


export CPATH = $(CommonPath)

links: ./common ctq61.00.tbl 

./common:
	ln -s $(CommonPath) ./common

ctq61.00.tbl:
	ln -s $(CommonPath)/ctq61.00.tbl ctq61.00.tbl

#
# ======================================================================
# =====         FNT1001             ====================================
# ======================================================================
# -------------------------------- with Rsep=1.3
# Tevatron - Scenario FNT1001 (with Rsep=1.3)
fnt1001rsep: links fnt1001rsep.cc $(CommonPath)/rsep-et-07.cc \
	$(cteq5pdf) $(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1001rsep.cc \
	$(CommonPath)/rsep-et-07.cc $(cteq5pdf) $(cteq6pdf) \
	-o fnt1001rsep -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# Tevatron - Scenario FNT1001 (w/o Rsep)
fnt1001midp: links fnt1001midp.cc $(CommonPath)/cone-et-07.cc \
	$(cteq5pdf) $(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1001midp.cc \
	$(CommonPath)/cone-et-07.cc $(cteq5pdf) $(cteq6pdf) \
	-o fnt1001midp -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# Tevatron - Scenario FNT1001 (w/o Rsep) withut CDF bin, but with and without pdfs
fnt1001midp-ref: links fnt1001midp-ref.cc $(CommonPath)/cone-et-07.cc \
	$(cteq5pdf) $(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1001midp-ref.cc \
	$(CommonPath)/cone-et-07.cc $(cteq5pdf) $(cteq6pdf) \
	-o fnt1001midp-ref -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3


#
# ======================================================================
# =====         FNT1002             ====================================
# ======================================================================
# D0 - Run I incl jets / ratio 630GeV/1800GeV

fnt1002-630: links fnt1002-630.cc $(CommonPath)/cone-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1002-630.cc \
	$(CommonPath)/cone-et-07.cc $(cteq6pdf) \
	-o fnt1002-630 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnt1002-1800: links fnt1002-1800.cc $(CommonPath)/cone-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1002-1800.cc \
	$(CommonPath)/cone-et-07.cc $(cteq6pdf) \
	-o fnt1002-1800 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3


#
# ======================================================================
# =====         FNT200a             ====================================
# ======================================================================
# D0 - prel Run II incl jets (flexible pT binning)
fnt200a: links fnt200a.cc $(CommonPath)/cone-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt200a.cc \
	$(CommonPath)/cone-e-07.cc $(cteq6pdf) \
	-o fnt200a -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

#
# ======================================================================
# =====         FNR0001             ====================================
# ======================================================================
# RHIC - Scenario FNR0001 (3 algos)
fnr0001: links fnr0001.cc $(CommonPath)/kt2jet-et-10.cc \
	$(CommonPath)/cone-et-10.cc $(CommonPath)/cone-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnr0001.cc \
	$(CommonPath)/kt2jet-et-10.cc $(CommonPath)/cone-et-10.cc\
	$(CommonPath)/cone-et-07.cc $(cteq6pdf) \
	-o fnr0001 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

#
# ======================================================================
# =====         nlofast-add             ================================
# ======================================================================
#
#

# nlofast add for table addition
nlofast-add: nlofast-add.o $(CommonPath)/gzstream.o
	g++ -o $(NLOJET)/bin/nlofast-add nlofast-add.o \
	$(CommonPath)/gzstream.o -lz
#
# ======================================================================
# =====         txt2raw                 ================================
# ======================================================================
#
#

# convert ascii to raw format
txt2raw: txt2raw.o 
	g++ -o $(NLOJET)/bin/txt2raw txt2raw.o 

clean:
	rm -rf .libs .obj *~ *.la *.o $(CommonPath)/*.o ctq61.00.tbl ./common
