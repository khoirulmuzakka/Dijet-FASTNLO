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

all: fnt1001midp

links: ./common ctq61.00.tbl 

./common:
	ln -s $(CommonPath) ./common

ctq61.00.tbl:
	ln -s $(CommonPath)/ctq61.00.tbl ctq61.00.tbl

#
# ======================================================================
# =====         FNT1001             ====================================
# ======================================================================
# -------------------------------- 

# D0 - Scenario FNT1001 (w/o Rsep)
fnt1001midp: links fnt1001midp.cc $(CommonPath)/cone-et-07.cc
	$(NLOJET)/bin/create-nlojet-user fnt1001midp.cc \
	$(CommonPath)/cone-et-07.cc \
	-o fnt1001midp -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

clean:
	rm -rf .libs .obj *~ *.la *.o $(CommonPath)/*.o ctq61.00.tbl ./common
