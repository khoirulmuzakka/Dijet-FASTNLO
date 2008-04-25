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

midpt_e_05 = $(CommonPath)/midpt-e-05.cc \
             $(CommonPath)/px_interface.f \
             $(CommonPath)/pxcone_mod_new.f

midpt_e_07 = $(CommonPath)/midpt-e-07.cc \
             $(CommonPath)/px_interface.f \
             $(CommonPath)/pxcone_mod_new.f

export CPATH = $(CommonPath)

links: ./common ctq61.00.tbl 

./common:
	ln -s $(CommonPath) ./common

ctq61.00.tbl:
	ln -s $(CommonPath)/ctq61.00.tbl ctq61.00.tbl


#
# ======================================================================
# =====         FNT1001  v1c      CDF incl jets @ 1800GeV     ==========
# ======================================================================
# 
fnt1001midp: links fnt1001midp.cc $(CommonPath)/cone-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1001midp.cc \
	$(CommonPath)/cone-et-07.cc $(cteq6pdf) \
	-o fnt1001midp -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnt1001rsep: links fnt1001rsep.cc $(CommonPath)/rsep-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1001rsep.cc \
	$(CommonPath)/rsep-et-07.cc $(cteq6pdf) \
	-o fnt1001rsep -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3


#
# ======================================================================
# =====         FNT1002  v1c     D0 incl jets @ 1800GeV     ============
# ======================================================================
# 
fnt1002midp: links fnt1002midp.cc $(CommonPath)/cone-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1002midp.cc \
	$(CommonPath)/cone-et-07.cc $(cteq6pdf) \
	-o fnt1002midp -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnt1002rsep: links fnt1002rsep.cc $(CommonPath)/rsep-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1002rsep.cc \
	$(CommonPath)/rsep-et-07.cc $(cteq6pdf) \
	-o fnt1002rsep -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

#
# ======================================================================
# =====         FNT1003  v1c   CDF dijets (ET,eta1,eta2) @ 1800GeV  ====
# ======================================================================
# 
fnt1003midp: links fnt1003midp.cc $(CommonPath)/cone-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1003midp.cc \
	$(CommonPath)/cone-et-07.cc $(cteq6pdf) \
	-o fnt1003midp -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnt1003rsep: links fnt1003rsep.cc $(CommonPath)/rsep-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1003rsep.cc \
	$(CommonPath)/rsep-et-07.cc $(cteq6pdf) \
	-o fnt1003rsep -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

#
# ======================================================================
# =====         FNT1004/5  D0 incl jets ratio 630/1800GeV           ====
# ======================================================================
# 
fnt1004midp: links fnt1004midp.cc $(CommonPath)/cone-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1004midp.cc \
	$(CommonPath)/cone-et-07.cc $(cteq6pdf) \
	-o fnt1004midp -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnt1004rsep: links fnt1004rsep.cc $(CommonPath)/rsep-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1004rsep.cc \
	$(CommonPath)/rsep-et-07.cc $(cteq6pdf) \
	-o fnt1004rsep -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnt1005midp: links fnt1005midp.cc $(CommonPath)/cone-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1005midp.cc \
	$(CommonPath)/cone-et-07.cc $(cteq6pdf) \
	-o fnt1005midp -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnt1005rsep: links fnt1005rsep.cc $(CommonPath)/rsep-et-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1005rsep.cc \
	$(CommonPath)/rsep-et-07.cc $(cteq6pdf) \
	-o fnt1005rsep -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

#
# ======================================================================
# =====         FNT1007  CDF dijet-mass
# ======================================================================
# 
fnt1007midp: links fnt1007midp.cc $(CommonPath)/cone-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1007midp.cc \
	$(CommonPath)/cone-e-07.cc $(cteq6pdf) \
	-o fnt1007midp -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnt1007rsep: links fnt1007rsep.cc $(CommonPath)/rsep-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1007rsep.cc \
	$(CommonPath)/rsep-e-07.cc $(cteq6pdf) \
	-o fnt1007rsep -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

#
# ======================================================================
# =====      FNT2001  v1c   D0 Dijet Delta Phi @ 1960GeV  ======
# ======================================================================
# 
fnt2001norm: links fnt2001norm.cc $(CommonPath)/cone-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2001norm.cc \
	$(CommonPath)/cone-e-07.cc $(cteq6pdf) \
	-o fnt2001norm -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnt2001diff: links fnt2001diff.cc $(midpt_e_07) \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2001diff.cc \
	$(midpt_e_07) $(cteq6pdf) \
	-o fnt2001diff -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

#
# ======================================================================
# =====      FNT2002  v1c   CDF incl jets cone algo (pT) @ 1960GeV  ====
# ======================================================================
# 
fnt2002midp: links fnt2002midp.cc $(CommonPath)/cone-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2002midp.cc \
	$(CommonPath)/cone-e-07.cc $(cteq6pdf) \
	-o fnt2002midp -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnt2002rsep: links fnt2002rsep.cc $(CommonPath)/rsep-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2002rsep.cc \
	$(CommonPath)/rsep-e-07.cc $(cteq6pdf) \
	-o fnt2002rsep -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

#
# ======================================================================
# =====      FNT2003  v1c   CDF incl jets kT algo (pT) @ 1960GeV  ======
# ======================================================================
# 
fnt2003: links fnt2003.cc $(CommonPath)/kt-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2003.cc \
	$(CommonPath)/kt-e-07.cc $(cteq6pdf) \
	-o fnt2003 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
#
# ===  --> extended version - using multiple y ranges 
fnt2003plus: links fnt2003plus.cc $(CommonPath)/kt-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2003plus.cc \
	$(CommonPath)/kt-e-07.cc $(cteq6pdf) \
	-o fnt2003plus -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# --- ... as above but for D=0.5
fnt2003plusD05: links fnt2003plusD05.cc $(CommonPath)/kt-e-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2003plusD05.cc \
	$(CommonPath)/kt-e-05.cc $(cteq6pdf) \
	-o fnt2003plusD05 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# --- ... as above but for D=1.0
fnt2003plusD10: links fnt2003plusD10.cc $(CommonPath)/kt-e-10.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2003plusD10cc \
	$(CommonPath)/kt-e-10.cc $(cteq6pdf) \
	-o fnt2003plusD10 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3


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


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
#
# ======================================================================
# =====   FNL0001: LHC test kT D=0.7 ===================================
# ======================================================================
# first LHC scenario with fine binning
#
fnl0001: links fnl0001.cc $(CommonPath)/kt2jet-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0001.cc \
	$(CommonPath)/kt2jet-e-07.cc $(cteq6pdf) \
	-o fnl0001 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0002: LHC CMS incl. jets with kT D=1.0 ====================
# ======================================================================
# LHC scenario for CMS @ 14 TeV
#
fnl0002: links fnl0002.cc $(CommonPath)/kt-e-10.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0002.cc \
	$(CommonPath)/kt-e-10.cc $(cteq6pdf) \
	-o fnl0002 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0003: LHC CMS incl. jets with Midpoint R=0.7 ==============
# ======================================================================
# LHC scenario for CMS @ 14 TeV
#
fnl0003: links fnl0003.cc $(CommonPath)/cone-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0003.cc \
	$(CommonPath)/cone-e-07.cc $(cteq6pdf) \
	-o fnl0003 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0004: LHC ATLAS test kT D=0.7 =============================
# ======================================================================
# ATLAS test scenario
#
fnl0004: links fnl0004.cc $(CommonPath)/kt2jet-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0004.cc \
	$(CommonPath)/kt2jet-e-07.cc $(cteq6pdf) \
	-o fnl0004 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0005: LHC ATLAS test kT D=0.7 dijet mass ==================
# ======================================================================
# ATLAS test scenario for dijet mass
#
fnl0005: links fnl0005.cc $(CommonPath)/kt2jet-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0005.cc \
	$(CommonPath)/kt2jet-e-07.cc $(cteq6pdf) \
	-o fnl0005ref -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0006: LHC ATLAS test kT D=0.7 dijet chi ===================
# ======================================================================
# ATLAS test scenario with realistic binning for chi distribution
#
fnl0006: links fnl0006.cc $(CommonPath)/kt2jet-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0006.cc \
	$(CommonPath)/kt2jet-e-07.cc $(cteq6pdf) \
	-o fnl0006 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0007: LHC CMS incl. jets with kT D=0.6 ====================
# ======================================================================
# LHC scenario for CMS @ 14 TeV
#
fnl0007: links fnl0007.cc $(CommonPath)/kt-e-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0007.cc \
	$(CommonPath)/kt-e-06.cc $(cteq6pdf) \
	-o fnl0007 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0008: LHC CMS incl. jets with fastjet SISCone R=0.7 =======
# ======================================================================
# LHC scenario for CMS @ 14 TeV
#
fnl0008: links fnl0008.cc $(CommonPath)/fj-sc-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0008.cc \
	$(CommonPath)/fj-sc-07.cc $(cteq6pdf) \
	-o fnl0008 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL0009: LHC CMS incl. jets with Midpoint R=0.7 ==============
# ======================================================================
# LHC scenario for CMS @ 14 TeV
#
fnl0009: links fnl0009.cc $(CommonPath)/cone-e-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0009.cc \
	$(CommonPath)/cone-e-05.cc $(cteq6pdf) \
	-o fnl0009 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0010: LHC CMS incl. jets with fastjet kT D=0.6 ============
# ======================================================================
# LHC scenario for CMS @ 14 TeV
#
fnl0010: links fnl0010.cc $(CommonPath)/fj-kt-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0010.cc \
	$(CommonPath)/fj-kt-06.cc $(cteq6pdf) \
	-o fnl0010 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL0011: LHC CMS incl. jets with fastjet Midpoint R=0.7 ======
# ======================================================================
# LHC scenario for CMS @ 14 TeV
#
fnl0011: links fnl0011.cc $(CommonPath)/fj-mc-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0011.cc \
	$(CommonPath)/fj-mc-07.cc $(cteq6pdf) \
	-o fnl0011 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL1007: LHC CMS incl. jets with kT D=0.6 ====================
# ======================================================================
# LHC scenario for CMS like 0007 @ 10 TeV
#
fnl1007: links fnl1007.cc $(CommonPath)/kt-e-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1007.cc \
	$(CommonPath)/kt-e-06.cc $(cteq6pdf) \
	-o fnl1007 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL1010: LHC CMS incl. jets with fastjet kT D=0.6 ============
# ======================================================================
# LHC scenario for CMS like 0010 @ 10 TeV
#
fnl1010: links fnl1010.cc $(CommonPath)/fj-kt-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1010.cc \
	$(CommonPath)/fj-kt-06.cc $(cteq6pdf) \
	-o fnl1010 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
#
# ======================================================================
# =====         FNR0001    RHIC test   kT D=0.7       =================
# ======================================================================
# first RHIC scenario with fine binning
#
fnr0001: links fnr0001.cc $(CommonPath)/kt2jet-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnr0001.cc \
	$(CommonPath)/kt2jet-e-07.cc $(cteq6pdf) \
	-o fnr0001 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 






# ======================================================================
# =====         nlofast-add             ================================
# ======================================================================
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
