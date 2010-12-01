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

# --- ... as original version with 12 x bins, but different x binning
fnt2003.x_06_2: links fnt2003.x_06_2.cc $(CommonPath)/kt-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2003.x_06_2.cc \
	$(CommonPath)/kt-e-07.cc $(cteq6pdf) \
	-o fnt2003.x_06_2 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
fnt2003.x_24_2: links fnt2003.x_24_2.cc $(CommonPath)/kt-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2003.x_24_2.cc \
	$(CommonPath)/kt-e-07.cc $(cteq6pdf) \
	-o fnt2003.x_24_2 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
fnt2003.x_48_2: links fnt2003.x_48_2.cc $(CommonPath)/kt-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2003.x_48_2.cc \
	$(CommonPath)/kt-e-07.cc $(cteq6pdf) \
	-o fnt2003.x_48_2 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
fnt2003.x_12_1: links fnt2003.x_12_1.cc $(CommonPath)/kt-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2003.x_12_1.cc \
	$(CommonPath)/kt-e-07.cc $(cteq6pdf) \
	-o fnt2003.x_12_1 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3



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

# ======================================================================
# =====   FNT2003: Tevatron Run II CDF incl. jets with kT D=0.7   ======
# ======================================================================
# hep-ex/0512062
#
#fnt2003: links fnt2003.cc $(CommonPath)/kt-e-07.cc \
#	$(cteq6pdf)
#	$(NLOJET)/bin/create-nlojet-user fnt2003.cc \
#	$(CommonPath)/kt-e-07.cc $(cteq6pdf) \
#	-o fnt2003 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNT2005: Tevatron Run II CDF incl. jets with kT D=0.5   ======
# ======================================================================
# hep-ex/0701051
#
fnt2005: links fnt2005.cc $(CommonPath)/kt-e-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2005.cc \
	$(CommonPath)/kt-e-05.cc $(cteq6pdf) \
	-o fnt2005 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNT2006: Tevatron Run II CDF incl. jets with kT D=1.0   ======
# ======================================================================
# hep-ex/0701051
#
fnt2006: links fnt2006.cc $(CommonPath)/kt-e-10.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2006.cc \
	$(CommonPath)/kt-e-10.cc $(cteq6pdf) \
	-o fnt2006 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNT2310: CMS incl. jets with fastjet kT D=0.6 ============
# ======================================================================
# LHC scenario for CMS, pt binning acc. to expected resolution
# For comparison done at Tevatron Run II energies in p\bar{p} collisions
#
fnt2310: links fnt2310.cc $(CommonPath)/fj-kt-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2310.cc \
	$(CommonPath)/fj-kt-06.cc $(cteq6pdf) \
	-o fnt2310 -Wno-deprecated -Wf,-O3 \
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
	-o fnl0005 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

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
# =====   FNL0007: LHC CMS incl. jets with kT R=0.6 ====================
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
fnl0009: links fnl0009.cc $(CommonPath)/cone-e-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0009.cc \
	$(CommonPath)/cone-e-07.cc $(cteq6pdf) \
	-o fnl0009 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0010: LHC CMS incl. jets with fastjet kT R=0.6 ============
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
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/CDFCones \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/CDFCones,-lCDFConesPlugin,-L$(FASTJET)/plugins/CDFCones/CDFcode,-lCDFcode

# ======================================================================
# =====   FNL0017: LHC CMS incl. jets with kT R=0.4 ====================
# ======================================================================
# LHC scenario for CMS @ 14 TeV
#
fnl0017: links fnl0017.cc $(CommonPath)/kt-e-04.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0017.cc \
	$(CommonPath)/kt-e-04.cc $(cteq6pdf) \
	-o fnl0017 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0018: LHC CMS incl. jets with fastjet SISCone R=0.5 =======
# ======================================================================
# LHC scenario for CMS @ 14 TeV
#
fnl0018: links fnl0018.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0018.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl0018 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL0019: LHC CMS incl. jets with Midpoint R=0.5 ==============
# ======================================================================
# LHC scenario for CMS @ 14 TeV
#
fnl0019: links fnl0019.cc $(CommonPath)/cone-e-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0019.cc \
	$(CommonPath)/cone-e-05.cc $(cteq6pdf) \
	-o fnl0019 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0020: LHC CMS incl. jets with fastjet kT R=0.4 ============
# ======================================================================
# LHC scenario for CMS @ 14 TeV
#
fnl0020: links fnl0020.cc $(CommonPath)/fj-kt-04.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0020.cc \
	$(CommonPath)/fj-kt-04.cc $(cteq6pdf) \
	-o fnl0020 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL0021: LHC CMS incl. jets with fastjet Midpoint R=0.5 ======
# ======================================================================
# LHC scenario for CMS @ 14 TeV
#
fnl0021: links fnl0021.cc $(CommonPath)/fj-mc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0021.cc \
	$(CommonPath)/fj-mc-05.cc $(cteq6pdf) \
	-o fnl0021 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/CDFCones \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/CDFCones,-lCDFConesPlugin,-L$(FASTJET)/plugins/CDFCones/CDFcode,-lCDFcode

# ======================================================================
# =====   FNL0117: LHC CMS incl. jets with kT R=0.4 ====================
# ======================================================================
# LHC scenario for CMS like 0017 but for forward jets
#
fnl0117: links fnl0117.cc $(CommonPath)/kt-e-04.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0117.cc \
	$(CommonPath)/kt-e-04.cc $(cteq6pdf) \
	-o fnl0117 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0118: LHC CMS incl. jets with fastjet SISCone R=0.5 =======
# ======================================================================
# LHC scenario for CMS like 1118 but for 14 TeV
#
fnl0118: links fnl0118.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0118.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl0118 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl0118.x_06_2: links fnl0118.x_06_2.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0118.x_06_2.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl0118.x_06_2 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl0118.x_24_2: links fnl0118.x_24_2.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0118.x_24_2.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl0118.x_24_2 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl0118.x_48_2: links fnl0118.x_48_2.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0118.x_48_2.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl0118.x_48_2 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl0118.x_12_1: links fnl0118.x_12_1.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0118.x_12_1.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl0118.x_12_1 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL0217: LHC CMS incl. jets with kT R=0.4 ====================
# ======================================================================
# LHC scenario for CMS like 0117 but for forward jets with pseudo-rapidity
#
fnl0217: links fnl0217.cc $(CommonPath)/kt-e-04.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0217.cc \
	$(CommonPath)/kt-e-04.cc $(cteq6pdf) \
	-o fnl0217 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL0218: LHC CMS incl. jets with fastjet SISCone R=0.5 =======
# ======================================================================
# LHC scenario for CMS like 1218 but for 14 TeV
#
fnl0218: links fnl0218.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0218.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl0218 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL0222: LHC CMS forward jets with fastjet anti-kT R=0.5 =====
# ======================================================================
# LHC scenario for CMS at 14 TeV
#
fnl0222: links fnl0222.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0222.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl0222 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL0310: LHC CMS incl. jets with fastjet kT R=0.6 ============
# ======================================================================
# LHC scenario for CMS, pt binning acc. to expected resolution, @ 14 TeV 
#
fnl0310: links fnl0310.cc $(CommonPath)/fj-kt-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0310.cc \
	$(CommonPath)/fj-kt-06.cc $(cteq6pdf) \
	-o fnl0310 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL0408: LHC CMS dijet mass with fastjet SISCone R=0.7 =======
# ======================================================================
# LHC scenario for CMS, @ 14 TeV
#
fnl0408: links fnl0408.cc $(CommonPath)/fj-sc-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0408.cc \
	$(CommonPath)/fj-sc-07.cc $(cteq6pdf) \
	-o fnl0408 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL0602: LHC ATLAS dchi dijets with fastjet kT R=1.0 ============
# ======================================================================
# LHC scenario for ATLAS, @ 14 TeV 
#
fnl0602: links fnl0602.cc $(CommonPath)/fj-kt-10.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0602.cc \
	$(CommonPath)/fj-kt-10.cc $(cteq6pdf) \
	-o fnl0602 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL1007: LHC CMS incl. jets with kT R=0.6 ====================
# ======================================================================
# LHC scenario for CMS like 0007 @ 10 TeV
#
fnl1007: links fnl1007.cc $(CommonPath)/kt-e-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1007.cc \
	$(CommonPath)/kt-e-06.cc $(cteq6pdf) \
	-o fnl1007 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnl1007.y_1_x_06_2: links fnl1007.y_1_x_06_2.cc $(CommonPath)/kt-e-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1007.y_1_x_06_2.cc \
	$(CommonPath)/kt-e-06.cc $(cteq6pdf) \
	-o fnl1007.y_1_x_06_2 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnl1007.y_1_x_12_2: links fnl1007.y_1_x_12_2.cc $(CommonPath)/kt-e-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1007.y_1_x_12_2.cc \
	$(CommonPath)/kt-e-06.cc $(cteq6pdf) \
	-o fnl1007.y_1_x_12_2 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnl1007.y_1_x_24_2: links fnl1007.y_1_x_24_2.cc $(CommonPath)/kt-e-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1007.y_1_x_24_2.cc \
	$(CommonPath)/kt-e-06.cc $(cteq6pdf) \
	-o fnl1007.y_1_x_24_2 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnl1007.y_1_x_48_2: links fnl1007.y_1_x_48_2.cc $(CommonPath)/kt-e-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1007.y_1_x_48_2.cc \
	$(CommonPath)/kt-e-06.cc $(cteq6pdf) \
	-o fnl1007.y_1_x_48_2 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnl1007.y_1_x_12_1: links fnl1007.y_1_x_12_1.cc $(CommonPath)/kt-e-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1007.y_1_x_12_1.cc \
	$(CommonPath)/kt-e-06.cc $(cteq6pdf) \
	-o fnl1007.y_1_x_12_1 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

fnl1007.y_1_x_48_1: links fnl1007.y_1_x_48_1.cc $(CommonPath)/kt-e-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1007.y_1_x_48_1.cc \
	$(CommonPath)/kt-e-06.cc $(cteq6pdf) \
	-o fnl1007.y_1_x_48_1 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

# ======================================================================
# =====   FNL1008: LHC CMS incl. jets with fastjet SISCone R=0.7 =======
# ======================================================================
# LHC scenario for CMS like 0008 @ 10 TeV
#
fnl1008: links fnl1008.cc $(CommonPath)/fj-sc-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1008.cc \
	$(CommonPath)/fj-sc-07.cc $(cteq6pdf) \
	-o fnl1008 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL1010: LHC CMS incl. jets with fastjet kT R=0.6 ============
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

# ======================================================================
# =====   FNL1018: LHC CMS incl. jets with fastjet SISCone R=0.5 =======
# ======================================================================
# LHC scenario for CMS like 0018 @ 10 TeV
#
fnl1018: links fnl1018.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1018.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1018 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL1118: LHC CMS incl. jets with fastjet SISCone R=0.5 =======
# ======================================================================
# LHC scenario for CMS like 1018 but for forward jets |y| = 3 - 5
#
fnl1118: links fnl1118.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1118.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1118 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1118.x_06_2: links fnl1118.x_06_2.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1118.x_06_2.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1118.x_06_2 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1118.x_12_2: links fnl1118.x_12_2.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1118.x_12_2.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1118.x_12_2 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1118.x_24_2: links fnl1118.x_24_2.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1118.x_24_2.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1118.x_24_2 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1118.x_48_2: links fnl1118.x_48_2.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1118.x_48_2.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1118.x_48_2 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1118.x_12_1: links fnl1118.x_12_1.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1118.x_12_1.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1118.x_12_1 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1118.x_48_1: links fnl1118.x_48_1.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1118.x_48_1.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1118.x_48_1 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL1218: LHC CMS incl. jets with fastjet SISCone R=0.5 =======
# ======================================================================
# LHC scenario for CMS like 1118 but using pseudo-rapidity
#
fnl1218: links fnl1218.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1218.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1218 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL1222: LHC CMS forward jets with fastjet anti-kT R=0.5 =====
# ======================================================================
# LHC scenario for CMS at 10 TeV
#
fnl1222: links fnl1222.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1222.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl1222 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL1308: LHC CMS incl. jets with fastjet SISCone R=0.7 =======
# ======================================================================
# LHC scenario for CMS,  pt binning acc. to expected resolution
#
fnl1308: links fnl1308.cc $(CommonPath)/fj-sc-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1308.cc \
	$(CommonPath)/fj-sc-07.cc $(cteq6pdf) \
	-o fnl1308 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL1310: LHC CMS incl. jets with fastjet kT R=0.6 ============
# ======================================================================
# LHC scenario for CMS, pt binning acc. to expected resolution
#
fnl1310: links fnl1310.cc $(CommonPath)/fj-kt-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1310.cc \
	$(CommonPath)/fj-kt-06.cc $(cteq6pdf) \
	-o fnl1310 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL1408: LHC CMS Di-jet mass with fastjet SISCone R=0.7 ======
# ======================================================================
# LHC scenario for CMS
#
fnl1408: links fnl1408.cc $(CommonPath)/fj-sc-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1408.cc \
	$(CommonPath)/fj-sc-07.cc $(cteq6pdf) \
	-o fnl1408 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL1518: LHC CMS Az. decorrelation with fastjet SISCone R=0.5 
# ======================================================================
# LHC scenario for CMS
#
fnl1518norm: links fnl1518norm.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1518norm.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1518norm -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1518diff: links fnl1518diff.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1518diff.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1518diff -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1518diffpt1: links fnl1518diffpt1.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1518diffpt1.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1518diffpt1 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1518diffpt2: links fnl1518diffpt2.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1518diffpt2.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1518diffpt2 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1518diffpt3: links fnl1518diffpt3.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1518diffpt3.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1518diffpt3 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1518diffpt4: links fnl1518diffpt4.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1518diffpt4.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1518diffpt4 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1518diffpt5: links fnl1518diffpt5.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1518diffpt5.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1518diffpt5 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl1518diffpt6: links fnl1518diffpt6.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl1518diffpt6.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl1518diffpt6 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL2218: LHC CMS forward jets with fastjet SISCone R=0.5 =====
# ======================================================================
# LHC scenario for CMS like 1218 but for 7 TeV
#
fnl2218: links fnl2218.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2218.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl2218 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL2222: LHC CMS forward jets with fastjet anti-kT R=0.5 =====
# ======================================================================
# LHC scenario for CMS at 7 TeV
#
fnl2222: links fnl2222.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2222.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2222 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL2222a: LHC CMS forward jets with fastjet anti-kT R=0.5 ====
# ======================================================================
# LHC scenario for CMS at 7 TeV
#
fnl2222a: links fnl2222a.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2222a.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2222a -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

fnl2222amu8: links fnl2222amu8.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2222amu8.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2222amu8 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL2308: LHC CMS incl. jets with fastjet SISCone R=0.7 =======
# ======================================================================
# LHC scenario for CMS,  pt binning acc. to expected resolution, @ 7 TeV
#
fnl2308: links fnl2308.cc $(CommonPath)/fj-sc-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2308.cc \
	$(CommonPath)/fj-sc-07.cc $(cteq6pdf) \
	-o fnl2308 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL2310: LHC CMS incl. jets with fastjet kT R=0.6 ============
# ======================================================================
# LHC scenario for CMS, pt binning acc. to expected resolution, @ 7 TeV 
#
fnl2310: links fnl2310.cc $(CommonPath)/fj-kt-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2310.cc \
	$(CommonPath)/fj-kt-06.cc $(cteq6pdf) \
	-o fnl2310 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL2318: LHC CMS incl. jets with fastjet SISCone R=0.5 =======
# ======================================================================
# LHC scenario for CMS,  pt binning acc. to expected resolution, @ 7 TeV
#
fnl2318: links fnl2318.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2318.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl2318 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL2320: LHC CMS incl. jets with fastjet kT R=0.5 ============
# ======================================================================
# LHC scenario for CMS, pt binning acc. to expected resolution, @ 7 TeV 
#
fnl2320: links fnl2320.cc $(CommonPath)/fj-kt-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2320.cc \
	$(CommonPath)/fj-kt-05.cc $(cteq6pdf) \
	-o fnl2320 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL2322: LHC CMS incl. jets with fastjet anti-kT R=0.5 =======
# ======================================================================
# LHC scenario for CMS, pt binning acc. to expected resolution, @ 7 TeV 
#
fnl2322: links fnl2322.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2322.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2322 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL2323: LHC CMS incl. jets with fastjet CA R=0.5 ============
# ======================================================================
# LHC scenario for CMS, pt binning acc. to expected resolution, @ 7 TeV 
#
fnl2323: links fnl2323.cc $(CommonPath)/fj-ca-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2323.cc \
	$(CommonPath)/fj-ca-05.cc $(cteq6pdf) \
	-o fnl2323 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL2332: LHC CMS incl. jets with fastjet anti-kT R=0.7 =======
# ======================================================================
# LHC scenario for CMS, pt binning acc. to expected resolution, @ 7 TeV 
#
fnl2332a: links fnl2332a.cc $(CommonPath)/fj-ak-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2332a.cc \
	$(CommonPath)/fj-ak-07.cc $(cteq6pdf) \
	-o fnl2332a -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL2342: LHC CMS incl. jets with fastjet anti-kT R=0.5 =======
# ======================================================================
# LHC scenario for CMS, pt binning acc. to expected resolution, @ 7 TeV 
#
fnl2342: links fnl2342.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2342.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2342 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

fnl2342a: links fnl2342a.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2342a.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2342a -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

fnl2342arap1: links fnl2342arap1.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2342arap1.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2342arap1 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL2408: LHC CMS dijet mass with fastjet SISCone R=0.7 =======
# ======================================================================
# LHC scenario for CMS, @ 7 TeV
#
fnl2408: links fnl2408.cc $(CommonPath)/fj-sc-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2408.cc \
	$(CommonPath)/fj-sc-07.cc $(cteq6pdf) \
	-o fnl2408 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL2412: LHC CMS dijet mass with fastjet anti-kT R=0.7 =======
# ======================================================================
# LHC scenario for CMS, @ 7 TeV
#
fnl2412: links fnl2412.cc $(CommonPath)/fj-ak-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2412.cc \
	$(CommonPath)/fj-ak-07.cc $(cteq6pdf) \
	-o fnl2412 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2412a: links fnl2412a.cc $(CommonPath)/fj-ak-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2412a.cc \
	$(CommonPath)/fj-ak-07.cc $(cteq6pdf) \
	-o fnl2412a -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

# ======================================================================
# =====   FNL2442: LHC CMS dijet mass ratio with fastjet anti-kT R=0.7 =
# ======================================================================
# LHC scenario for CMS, @ 7 TeV
#
fnl2442: links fnl2442.cc $(CommonPath)/fj-ak-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2442.cc \
	$(CommonPath)/fj-ak-07.cc $(cteq6pdf) \
	-o fnl2442 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2442x40: links fnl2442x40.cc $(CommonPath)/fj-ak-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2442x40.cc \
	$(CommonPath)/fj-ak-07.cc $(cteq6pdf) \
	-o fnl2442x40 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2442a: links fnl2442a.cc $(CommonPath)/fj-ak-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2442a.cc \
	$(CommonPath)/fj-ak-07.cc $(cteq6pdf) \
	-o fnl2442a -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

# ======================================================================
# =====   FNL2518: LHC CMS Az. decorrelation with fastjet SISCone R=0.5 
# ======================================================================
# LHC scenario for CMS, @ 7 TeV
#
fnl2518norm: links fnl2518norm.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2518norm.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl2518norm -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl2518diff: links fnl2518diff.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2518diff.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl2518diff -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl2518diffpt1: links fnl2518diffpt1.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2518diffpt1.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl2518diffpt1 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl2518diffpt2: links fnl2518diffpt2.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2518diffpt2.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl2518diffpt2 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl2518diffpt3: links fnl2518diffpt3.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2518diffpt3.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl2518diffpt3 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl2518diffpt4: links fnl2518diffpt4.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2518diffpt4.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl2518diffpt4 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl2518diffpt5: links fnl2518diffpt5.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2518diffpt5.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl2518diffpt5 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

fnl2518diffpt6: links fnl2518diffpt6.cc $(CommonPath)/fj-sc-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2518diffpt6.cc \
	$(CommonPath)/fj-sc-05.cc $(cteq6pdf) \
	-o fnl2518diffpt6 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

# ======================================================================
# =====   FNL2522: LHC CMS Az. decorrelation with fastjet anti-kT R=0.5 
# ======================================================================
# LHC scenario for CMS, @ 7 TeV
#
fnl2522norm: links fnl2522norm.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2522norm.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2522norm -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2522norma: links fnl2522norma.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2522norma.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2522norma -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2522normb: links fnl2522normb.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2522normb.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2522normb -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2522diffpt1: links fnl2522diffpt1.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2522diffpt1.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2522diffpt1 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2522diffpt2: links fnl2522diffpt2.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2522diffpt2.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2522diffpt2 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2522diffpt3: links fnl2522diffpt3.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2522diffpt3.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2522diffpt3 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2522diffpt4: links fnl2522diffpt4.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2522diffpt4.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2522diffpt4 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2522diffpt5: links fnl2522diffpt5.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2522diffpt5.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2522diffpt5 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2522diffpt5b: links fnl2522diffpt5b.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2522diffpt5b.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2522diffpt5b -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

# ======================================================================
# =====   FNL2612: LHC CMS dchi dijets with fastjet ak R=0.7 ===========
# ======================================================================
# LHC scenario for CMS, @ 7 TeV 
#
# Like fnl2622d with ak07
fnl2612d: links fnl2612d.cc $(CommonPath)/fj-ak-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2612d.cc \
	$(CommonPath)/fj-ak-07.cc $(cteq6pdf) \
	-o fnl2612d -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL2622: LHC CMS dchi dijets with fastjet ak R=0.5 ===========
# ======================================================================
# LHC scenario for CMS, @ 7 TeV 
#
fnl2622a: links fnl2622a.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2622a.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2622a -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

fnl2622b: links fnl2622b.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2622b.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2622b -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# Like fnl2622a, dijet mass bins changed 
fnl2622c: links fnl2622c.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2622c.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2622c -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# Like fnl2622b, again dijet mass bins changed 
fnl2622d: links fnl2622d.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2622d.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2622d -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL2722: LHC CMS 3-jet ratio with fastjet anti-kT R=0.5 =
# ======================================================================
# LHC scenario for CMS, @ 7 TeV
#
fnl2722norm: links fnl2722norm.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2722norm.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2722norm -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2722: links fnl2722.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2722.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2722 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2722num: links fnl2722num.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2722num.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2722num -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

fnl2722den: links fnl2722den.cc $(CommonPath)/fj-ak-05.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl2722den.cc \
	$(CommonPath)/fj-ak-05.cc $(cteq6pdf) \
	-o fnl2722den -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet

# ======================================================================
# =====   FNL3310: LHC CMS incl. jets with fastjet kT R=0.6 ============
# ======================================================================
# LHC scenario for CMS, pt binning acc. to expected resolution, @ 2.36 TeV 
#
fnl3310: links fnl3310.cc $(CommonPath)/fj-kt-06.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl3310.cc \
	$(CommonPath)/fj-kt-06.cc $(cteq6pdf) \
	-o fnl3310 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include \
	-Wl,-L$(FASTJET)/lib,-lfastjet 

# ======================================================================
# =====   FNL3408: LHC CMS dijet mass with fastjet SISCone R=0.7 =======
# ======================================================================
# LHC scenario for CMS, @ 2.36 TeV
#
fnl3408: links fnl3408.cc $(CommonPath)/fj-sc-07.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl3408.cc \
	$(CommonPath)/fj-sc-07.cc $(cteq6pdf) \
	-o fnl3408 -Wno-deprecated -Wf,-O3 \
	-Wx,-O3,-Wall,-I,$(FASTJET)/include,-I,$(FASTJET)/plugins/SISCone \
	-Wl,-L$(FASTJET)/lib,-lfastjet,-L$(FASTJET)/plugins/SISCone,-lSISConePlugin,-L$(FASTJET)/plugins/SISCone/siscone/siscone,-lsiscone

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
	rm -rf .libs .obj *~ *.la *.lo *.o $(CommonPath)/*.o $(CommonPath)/*.lo ctq61.00.tbl ./common
