#
# ***  create fastNLO author or standalone job 
#

CommonPath        = ../common

cteq5pdf = $(CommonPath)/Ctq5Par.f \
	   $(CommonPath)/cteq5mdis.cc

cteq6pdf = $(CommonPath)/Cteq61Pdf.f \
	   $(CommonPath)/cteq6dis.cc

export CPATH = $(CommonPath)

links: ./common ctq61.00.tbl 

./common:
	ln -s $(CommonPath) ./common

ctq61.00.tbl:
	ln -s $(CommonPath)/ctq61.00.tbl ctq61.00.tbl


#
# ======================================================================
# =====         FNH1001        H1 incl jets @ 300GeV     ==========
# ======================================================================
# 
fnh1001: links fnh1001.cc $(CommonPath)/kt-et-10-dis.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnh1001.cc \
	$(CommonPath)/kt-et-10-dis.cc $(CommonPath)/alew.f $(cteq6pdf) \
	-o fnh1001 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
#
# ======================================================================
# =====         FNH1002        ZEUS incl jets @ 300GeV     ==========
# ======================================================================
# 
fnh1002: links fnh1002.cc $(CommonPath)/kt-et-10-dis.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnh1002.cc \
	$(CommonPath)/kt-et-10-dis.cc $(CommonPath)/alew.f $(cteq6pdf) \
	-o fnh1002 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
#
# ======================================================================
# =====         FNH1003     H1 fwd. incl jets @ 300GeV, lowQ2  ==========
# ======================================================================
# 
fnh1003: links fnh1003.cc $(CommonPath)/kt-et-10-dis.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnh1003.cc \
	$(CommonPath)/kt-et-10-dis.cc $(CommonPath)/alew.f $(cteq6pdf) \
	-o fnh1003 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
# ======================================================================
# =====         FNH1004        H1 dijets         @ 300GeV     ==========
# ======================================================================
# 
fnh1004: links fnh1004.cc $(CommonPath)/kt-et-10-dis.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnh1004.cc \
	$(CommonPath)/kt-et-10-dis.cc $(CommonPath)/alew.f $(cteq6pdf) \
	-o fnh1004 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
# ======================================================================
# =====         FNH1005        H1 fwd+dijets   @ 300GeV     ==========
# ======================================================================
# 
fnh1005: links fnh1005.cc $(CommonPath)/kt-et-10-dis.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnh1005.cc \
	$(CommonPath)/kt-et-10-dis.cc $(CommonPath)/alew.f $(cteq6pdf) \
	-o fnh1005 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
# ======================================================================
# =====         FNH1006        ZEUS fwd+dijets   @ 319GeV     ==========
# ======================================================================
# 
fnh1006: links fnh1006.cc $(CommonPath)/kt-et-10-dis.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnh1006.cc \
	$(CommonPath)/kt-et-10-dis.cc $(CommonPath)/alew.f $(cteq6pdf) \
	-o fnh1006 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
#
# ======================================================================
# =====         FNH2001        H1 incl jets @ 318GeV     ==========
# ======================================================================
# 
fnh2001: links fnh2001.cc $(CommonPath)/kt-et-10-dis.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnh2001.cc \
	$(CommonPath)/kt-et-10-dis.cc $(CommonPath)/alew.f $(cteq6pdf) \
	-o fnh2001 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
#
# ======================================================================
# =====         FNH2002        ZEUS incl jets @ 318GeV     ==========
# ======================================================================
# 
fnh1002: links fnh2002.cc $(CommonPath)/kt-et-10-dis.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnh2002.cc \
	$(CommonPath)/kt-et-10-dis.cc $(CommonPath)/alew.f $(cteq6pdf) \
	-o fnh2002 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
#
# ======================================================================
# =====         FNHXXXX       testing                         ==========
# ======================================================================
# 
fnhxxxx: links fnhxxxx.cc $(CommonPath)/kt-et-10-dis.cc \
	$(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnhxxxx.cc \
	$(CommonPath)/kt-et-10-dis.cc $(CommonPath)/alew.f $(cteq6pdf) \
	-o fnhxxxx -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3

clean:
	rm -rf .libs .obj *~ *.la *.o $(CommonPath)/*.o ctq61.00.tbl ./common
