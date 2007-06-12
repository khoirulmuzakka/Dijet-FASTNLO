#
# ***  create fastNLO author job
#

CommonPath        = ../common

cteq5pdf = $(CommonPath)/Ctq5Par.f \
	   $(CommonPath)/cteq5m.cc

cteq6pdf = $(CommonPath)/Cteq61Pdf.f \
	   $(CommonPath)/cteq6.cc

midpt_e_07 = $(CommonPath)/midpt-e-07.cc \
	     $(CommonPath)/px_interface.f \
	     $(CommonPath)/pxcone_mod_new.f

export CPATH = $(CommonPath)

CommonDir: ./common

./common: 
	ln -s $(CommonPath) ./common


# ======================================================================
# ==========  FNT1002
# ===== Tevatron -   D0 630/1800GeV  single and ratio
# ------------------------------------------------Tevatron Run I 630GeV
fnt1002lo:  CommonDir fnt1002lo.cc $(CommonPath)/cone-et-07.cc $(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1002lo.cc \
	 $(CommonPath)/cone-et-07.cc  $(cteq6pdf) \
	-o fnt1002lo -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
# -----------------------------------------------Tevatron Run I 1800GeV
# Tevatron Run I 1800GeV - to normalize 630GeV results
fnt1002hi:  CommonDir fnt1002hi.cc $(CommonPath)/cone-et-07.cc $(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt1002hi.cc \
	 $(CommonPath)/cone-et-07.cc  $(cteq6pdf) \
	-o fnt1002hi -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3


# ======================================================================
# =========    FNL0001
# ===== LHC - kT algo, D=0.7, E-scheme - 3 rapid. regions
fnl0001:  CommonDir fnl0001.cc $(CommonPath)/kt2jet-e-07.cc $(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnl0001.cc \
	 $(CommonPath)/kt2jet-e-07.cc  $(cteq6pdf) \
	-o fnl0001 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3


# ======================================================================
# =========    FNR0001
# === RHIC - 3 algos at once / all ET scheme:  kT-D10, Cone-D07, Cone-D10 
fnr0001:  CommonDir fnr0001.cc $(CommonPath)/kt2jet-et-10.cc \
	 $(CommonPath)/cone-et-10.cc $(CommonPath)/cone-et-07.cc  $(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnr0001.cc \
	 $(CommonPath)/kt2jet-et-10.cc \
	 $(CommonPath)/cone-et-10.cc $(CommonPath)/cone-et-07.cc  $(cteq6pdf) \
	-o fnr0001 -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3


# ======================================================================
# ==========  FNT2001
# =========   Tevatron Run II 1960GeV - D0  DeltaPhi
# --------------------------------- normalization
fnt2001norm:  CommonDir fnt2001norm.cc $(midpt_e_07) $(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2001norm.cc \
	 $(midpt_e_07)  $(cteq6pdf) \
	-o fnt2001norm -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3
# --------------------------------- differential distribution
fnt2001diff:  CommonDir fnt2001diff.cc $(midpt_e_07) $(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt2001diff.cc \
	 $(midpt_e_07)  $(cteq6pdf) \
	-o fnt2001diff -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3


# ======================================================================
# ==========  FNT200a
# -----------   Tevatron Run II 1960GeV - D0 preliminary incl jets
fnt200a:  CommonDir fnt200a.cc $(CommonPath)/cone-e-07.cc $(cteq6pdf)
	$(NLOJET)/bin/create-nlojet-user fnt200a.cc \
	 $(CommonPath)/cone-e-07.cc  $(cteq6pdf) \
	-o fnt200a -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3


# ======================================================================
# HERA jets
fnhera920incl01a: CommonDir fnhera920incl01a.cc $(CommonPath)/kt-et-10.cc $(cteq5pdf)
	$(NLOJET)/bin/create-nlojet-user fnhera920incl01a.cc \
	$(CommonPath)/kt-et-10.cc $(CommonPath)/alew.f $(cteq5pdf)\
	-o fnhera920incl01a -Wno-deprecated -Wx,-O3,-Wall -Wf,-O3


# ======================================================================
# ======================================================================
nlofast-add: nlofast-add.o $(CommonPath)/gzstream.o
	g++ -o $(NLOJET)/bin/nlofast-add nlofast-add.o $(CommonPath)/gzstream.o -lz


clean:
	rm -rf .libs .obj *~ *.la *.o $(CommonPath)/*.o
