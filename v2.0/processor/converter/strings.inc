*     --- force emacs to fortran mode -*-fortran-*-
************************************************************************
*     
*     Define character strings used in fastNLO v1.4
*     
************************************************************************
      CHARACTER*17 CIREACTION(4)
      CHARACTER*17 CIPROC(4)
      CHARACTER*28 CIALGO(6)
      CHARACTER*5  CJETRES1(6)
      CHARACTER*12 CJETRES2(6)
      CHARACTER*64 CNLOJET(3)
      CHARACTER*64 CTHRCOR(3)
      
      DATA CIREACTION/"electron-proton","proton-proton",
     >     "proton-antiproton","photon-proton"/
      DATA CIPROC/"inclusive jets","dijets","three-jets",
     >     "four-jets"/
      DATA CIALGO/"kT","MP","MP_RSEP","SISCone","anti-kT","C/A"/
      DATA CJETRES1/"R","R","R","R","R","R"/
      DATA CJETRES2/"n/a","f_overlap","Rsep","f_overlap","n/a","n/a"/
      
      DATA CNLOJET/
     >     "NLOJet++_2.0.1",
     >     "Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),",
     >     "Z. Nagy, Phys. Rev. D68, 094002 (2003)."/
      DATA CTHRCOR/
     >     "2-loop threshold corrections for the inclusive ",
     >     "jet cross section in pp and ppbar",
     >     "N. Kidonakis, J.F. Owens, Phys. Rev. D63, 054019 (2001)."/
