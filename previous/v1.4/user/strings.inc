* --- force emacs to fortran mode -*-fortran-*-
C
C ----- define character strings
C
      CHARACTER*17 CIREACTION(4)
      CHARACTER*17 CIPROC(4)
      CHARACTER*28 CIALGO(6)
      CHARACTER*5  CJETRES1(6)
      CHARACTER*12 CJETRES2(6)
      CHARACTER*72 CHEADER(13)
      CHARACTER*72 CNLOJET(4)
      CHARACTER*72 CTHRCOR(5)

      DATA CIREACTION/"electron-proton","proton-proton",
     >     "proton-antiproton","photon-proton"/
      DATA CIPROC/"inclusive jets","dijets","three-jets",
     >     "four-jets"/
      DATA CIALGO/"kT","MP","MP_RSEP","SISCone","anti-kT","C/A"/
      DATA CJETRES1/"R","R","R","R","R","R"/
      DATA CJETRES2/"n/a","f_overlap","Rsep","f_overlap","n/a","n/a"/

      DATA CHEADER/
     >     "##########################################################",
     >     " ",
     >     "  fastNLO - version 1.4         Jan. 31, 2006",
     >     " ",
     >     "  Thomas Kluge, Klaus Rabbertz, Markus Wobisch",
     >     " ",
     >     "  if you use this code, please cite as reference:",
     >     "      T. Kluge, K. Rabbertz and M. Wobisch,",
     >     "      hep-ph/0609285",
     >     " ",
     >     "      ... now reading the coefficient table:",
     >     " ",
     >     "##########################################################"/

      DATA CNLOJET/
     >     "For NLOJET++ please cite:",
     >     "  Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),",
     >     "  Z. Nagy, Phys. Rev. D68, 094002 (2003).",
     >     " "/

      DATA CTHRCOR/
     >     "The 2-loop threshold corrections for the inclusive jet",
     >     "cross section in pp and ppbar have been computed by:",
     >     "  N. Kidonakis and J.F. Owens - please cite",
     >     "  N. Kidonakis, J.F. Owens, Phys. Rev. D63, 054019 (2001).",
     >     " "/
