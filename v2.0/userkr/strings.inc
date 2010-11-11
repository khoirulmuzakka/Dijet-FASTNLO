* --- force emacs to fortran mode -*-fortran-*-
C
C ----- define character strings
C
      CHARACTER*72 CHEADER(13)
      CHARACTER*72 CNLOJET(4)
      CHARACTER*72 CTHRCOR(5)
      DATA CHEADER/
     >     "##########################################################",
     >     " ",
     >     "  fastNLO - version 2.0",
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
