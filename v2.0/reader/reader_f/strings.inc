* --- force emacs to fortran mode -*-fortran-*-
C
C ----- define character strings
C
      INTEGER NHEADER
      PARAMETER (NHEADER=7)
      CHARACTER*72 CHEADER(NHEADER)
ckr      CHARACTER*72 CNLOJET(4)
ckr      CHARACTER*72 CTHRCOR(5)

      DATA CHEADER/
     >     " fastNLO_reader_2.1.0",
     >     "",
     >     " D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch",
     >     "",
     >     " If you use this code, please cite:",
     >     "   T. Kluge, K. Rabbertz, M. Wobisch, hep-ph/0609285",
     >     "   D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch, TBD"/

Comment:       DATA CNLOJET/
Comment:      >     "For NLOJET++ please cite:",
Comment:      >     "  Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),",
Comment:      >     "  Z. Nagy, Phys. Rev. D68, 094002 (2003).",
Comment:      >     " "/
Comment: 
Comment:       DATA CTHRCOR/
Comment:      >     "The 2-loop threshold corrections for the inclusive jet",
Comment:      >     "cross section in pp and ppbar have been computed by:",
Comment:      >     "  N. Kidonakis and J.F. Owens - please cite",
Comment:      >     "  N. Kidonakis, J.F. Owens, Phys. Rev. D63, 054019 (2001).",
Comment:      >     " "/
