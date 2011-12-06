* --- force emacs to fortran mode -*-fortran-*-
C
C ----- define character strings
C
      INTEGER NHEADER
      PARAMETER (NHEADER=8)
      CHARACTER*72 CHEADER(NHEADER)

      DATA CHEADER/
     >     " fastNLO_reader_2.1.0",
     >     "",
     >     " D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch",
     >     "",
     >     " If you use this code, please cite:",
     >     "   T. Kluge, K. Rabbertz, M. Wobisch, hep-ph/0609285",
     >     "   D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch,",
     >     "      arXiv:nnnn.mmmm"/
