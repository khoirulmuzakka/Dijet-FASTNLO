* --- force emacs to fortran mode -*-fortran-*-
***********************************************************************
*     
*     Define character strings for print out (from v2.0)
*     
***********************************************************************
      CHARACTER*41  CSEP0,DSEP0,LSEP0,SSEP0
      PARAMETER(CSEP0='#########################################')
      PARAMETER(DSEP0='#========================================')
      PARAMETER(LSEP0='#----------------------------------------')
      PARAMETER(SSEP0='*****************************************')
      CHARACTER*82  CSEPS,DSEPS,LSEPS,SSEPS
      PARAMETER(CSEPS='#########################################'//
     >     '#########################################')
      PARAMETER(DSEPS='#========================================'//
     >     '=========================================')
      PARAMETER(LSEPS='#----------------------------------------'//
     >     '-----------------------------------------')
      PARAMETER(SSEPS='*****************************************'//
     >     '*****************************************')
      CHARACTER*164 CSEPL,DSEPL,LSEPL,SSEPL
      PARAMETER(CSEPL='#########################################'//
     >     '#########################################'//
     >     '#########################################'//
     >     '#########################################')
      PARAMETER(DSEPL='========================================='//
     >     '========================================='//
     >     '========================================='//
     >     '=========================================')
      PARAMETER(LSEPL='-----------------------------------------'//
     >     '-----------------------------------------'//
     >     '-----------------------------------------'//
     >     '-----------------------------------------')
      PARAMETER(SSEPL='*****************************************'//
     >     '*****************************************'//
     >     '*****************************************'//
     >     '*****************************************')
C
C ----- define character strings
C
      INTEGER NREAC,NPROC,NALGO,NNLOJET,NTHRCOR
      PARAMETER(NREAC=4,NPROC=4,NALGO=6,NNLOJET=2,NTHRCOR=2)
      CHARACTER*17 CIREACTION(NREAC)
      CHARACTER*17 CIPROC(NPROC)
      CHARACTER*28 CIALGO(NALGO)
      CHARACTER*5  CJETRES1(NALGO)
      CHARACTER*12 CJETRES2(NALGO)
      CHARACTER*72 CNLOJET(NNLOJET)
      CHARACTER*72 CTHRCOR(NTHRCOR)

      DATA CIREACTION/"electron-proton","proton-proton",
     >     "proton-antiproton","photon-proton"/

      DATA CIPROC/"inclusive jets","dijets","three-jets",
     >     "four-jets"/

      DATA CIALGO/"kT","MP","MP_RSEP","SISCone","anti-kT","C/A"/
      DATA CJETRES1/"R","R","R","R","R","R"/
      DATA CJETRES2/"n/a","f_overlap","Rsep","f_overlap","n/a","n/a"/
      
      DATA CNLOJET/
     >     "  Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),",
     >     "  Z. Nagy, Phys. Rev. D68, 094002 (2003)."/

      DATA CTHRCOR/
     >     "  N. Kidonakis and J.F. Owens - please cite",
     >     "  N. Kidonakis, J.F. Owens, Phys. Rev. D63, 054019 (2001)."/
