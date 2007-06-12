* --- force emacs to fortran mode -*-fortran-*-
c ------------------------------------------------------------
c
c  Commonblock definition for fastNLO
c 
c TK 11/11/05 Changed dimensions of nevt and npow to 6
c kr 30/08/05 Added trivial comment to start emacs fortran mode 
c MW 08/26/05 restructure "array" - remove two dimensions
c TK 12/07/05 NPOW added
c MW 04/15/05
c ------------------------------------------------------------
c      INTEGER NWPAWC
c      PARAMETER (NWPAWC=900000)
c      REAL HMEMOR(NWPAWC)

      INTEGER ISEPARATOR, NRAPMAX,NPTMAX,NRAPPT,NXMAX,NSCALEMAX
c - needed for bin boundaries, xlimit, scales
      PARAMETER (NRAPMAX=6,NPTMAX=33,NSCALEMAX=6)
      PARAMETER (ISEPARATOR=1234567890,NRAPPT=123,NXMAX=1275)
c      PARAMETER (ISEPARATOR=1234567890,NRAPMAX=6,NPTMAX=75,NXMAX=50)

      INTEGER*8 NEVT(6)
      INTEGER IREACTION, IPROC, IALGO, NORD, 
     +        NXTOT, IXSCHEME,
     +        IPDFWGT,NRAPIDITY,NPT(NRAPMAX),NSCALEVAR
      INTEGER NPOW(6)
      DOUBLE PRECISION ECMS, JETRES1, JETRES2, 
     +        RAPBIN(NRAPMAX+1), PTBIN(NRAPMAX,NPTMAX+1), 
     +        XLIMIT(NRAPMAX,NPTMAX),
     +        MURSCALE(NSCALEMAX), MURVAL(NRAPMAX,NPTMAX), 
     +        MUFSCALE(NSCALEMAX), MUFVAL(NRAPMAX,NPTMAX),
     +        array(NRAPPT,NXMAX,7,1+NSCALEMAX),
     +        pdf(NRAPPT,NXMAX,3),    ! new array for PDFs
     +        result(NRAPPT,4,2,NSCALEMAX),
c     +        result(NRAPMAX*NPTMAX,8,6)
     +        xsect(NRAPPT,2,NSCALEMAX)    ! new array directly returned as argument

      COMMON /CFASTNLO/ ECMS, JETRES1, JETRES2, 
     +        RAPBIN, PTBIN,  XLIMIT,
     +        MURSCALE, MURVAL, MUFSCALE, MUFVAL,
     +        array, pdf, result,
     +        NEVT, IREACTION, IPROC, IALGO, NORD, NPOW,
     +        NXTOT, IXSCHEME,
     +        IPDFWGT,NRAPIDITY,NPT,NSCALEVAR


      
c
c ----- define character strings
c
      CHARACTER*17 CIREACTION(3)
      CHARACTER*17 CIPROC(4)
      CHARACTER*28 CIALGO(4)
      CHARACTER*5 CJETRES1(4)
      CHARACTER*12 CJETRES2(4)

      DATA CIREACTION/'electron-proton  ','proton-proton    ',
     +     'proton-antiproton'/
      DATA CIPROC/'inclusive jets   ','dijets','three-jets','four-jets'/
      DATA CIALGO/'kT algorithm','midpoint cone algorithm',
     +     'rsep cone algorithm','searchcone algorithm'/
      DATA CJETRES1/'D','Rcone','Rcone','Rcone'/
      DATA CJETRES2/'n/a','f_overlap','Rsep','R_searchcone'/
