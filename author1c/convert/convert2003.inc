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
c MW 01/17/06 new variables for tableformat version 1c
c ------------------------------------------------------------
c      INTEGER NWPAWC
c      PARAMETER (NWPAWC=900000)
c      REAL HMEMOR(NWPAWC)

      INTEGER ISEPARATOR, NRAPMAX,NPTMAX,NBINTOTMAX,NXMAX,NSCALEMAX,
     +     NMAXSUBPROC,NSCALEBINMAX
c - needed for bin boundaries, xlimit, scales
c
c - new 'good default' settings (o.k. for everything)
      PARAMETER (NRAPMAX=50,NPTMAX=160,NSCALEMAX=6)
      PARAMETER (ISEPARATOR=1234567890,NBINTOTMAX=900,NXMAX=300)
      PARAMETER (NMAXSUBPROC=7,NSCALEBINMAX=2)

      INTEGER*8 NEVT(6)  
      INTEGER ITABVERSION, IREACTION, IXSECTUNITS, IPROC, IALGO, NORD, 
     +        NXTOT, IXSCHEME,IPDFWGT,IREF, NBINTOT,NDIMENSION,
     +        NRAPIDITY,NPT(NRAPMAX),NSCALEVAR,NSCALEBIN, NXSUM,NSUBPROC
      INTEGER NPOW(6)
      CHARACTER*64 NAMELABEL(5)
      CHARACTER*64 POWLABEL(6)
      CHARACTER*64 CODELABEL(6)
      CHARACTER*24 DIMLABEL(2)
      CHARACTER*64 SCALELABEL
      DOUBLE PRECISION ECMS, JETRES1, JETRES2, 
     +        RAPBIN(NRAPMAX+1), PTBIN(NRAPMAX,NPTMAX+1), 
     +        XLIMIT(NRAPMAX,NPTMAX),
     +        MURSCALE(NSCALEMAX), MURVAL(NRAPMAX,NPTMAX,NSCALEBINMAX), 
     +        MUFSCALE(NSCALEMAX), MUFVAL(NRAPMAX,NPTMAX,NSCALEBINMAX),
     +        array(NBINTOTMAX,NXMAX,NMAXSUBPROC,1+NSCALEMAX*2,NSCALEBINMAX), ! 
*(Nord-1)
     +        pdf(NBINTOTMAX,NXMAX,NMAXSUBPROC,NSCALEBINMAX),    ! array for PDFs
     +        result(NBINTOTMAX,NMAXSUBPROC+1,3),
     +        xsect(NBINTOTMAX,3)    ! array returned as argument
c -> remove multiple scales in arrays
c     +        result(NBINTOTMAX,NMAXSUBPROC+1,3,NSCALEMAX),
c     +        xsect(NBINTOTMAX,3,NSCALEMAX)    ! array returned as argument

      COMMON /CFNX9999/ ECMS, JETRES1, JETRES2, 
     +        RAPBIN, PTBIN,  XLIMIT,
     +        MURSCALE, MURVAL, MUFSCALE, MUFVAL,
     +        array, pdf, result,
     +        NEVT, ITABVERSION, IREACTION, IXSECTUNITS, IPROC, IALGO, 
     +        NORD, NPOW,
     +        NAMELABEL,POWLABEL,CODELABEL,DIMLABEL,SCALELABEL, 
     +        NXTOT, IXSCHEME,IPDFWGT,IREF,NBINTOT,NDIMENSION,NRAPIDITY,
     +        NPT,NSCALEVAR,NSCALEBIN, NXSUM,NSUBPROC

