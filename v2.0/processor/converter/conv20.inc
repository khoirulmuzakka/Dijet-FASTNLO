* --- force emacs to fortran mode -*-fortran-*-
c ------------------------------------------------------------
c
c  Commonblock definition for fastNLO v1.5
c
c   >> this code should not be edited
c   >> if you have questions, please contact 
c      the authors at:  fastnlo@cedar.ac.uk
c 
c MW 11/20/06 implement tableformat version 1.5 
c ------------------------------------------------------------

      INTEGER ISEPARATOR, NRAPMAX,NPTMAX,NBINTOTMAX,NXMAX,NSCALEMAX,
     +     NMAXSUBPROC,NSCALEBINMAX
c - needed for bin boundaries, xlimit, scales
c
c - scenario-dependent limits
      PARAMETER (ISEPARATOR=1234567890)
      PARAMETER (NRAPMAX=50,NPTMAX=160,NBINTOTMAX=800,NXMAX=200) ! 50xBins
      PARAMETER (NMAXSUBPROC=7,NSCALEBINMAX=4,NSCALEMAX=6)

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
     +        array(NBINTOTMAX,NXMAX,NMAXSUBPROC,1+NSCALEMAX*2,NSCALEBINMAX),
     +        pdf(NBINTOTMAX,NXMAX,NMAXSUBPROC,NSCALEBINMAX), ! PDF array
     +        result(NBINTOTMAX,NMAXSUBPROC+1,3),
     +        xsect(NBINTOTMAX,3)    ! array returned as argument

      COMMON /CFNX9999/ ECMS, JETRES1, JETRES2, 
     +        RAPBIN, PTBIN,  XLIMIT,
     +        MURSCALE, MURVAL, MUFSCALE, MUFVAL,
     +        array, pdf, result,
     +        NEVT, ITABVERSION, IREACTION, IXSECTUNITS, IPROC, IALGO, 
     +        NORD, NPOW,
     +        NAMELABEL,POWLABEL,CODELABEL,DIMLABEL,SCALELABEL, 
     +        NXTOT, IXSCHEME,IPDFWGT,IREF,NBINTOT,NDIMENSION,NRAPIDITY,
     +        NPT,NSCALEVAR,NSCALEBIN, NXSUM,NSUBPROC
      
c
c ----- define character strings
c
      CHARACTER*17 CIREACTION(4)
      CHARACTER*17 CIPROC(4)
      CHARACTER*28 CIALGO(4)
      CHARACTER*5 CJETRES1(4)
      CHARACTER*12 CJETRES2(4)

      DATA CIREACTION/'electron-proton  ','proton-proton    ',
     +     'proton-antiproton','photon-proton'/
      DATA CIPROC/'inclusive jets   ','dijets','three-jets','four-jets'/
      DATA CIALGO/'kT algorithm','midpoint cone algorithm',
     +     'rsep cone algorithm','searchcone algorithm'/
      DATA CJETRES1/'D','Rcone','Rcone','Rcone'/
      DATA CJETRES2/'n/a','f_overlap','Rsep','R_searchcone'/

      CHARACTER*64 CHEADER(22)
      DATA CHEADER/
     +     '################################################################',
     +     ' ',
     +     '  fastNLO - version 1.5         May 30, 2007',
     +     ' ',
     +     '  Thomas Kluge, Klaus Rabbertz, Markus Wobisch',
     +     ' ',
     +     '  if you use this code, please cite as reference:',
     +     '      T. Kluge, K. Rabbertz and M. Wobisch,',
     +     '         publication in preparation,',
     +     '         (hep-ex No. to be added)',
     +     ' ',
     +     '     .... now reading the coefficient table:',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' '/

      CHARACTER*64 CNLOJET(4)
      DATA CNLOJET/
     +     '  for NLOJET++ please cite:',
     +     '      Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),',
     +     '      Z. Nagy, Phys. Rev. D68, 094002 (2003).',
     +     ' '/

      CHARACTER*64 CTHRCOR(5)
      DATA CTHRCOR/
     +     '  the 2-loop threshold corrections for the inclusive jet',
     +     '  cross section in pp and ppbar have been computed by:',
     +     '  N. Kidonakis and J.F. Owens - please cite',
     +     '     N. Kidonakis, J.F. Owens, Phys. Rev. D63, 054019 (2001).',
     +     ' '/
