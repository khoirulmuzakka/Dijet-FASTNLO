*     --- force emacs to fortran mode -*-fortran-*-
***********************************************************************
*
*     Common block definition for v14 style binning
*     Needs fnx9999.inc and uncert.inc to be included!
*
***********************************************************************
      INTEGER NPTMAX,NRAPMAX,NSCALEMAX
      PARAMETER(NPTMAX=MXOBSBIN,NRAPMAX=MXOBSBIN,NSCALEMAX=8)
      INTEGER NRAPIDITY,NPT(NRAPMAX),NORD,NSCLS
      DOUBLE PRECISION RAPBIN(NRAPMAX+1),PTBIN(NRAPMAX,NPTMAX+1)
      COMMON/V14/RAPBIN,PTBIN,NRAPIDITY,NPT,NORD,NSCLS
