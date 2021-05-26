* --- force emacs to fortran mode -*-fortran-*-
C ---------------------------------------------------------------------
C
C  Commonblock definition for fastNLO v1.4
C
C   >> this code usually should not be edited
C   >> if you have questions, please contact 
C      the authors at:  fastnlo@hepforge.org
C 
C KR 01/07/11 Final version for v14 fastnloweb code
C ---------------------------------------------------------------------
C - Needed for bin boundaries, xlimit, scales
      INTEGER ISEPARATOR, NRAPMAX,NPTMAX,NBINTOTMAX,NXMAX,NSCALEMAX,
     >     NMAXSUBPROC,NSCALEBINMAX
      PARAMETER(ISEPARATOR=1234567890)
C - Attention! Scenario dependend limits! Will be checked when reading table!
C - Set to LARGE values, but need to be adapted eventually ...
      PARAMETER(NRAPMAX=50,NPTMAX=50,NBINTOTMAX=500,NXMAX=250)
      PARAMETER(NMAXSUBPROC=7,NSCALEBINMAX=4,NSCALEMAX=4)
      
      INTEGER*8 NEVT(6)  
      INTEGER ITABVERSION, IREACTION, IXSECTUNITS, IPROC, IALGO, NORD, 
     >     NXTOT, IXSCHEME,IPDFWGT,IREF, NBINTOT,NDIMENSION,
     >     NRAPIDITY,NPT(NRAPMAX),NSCALEVAR,NSCALEBIN, NXSUM,NSUBPROC
      INTEGER NPOW(6)
      DOUBLE PRECISION ECMS,JETRES1,JETRES2, 
     >     RAPBIN(NRAPMAX+1),PTBIN(NRAPMAX,NPTMAX+1), 
     >     XLIMIT(NRAPMAX,NPTMAX),
     >     MURSCALE(NSCALEMAX),MURVAL(NRAPMAX,NPTMAX,NSCALEBINMAX), 
     >     MUFSCALE(NSCALEMAX),MUFVAL(NRAPMAX,NPTMAX,NSCALEBINMAX),
     >     ARRAY(NBINTOTMAX,NXMAX,NMAXSUBPROC,
     >     1+NSCALEMAX*2,NSCALEBINMAX),
     >     PDF(NBINTOTMAX,NXMAX,NMAXSUBPROC,NSCALEBINMAX), ! PDF Array
     >     RESULT(NBINTOTMAX,NMAXSUBPROC+1,3),
     >     XSECT(NBINTOTMAX,3)  ! Array returned as argument
      
      COMMON/CFNX9999/ECMS,JETRES1,JETRES2, 
     >     RAPBIN,PTBIN,XLIMIT,
     >     MURSCALE,MURVAL,MUFSCALE,MUFVAL,
     >     ARRAY,PDF,RESULT,
     >     NEVT,ITABVERSION,IREACTION,IXSECTUNITS,IPROC,IALGO, 
     >     NORD,NPOW,
     >     NXTOT,IXSCHEME,IPDFWGT,IREF,NBINTOT,NDIMENSION,NRAPIDITY,
     >     NPT,NSCALEVAR,NSCALEBIN,NXSUM,NSUBPROC

      CHARACTER*64 NAMELABEL(5)
      CHARACTER*64 POWLABEL(6)
      CHARACTER*64 CODELABEL(6)
      CHARACTER*24 DIMLABEL(2)
      CHARACTER*64 SCALELABEL

      COMMON/SFNX9999/
     >     NAMELABEL,POWLABEL,CODELABEL,DIMLABEL,SCALELABEL
