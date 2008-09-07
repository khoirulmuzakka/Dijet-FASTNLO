* --- force emacs to fortran mode -*-fortran-*-
C ---------------------------------------------------------------------
C
C  Commonblock definition for fastNLO
C
C   >> this code should not be edited
C   >> if you have questions, please contact 
C      the authors at:  fastnlo@cedar.ac.uk
C 
C TK 11/11/05 Changed dimensions of nevt and npow to 6
C KR 30/08/05 Added trivial comment to start emacs fortran mode 
C MW 08/26/05 restructure "array" - remove two dimensions
C TK 12/07/05 NPOW added
C MW 04/15/05
C MW 01/17/06 new variables for tableformat version 1c
C MW 02/01/06 implement tableformat version 1.4 
C KR 01/06/08 on request apply some f77 rules
C ---------------------------------------------------------------------
C - Needed for bin boundaries, xlimit, scales
      INTEGER ISEPARATOR, NRAPMAX,NPTMAX,NBINTOTMAX,NXMAX,NSCALEMAX,
     >     NMAXSUBPROC,NSCALEBINMAX
C - Scenario-dependent limits
      PARAMETER(ISEPARATOR=1234567890)
      PARAMETER(NRAPMAX=50,NPTMAX=50,NBINTOTMAX=400,NXMAX=200) 
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
