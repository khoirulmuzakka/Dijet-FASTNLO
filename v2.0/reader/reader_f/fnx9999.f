***********************************************************************
*     
*     fastNLO reader code
*     
*     Initial version: T. Kluge, M. Wobisch, 2006
*     Updated for v2:  K. Rabbertz, M. Wobisch, 2011
*     
*     This code usually should not be edited.
*     
*     Contains:
*     ---------
*     FX9999CC  main routine
*     FX9999IN  initialize fastNLO table
*     FX9999PT  determine pointers to contributions and scales
*     FX9999PM  get PDFs, multiply with perturbative coefficients and
*     >         alpha_s
*     FX9999MT  multiply coefficients and PDFs
*     FX9999GP  read the PDFs in all bins
*     FX9999PL  compute PDF linear combinations
*     FX9999PR  print results
*     FX9999NF  print scenario information (physics & technical)
*     FX9999TB  returns information on variables from the table
*     FX9999NM  NOT IMPLEMENTED YET: norm. distr. by its own integral
*     
*     Requires:
*     ---------
*     fnx9999.inc common block definitions 
*     fx9999rw.f  read (or write) tables
*     
*     Scenario independent routines in :
*     ----------------------------------
*     fnio.f      contains routines for table I/O
*     fnset.f     set flags which define contributions to be considered
*     fn-interface.f:
*     FNALPHAS    alpha_s interface (double precision function)
*     FNPDF       PDF interface
*     
***********************************************************************

      SUBROUTINE FX9999CC(XMUR,XMUF,XSECT,XSUNCOR,XSCOR)
***********************************************************************
*     
*     fastNLO user code v2 - main routine 
*     
*     Input:
*     ------
*     XMUR  pre-factor for nominal renormalization scale     
*     >     any choice is possible, but please note 
*     >     that 2-loop threshold corrections work
*     >     only for xmur=xmuf
*     XMUF  pre-factor for nominal factorization scale
*     >     available choices depend on scenario
*     >     (see scenario information)
*     
*     Output:
*     -------
*     XSECT(NOBSBIN)     array of cross sections 
*     XSUNCOR(NOBSBIN,2) quadr. summed uncorr. uncertainty,
*     >                  lower (.,1), upper (.,2)
*     XSCOR(NOBSBIN,2)   quadr. summed corr. uncertainty,
*     >                  lower (.,1), upper (.,2)
*     
*     Current restrictions:
*     ---------------------
*     - can deal with single scale dim. only (no DIS with muf=Q, mur=pT)
*     - assumes same x section units (nb, pb, ...) for all contributions
*     - assumes that max. one contribution per allowed type is present
*     
***********************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION XMUR, XMUF
      INCLUDE 'fnx9999.inc'
      INTEGER IPOINT, I,J,K

*---  Reset output arrays
      DO I=1,MXOBSBIN
         XSECT(I)      = 0D0
         XSNORM(I)     = 0D0
         XSUNCOR(I,1)  = 0D0
         XSUNCOR(I,2)  = 0D0
         XSCOR(I,1)    = 0D0
         XSCOR(I,2)    = 0D0
*---  Also reset result array (in common block)
         DO J=0,MXSUBPROC
            DO K=0,MXCTRB
               RESULT(I,J,K) = 0D0
            ENDDO
         ENDDO
      ENDDO
      
*---  Determine pointers to access contributions and scales
      CALL FX9999PT(XMUR,XMUF,1)

*---  Differentiate between input requiring scales and alpha_s or not
*---  Start with perturbative contributions
      IF (ICONTRSELECTOR(ILO).EQ.1.OR.
     >     ICONTRSELECTOR(INLO).EQ.1.OR.
     >     ICONTRSELECTOR(ITHC1L).EQ.1.OR.
     >     ICONTRSELECTOR(ITHC2L).EQ.1) THEN

*---  Loop over logical contributions
         DO I=ILO,ITHC2L
            IF (ICONTRSELECTOR(I).EQ.1.AND.ICONTRPOINTER(I).NE.-1) THEN
               
*---  Get PDFs - multiply with perturbative coefficients and alpha_s
               CALL FX9999PM(I,XMUR,XMUF)
            ENDIF
         ENDDO

*---  Add results in output array
         DO I=ILO,ITHC2L
            IF (ICONTRSELECTOR(I).EQ.1.AND.ICONTRPOINTER(I).NE.-1) THEN
               IPOINT = ICONTRPOINTER(I)
               DO J=1,NOBSBIN
                  DO K=1,NSUBPROC(IPOINT)
C---  DO K=1,1 ! Test - only gg 
C---  DO K=2,2 ! Test - only g (DIS)
C---  DO K=2,5 ! Test - only qq
                     XSECT(J) = XSECT(J) + RESULT(J,K,I)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
*---  Either multiply x section by multiplicative correction or ...
         I = INPC1
         IF (ICONTRSELECTOR(I).EQ.1.AND.ICONTRPOINTER(I).NE.-1) THEN
            DO J=1,NOBSBIN
               XSECT(J) = XSECT(J) * MFACT(1,J)
            ENDDO
         ENDIF
*---  ... fill mult. correction and uncertainties into output arrays
*---  (Only works with one correction of this type for now)
      ELSEIF (ICONTRSELECTOR(INPC1).EQ.1.AND.
     >        ICONTRPOINTER(INPC1).NE.-1) THEN
         DO J=1,NOBSBIN
            XSECT(J) = MFACT(1,J)
*---  For now add sources of uncertainty quadratically
            DO K=1,NMUNCORREL(1)
               XSUNCOR(J,1) = XSUNCOR(J,1) +
     >              MUNCORLO(1,J,K)*MUNCORLO(1,J,K)
               XSUNCOR(J,2) = XSUNCOR(J,2) +
     >              MUNCORUP(1,J,K)*MUNCORUP(1,J,K)
            ENDDO
            XSUNCOR(J,1) = - SQRT(MIN(0D0,XSUNCOR(J,1)))
            XSUNCOR(J,2) = + SQRT(MIN(0D0,XSUNCOR(J,2)))
            DO K=1,NMCORREL(1)
               XSCOR(J,1) = XSCOR(J,1) +
     >              MCORLO(1,J,K)*MCORLO(1,J,K)
               XSCOR(J,2) = XSCOR(J,2) +
     >              MCORUP(1,J,K)*MCORUP(1,J,K)
            ENDDO
            XSCOR(J,1) = - SQRT(MIN(0D0,XSCOR(J,1)))
            XSCOR(J,2) = + SQRT(MIN(0D0,XSCOR(J,2)))
         ENDDO

*---  Fill data points and uncertainties into output arrays
      ELSEIF (ICONTRSELECTOR(IDATA).EQ.1.AND.
     >        ICONTRPOINTER(IDATA).NE.-1) THEN
         DO J=1,NOBSBIN
            XSECT(J) = DYVAL(J)
*---  For now add sources of uncertainty quadratically
            DO K=1,NDUNCORREL
               XSUNCOR(J,1) = XSUNCOR(J,1) +
     >              DUNCORLO(J,K)*DUNCORLO(J,K)
               XSUNCOR(J,2) = XSUNCOR(J,2) +
     >              DUNCORUP(J,K)*DUNCORUP(J,K)
            ENDDO
            XSUNCOR(J,1) = - SQRT(MAX(0D0,XSUNCOR(J,1)))
            XSUNCOR(J,2) = + SQRT(MAX(0D0,XSUNCOR(J,2)))
            DO K=1,NDCORREL
               XSCOR(J,1) = XSCOR(J,1) +
     >              DCORLO(J,K)*DCORLO(J,K)
               XSCOR(J,2) = XSCOR(J,2) +
     >              DCORUP(J,K)*DCORUP(J,K)
            ENDDO
            XSCOR(J,1) = - SQRT(MAX(0D0,XSCOR(J,1)))
            XSCOR(J,2) = + SQRT(MAX(0D0,XSCOR(J,2)))
         ENDDO
      ELSE
         WRITE(*,*)"FX9999CC: ERROR! Contribution selected "//
     >        "for output not defined, aborted."
         Do I=ILO,IDATA
            Write(*,*)"          IContr, NContr, Pointer, Selector",
     >           I,NCONTRCOUNTER(I),ICONTRPOINTER(I),ICONTRSELECTOR(I)
         ENDDO
         STOP
      ENDIF

*---  Normalization: Todo
      IF (INORMFLAG.EQ.0) THEN
         CONTINUE
      ELSEIF (INORMFLAG.EQ.1) THEN
*---  CALL FX9999NM          ! normalize by own integral
         WRITE(*,*)"FX9999CC: ERROR! Cross section normalization "//
     >        "not implemented yet, aborted. INormFlag = ",INORMFLAG
         STOP
      ELSEIF (INORMFLAG.EQ.2 .OR. INORMFLAG.EQ.3) THEN ! other normalization
         WRITE(*,*)"FX9999CC: ERROR! Cross section normalization "//
     >        "not implemented yet, aborted. INormFlag = ",INORMFLAG
      ELSE
         WRITE(*,*)"FX9999CC: ERROR! Illegal normalization "//
     >        "flag, aborted. INormFlag = ",INORMFLAG
      ENDIF

      RETURN
      END
***********************************************************************

      SUBROUTINE FX9999IN(FILENAME)
***********************************************************************
*     
*     Initialize fastNLO table
*     
*     Input: Filename of fastNLO table
*     
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INCLUDE 'strings.inc'
      INTEGER I,J,K, ICTRB(30), ISTAT,NSEP
      CHARACTER*(*) FILENAME
      CHARACTER*255 OLDFILENAME
      CHARACTER*72 CHEAD
      CHARACTER*80 CHTMP
      DATA OLDFILENAME/'xxxx'/

*---  Initialization
*---  In 1st scenario call: read fastNLO coefficient table
      IF (FILENAME.NE.OLDFILENAME) THEN
         CALL FX9999RW('read',FILENAME,ICTRB)
         OLDFILENAME = FILENAME
      ENDIF

      RETURN 
      END
***********************************************************************

      SUBROUTINE FX9999PT(XMUR,XMUF,IPRINT)
***********************************************************************
*     
*     Determine pointers to contributions and scales
*     
*     Input:
*     ------    
*     XMUR  pre-factor for nominal renormalization scale     
*     >     any choice is possible, but please note 
*     >     that 2-loop threshold corrections work
*     >     only for xmur=xmuf
*     XMUF  pre-factor for nominal factorization scale
*     >     available choices depend on scenario
*     >     (see scenario information)
*     IPRINT  verbosity
*     
*     Current restrictions:
*     ---------------------
*     - so far it is assumed that all scales are stored in first
*     > dimension, multiple scale dimensions are not yet implemented
*     > here
*     - so far it is assumed that (in one dimension) the scale factors
*     > for mur and muf are the same (i.e. can not have fixed muf, but
*     > different mur values
*     > (which, in principle would work with a posteriori variations)
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INTEGER IPRINT,I,J,K,I1
      DOUBLE PRECISION XMUR,XMUF

*---  Initialize counters and pointers for contribution types 
      DO I=1,MXCTRB
         NCONTRCOUNTER(I)  =  0
         ICONTRPOINTER(I)  = -1
         ICONTRSELECTOR(I) =  0
         ISCALEPOINTER(I)  = -1
      ENDDO

*---  Check input
      IF (XMUR.LT.1D-3.OR.XMUF.LT.1D-3) THEN
         WRITE(*,*)"FX9999PT: ERROR! Scale factors smaller than "//
     >        "0.001 are not allowed, stopped! xmur, xmuf = ",xmur,xmuf
         STOP
      ELSEIF (XMUR.GT.1D3.OR.XMUF.GT.1D3) THEN
         WRITE(*,*)"FX9999PT: ERROR! Scale factors largeer than "//
     >        "1000. are not allowed, stopped! xmur, xmuf = ",xmur,xmuf
         STOP
      ENDIF

*---  Loop once over all contributions and register types
      DO I=1,NCONTRIB
         IF (ICONTRFLAG1(I).EQ.1.AND.ICONTRFLAG2(I).EQ.1.AND.
     >        IREF(I).EQ.PREFTAB) THEN
            NCONTRCOUNTER(ILO) = NCONTRCOUNTER(ILO) + 1
            ICONTRPOINTER(ILO) = I
         ELSEIF (ICONTRFLAG1(I).EQ.1.AND.ICONTRFLAG2(I).EQ.2.AND.
     >           IREF(I).EQ.PREFTAB) THEN
            NCONTRCOUNTER(INLO) = NCONTRCOUNTER(INLO) + 1
            ICONTRPOINTER(INLO) = I
         ELSEIF (ICONTRFLAG1(I).EQ.2.AND.ICONTRFLAG2(I).EQ.1.AND.
     >           IREF(I).EQ.PREFTAB) THEN
            NCONTRCOUNTER(ITHC1L) = NCONTRCOUNTER(ITHC1L) + 1
            ICONTRPOINTER(ITHC1L) = I
         ELSEIF (ICONTRFLAG1(I).EQ.2.AND.ICONTRFLAG2(I).EQ.2.AND.
     >           IREF(I).EQ.PREFTAB) THEN
            NCONTRCOUNTER(ITHC2L) = NCONTRCOUNTER(ITHC2L) + 1
            ICONTRPOINTER(ITHC2L) = I
         ELSEIF (ICONTRFLAG1(I).EQ.4.AND.ICONTRFLAG2(I).EQ.1.AND.
     >           IADDMULTFLAG(I).EQ.1.AND.
     >           IREF(I).EQ.PREFTAB) THEN
            NCONTRCOUNTER(INPC1) = NCONTRCOUNTER(INPC1) + 1
            ICONTRPOINTER(INPC1) = I
         ELSEIF (ICONTRFLAG1(I).EQ.0.AND.ICONTRFLAG2(I).EQ.0.AND.
     >           IDATAFLAG(I).EQ.1.AND.
     >           IREF(I).EQ.PREFTAB) THEN
            NCONTRCOUNTER(IDATA) = NCONTRCOUNTER(IDATA) + 1
            ICONTRPOINTER(IDATA) = I
         ELSE
            WRITE(*,*)"FX9999PT: ERROR! Unknown contribution type "//
     >           "encountered, stopped!"
            WRITE(*,*)"          IContrFlag1, IContrFlag2 = ",
     >           ICONTRFLAG1(I),ICONTRFLAG2(I)
            WRITE(*,*)"          IAddMultFlag,IDataFlag = ",
     >           IADDMULTFLAG(I),IDATAFLAG(I)
            WRITE(*,*)"          Iref = ",
     >           IREF(I)
            STOP
         ENDIF
      ENDDO
      
*---  Find particular contributions
*---  Find LO contribution
      IF (PORDPTHY.GE.1) THEN
         IF (ICONTRPOINTER(ILO).EQ.-1) THEN
            WRITE(*,*)"FX9999PT: ERROR! Requested contribution "//
     >           "not available, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY
            STOP
         ELSEIF (NCONTRCOUNTER(ILO).NE.1) THEN
            WRITE(*,*)"FX9999PT: ERROR! More than one contribution "//
     >           "available, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY
            WRITE(*,*)"          icontr = ",ICONTRPOINTER(ILO),
     >           "NContrCounter = ",NCONTRCOUNTER(ILO)
            STOP
         ELSE
            ICONTRSELECTOR(ILO) = 1
         ENDIF
      ENDIF
      
*---  Find NLO contribution
      IF (PORDPTHY.GE.2) THEN
         IF (ICONTRPOINTER(INLO).EQ.-1) THEN
            WRITE(*,*)"FX9999PT: ERROR! Requested contribution "//
     >           "not available, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY
            STOP
         ELSEIF (NCONTRCOUNTER(INLO).NE.1) THEN
            WRITE(*,*)"FX9999PT: ERROR! More than one contribution "//
     >           "available, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY
            WRITE(*,*)"          icontr = ",ICONTRPOINTER(INLO),
     >           "NContrCounter = ",NCONTRCOUNTER(INLO)
            STOP
         ELSE
            ICONTRSELECTOR(INLO) = 1
         ENDIF
      ENDIF
      
*---  Find 1-loop TC
      If (PTHRESHCOR.EQ.1) THEN
         If (PORDPTHY.LT.0) THEN
            WRITE(*,*)"FX9999PT: ERROR! Illegal choice of "//
     >           "perturbative order, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY,
     >           ", PTHRESHCOR = ",PTHRESHCOR
            STOP
         ELSEIf (PORDPTHY.EQ.0) THEN
            WRITE(*,*)"FX9999PT: Warning! 1-loop TC returned "//
     >           "separately, but need to be added to LO!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY,
     >           ", PTHRESHCOR = ",PTHRESHCOR
         ELSEIf (PORDPTHY.GE.2) THEN
            WRITE(*,*)"FX9999PT: ERROR! Inconsistent choice of "//
     >           "1-loop TC with NLO, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY,
     >           ", PTHRESHCOR = ",PTHRESHCOR
            STOP
         ENDIF
         IF (ICONTRPOINTER(ITHC1L).EQ.-1) THEN
            WRITE(*,*)"FX9999PT: ERROR! Requested contribution "//
     >           "not available, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY,
     >           ", PTHRESHCOR = ",PTHRESHCOR
            STOP
         ELSEIF (NCONTRCOUNTER(ITHC1L).NE.1) THEN
            WRITE(*,*)"FX9999PT: ERROR! More than one contribution "//
     >           "available, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY,
     >           ", PTHRESHCOR = ",PTHRESHCOR
            WRITE(*,*)"          icontr = ",ICONTRPOINTER(ITHC1L),
     >           "NContrCounter = ",NCONTRCOUNTER(ITHC1L)
            STOP
         ELSE
            ICONTRSELECTOR(ITHC1L) = 1
         ENDIF
      ENDIF

*---  Find 2-loop TC
      IF (PTHRESHCOR.EQ.2) THEN
         IF (PORDPTHY.LT.0) THEN
            WRITE(*,*)"FX9999PT: ERROR! Illegal choice of "//
     >           "perturbative order, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY,
     >           ", PTHRESHCOR = ",PTHRESHCOR
            STOP
         ELSEIf (PORDPTHY.EQ.0) THEN
            WRITE(*,*)"FX9999PT: Warning! 2-loop TC returned "//
     >           "separately, but need to be added to NLO!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY,
     >           ", PTHRESHCOR = ",PTHRESHCOR
         ELSEIf (PORDPTHY.EQ.1) THEN
            WRITE(*,*)"FX9999PT: ERROR! Inconsistent choice of "//
     >           "2-loop TC with LO, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY,
     >           ", PTHRESHCOR = ",PTHRESHCOR
            STOP
         ELSEIf (PORDPTHY.GE.3) THEN
            WRITE(*,*)"FX9999PT: ERROR! Inconsistent choice of "//
     >           "2-loop TC with N?LO, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY,
     >           ", PTHRESHCOR = ",PTHRESHCOR
            STOP
         ENDIF
         IF (ICONTRPOINTER(ITHC2L).EQ.-1) THEN
            WRITE(*,*)"FX9999PT: ERROR! Requested contribution "//
     >           "not available, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY,", PTHRESHCOR = "
     >           ,PTHRESHCOR
            STOP
         ELSEIF (NCONTRCOUNTER(ITHC2L).NE.1) THEN
            WRITE(*,*)"FX9999PT: ERROR! More than one contribution "//
     >           "available, stopped!"
            WRITE(*,*)"          PORDPTHY = ",PORDPTHY,", PTHRESHCOR = "
     >           ,PTHRESHCOR
            WRITE(*,*)"          icontr = ",ICONTRPOINTER(ITHC2L),
     >           "NContrCounter = ",NCONTRCOUNTER(ITHC2L)
            STOP
         ELSE
            ICONTRSELECTOR(ITHC2L) = 1
         ENDIF
      ENDIF
      
*---  Find NP corrections
      IF (PNPCOR.EQ.1) THEN
         IF (ICONTRPOINTER(INPC1).EQ.-1) THEN
            WRITE(*,*)"FX9999PT: ERROR! Requested contribution "//
     >           "not available, stopped!"
            WRITE(*,*)"          PNPCOR = ",PNPCOR
            STOP
         ELSEIF (NContrCounter(INPC1).ne.1) THEN
            WRITE(*,*)"FX9999PT: ERROR! More than one contribution "//
     >           "available, stopped!"
            WRITE(*,*)"          PNPCOR = ",PNPCOR
            WRITE(*,*)"          icontr = ",ICONTRPOINTER(INPC1),
     >           "NContrCounter = ",NCONTRCOUNTER(INPC1)
            STOP
         ELSE
            ICONTRSELECTOR(INPC1) = 1
         ENDIF
      ENDIF

*---  Find Data
      IF (PDATA.EQ.1) THEN
         IF (ICONTRPOINTER(IDATA).EQ.-1) THEN
            WRITE(*,*)"FX9999PT: ERROR! Requested contribution "//
     >           "not available, stopped!"
            WRITE(*,*)"          PDATA = ",PDATA
            STOP
         ELSEIF (NContrCounter(idata).ne.1) THEN
            WRITE(*,*)"FX9999PT: ERROR! More than one contribution "//
     >           "available, stopped!"
            WRITE(*,*)"          PDATA = ",PDATA
            WRITE(*,*)"          icontr = ",ICONTRPOINTER(IDATA),
     >           "NContrCounter = ",NCONTRCOUNTER(IDATA)
            STOP
         ELSE
            ICONTRSELECTOR(IDATA) = 1
         ENDIF
      ENDIF

*---  Print contribution summary (if requested)
      IF (IPRINT.GT.1) THEN
         WRITE(*,'(A)')
     >        " FX9999PT: Contribution summary:"
         WRITE(*,'(A)')
     >        " FX9999PT: Name  ,  Enum, IPoint, NCount, ISelct"
         DO I=ILO,IDATA
            WRITE(*,'(A,I6,I8,I8,I8)')" FX9999PT: "//CNAME(I)//":",
     >           I,ICONTRPOINTER(I),NCONTRCOUNTER(I),
     >           ICONTRSELECTOR(I)
         ENDDO
      ENDIF

*---  Check availability of scale choices and assign pointers
*---  The current treatment works only for a single scale dimension!
*---  For 2nd scale dimension need 2nd scale pointer!
*     (Selection is based on factorization scale since renormalization
*     scale is flexible - except for threshold corrections with
*     IScaleDep = 2)
      DO I=ILO,ITHC2L
         IF (ICONTRSELECTOR(I).EQ.1.AND.ICONTRPOINTER(I).NE.-1) THEN
            I1 = ICONTRPOINTER(I)
            ISCALEPOINTER(I) = 0
            IF (ISCALEDEP(I1).EQ.0) THEN ! Born-type w/o scale dep - use any scale
               IF (NSCALEVAR(I1,1).GE.1) THEN
                  ISCALEPOINTER(I) = 1 ! Use 1st scale variation
               ELSE
                  WRITE(*,*)"FX9999PT: ERROR! Not a single scale "//
     >                 "available in contribution, stopped!"
                  WRITE(*,*)"          IContr = ",ICONTRPOINTER(I),
     >                 ", NScaleVar(.,1) = ",NSCALEVAR(I1,1)
                  STOP
               ENDIF
               IF (NSCALEVAR(I1,1).GT.1) THEN
                  WRITE(*,*)"FX9999PT: WARNING! Why more than one "//
     >                 " scale variation for contribution"
                  WRITE(*,*)
     >                 "        with scale-independent coefficients?"
               ENDIF
*---  Check for presence of more restrictive factorization scale
            ELSE
               DO J=1,NSCALEVAR(I1,1)
                  IF (ABS(SCALEFAC(I1,1,J)-XMUF).LT.TINY(1.D0)) THEN
                     ISCALEPOINTER(I)=J
                     EXIT       ! Quit DO loop at first match
                  ENDIF
               ENDDO
            ENDIF
            IF (ISCALEPOINTER(I).EQ.0) THEN
               IF (IPRINT.GT.0) THEN
                  WRITE(*,*)"FX9999PT: WARNING! The requested "//
     >                 " factorization scale xmuf = ",xmuf
                  WRITE(*,*)"          is not available!"
                  DO J=1,NCONTRDESCR(I1)
                     WRITE(*,*)'  ',CTRBDESCRIPT(I1,J)
                  ENDDO
                  STOP
               ENDIF
            ENDIF 
         ENDIF
      ENDDO
      
*---  For IScaleDep = 2 check if renormalization scale is directly
*---  available, i.e. identical to factorization scale
*---  Can be provided a posteriori otherwise
      IF (ABS(XMUR-XMUF).GT.TINY(1D0)) THEN
         DO I=ILO,ITHC2L
            IF (ICONTRSELECTOR(I).EQ.1.AND.ICONTRPOINTER(I).NE.-1) THEN
               I1 = ICONTRPOINTER(I)
               IF (ISCALEDEP(I1).EQ.2) THEN
                  ISCALEPOINTER(I) = -1
                  IF (IPRINT.GT.0) THEN
                     WRITE(*,*)"FX9999PT: WARNING! The requested "//
     >                    " renormalization scale xmur = ",XMUR
                     WRITE(*,*)"          is not available, stopped!"
                     WRITE(*,*)"          Only xmur=xmuf is possible."
                     DO J=1,NCONTRDESCR(I1)
                        WRITE(*,*)'  ',CTRBDESCRIPT(I1,J)
                     ENDDO
                     STOP
                  ENDIF
               ENDIF
            ENDIF
            IF (IPRINT.GT.1) THEN
               WRITE(*,*) "FX9999PT: Pointer number ",I," to "//
     >              "IContr, IScale:",ICONTRPOINTER(I),ISCALEPOINTER(I)
            ENDIF
         ENDDO
      ENDIF
      
      RETURN
      END
***********************************************************************

      SUBROUTINE FX9999PM(ICTRB,XMUR,XMUF)
***********************************************************************
*
*     Get PDFs and multiply with alpha_s and coefficients
*     
*     Input:
*     ------    
*     ICTRB  number of contribution (logical - not in table!)
*     XMUR  pre-factor for nominal renormalization scale     
*     >     any choice is possible, but please note 
*     >     that 2-loop threshold corrections work
*     >     only for xmur=xmuf
*     XMUF  pre-factor for nominal factorization scale
*     >     available choices depend on scenario
*     >     (see scenario information)
*     
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INTEGER ICTRB, IADDPOW, IC, IS
      DOUBLE PRECISION XMUR, XMUF
      DOUBLE PRECISION FACTOR, LOGMUR
      DOUBLE PRECISION NF, CA, CF, BETA0, BETA1
      PARAMETER (NF=5D0, CA=3D0, CF=4D0/3D0)
      PARAMETER (BETA0=(11D0*CA-2D0*NF)/3D0) 
      PARAMETER (BETA1=34D0*CA*CA/3D0-2D0*NF*(CF+5D0*CA/3D0))

*---  Compute 'standard' contribution at chosen scales
      IADDPOW = 0
      FACTOR  = 1D0
      CALL FX9999GP(ICTRB,XMUF)
      CALL FX9999MT(ICTRB,XMUR,XMUF,IADDPOW,FACTOR)
      
*---  A posteriori MuR variation for scale dependent parts
      IC = ICONTRPOINTER(ICTRB)
      IS = ISCALEPOINTER(ICTRB)
      IF (ISCALEDEP(IC).NE.0.AND.ABS(XMUR-XMUF).GT.TINY(1D0)) THEN
         IF (ISCALEDEP(IC).GT.1) THEN
            WRITE(*,*)'FX9999PM: ERROR! A posteriori variation '//
     >           'of MuR for contribution type',ICTRB,
     >           'with scale dependence',ISCALEDEP(IC),
     >           'not implemented! Stopped.'
            STOP
         ENDIF
   
*---  Abs. order of LO:     ILOORD
*---  Abs. order in alphas: NPOW(IC)
*---  >  NLO if (NPOW-ILOORD) = 1, NNLO if (NPOW(IC)-ILOORD) = 2

*---  If LO
         IF ( (NPOW(IC)-ILOORD).EQ.0 ) THEN
*---  Nothing to do, direct MuR dependence via alpha_s already taken
*---  into account by first call to FX9999MT. Also, MuRcache
*---  was filled in FX9999MT for a posteriori use.
*---  If NLO 
         ELSEIF ( (NPOW(IC)-ILOORD).EQ.1 ) THEN
            IADDPOW = 1
            LOGMUR  = LOG(XMUR/SCALEFAC(IC,1,IS))
            FACTOR  = DBLE(ILOORD)*BETA0*LOGMUR
            CALL FX9999GP(ICTRB-1,XMUF)
*---  For NLO, ICTRB-1 is LO as required to store c_1 dependent
*---  NLO scale modification. However, this has to go one order higher,
*---  that is IADDPOW = 1, and also use the cached MuR values of the NLO
*---  scale nodes.
            CALL FX9999MT(ICTRB-1,XMUR,XMUF,IADDPOW,FACTOR) ! 1: mod NLO
*---  If NNLO
         ELSEIF ( (NPOW(IC)-ILOORD).EQ.2 ) THEN
            LOGMUR = LOG(XMUR/SCALEFAC(IC,1,IS))
            WRITE(*,*)'FX9999PM: ERROR! A posteriori scale variation'
     >           //' beyond NLO not yet implemented! Stopped.'
            STOP
C---  Factor = ...dble(ILOord)*beta0*logmur ! n beta0 logmu
            CALL FX9999GP(ICTRB-2,XMUF)
            CALL FX9999MT(ICTRB-2,XMUR,XMUF,2,FACTOR) ! 2: mod LO
            
C---  Factor = ...dble(ILOord)*beta0*logmur ! n beta0 logmu
            CALL FX9999GP(ICTRB-1,XMUF)
            CALL FX9999MT(ICTRB-1,XMUR,XMUF,1,FACTOR) ! 1: mod NLO
         ELSE
            WRITE(*,*)'FX9999PM: ERROR! A posteriori scale variation'
     >           //' beyond NNLO not implemented! Stopped.'
            STOP
         ENDIF
      ENDIF

      RETURN
      END
***********************************************************************

      SUBROUTINE FX9999MT(ICTRB,XMUR,XMUF,IADDPOW,FACTOR)
***********************************************************************
*     
*     Multiply the PDFs and the perturbative coefficients 
*     
*     Input:
*     ------    
*     ICTRB number of logical contribution 
*     XMUR  pre-factor for nominal renormalization scale     
*     >     any choice is possible, but please note 
*     >     that 2-loop threshold corrections work
*     >     only for xmur=xmuf
*     XMUF  pre-factor for nominal factorization scale
*     >     available choices depend on scenario
*     >     (see scenario information)
*     IADDPOW add. power in alphas - for a posteriori mur variation
*     >       also: increment in result-contribution index where 
*     >       result is stored
*     FACTOR  =! 1 for a posteriori mur variations
*     
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INTEGER IADDPOW
      INTEGER IXMUR,IXMUF, ICTRB,IC,IS, I,J,K,L,M,N, NXMAX,
     >     ISTORE
      DOUBLE PRECISION XMUR, XMUF, FACTOR, MUR, FNALPHAS, SCF
      DOUBLE PRECISION AS, ASPOW
      DOUBLE PRECISION MURCACHE(MXOBSBIN,MXSCALENODE)

*---  Set pointers to contribution in table and to scale variation
      IC = ICONTRPOINTER(ICTRB)
      IS = ISCALEPOINTER(ICTRB)

*---  If regular contribution (i.e. no a posteriori scale variation)
*---  then store in slot number of logical contribution,
*---  otherwise increment by no. of add powers in alphas
      ISTORE = ICTRB + IADDPOW

*---  Loop over coefficient array - compare with order in table !
*---  Storage loop: observable, scalebins, (get alphas), xbins, subproc
      DO J=1,NOBSBIN
         IF (NPDFDIM(IC).EQ.0) THEN ! 1-d case (DIS)
            NXMAX =  NXTOT(IC,1,J) 
         ELSEIF (NPDFDIM(IC).EQ.1) THEN ! 2-d half-matrix
            NXMAX = (NXTOT(IC,1,J)*NXTOT(IC,1,J)+NXTOT(IC,1,J))/2 
         ELSEIF (NPDFDIM(IC).EQ.2) THEN ! 2-d full matrix
            WRITE(*,*)'FX9999MT: ERROR! 2-d case not yet implemented!'
     >           //' Stopped.'
            STOP
         ENDIF
         DO K=1,NSCALENODE(IC,1)
*---  If regular contribution (i.e. no a posteriori scale variation)
*---  then store MuR in cache for use with a posteriori scale variation
            IF (IADDPOW.EQ.0) THEN
               MUR = XMUR / SCALEFAC(IC,1,IS) * SCALENODE(IC,J,1,IS,K) 
               MURCACHE(J,K) = MUR
*---  Otherwise retrieve previously (Check! TBD) cached MuR values
            ELSE
               MUR = MURCACHE(J,K)
            ENDIF
*---  Get alpha_s
            AS =  FNALPHAS(MUR)
            ASPOW = AS**DBLE(NPOW(IC)+IADDPOW)
            DO L=1,NXMAX
               DO M=1,NSUBPROC(IC)
                  IF (PREFTAB.EQ.0) THEN
                     RESULT(J,M,ISTORE) = RESULT(J,M,ISTORE) + 
     >                    SIGMATILDE(IC,J,1,IS,K,L,M)
     >                    * ASPOW
     >                    * PDF(J,K,L,M)
     >                    * FACTOR
                  ELSE          ! Reference Table
                     RESULT(J,M,ISTORE) = RESULT(J,M,ISTORE) + 
     >                    SIGMATILDE(IC,J,1,IS,K,L,M)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
      RETURN 
      END
***********************************************************************

      SUBROUTINE FX9999GP(ICTRB,XMUF)
***********************************************************************
*     
*     Read the PDFs in all bins - at all xmax, xmin bins
*     the default factorization scale (in GeV) is multiplied by
*     muffactor
*     
*     Input:
*     ------    
*     ICTRB number of logical contribution 
*     XMUF  pre-factor for nominal factorization scale
*     >     available choices depend on scenario
*     >     (see scenario information)
*
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      DOUBLE PRECISION XMUF, X, MUF, 
     >     TMPPDF(-6:6), XPDF1(MXNXTOT,-6:6),XPDF2(MXNXTOT,-6:6), H(10)
      INTEGER ICTRB, IC,IS, I,J,K,L,M, NX, NX2LIMIT

*---  Set pointers to contribution in table and to scale variation
      IC = ICONTRPOINTER(ICTRB)
      IS = ISCALEPOINTER(ICTRB)
      
      DO I=1,NOBSBIN
         DO J=1,NSCALENODE(IC,1)
*---  For MuF independent contributions (e.g. LO) the scale factor XMUF
*---  has to be applied to the scale IS stored within the table
*---  contribution IC.
*---  For MuF dependent contributions the scale factor XMUF is
*---  integrated already in the properly chosen scale IS stored within
*---  the table contribution IC.
*---  Note:
*---  For threshold corrections MuR and MuF must be identical!
*---  The LO jet cross sections DO change with MuF for identical
*---  MuR because of the PDF's. Only the first perturbative coefficient
*---  is independent of MuF.
            IF (ISCALEDEP(IC).EQ.0) THEN
               MUF = XMUF*SCALENODE(IC,I,1,IS,J) ! 1: ScaleDim  IS: ScaleVar
            ELSE
               MUF = SCALENODE(IC,I,1,IS,J)
            ENDIF
            NX =0
            DO K=1,NXTOT(IC,1,I) ! Fill first PDF
               X = XNODE1(IC,I,K) 
               CALL FNPDF(X, MUF, TMPPDF)
               DO M=-6,6
                  XPDF1(K,M) = TMPPDF(M)
                  IF (NPDF(IC).EQ.1) THEN ! DIS
                     CONTINUE   
                  ELSEIF (NPDF(IC).EQ.2) THEN ! Two hadrons
                     IF (NPDFPDG(IC,1).EQ.NPDFPDG(IC,2)) THEN ! Identical hh
                        XPDF2(K,M) = TMPPDF(M)
                     ELSEIF (NPDFPDG(IC,1).EQ.-NPDFPDG(IC,2)) THEN ! hh-bar
                        XPDF2(K,-M) = TMPPDF(M)
                     ELSE
                        WRITE(*,*)
     >                       'FX9999GP: ERROR! So far only the '//
     >                       'scattering of identical hadrons or '
                        WRITE(*,*)'          hadron anti-hadron is'//
     >                       ' implemented. gamma-p to be done...'
                        WRITE(*,*)'          Stopped'
                        STOP
                     ENDIF
                  ELSE
                     WRITE(*,*)'FX9999GP: ERROR! Neither one nor '//
     >                    'two hadrons ... ? NPDF = ',NPDF(ic)
                     STOP
                  ENDIF
               ENDDO            ! M
            ENDDO               ! K 
            IF (NPDFDIM(IC).EQ.2) THEN ! Fill second PDF
               DO K=1,NXTOT(IC,2,I)
                  X = XNODE2(IC,I,L)
                  CALL FNPDF(X, MUF, TMPPDF)
                  DO M=-6,6
                     XPDF2(J,M) = TMPPDF(M)
                  ENDDO
               ENDDO
            ENDIF
            DO K=1,NXTOT(IC,1,I) ! Build PDF 1,2 linear combinations
               IF (NPDF(IC).EQ.1) NX2LIMIT = 1 ! 1-d case (DIS)
               IF (NPDFDIM(IC).EQ.1) NX2LIMIT = K ! 2-d half matrix
               IF (NPDFDIM(IC).EQ.2) NX2LIMIT = NXTOT(IC,2,I) ! 2-d full matrix
               DO L=1,NX2LIMIT
                  NX = NX+1
                  CALL FX9999PL(IPDFDEF(IC,1),IPDFDEF(IC,2),
     >                 IPDFDEF(IC,3),K,L,XPDF1,XPDF2,H)
                  DO M=1,NSUBPROC(IC)
                     PDF(I,J,NX,M) = H(M)
                  ENDDO         ! m subproc
               ENDDO            ! x2-loop
            ENDDO               ! x1-loop
         ENDDO                  ! ScaleNode-loop
      ENDDO                     ! ObsBin-loop

      RETURN
      END
***********************************************************************

      SUBROUTINE FX9999PL(ICF1,ICF2,ICF3,I,J,XPDF1,XPDF2,H)
***********************************************************************
*     
*     Compute PDF linear combinations - for different sets of
*     subprocesses
*     
*     Depending on icf1,2,3, the product of the i-th and j-th entries
*     of the PDF array XPDF are multiplied into the relevant linear 
*     combinations in the array H
*     
*     Input:
*     ------    
*     IREACT             Flag for reaction (1:DIS, 2:pp, 3:ppbar)
*     I                  x-index of first hadron        
*     J                  x-index of second hadron (if pp, ppbar)
*     XPDF1(NXMAX,-6:6)  PDF array for all x-bins
*     XPDF2(NXMAX,-6:6)  PDF array for all x-bins
*     
*     Output:
*     -------    
*     H(10)              PDF linear combinations
*     
*     Current restrictions:
*     - Check order of DIS subprocesses ...!
*
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INTEGER ICF1,ICF2,ICF3, I,J,K
      DOUBLE PRECISION XPDF1(MXNXTOT,-6:6),XPDF2(MXNXTOT,-6:6), H(10),
     >     G1, G2,              ! gluon densities from both hadrons
     >     SUMQ1, SUMQ2,        ! sum of quark densities
     >     SUMQB1, SUMQB2,      ! sum of anti-quark densities
     >     Q1(6), Q2(6), QB1(6), QB2(6), ! arrays of 6 (anti-)quark densities
     >     S,A                  ! products S,A

*---  DIS: Inclusive and jets, gammaP direct. TBD
      IF (ICF1.EQ.2 .AND. (ICF2.EQ.0 .OR. ICF2.EQ.1)) THEN 
C---  H(3) = 0d0             ! Delta  at O(as^0)
C---  Do k=1,5,2
C---  H(3) = H(3) + (XPDF1(i,k)+XPDF1(i,-k)+
C---  >           4d0*(XPDF1(i,k+1)+XPDF1(i,-k-1)))/9d0
C---  ENDDO
C---  H(1) = XPDF1(i,0)      ! Gluon  at O(as^1)
C---  H(2) = 0d0             ! Sigma  at O(as^2)
C---  Do k=1,6
C---  H(2) = H(2)+XPDF1(i,k)+XPDF1(i,-k)
C---  ENDDO

*---  This is not final - just for an early DIS table. TBD
C---  H(2) = 0d0             ! Delta  at O(as^0)
C---  Do k=1,5,2
C---  H(2) = H(2) + (XPDF1(i,k)+XPDF1(i,-k)+
C---  >           4d0*(XPDF1(i,k+1)+XPDF1(i,-k-1)))/9d0
C---  ENDDO
C---  H(1) = XPDF1(i,0)      ! Gluon  at O(as^1)
C---  H(3) = 0d0             ! Sigma  at O(as^2)
C---  Do k=1,6
C---  H(3) = H(3)+XPDF1(i,k)+XPDF1(i,-k)
C---  ENDDO

*---  Final: 1 Delta, 2 Gluon, 3 Sigma.
         H(1) = 0D0             ! Delta  at O(as^0)
         DO K=1,5,2
            H(1) = H(1) + (XPDF1(I,K)+XPDF1(I,-K)+
     >           4D0*(XPDF1(I,K+1)+XPDF1(I,-K-1)))/9D0
         ENDDO
         H(2) = XPDF1(I,0)      ! Gluon  at O(as^1)
         H(3) = 0D0             ! Sigma  at O(as^2)
         DO K=1,6
            H(3) = H(3)+XPDF1(I,K)+XPDF1(I,-K)
         ENDDO

*---  hadron-hadron: jets
      ELSEIF (ICF1.EQ.3.AND.ICF2.EQ.1.AND.(ICF3.GE.1.AND.ICF3.LE.2))THEN
         SUMQ1  = 0D0
         SUMQB1 = 0D0
         SUMQ2  = 0D0
         SUMQB2 = 0D0
         DO K=1,6
            Q1(K)  = XPDF1(I,K) ! Read 1st PDF at x1
            QB1(K) = XPDF1(I,-K)
            SUMQ1  = SUMQ1  + Q1(K)
            SUMQB1 = SUMQB1 + QB1(K)
            Q2(K)  = XPDF2(J,K) ! Read 2nd PDF at x2
            QB2(K) = XPDF2(J,-K)
            SUMQ2  = SUMQ2  + Q2(K)
            SUMQB2 = SUMQB2 + QB2(K)
         ENDDO
         G1     = XPDF1(I,0)
         G2     = XPDF1(J,0)
*---  Compute S, A
         S = 0D0
         A = 0D0
         DO K=1,6
            S = S + (Q1(K)*Q2(K)) + (QB1(K)*QB2(K)) 
            A = A + (Q1(K)*QB2(K)) + (QB1(K)*Q2(K)) 
         ENDDO
*---  Compute seven combinations
         H(1) = G1*G2
         H(2) = SUMQ1*SUMQ2 + SUMQB1*SUMQB2 - S
         H(3) = S
         H(4) = A
         H(5) = SUMQ1*SUMQB2 + SUMQB1*SUMQ2 - A
         H(6) = (SUMQ1+SUMQB1)*G2
         H(7) = G1*(SUMQ2+SUMQB2)
         IF (ICF3.EQ.1) H(6) = H(6)+H(7) ! Case: 6 subprocesses

*---  gammaP: direct, jets
      ELSEIF (ICF1.EQ.2 .AND. ICF2.EQ.2) THEN 
         WRITE(*,*)'FX9999PL: gammaP to be implemented, stopped!'
         STOP
      ELSE
         WRITE(*,*)'FX9999PL: Combination not defined:'
         WRITE(*,*)'          icf1,2,3 = ',ICF1,ICF2,ICF3
         WRITE(*,*)'          Stopped.'
         STOP
      ENDIF

      RETURN
      END
***********************************************************************

      SUBROUTINE FX9999PR(XSECT)
***********************************************************************
*
*     Print fastNLO cross section array
*     
*     Input:
*     ------
*     XSECT  Result array
*     
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INCLUDE 'strings.inc'
      INTEGER I,J
      DOUBLE PRECISION MIN2OLD,MAX2OLD

      WRITE(*,*)"fastNLO results:     ndim=",NDIM
      WRITE(*,9003)IXSECTUNITS(1),CUNITS(IXSECTUNITS(1))

      MIN2OLD=-99999D0
      MAX2OLD=-99999D0

      IF (NDIM.EQ.1) 
     >     WRITE(*,*)' BINS IN ',DIMLABEL(NDIM)

      DO I=1,NOBSBIN
         IF (NDIM.GT.1) THEN
            IF ((MIN2OLD.NE.LOBIN(I,2)).OR.
     >           (MAX2OLD.NE.UPBIN(I,2))) THEN
               MIN2OLD = LOBIN(I,2)
               MAX2OLD = UPBIN(I,2)
               WRITE(*,*)'------------------------------------'
               WRITE(*,9005) MIN2OLD,MAX2OLD,DIMLABEL(NDIM)
               WRITE(*,*)' bins in ',DIMLABEL(1)
            ENDIF
         ENDIF
         WRITE(*,9010)I,LOBIN(I,1),UPBIN(I,1),XSECT(I)
      ENDDO
 9003 FORMAT(' results in units of 10^-',I2,' barn (',A2,')')
 9005 FORMAT('  range ',F10.3,' -',F10.3,' in ',A64)
 9010 FORMAT(I4,F9.3,' -',F9.3,' :',E12.4)

      RETURN
      END
***********************************************************************

      SUBROUTINE FX9999NF
***********************************************************************
*
*     fastNLO user code v2.0 - print scenario information
*
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INCLUDE 'strings.inc'
      INTEGER I,J,K
      CHARACTER*80 CHTMP

      WRITE(*,'(A)')
      WRITE(*,*)CSEPS
      WRITE(*,*)"# Information on fastNLO scenario: ",
     >     SCENNAME(1:LEN_TRIM(SCENNAME))
      WRITE(*,*)LSEPS
      WRITE(*,*)"# Description:"
      DO I=1,NSCDESCRIPT
         CHTMP = SCDESCRIPT(I)
         WRITE(*,*)"#   ",CHTMP(1:LEN_TRIM(CHTMP))  
      ENDDO
      WRITE(*,*)"#"
      WRITE(*,'(A,G10.4,A)')
     >     " # Centre-of-mass energy Ecms: ",
     >     ECMS," GeV"
      WRITE(*,*)"#"
      WRITE(*,'(A,I3,A,I1,A)')
     >     " # Tot. no. of observable bins: ",NOBSBIN," in ",NDIM,
     >     " dimensions:"
      WRITE(*,*)"#"
      WRITE(*,'(A,I1)')" # No. of contributions: ",NCONTRIB
      DO I=1,NCONTRIB         
         WRITE(*,'(A,I1,A)')" # Contribution ",I,":"
         DO J=1,NCONTRDESCR(I)
            CHTMP = CTRBDESCRIPT(I,J)
            WRITE(*,*)"#   ",CHTMP(1:LEN_TRIM(CHTMP))  
         ENDDO
         WRITE(*,'(A,I16)')" #   No. of events: ",NEVT(I)
         WRITE(*,*)"#   provided by:"
         DO J=1,NCODEDESCR(I)
            CHTMP = CODEDESCRIPT(I,J)
            WRITE(*,*)"#   ",CHTMP(1:LEN_TRIM(CHTMP))  
         ENDDO
         WRITE(*,'(A,I1)')" #   Scale dimensions: ",
     >        NSCALEDIM(I)
         DO J=1,NSCALEDIM(I)
            DO K=1,NSCALEDESCRIPT(I,J)
               CHTMP = SCALEDESCRIPT(I,J,K)
               WRITE(*,'(A,I1,A,A)')
     >              " #     Scale description for dimension ",
     >              J,":          ",CHTMP(1:LEN_TRIM(CHTMP))  
            ENDDO
            WRITE(*,'(A,I1,A,I1)')
     >           " #     Number of scale variations for dimension ",
     >           J,": ",NSCALEVAR(I,J)
            WRITE(*,'(A,I1,A)')
     >           " #     Available scale settings for dimension ",
     >           J,":"
            DO K=1,NSCALEVAR(I,J)
               WRITE(*,'(A,I1,A,F10.4)')
     >              " #       Scale factor number ",
     >              K,":                   ",SCALEFAC(I,J,K)
            ENDDO
            WRITE(*,'(A,I1,A,I1)')
     >           " #     Number of scale nodes for dimension ",
     >           J,":      ",NSCALENODE(I,J)
         ENDDO
      ENDDO
      WRITE(*,*)"#"
      WRITE(*,*)CSEPS

      RETURN
      END
***********************************************************************

      SUBROUTINE FX9999TB(FNSTRING,IVAR,DVAR)
***********************************************************************
*
*     fastNLO user code v2.0 - returns information from table variables
*
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      CHARACTER*(*) FNSTRING
      INTEGER IVAR
      DOUBLE PRECISION DVAR

      IVAR = 0
      DVAR = 0D0
      IF (FNSTRING .EQ. 'NCONTRIB') THEN
         IVAR = NCONTRIB
      ELSEIF (FNSTRING .EQ. 'NOBSBIN') THEN
         IVAR = NOBSBIN
      ELSEIF (FNSTRING .EQ. 'NEVT1') THEN
         IVAR = NEVT(1)
      ELSEIF (FNSTRING .EQ. 'NEVT2') THEN
         IVAR = NEVT(2)
      ELSE
         WRITE(*,*)'FX9999TB: Unknown input variable:',fnstring
         STOP
      ENDIF

      RETURN
      END
***********************************************************************



      SUBROUTINE FX9999CL
***********************************************************************
*
*     fastNLO user code v2.0 - print contribution list
*
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INCLUDE 'strings.inc'
      INTEGER I,J,K
      CHARACTER*80 CHTMP

      WRITE(*,'(A)')
      WRITE(*,*)CSEPS
      WRITE(*,*)"# Overview on contribution types and "//
     >     "numbers contained in table:"
      WRITE(*,*)LSEPS
      WRITE(*,'(A,I2)')" # Number of contributions: ",Ncontrib
      DO I=ILO,IDATA
         WRITE(*,'(A,I2)')" #   No.: ",I
ckr TBD     >        I,ICONTRPOINTER(I),NCONTRCOUNTER(I),
ckr TBD     >        ICONTRSELECTOR(I)
      ENDDO
      WRITE(*,*)CSEPS

      RETURN
      END
***********************************************************************
