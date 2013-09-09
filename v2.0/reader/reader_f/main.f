      PROGRAM MAIN
***********************************************************************
*
*     fastNLO_reader: MAIN
*
*     M. Wobisch, K. Rabbertz, 2011
*
*     Program to read fastNLO v2 tables and derive
*     QCD cross sections using PDFs from LHAPDF
*
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INCLUDE 'strings.inc'
      CHARACTER*2  CH2TMP
      CHARACTER*14 CHTMP3
      CHARACTER*16 CHTMP1
      CHARACTER*18 CHTMP2
      CHARACTER*255 FILENAME,PDFSET,CHRES,CHFRM
      INTEGER I, J, IS, IPRINT, NDIMBINS(MXDIM), NBLANK
      INTEGER NSCDM, NSCLS
      INTEGER MXSCALECOMB
      PARAMETER (MXSCALECOMB=2*MXSCALEVAR)
      LOGICAL LTHC,LTHC1SEP,LTHC2SEP,LNPC1SEP
      DOUBLE PRECISION ALPS,FNALPHAS,SCALER,SCALEF
      DOUBLE PRECISION XMURS(MXSCALECOMB),XMUFS(MXSCALECOMB)
      DATA IPRINT/1/
*---  Define series of scale factor settings to test. Last and 8th entry
*---  is (0,0) and is not to be used!
      DATA XMURS/1.0D0,0.5D0,2.0D0,0.5D0,1.0D0,1.D0,2.D0,0.0D0/
      DATA XMUFS/1.0D0,0.5D0,2.0D0,1.0D0,0.5D0,2.D0,1.D0,0.0D0/

*---  ATTENTION: This is the most likely source of Fortran problems!
*---  For each scenario, the result array must be declared at least
*---  as large as in the definition in the common block of the
*---  corresponding scenario.
*---  See the scenario info printout "Tot. no. of observable bins"
*---  and compare to the value of the parameter MxObsBin
*---  in the file fnx9999.inc or [scenario].inc.
*---  Adapt the following to your scenario:
*---  Integer MxObsBin
*---  Parameter (MxObsBin = nnn)
*---  We recommend to name the array according to the scenario.
      DOUBLE PRECISION XSLO(MXOBSBIN),XSCLLO(MXOBSBIN)
      DOUBLE PRECISION XSNLO(MXOBSBIN),XSCLNLO(MXOBSBIN)
      DOUBLE PRECISION XSTHC1(MXOBSBIN),XSCLTHC1(MXOBSBIN)
      DOUBLE PRECISION XSTHC2(MXOBSBIN),XSCLTHC2(MXOBSBIN)
      DOUBLE PRECISION DXSUCTMP(MXOBSBIN,2),DXSCORTMP(MXOBSBIN,2)
      DOUBLE PRECISION XSNPC1(MXOBSBIN),XSCLNPC1(MXOBSBIN)
      DOUBLE PRECISION DXSUCNPC1(MXOBSBIN,2),DXSCORNPC1(MXOBSBIN,2)
      DOUBLE PRECISION XSDAT(MXOBSBIN),XSCLDAT(MXOBSBIN)
      DOUBLE PRECISION DXSUCDATA(MXOBSBIN,2),DXSCORDATA(MXOBSBIN,2)
      DOUBLE PRECISION KFAC(MXOBSBIN),KTHC1(MXOBSBIN)
      DOUBLE PRECISION KTHC2(MXOBSBIN),KNPC1(MXOBSBIN)

*---  Initialization
C---  *---  Fortran 90 functions for computing precision:
C---  *---  tiny(x), huge(x), precision(x)
C---  Write(*,*)"fnlo-fread: F90 double tiny, huge, precision = ",
C---  >     tiny(1d0),huge(1d0),precision(1d0)
      DO I=1,NOBSBIN
         XSLO(I)   = -1.D0
         XSNLO(I)  = -1.D0
         XSTHC1(I) = -1.D0
         XSTHC2(I) = -1.D0
         XSNPC1(I) = -1.D0
         KFAC(I)   =  0.D0
         KTHC1(I)  =  0.D0
         KTHC2(I)  =  0.D0
         KNPC1(I)  =  0.D0
      ENDDO
      LTHC1SEP = .FALSE.
      LTHC2SEP = .FALSE.
      LNPC1SEP = .FALSE.

*---  Parse command line
      WRITE(*,'(A)')
      WRITE(*,*)CSEPS
      WRITE(*,*)"# fnlo-read: Program Steering"
      WRITE(*,*)LSEPS
      IF (IARGC().LT.1) THEN
         FILENAME = "table.tab"
         WRITE(*,*)"# fnlo-read: WARNING! No table name given,"
         WRITE(*,*)"# taking the default table.tab instead!"
         WRITE(*,*)"#   For an explanation of command line "//
     >        "arguments type:"
         WRITE(*,*)"#   ./fnlo-fread -h"
      ELSE
         CALL GETARG(1,FILENAME)
         IF (FILENAME(1:LEN_TRIM(FILENAME)).EQ."-h") THEN
            WRITE(*,*)'#'
            WRITE(*,*)'# Usage: ./fnlo-fread [arguments]'
            WRITE(*,*)'# Table input file, def. = table.tab'
            WRITE(*,*)'# PDF set, def. = cteq6m.LHpdf'
            WRITE(*,*)'# Number of mu_r, mu_f scale settings to '//
     >           'investigate, if possible, def. = 1, max. = 7'
            WRITE(*,*)'#'
            WRITE(*,*)'# Give full path(s) if these are not in the cwd.'
            WRITE(*,*)'# Use "_" to skip changing a default argument.'
            WRITE(*,*)'#'
            STOP
         ELSEIF (FILENAME(1:1).EQ."_") THEN
            FILENAME = "table.tab"
            WRITE(*,*)
            WRITE(*,*)"# fnlo-read: WARNING! No table name given,"
            WRITE(*,*)"# taking the default table.tab instead!"
         ELSE
            WRITE(*,*)"# fnlo-read: Evaluating table: ",
     >           FILENAME(1:LEN_TRIM(FILENAME))
         ENDIF
      ENDIF

*---  PDF set
      PDFSET = "X"
      IF (IARGC().GE.2) THEN
         CALL GETARG(2,PDFSET)
      ENDIF
      IF (IARGC().LT.2.OR.PDFSET(1:1).EQ."_") THEN
         PDFSET = "cteq6m.LHpdf"
         WRITE(*,*)"# fnlo-read: WARNING! No PDF set given,"
         WRITE(*,*)"# taking cteq6m.LHpdf instead!"
      ELSE
         WRITE(*,*)"# fnlo-read: Using PDF set   : ",
     >        PDFSET(1:LEN_TRIM(PDFSET))
      ENDIF

*---  Number of scale settings
      CH2TMP = "X"
      IF (IARGC().GE.3) THEN
         CALL GETARG(3,CH2TMP)
      ENDIF
      IF (IARGC().LT.3.OR.CH2TMP(1:1).EQ."_") THEN
         NSCLS = 1
         WRITE(*,*)
     >        "# fnlo-read: No request given for number of "//
     >        "scale settings,"
         WRITE(*,*)
     >        "#            investigating primary scale only."
      ELSE
         READ(CH2TMP,'(I2)'),NSCLS
         IF (NSCLS.LT.1) THEN
            WRITE(*,*)
     >           "# fnlo-read: ERROR! No scale setting "//
     >           "or even less??? Aborting! NSCLS = ",
     >           NSCLS
            STOP
         ELSEIF (NSCLS.GT.MXSCALECOMB-1) THEN
            WRITE(*,*)
     >           "# fnlo-read: ERROR! Too many scale settings "//
     >           "requested, aborting! NSCLS = ",
     >           NSCLS
            STOP
         ELSE
            WRITE(*,'(A,I1,A)')
     >           " # fnlo-read: If possible, will try to do ",
     >           NSCLS, " scale setting(s)."
         ENDIF
      ENDIF

*---  Too many arguments
      IF (IARGC().GT.3) THEN
         WRITE(*,*)
     >        "fnlo-read: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF
      WRITE(*,*)CSEPS

*---  Print header
      Call FNLOFREAD()

*---  Initialize table
      Call FX9999IN(FILENAME)

*---  Determine pointers to contributions and print out contribution list
      Call FX9999PT(1D0,1D0,1)
      Call FX9999CL

*---  Print out scenario information
      Call FX9999NF

*---  Initialize LHAPDF
C---  CALL SETLHAPARM('SILENT')
      CALL SETLHAPARM('LOWKEY')
      CALL INITPDFSET(PDFSET(1:LEN_TRIM(PDFSET)))

*---  Initialize one member, 0=best fit member
      CALL INITPDF(0)

*---  Initial call to alpha_s interface
      ALPS = FNALPHAS(91.1876D0)

*---  Compute the cross sections
      WRITE(*,'(A)')""
      WRITE(*,'(A)')CSEPL
      WRITE(*,'(A)')"fnlo-read: Calculate my cross sections"
      WRITE(*,'(A)')CSEPL

*---  Loop over pre-defined scale settings in XMURS(F)
*---  For now assume only one scale dimension, since (MxScaleDim=1)!
      DO IS=1,NSCLS
         SCALER = XMURS(IS)
         SCALEF = XMUFS(IS)

*---  Initial settings
         Call FNSET("P_RESET",0) ! Reset all selections to zero

*---  Check on pointers to access contributions and scales
*---  If THC 1-loop not available --> only do LO
*---  If THC 2-loop not available --> only do NLO
*---  If NLO not available --> skip this scale setting
         IF (LCONTR(ILO)) CALL FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
         IF (LCONTR(ITHC1L)) THEN
            LTHC = .TRUE.
            CALL FNSET("P_THRESHCOR",1) ! select no. of loops in threshold correction
         ENDIF
         CALL FX9999PT(SCALER,SCALEF,IPRINT)
         IF (LCONTR(ITHC1L).AND..NOT.ISCALEPOINTER(ITHC1L).GT.0) THEN
            LTHC = .FALSE.
            CALL FNSET("P_THRESHCOR",0)
         ENDIF
         CALL FNSET("P_RESET",0) ! Reset all selections to zero
         IF (LCONTR(INLO)) CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
         IF (LCONTR(ITHC2L)) THEN
            LTHC = .TRUE.
            CALL FNSET("P_THRESHCOR",2) ! select no. of loops in threshold correction
         ENDIF
         CALL FX9999PT(SCALER,SCALEF,IPRINT)
         IF (LCONTR(INLO).AND..NOT.ISCALEPOINTER(INLO).GT.0) CYCLE
         IF (LCONTR(ITHC2L).AND..NOT.ISCALEPOINTER(ITHC2L).GT.0) THEN
            LTHC = .FALSE.
            CALL FNSET("P_THRESHCOR",0)
         ENDIF

*---  Reset all selections
         CALL FNSET("P_RESET",0) ! Reset all selections to zero

*---  In the following, set IPRINT to 1 for more verbose output
*---  Calculate LO cross sections
         IF (LCONTR(ILO)) THEN
            CALL FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
            CALL FX9999CC(SCALER, SCALEF, XSLO, XSCLLO,
     >           DXSUCTMP, DXSCORTMP)
         ENDIF

*---  Calculate NLO cross sections
         IF (LCONTR(ILO).AND.LCONTR(INLO)) THEN
            CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            CALL FX9999CC(SCALER, SCALEF, XSNLO, XSCLNLO,
     >           DXSUCTMP, DXSCORTMP)
         ENDIF

*---  Calculate LO cross section incl. 1-loop threshold corrections
*---  Only xmur = xmuf allowed here!
         IF (LCONTR(ILO).AND.LCONTR(ITHC1L)
     >        .AND.LTHC
     >        .AND..NOT.LTHC1SEP) THEN
            CALL FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
            CALL FNSET("P_THRESHCOR",1) ! select no. of loops in threshold correction
            CALL FX9999CC(SCALEF, SCALEF, XSTHC1, XSCLTHC1,
     >           DXSUCTMP, DXSCORTMP)
         ENDIF

*---  Print out 1-loop threshold corrections separately
         IF (LCONTR(ITHC1L)
     >        .AND.LTHC
     >        .AND.LTHC1SEP) THEN
            CALL FNSET("P_ORDPTHY",0) ! select order pert. theory: 1=LO, 2=NLO
            CALL FNSET("P_THRESHCOR",1) ! select no. of loops in thr. corr.
            CALL FX9999CC(SCALEF, SCALEF, XSTHC1, XSCLTHC1,
     >           DXSUCTMP,DXSCORTMP)
         ENDIF

*---  Calculate NLO cross section incl. 2-loop threshold corrections
*---  Only xmur = xmuf allowed here!
         IF (LCONTR(ILO).AND.LCONTR(INLO).AND.LCONTR(ITHC2L)
     >        .AND.LTHC
     >        .AND..NOT.LTHC2SEP) THEN
            CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            CALL FNSET("P_THRESHCOR",2) ! select no. of loops in threshold correction
            CALL FX9999CC(SCALEF, SCALEF, XSTHC2, XSCLTHC2,
     >           DXSUCTMP, DXSCORTMP)
         ENDIF

*---  Print out 2-loop threshold corrections separately
         IF (LCONTR(ITHC2L).AND.LTHC.AND.LTHC2SEP) THEN
            CALL FNSET("P_ORDPTHY",0) ! select order pert. theory: 1=LO, 2=NLO
            CALL FNSET("P_THRESHCOR",2) ! select no. of loops in thr. corr.
            CALL FX9999CC(SCALEF, SCALEF, XSTHC2, XSCLTHC2,
     >           DXSUCTMP,DXSCORTMP)
         ENDIF

*---  Apply non-perturbative corrections to NLO cross section
         IF (LCONTR(ILO).AND.LCONTR(INLO).AND.LCONTR(INPC1)
     >        .AND..NOT.LNPC1SEP) THEN
            CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            CALL FNSET("P_THRESHCOR",0) ! deselect threshold corrections
            CALL FNSET("P_NPCOR",1) ! select non-perturbative corrections
            CALL FX9999CC(SCALER, SCALEF, XSNPC1, XSCLNPC1,
     >           DXSUCNPC1, DXSCORNPC1)
         ENDIF

*---  Print out non-perturbative corrections separately
         IF (LCONTR(INPC1).AND.LNPC1SEP) THEN
            CALL FNSET("P_ORDPTHY",0) ! select order pert. theory: 1=LO, 2=NLO
            CALL FNSET("P_THRESHCOR",0) ! deselect threshold corrections
            CALL FNSET("P_NPCOR",1) ! select non-perturbative corrections
            CALL FX9999CC(SCALER, SCALEF, XSNPC1, XSCLNPC1,
     >           DXSUCNPC1,DXSCORNPC1)
            LNPC1SEP = .TRUE.
         ENDIF

*---  Reset selection for within scale loop
         CALL FNSET("P_RESET",0) ! Reset all selections to zero

*---  Cross section printout
         DO I=1,NOBSBIN
            IF (LCONTR(ILO).AND.LCONTR(INLO).AND.
     >           (ABS(XSLO(I)).GT.TINY(1D0))) THEN
               KFAC(I) = XSNLO(I) / XSLO(I)
            ENDIF
            IF (LCONTR(ITHC1L).AND.LTHC) THEN
               IF (LTHC1SEP) THEN
                  KTHC1(I) = XSTHC1(I)
               ELSE
                  IF ((ABS(XSLO(I)).GT.TINY(1D0))) THEN
                     KTHC1(I) = XSTHC1(I) / XSLO(I)
                  ENDIF
               ENDIF
            ENDIF
            IF (LCONTR(ITHC2L).AND.LTHC) THEN
               IF (LTHC2SEP) THEN
                  KTHC2(I) = XSTHC2(I)
               ELSE
                  IF ((ABS(XSNLO(I)).GT.TINY(1D0))) THEN
                     KTHC2(I) = XSTHC2(I) / XSNLO(I)
                  ENDIF
               ENDIF
            ENDIF
            IF (LNPC1SEP) THEN
               KNPC1(I) = XSNPC1(I)
            ELSE
               IF ((ABS(XSNLO(I)).GT.TINY(1D0))) THEN
                  KNPC1(I) = XSNPC1(I) / XSNLO(I)
               ENDIF
            ENDIF
         ENDDO
         WRITE(*,'(A)')DSEPL
         WRITE(*,'(A)')" My Cross Sections"
         WRITE(*,'(2(A,F10.3))')" The scale factors xmur, "//
     >        "xmuf chosen here are: ",SCALER,", ",SCALEF
         WRITE(*,'(A)')LSEPL
         CHTMP1 = DIMLABEL(1)
         CHTMP1 = "[ "//CHTMP1(1:12)//" ]"
         CHTMP2 = DIMLABEL(2)
         CHTMP2 = "[  "//CHTMP2(1:12)//"  ]"
         CHRES  = ""
         CHFRM  =
     >        "(1P,X,I5,X,G10.4,(X,I5,2(X,G10.4)),(X,I5,2(2X,E8.2)),0P"
         IF (LCONTR(ILO)) THEN
            CHTMP3  = SCALEDESCRIPT(ILO,1,1)
            CHTMP3  = "<"//CHTMP3(1:12)//">"
            CHRES = CHRES(1:LEN_TRIM(CHRES))//
     >           CHTMP3//"  "
            CHRES = CHRES(1:LEN_TRIM(CHRES))//
     >           "  LO cross section"
            CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//
     >           ",1P,(X,G10.4),5X,(X,E18.11),0P"
            IF (LCONTR(ITHC1L).AND.LTHC) THEN
               CHRES =
     >              CHRES(1:LEN_TRIM(CHRES))//
     >              "    KTHC1"
               CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//
     >              ",(X,F9.5)"
            ENDIF
            IF (LCONTR(INLO)) THEN
               CHRES = CHRES(1:LEN_TRIM(CHRES))/
     >              /"   NLO cross section"
               CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//
     >              ",1P,(X,E18.11),0P"
               CHRES = CHRES(1:LEN_TRIM(CHRES))//
     >              "   KNLO"
               CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//
     >              ",(X,F9.5)"
               CH2TMP = "  "
               NBLANK = 2
               IF (LCONTR(ITHC2L).AND.LTHC) THEN
                  CHRES =
     >                 CHRES(1:LEN_TRIM(CHRES))//
     >                 CH2TMP(1:NBLANK)//"    KTHC2"
                  NBLANK = 1
                  CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//
     >                 ",(X,F9.5)"
               ENDIF
               IF (LCONTR(INPC1)) THEN
                  CHRES =
     >                 CHRES(1:LEN_TRIM(CHRES))//
     >                 CH2TMP(1:NBLANK)//"    KNPC1"
                  CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//
     >                 ",(X,F9.5)"
               ENDIF
            ENDIF
         ENDIF
         CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//",X)"

         WRITE(*,'(A)')"  IObs  Bin Size "//
     >        "IODim1  "//
     >        CHTMP1//"    "//
     >        "IODim2  "//
     >        CHTMP2//"  "//CHRES(1:LEN_TRIM(CHRES))
         WRITE(*,'(A)')LSEPL

         DO I=1,NOBSBIN
            DO J=1,NDIM
               IF (I.EQ.1) THEN
                  NDIMBINS(J) = 1
               ELSEIF (LOBIN(I-1,J).LT.LOBIN(I,J)) THEN
                  NDIMBINS(J) = NDIMBINS(J) + 1
               ELSEIF (LOBIN(I,J).LT.LOBIN(I-1,J)) THEN
                  NDIMBINS(J) = 1
               ENDIF
            ENDDO

            IF (LCONTR(ILO).AND.LCONTR(INLO).AND.LCONTR(ITHC2L).AND.LTHC
     >           .AND.LCONTR(INPC1)) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
ckr     >              XSCLNPC1(I),
     >              XSCLNLO(I),
     >              XSLO(I),XSNLO(I),KFAC(I),KTHC2(I),KNPC1(I)
            ELSEIF (LCONTR(ILO).AND.LCONTR(INLO).AND.LCONTR(ITHC2L)
     >              .AND.LTHC) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
ckr     >              XSCLTHC2(I),
     >              XSCLNLO(I),
     >              XSLO(I),XSNLO(I),KFAC(I),KTHC2(I)
            ELSEIF (LCONTR(ILO).AND.LCONTR(INLO).AND.LCONTR(INPC1)) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
ckr     >              XSCLNPC1(I),
     >              XSCLNLO(I),
     >              XSLO(I),XSNLO(I),KFAC(I),KNPC1(I)
            ELSEIF (LCONTR(ILO).AND.LCONTR(INLO).AND.LCONTR(ITHC1L)
     >              .AND.LTHC) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
     >              XSCLNLO(I),
     >              XSLO(I),XSNLO(I),KFAC(I),KTHC1(I)
            ELSEIF (LCONTR(ILO).AND.LCONTR(INLO)) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
     >              XSCLNLO(I),
     >              XSLO(I),XSNLO(I),KFAC(I)
            ELSEIF (LCONTR(ILO).AND.LCONTR(ITHC1L).AND.LTHC) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
     >              XSCLLO(I),
     >              XSLO(I),KTHC1(I)
            ELSEIF (LCONTR(ILO)) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
     >              XSCLLO(I),
     >              XSLO(I)
            ELSE
               WRITE(*,*)
     >              "fnlo-read: Nothing to report!"
            ENDIF
         ENDDO
      ENDDO

*---  Print out data (set IPRINT to 1 for more verbose output)
      IF (LCONTR(IDATA)) THEN
         CALL FNSET("P_RESET",0) ! Reset all selections to zero
         CALL FNSET("P_DATA",1) ! Select data
         CALL FX9999CC(SCALER, SCALEF, XSDAT, XSCLDAT,
     >        DXSUCDATA, DXSCORDATA)

*---  Data section printout
         CHFRM  =
     >        "(1P,X,I5,X,G10.4,(X,I5,2(X,G10.4)),"//
     >        "(X,I5,2(2X,E8.2)),(3X,E10.3),SP,4(X,E9.2),X)"
         WRITE(*,'(A)')DSEPL
         WRITE(*,'(A)')" Measurement"
         WRITE(*,'(A)')LSEPL
         CHTMP1 = DIMLABEL(1)
         CHTMP1 = "[ "//CHTMP1(1:12)//" ]"
         CHTMP2 = DIMLABEL(2)
         CHTMP2 = "[  "//CHTMP2(1:12)//"  ]"
         WRITE(*,'(A)')"  IObs  Bin Size "//
     >        "IODim1  "//
     >        CHTMP1//"    "//
     >        "IODim2  "//
     >        CHTMP2//"    "//
     >        "X Section "//
     >        "QSumUnc.+ "//
     >        "QSumUnc.- "//
     >        "QSumCor.+ "//
     >        "QSumCor.-"
         WRITE(*,'(A)')LSEPL
         DO I=1,NOBSBIN
            DO J=1,NDIM
               IF (I.EQ.1) THEN
                  NDIMBINS(J) = 1
               ELSEIF (LOBIN(I-1,J).LT.LOBIN(I,J)) THEN
                  NDIMBINS(J) = NDIMBINS(J) + 1
               ELSEIF (LOBIN(I,J).LT.LOBIN(I-1,J)) THEN
                  NDIMBINS(J) = 1
               ENDIF
            ENDDO
            WRITE(*,CHFRM)I,BINSIZE(I),
     >           NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >           NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
     >           XSDAT(I),
     >           DXSUCDATA(I,2),DXSUCDATA(I,1),
     >           DXSCORDATA(I,2),DXSCORDATA(I,1)
         ENDDO
      ENDIF

      END
