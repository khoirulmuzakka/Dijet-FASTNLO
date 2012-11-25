      SUBROUTINE FNLOFREAD
***********************************************************************
*     
*     fastNLO_reader: FNLOFREAD
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
      CHARACTER*12 CHTMP3
      CHARACTER*16 CHTMP1
      CHARACTER*18 CHTMP2
      CHARACTER*255 FILENAME,PDFSET,CHRES,CHFRM
      INTEGER I, J, IS, IPRINT, NDIMBINS(MXDIM)
      INTEGER NSCDM, NSCLS
      INTEGER MXSCALECOMB
      INTEGER NCHPLUS
      PARAMETER (MXSCALECOMB=2*MXSCALEVAR)
      LOGICAL LLO,LNLO,LTHC1L,LTHC2L,LNPC1,LDATA
      LOGICAL LTHCSEP,LNPCSEP
      DOUBLE PRECISION ALPS,FNALPHAS,SCALER,SCALEF
      DOUBLE PRECISION XMURS(MXSCALECOMB),XMUFS(MXSCALECOMB)
      DATA IPRINT/0/
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
      DOUBLE PRECISION XSTHC(MXOBSBIN),XSCLTHC(MXOBSBIN)
      DOUBLE PRECISION DXSUCTMP(MXOBSBIN,2),DXSCORTMP(MXOBSBIN,2) 
      DOUBLE PRECISION XSNPC(MXOBSBIN),XSCLNPC(MXOBSBIN)
      DOUBLE PRECISION DXSUCNPC(MXOBSBIN,2),DXSCORNPC(MXOBSBIN,2) 
      DOUBLE PRECISION XSDAT(MXOBSBIN),XSCLDAT(MXOBSBIN)
      DOUBLE PRECISION DXSUCDATA(MXOBSBIN,2),DXSCORDATA(MXOBSBIN,2) 
      DOUBLE PRECISION KFAC(MXOBSBIN),KTHC(MXOBSBIN),KNPC(MXOBSBIN) 

*---  Initialization
C---  *---  Fortran 90 functions for computing precision:
C---  *---  tiny(x), huge(x), precision(x)
C---  Write(*,*)"fnlo-fread: F90 double tiny, huge, precision = ",
C---  >     tiny(1d0),huge(1d0),precision(1d0)
      DO I=1,NOBSBIN
         XSLO(I)  = -1.D0
         XSNLO(I) = -1.D0
         XSTHC(I) = -1.D0
         KFAC(I)  =  0.D0
         KTHC(I)  =  0.D0
         KNPC(I)  =  0.D0
      ENDDO
      LTHCSEP = .FALSE.
      LNPCSEP = .FALSE.

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

*---  Initialize table
      Call FX9999IN(FILENAME)

*---  Initial call to alpha_s interface
      ALPS = FNALPHAS(91.1876D0)

*---  Print out contribution list
      Call FX9999CL

*---  Print out scenario information
      Call FX9999NF

*---  Check on existence of LO, NLO, THC, NPC and DATA contributions
      LLO    = .FALSE.
      LNLO   = .FALSE.
      LTHC1L = .FALSE.
      LTHC2L = .FALSE.
      LNPC1  = .FALSE.
      LDATA  = .FALSE.
      DO I=1,NCONTRIB
         IF (ICONTRFLAG1(I).EQ.1.AND.ICONTRFLAG2(I).EQ.1)LLO    = .TRUE.
         IF (ICONTRFLAG1(I).EQ.1.AND.ICONTRFLAG2(I).EQ.2)LNLO   = .TRUE.
         IF (ICONTRFLAG1(I).EQ.2.AND.ICONTRFLAG2(I).EQ.1)LTHC1L = .TRUE.
         IF (ICONTRFLAG1(I).EQ.2.AND.ICONTRFLAG2(I).EQ.2)LTHC2L = .TRUE.
         IF (ICONTRFLAG1(I).EQ.4.AND.ICONTRFLAG2(I).EQ.1.AND.
     >        IADDMULTFLAG(I).EQ.1)
     >        LNPC1 = .TRUE.
         IF (ICONTRFLAG1(I).EQ.0.AND.ICONTRFLAG2(I).EQ.0.AND.
     >        IDATAFLAG(I).EQ.1)
     >        LDATA = .TRUE.
      ENDDO
ckr TBD
      lthc1l = .false.
      lthc2l = .false.
      lnpc1 = .false.
      ldata = .false.

*---  Initialize LHAPDF  
C---  CALL SETLHAPARM('SILENT')
      CALL SETLHAPARM('LOWKEY')
      CALL INITPDFSET(PDFSET(1:LEN_TRIM(PDFSET)))

*---  Initialize one member, 0=best fit member
      CALL INITPDF(0)

*---  Compute the cross sections
      WRITE(*,'(A)')""
      WRITE(*,'(A)')CSEPL
      WRITE(*,'(A)')"fnlo-read: Calculate my cross sections"
      WRITE(*,'(A)')CSEPL

*---  Initial settings
      Call FNSET("P_RESET",0)   ! Reset all selections to zero
      
*---  Loop over pre-defined scale settings in XMURS(F)
*---  For now assume only one scale dimension, since (MxScaleDim=1)!
      DO IS=1,NSCLS
         SCALER = XMURS(IS)
         SCALEF = XMUFS(IS)

*---  Check on pointers to access contributions and scales for NLO and
*---  2-loop --> otherwise skip this scale setting 
         IF (LNLO) CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
         IF (LTHC2L) CALL FNSET("P_THRESHCOR",2) ! select no. of loops in threshold correction
         CALL FX9999PT(SCALER,SCALEF,IPRINT)
         IF (LNLO.AND..NOT.ISCALEPOINTER(INLO).GT.0) CYCLE
         IF (LTHC2L.AND..NOT.ISCALEPOINTER(ITHC2L).GT.0) CYCLE
         CALL FNSET("P_RESET",0) ! Reset all selections to zero

*---  Calculate LO cross sections (set IPRINT to 1 for more verbose
*---  output)
         IF (LLO) THEN
            CALL FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
            CALL FX9999CC(SCALER, SCALEF, XSLO, XSCLLO,
     >           DXSUCTMP, DXSCORTMP)
         ENDIF
         
*---  Calculate NLO cross sections (set IPRINT to 1 for more verbose
*---  output)
         IF (LLO.AND.LNLO) THEN
            CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            CALL FX9999CC(SCALER, SCALEF, XSNLO, XSCLNLO,
     >           DXSUCTMP, DXSCORTMP)
         ENDIF

*---  Calculate NLO cross section incl. 2-loop threshold corrections
*---  Only xmur = xmuf allowed here!
*---  (set IPRINT to 1 for more verbose output)
         IF (LLO.AND.LNLO.AND.LTHC2L) THEN
            CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            CALL FNSET("P_THRESHCOR",2) ! select no. of loops in threshold correction
            CALL FX9999CC(SCALEF, SCALEF, XSTHC, XSCLTHC,
     >           DXSUCTMP, DXSCORTMP)
         ENDIF

*---  Print out 2-loop threshold corrections (set IPRINT to 1 for more
*---  verbose output)
C---  IF (LTHC2L) THEN
C---  CALL FNSET("P_ORDPTHY",0) ! select order pert. theory: 1=LO, 2=NLO
C---  CALL FNSET("P_THRESHCOR",2) ! select no. of loops in thr. corr.
C---  CALL FX9999CC(SCALEF, SCALEF, XSTHC, XSCLTHC, DXSUCTMP, DXSCORTMP)
C---  LTHCSEP = .TRUE.
C---  ENDIF

*---  Apply non-perturbative corrections to NLO cross section (set
*---  IPRINT to 1 for more verbose output)
         IF (LLO.AND.LNLO.AND.LNPC1) THEN
            CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            CALL FNSET("P_THRESHCOR",0) ! deselect threshold corrections
            CALL FNSET("P_NPCOR",1) ! select non-perturbative corrections
            CALL FX9999CC(SCALER, SCALEF, XSNPC, XSCLNPC,
     >           DXSUCNPC, DXSCORNPC)
         ENDIF

*---  Print out non-perturbative corrections (set IPRINT to 1 for more
*---  verbose output)
C---  IF (LNPC1) THEN
C---  CALL FNSET("P_ORDPTHY",0) ! select order pert. theory: 1=LO, 2=NLO
C---  CALL FNSET("P_THRESHCOR",0) ! deselect threshold corrections
C---  CALL FNSET("P_NPCOR",1) ! select non-perturbative corrections
C---  CALL FX9999CC(SCALER, SCALEF, XSNPC, XSCLNPC, DXSUCNPC, DXSCORNPC)
C---  LNPCSEP = .TRUE.
C---  ENDIF

*---  Reset selection for within scale loop
         CALL FNSET("P_RESET",0) ! Reset all selections to zero

*---  Cross section printout
         DO I=1,NOBSBIN
            IF ((ABS(XSLO(I)).GT.TINY(1D0))) THEN
               KFAC(I) = XSNLO(I) / XSLO(I)
            ENDIF
            IF (LTHCSEP) THEN
               KTHC(I) = XSTHC(I)
            ELSE
               IF ((ABS(XSNLO(I)).GT.TINY(1D0))) THEN
                  KTHC(I) = XSTHC(I) / XSNLO(I)
               ENDIF
            ENDIF
            IF (LNPCSEP) THEN
               KNPC(I) = XSNPC(I)
            ELSE
               IF ((ABS(XSNLO(I)).GT.TINY(1D0))) THEN
                  KNPC(I) = XSNPC(I) / XSNLO(I)
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
     >        "(1P,X,I5,X,G10.4,(X,I5,2(X,G10.4)),(X,I5,2(2X,E8.2))"
         IF (LLO) THEN
            CHTMP3  = SCALEDESCRIPT(ILO,1,1)
            CHTMP3  = "<"//CHTMP3(1:10)//">"
            CHRES = CHRES(1:LEN_TRIM(CHRES))//
     >           CHTMP3//"  "
            CHRES = CHRES(1:LEN_TRIM(CHRES))//"  LO cross section"
            CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//"(X,G10.4),3X,(X,E18.11)"
            IF (LNLO) THEN
               CHRES = CHRES(1:LEN_TRIM(CHRES))/
     >              /"   NLO cross section"
               CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//",(X,E18.11)"
               CHRES = CHRES(1:LEN_TRIM(CHRES))//"   K NLO"
               CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//",0P,(X,F9.5)"
               IF (LTHC2L) THEN
                  CHRES =
     >                 CHRES(1:LEN_TRIM(CHRES))//"     K THC"
                  CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//",(X,F9.5)"
               ENDIF
               IF (LNPC1) THEN
                  CHRES =
     >                 CHRES(1:LEN_TRIM(CHRES))//"     K NPC"
                  CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//",(X,F9.5)"
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

            IF (LLO.AND.LNLO.AND.LTHC2L.AND.LNPC1) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
     >              XSCLNPC(I),
     >              XSLO(I),XSNLO(I),KFAC(I),KTHC(I),KNPC(I)
            ELSEIF (LLO.AND.LNLO.AND.LTHC2L) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
     >              XSCLTHC(I),
     >              XSLO(I),XSNLO(I),KFAC(I),KTHC(I)
            ELSEIF (LLO.AND.LNLO.AND.LNPC1) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
     >              XSCLNPC(I),
     >              XSLO(I),XSNLO(I),KFAC(I),KNPC(I)
            ELSEIF (LLO.AND.LNLO) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
     >              XSCLNLO(I),
     >              XSLO(I),XSNLO(I),KFAC(I)
            ELSEIF (LLO) THEN
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
      IF (LDATA) THEN
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
      
      RETURN
      END
