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
      INTEGER NXMU
      PARAMETER(NXMU = 4)
      CHARACTER*16 CHTMP1
      CHARACTER*18 CHTMP2
      CHARACTER*255 FILENAME,PDFSET,CHRES,CHFRM
      INTEGER I, J, IS, IPRINT, NDIMBINS(MXDIM)
      INTEGER NSCDM, NSCLS, ISCLPNT(NXMU)
      LOGICAL LLO,LNLO,LTHC1L,LTHC2L,LNPC1,LDATA
      LOGICAL LTHCSEP,LNPCSEP
      DOUBLE PRECISION ALPS,FNALPHAS,SCALEF,XMU(NXMU)
      DATA ISCLPNT/0,0,0,0/
      DATA IPRINT/0/
*---  Define series of scale factor settings to test
      DATA XMU/1D0,0.25D0,0.5D0,2D0/

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
      DOUBLE PRECISION XSLO(MXOBSBIN) 
      DOUBLE PRECISION XSNLO(MXOBSBIN) 
      DOUBLE PRECISION XSTHC(MXOBSBIN) 
      DOUBLE PRECISION DXSUCTMP(MXOBSBIN,2),DXSCORTMP(MXOBSBIN,2) 
      DOUBLE PRECISION XSNPC(MXOBSBIN)
      DOUBLE PRECISION DXSUCNPC(MXOBSBIN,2),DXSCORNPC(MXOBSBIN,2) 
      DOUBLE PRECISION XSDAT(MXOBSBIN) 
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
            WRITE(*,*)'# PDF set, def. = cteq6mE.LHgrid'
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
         PDFSET = "cteq6mE.LHgrid"
         WRITE(*,*)"# fnlo-read: WARNING! No PDF set given,"
         WRITE(*,*)"# taking cteq6mE.LHgrid instead!"
      ELSE
         WRITE(*,*)"# fnlo-read: Using PDF set   : ",
     >        PDFSET(1:LEN_TRIM(PDFSET))
      ENDIF

*---  Too many arguments
      IF (IARGC().GT.2) THEN
         WRITE(*,*)
     >        "fnlo-read: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF
      WRITE(*,*)CSEPS

*---  Initialize table
      Call FX9999IN(FILENAME)

*---  Initial call to alpha_s interface
      ALPS = FNALPHAS(91.1876D0)

*---  Print out scenario information
      Call FX9999NF

*---  Check on existence of LO, NLO, and THC contributions
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

*---  Initialize LHAPDF  
      CALL SETLHAPARM('SILENT')
      CALL INITPDFSET(PDFSET(1:LEN_TRIM(PDFSET)))

*---  Initialize one member, 0=best fit member
      CALL INITPDF(0)

*---  Compute the cross sections
      WRITE(*,'(A)')""
      WRITE(*,'(A)')CSEPL
      WRITE(*,'(A)')"fnlo-read: Calculate cross sections"
      WRITE(*,'(A)')CSEPL

*---  Initial settings
      Call FNSET("P_RESET",0)   ! Reset all selections to zero
      
*---  Loop over allowed scale settings for values pre-set in XMU
*---  For now assume only one scale dimension, since (MxScaleDim=1)!
      DO IS=1,NXMU

*---  Check on pointers to access contributions and scales for NLO and
*---  2-loop --> otherwise skip this scale setting 
         IF (LNLO) CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
         IF (LTHC2L) CALL FNSET("P_THRESHCOR",2) ! select no. of loops in threshold correction
         CALL FX9999PT(XMU(IS),XMU(IS),0)
         IF (LNLO.AND..NOT.ISCALEPOINTER(INLO).GT.0) CYCLE
         IF (LTHC2L.AND..NOT.ISCALEPOINTER(ITHC2L).GT.0) CYCLE
         SCALEF = XMU(IS)
         CALL FNSET("P_RESET",0) ! Reset all selections to zero

*---  Calculate LO cross sections (set IPRINT to 1 for more verbose
*---  output)
         IF (LLO) THEN
            CALL FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
            CALL FX9999CC(SCALEF, SCALEF, XSLO, DXSUCTMP, DXSCORTMP)
         ENDIF
         
*---  Calculate NLO cross sections (set IPRINT to 1 for more verbose
*---  output)
         IF (LLO.AND.LNLO) THEN
            CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            CALL FX9999CC(SCALEF, SCALEF, XSNLO, DXSUCTMP, DXSCORTMP)
         ENDIF

*---  Calculate NLO cross section incl. 2-loop threshold corrections
*---  (set IPRINT to 1 for more verbose output)
         IF (LLO.AND.LNLO.AND.LTHC2L) THEN
            CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            CALL FNSET("P_THRESHCOR",2) ! select no. of loops in threshold correction
            CALL FX9999CC(SCALEF, SCALEF, XSTHC, DXSUCTMP, DXSCORTMP)
         ENDIF

*---  Print out 2-loop threshold corrections (set IPRINT to 1 for more
*---  verbose output)
C---  IF (LTHC2L) THEN
C---  CALL FNSET("P_ORDPTHY",0) ! select order pert. theory: 1=LO, 2=NLO
C---  CALL FNSET("P_THRESHCOR",2) ! select no. of loops in thr. corr.
C---  CALL FX9999CC(SCALEF, SCALEF, XSTHC, DXSUCTMP, DXSCORTMP)
C---  LTHCSEP = .TRUE.
C---  ENDIF

*---  Apply non-perturbative corrections to NLO cross section (set
*---  IPRINT to 1 for more verbose output)
         IF (LLO.AND.LNLO.AND.LNPC1) THEN
            CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            CALL FNSET("P_THRESHCOR",0) ! deselect threshold corrections
            CALL FNSET("P_NPCOR",1) ! select non-perturbative corrections
            CALL FX9999CC(SCALEF, SCALEF, XSNPC, DXSUCNPC, DXSCORNPC)
         ENDIF

*---  Print out non-perturbative corrections (set IPRINT to 1 for more
*---  verbose output)
C---  IF (LNPC1) THEN
C---  CALL FNSET("P_ORDPTHY",0) ! select order pert. theory: 1=LO, 2=NLO
C---  CALL FNSET("P_THRESHCOR",0) ! deselect threshold corrections
C---  CALL FNSET("P_NPCOR",1) ! select non-perturbative corrections
C---  CALL FX9999CC(SCALEF, SCALEF, XSNPC, DXSUCNPC, DXSCORNPC)
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
         WRITE(*,'(A)')" Cross Sections"
         WRITE(*,"(A,F10.3)")" The scale factor chosen here is: ",
     >        SCALEF
         WRITE(*,'(A)')LSEPL
         CHTMP1 = DIMLABEL(1)
         CHTMP1 = "[ "//CHTMP1(1:12)//" ]"
         CHTMP2 = DIMLABEL(2)
         CHTMP2 = "[  "//CHTMP2(1:12)//"  ]"
         CHRES  = ""
         CHFRM  =
     >        "(1P,X,I5,X,G10.4,(X,I5,2(X,G10.4)),(X,I5,2(2X,E8.2))"
         IF (LLO) THEN
            CHRES = CHRES(1:LEN_TRIM(CHRES))//"LO cross section"
            CHFRM = CHFRM(1:LEN_TRIM(CHFRM))//",(X,E18.11)"
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
     >              XSLO(I),XSNLO(I),KFAC(I),KTHC(I),KNPC(I)
            ELSEIF (LLO.AND.LNLO.AND.LTHC2L) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
     >              XSLO(I),XSNLO(I),KFAC(I),KTHC(I)
            ELSEIF (LLO.AND.LNLO.AND.LNPC1) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
     >              XSLO(I),XSNLO(I),KFAC(I),KNPC(I)
            ELSEIF (LLO.AND.LNLO) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
     >              XSLO(I),XSNLO(I),KFAC(I)
            ELSEIF (LLO) THEN
               WRITE(*,CHFRM)I,BINSIZE(I),
     >              NDIMBINS(1),LOBIN(I,1),UPBIN(I,1),
     >              NDIMBINS(2),LOBIN(I,2),UPBIN(I,2),
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
         CALL FX9999CC(SCALEF, SCALEF, XSDAT, DXSUCDATA, DXSCORDATA)

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
