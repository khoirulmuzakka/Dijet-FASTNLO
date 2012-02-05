      SUBROUTINE FNLOFREAD
*********************************************************************
*     
*     fastNLO_reader: FNLOFREAD
*     Program to read fastNLO v2 tables and derive
*     QCD cross sections using PDFs from LHAPDF
*     
*     M. Wobisch, K. Rabbertz
*     
*********************************************************************
      Implicit None
      Include 'fnx9999.inc'
      Include 'strings.inc'
      Integer NXMU
      Parameter(NXMU = 4)
      Character*16 CHTMP1,CHTMP2
      Character*255 FILENAME,PDFSET,CHRES,CHFRM
      Integer i, j, IS, IPRINT, NDimBins(MxDim)
      Integer NSCDM, NSCLS, ISCLPNT(NXMU)
      Logical LLO,LNLO,LTHC1L,LTHC2L,LNPC1,LDATA
      Logical LTHCSEP,LNPCSEP
      Double Precision ALPS,FNALPHAS,SCALEF,XMU(NXMU)
      Data ISCLPNT/0,0,0,0/
      Data XMU/1D0,2D0,0.5D0,0.25D0/
      Data IPRINT/0/

*---  Attention - this is the most likely source of Fortran errors!
*---  For each scenario, the result array must be declared at least  
*---  as large as in the definition in the common block of the
*---  corresponding scenario. 
*---  --> See the value of the parameter MxObsBin
*---  in the file [scenario].inc or fnx9999.inc
*---  We recommend to name the array according to the scenario
*---  Adapt the following to your scenario!
*---  Integer MxObsBin
*---  Parameter (MxObsBin = 200)
      Double Precision xslo(MxObsBin) 
      Double Precision xsnlo(MxObsBin) 
      Double Precision xsthc(MxObsBin) 
      Double Precision dxsuctmp(MxObsBin,2),dxscortmp(MxObsBin,2) 
      Double Precision xsnpc(MxObsBin)
      Double Precision dxsucnpc(MxObsBin,2),dxscornpc(MxObsBin,2) 
      Double Precision xsdat(MxObsBin) 
      Double Precision dxsucdata(MxObsBin,2),dxscordata(MxObsBin,2) 
      Double Precision kfac(MxObsBin),kthc(MxObsBin),knpc(MxObsBin) 

*---  Initialization
Comment: *---  Fortran 90 functions for computing precision:
Comment: *---  tiny(x), huge(x), precision(x)
Comment:       Write(*,*)"fnlo-fread: F90 double tiny, huge, precision = ",
Comment:      >     tiny(1d0),huge(1d0),precision(1d0)
      DO I=1,NObsBin
         xslo(I)  = -1.d0
         xsnlo(I) = -1.d0
         xsthc(I) = -1.d0
         kfac(I)  =  0.d0
         kthc(I)  =  0.d0
         knpc(I)  =  0.d0
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
      DO I=1,NContrib
         IF (IContrFlag1(I).EQ.1.AND.IContrFlag2(I).EQ.1)LLO    = .TRUE.
         IF (IContrFlag1(I).EQ.1.AND.IContrFlag2(I).EQ.2)LNLO   = .TRUE.
         IF (IContrFlag1(I).EQ.2.AND.IContrFlag2(I).EQ.1)LTHC1L = .TRUE.
         IF (IContrFlag1(I).EQ.2.AND.IContrFlag2(I).EQ.2)LTHC2L = .TRUE.
         IF (IContrFlag1(I).EQ.4.AND.IContrFlag2(I).EQ.1.AND.
     >        IAddMultFlag(I).EQ.1)
     >        LNPC1 = .TRUE.
         IF (IContrFlag1(I).EQ.0.AND.IContrFlag2(I).EQ.0.AND.
     >        IDataFlag(I).EQ.1)
     >        LDATA = .TRUE.
      ENDDO

*---  Initialize LHAPDF  
      Call SetLHAPARM('SILENT')
      Call InitPDFset(PDFSET(1:LEN_TRIM(PDFSET)))

*---  Initialize one member, 0=best fit member
      call InitPDF(0)

*---  Compute the cross sections
      WRITE(*,'(A)')""
      WRITE(*,'(A)')CSEPL
      WRITE(*,'(A)')"fnlo-read: Calculate cross sections"
      WRITE(*,'(A)')CSEPL

*---  Initial settings
      Call FNSET("P_RESET",0)   ! Reset all selections to zero
      
*---  Loop over allowed scale settings for values required in XMU
*---  For now assume only one scale dimension, since (MxScaleDim=1)!
      DO IS=1,NXMU

*---  Check on pointers to access contributions and scales for NLO and
*---  2-loop --> otherwise skip this scale setting 
         If (LNLO) Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
         If (LTHC2L) Call FNSET("P_THRESHCOR",2) ! select no. of loops in threshold correction
         Call FX9999PT(XMU(IS),XMU(IS),0)
         IF (LNLO.AND..NOT.IScalePointer(INLO).GT.0) CYCLE
         If (LTHC2L.AND..NOT.IScalePointer(ITHC2L).GT.0) CYCLE
         SCALEF = XMU(IS)
         Call FNSET("P_RESET",0) ! Reset all selections to zero

*---  Calculate LO cross sections (set IPRINT to 1 for more verbose
*---  output)
         IF (LLO) THEN
            Call FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
            Call FX9999CC(SCALEF, SCALEF, XSLO, DXSUCTMP, DXSCORTMP)
         ENDIF
         
*---  Calculate NLO cross sections (set IPRINT to 1 for more verbose
*---  output)
         IF (LLO.AND.LNLO) THEN
            Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            Call FX9999CC(SCALEF, SCALEF, XSNLO, DXSUCTMP, DXSCORTMP)
         ENDIF

*---  Calculate NLO cross section incl. 2-loop threshold corrections
*---  (set IPRINT to 1 for more verbose output)
         IF (LLO.AND.LNLO.AND.LTHC2L) THEN
            Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            Call FNSET("P_THRESHCOR",2) ! select no. of loops in threshold correction
            Call FX9999CC(SCALEF, SCALEF, XSTHC, DXSUCTMP, DXSCORTMP)
         ENDIF

*---  Print out 2-loop threshold corrections (set IPRINT to 1 for more
*---  verbose output)
Comment:          IF (LTHC2L) THEN
Comment:             Call FNSET("P_ORDPTHY",0) ! select order pert. theory: 1=LO, 2=NLO
Comment:             Call FNSET("P_THRESHCOR",2) ! select no. of loops in threshold correction
Comment:             Call FX9999CC(SCALEF, SCALEF, XSTHC, DXSUCTMP, DXSCORTMP)
Comment:             LTHCSEP = .TRUE.
Comment:          ENDIF

*---  Apply non-perturbative corrections to NLO cross section (set
*---  IPRINT to 1 for more verbose output)
         IF (LLO.AND.LNLO.AND.LNPC1) THEN
            Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            Call FNSET("P_THRESHCOR",0) ! deselect threshold corrections
            Call FNSET("P_NPCOR",1) ! select non-perturbative corrections
            Call FX9999CC(SCALEF, SCALEF, XSNPC, DXSUCNPC, DXSCORNPC)
         ENDIF

*---  Print out non-perturbative corrections (set IPRINT to 1 for more
*---  verbose output)
Comment:          IF (LNPC1) THEN
Comment:             Call FNSET("P_ORDPTHY",0) ! select order pert. theory: 1=LO, 2=NLO
Comment:             Call FNSET("P_THRESHCOR",0) ! deselect threshold corrections
Comment:             Call FNSET("P_NPCOR",1) ! select non-perturbative corrections
Comment:             Call FX9999CC(SCALEF, SCALEF, XSNPC, DXSUCNPC, DXSCORNPC)
Comment:             LNPCSEP = .TRUE.
Comment:          ENDIF

*---  Reset selection for within scale loop
         Call FNSET("P_RESET",0) ! Reset all selections to zero

*---  Cross section printout
         DO I=1,NObsBin
            IF ((ABS(xslo(I)).GT.Tiny(1D0))) THEN
               kfac(I) = xsnlo(I) / xslo(I)
            ENDIF
            If (LTHCSEP) Then
               kthc(I) = xsthc(I)
            Else
               IF ((ABS(xsnlo(I)).GT.Tiny(1D0))) THEN
                  kthc(I) = xsthc(I) / xsnlo(I)
               ENDIF
            Endif
            If (LNPCSEP) Then
               knpc(I) = xsnpc(I)
            Else
               IF ((ABS(xsnlo(I)).GT.Tiny(1D0))) THEN
                  knpc(I) = xsnpc(I) / xsnlo(I)
               ENDIF
            Endif
         ENDDO
         WRITE(*,'(A)')DSEPL
         WRITE(*,'(A)')" Cross Sections"
         WRITE(*,"(A,F10.3)")" The scale factor chosen here is: ",
     >        SCALEF
         WRITE(*,'(A)')LSEPL
         CHTMP1 = DimLabel(1)
         CHTMP1 = "[ "//CHTMP1(1:12)//" ]"
         CHTMP2 = DimLabel(2)
         CHTMP2 = "[ "//CHTMP2(1:12)//" ]"
         CHRES  = ""
         CHFRM  =
     >        "(1P,X,I5,X,G10.4,(X,I5,2(X,G10.4)),(X,I5,2(2X,E7.1))"
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

         DO I=1,NObsBin
            DO J=1,NDim
               IF (I.EQ.1) THEN
                  NDimBins(J) = 1 
               ELSEIF (LoBin(I-1,J).LT.LoBin(I,J)) THEN
                  NDimBins(J) = NDimBins(J) + 1 
               ELSEIF (LoBin(I,J).LT.LoBin(I-1,J)) THEN
                  NDimBins(J) = 1 
               ENDIF
            ENDDO

            IF (LLO.AND.LNLO.AND.LTHC2L.AND.LNPC1) THEN
               WRITE(*,CHFRM)I,BinSize(I),
     >              NDimBins(1),LoBin(I,1),UpBin(I,1),
     >              NDimBins(2),LoBin(I,2),UpBin(I,2),
     >              XSLO(I),XSNLO(I),KFAC(I),KTHC(I),KNPC(I)
            ELSEIF (LLO.AND.LNLO.AND.LTHC2L) THEN
               WRITE(*,CHFRM)I,BinSize(I),
     >              NDimBins(1),LoBin(I,1),UpBin(I,1),
     >              NDimBins(2),LoBin(I,2),UpBin(I,2),
     >              XSLO(I),XSNLO(I),KFAC(I),KTHC(I)
            ELSEIF (LLO.AND.LNLO.AND.LNPC1) THEN
               WRITE(*,CHFRM)I,BinSize(I),
     >              NDimBins(1),LoBin(I,1),UpBin(I,1),
     >              NDimBins(2),LoBin(I,2),UpBin(I,2),
     >              XSLO(I),XSNLO(I),KFAC(I),KNPC(I)
            ELSEIF (LLO.AND.LNLO) THEN
               WRITE(*,CHFRM)I,BinSize(I),
     >              NDimBins(1),LoBin(I,1),UpBin(I,1),
     >              NDimBins(2),LoBin(I,2),UpBin(I,2),
     >              XSLO(I),XSNLO(I),KFAC(I)
            ELSEIF (LLO) THEN
               WRITE(*,CHFRM)I,BinSize(I),
     >              NDimBins(1),LoBin(I,1),UpBin(I,1),
     >              NDimBins(2),LoBin(I,2),UpBin(I,2),
     >              XSLO(I)
            ELSE
               WRITE(*,*)
     >              "fnlo-read: Nothing to report!"
            ENDIF
            
         ENDDO
      ENDDO

*---  Print out data (set IPRINT to 1 for more verbose output)
      IF (LDATA) THEN
         Call FNSET("P_RESET",0) ! Reset all selections to zero
         Call FNSET("P_DATA",1) ! Select data
         Call FX9999CC(SCALEF, SCALEF, XSDAT, DXSUCDATA, DXSCORDATA)

*---  Data section printout
         WRITE(*,'(A)')DSEPL
         WRITE(*,'(A)')" Measurement"
         WRITE(*,'(A)')LSEPL
         CHTMP1 = DimLabel(1)
         CHTMP1 = "[ "//CHTMP1(1:12)//" ]"
         CHTMP2 = DimLabel(2)
         CHTMP2 = "[ "//CHTMP2(1:12)//" ]"
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
         DO I=1,NObsBin
            DO J=1,NDim
               IF (I.EQ.1) THEN
                  NDimBins(J) = 1 
               ELSEIF (LoBin(I-1,J).LT.LoBin(I,J)) THEN
                  NDimBins(J) = NDimBins(J) + 1 
               ELSEIF (LoBin(I,J).LT.LoBin(I-1,J)) THEN
                  NDimBins(J) = 1 
               ENDIF
            ENDDO
            WRITE(*,999)I,BinSize(I),
     >           NDimBins(1),LoBin(I,1),UpBin(I,1),
     >           NDimBins(2),LoBin(I,2),UpBin(I,2),
     >           XSDAT(I),
     >           DXSUCDATA(I,2),DXSUCDATA(I,1),
     >           DXSCORDATA(I,2),DXSCORDATA(I,1)
         ENDDO
      Endif
      
 999  FORMAT(1P,X,I5,X,G10.4,(X,I5,2(X,G10.4)),
     >     (X,I5,2(2X,E7.1)),(3X,E10.3),SP,4(X,E9.2),X)
      
      Return
      End
