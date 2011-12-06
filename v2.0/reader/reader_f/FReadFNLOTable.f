      PROGRAM FReadFNLOTable
*     ------------------------------------------------------------------
*     
*     M. Wobisch                                08/10/2010
*     
*     fastNLO_reader:
*     Program to read fastNLO v2 tables and derive
*     QCD cross sections using PDFs from LHAPDF
*     
*     K. Rabbertz                               16.11.2011
*     
*     First version common with C++ reader
*     
*     ------------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Character*16 CHTMP1,CHTMP2
      Character*37  CSEP37,DSEP37,LSEP37,SSEP37
      Character*74  CSEPS,DSEPS,LSEPS,SSEPS
      Character*148 CSEPL,DSEPL,LSEPL,SSEPL
      Character*255 FILENAME,PDFSET
      Integer i, j, IS, IPRINT, NDimBins(MxDim)
      Logical LLO,LNLO,LTHC1L,LTHC2L,LMULT,LDATA
      Double Precision SCALEF
      Data IPRINT/0/
      Data CSEP37,DSEP37,LSEP37,SSEP37/
     >     '#####################################',
     >     '=====================================',
     >     "-------------------------------------",
     >     "*************************************"/

c     - Attention - this is the most likely source of Fortran errors in !!!
c     fastNLO
c     For each scenario, the result array must be declared at least  
c     as large as in the definition in the common block of the
c     corresponding scenario. 
c     -> See the value of the parameter MxObsBin
c     in the file [scenario].inc
c     We recommend to name the array according to the scenario
c     Adapt the following to your scenario!
c     Integer MxObsBin
c     Parameter (MxObsBin = 200)
      Double Precision xslo(MxObsBin) 
      Double Precision xsnlo(MxObsBin) 
      Double Precision xsthc(MxObsBin) 
      Double Precision xsnpc(MxObsBin)
      Double Precision dxsucnpc(MxObsBin,2),dxscornpc(MxObsBin,2) 
      Double Precision xsdat(MxObsBin) 
      Double Precision dxsucdata(MxObsBin,2),dxscordata(MxObsBin,2) 
      Double Precision kfac(MxObsBin),kthc(MxObsBin),knpc(MxObsBin) 

*---  Initialization
      CSEPL = CSEP37//CSEP37//CSEP37//CSEP37
      DSEPL = DSEP37//DSEP37//DSEP37//DSEP37
      LSEPL = LSEP37//LSEP37//LSEP37//LSEP37
      SSEPL = SSEP37//SSEP37//SSEP37//SSEP37
      CSEPS = CSEP37//CSEP37
      DSEPS = DSEP37//DSEP37
      LSEPS = LSEP37//LSEP37
      SSEPS = SSEP37//SSEP37
Comment:       CSEPS(1:2) = "# "
      DSEPS(1:2) = "# "
      LSEPS(1:2) = "# "
      SSEPS(1:2) = "# "
Comment:       CSEPL(1:2) = "# "
Comment:       DSEPL(1:2) = "# "
Comment:       LSEPL(1:2) = "# "
Comment:       SSEPL(1:2) = "# "
ckr CHeck Fortran 90 functions tiny(x),huge(x),precision(x) 
      write(*,*)"tiny,huge,precision",tiny(1d0),huge(1d0),precision(1d0)
      DO I=1,NObsBin
         xslo(I)  = -1.d0
         xsnlo(I) = -1.d0
         xsthc(I) = -1.d0
         kfac(I)  =  0.d0
         kthc(I)  =  0.d0
      ENDDO

*---  Initial output
      WRITE(*,*)""
      WRITE(*,*)CSEPS
      WRITE(*,*)"# ReadFNLOTable"
      WRITE(*,*)CSEPS
      WRITE(*,*)"# Program to read fastNLO v2 tables and derive"
      WRITE(*,*)"# QCD cross sections using PDFs from LHAPDF"
      WRITE(*,*)CSEPS

*---  Parse command line
      WRITE(*,*)"# ReadFNLOTable: Program Steering"
      WRITE(*,*)LSEPS
      IF (IARGC().LT.1) THEN
         FILENAME = "table.tab"
         WRITE(*,*)"# ReadFNLOTable: WARNING! No table name given,"
         WRITE(*,*)"# taking the default table.tab instead!"
         WRITE(*,*)"#   For an explanation of command line "//
     >        "arguments type:"
         WRITE(*,*)"#   ./FReadFNLOTable -h"
      ELSE
         CALL GETARG(1,FILENAME)
         IF (FILENAME(1:LEN_TRIM(FILENAME)).EQ."-h") THEN
            WRITE(*,*)'#'
            WRITE(*,*)'# Usage: ./FReadFNLOTable [arguments]'
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
            WRITE(*,*)"# ReadFNLOTable: WARNING! No table name given,"
            WRITE(*,*)"# taking the default table.tab instead!"
         ELSE
            WRITE(*,*)"# ReadFNLOTable: Evaluating table: ",
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
         WRITE(*,*)"# ReadFNLOTable: WARNING! No PDF set given,"
         WRITE(*,*)"# taking cteq6mE.LHgrid instead!"
      ELSE
         WRITE(*,*)"# ReadFNLOTable: Using PDF set   : ",
     >        PDFSET(1:LEN_TRIM(PDFSET))
      ENDIF

*---  Too many arguments
      IF (IARGC().GT.2) THEN
         WRITE(*,*)
     >        "ReadFNLOTable: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF

*---  Initialize table
      Call FX9999IN(FILENAME)

*---  Print out scenario information
      Call FX9999NF

*---  Check on existence of LO, NLO, and THC contributions
      LLO    = .FALSE.
      LNLO   = .FALSE.
      LTHC1L = .FALSE.
      LTHC2L = .FALSE.
      LMULT  = .FALSE.
      LDATA  = .FALSE.
      DO I=1,NContrib
         IF (IContrFlag1(I).EQ.1.AND.IContrFlag2(I).EQ.1) LLO  = .TRUE.
         IF (IContrFlag1(I).EQ.1.AND.IContrFlag2(I).EQ.2) LNLO = .TRUE.
         IF (IContrFlag1(I).EQ.2.AND.IContrFlag2(I).EQ.1.AND.
     >        IContrFlag3(I).EQ.1) LTHC1L = .TRUE.
         IF (IContrFlag1(I).EQ.2.AND.IContrFlag2(I).EQ.1.AND.
     >        IContrFlag3(I).EQ.2) LTHC2L = .TRUE.
         IF (IContrFlag1(I).EQ.0.AND.IContrFlag2(I).EQ.0.AND.
     >        IContrFlag3(I).EQ.0.AND.IAddMultFlag(I).EQ.1)
     >        LMULT = .TRUE.
         IF (IContrFlag1(I).EQ.0.AND.IContrFlag2(I).EQ.0.AND.
     >        IContrFlag3(I).EQ.0.AND.IDataFlag(I).EQ.1)
     >        LDATA = .TRUE.
      ENDDO

*---  Initialize LHAPDF  
      Call SetLHAPARM('SILENT')
      Call InitPDFset(PDFSET(1:LEN_TRIM(PDFSET)))

*---  Initialize one member, 0=best fit member
      call InitPDF(0)

*---  Compute the cross sections
      WRITE(*,*)""
      WRITE(*,'(A)')CSEPL
      WRITE(*,'(A)')"ReadFNLOTable: Calculate cross sections"
      WRITE(*,'(A)')CSEPL

*---  Initial settings
      Call FNSET("P_REFTAB",0)  ! evaluate standard table: 0, or reference: 1
      Call FNSET("P_ORDPTHY",0) ! select order pert. theory: 1=LO, 2=NLO
      Call FNSET("P_THRESHCOR",0) ! select no. loops in threshold
      Call FNSET("P_NPCOR",0)   ! deselect non-perturbative corrections
      Call FNSET("P_DATA",0)    ! deselect data
         
*---  Loop over scale settings in order of appearance in the table
*---  Assume for now that NLO = ctrb. 2 and scale dimension = 1
      DO IS=1,NScaleVar(2,1)
         SCALEF = ScaleFac(2,1,IS)
         
*---  Calculate LO cross sections (set IPRINT to 1 for more verbose
*---  output)
         IF (LLO) THEN
            Call FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
            Call FX9999CC(SCALEF, SCALEF, XSLO, DXSUCNPC, DXSCORNPC)
         ENDIF
         
*---  Calculate NLO cross sections (set IPRINT to 1 for more verbose
*---  output)
         IF (LLO.AND.LNLO) THEN
            Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
            Call FX9999CC(SCALEF, SCALEF, XSNLO, DXSUCNPC, DXSCORNPC)
         ENDIF

*---  Calculate NLO cross section incl. 2-loop threshold corrections
*---  (set IPRINT to 1 for more verbose output)
Comment:          IF (LLO.AND.LNLO.AND.LTHC2L) THEN
Comment:             Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
Comment:             Call FNSET("P_THRESHCOR",2) ! select no. of loops in threshold correction
Comment:             Call FX9999CC(SCALEF, SCALEF, XSTHC, DXSUCNPC, DXSCORNPC)
Comment:          ENDIF

*---  Print out 2-loop threshold corrections (set IPRINT to 1 for more
*---  verbose output)
         IF (LTHC2L) THEN
            Call FNSET("P_ORDPTHY",0) ! select order pert. theory: 1=LO, 2=NLO
            Call FNSET("P_THRESHCOR",2) ! select no. of loops in threshold correction
            Call FX9999CC(SCALEF, SCALEF, XSTHC, DXSUCNPC, DXSCORNPC)
         ENDIF

*---  Apply non-perturbative corrections to NLO cross section (set
*---  IPRINT to 1 for more verbose output)
Comment:          IF (LLO.AND.LNLO.AND.LMULT) THEN
Comment:             Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
Comment:             Call FNSET("P_THRESHCOR",0) ! deselect threshold corrections
Comment:             Call FNSET("P_NPCOR",1) ! select non-perturbative corrections
Comment:             Call FX9999CC(SCALEF, SCALEF, XSNPC, DXSUCNPC, DXSCORNPC)
Comment:          ENDIF

*---  Print out non-perturbative corrections (set IPRINT to 1 for more
*---  verbose output)
         IF (LMULT) THEN
            Call FNSET("P_ORDPTHY",0) ! select order pert. theory: 1=LO, 2=NLO
            Call FNSET("P_THRESHCOR",0) ! deselect threshold corrections
            Call FNSET("P_NPCOR",1) ! select non-perturbative corrections
            Call FX9999CC(SCALEF, SCALEF, XSNPC, DXSUCNPC, DXSCORNPC)
         ENDIF

*---  Cross section printout
         DO I=1,NObsBin
            IF ((ABS(xslo(I)).GT.1.D-99)) THEN
               kfac(I) = xsnlo(I) / xslo(I)
            ENDIF
            IF ((ABS(xsnlo(I)).GT.1.D-99)) THEN
               kthc(I) = xsthc(I) / xsnlo(I)
            ENDIF
         ENDDO
         WRITE(*,'(A)')DSEPL
         WRITE(*,'(A)')" Cross Sections"
         WRITE(*,"(A,G10.2)")" The scale factor chosen here is: ",
     >        SCALEF
         WRITE(*,'(A)')LSEPL
         CHTMP1 = DimLabel(1)
         CHTMP1 = "[ "//CHTMP1(1:12)//" ]"
         CHTMP2 = DimLabel(2)
         CHTMP2 = "[ "//CHTMP2(1:12)//" ]"
         WRITE(*,'(A)')"  IObs  Bin Size "//
     >        "IODim1  "//
     >        CHTMP1//"    "//
     >        "IODim2  "//
     >        CHTMP2//"  "//
     >        "LO cross section   "//
     >        "NLO cross section  "//
     >        "K factor           "//
     >        "K thr. corr."
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
     >           XSLO(I),XSNLO(I),KFAC(I),KTHC(I),KNPC(I)
         ENDDO
      ENDDO

*---  Print out data (set IPRINT to 1 for more verbose output)
      IF (LDATA) THEN
         Call FNSET("P_ORDPTHY",0) ! deselect pert. theory
         Call FNSET("P_THRESHCOR",0) ! deselect threshold corrections
         Call FNSET("P_NPCOR",0) ! deselect non-perturbative corrections
         Call FNSET("P_DATA",1) ! select data
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
     >        CHTMP2//"  "//
     >        "X Section          "//
     >        "Lower unc. uncert. "//
     >        "Upper unc. uncert. "//
     >        "Lower corr. uncert."//
     >        "Upper corr. uncert."
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
     >           DXSUCDATA(I,1),DXSUCDATA(I,2),
     >           DXSCORDATA(I,1),DXSCORDATA(I,2)
         ENDDO
      Endif
      
 999  FORMAT(1P,X,I5,X,G10.4,(X,I5,2(X,G10.4)),
     >     (X,I5,2(2X,E7.1)),5(X,E18.11),X)
      
      End
