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
      Character*33  CSEP33,DSEP33,LSEP33,SSEP33
      Character*66  CSEPS,DSEPS,LSEPS,SSEPS
      Character*132 CSEPL,DSEPL,LSEPL,SSEPL
      Character*255 FILENAME,PDFSET
      Integer i, j, IS, IPRINT, NDimBins(MxDim)
      Double Precision SCALEF(4)
      Data IPRINT/0/
      Data SCALEF/0.25D0,0.5D0,1.0D0,2.0D0/
      Data CSEP33,DSEP33,LSEP33,SSEP33/
     >     '#################################',
     >     '=================================',
     >     "---------------------------------",
     >     "*********************************"/

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
      Double Precision kfac(MxObsBin) 

*---  Initialization
      CSEPL = CSEP33//CSEP33//CSEP33//CSEP33
      DSEPL = DSEP33//DSEP33//DSEP33//DSEP33
      LSEPL = LSEP33//LSEP33//LSEP33//LSEP33
      SSEPL = SSEP33//SSEP33//SSEP33//SSEP33
      CSEPS = CSEP33//CSEP33
      DSEPS = DSEP33//DSEP33
      LSEPS = LSEP33//LSEP33
      SSEPS = SSEP33//SSEP33
Comment:       CSEPS(1:2) = "# "
      DSEPS(1:2) = "# "
      LSEPS(1:2) = "# "
      SSEPS(1:2) = "# "
Comment:       CSEPL(1:2) = "# "
Comment:       DSEPL(1:2) = "# "
Comment:       LSEPL(1:2) = "# "
Comment:       SSEPL(1:2) = "# "
      DO I=1,MxObsBin
         xslo(I)  = -1.d0
         xsnlo(I) = -1.d0
         kfac(I)  =  0.d0
      ENDDO

*---  Initial output
      WRITE(*,*)
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
         WRITE(*,*)"# ReadFNLOTable: Using PDF set: ",
     >        PDFSET(1:LEN_TRIM(PDFSET))
      ENDIF

*---  Too many arguments
      IF (IARGC().GT.2) THEN
         WRITE(*,*)
     >        "ReadFNLOTable: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF

*---  Read table once to get scenario information
      Call FNSET("P_REFTAB",0)  ! evaluate standard table:0 or reference:1
      Call FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
      Call FX9999CC(FILENAME, 1.D0, 1.D0, IPRINT, XSLO)

*---  Print out scenario information
      Call FX9999NF

*---  Initialize LHAPDF  
      call SetLHAPARM('SILENT')
      call InitPDFset(PDFSET(1:LEN_TRIM(PDFSET)))

*---  Initialize one member, 0=best fit member
      call InitPDF(0)

*---  Compute the cross sections
      WRITE(*,'(A)')CSEPL
      WRITE(*,'(A)')"ReadFNLOTable: Calculate cross sections"
      WRITE(*,'(A)')CSEPL

c     - new call: a single call for each scale
c     1st argument:  name of table
c     2nd argument:  xmur  prefactor for nominal ren-scale
c     any choice is possible, but please note 
c     that NNLO-NLL works only for xmur=xmuf
c     3rd argument:  xmuf  prefactor for nominal fact-scale
c     only a few choices are possible
c     (see output or table documentation)
c     4th argument:  0: no ascii output       1: print results
c     5th argument:  array to return results

*---  Loop over scale settings
      DO IS=1,4

*---  Evaluate table
         Call FNSET("P_REFTAB",0) ! evaluate standard table:0 or reference:1
         
*---  Calculate LO cross sections (set IPRINT to 1 for more verbose
*---  output) 
         Call FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
         Call FX9999CC(FILENAME, SCALEF(IS), SCALEF(IS), IPRINT, XSLO)
         
*---  Calculate NLO cross sections (set IPRINT to 1 for more verbose
*---  output) 
         Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
         Call FX9999CC(FILENAME, SCALEF(IS), SCALEF(IS), IPRINT, XSNLO)

*---  Calculate and print threshold corrections
c     Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
c     Call FNSET("P_THRESHCOR",2) ! select No. loops in threshold
c     corrections

*---  Cross section printout
         DO I=1,MxObsBin
            IF ((ABS(1.D0 + xslo(I)).GT.1.D-99) .OR.
     >           (ABS(1.D0 + xsnlo(I)).GT.1.D-99)) THEN
            ENDIF
            IF ((ABS(xslo(I)).GT.1.D-99)) THEN
               kfac(I) = xsnlo(I) / xslo(I)
            ENDIF
         ENDDO
         WRITE(*,'(A)')DSEPL
         WRITE(*,'(A)')" Cross Sections"
         WRITE(*,"(A,G10.2)")" The scale factor chosen here is: ",
     >        SCALEF(IS)
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
     >        "K factor"
         WRITE(*,'(A)')LSEPL
 900     FORMAT(1P,X,I5,X,G10.4,(X,I5,2(X,G10.4)),
     >        (X,I5,2(2X,E7.1)),3(X,E18.11),X)
         DO I=1,MxObsBin
            DO J=1,NDim
               IF (I.EQ.1) THEN
                  NDimBins(J) = 1 
               ELSEIF (LoBin(I-1,J).LT.LoBin(I,J)) THEN
                  NDimBins(J) = NDimBins(J) + 1 
               ELSEIF (LoBin(I,J).LT.LoBin(I-1,J)) THEN
                  NDimBins(J) = 1 
               ENDIF
            ENDDO
            IF ((ABS(1.D0 + xslo(I)).GT.1.D-99) .OR.
     >           (ABS(1.D0 + xsNlo(I)).GT.1.D-99)) THEN
               WRITE(*,900)I,BinSize(I),
     >              NDimBins(1),LoBin(I,1),UpBin(I,1),
     >              NDimBins(2),LoBin(I,2),UpBin(I,2),
     >              XSLO(I),XSNLO(I),KFAC(I)
            ENDIF
         ENDDO
      ENDDO

      End
