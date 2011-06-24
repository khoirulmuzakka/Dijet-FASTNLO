      PROGRAM EXAMPLE
* -------------------------------------------------------------------
* M. Wobisch                                08/10/2010
*
* fastNLO - example program to derive QCD cross sections
*           from fastNLO v2 tables using PDFs from LHAPDF
*
* KR 07112010: Improved commandline steering
* KR 16052011: Simplify and improve pure text output
* -------------------------------------------------------------------
      Implicit None
      Character*255 FILENAME,PDFSET
      Integer i, j, IPRINT
      Double Precision SCALEF(4)
      Data IPRINT/0/
      Data SCALEF/0.25D0,0.5D0,1.0D0,2.0D0/

c - Attention - this is the most likely source of Fortran errors in fastNLO!!!
c        For each scenario, the result array must be declared at least  
c        as large as in the definition in the common block of the
c        corresponding scenario. 
c           -> See the value of the parameter MxObsBin
c                           in the file [scenario].inc
c        We recommend to name the array according to the scenario
c        Adapt the following to your scenario!
      Integer MxObsBin
      Parameter (MxObsBin = 200)
      Double Precision xslo(MxObsBin) 
      Double Precision xsnlo(MxObsBin) 
      Double Precision kfac(MxObsBin) 

*---Initialization
      DO I=1,MxObsBin
         xslo(I)  = -1.d0
         xsnlo(I) = -1.d0
         kfac(I)  =  0.d0
      ENDDO


*---Parse command line
      WRITE(*,*)" "
      WRITE(*,*)"########################################"//
     >     "################################"
      WRITE(*,*)"# EXAMPLE"
      WRITE(*,*)"########################################"//
     >     "################################"
      WRITE(*,*)"# Example program to derive QCD cross sections"
      WRITE(*,*)"# from fastNLO v2 tables using PDFs from LHAPDF"
      WRITE(*,*)"########################################"//
     >     "################################"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      WRITE(*,*)"EXAMPLE: Program Steering"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      IF (IARGC().LT.1) THEN
         FILENAME = "table.tab"
         WRITE(*,*)
     >        "EXAMPLE: WARNING! No table name given, "//
     >        "taking the default table.tab instead!"
         WRITE(*,*)"      For an explanation of command line "//
     >        "arguments type:"
         WRITE(*,*)"      ./example -h"
      ELSE
         CALL GETARG(1,FILENAME)
         IF (FILENAME(1:LEN_TRIM(FILENAME)).EQ."-h") THEN
            WRITE(*,*)' '
            WRITE(*,*)'Usage: ./example [arguments]'
            WRITE(*,*)'  Table input file, def. = table.tab'
            WRITE(*,*)'  PDF set, def. = cteq6mE.LHgrid'
            WRITE(*,*)' '
            WRITE(*,*)'  Give full path(s) if these are not in the cwd.'
            WRITE(*,*)'  Use \"_\" to skip changing a default argument.'
            WRITE(*,*)' '
            STOP
         ELSEIF (FILENAME(1:1).EQ."_") THEN
            FILENAME = "table.tab"
            WRITE(*,*)
     >           "EXAMPLE: WARNING! No table name given, "//
     >           "taking the default table.tab instead!"
         ELSE
            WRITE(*,*)"EXAMPLE: Evaluating table: ",
     >           FILENAME(1:LEN_TRIM(FILENAME))
         ENDIF
      ENDIF

*---PDF set
      PDFSET = "X"
      IF (IARGC().GE.2) THEN
         CALL GETARG(2,PDFSET)
      ENDIF
      IF (IARGC().LT.2.OR.PDFSET(1:1).EQ."_") THEN
         PDFSET = "cteq6mE.LHgrid"
         WRITE(*,*)
     >        "EXAMPLE: WARNING! No PDF set given, "//
     >        "taking cteq6mE.LHgrid instead!"
      ELSE
         WRITE(*,*)"EXAMPLE: Using PDF set: ",
     >        PDFSET(1:LEN_TRIM(PDFSET))
      ENDIF

*---Too many arguments
      IF (IARGC().GT.2) THEN
         WRITE(*,*)"EXAMPLE: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF

*---Read table once to get scenario information
      Call FNSET("P_REFTAB",0)  ! evaluate standard table:0 or reference:1
      Call FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
      Call FX9999CC(FILENAME, 1.D0, 1.D0, IPRINT, XSLO)
*---Print out scenario information
      Call FX9999NF

*---Initialize LHAPDF  
      call InitPDFset(PDFSET(1:LEN_TRIM(PDFSET)))

*---Initialize one member, 0=best fit member
      call InitPDF(0)

*---Compute the cross sections
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      WRITE(*,*)"EXAMPLE: Calculate cross sections"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"

c- new call: a single call for each scale
c         1st argument:  name of table
c         2nd argument:  xmur  prefactor for nominal ren-scale
c                              any choice is possible, but please note 
c                              that NNLO-NLL works only for xmur=xmuf
c         3rd argument:  xmuf  prefactor for nominal fact-scale
c                              only a few choices are possible
c                              (see output or table documentation)
c         4th argument:  0: no ascii output       1: print results
c         5th argument:  array to return results

*---Loop over scale settings
      DO J=1,4

*---Evaluate table
         Call FNSET("P_REFTAB",0) ! evaluate standard table:0 or reference:1
         
*---Calculate LO cross sections (set IPRINT to 1 for more verbose output) 
         Call FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
         Call FX9999CC(FILENAME, SCALEF(J), SCALEF(J), IPRINT, XSLO)
      
*---Calculate NLO cross sections (set IPRINT to 1 for more verbose output) 
         Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
         Call FX9999CC(FILENAME, SCALEF(J), SCALEF(J), IPRINT, XSNLO)

*---Calculate and print threshold corrections
c      Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
c      Call FNSET("P_THRESHCOR",2) ! select No. loops in threshold corrections

*---Cross section printout
         DO I=1,MxObsBin
            IF ((ABS(1.D0 + xslo(I)).GT.1.D-99) .OR.
     >           (ABS(1.D0 + xsnlo(I)).GT.1.D-99)) THEN
            ENDIF
            IF ((ABS(xslo(I)).GT.1.D-99)) THEN
               kfac(I) = xsnlo(I) / xslo(I)
            ENDIF
         ENDDO
         WRITE(*,*)"========================================"//
     >        "================================"
         WRITE(*,*)" Cross Sections"
         WRITE(*,*)" The scale factor no. ",J," is: ",SCALEF(J)
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         WRITE(*,*)" bin     LO cross section      "//
     >        "NLO cross section     K factor"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
 900     FORMAT(1P,I5,3(4X,E18.11))
         DO I=1,MxObsBin
            IF ((ABS(1.D0 + xslo(I)).GT.1.D-99) .OR.
     >           (ABS(1.D0 + xsNlo(I)).GT.1.D-99)) THEN
               WRITE(*,900)I,XSLO(I),XSNLO(I),KFAC(I)
            ENDIF
         ENDDO
      ENDDO

      End
