      PROGRAM EXAMPLEREF
* -------------------------------------------------------------------
* M. Wobisch                                08/10/2010
*
* fastNLO - example program to compare fastNLO v2 results with
*           the results stored in the fastNLO reference table
*
* KR 07112010: Improved commandline steering
* KR 16052011: Simplify and improve pure text output
* -------------------------------------------------------------------
      Implicit None
      Character*255 FILENAME,PDFSET
      Integer i, j, IPRINT
      Double Precision DREF
      Data IPRINT/0/
      
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
      Double Precision xslo(MxObsBin),reflo(MxObsBin)
      Double Precision xsnlo(MxObsBin),refnlo(MxObsBin)

*---Initialization
      DO I=1,MxObsBin
         xslo(I)   = -1.d0
         xsnlo(I)  = -1.d0
         reflo(I)  = -1.d0
         refnlo(I) = -1.d0
      ENDDO



*---Parse command line
      WRITE(*,*)" "
      WRITE(*,*)"########################################"//
     >     "################################"
      WRITE(*,*)"# EXAMPLE-REF"
      WRITE(*,*)"########################################"//
     >     "################################"
      WRITE(*,*)"# Example program to estimate the accuracy"
      WRITE(*,*)"# of the fastNLO v2 approximation table"
      WRITE(*,*)"########################################"//
     >     "################################"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      WRITE(*,*)"EXAMPLE-REF: Program Steering"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      IF (IARGC().LT.1) THEN
         FILENAME = "table.tab"
         WRITE(*,*)
     >        "EXAMPLE-REF: WARNING! No table name given, "//
     >        "taking the default table.tab instead!"
         WRITE(*,*)"      For an explanation of command line "//
     >        "arguments type:"
         WRITE(*,*)"      ./example-ref -h"
      ELSE
         CALL GETARG(1,FILENAME)
         IF (FILENAME(1:LEN_TRIM(FILENAME)).EQ."-h") THEN
            WRITE(*,*)' '
            WRITE(*,*)'Usage: ./example-ref [arguments]'
            WRITE(*,*)'  Table input file, def. = table.tab'
            WRITE(*,*)'  PDF set, def. = cteq6mE.LHgrid'
            WRITE(*,*)'  Warning! The ref. PDF should not be changed!'
            WRITE(*,*)' '
            WRITE(*,*)'  Give full path(s) if these are not in the cwd.'
            WRITE(*,*)'  Use \"_\" to skip changing a default argument.'
            WRITE(*,*)' '
            STOP
         ELSEIF (FILENAME(1:1).EQ."_") THEN
            FILENAME = "table.tab"
            WRITE(*,*)
     >           "EXAMPLE-REF: WARNING! No table name given, "//
     >           "taking the default table.tab instead!"
         ELSE
            WRITE(*,*)"EXAMPLE-REF: Evaluating reference table: ",
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
     >        "EXAMPLE-REF: No PDF set given, "//
     >        "taking the default of cteq6mE.LHgrid!"
      ELSE
         WRITE(*,*)"EXAMPLE-REF: Using PDF set: ",
     >        PDFSET(1:LEN_TRIM(PDFSET))
      ENDIF

*---Too many arguments
      IF (IARGC().GT.2) THEN
         WRITE(*,*)"EXAMPLE-REF: ERROR! Too many arguments, aborting!"
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
      WRITE(*,*)"EXAMPLE-REF: Calculate cross sections"
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

*---Evaluate reference table
      Call FNSET("P_REFTAB",1)  ! evaluate standard table:0 or reference:1

*---Calculate reference LO cross sections (set IPRINT to 1 for more verbose output) 
      Call FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
      Call FX9999CC(FILENAME, 1.d0, 1.d0, IPRINT, reflo)

*---Evaluate normal table, calculate approximated LO cross sections
      Call FNSET("P_REFTAB",0)  ! evaluate standard table:0 or reference:1
      Call FX9999CC(FILENAME, 1.d0, 1.d0, IPRINT, xslo)

*---Evaluate reference table
      Call FNSET("P_REFTAB",1)  ! evaluate standard table:0 or reference:1

*---Calculate reference NLO cross sections (set IPRINT to 1 for more verbose output) 
      Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
      Call FX9999CC(FILENAME, 1.d0, 1.d0, IPRINT, refnlo)

*---Evaluate normal table, calculate approximated NLO cross sections
      Call FNSET("P_REFTAB",0)  ! evaluate standard table:0 or reference:1
      Call FX9999CC(FILENAME, 1.d0, 1.d0, IPRINT, xsnlo)

*---Cross section printout
      WRITE(*,*)"========================================"//
     >     "================================"
      WRITE(*,*)" Relative Algorithmic Errors"
      WRITE(*,*)" The scale factor is ",1.d0
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      WRITE(*,*)" bin     LO cross section      "//
     >     "LO reference          algorithmic deviation"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
 900  FORMAT(1P,I5,3(4X,E18.11))
      Do I=1,MxObsBin
         IF (ABS(REFLO(I)).GT.1.D-99) THEN
            DREF = (XSLO(I)-REFLO(I)) / REFLO(I)
            IF ((ABS(1.D0 + xslo(I)).GT.1.D-99)) THEN
               WRITE(*,900)I,XSLO(I),REFLO(I),DREF
            ENDIF
         ENDIF
      Enddo
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      WRITE(*,*)" bin     NLO cross section     "//
     >     "NLO reference         algorithmic deviation"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      Do i=1,MxObsBin
         IF (ABS(REFNLO(I)).GT.1.D-99) THEN
            DREF = (XSNLO(I)-REFNLO(I)) / REFNLO(I)
            IF ((ABS(1.D0 + xsNlo(I)).GT.1.D-99)) THEN
               WRITE(*,900)I,XSNLO(I),REFNLO(I),DREF
            ENDIF
         ENDIF
      Enddo
      
      End
