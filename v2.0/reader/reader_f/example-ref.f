      PROGRAM EXAMPLEREF
* -------------------------------------------------------------------
* M. Wobisch                                08/10/2010
*
* fastNLO - example program to compare fastNLO v2 results with
*           the results stored in the fastNLO reference table
*
* KR 07112010: Improved commandline steering
* -------------------------------------------------------------------
      Implicit None
      Character*255 FILENAME,HISTOFILE,PDFSET
      Integer i,j
      Double Precision scale(4)
      
c - Attention - this is the most likely source of errors in fastNLO!!!
c        For each scenario, the result array must be declared 
c        consistent with its definition in the commonblock of 
c        the corresponding scenario. 
c           -> See the value of the parameter MxObsBin
c                           in the file [scenario].inc
c        We recommend to name the array according to the scenario
      Double Precision xs(200),ref(200),xsnlo(200),refnlo(200)

      Scale(1) = 1d0
      Scale(2) = 2d0
      Scale(3) = 0.5d0
      Scale(4) = 0.25d0

c --- parse command line
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
            WRITE(*,*)'  HBOOK output file, def. = fastnlo.hbk'
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
     >           "EXAMPLE: WARNING! No table name given, "//
     >           "taking the default table.tab instead!"
            WRITE(*,*)"      For an explanation of command line "//
     >           "arguments type:"
            WRITE(*,*)"      ./example -h"
         ELSE
            WRITE(*,*)"ALLUNC: Evaluating table: ",
     >           FILENAME(1:LEN_TRIM(FILENAME))
         ENDIF
      ENDIF

*---HBOOK filename
      HISTOFILE = "X"
      IF (IARGC().GE.2) THEN
         CALL GETARG(2,HISTOFILE)
      ENDIF
      IF (IARGC().LT.2.OR.HISTOFILE(1:1).EQ."_") THEN
         HISTOFILE = "fastnlo.hbk"
         WRITE(*,*)
     >        "EXAMPLE: WARNING! No output filename given, "//
     >        "taking fastnlo.hbk instead!"
      ELSE
         WRITE(*,*)"EXAMPLE: Creating output file: ",
     >        HISTOFILE(1:LEN_TRIM(HISTOFILE))
      ENDIF
      
*---PDF set
      PDFSET = "X"
      IF (IARGC().GE.3) THEN
         CALL GETARG(3,PDFSET)
      ENDIF
      IF (IARGC().LT.3.OR.PDFSET(1:1).EQ."_") THEN
         PDFSET = "cteq6mE.LHgrid"
         WRITE(*,*)
     >        "EXAMPLE: WARNING! No PDF set given, "//
     >        "taking cteq66.LHgrid instead!"
      ELSE
         WRITE(*,*)"EXAMPLE: Using PDF set: ",
     >        PDFSET(1:LEN_TRIM(PDFSET))
      ENDIF

*---Too many arguments
      IF (IARGC().GT.3) THEN
         WRITE(*,*)"\nEXAMPLE: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF

ckr      If ( IARGC().lt.1)  FILENAME = 'table.tab'
ckr      If ( IARGC().lt.2)  HISTOFILE= 'fastnlo.hbk'
ckr      If ( IARGC().gt.0)  Call GETARG(1,FILENAME)
ckr      If ( IARGC().gt.1)  Call GETARG(2,HISTOFILE)
ckr      If ( IARGC().gt.2)  Then
ckr         Write(*,*) 'fastNLO: too many arguments given. Stopping'
ckr         Return
ckr      Endif

c - Initialize LHAPDF  
c - for Reference tables - use CTEQ6M
      call InitPDFset(PDFSET(1:LEN_TRIM(PDFSET)))
c - initialize one member, 0=best fit member
      call InitPDF(0)

c - compute the cross sections
      Write(*,*) '     fastNLO: call user code'

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



      Write(*,*) ' evaluate relative differences to reference results ',
     +     'for table ',FILENAME

      Do i=1,4                  ! scale
         Call FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
         Call FNSET("P_REFTAB",1) ! evaluate standard table:0 or reference:1
         Call FX9999CC(FILENAME, scale(i), scale(i), 1, ref)
         Call FNSET("P_REFTAB",0) ! evaluate standard table:0 or reference:1
         Call FX9999CC(FILENAME, scale(i), scale(i), 1, xs)

         Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
         Call FNSET("P_REFTAB",1) ! evaluate standard table:0 or reference:1
         Call FX9999CC(FILENAME, scale(i), scale(i), 1, refnlo)
         Call FNSET("P_REFTAB",0) ! evaluate standard table:0 or reference:1
         Call FX9999CC(FILENAME, scale(i), scale(i), 1, xsnlo)

         Write(*,*)
         Write(*,*) '    now for scale factor ',i,' =',scale(i)
         Do j=1,110
            Write(*,*) j, 100d0*(xs(j)-ref(j))/ref(j), 
     +           100d0*(xsnlo(j)-refnlo(j))/refnlo(j)
         Enddo
      Enddo

      End
