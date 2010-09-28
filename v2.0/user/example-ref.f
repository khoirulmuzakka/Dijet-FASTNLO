      PROGRAM EXAMPLEREF
* -------------------------------------------------------------------
* M. Wobisch                                08/10/2010
*
* fastNLO - example program to compare fastNLO v2 results with
*           the results stored in the fastNLO reference table
*
* -------------------------------------------------------------------
      Implicit None
      Character*255 FILENAME
      Character*255 HISTOFILE
      Integer i,j          !, LENNOC
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
      If ( IARGC().lt.1)  FILENAME = 'table.tab'
      If ( IARGC().lt.2)  HISTOFILE= 'fastnlo.hbk'
      If ( IARGC().gt.0)  Call GETARG(1,FILENAME)
      If ( IARGC().gt.1)  Call GETARG(2,HISTOFILE)
      If ( IARGC().gt.2)  Then
         Write(*,*) 'fastNLO: too many arguments given. Stopping'
         Return
      Endif

c - Initialize LHAPDF  
c - for Reference tables - use CTEQ6M
      call InitPDFset('/usr/local/share/lhapdf/PDFsets/cteq6mE.LHgrid')
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
