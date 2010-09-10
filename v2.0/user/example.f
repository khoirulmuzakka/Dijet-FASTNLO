      PROGRAM EXAMPLE
* -------------------------------------------------------------------
* M. Wobisch                                08/10/2010
*
* fastNLO - example program to compute Run I cross section
*           using PDFs from LHAPDF
*
* -------------------------------------------------------------------
      Implicit None
      Character*255 FILENAME
      Character*255 HISTOFILE
      Integer i          !, LENNOC

c - Attention - this is the most likely source of errors in fastNLO!!!
c        For each scenario, the result array must be declared 
c        consistent with its definition in the commonblock of 
c        the corresponding scenario. 
c           -> See the value of the parameter MxObsBin
c                           in the file [scenario].inc
c        We recommend to name the array according to the scenario
      Double Precision xst1001(200) 

c --- parse command line
      If ( IARGC().lt.1)  FILENAME = 'table.tab'
      If ( IARGC().lt.2)  HISTOFILE= 'fastnlo.hbk'
      If ( IARGC().gt.0)  Call GETARG(1,FILENAME)
      If ( IARGC().gt.1)  Call GETARG(2,HISTOFILE)
      If ( IARGC().gt.2)  Then
         Write(*,*) 'fastNLO: too many arguments given. Stopping'
         Return
      Endif

c - Initialize LHAPDF    - for CTEQ6.5M 
      call InitPDFset('/usr/local/share/lhapdf/PDFsets/cteq66.LHgrid')
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


c      Call FNSET("P_REFTAB",1) ! evaluate standard table:0 or reference:1
      Call FNSET("P_REFTAB",0) ! evaluate standard table:0 or reference:1

      Call FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
      Call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
c      Call FNSET("P_THRESHCOR",2) ! select No. loops in threshold corrections


      Call FX9999CC(FILENAME, 1.0d0, 1.0d0, 1, XST1001)
c      Call FX9999CC(FILENAME, 2.0d0, 2.0d0, 1, XST1001)
c      Call FX9999CC(FILENAME, 0.5d0, 0.5d0, 1, XST1001)
c      Call FX9999CC(FILENAME, 0.25d0, 0.25d0, 1, XST1001)
c      Call FX9999CC(FILENAME, 1.0d0, 0.5d0, 1, XST1001)
c      Call FX9999CC(FILENAME, 0.5d0, 1.0d0, 1, XST1001)

c - test aposteriori mur variation
c      Call FX9999CC(FILENAME, 1d0, 0.25d0, 1, XST1001)
c      Call FX9999CC(FILENAME, 2d0, 0.25d0, 1, XST1001)
c      Call FX9999CC(FILENAME, 1d0, 0.5d0, 1, XST1001)

c - the results of the last call can be accessed in the array:  xst1001(n)
c      n: continuous bin number for all observable bins


c      Call FX9999NF ! print scenario info


      End
