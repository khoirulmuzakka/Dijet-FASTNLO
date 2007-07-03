      PROGRAM EXAMPLE
* -------------------------------------------------------------------
* M. Wobisch 06/29/2005
*
* fastNLO - example program to compute Run I cross section
*           using PDFs from LHAPDF
*
* -------------------------------------------------------------------
      Implicit None
      CHARACTER*255 FILENAME, HISTOFILE
      Integer i

c - Attention!!! - this mus be declared consistent with its 
c                  definition in the commonblock!!!!!
      double precision xst1001(200) 

c --- parse command line
      IF ( IARGC().LT.1)  FILENAME = 'table.tab'
      IF ( IARGC().LT.2)  HISTOFILE= 'fastnlo.hbk'
      IF ( IARGC().GT.0)  CALL GETARG(1,FILENAME)
      IF ( IARGC().GT.1)  CALL GETARG(2,HISTOFILE)
      IF ( IARGC().GT.2)  THEN
         write(*,*) 'fastNLO: too many arguments given. Stopping'
         RETURN
      ENDIF

c - Initialize LHAPDF    - for CTEQ6.1M      MRST2004NLO:0.1205
      call InitPDFset('/usr/local/share/lhapdf/PDFsets/cteq65.LHgrid')
c - initialize one member, 0=best fit member
      call InitPDF(0)

c - compute the cross sections
      write(*,*) '     fastNLO: call user code'

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

      call FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
      call FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
c      call FNSET("P_THRESHCOR",2) ! select No. loops in threshold corrections
      call FX9999CC(FILENAME, 1.0d0 , 1.0d0, 1 , XST1001)
c      call FX9999NF             ! print scenario info

c - the results of the last call can be accessed in the array:  xst1001(n,iord)
c      n: continuous bin number for all D0/CDF bins (see documentation)
c   iord: order  1 LO, 
c                2 NLO correction (add 1,2 to get the NLO x-section)
c                3 NNLO-NLL correction (add 1,2,3 to get full prediction)


      END
