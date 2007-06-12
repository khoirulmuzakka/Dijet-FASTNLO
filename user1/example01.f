      PROGRAM EXAMPLE
* -------------------------------------------------------------------
* M. Wobisch 06/29/2005
*
* fastNLO - example program to compute Run I cross section
*           using PDFs from LHAPDF
*
* -------------------------------------------------------------------
      implicit none
      CHARACTER*255 FILENAME
      CHARACTER*255 HISTOFILE

      double precision h1incl(4,4) , h1dijet(4,5) ! for FNH0001

c - Attention!!! - this mus be declared consistent with its 
c                  definition in the commonblock!!!!!
      double precision xst1001(286,3,4) 


c --- parse command line
      IF ( IARGC().LT.1)  FILENAME = 'table.txt.gz'
      IF ( IARGC().LT.2)  HISTOFILE= 'fastnlo.hbk'
      IF ( IARGC().GT.0)  CALL GETARG(1,FILENAME)
      IF ( IARGC().GT.1)  CALL GETARG(2,HISTOFILE)
      IF ( IARGC().GT.2)  THEN
         write(*,*) 'fastNLO: too many arguments given. Stopping'
         RETURN
      ENDIF

      
c - Initialize LHAPDF    - for CTEQ6.1M   
      call InitPDFset('/disk2/work/wobisch/lhapdf-4.1/PDFsets/cteq61.LHgrid')
c
c      call InitPDFset('/h1/h1gen/lhapdf/LHAPDFv4/PDFsets/cteq61.LHgrid')
c - initialize one member, 0=best fit member
      call InitPDF(0)

c - .... or initialize CTEQ code - needs table file:  ctq61.00.tbl
c      Call SetCtq6(200)         ! for CTEQ6.1M


c - compute the cross sections
      write(*,*) '     fastNLO: compute the cross section'

c- new call
c         1st argument:  name of table
c         2nd argument:  xmur  negative: mur=-xmur*muf
c                              positive: mur= xmur*mur_0
c                              -> NNLO-NLL works only for xmur=-1.0
c         3rd argument:  0: only central fact scale    1: also variations
c         4th argument:  0: no ascii output            1: print results
c         5th argument:  array to return results
      call FX9999CC(FILENAME, -1.0d0 , 1, 1 , XST1001)


c - now the results can be accessed in the array:  xst1001(n,iord,iscale)
c      n: continuous bin number for all D0/CDF bins (see documentation)
c   iord: order  1 LO, 
c                2 NLO correction (add 1,2 to get the NLO x-section)
c                3 NNLO-NLL correction (add 1,2,3 to get full prediction)
c iscale: scale setting for mu_r,mu_f (see output for values)



c - book, fill and store histograms - works now! (Sept 15, 2005 MW)
c - >> works now also for threshold corrections
c      call FNHIST(HISTOFILE,1)  ! to plot x-sect histos
c      call FNHIST(HISTOFILE,2)  ! to plot x-sect histos + x-histos



C=========== and now do the H1 jets:   FNH0001  >>> very old MW code
c      call tsh0001(h1incl,h1dijet)
c      call tsh0001out(h1incl,h1dijet)


      END
