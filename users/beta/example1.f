      PROGRAM EXAMPLE
* -------------------------------------------------------------------
* M. Wobisch 06/29/2005
*
* fastNLO - example program for the fastNLO beta-release
*           -> computation of Run I cross jet sections
*              using PDFs from LHAPDF
*
* -------------------------------------------------------------------
      implicit none
      CHARACTER*255 FILENAME
      double precision xst0001(123,2,4) ! for FNB0001
      integer iset

      
c - Initialize LHAPDF    - for CTEQ6.1M   
      call InitPDFset('/disk2/work/wobisch/lhapdf-4.1/PDFsets/cteq61.LHgrid')
c      call InitPDFset('/h1/h1gen/lhapdf/LHAPDFv4/PDFsets/cteq61.LHgrid')

c - initialize one member, 0=best fit member
      call InitPDF(0)




c - please select the table:
c     the gzipped file is much smaller to store ...
c         ... but the ASCII file is much faster to read!
c      FILENAME = 'fbt1001.tab.gz'
      FILENAME = 'fbt1001.tab'




      write(*,*) '     fastNLO: compute the cross section'

c - compute the cross sections for the central CTEQ6.1M PDFs
c           for different ren./fact. scales

c         1st argument:  name of table
c         2nd argument:  0: only central scale       1: also scale variations
c         3rd argument:  0: no ascii output          1: print results
c         4th argument:  array to return results
      call FB0001CC(FILENAME, 1, 1, XST0001)
c
c - LHAPDF needs to be initialized before FB0001CC is called!
c


c - now the results can be accessed in the array:  xst1001(n,iord,iscale)
c      n: continuous bin number for all D0/CDF bins (see documentation)
c   iord: order  1 LO, 2 NLO correction (add both to get the NLO x-section)
c iscale: scale setting for mu_r,mu_f (see output for values)




c - now compute the cross sections for the 40 CTEQ6.1M error PDFs
c             (but only for central scale  mu=pT/2)

      write(*,*) '   ------------------------------------------'
      write(*,*) '   ----  and now the PDF uncertainties...'
      write(*,*) '   ------------------------------------------'
      do iset =1,40
         write(*,*) '   compute cross sections for CTEQ6.1 set No. ',iset
         call InitPDF(iset)
         call FB0001CC(FILENAME, 0, 1, XST0001)
      enddo

      END
