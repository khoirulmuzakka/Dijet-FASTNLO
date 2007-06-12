      PROGRAM EXAMPLE01
* -------------------------------------------------------------------
* M. Wobisch 06/29/2005
*
* fastNLO - example program to compute Run I cross section
*           using PDFs from LHAPDF
*
* -------------------------------------------------------------------
      implicit none
      CHARACTER*255 FILENAME
      double precision h1incl(4,4) , h1dijet(4,5) ! for FNH0001
      double precision xst1001(123,2,5) ! for FNT1001


c --- parse command line
      IF ( IARGC().EQ.0) THEN
c default filename
         FILENAME = 'table.txt.gz'
      ELSEIF ( IARGC().EQ.1) THEN
         CALL GETARG(1,FILENAME)
      ELSE
         write(*,*) 'fastNLO: more than one argument given. Stopping'
         RETURN 
      ENDIF
      
c - Initialize LHAPDF    - for CTEQ6.1M   
c      call InitPDFset('/disk2/work/wobisch/LHAPDFv4/PDFsets/cteq61.LHgrid')
c      call InitPDFset('/disk2/work/wobisch/LHAPDFv4/PDFsets/cteq4m.LHgrid')
c      call InitPDFset('/disk2/work/wobisch/LHAPDFv4/PDFsets/MRST2004nlo.LHgrid')
c      call InitPDFset('/disk2/work/wobisch/LHAPDFv4/PDFsets/H12000ms.LHgrid')
c      call InitPDFset('/disk2/work/wobisch/LHAPDFv4/PDFsets/ZEUS2002_TR.LHpdf')
c
c     call InitPDFset('/h1/h1gen/lhapdf/LHAPDFv4/PDFsets/cteq61.LHgrid')
      call InitPDFset('/h1/h1gen/lhapdf/LHAPDFv4/PDFsets/cteq5m1.LHgrid')
c - initialize one member, 0=best fit member
      call InitPDF(0)




c - compute the cross sections
      write(*,*) '     fastNLO: compute the cross section'
      call FNCALC(FILENAME, 1, XST1001)
c - now the results can be accessed in the array:  result(n,m,iord)
c      n: continuous bin number for all D0/CDF bins (see documentation)
c      m: contrib. of partonic subprocess (1-7) or total x-section (8)
c   iord: order  1 LO, 2 NLO correction (add both to get the NLO x-section)



c - if desired:  ASCII output of results
      call FNOUT


c - book, fill and store histograms
c      call FNHIST('path/fastnlo-mw-test.hbk')
c - temporarily disabled while changing array structure

C=========== and now do the H1 jets:   FNH0001
c      call tsh0001(h1incl,h1dijet)
c      call tsh0001out(h1incl,h1dijet)


      RETURN
      END
