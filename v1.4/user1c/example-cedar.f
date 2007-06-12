      PROGRAM EXAMPLE
* -------------------------------------------------------------------
* M. Wobisch 06/29/2005
*
* fastNLO - example program to compute Run I cross section
*           using PDFs from LHAPDF
*
* input:  Table
*         PDFset (LHAPDF string)
*         alphas(Mz)
*         mu_f
*         mu_r  (if mur==0 then mu_r=mu_f)
*
* -------------------------------------------------------------------
      implicit none
      CHARACTER*255 FILENAME, PDFSET, arg3, arg4,arg5
      CHARACTER*255 HISTOFILE
      integer i
      Double Precision mur,muf, alphas, ASMZ
      COMMON /FNCEDAR/ ASMZ


c - Attention!!! - this mus be declared consistent with its 
c                  definition in the commonblock!!!!!
      double precision xst1001(900,3) 
 

c --- parse command line
      IF ( IARGC().LT.1)  FILENAME = 'table.txt.gz'
      IF ( IARGC().LT.2)  HISTOFILE= 'fastnlo.hbk'
      IF ( IARGC().GT.0)  CALL GETARG(1,FILENAME)
      IF ( IARGC().GT.1)  CALL GETARG(2,PDFSET)
      IF ( IARGC().GT.2)  CALL GETARG(3,arg3)
      read(arg3,*) alphas
      IF ( IARGC().GT.3)  CALL GETARG(4,arg4)
      read(arg4,*) muf
      IF ( IARGC().GT.4)  CALL GETARG(5,arg5)
      read(arg5,*) mur
      IF ( IARGC().GT.5)  THEN
         write(*,*) 'fastNLO: too many arguments given. Stopping'
         RETURN
      ENDIF


      if (abs(mur).lt.1E-07) mur = muf
      write(*,*) "  "
      write(*,*) " in the matrixelements you use  alpha_s(Mz)=",alphas
      write(*,*) "  "

c      write(*,*) 'fastNLO input parameters',filename,pdfset,alphas,muf,mur

c - Initialize LHAPDF    - for CTEQ6.1M   
c      call InitPDFset('/disk2/work/wobisch/lhapdf-4.1/PDFsets/cteq61.LHgrid')
c
c      call InitPDFset('/h1/h1gen/lhapdf/LHAPDFv4/PDFsets/cteq61.LHgrid')
c - initialize one member, 0=best fit member
c      call InitPDF(0)

c - .... or initialize CTEQ code - needs table file:  ctq61.00.tbl
      Call SetCtq6(200)         ! for CTEQ6.1M


c - compute the cross sections
c      write(*,*) '     fastNLO: compute the cross section'

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

      call FX9999CC(FILENAME, mur , muf, 1 , XST1001)


c - the results of the last call can be accessed in the array:  xst1001(n,iord)
c      n: continuous bin number for all D0/CDF bins (see documentation)
c   iord: order  1 LO, 
c                2 NLO correction (add 1,2 to get the NLO x-section)
c                3 NNLO-NLL correction (add 1,2,3 to get full prediction)


      write(*,*) "  "
      write(*,*) " in the matrixelements you use  alpha_s(Mz)=",alphas
      write(*,*) "  "


      END
