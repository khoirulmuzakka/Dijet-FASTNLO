      PROGRAM FASTNLOWEB
* -------------------------------------------------------------------
* M. Wobisch 06/29/2005
*
* fastNLO - program to compute jet cross sections
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
      CHARACTER*255 FILENAME, PDFSET, arg3, arg4,arg5,arg6
      CHARACTER*255 HISTOFILE
      INTEGER PDFMEMBER
      DOUBLE PRECISION alphasPDF,aspdf
      integer i

      Double Precision mur,muf, alphas, ASMZ
      COMMON /FNCEDAR/ ASMZ

c - Attention!!! - this mus be declared consistent with its 
c                  definition in the commonblock!!!!!
      double precision xst1001(900,3) 
 
      PDFMEMBER = 0
c --- parse command line
      IF ( IARGC().LT.1)  FILENAME = 'table.txt.gz'
      IF ( IARGC().LT.2)  HISTOFILE= 'fastnlo.hbk'
      IF ( IARGC().GT.0)  CALL GETARG(1,FILENAME)
      IF ( IARGC().GT.1)  CALL GETARG(2,PDFSET)
      IF ( IARGC().GT.2) THEN
         CALL GETARG(3,arg3)
         read(arg3,*) PDFMEMBER
      ENDIF
      IF ( IARGC().GT.3)  THEN 
         CALL GETARG(4,arg4)
         read(arg4,*) alphas
      ENDIF
      IF ( IARGC().GT.4)  THEN 
         CALL GETARG(5,arg5)
         read(arg5,*) muf
      ENDIF
      IF ( IARGC().GT.5) THEN
         CALL GETARG(6,arg6)
         read(arg6,*) mur
      ENDIF
      IF ( IARGC().GT.6)  THEN
         write(*,*) 'fastNLO: too many arguments given. Stopping'
         RETURN
      ENDIF


      if (abs(mur).lt.1E-07) mur = muf
      ASMZ = alphas

      write(*,*) "  "
      write(*,*) " in the matrixelements you use  alpha_s(Mz)=",alphas
      write(*,*) "  "

c      write(*,*) 'fastNLO input parameters',filename,pdfset,pdfmember,alphas,muf,mur

c - Initialize LHAPDF
      call InitPDFset(PDFSET)

c - initialize one member, 0=best fit member
      call InitPDF(PDFMEMBER)


      aspdf = alphasPDF(91.1876d0)
      write(*,*) " in the pdf was used: alpha_s(Mz)=",aspdf

c - compute the cross sections
      write(*,*) '     fastNLO: compute the cross section'

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
