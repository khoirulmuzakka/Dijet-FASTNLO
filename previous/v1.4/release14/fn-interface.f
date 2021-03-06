******************************************************************
* M. Wobisch - July 26, 2005           fn-interface.f
*
*   fastNLO user interface to PDF and alpha_s code
*   --> to be edited by user 
*       to interface own PDF/alphas code
*
* included:    
*      DOUBLE PRECISION FUNCTION FNALPHAS(MUR)
*      SUBROUTINE FNPDF(X,MUF,XPDF)
*
*  in the default version the PDF interface FNPDF is set up
*  to access LHAPDF - the alpha_s routine FNALPHAS calls
*  alphas-demo.f which is an iterative solution of the 2-loop RGE
*******************************************************************


      DOUBLE PRECISION FUNCTION FNALPHAS(MUR)
*-----------------------------------------------------------------
* MW 06/29/2005  - alphas interface to the fastNLO code
*
* alpha_s computation
*   input:   MUR     renormalization scale in GeV
*   output:  value of (alpha_s/2Pi) at scale MUR
*
*  !!!! again: this function must return  alpha_s/(2Pi)  !!!!!
*
*-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER nloop
      DOUBLE PRECISION MUR, ALPSMZ, PI
      DOUBLE PRECISION ALPS_IT 
c      DOUBLE PRECISION MY_FAVOURITE-ALPHAS
      PARAMETER (PI=3.1415927d0)

c === example: exact, iterative solution of the 2-loop RGE 
      nloop=2
      alpsmz=0.118              ! set here the value of alpha_s(Mz)
c      alpsmz=0.1185              ! for H1-2000 MSbar
c      alpsmz=0.1205             ! for MRST2004
      FNALPHAS = ALPS_IT(MUR,ALPSMZ,NLOOP)/2d0/PI


c - here you can call your own alpha_s code
c           -> remember to divide by 2Pi
c
c     FNALPHAS = MY_FAVOURITE-ALPHAS(MUR)
c

      RETURN
      END

C *****************************************************************

      SUBROUTINE FNPDF(X,MUF,XPDF)
*-----------------------------------------------------------------
* MW 06/29/2005 
*
* PDF interface to the fastNLO usercode
*
*   input   X       parton momentum fraction 
*           MUF     factorization scale in GeV
*
*   output  XPDF(-6:6) array of PDF momentum densities i.e. x*pdf!
*                      using the LHAPDF numbering convention:
*        tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t
*         -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 , 0, 1, 2, 3, 4, 5, 6
*
*-----------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION X, MUF, XPDF(-6:6)


c ======= example: interface to LHAPDF ============================
c remember to initialize the LHAPDF set first (in the main routine)
c                                       -> see the fastNLO example
c             call InitPDFset("cteq61.LHgrid")
c             call InitPDF(0)
      call evolvePDF(X,MUF,XPDF)
c


c === here you can call your own PDF code
c  -> remember that these are *momentum densities* i.e. x * PDF
c     and the scale is in GeV
c
c     call MY-FAVORITE-PDFS(....)
c
      RETURN
      END
