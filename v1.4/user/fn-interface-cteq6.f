******************************************************************
* M. Wobisch - July 26, 2005           fn-interface.f
*
*   fastNLO user interface to PDF and alpha_s code
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
*   input   MUR     renormalization scale in GeV
*   output  value of alpha_s at scale MUR divided by 2Pi
*
*  !!!! again: this function must return  alpha_s/(2Pi)  !!!!!
*
*-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER nloop
      DOUBLE PRECISION MUR, ALPSMZ, PI
c      DOUBLE PRECISION MY_FAVOURITE-ALPHAS
      DOUBLE PRECISION ALPS_IT
      PARAMETER (PI=3.1415927d0)

c - here you can call your own alpha_s code
c           -> remember to divide by 2Pi
c
c     FNALPHAS = MY_FAVOURITE-ALPHAS(MUR)
c

c === example: an exact, iterative solution of the 2-loop RGE
      nloop=2
      alpsmz=0.118              ! set here the value of alpha_s(Mz)
      FNALPHAS = ALPS_IT(MUR,ALPSMZ,NLOOP)/2d0/PI

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
*   output  XPDF(-6:6) array of PDF momentum densities i.e. x*pdf!
*                      using the LHAPDF numbering convention:
*        tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t
*         -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 , 0, 1, 2, 3, 4, 5, 6
*
*-----------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION X, MUF, XPDF(-6:6)
      integer i

c === here you can call your own PDF code
c  -> remember that these are *momentum densities* i.e. x*PDF
c    and the scale is in GeV
c
c     call MY-FAVORITE-PDFS(....)
c


c - use the CTEQ6 original code
      DOUBLE PRECISION Ctq6Pdf
      XPDF(-6) =0d0
      XPDF(-5) = X * Ctq6Pdf (-5, X, muf)
      XPDF(-4) = X * Ctq6Pdf (-4, X, muf)
      XPDF(-3) = X * Ctq6Pdf (-3, X, muf)
      XPDF(-2) = X * Ctq6Pdf (-1, X, muf)  ! swapped!!!
      XPDF(-1) = X * Ctq6Pdf (-2, X, muf)  ! swapped
      XPDF(0) = X * Ctq6Pdf (0, X, muf)
      XPDF(1) = X * Ctq6Pdf (2, X, muf)     ! swapped
      XPDF(2) = X * Ctq6Pdf (1, X, muf)     ! swapped
      XPDF(3) = X * Ctq6Pdf (3, X, muf)
      XPDF(4) = X * Ctq6Pdf (4, X, muf)
      XPDF(5) = X * Ctq6Pdf (5, X, muf)
      XPDF(6) =0d0

c ======= as example this is the interface to LHAPDF ==================
c
c remember to initialize the LHAPDF set first (in the main routine)
c                                       -> see the fastNLO example
c             call InitPDFset("cteq61.LHgrid")
c             call InitPDF(0)
c
c - get ALL pdfs in one call
c - note: the momentum densities are returned, i.e. x X pdf!
c - the LHAPDF numbering convention:
c - tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t
c -  -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 , 0, 1, 2, 3, 4, 5, 6
c
c      call evolvePDF(X,MUF,XPDF)
c
      RETURN
      END
