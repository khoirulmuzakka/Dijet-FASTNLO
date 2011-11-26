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
      DOUBLE PRECISION MUR
      
      DOUBLE PRECISION ALPSMZ, PI4, PI
      DOUBLE PRECISION ALPS_IT 
      INTEGER IFIRST, NLOOP

      DATA IFIRST/0/      

      IF (IFIRST.EQ.0) THEN
         IFIRST = 1
         PI4 = 4D0 * 4D0 * ATAN(1D0)
         PI  = PI4/4D0
      ENDIF

*     Example:
*     Calculation of alpha_s in the MSbar scheme for given alpha_s(Mz)
*     using exact, iterative solution of 2-/3-/4-loop formulas
*     as used by GRV hep-ph/9806404
      NLOOP  = 2
      ALPSMZ = 0.1185           ! Bethke 2011 (Ref. ???)

*     One can also call ones own alpha_s code here
*     --> Only one has to remember to divide by 2Pi!
      FNALPHAS = ALPS_IT(MUR,ALPSMZ,NLOOP)/2D0/PI
ckr Fix alpha_s for debugging
ckr      FNALPHAS = ALPSMZ/2D0/PI

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
