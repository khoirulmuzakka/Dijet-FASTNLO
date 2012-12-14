***********************************************************************
*
*     fastNLO user interface to PDF and alpha_s code
*
*     Initial version: M. Wobisch, 2005
*     Updated for v2:  K. Rabbertz, M. Wobisch, 2011
*
*     Contains:
*     ---------
*     DOUBLE PRECISION FUNCTION FNALPHAS(MUR)
*     SUBROUTINE FNPDF(X,MUF,XPDF)
*
*     By default:
*     -----------
*     - The PDF interface FNPDF is set up to access LHAPDF
*     - The alpha_s routine FNALPHAS calls a calculation of
*     > alpha_s in the MSbar scheme for given alpha_s(Mz)
*     > using an exact, iterative solution of 2-/3-/4-loop
*     > formulas as used by GRV hep-ph/9806404
*
***********************************************************************

      DOUBLE PRECISION FUNCTION FNALPHAS(MUR)
***********************************************************************
*
*     Interface for alpha_s computation
*
*     Input:
*     ------
*     MUR       renormalization scale in GeV
*
*     Output:
*     -------
*     FNALPHAS  value of (alpha_s/2Pi) at scale MUR
*
*     ATTENTION: This function MUST return alpha_s/(2Pi) !
*
***********************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION MUR

      DOUBLE PRECISION ALPSMZ, PI
      DOUBLE PRECISION ALPS_IT,ALPS_IT_FNLO14
      INTEGER IFIRST, NLOOP

      DATA IFIRST/0/

      IF (IFIRST.EQ.0) THEN
         IFIRST = 1
         PI = 4D0 * ATAN(1D0)
      ENDIF

*---  Example:
*     Calculation of alpha_s in the MSbar scheme for given alpha_s(Mz)
*     using an exact, iterative solution of 2-/3-/4-loop formulas
*     as used by GRV hep-ph/9806404
      NLOOP  = 2
      ALPSMZ = 0.1184D0         ! PDG 2012; previous was 0.1185D0, Bethke 2011.

*---  One can also call ones own alpha_s code here ...
*---  --> Only one has to remember to divide by 2Pi!
C---  New code FNLO v2
      FNALPHAS = ALPS_IT(MUR,ALPSMZ,NLOOP)/2D0/PI
C---  Old code FNLO v14
C---  FNALPHAS = ALPS_IT_FNLO14(MUR,ALPSMZ,NLOOP)/2D0/PI
C---  Fix alpha_s for debugging
C---  FNALPHAS = ALPSMZ/2D0/PI

      RETURN
      END
***********************************************************************

      SUBROUTINE FNPDF(X,MUF,XPDF)
***********************************************************************
*
*     PDF interface to the fastNLO usercode
*
*     Input:
*     ------
*     X    parton momentum fraction
*     MUF  factorization scale in GeV
*
*     Output:
*     -------
*     XPDF(-6:6)  array of PDF momentum densities i.e. x*pdf !
*     >           using the LHAPDF numbering convention:
*     >           tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t
*     >            -6 , -5 , -4 , -3 , -2 , -1 ,0,1,2,3,4,5,6
*
***********************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION X, MUF, XPDF(-6:6)
*---  Example:
*     Interface to LHAPDF
*     Remember to initialize the LHAPDF set first (in the main routine):
*     >  CALL INITPDFSET("cteq61.LHgrid")
*     >  CALL INITPDF(0)
      CALL EVOLVEPDF(X,MUF,XPDF)
ckr Temporary
*---  Temporary fix for lhapdf-5.8.7 and lhapdf-5.8.8 bug with ABKM09 or ABM11
*---  Set top-antitop PDF to zero
      XPDF(-6) = 0D0
      XPDF( 6) = 0D0
ckr Temporary!

*---  Here one can also call ones own PDF code
*---  --> Only remember that these are MOMENTUM DENSITIES, i.e. x * PDF
*---  --> and the scale is in GeV
C---  CALL MY-FAVORITE-PDFS(....)

      RETURN
      END
