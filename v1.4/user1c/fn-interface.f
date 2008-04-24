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
      DOUBLE PRECISION ALPS_IT,ALPHASPDF,RALPSL,RALPSM
      DOUBLE PRECISION ASMZ,ASMZPDF,A2PI1,A2PI2
      INTEGER IOAS,NLOOP

      DOUBLE PRECISION ZMASS
ckr Old Z mass
ckr      PARAMETER (ZMASS = 91.187D0)
ckr Z mass from PDG 2006
      PARAMETER (ZMASS = 91.1876D0)

ckr PDG 2006
      DOUBLE PRECISION ASMZPDG
      PARAMETER (ASMZPDG = 0.1176D0)

ckr For LHC energies assume always 5 flavours
      INTEGER NF
      PARAMETER (NF = 5)

ckr 30.01.2008: Initialize pi in double precision at first call like elsewhere
ckr      PARAMETER (PI=3.1415927D0)
      INTEGER IFIRST
      DOUBLE PRECISION PI
      DATA IFIRST,PI/0,0.D0/
      SAVE IFIRST,PI
      
      call GetOrderAs(IOAS)
      ASMZPDF = ALPHASPDF(ZMASS)
ckr Round value to 6 digits only
      ASMZPDF = ANINT(ASMZPDF*1D6)/1D6
      IF (IFIRST.EQ.0) THEN
         IFIRST = 1
         PI = 4D0 * ATAN(1D0)
         WRITE(*,*)"FNALPHAS: PI =",PI 
         WRITE(*,*)"FNALPHAS: M_Z/GeV =",ZMASS 
         WRITE(*,*)"FNALPHAS: a_s(M_Z)_PDG (parameter) =",ASMZPDG 
         WRITE(*,*)"FNALPHAS: a_s(M_Z)_PDF (1st call)  =",ASMZPDF 
         WRITE(*,*)"FNALPHAS: a_s was used in",IOAS+1,
     >        "-loop order in PDF"
      ENDIF

c === example: exact, iterative solution of the 2-loop RGE 
c      nloop=2
c      ASMZ=0.118              ! set here the value of alpha_s(Mz)
c      ASMZ=0.1185             ! for H1-2000 MSbar
c      ASMZ=0.1205             ! for MRST2004
ckr 30.01.2008: Initialize alphas(M_Z) in double precision!
ckr      ASMZ=0.118D0              ! set here the value of alpha_s(Mz)
ckr      FNALPHAS = ALPS_IT(MUR,ASMZ,NLOOP)/2d0/PI
ckr      A2PI1 = ALPS_IT(MUR,ASMZ,NLOOP)/2d0/PI

c === here you can call your own alpha_s code
c           -> remember to divide by 2Pi

      NLOOP = IOAS+1
      ASMZ = ASMZPDF
ckr      ASMZ = ASMZPDG
      FNALPHAS = ALPHASPDF(MUR)/2D0/PI
ckr      FNALPHAS = RALPSL(MUR,ZMASS,ASMZ,NF,NLOOP)/2D0/PI
ckr      FNALPHAS = RALPSM(MUR,ZMASS,ASMZ,NF,NLOOP)/2D0/PI
ckr      A2PI2 = RALPSM(MUR,ZMASS,ASMZ,NF,NLOOP)/2D0/PI
ckr      write(*,*)"d(A2PI) =",abs(a2pi2-a2pi1)/a2pi1
ckr      FNALPHAS = A2PI2

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
