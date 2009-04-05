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

      CHARACTER*255 ASMODE
      DOUBLE PRECISION ASMZVAL
      INTEGER IASLOOP
      COMMON/STEER/ASMZVAL,IASLOOP,ASMODE
      
      DOUBLE PRECISION ALPS_IT,ALPHASPDF,RALPSM
      DOUBLE PRECISION AS,ASMZ,ASMZPDF
      INTEGER IOAS,NLOOP

ckr Set default values
      DOUBLE PRECISION ZMASS
ckr Old Z mass
ckr      PARAMETER (ZMASS = 91.187D0)
ckr Z mass from PDG 2006
      PARAMETER (ZMASS = 91.1876D0)

ckr PDG 2006
      DOUBLE PRECISION ASMZPDG
      PARAMETER (ASMZPDG = 0.1176D0)

ckr For LHC energies assume 5 flavours
      INTEGER NF
      PARAMETER (NF = 5)

ckr 30.01.2008: Initialize pi in double precision at first call
ckr      PARAMETER (PI=3.1415927D0)
      INTEGER IFIRST
      DOUBLE PRECISION PI
      SAVE IFIRST,PI
      DATA IFIRST,PI/0,0.D0/

c - Get info from PDF set      
      CALL GETORDERAS(IOAS)
      NLOOP = IOAS+1
      ASMZPDF = ALPHASPDF(ZMASS)
ckr Round value to 6 digits only
      ASMZPDF = ANINT(ASMZPDF*1D6)/1D6

c - Print info
      IF (IFIRST.EQ.0) THEN
         IFIRST = 1
         PI = 4D0 * ATAN(1D0)
         WRITE(*,*)"FNALPHAS: PI =",PI 
         WRITE(*,*)"FNALPHAS: M_Z_PDG/GeV =",ZMASS 
         WRITE(*,*)"FNALPHAS: a_s(M_Z)_PDG (parameter) =",ASMZPDG 
         WRITE(*,*)"FNALPHAS: a_s(M_Z)_PDF (1st call)  =",ASMZPDF 
         IF (ASMZVAL.GT.0.D0) THEN
            WRITE(*,*)"FNALPHAS: a_s(M_Z) requested       =",ASMZVAL 
         ENDIF
         WRITE(*,*)"FNALPHAS: a_s was used in",NLOOP,
     >        "-loop order in PDF"
         IF (IASLOOP.GT.0) THEN
            WRITE(*,*)"FNALPHAS: a_s requested to be in",IASLOOP,
     >           "-loop order" 
         ENDIF
      ENDIF

c - ASMODE "PDF": Take alpha_s etc. from PDF set
      IF (ASMODE.EQ."PDF") THEN 
         IF (ASMZVAL.GT.0.D0.OR.IASLOOP.GT.0) THEN
            WRITE(*,*)"fastNLO: ERROR! In a_s mode PDF "//
     >           "alpha_s(M_Z) and the loop order cannot be changed!"
            STOP
         ENDIF
         AS = ALPHASPDF(MUR)

c - ASMODE "MW": Use Markus original example code
c -              Exact, iterative solution of the 2-loop RGE 
      ELSEIF (ASMODE.EQ."MW") THEN
c   ASMZ = 0.1185             ! for H1-2000 MSbar
c   ASMZ = 0.1205             ! for MRST2004
         ASMZ = ASMZPDG
         IF (ASMZVAL.GT.0.D0) THEN
            ASMZ = ASMZVAL
         ENDIF
c - Only NLOOP=2 or 4 allowed
         NLOOP = 2
         IF (IASLOOP.EQ.2.OR.IASLOOP.EQ.4) THEN
            NLOOP = IASLOOP
         ELSEIF (IASLOOP.GT.0) THEN
            WRITE(*,*)"fastNLO: ERROR! In a_s mode MW "//
     >           "only loop orders 2 and 4 are possible!"
            STOP
         ENDIF
         AS = ALPS_IT(MUR,ASMZ,NLOOP)

c - ASMODE "KR":
c Calculation of the running strong coupling up to 3-loop order
c in the MSbar scheme for given alpha_s(M_Z) according to Giele,
c Glover, Yu: hep-ph/9506442.
      ELSEIF (ASMODE.EQ."KR") THEN
         ASMZ = ASMZPDG
         IF (ASMZVAL.GT.0.D0) THEN
            ASMZ = ASMZVAL
         ENDIF
c - Only NLOOP=1,2,3 allowed
         NLOOP = 2
         IF (IASLOOP.EQ.1.OR.IASLOOP.EQ.2.OR.IASLOOP.EQ.3) THEN
            NLOOP = IASLOOP
         ELSEIF (IASLOOP.GT.0) THEN
            WRITE(*,*)"fastNLO: ERROR! In a_s mode KR "//
     >           "only loop orders 1,2 and 3 are possible!"
            STOP
         ENDIF
         AS = RALPSM(MUR,ZMASS,ASMZ,NF,NLOOP)
      ELSE
c - Call your own alpha_s code
         WRITE(*,*)"fastNLO: ERROR! Own alpha_s code called, "//
     >        "but not implemented! ASMODE =", ASMODE
         STOP
      ENDIF
      FNALPHAS = AS/2D0/PI

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
