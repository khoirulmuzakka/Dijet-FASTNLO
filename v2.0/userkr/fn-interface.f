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

ckr
      CHARACTER*255 ASMODE
      DOUBLE PRECISION ASMZVAL
      INTEGER IASLOOP
      COMMON/STEER/ASMZVAL,IASLOOP,ASMODE
ckr

ckr
      DOUBLE PRECISION ALPHASPDF,RALPSM,PYALPS
      DOUBLE PRECISION AS,ASMZ,ASMZPDF,QLAM4,QLAM5
      INTEGER IOAS
ckr

ckr Z mass and alpha_s from PDG 2012
      DOUBLE PRECISION ZMASS
      PARAMETER (ZMASS = 91.1876D0)
      DOUBLE PRECISION ASMZPDG
      PARAMETER (ASMZPDG = 0.1184D0)
ckr

ckr For LHC energies assume 5 flavours
      INTEGER NF
      PARAMETER (NF = 5)
ckr

      DATA IFIRST/0/

*---  Get info from PDF set
      CALL GETORDERAS(IOAS)
      CALL GETLAM4(0,QLAM4)
      CALL GETLAM5(0,QLAM5)
      NLOOP = IOAS+1
      ASMZPDF = ALPHASPDF(ZMASS)
ckr Round value to 6 digits only for comparisons
ckr      ASMZPDF = ANINT(ASMZPDF*1D6)/1D6

*---  Print info
      IF (IFIRST.EQ.0) THEN
         IFIRST = 1
         PI = 4D0 * ATAN(1D0)
         WRITE(*,*)"FNALPHAS: PI =",PI
         WRITE(*,*)"FNALPHAS: M_Z_PDG/GeV =",ZMASS
         WRITE(*,*)"FNALPHAS: a_s(M_Z)_PDG (parameter) =",ASMZPDG
         WRITE(*,*)"FNALPHAS: a_s(M_Z)_PDF (1st call)  =",ASMZPDF
         WRITE(*,*)"FNALPHAS: Lambda_4_PDF (1st call)  =",QLAM4
         WRITE(*,*)"FNALPHAS: Lambda_5_PDF (1st call)  =",QLAM5
         IF (ASMZVAL.GT.0.D0) THEN
            IF (ASMODE.EQ."PY") THEN
               WRITE(*,*)"FNALPHAS: Lambda_4 requested       =",ASMZVAL
            ELSE
               WRITE(*,*)"FNALPHAS: a_s(M_Z) requested       =",ASMZVAL
            ENDIF
         ENDIF
         WRITE(*,*)"FNALPHAS: a_s was used in",NLOOP,
     >        "-loop order in PDF"
         IF (IASLOOP.GT.0) THEN
            WRITE(*,*)"FNALPHAS: a_s requested to be in",IASLOOP,
     >           "-loop order"
         ENDIF
      ENDIF

*---  ASMODE "PDF": Take alpha_s etc. from PDF set
      IF (ASMODE.EQ."PDF") THEN
         IF (ASMZVAL.GT.0.D0.OR.IASLOOP.GT.0) THEN
            WRITE(*,*)"fastNLO: ERROR! In a_s mode PDF "//
     >           "alpha_s(M_Z) and the loop order cannot be changed!"
            STOP
         ENDIF
         AS = ALPHASPDF(MUR)

*---  ASMODE "PY":
*     Calculation of the running strong coupling up to 2-loop order
*     as in PYTHIA 6.4 using Lambda ...
      ELSEIF (ASMODE.EQ."PY") THEN
         ASMZ = QLAM4
         IF (ASMZVAL.GT.0.D0) THEN
            ASMZ = ASMZVAL
         ENDIF
*---  Only NLOOP=0,1,2 allowed (0 = fixed value)
         IF (IASLOOP.EQ.0.OR.IASLOOP.EQ.1.OR.IASLOOP.EQ.2) THEN
            NLOOP = IASLOOP
         ELSEIF (IASLOOP.GT.0) THEN
            WRITE(*,*)"fastNLO: ERROR! In a_s mode PY "//
     >           "only loop orders 1 and 2 are possible!"
            STOP
         ENDIF
         AS = PYALPS(MUR,ZMASS,ASMZ,4,NLOOP)

*---  ASMODE "KR":
*     Calculation of the running strong coupling up to 3-loop order
*     in the MSbar scheme for given alpha_s(M_Z) according to Giele,
*     Glover, Yu: hep-ph/9506442.
      ELSEIF (ASMODE.EQ."KR") THEN
         ASMZ = ASMZPDF
         IF (ASMZVAL.GT.0.D0) THEN
            ASMZ = ASMZVAL
         ENDIF
*---  Only NLOOP=1,2,3 allowed
         IF (IASLOOP.EQ.1.OR.IASLOOP.EQ.2.OR.IASLOOP.EQ.3) THEN
            NLOOP = IASLOOP
         ELSEIF (IASLOOP.GT.0) THEN
            WRITE(*,*)"fastNLO: ERROR! In a_s mode KR "//
     >           "only loop orders 1,2 and 3 are possible!"
            STOP
         ENDIF
         AS = RALPSM(MUR,ZMASS,ASMZ,NF,NLOOP)

*---  ASMODE "MW14": Use Markus original example code v14
*---  Exact, iterative solution of the RGE at 2-loop or 4-loop order
*---  following GRV hep-ph/9806404
      ELSEIF (ASMODE.EQ."MW14") THEN
         ASMZ = ASMZPDF
         IF (ASMZVAL.GT.0.D0) THEN
            ASMZ = ASMZVAL
         ENDIF
*---  Only NLOOP=2 or 4 allowed
         IF (IASLOOP.EQ.2.OR.IASLOOP.EQ.4) THEN
            NLOOP = IASLOOP
         ELSEIF (IASLOOP.GT.0) THEN
            WRITE(*,*)"fastNLO: ERROR! In a_s mode MW14 "//
     >           "only loop orders 2 and 4 are possible!"
            STOP
         ENDIF
         AS = ALPS_IT_FNLO14(MUR,ASMZ,NLOOP)

*---  ASMODE "MW": Use Markus example code v20
*---  Exact, iterative solution of the RGE at 2-, 3- or 4-loop order
*---  following GRV hep-ph/9806404
      ELSEIF (ASMODE.EQ."MW") THEN
         ASMZ = ASMZPDF
         IF (ASMZVAL.GT.0.D0) THEN
            ASMZ = ASMZVAL
         ENDIF
*---  NLOOP=2, 3 or 4 allowed
         IF (IASLOOP.GE.2.AND.IASLOOP.LE.4) THEN
            NLOOP = IASLOOP
         ELSEIF (IASLOOP.GT.0) THEN
            WRITE(*,*)"fastNLO: ERROR! In a_s mode MW "//
     >           "only loop orders 2, 3 or 4 are possible!"
            STOP
         ENDIF
         AS = ALPS_IT(MUR,ASMZ,NLOOP)

      ELSE
*---  Call your own alpha_s code
         WRITE(*,*)"fastNLO: ERROR! Own alpha_s code called, "//
     >        "but not implemented! ASMODE =", ASMODE
         STOP
      ENDIF
      FNALPHAS = AS/2D0/PI

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

*---  Temporary fix for lhapdf-5.8.7 and lhapdf-5.8.8 bug with ABKM09 or ABM11
*---  Set top-antitop PDF to zero
*---  Fixed in LHAPDF-5.8.9b1
C      XPDF(-6) = 0D0
C      XPDF( 6) = 0D0


*---  Here one can also call ones own PDF code
*---  --> Only remember that these are MOMENTUM DENSITIES, i.e. x * PDF
*---  --> and the scale is in GeV
C---  CALL MY-FAVORITE-PDFS(....)

      RETURN
      END
