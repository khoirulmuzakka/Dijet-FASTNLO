***********************************************************************
*     
*     Calculation of alpha_s in the MSbar scheme for given alpha_s(Mz)
*     using an exact, iterative solution of 2-/3-/4-loop formulas
*     as used by GRV hep-ph/9806404
*     
*     Initial version: M. Wobisch, 1999
*     Added 3-loop:    M. Wobisch, 2009
*     Updated for v2:  K. Rabbertz, M Wobisch, 2011
*     
*     Current restrictions:
*     ---------------------
*     So far, the code is written for nf=5, and no flavor thresholds 
*     have been implemented. Therefore the code should be used only 
*     for mu_r > m_bottom (which is safe for jet observables)
*     
***********************************************************************
      DOUBLE PRECISION FUNCTION ALPS_IT(MU,ALPSMZ,NLOOP)
      IMPLICIT NONE
      DOUBLE PRECISION MU, ALPSMZ,  FBETA
      DOUBLE PRECISION B0, B1, B2,B3, B10 , ZETA3,
     >     PI4, F, FP,FM, LL2,
     >     ONED, TWOD, ZMASS, ZMASS2, ALPHAS, ASAPPROX, Q2, LAM2, LQ2,
     >     MUCACHE,ASCACHE,ASMZCACHE
      INTEGER  NLOOP, NF, IFIRST, I, NLOOPCACHE
      PARAMETER (ZMASS = 91.1876D0) ! PDG 2011
      DATA IFIRST/0/, ONED/1D0/, TWOD/2D0/
      SAVE IFIRST, NF, ONED, TWOD, PI4, B0, B1, B10, ZMASS2, ASCACHE,
     >     MUCACHE
      INCLUDE 'strings.inc'

*---  Initialize pi and beta functions
      IF (IFIRST.EQ.0) THEN
         IFIRST = 1
         NF = 5
         PI4 = 4D0 * 4D0 * ATAN(1D0)
         B0  = 11D0 - 2D0/3D0 * DBLE(NF)
         B1  = 102D0 - 38D0 / 3D0 * DBLE(NF)
         B10 = B1 / B0 / B0
         ZMASS2 = ZMASS*ZMASS
         ASMZCACHE = 0D0
         MUCACHE = 0D0
         NLOOPCACHE = 0
*---  Print info
         WRITE(*,'(A)')""
         WRITE(*,'(X,A)')CSEPS
         WRITE(*,'(X,A)')"# alphas-grv: First call:"
         WRITE(*,'(X,A)')CSEPS
         WRITE(*,'(X,A,F18.15)')"# ALPHAS-GRV: PI              = ",
     >        PI4/4D0 
         WRITE(*,'(X,A,F9.6)')"# ALPHAS-GRV: M_Z/GeV         = ",ZMASS 
         WRITE(*,'(X,A,F9.6)')"# ALPHAS-GRV: a_s(M_Z)        = ",ALPSMZ 
         WRITE(*,'(X,A,I2)')"# APLHAS-GRV: a_s loop        = ",NLOOP
         WRITE(*,'(X,A,L4)')"# APLHAS-GRV: flavor-matching = ",.FALSE.
         WRITE(*,'(X,A,I2)')"# APLHAS-GRV: nf (M_Z)        = ",NF
         WRITE(*,'(X,A)')CSEPS
      ENDIF

      IF (MU.EQ.MUCACHE .AND. ALPSMZ.EQ.ASMZCACHE 
     >     .AND. NLOOP.EQ.NLOOPCACHE) THEN
         ALPS_IT = ASCACHE
         RETURN
      ENDIF

      Q2 = MU**2

*---  Exact formula -> extract Lambda from alpha_s(Mz)
      LAM2 = ZMASS2 / DEXP(FBETA(ALPSMZ,NLOOP))

*---  Extract approx alpha_s(mu) value - 2 loop approx is fine
      LL2 = ZMASS2 * DEXP( -PI4/B0/ALPSMZ + 
     >     B10 * DLOG( PI4/B0/ALPSMZ + B10) )
      LQ2 = DLOG( Q2 / LL2 )
      ASAPPROX = PI4/B0/LQ2 * (ONED - B10*DLOG(LQ2)/LQ2)
      ALPHAS = ASAPPROX

*---  Exact 4-loop value by Newton procedure
      DO I=1,3
         F  = DLOG(Q2/LAM2) - FBETA(ALPHAS,NLOOP)
         FP = - FBETA(ALPHAS*1.01D0,NLOOP)
         FM = - FBETA(ALPHAS*0.99D0,NLOOP)
         ALPHAS = ALPHAS - F/(FP-FM)*0.02D0*ALPHAS 
      ENDDO

*---  That's it - modify cache - set function - return
      MUCACHE = MU
      ASMZCACHE = ALPSMZ
      ASCACHE = ALPHAS
      ALPS_IT = ALPHAS

      RETURN 
      END
***********************************************************************

      DOUBLE PRECISION FUNCTION FBETA(ALPHAS,NLOOP)
      IMPLICIT NONE
      DOUBLE PRECISION ALPHAS, ALPI, B0, B1, B2,B3, B10,B20,B30,C,
     >     ZETA3, PI, ONED, TWOD
      INTEGER  NLOOP, IFIRST, NF
      DATA IFIRST/0/
      SAVE IFIRST,NF,PI,ONED,TWOD,ZETA3,B0,B1,B2,B3,B10,B20,B30,C

*---  Initialize pi and beta functions
      IF (IFIRST.EQ.0) THEN
         IFIRST = 1
         NF = 5
         PI = 4D0 * ATAN(1D0)
         ONED = 1.0D0
         TWOD = 2.0D0

         ZETA3 = 1.202056903D0
         B0  = (11D0 - 2D0/3D0 * DBLE(NF)) / 4D0
         B1  = (102D0 - 38D0 / 3D0 * DBLE(NF)) /16D0
         B2  = (2857D0/2D0 - 5033D0/18D0*DBLE(NF) + 
     >        325D0/54D0*DBLE(NF)**2) / 64D0
         B3  = ((149753D0/6D0 + 3564D0*ZETA3) - 
     >        (1078361D0/162D0 + 6508D0/27D0*ZETA3) * DBLE(NF) +
     >        (50065D0/162D0 + 6472D0/81D0*ZETA3) * DBLE(NF)**2 +
     >        1093D0/729D0 * DBLE(NF)**3 ) /256D0

         B10 = B1 / B0
         B20 = B2 / B0
         B30 = B3 / B0
         C = B10 / B0 * DLOG(B0) 
      ENDIF

      ALPI = ALPHAS / PI
      IF (NLOOP.EQ.2) THEN      ! 2-loop RGE
         FBETA = C + ONED/B0 * ( 
     >        ONED/ALPI + B10 * DLOG(ALPI) + (-B10**2) * ALPI 
     >        + (B10**3/TWOD)*ALPI**2 )
      ELSEIF (NLOOP.EQ.3) THEN  ! 3-loop RGE
         FBETA = C + ONED/B0 * ( 
     >        ONED/ALPI + B10 * DLOG(ALPI) + (B20-B10**2) * ALPI 
     >        + (B10*B20 + B10**3/TWOD)*ALPI**2 )
      ELSEIF (NLOOP.EQ.4) THEN  ! 4-loop RGE
         FBETA = C + ONED/B0 * ( 
     >        ONED/ALPI + B10 * DLOG(ALPI) + (B20-B10**2) * ALPI 
     >        + (B30/TWOD - B10*B20 + B10**3/TWOD)*ALPI**2 )
      ELSE
         WRITE(*,*)'FBETA: ERROR! Nloop =',nloop,' not implemented.'//
     >        ' Stopped.'
         STOP
      ENDIF
      
      RETURN
      END
