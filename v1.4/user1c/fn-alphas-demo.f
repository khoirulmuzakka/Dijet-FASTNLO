*******************************************************************
* M. Wobisch  25/05/99
*
* calculation of alpha_s in the MSbar scheme
* for given alpha_s(Mz)
*
* using exact / iterative solution of 2-loop formula
*
* as GRV hep-ph/9806404
*
*******************************************************************
      DOUBLE PRECISION FUNCTION ALPS_IT(MU,ALPSMZ,NLOOP)
*
*  1st version - only for 2-loop - as GRV code
*
      IMPLICIT NONE
      DOUBLE PRECISION MU, ALPSMZ, ALPS4_IT
      DOUBLE PRECISION B0, B1, B10 , PI4, F, FP,FM,
     +     ZMASS, ZMASS2, ALPHAS, ASAPPROX, Q2, LAM2, LQ2
      INTEGER  NLOOP, NF, IFIRST, I
ckr 30.01.2008: Initialize Z mass in double precision 
ckr      PARAMETER (ZMASS = 91.187)        ! PDG data book '98
      PARAMETER (ZMASS = 91.187D0)
      SAVE IFIRST, NF, PI4, B0, B1, B10, ZMASS2
      DATA IFIRST/0/

c - other routine for 4-loop
      IF (NLOOP .eq. 4) THEN
         ALPS_IT =  ALPS4_IT(MU,ALPSMZ,NLOOP)
         RETURN
      ELSEIF (NLOOP .eq. 3) THEN
         WRITE(*,*) ' 3-loop alpha_s not available !!!  '
         STOP
      ENDIF

c - initialize pi and beta functions
      IF (IFIRST.eq.0) THEN
         IFIRST = 1
c         WRITE(*,*) '  *   ALPS_IT:  exact 2-loop result for alpha_s'
         NF = 5
         PI4 = 4D0 * 4D0 * ATAN(1D0)
         B0  = 11D0 - 2D0/3D0 * DBLE(NF)
         B1  = 102D0 - 38D0 / 3D0 * DBLE(NF)
         B10 = B1 / B0 / B0
         ZMASS2 = ZMASS**2
         IF (NLOOP .ne. 2) WRITE(*,*) 'ALPS_IT:  only for 2-loop!!'
      ENDIF

c - exact formula to extract Lambda from alpha_s(Mz)
      Q2 = MU**2
      LAM2 = ZMASS2 * EXP( -PI4/B0/ALPSMZ + 
     +     B10 * DLOG( PI4/B0/ALPSMZ + B10) )

c - extract approx. alpha_s(mu) value 
      LQ2 = DLOG( Q2 / LAM2 ) 
      ASAPPROX = PI4/B0/LQ2 * (1D0 - B10*DLOG(LQ2)/LQ2)
      ALPHAS = ASAPPROX

c - exact 2loop value by Newton procedure
      DO I=1,6
         F  = LQ2 - PI4/B0/ALPHAS + B10*DLOG(PI4/B0/ALPHAS + B10)
         FP = - PI4/B0/(ALPHAS*1.01D0) + 
     +        B10 * DLOG(PI4/B0/(ALPHAS*1.01D0) + B10)
         FM = - PI4/B0/(ALPHAS*0.99D0) + 
     +        B10 * DLOG(PI4/B0/(ALPHAS*0.99D0) + B10)
         ALPHAS = ALPHAS - F/(FP-FM)*0.02D0*ALPHAS 
c      WRITE(*,*) ' LAMDA/a_s_approx/a_s = ',sqrt(lam2),ASAPPROX,ALPHAS
      ENDDO

c - that's it!
      ALPS_IT = ALPHAS

      RETURN 
      END

C --------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION ALPS4_IT(MU,ALPSMZ,NLOOP)
*
*  2nd version - for 4-loop RGE
*
      IMPLICIT NONE
      DOUBLE PRECISION MU, ALPSMZ,  FBETA
      DOUBLE PRECISION B0, B1, B10,
     +     PI4, F, FP,FM, LL2,
     +     ONED, ZMASS, ZMASS2, ALPHAS, ASAPPROX, Q2, LAM2, LQ2,
     +     MUCACHE,ASCACHE,ASMZCACHE
      INTEGER  NLOOP, NF, IFIRST, I
      PARAMETER (ZMASS = 91.187)        ! PDG data book '98
      SAVE IFIRST, NF, ONED, PI4, B0, B1, B10, ZMASS2, ASCACHE,
     +     MUCACHE
      DATA IFIRST/0/, ONED/1.D0/

c - initialize pi and beta functions
      IF (IFIRST.eq.0) THEN
         IFIRST = 1
         WRITE(*,*) '  *   ALPS_IT:  exact 4-loop result for alpha_s'
         WRITE(*,*) '  *                              in 3 iterations'
         NF = 5
         PI4 = 4D0 * 4D0 * ATAN(1D0)
         B0  = 11D0 - 2D0/3D0 * DBLE(NF)
         B1  = 102D0 - 38D0 / 3D0 * DBLE(NF)
         B10 = B1 / B0 / B0
         ZMASS2 = ZMASS**2
         IF (NLOOP .ne. 4) WRITE(*,*) 'ALPS_IT:  only for 4-loop!!'
         WRITE(*,*) '  *             1st call   mu , alpha_s(Mz)'
         WRITE(*,*) '  *                   ',real(mu),real(alpsmz)
         ASMZCACHE = 0d0
         MUCACHE = 0d0
      ENDIF

      IF (MU .eq. MUCACHE  .and.  ALPSMZ .eq. ASMZCACHE ) THEN
         ALPS4_IT = ASCACHE
         RETURN
      ENDIF

      Q2 = MU**2

c - exact formula -> extract Lambda from alpha_s(Mz)
      LAM2 = ZMASS2 / DEXP(FBETA(ALPSMZ))

c - extract approx alpha_s(mu) value - 2 loop approx is fine
      LL2 = ZMASS2 * DEXP( -PI4/B0/ALPSMZ + 
     +     B10 * DLOG( PI4/B0/ALPSMZ + B10) )
      LQ2 = DLOG( Q2 / LL2 )
      ASAPPROX = PI4/B0/LQ2 * (ONED - B10*DLOG(LQ2)/LQ2)
      ALPHAS = ASAPPROX

c - exact 4-loop value by Newton procedure
      DO I=1,3
         F  = DLOG(Q2/LAM2) - FBETA(ALPHAS)
         FP = - FBETA(ALPHAS*1.01D0)
         FM = - FBETA(ALPHAS*0.99D0)
         ALPHAS = ALPHAS - F/(FP-FM)*0.02D0*ALPHAS 
c         WRITE(*,*) ' i,alphas,q2 = ',i,alphas,real(q2)
      ENDDO

c - that's it - modify cache - set function - return
      MUCACHE = MU
      ASMZCACHE = ALPSMZ
      ASCACHE = ALPHAS
      ALPS4_IT = ALPHAS

      RETURN 
      END

C ------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION FBETA(ALPHAS)
      IMPLICIT NONE
      DOUBLE PRECISION ALPHAS, ALPI, B0, B1, B2,B3, B10,B20,B30,C,
     +     ZETA3, PI, ONED, TWOD
      INTEGER  IFIRST, NF
      SAVE IFIRST, NF, PI, ONED, TWOD, ZETA3,B0,B1,B2,B3,B10,B20,B30,C
      DATA IFIRST/0/

c - initialize pi and beta functions
      IF (IFIRST.eq.0) THEN
         IFIRST = 1
         NF = 5
         PI = 4D0 * ATAN(1D0)
         ONED = 1.0D0
         TWOD = 2.0D0

         ZETA3 = 1.202056903D0
         B0  = (11D0 - 2D0/3D0 * DBLE(NF)) / 4D0
         B1  = (102D0 - 38D0 / 3D0 * DBLE(NF)) /16D0
         B2  = (2857D0/2D0 - 5033D0/18D0*DBLE(NF) + 
     +        325D0/54D0*DBLE(NF)**2) / 64D0
         B3  = ((149753D0/6D0 + 3564D0*ZETA3) - 
     +        (1078361D0/162D0 + 6508D0/27D0*ZETA3) * DBLE(NF) +
     +        (50065D0/162D0 + 6472D0/81D0*ZETA3) * DBLE(NF)**2 +
     +        1093D0/729D0 * DBLE(NF)**3 ) /256D0

         B10 = B1 / B0
         B20 = B2 / B0
         B30 = B3 / B0
         C = B10 / B0 * DLOG(B0) 
      ENDIF

      ALPI = ALPHAS / PI
      FBETA = C + ONED/B0 * ( 
     +     ONED/ALPI + B10 * DLOG(ALPI) + (B20-B10**2) * ALPI 
     +     + (B30/TWOD - B10*B20 + B10**3/TWOD)*ALPI**2 )

      RETURN
      END
