CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  FUNCTION TO RETURN STRONG COUPLING ALPHA_S=G^2/4/PI         C
C  ---->  WITH INTERPOLATION ACROSS THRESHOLDS  <----          C
C  INPUTS:                                                     C
C         IALPHAS                                              C
C                 = 0: CONSTANT (0.2)                          C
C                 = 1: ONE LOOP                                C
C                 = 2: TWO LOOP                                C
C         MU2                                                  C
C                 = RENORMALIZATION SCALE SQUARED              C
C         LQCD5                                                C
C                 = QCD SCALE FOR NF=5                         C
C  OUTPUT:                                                     C
C         DOUBLE PRECISION ALPHAS                              C
C         NF = EFFECTIVE NUMBER OF FLAVORS                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION ALPHA_S(IALPHAS,XMU2,XLQCD5,NEFF)
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z )
      ALPHA_S = 0.2
      NEFF=5
      IF( IALPHAS.EQ.0 )RETURN
      ALPHA_S = ALFAS5(XMU2,XLQCD5,IALPHAS,NEFF)
      RETURN
      END
C.------------------------------------------------------------.C
C. FUNCTION TO RETURN ALPHA_S GIVEN LAMBDA_QCD^5 AND SCALE Q2 .C
C. AND ILOOP=1 (ONE LOOP) OR ILOOP=2 (TWO LOOP)               .C
C. SEE G. ALTARELLI ET AL. NUCL. PHYS. B308 (1988) 724.       .C
C. EQNS (7)->(13)                                             .C
C.------------------------------------------------------------.C
      FUNCTION ALFAS5(Q2,ALQCD5,ILOOP,NEFF)
      IMPLICIT DOUBLE PRECISION(A-Z)
      INTEGER ILOOP,NEFF
      A  = 1.00
c      MB = 4.750
c      mb=5.0
       mb=4.5
c      MC = 1.50
c      mc=1.6
      mc=1.3
      AMB=A*MB
      AMC=A*MC
      MU=SQRT(Q2)
      ALQCD52=ALQCD5*ALQCD5
      IF      (MU.GT.AMB) THEN
         ALFAS5 = ASLOOP(Q2,ALQCD52,5,ILOOP)
C. MUST RETURN NF=5 FOR RENORMALIZATION GROUP TO WORK
         NEFF=5
      ELSE IF (MU.GT.AMC) THEN
         AMB2=AMB*AMB
         ALFAS5 = 1 / ( 1 / ASLOOP(Q2  ,ALQCD52,4,ILOOP) + 
     &                    1 / ASLOOP(AMB2,ALQCD52,5,ILOOP) -
     &                    1 / ASLOOP(AMB2,ALQCD52,4,ILOOP) )
C. MUST RETURN NF=4 FOR RENORMALIZATION GROUP TO WORK
         NEFF=4
      ELSE
         AMB2=AMB*AMB
         AMC2=AMC*AMC
         ALFAS5 = 1 / ( 1 / ASLOOP(Q2  ,ALQCD52,3,ILOOP) + 
     &                    1 / ASLOOP(AMC2,ALQCD52,4,ILOOP) +
     &                    1 / ASLOOP(AMB2,ALQCD52,5,ILOOP) -
     &                    1 / ASLOOP(AMB2,ALQCD52,4,ILOOP) - 
     &                    1 / ASLOOP(AMC2,ALQCD52,3,ILOOP) )
C. MUST RETURN NF=3 FOR RENORMALIZATION GROUP TO WORK
         NEFF=3
      ENDIF
      RETURN
      END
C.---------------------------------------------------------------.C
C. FUNCTION TO RETURN ALPHAS GIVEN MU2,LQCD2,NF AND NORDER       .C
C.---------------------------------------------------------------.C
      FUNCTION ASLOOP(XMU2,XLQCD2,NF,NORDER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.14159)
      IF    ( NORDER.EQ.1 ) THEN
         B1=(33.-2.*NF)/(12.*PI)
         T=XMU2/XLQCD2
         ASLOOP=1/(B1*LOG(T))
      ELSEIF( NORDER.EQ.2 ) THEN
         B1=(33.-2.*NF)/(12.*PI)
         B2=(153.-19.*NF)/(2.*PI*(33.-2.*NF))
         T=XMU2/XLQCD2
         F1=B1*LOG(T)
         ASLOOP=(1.-B2*LOG(LOG(T))/F1)/F1
      ELSE
         WRITE(*,*) 'ERROR IN ASLOOP: NORDER = ',NORDER
         STOP
      ENDIF
      RETURN
      END
