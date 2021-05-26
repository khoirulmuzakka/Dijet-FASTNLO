ckr Based on copy from Pythia Version 6.4.2, 13.05.2009
C...PYALPS
C...Gives the value of alpha_strong.

C... ALPS = PYALPS(Q2)
C... Purpose: to calculate the running strong coupling constant αs,
C...          e.g. in matrix elements and resonance decay widths.
C...          (The function is not used in parton showers, how-
C...           ever, where formulae rather are written in terms of
C...           the relevant Λ values.)
C...          The ﬁrst- and second-order expressions are given by
C...          eqs. (27) and (32).
C...          See MSTU(111) - MSTU(118) and PARU(111) - PARU(118) for options.
C... Q2 :    the momentum transfer scale Q2 at which to evaluate αs.

C... MSTU(111) : (I, D=1) order of αs evaluation in the PYALPS function.
C...             Is overwritten in PYEEVT, PYONIA or PYINIT calls with
C...             the value desired for the process under study.
C...   = 0 : αs is ﬁxed at the value PARU(111).
C...         As extra safety, Λ =PARU(117) is set in PYALPS so that
C...         the ﬁrst-order running αs agrees with the desired ﬁxed αs
C...         for the Q2 value used.
C...   = 1 : ﬁrst-order running αs is used.
C...   = 2 : second-order running αs is used.
C... MSTU(112) : (D = 5) the nominal number of ﬂavours assumed in
C...             the αs expression, with respect to which Λ is deﬁned.
C... MSTU(113) : (D = 3) minimum number of ﬂavours that may be assumed in
C...             αs expression, see MSTU(112).
C... MSTU(114) : (D = 5) maximum number of ﬂavours that may be assumed in
C...             αs expression, see MSTU(112).
C... MSTU(115) : (D = 0) treatment of αs singularity for Q2 → 0 in
C...             PYALPS calls. (Relevant e.g. for QCD 2 → 2 matrix elements
C...             in the p⊥ → 0 limit, but not for showers,
C...             where PYALPS is not called.)
C...   = 0 : allow it to diverge like 1/ ln(Q2 /Λ2 ).
C...   = 1 : soften the divergence to 1/ ln(1 + Q2 /Λ2 ).
C...   = 2 : freeze Q2 evolution below PARU(114), i.e. the eﬀective argument is
C...            max(Q2 ,PARU(114)).
C... MSTU(118) : (I) number of ﬂavours nf found and used in latest PYALPS call.

C... PARU(111) : (D = 0.20) ﬁx αs value assumed in PYALPS when MSTU(111) = 0
C...             (and also in parton showers when αs is assumed ﬁx there).
C... PARU(112) : (I, D=0.25 GeV) Λ used in running αs expression in PYALPS.
C...             Like MSTU(111), this value is overwritten by the calling
C...             physics routines, and is therefore purely nominal.
C... PARU(113) : (D = 1.) the ﬂavour thresholds, for the eﬀective number of
C...             ﬂavours nf to use in the αs expression, are assumed to sit
C...             at Q2 =PARU(113)×m2 , where mqq is the quark mass.
C...             May be overwritten from the calling physics routine.
C... PARU(114) : (D = 4 GeV2 ) Q2 value below which the αs value is
C...             assumed constant for MSTU(115) = 2.
C... PARU(115) : (D = 10.) maximum αs value that PYALPS will ever return;
C...             is used as a last resort to avoid singularities.
C... PARU(117) : (I) Λ value (associated with MSTU(118) eﬀective ﬂavours)
C...             obtained in latest PYALPS call.
C... PARU(118) : (I) αs value obtained in latest PYALPS call.

C... Scale is given by MSTP(32)



ckr      FUNCTION PYALPS(Q2)
      FUNCTION PYALPS(MUR,ZMASS,LAM5,MSTU112,NLOOP)

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      DOUBLE PRECISION MUR,ZMASS,LAM5
      INTEGER NF,NLOOP
C...Commonblocks.
ckr      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
ckr      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
ckr      SAVE /PYDAT1/,/PYDAT2/
C...Coefficients for second-order threshold matching.
C...From W.J. Marciano, Phys. Rev. D29 (1984) 580.
      DIMENSION STEPDN(6),STEPUP(6)
C...     DATA STEPDN/0D0,0D0,(2D0*107D0/2025D0),(2D0*963D0/14375D0),
C...    &(2D0*321D0/3703D0),0D0/
C...     DATA STEPUP/0D0,0D0,0D0,(-2D0*107D0/1875D0),
C...    &(-2D0*963D0/13225D0),(-2D0*321D0/3381D0)/
      DATA STEPDN/0D0,0D0,0.10568D0,0.13398D0,0.17337D0,0D0/
      DATA STEPUP/0D0,0D0,0D0,-0.11413D0,-0.14563D0,-0.18988D0/

ckr Translate my parms into PYTHIA parms
ckr Initialize pi in double precision at first call
      INTEGER IFIRST
      DOUBLE PRECISION PARU1
      SAVE IFIRST,PARU1,PARU2
      DATA IFIRST,PARU1,PARU2/0,0D0,0D0/

      INTEGER MSTU111,MSTU112,MSTU113,MSTU114,MSTU115,MSTU118
      DOUBLE PRECISION PARU111,PARU112,PARU114,PARU115,PARU117,PARU118
      DOUBLE PRECISION PMAS(6,1)
      DATA MSTU113,MSTU114,MSTU115/3,5,0/
      DATA PARU111,PARU113,PARU114,PARU115/0.2D0,1.D0,4.D0,10.D0/
      DATA (PMAS(I,1),I=1,6)/2*0.33D0,0.5D0,1.5D0,4.8D0,175D0/

      IF (IFIRST.EQ.0) THEN
         IFIRST = 1
         PARU1 = 4D0 * ATAN(1D0)
         PARU2 = 2D0 * PARU1
      ENDIF
      Q2 = MUR*MUR
      PARU112 = LAM5
      MSTU111 = NLOOP
ckr Avoid overwriting input parameter
ckr      MSTU112 = NF
      NF = MSTU112

C...Constant alpha_strong trivial. Pick artificial Lambda.
      IF (MSTU111.LE.0) THEN
         PYALPS=PARU111
         MSTU118=MSTU112
         PARU117=0.2D0
         IF(Q2.GT.0.04D0) PARU117=SQRT(Q2)*EXP(-6D0*PARU1/
     &        ((33D0-2D0*MSTU112)*PARU111))
         PARU118=PARU111
         RETURN
      ENDIF

C...Find effective Q2, number of flavours and Lambda.
      Q2EFF=Q2
      IF(MSTU115.GE.2) Q2EFF=MAX(Q2,PARU114)
ckr      NF=MSTU112
      ALAM2=PARU112**2
 100  IF(NF.GT.MAX(3,MSTU113)) THEN
         Q2THR=PARU113*PMAS(NF,1)**2
         IF(Q2EFF.LT.Q2THR) THEN
            NF=NF-1
            Q2RAT=Q2THR/ALAM2
            ALAM2=ALAM2*Q2RAT**(2D0/(33D0-2D0*NF))
            IF(MSTU111.EQ.2) ALAM2=ALAM2*LOG(Q2RAT)**STEPDN(NF)
            GOTO 100
         ENDIF
      ENDIF
 110  IF(NF.LT.MIN(6,MSTU114)) THEN
         Q2THR=PARU113*PMAS(NF+1,1)**2
         IF(Q2EFF.GT.Q2THR) THEN
            NF=NF+1
            Q2RAT=Q2THR/ALAM2
            ALAM2=ALAM2*Q2RAT**(-2D0/(33D0-2D0*NF))
            IF(MSTU111.EQ.2) ALAM2=ALAM2*LOG(Q2RAT)**STEPUP(NF)
            GOTO 110
         ENDIF
      ENDIF
      IF(MSTU115.EQ.1) Q2EFF=Q2EFF+ALAM2
      PARU117=SQRT(ALAM2)

C...Evaluate first or second order alpha_strong.
      B0=(33D0-2D0*NF)/6D0
      ALGQ=LOG(MAX(1.0001D0,Q2EFF/ALAM2))
      IF(MSTU111.EQ.1) THEN
         PYALPS=MIN(PARU115,PARU2/(B0*ALGQ))
      ELSE
         B1=(153D0-19D0*NF)/6D0
         PYALPS=MIN(PARU115,PARU2/(B0*ALGQ)*(1D0-B1*LOG(ALGQ)/
     &        (B0**2*ALGQ)))
      ENDIF
      MSTU118=NF
      PARU118=PYALPS

      RETURN
      END
