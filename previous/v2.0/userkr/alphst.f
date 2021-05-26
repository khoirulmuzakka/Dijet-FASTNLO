      DOUBLE PRECISION FUNCTION RALPSM(MUR,ZMASS,ALPSMZ,NF,NLOOP)
************************************************************************
*
* RALPSM
* Calculation of the running strong coupling up to 3-loop order
* in the MSbar scheme for given alpha_s(M_Z) according to Giele,
* Glover, Yu: hep-ph/9506442.
*
* INPUT:    MUR   : scale at which to evaluate alphas
*           ZMASS : Z boson mass
*           ALPSMZ: alphas at the Z mass
*           NF    : no. of active flavours
*           NLOOP : loop order
*
* Original version by M. Wobisch, 21/06/97
* Slightly modified by K. Rabbertz, 31/03/99
* Modified and included in NLOLIB by K. Rabbertz, 30/03/2001
* Updated Z mass by K. Rabbertz, 26/03/2008
************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION MUR,ZMASS,ALPSMZ
      INTEGER NF,NLOOP

      DOUBLE PRECISION BETA

*---PDG 2012
      DOUBLE PRECISION ZMPDG
      PARAMETER (ZMPDG = 91.1876D0)

      DOUBLE PRECISION BPI(0:2),LFAC,LSCALE,PI

      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST/.TRUE./

*---Initialization
      IF (FIRST) THEN
         PI = 4.D0*ATAN(1.D0)
         BPI(0) = 1.D0/2.D0/PI
         BPI(1) = BPI(0)/4.D0/PI
         BPI(2) = BPI(1)/4.D0/PI
      ENDIF

*---alpha_s
      LSCALE = LOG(MUR/ZMASS)
      IF (NLOOP.EQ.1) THEN
         LFAC = BETA(0,NF)*BPI(0) * LSCALE
      ELSEIF (NLOOP.EQ.2) THEN
         LFAC = (BETA(0,NF)*BPI(0) + BETA(1,NF)*BPI(1)*ALPSMZ) * LSCALE
      ELSEIF (NLOOP.EQ.3) THEN
         LFAC = (BETA(0,NF)*BPI(0) + BETA(1,NF)*BPI(1)*ALPSMZ +
     &        BETA(2,NF)*BPI(2)*ALPSMZ*ALPSMZ) * LSCALE -
     &        BETA(0,NF)*BPI(0)*BETA(1,NF)*BPI(1)/2.D0*ALPSMZ*ALPSMZ *
     &        LSCALE*LSCALE
      ELSE
         WRITE(*,*)'RALPSM: Illegal loop order, stopped:',NLOOP
         STOP
      ENDIF

      RALPSM = ALPSMZ / (1.D0 + ALPSMZ * LFAC)

      IF (FIRST) THEN
         FIRST = .FALSE.
         WRITE(*,*)' '
         WRITE(*,*)'********************************************'
         WRITE(*,*)' RALPSM: First call to alpha_s calculation: '
         WRITE(*,*)' '
         WRITE(*,FMT='(A,F8.4,A,F7.4,A)')
     &        '  MZ =',ZMASS,' GeV, alpha_s(MZ) =',ALPSMZ
         WRITE(*,*)' NF =',NF,', NLOOP =',NLOOP
         WRITE(*,FMT='(A,F5.1,A,F11.8)')
     &        '  alpha_s(MUR =',MUR,') =',RALPSM
         WRITE(*,FMT='(A,F8.4,A)')
     &        '  Note: MZ PDG 2012 =',ZMPDG,' GeV'
         WRITE(*,*)'********************************************'
         WRITE(*,*)' '



      ENDIF

      RETURN
      END



      DOUBLE PRECISION FUNCTION RALPSL(MUR,ZMASS,LAMBDA,NF,NLOOP)
************************************************************************
*
* RALPSL
* Calculation of the running strong coupling up to 2-loop order
* in the MSbar scheme for given Lambda according to Ellis,
* Sterling, Webber: QCD and Collider Physics.
*
* INPUT:    MUR   : scale at which to evaluate alphas
*           ZMASS : Z boson mass
*           LAMBDA: Lambda parameter of QCD (GeV!)
*           NF    : no. of active flavours
*           NLOOP : loop order
*
* Included in NLOLIB by K. Rabbertz, 31/03/2001
* Updated Z mass by K. Rabbertz, 26/03/2008
************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION MUR,ZMASS,LAMBDA
      INTEGER NF,NLOOP

      DOUBLE PRECISION BETA

*---PDG 2012
      DOUBLE PRECISION ZMPDG
      PARAMETER (ZMPDG = 91.1876D0)

      DOUBLE PRECISION LSCALE2,PI

      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST/.TRUE./

*---Initialization
      IF (FIRST) THEN
         PI = 4.D0*ATAN(1.D0)
      ENDIF

*---alpha_s
      LSCALE2 = LOG(MUR*MUR/LAMBDA/LAMBDA)
      RALPSL  = 4.D0*PI/BETA(0,NF)/LSCALE2
      IF (NLOOP.EQ.1) THEN
      ELSEIF (NLOOP.EQ.2) THEN
         RALPSL = RALPSL*(1.D0-BETA(1,NF)/BETA(0,NF)/BETA(0,NF) *
     &        LOG(LSCALE2)/LSCALE2)
      ELSEIF (NLOOP.EQ.3) THEN
         WRITE(*,*)'RALPSL: 3-loop not implemented yet, stopped!'
         STOP
      ELSE
         WRITE(*,*)'RALPSL: Illegal loop order, stopped!'
         STOP
      ENDIF

      IF (FIRST) THEN
         FIRST = .FALSE.
         WRITE(*,*)' '
         WRITE(*,*)'********************************************'
         WRITE(*,*)' RALPSL: First call to alpha_s calculation: '
         WRITE(*,*)' '
         WRITE(*,FMT='(A,F8.4,A,F8.4,A)')
     &        ' MZ =',ZMASS,' GeV, Lambda =',LAMBDA,' GeV'
         WRITE(*,*)' NF =',NF,', NLOOP =',NLOOP
         WRITE(*,*)' alpha_s(MUR =',MUR,') =',RALPSL
         WRITE(*,FMT='(A,F8.4,A)')
     &        '  Note: MZ PDG 2012 =',ZMPDG,' GeV'
         WRITE(*,*)'********************************************'
         WRITE(*,*)' '
      ENDIF

      RETURN
      END



      DOUBLE PRECISION FUNCTION LAMASZ(ZMASS,ALPSMZ,NF,NLOOP)
************************************************************************
*
* LAMASZ
* Calculation of Lambda_5 up to 2-loop order
* in the MSbar scheme for given alpha_s(MZ) according to Ellis,
* Sterling, Webber: QCD and Collider Physics.
*
* INPUT:    ZMASS : Z boson mass
*           ALPSMZ: alphas at the Z mass
*           NF    : no. of active flavours
*           NLOOP : loop order
*
* Included in NLOLIB by K. Rabbertz, 02/04/2001
* Updated Z mass by K. Rabbertz, 26/03/2008
************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION ZMASS,ALPSMZ
      INTEGER NF,NLOOP

      DOUBLE PRECISION BETA

*---PDG 2012
      DOUBLE PRECISION ZMPDG
      PARAMETER (ZMPDG = 91.1876D0)

      DOUBLE PRECISION BFAC,BPRIM,EXPON,LFAC,PI

      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST/.TRUE./

*---Initialization
      IF (FIRST) THEN
         PI = 4.D0*ATAN(1.D0)
      ENDIF

*---Lambda_5
      BFAC  = BETA(0,NF)/4.D0/PI
      EXPON = 1.D0/BFAC/ALPSMZ
      LFAC  = 1.D0
      IF (NLOOP.EQ.1) THEN
      ELSEIF (NLOOP.EQ.2) THEN
         BPRIM = BETA(1,NF)/BETA(0,NF)/4.D0/PI
         EXPON = EXPON + BPRIM/BFAC*
     &        LOG(BPRIM*ALPSMZ/(1.D0 + BPRIM*ALPSMZ))
         LFAC  = 1.D0/(BFAC/BPRIM)**(BPRIM/BFAC/2.D0)
      ELSEIF (NLOOP.EQ.3) THEN
         WRITE(*,*)'LAMASZ: 3-loop not implemented yet, stopped!'
         STOP
      ELSE
         WRITE(*,*)'LAMASZ: Illegal loop order, stopped!'
         STOP
      ENDIF
      LAMASZ = LFAC*SQRT(EXP(-EXPON))*ZMASS

      IF (FIRST) THEN
         FIRST = .FALSE.
 1001    FORMAT(A,F8.4,A)
         WRITE(*,*)' '
         WRITE(*,*)'********************************************'
         WRITE(*,*)' LAMASZ: First call to Lambda(alpha_s(MZ)): '
         WRITE(*,*)' '
         WRITE(*,*)' MZ =',ZMASS,' GeV, alpha_s(MZ) =',ALPSMZ
         WRITE(*,*)' NF =',NF,', NLOOP =',NLOOP
         WRITE(*,*)' Lambda_5 =',LAMASZ
         WRITE(*,1001)'  Note: MZ PDG 2012 =',ZMPDG,' GeV'
         WRITE(*,*)'********************************************'
         WRITE(*,*)' '
      ENDIF

      RETURN
      END



      DOUBLE PRECISION FUNCTION BETA(NORD,NF)
************************************************************************
*
* BETA
* Calculation of the beta function of QCD up to 4-loop order
* in the MSbar scheme according to van Ritbergen, Vermaseren, Larin:
* hep-ph/9701390.
* Beware: Definitions of BETA may differ and from 2-loop on depend on
*         the renormalization scheme!
* Here: 1/(4\pi)*\beta(\alpha_s) :=
*       1/(4\pi)*\mu^2*(d\alpha_s)/(d\mu^2) :=
*       - \sum_{n=0}^{\infty} \beta_n * (\alpha_s/(4\pi))^{n+2}
*
* INPUT:    NORD  : order = (loop order - 1)
*           NF    : no. of active flavours
*
* K. Rabbertz, 30/03/2001
************************************************************************
      IMPLICIT NONE

      INTEGER NORD,NF

      DOUBLE PRECISION ZETA3
      PARAMETER (ZETA3 = 1.202056903D0)

      INTEGER NLOOP

      NLOOP = NORD + 1
      IF (NLOOP.EQ.1) THEN
         BETA = 11.D0 - 2.D0/3.D0*NF
      ELSEIF (NLOOP.EQ.2) THEN
         BETA = 102.D0 - 38.D0/3.D0*NF
      ELSEIF (NLOOP.EQ.3) THEN
         BETA = 2857.D0/2.D0 - 5033.D0/18.D0*NF + 325.D0/54.D0*NF*NF
      ELSEIF (NLOOP.EQ.4) THEN
         BETA = (149753.D0/6.D0 + 3564.D0*ZETA3) -
     &        (1078361.D0/162.D0 + 6508.D0/27.D0*ZETA3) * NF +
     &        (50065.D0/162.D0 + 6472.D0/81.D0*ZETA3) *NF*NF +
     &        (1093.D0/729.D0) * NF*NF*NF
      ELSE
         WRITE(*,*)'BETA: Illegal order, stopped:',NORD
         STOP
      ENDIF

      RETURN
      END
