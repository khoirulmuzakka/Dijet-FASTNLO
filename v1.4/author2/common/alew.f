
      function xalfaem(q2)
c Alpha_em(MSbar) at the scale q2 = q^2.
c Uses alpha_Thomson below the electron mass, an interpolation between
c m_e and m_tau, and the evolution equation above m_b
c-------------------------------------------------------------------------
      implicit real*8 (a-z)
      integer npoints,ideg
      parameter (npoints=3,ideg=3)
      real*4 ooa(npoints),xlogmu(npoints),divdif
c 1/alpha_em at m_e=0.000511,m_mu=0.1056,m_tau=1.777
      data ooa     / 137.036, 135.95, 133.513 /
c logs of sqrt(q2) at m_e=0.000511,m_mu=0.1056,m_tau=1.777
      data xlogmu  / -7.57914, -2.2481, 0.574927 /
      data zm/91.2d0/,ooaz/127.934d0/,pi/3.141592653589793D0/,nc/3/
c
      if(q2.lt.exp(2.*xlogmu(1))) then
         xalfaem = 1.d0/ooa(1)
      elseif(q2.lt.exp(2.*xlogmu(3))) then
         xlogq = log(q2)/2.d0
         xalfaem = 1.d0/divdif(ooa,xlogmu,npoints,sngl(xlogq),ideg)
      elseif(q2.lt.5.**2) then
         b = 3 + 2*nc*(1d0/3d0)**2 + 2*nc*(2d0/3d0)**2
         xlq = log(q2) - 2.*xlogmu(3)
         xalfaem = 1d0/ooa(3)/(1.d0 - 1.d0/3.d0/pi/ooa(3)*b*xlq)
      else
         b = 3 + 3*nc*(1d0/3d0)**2 + 2*nc*(2d0/3d0)**2
         xlq = log(q2/zm**2)
         xalfaem = 1d0/ooaz/(1.d0 - 1.d0/3.d0/pi/ooaz*b*xlq)
      endif
      return
      end


c From CERN libraries, version 2000. A call to the CERN routine ABEND
c has been substituted with stop (ABEND stops the program gracefully)
c #include "kernnum/pilot.h" --> this is only needed in libraries
      FUNCTION DIVDIF(F,A,NN,X,MM)
      DIMENSION A(NN),F(NN),T(20),D(20)
      LOGICAL EXTRA
      LOGICAL MFLAG,RFLAG
      DATA MMAX/10/
C
C  TABULAR INTERPOLATION USING SYMMETRICALLY PLACED ARGUMENT POINTS.
C
C  START.  FIND SUBSCRIPT IX OF X IN ARRAY A.
      IF( (NN.LT.2) .OR. (MM.LT.1) ) GO TO 20
      N=NN
      M=MIN0(MM,MMAX,N-1)
      MPLUS=M+1
      IX=0
      IY=N+1
      IF(A(1).GT.A(N)) GO TO 4
C     (SEARCH INCREASING ARGUMENTS.)
    1    MID=(IX+IY)/2
         IF(X.GE.A(MID)) GO TO 2
            IY=MID
            GO TO 3
C        (IF TRUE.)
    2       IX=MID
    3    IF(IY-IX.GT.1) GO TO 1
         GO TO 7
C     (SEARCH DECREASING ARGUMENTS.)
    4    MID=(IX+IY)/2
         IF(X.LE.A(MID)) GO TO 5
            IY=MID
            GO TO 6
C        (IF TRUE.)
    5       IX=MID
    6    IF(IY-IX.GT.1) GO TO 4
C
C  COPY REORDERED INTERPOLATION POINTS INTO (T(I),D(I)), SETTING
C  *EXTRA* TO TRUE IF M+2 POINTS TO BE USED.
    7 NPTS=M+2-MOD(M,2)
      IP=0
      L=0
      GO TO 9
    8    L=-L
         IF(L.GE.0) L=L+1
    9    ISUB=IX+L
         IF((1.LE.ISUB).AND.(ISUB.LE.N)) GO TO 10
C        (SKIP POINT.)
            NPTS=MPLUS
            GO TO 11
C        (INSERT POINT.)
   10       IP=IP+1
            T(IP)=A(ISUB)
            D(IP)=F(ISUB)
   11    IF(IP.LT.NPTS) GO TO 8
      EXTRA=NPTS.NE.MPLUS
C
C  REPLACE D BY THE LEADING DIAGONAL OF A DIVIDED-DIFFERENCE TABLE, SUP-
C  PLEMENTED BY AN EXTRA LINE IF *EXTRA* IS TRUE.
      DO 14 L=1,M
         IF(.NOT.EXTRA) GO TO 12
            ISUB=MPLUS-L
            D(M+2)=(D(M+2)-D(M))/(T(M+2)-T(ISUB))
   12    I=MPLUS
         DO 13 J=L,M
            ISUB=I-L
            D(I)=(D(I)-D(I-1))/(T(I)-T(ISUB))
            I=I-1
   13    CONTINUE
   14 CONTINUE
C
C  EVALUATE THE NEWTON INTERPOLATION FORMULA AT X, AVERAGING TWO VALUES
C  OF LAST DIFFERENCE IF *EXTRA* IS TRUE.
      SUM=D(MPLUS)
      IF(EXTRA) SUM=0.5*(SUM+D(M+2))
      J=M
      DO 15 L=1,M
         SUM=D(J)+(X-T(J))*SUM
         J=J-1
   15 CONTINUE
      DIVDIF=SUM
      RETURN
C
   20 CALL KERMTR('E105.1',LGFILE,MFLAG,RFLAG)
      DIVDIF=0
      IF(MFLAG) THEN
         IF(LGFILE.EQ.0) THEN
            IF(MM.LT.1) WRITE(*,101) MM
            IF(NN.LT.2) WRITE(*,102) NN
         ELSE
            IF(MM.LT.1) WRITE(LGFILE,101) MM
            IF(NN.LT.2) WRITE(LGFILE,102) NN
         ENDIF
      ENDIF
      IF(.NOT.RFLAG) STOP
      RETURN
  101 FORMAT( 7X, 'FUNCTION DIVDIF ... M =',I6,' IS LESS THAN 1')
  102 FORMAT( 7X, 'FUNCTION DIVDIF ... N =',I6,' IS LESS THAN 2')
      END


      SUBROUTINE KERSET(ERCODE,LGFILE,LIMITM,LIMITR)
      PARAMETER(KOUNTE  =  27)
      CHARACTER*6         ERCODE,   CODE(KOUNTE)
      LOGICAL             MFLAG,    RFLAG
      INTEGER             KNTM(KOUNTE),       KNTR(KOUNTE)
      DATA      LOGF      /  0  /
      DATA      CODE(1), KNTM(1), KNTR(1)  / 'C204.1', 255, 255 /
      DATA      CODE(2), KNTM(2), KNTR(2)  / 'C204.2', 255, 255 /
      DATA      CODE(3), KNTM(3), KNTR(3)  / 'C204.3', 255, 255 /
      DATA      CODE(4), KNTM(4), KNTR(4)  / 'C205.1', 255, 255 /
      DATA      CODE(5), KNTM(5), KNTR(5)  / 'C205.2', 255, 255 /
      DATA      CODE(6), KNTM(6), KNTR(6)  / 'C305.1', 255, 255 /
      DATA      CODE(7), KNTM(7), KNTR(7)  / 'C308.1', 255, 255 /
      DATA      CODE(8), KNTM(8), KNTR(8)  / 'C312.1', 255, 255 /
      DATA      CODE(9), KNTM(9), KNTR(9)  / 'C313.1', 255, 255 /
      DATA      CODE(10),KNTM(10),KNTR(10) / 'C336.1', 255, 255 /
      DATA      CODE(11),KNTM(11),KNTR(11) / 'C337.1', 255, 255 /
      DATA      CODE(12),KNTM(12),KNTR(12) / 'C341.1', 255, 255 /
      DATA      CODE(13),KNTM(13),KNTR(13) / 'D103.1', 255, 255 /
      DATA      CODE(14),KNTM(14),KNTR(14) / 'D106.1', 255, 255 /
      DATA      CODE(15),KNTM(15),KNTR(15) / 'D209.1', 255, 255 /
      DATA      CODE(16),KNTM(16),KNTR(16) / 'D509.1', 255, 255 /
      DATA      CODE(17),KNTM(17),KNTR(17) / 'E100.1', 255, 255 /
      DATA      CODE(18),KNTM(18),KNTR(18) / 'E104.1', 255, 255 /
      DATA      CODE(19),KNTM(19),KNTR(19) / 'E105.1', 255, 255 /
      DATA      CODE(20),KNTM(20),KNTR(20) / 'E208.1', 255, 255 /
      DATA      CODE(21),KNTM(21),KNTR(21) / 'E208.2', 255, 255 /
      DATA      CODE(22),KNTM(22),KNTR(22) / 'F010.1', 255,   0 /
      DATA      CODE(23),KNTM(23),KNTR(23) / 'F011.1', 255,   0 /
      DATA      CODE(24),KNTM(24),KNTR(24) / 'F012.1', 255,   0 /
      DATA      CODE(25),KNTM(25),KNTR(25) / 'F406.1', 255,   0 /
      DATA      CODE(26),KNTM(26),KNTR(26) / 'G100.1', 255, 255 /
      DATA      CODE(27),KNTM(27),KNTR(27) / 'G100.2', 255, 255 /
      LOGF  =  LGFILE
      L  =  0
      IF(ERCODE .NE. ' ')  THEN
         DO 10  L = 1, 6
            IF(ERCODE(1:L) .EQ. ERCODE)  GOTO 12
 10      CONTINUE
 12      CONTINUE
      ENDIF
      DO 14     I  =  1, KOUNTE
         IF(L .EQ. 0)  GOTO 13
         IF(CODE(I)(1:L) .NE. ERCODE(1:L))  GOTO 14
 13      IF(LIMITM.GE.0) KNTM(I)  =  LIMITM
         IF(LIMITR.GE.0) KNTR(I)  =  LIMITR
 14   CONTINUE
      RETURN
      ENTRY KERMTR(ERCODE,LOG,MFLAG,RFLAG)
      LOG  =  LOGF
      DO 20     I  =  1, KOUNTE
         IF(ERCODE .EQ. CODE(I))  GOTO 21
 20   CONTINUE
      WRITE(*,1000)  ERCODE
      STOP
      RETURN
 21   RFLAG  =  KNTR(I) .GE. 1
      IF(RFLAG  .AND.  (KNTR(I) .LT. 255))  KNTR(I)  =  KNTR(I) - 1
      MFLAG  =  KNTM(I) .GE. 1
      IF(MFLAG  .AND.  (KNTM(I) .LT. 255))  KNTM(I)  =  KNTM(I) - 1
      IF(.NOT. RFLAG)  THEN
         IF(LOGF .LT. 1)  THEN
            WRITE(*,1001)  CODE(I)
         ELSE
            WRITE(LOGF,1001)  CODE(I)
         ENDIF
      ENDIF
      IF(MFLAG .AND. RFLAG)  THEN
         IF(LOGF .LT. 1)  THEN
            WRITE(*,1002)  CODE(I)
         ELSE
            WRITE(LOGF,1002)  CODE(I)
         ENDIF
      ENDIF
      RETURN
 1000 FORMAT(' KERNLIB LIBRARY ERROR. ' /
     +     ' ERROR CODE ',A6,' NOT RECOGNIZED BY KERMTR',
     +     ' ERROR MONITOR. RUN ABORTED.')
 1001 FORMAT(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',
     +     'CONDITION ',A6)
 1002 FORMAT(/' ***** CERN LIBRARY ERROR CONDITION ',A6)
      END
