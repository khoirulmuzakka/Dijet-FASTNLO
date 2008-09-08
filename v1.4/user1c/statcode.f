      SUBROUTINE STATCODE(BORNN,BORNNAME,NLON,NLONAME,HISTFILE)

      IMPLICIT NONE
c - Attention!!! - this must be declared consistent with its 
c                  definition in the commonblock!!!!!
c      DOUBLE PRECISION XSECT(900,3) 
      INCLUDE 'fnx9999.inc'
      
      CHARACTER*(*) BORNNAME,NLONAME,HISTFILE
      CHARACTER*255 FILENAME, LOFILE, NLOFILE
      CHARACTER*4 NO
      
      INTEGER BORNN,NLON,NBORN,NNLO
      INTEGER I,J,K,KHIST,IBIN,NBINS,NMAX,ICYCLE,ISTAT,LENOCC

      INTEGER IORD,ISUB,J1,J2,IHIST
      REAL PT(NPTMAX),SERR(NPTMAX,NRAPIDITY),HI
      
      DOUBLE PRECISION WGT(3*NBINTOTMAX,4),WGT2(3*NBINTOTMAX,4)
      DOUBLE PRECISION WGTX(3*NBINTOTMAX,4),WGTX2(3*NBINTOTMAX,4)
      DOUBLE PRECISION MEANL(3*NBINTOTMAX,4),SIGMAL(3*NBINTOTMAX,4)
      DOUBLE PRECISION MEANN(3*NBINTOTMAX,4),SIGMAN(3*NBINTOTMAX,4)
      DOUBLE PRECISION MINL(3*NBINTOTMAX,4),MAXL(3*NBINTOTMAX,4)
      DOUBLE PRECISION MINN(3*NBINTOTMAX,4),MAXN(3*NBINTOTMAX,4)
      DOUBLE PRECISION NEVTS,VAL
      DOUBLE PRECISION CBIN(1000,3*NBINTOTMAX,4),WTAB(1000)
      INTEGER NJMIN(3*NBINTOTMAX),NJMAX(3*NBINTOTMAX)
      DOUBLE PRECISION MU(4) ! Should be read from table ...
      DATA KHIST/3/
      DATA MU/0.25D0, 0.5D0, 1.0D0, 2.0D0/
 
c - HBOOK common
      INTEGER NWPAWC
      PARAMETER (NWPAWC=2500000)
      REAL HMEMOR(NWPAWC)
      COMMON /PAWC/ HMEMOR

c - Initialization
      DO J=1,3*NBINTOTMAX       ! BINS
         NJMIN(J) = -1
         NJMAX(J) = -1
         DO K=1,4               ! SCALES
            MEANL(J,K)  = 0D0
            SIGMAL(J,K) = 0D0
            MEANN(J,K)  = 0D0
            SIGMAN(J,K) = 0D0
            MAXL(J,K)   = -1D99
            MAXN(J,K)   = -1D99
            MINL(J,K)   = 1D99
            MINN(J,K)   = 1D99
         ENDDO
      ENDDO
      
c - Open & book
Comment:       CALL HLIMIT(NWPAWC)
Comment:       CALL HROPEN(11,'fastNLO',HISTFILE,'N',1024,ISTAT)
Comment:       IF (ISTAT.NE.0) THEN
Comment:          WRITE(*,*)"STATERR: ERROR! Could not open histofile: ",ISTAT
Comment:       ENDIF

c ==================================================================
c === For LO ======================================================
c ==================================================================
      WRITE(*,*)"\n *************************************************"
      WRITE(*,*)"STATERR: Compute statistical errors for LO tables"
      WRITE(*,*)"*************************************************"
      
      DO J=1,3*NBINTOTMAX
         DO K=1,4
            WGT(J,K)   = 0D0
            WGT2(J,K)  = 0D0
            WGTX(J,K)  = 0D0
            WGTX2(J,K) = 0D0
         ENDDO
      ENDDO

c - Loop over files
      NBORN = 0
      DO I=0,BORNN
         WRITE(*,*) " ###############################################"
         WRITE(*,*) " ############  NEXT LO TABLE  ##################"
         WRITE(*,*) " ###############################################"
         WRITE(NO,'(I4.4)'),I
         FILENAME = BORNNAME(1:LENOCC(BORNNAME))//NO//".tab"
         WRITE(*,*)"LO filename:",FILENAME
         OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=ISTAT)
         IF (ISTAT.NE.0) THEN
            WRITE(*,*)"STATERR: WARNING! Table file not found, "//
     >           "skipped! IOSTAT = ",ISTAT
ckr While using f90 DO-ENDDO ... use also the EXIT statement instead of GOTO
ckr            GOTO 10
            EXIT
         ELSE
            CLOSE(2)
         ENDIF
         NBORN = NBORN + 1

c - Loop over scales
         DO K=1,4
            CALL FX9999CC(FILENAME,MU(K),MU(K),0,XSECT)
c - Take LO file with 1D9 events as weight 1 
            IF (K.EQ.1) WTAB(NBORN) = NEVTS(1)/1D9
ckr Histo booking, for first LO table and first scale only
            IF (I.EQ.0 .AND. K.EQ.1) THEN 
               NBINS = NMAX(1)
               WRITE(*,*)"STATERR: This file has ",NBINS," bins"
               IF (NBINS.GT.NBINTOTMAX) THEN
                  WRITE(*,*)"STATERR: ERROR! NBINTOTMAX too small:",
     >                 NBINTOTMAX
                  STOP
               ENDIF
c - LO only at single scale
Comment:                CALL HBOOK1(1000001,"LO sigma/mu in % vs. bin number",
Comment:      >              NBINS,0.5,REAL(NBINS+0.5),0)
Comment:                CALL HBOOK1(1000002,"LO dmax/2/mu in % vs. bin number",
Comment:      >              NBINS,0.5,REAL(NBINS+0.5),0)
Comment:                DO IBIN=0,NBINS
Comment:                   WRITE(NO,'(I4)'),IBIN
Comment:                   CALL HBOOK1(1001000+IBIN,
Comment:      >                 "Distr. of LO x sections norm. to "//
Comment:      >                 "mu,sigma for bin no. "//NO,
Comment:      >                 21,-5.5,5.5,0)
Comment:                   CALL HBOOK1(1002000+IBIN,
Comment:      >                 "Distr. of LO x sections norm. to "//
Comment:      >                 "mu,dmax half for bin no. "//NO,
Comment:      >                 31,-1.5,1.5,0)
Comment:                ENDDO
            ENDIF
            
            DO J=1,NBINS
               VAL = XSECT(J,1)
               CBIN(NBORN,J,K) = VAL
               WGT(J,K)   = WGT(J,K)  + WTAB(NBORN)
               WGT2(J,K)  = WGT2(J,K) + WTAB(NBORN)*WTAB(NBORN)
               WGTX(J,K)  = WGTX(J,K) + VAL * WTAB(NBORN)
               WGTX2(J,K) = WGTX2(J,K)+ VAL*VAL * WTAB(NBORN)
               IF (VAL.LT.MINL(J,K)) MINL(J,K) = VAL
               IF (VAL.GT.MAXL(J,K)) MAXL(J,K) = VAL
            ENDDO
         ENDDO
         WRITE(*,*)"STATERR: Total weight:",WTAB(NBORN)
ckr 10   ENDDO
      ENDDO

c - Extract mean values and standard deviations
      WRITE(*,*)"\n *************************************************"
      WRITE(*,*)"STATERR: Looping over scales, LO ..."
      WRITE(*,*)"*************************************************"
      DO K=1,4
         WRITE(*,*)"STATERR: Next scale: ",K," ; weight: ",WTAB(NBORN)
         DO J=1,NBINS
ckr Markus: Entspricht sigma mit 1/sqrt(n) fuer ungewichtete Daten
c            MEANL(J,K)  = WGTX(J,K) / WGT(J,K)
c            SIGMAL(J,K) = (WGTX2(J,K)-(WGTX(J,K)*WGTX(J,K)/WGT(J,K)))
c     >           /WGT(J,K)
c            SIGMAL(J,K) = SQRT(SIGMAL(J,K)/NBORN)
ckr Klaus: Entspricht sigma mit 1/sqrt(n-1) fuer ungew. Daten
ckr        Ein Freiheitsgrad fuer Mittelwertsberechnung!
            MEANL(J,K)  = WGTX(J,K) / WGT(J,K)
            SIGMAL(J,K) = (WGTX2(J,K)/WGT(J,K)-MEANL(J,K)*MEANL(J,K)) /
     >           (1D0 - WGT2(J,K)/WGT(J,K)/WGT(J,K))
            SIGMAL(J,K) = SQRT(SIGMAL(J,K)/NBORN)
c            WRITE(*,*)"WGT,WGT2",WGT(J,K),WGT2(J,K)
c            WRITE(*,*)"MEANL,N",MEANL(J,K),MEANN(J,K)
c            WRITE(*,*)"SIGMAL,N",SIGMAL(J,K),SIGMAN(J,K)
            WRITE(*,900) J,MEANL(J,K),
     >           100D0*SIGMAL(J,K)/MEANL(J,K),
     >           100D0*(MINL(J,K)-MEANL(J,K))/MEANL(J,K),
     >           100D0*(MAXL(J,K)-MEANL(J,K))/MEANL(J,K)
         ENDDO
      ENDDO 
      
c Fill histos for scale KHIST, normally no. 3, mu = 1.0      
      IORD = 1
      K = KHIST
      ISUB = 0

ckr      DO J=1,NBINS
      J = 0
      DO J1=1,NRAPIDITY         ! RAPIDITY BINS
         IHIST = IORD*1000000 + K*100000 + ISUB*10000 + J1*100
         DO J2=1,NPT(J1)        ! PT BINS
            J = J + 1
            PT(J) = REAL(PTBIN(J1,J2))
ckr If histogram with central result (IHIST+0) filled, derive stat. unc.
ckr rel. to MEAN and then use central result. If not, just use SIGMA. 
            IF (HI(IHIST,J2).GT.0.) THEN
               SERR(J2,J1) = REAL(SIGMAL(J,K)/MEANL(J,K))*HI(IHIST,J2)
            ELSE
               SERR(J2,J1) = REAL(SIGMAL(J,K))
            ENDIF
            CALL HFILL(IHIST+3,PT(J),0.,
     >           MAX(0.,REAL(100D0*SIGMAL(J,K)/MEANL(J,K))))
            CALL HFILL(IHIST+4,PT(J),0.,
     >           MAX(0.,REAL(50D0*(MAXL(J,K)-MINL(J,K))/MEANL(J,K))))
Comment:         DO I=1,NBORN
CKR            CALL HFILL(1000000+J,REAL(CBIN(I,J,K)),0.,REAL(WTAB(I))
c            write(*,*)"hfill cbin,mean,sig,x",
c     >           CBIN(I,J,K),MEANL(J,K),SIGMAL(J,K),
c     >           (CBIN(I,J,K)-MEANL(J,K))/SIGMAL(J,K)
Comment:             CALL HFILL(1001000,
Comment:      >           REAL((CBIN(I,J,K)-MEANL(J,K))/SIGMAL(J,K)),
Comment:      >           0.,REAL(WTAB(I)))
Comment:             CALL HFILL(1001000+J,
Comment:      >           REAL((CBIN(I,J,K)-MEANL(J,K))/SIGMAL(J,K)),
Comment:      >           0.,REAL(WTAB(I)))
Comment:             CALL HFILL(1002000,
Comment:      >           REAL(2D0*(CBIN(I,J,K)-MEANL(J,K))/(MAXL(J,K)-MINL(J,K))
Comment:      >           ),0.,REAL(WTAB(I)))
Comment:             CALL HFILL(1002000+J,
Comment:      >           REAL(2D0*(CBIN(I,J,K)-MEANL(J,K))/(MAXL(J,K)-MINL(J,K))
Comment:      >           ),0.,REAL(WTAB(I)))
Comment:         ENDDO
         ENDDO
         CALL HPAKE(IHIST,SERR(1,J1))
      ENDDO

ckr      GOTO 999

c ==================================================================
c === for NLO ======================================================
c ==================================================================
      WRITE(*,*)"\n **************************************************"
      WRITE(*,*)"STATERR: Compute statistical errors for NLO tables"
      WRITE(*,*)"**************************************************"

      DO J=1,3*NBINTOTMAX
         DO K=1,4
            WGT(J,K)   = 0D0
            WGT2(J,K)  = 0D0
            WGTX(J,K)  = 0D0
            WGTX2(J,K) = 0D0
         ENDDO
      ENDDO

c - Loop over files
      NNLO = 0
      DO I=0,NLON
         WRITE(*,*) " ###############################################"
         WRITE(*,*) " ############  NEXT NLO TABLE  #################"
         WRITE(*,*) " ###############################################"
         WRITE(NO,'(I4.4)'),I
         FILENAME = NLONAME(1:LENOCC(NLONAME))//NO//".tab"
         WRITE(*,*)"NLO filename:",FILENAME
         OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=ISTAT)
         IF (ISTAT.NE.0) THEN
            WRITE(*,*)"STATERR: WARNING! Table file not found, "//
     >           "skipped! IOSTAT = ",ISTAT
ckr While using f90 DO-ENDDO ... use also the EXIT statement instead of GOTO
ckr            GOTO 20
            EXIT
         ENDIF
         NNLO = NNLO + 1

c - Loop over scales
         DO K=1,4
            CALL FX9999CC(FILENAME,MU(K),MU(K),0,XSECT)
c - Take NLO file with 1D8 events as weight 1 
            IF (K.EQ.1) WTAB(NNLO) = NEVTS(2)/1D8
ckr Histo booking, for first NLO table and first scale only
            IF (I.EQ.0 .AND. K.EQ.1) THEN 
               NBINS = NMAX(1)
               WRITE(*,*)"STATERR: This file has ",NBINS," bins"
               IF (NBINS.GT.NBINTOTMAX) THEN
                  WRITE(*,*)"STATERR: ERROR! NBINTOTMAX too small:",
     >                 NBINTOTMAX
                  STOP
               ENDIF
c - NLO for now only at single scale
Comment:                CALL HBOOK1(2000001,"NLO sigma/mu in % vs. bin number",
Comment:      >              NBINS,0.5,REAL(NBINS+0.5),0)
Comment:                CALL HBOOK1(2000002,"NLO dmax/2/mu in % vs. bin number",
Comment:      >              NBINS,0.5,REAL(NBINS+0.5),0)
Comment:                DO IBIN=0,NBINS
Comment:                   WRITE(NO,'(I4)'),IBIN
Comment:                   CALL HBOOK1(2001000+IBIN,
Comment:      >                 "Distr. of NLO x sections norm. to "//
Comment:      >                 "mu,sigma for bin no. "//NO,
Comment:      >                 21,-5.5,5.5,0)
Comment:                   CALL HBOOK1(2002000+IBIN,
Comment:      >                 "Distr. of NLO x sections norm. to "//
Comment:      >                 "mu,dmax half for bin no. "//NO,
Comment:      >                 31,-1.5,1.5,0)
Comment:                ENDDO
            ENDIF

            DO J=1,NBINS
               VAL = XSECT(J,1)+XSECT(J,2)
               CBIN(NNLO,J,K) = VAL
c               IF (VAL.LT.0D0) WRITE(*,*)"J,K,NNLO,VAL",
c     >              J,K,NNLO,VAL
               WGT(J,K)   = WGT(J,K)  + WTAB(NNLO)
               WGT2(J,K)  = WGT2(J,K) + WTAB(NNLO)*WTAB(NNLO)
               WGTX(J,K)  = WGTX(J,K) + VAL * WTAB(NNLO)
               WGTX2(J,K) = WGTX2(J,K)+ VAL*VAL * WTAB(NNLO)
               IF (VAL.LT.MINN(J,K)) THEN
                  MINN(J,K) = VAL
                  NJMIN(J)  = I
               ENDIF
               IF (VAL.GT.MAXN(J,K)) THEN
                  MAXN(J,K) = VAL
                  NJMAX(J)  = I
               ENDIF
            ENDDO
         ENDDO
         WRITE(*,*)"STATERR: Total weight:",WTAB(NNLO)
ckr 20   ENDDO
      ENDDO

c - Extract mean values and standard deviations
      WRITE(*,*)"\n *************************************************"
      WRITE(*,*)"STATERR: Looping over scales, NLO ..."
      WRITE(*,*)"*************************************************"
      DO K=1,4
         WRITE(*,*)"STATERR: Next scale: ",K," ; weight: ",WTAB(NNLO)
         DO J=1,NBINS
            MEANN(J,K)  = WGTX(J,K) / WGT(J,K)
            SIGMAN(J,K) = (WGTX2(J,K)/WGT(J,K)-MEANN(J,K)*MEANN(J,K)) /
     >           (1D0 - WGT2(J,K)/WGT(J,K)/WGT(J,K))
            SIGMAN(J,K) = SQRT(SIGMAN(J,K)/NNLO)
            WRITE(*,901) J,MEANN(J,K),
     >           100D0*SIGMAN(J,K)/MEANN(J,K),
     >           100D0*(MINN(J,K)-MEANN(J,K))/MEANN(J,K),
     >           100D0*(MAXN(J,K)-MEANN(J,K))/MEANN(J,K),
     >           NJMIN(J),NJMAX(J)
         ENDDO
      ENDDO
      
 900  FORMAT (I4,E16.5,"  in %:",3F10.3)
 901  FORMAT (I4,E16.5,"  in %:",3F10.3,2X,2I5)

c Fill histos for scale KHIST, normally no. 3, mu = 1.0      
      IORD = 0
      K = KHIST
      ISUB = 0

ckr      DO J=1,NBINS
      J = 0
      DO J1=1,NRAPIDITY         ! RAPIDITY BINS
         IHIST = IORD*1000000 + K*100000 + ISUB*10000 + J1*100
         DO J2=1,NPT(J1)        ! PT BINS
            J = J + 1
            PT(J) = REAL(PTBIN(J1,J2))
ckr If histogram with central result (IHIST+0) filled, derive stat. unc.
ckr rel. to MEAN and then use central result. If not, just use SIGMA. 
            IF (HI(IHIST,J2).GT.0.) THEN
               SERR(J2,J1) = REAL(SIGMAN(J,K)/MEANN(J,K))*HI(IHIST,J2)
            ELSE
               SERR(J2,J1) = REAL(SIGMAN(J,K))
            ENDIF
            CALL HFILL(IHIST+3,PT(J),0.,
     >           MAX(0.,REAL(100D0*SIGMAN(J,K)/MEANN(J,K))))
            CALL HFILL(IHIST+4,PT(J),0.,
     >           MAX(0.,REAL(50D0*(MAXN(J,K)-MINN(J,K))/MEANN(J,K))))
Comment:          DO I=1,NNLO
Comment:             CALL HFILL(2001000,
Comment:      >           REAL((CBIN(I,J,K)-MEANN(J,K))/SIGMAN(J,K)),
Comment:      >           0.,REAL(WTAB(I)))
Comment:             CALL HFILL(2001000+J,
Comment:      >           REAL((CBIN(I,J,K)-MEANN(J,K))/SIGMAN(J,K)),
Comment:      >           0.,REAL(WTAB(I)))
Comment:             CALL HFILL(2002000,
Comment:      >           REAL(2D0*(CBIN(I,J,K)-MEANN(J,K))/(MAXN(J,K)-MINN(J,K))
Comment:      >           ),0.,REAL(WTAB(I)))
Comment:             CALL HFILL(2002000+J,
Comment:      >           REAL(2D0*(CBIN(I,J,K)-MEANN(J,K))/(MAXN(J,K)-MINN(J,K))
Comment:      >           ),0.,REAL(WTAB(I)))
Comment:         ENDDO
         ENDDO
         CALL HPAKE(IHIST,SERR(1,J1))
      ENDDO
      
 999  CONTINUE
      
      WRITE(*,*)"\n **********************************************"
      WRITE(*,*)"STATERR: Job finished, storing histos in file: ",
     >     HISTFILE(1:LENOCC(HISTFILE))
      WRITE(*,*)"STATERR: Scale used for histo filling: ",KHIST
      WRITE(*,*)"**********************************************"
Comment:       CALL HROUT (0,ICYCLE,' ')
Comment:       CALL HREND ('fastNLO')
      
      END
      


C-------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION NEVTS(I)
      IMPLICIT NONE
      INTEGER I
      INCLUDE 'fnx9999.inc'

      NEVTS=DBLE(NEVT(I))

      RETURN
      END
C-------------------------------------------------------------------
      INTEGER FUNCTION NMAX(I)
      IMPLICIT NONE
      INTEGER I,J1,J2, NBINS
      INCLUDE 'fnx9999.inc'

      NBINS = 0
      DO J1=1,NRAPIDITY          ! RAPIDITY BINS
         DO J2=1,NPT(J1)          ! PT BINS
            NBINS = NBINS + 1
         ENDDO
      ENDDO
      NMAX= NBINS

      RETURN
      END
