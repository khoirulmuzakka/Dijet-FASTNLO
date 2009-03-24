      SUBROUTINE STATCODE(BORNN,BORNNAME,NLON,NLONAME)
      
      IMPLICIT NONE
c - Attention!!! - this must be declared consistent with its 
c                  definition in the commonblock!!!!!
c      DOUBLE PRECISION XSECT(900,3) 
      INCLUDE 'fnx9999.inc'
      
      CHARACTER*(*) BORNNAME,NLONAME
      CHARACTER*255 FILEBASE,FILENAME,LOFILE,NLOFILE
      CHARACTER*4 NO
      
ckr      INTEGER BORNN,NLON,NBORN,NNLO,NTAB,NCOUNT
      INTEGER BORNN,NLON,NTAB,ICOUNT,NCOUNT
      INTEGER I,J,IBIN,NBINS,NMAX,ICYCLE,ISTAT,LENOCC
      INTEGER IORD,ISCL,ISUB,IRAP,IPT,IHIST
      INTEGER NJMIN(NPTMAX,NRAPIDITY,NSCALEVAR)
      INTEGER NJMAX(NPTMAX,NRAPIDITY,NSCALEVAR)

      REAL PT(NPTMAX,NRAPIDITY),SERR(NPTMAX,NRAPIDITY),HI
      
      DOUBLE PRECISION WGT(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)
      DOUBLE PRECISION WGT2(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)
      DOUBLE PRECISION WGTX(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)
      DOUBLE PRECISION WGTX2(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)
      DOUBLE PRECISION CBIN(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD,
     >     NBINTOTMAX)
      DOUBLE PRECISION WTAB(3*NBINTOTMAX)
      DOUBLE PRECISION MINE(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)
      DOUBLE PRECISION MAXE(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)
      DOUBLE PRECISION MEAN(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)
      DOUBLE PRECISION SIGMA(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)

      DOUBLE PRECISION MUR(NSCALEVAR),MUF(NSCALEVAR)
      DOUBLE PRECISION BWGT,NEVTS,VAL


 
c - HBOOK common
      INTEGER NWPAWC
      PARAMETER (NWPAWC=2500000)
      REAL HMEMOR(NWPAWC)
      COMMON /PAWC/ HMEMOR



c - Initialization
      DO ISCL = 1,NSCALEVAR
         MUR(ISCL) = MURSCALE(ISCL)
         MUF(ISCL) = MUFSCALE(ISCL)
      ENDDO
      


c ==================================================================
c === Loop over all orders up to NLO ===============================
c ==================================================================
      DO IORD=0,MIN(2,NORD)
ckr      DO IORD=1,MIN(2,NORD)

c - Loop initialization
         DO ISCL=1,NSCALEVAR
            DO IRAP=1,NRAPIDITY
               DO IPT=1,NPT(IRAP)
                  NJMIN(IPT,IRAP,ISCL)      = -1
                  NJMAX(IPT,IRAP,ISCL)      = -1
                  WGT(IPT,IRAP,ISCL,IORD)   = 0D0
                  WGT2(IPT,IRAP,ISCL,IORD)  = 0D0
                  WGTX(IPT,IRAP,ISCL,IORD)  = 0D0
                  WGTX2(IPT,IRAP,ISCL,IORD) = 0D0
                  MINE(IPT,IRAP,ISCL,IORD)  =  1D99
                  MAXE(IPT,IRAP,ISCL,IORD)  = -1D99
                  MEAN(IPT,IRAP,ISCL,IORD)  =  0D0
                  SIGMA(IPT,IRAP,ISCL,IORD) =  0D0
               ENDDO
            ENDDO
         ENDDO

      WRITE(*,*)"\n *************************************************"
      WRITE(*,*)"STATERR: Compute statistical errors for order:", IORD
      WRITE(*,*)"*************************************************"

c - Loop over files
      IF (IORD.EQ.0) THEN
c - Total x section
         NTAB = NLON 
         BWGT = 1D0/1D8
         FILEBASE = NLONAME(1:LENOCC(NLONAME))
      ELSEIF (IORD.EQ.1) THEN
c - LO
         NTAB = BORNN
         BWGT = 1D0/1D9
         FILEBASE = BORNNAME(1:LENOCC(BORNNAME))
      ELSEIF (IORD.EQ.2) THEN
c - NLO
         NTAB = NLON
         BWGT = 1D0/1D8
         FILEBASE = NLONAME(1:LENOCC(NLONAME))
      ELSE
         WRITE(*,*)"STATERR: ERROR! Illegal order for stat. calc:",
     >        IORD
         STOP
      ENDIF
      NCOUNT = 0

      DO I=0,NTAB
         WRITE(*,*) " ###############################################"
         IF (IORD.EQ.1) THEN
            WRITE(*,*)
     >           " ############  NEXT LO TABLE  ##################"
         ELSE
            WRITE(*,*)
     >           " ############  NEXT NLO TABLE ##################"
         ENDIF
         WRITE(*,*) " ###############################################"
         WRITE(NO,'(I4.4)'),I
         FILENAME = FILEBASE(1:LENOCC(FILEBASE))//NO//".tab"
         WRITE(*,*)"Filename for order",IORD,":",FILENAME
         OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=ISTAT)
         IF (ISTAT.NE.0) THEN
            WRITE(*,*)"STATERR: WARNING! Table file not found, "//
     >           "skipped! IOSTAT = ",ISTAT
ckr While using f90 DO-ENDDO ... one could also use the EXIT statement
ckr instead of GOTO. However, EXIT leaves the loop!
ckr To continue the loop use CYCLE instead! 
            GOTO 10
ckr            CYCLE
         ELSE
            CLOSE(2)
         ENDIF
         NCOUNT = NCOUNT + 1
         
c - Loop over scales
         DO ISCL=1,NSCALEVAR
            CALL FX9999CC(FILENAME,MUR(ISCL),MUF(ISCL),0,XSECT)
c - Take LO file with 1D9 events as weight 1 
            WTAB(NCOUNT) = NEVTS(IORD)*BWGT
            IF (I.EQ.0 .AND. ISCL.EQ.1) THEN 
               NBINS = NMAX(1)
               WRITE(*,*)"STATERR: This file has ",NBINS," bins"
               IF (NBINS.GT.NBINTOTMAX) THEN
                  WRITE(*,*)"STATERR: ERROR! NBINTOTMAX too small:",
     >                 NBINTOTMAX
                  STOP
               ENDIF
            ENDIF

            J = 0
            DO IRAP=1,NRAPIDITY
               DO IPT=1,NPT(IRAP)
                  J = J+1
                  IF (IORD.EQ.0) THEN
                     VAL = XSECT(J,1)+XSECT(J,2)
                  ELSEIF (IORD.EQ.1) THEN
                     VAL = XSECT(J,1)
                  ELSEIF (IORD.EQ.2) THEN
                     VAL = XSECT(J,2)
                  ELSE
                     WRITE(*,*)
     >                    "STATERR: ERROR! Order too large:",
     >                    IORD
                     STOP
                  ENDIF
                  CBIN(IPT,IRAP,ISCL,IORD,NCOUNT) = VAL
                  WGT(IPT,IRAP,ISCL,IORD) =
     >                 WGT(IPT,IRAP,ISCL,IORD) +
     >                 WTAB(NCOUNT)
                  WGT2(IPT,IRAP,ISCL,IORD) =
     >                 WGT2(IPT,IRAP,ISCL,IORD) +
     >                 WTAB(NCOUNT)*WTAB(NCOUNT)
                  WGTX(IPT,IRAP,ISCL,IORD) =
     >                 WGTX(IPT,IRAP,ISCL,IORD)
     >                 + VAL * WTAB(NCOUNT)
                  WGTX2(IPT,IRAP,ISCL,IORD) =
     >                 WGTX2(IPT,IRAP,ISCL,IORD) +
     >                 VAL*VAL * WTAB(NCOUNT)
                  IF (VAL.LT.MINE(IPT,IRAP,ISCL,IORD)) THEN
                     MINE(IPT,IRAP,ISCL,IORD) = VAL
                     NJMIN(IPT,IRAP,ISCL) = I
                  ENDIF
                  IF (VAL.GT.MAXE(IPT,IRAP,ISCL,IORD)) THEN
                     MAXE(IPT,IRAP,ISCL,IORD) = VAL
                     NJMAX(IPT,IRAP,ISCL) = I
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         WRITE(*,*)"STATERR: NCOUNT, total weight:",NCOUNT,WTAB(NCOUNT)
 10   ENDDO
ckr      ENDDO

c - Extract mean values and standard deviations
      WRITE(*,*)"\n *************************************************"
      WRITE(*,*)"STATERR: Looping over scales for order:",IORD
      WRITE(*,*)"*************************************************"
      DO ISCL=1,NSCALEVAR
         WRITE(*,*)"STATERR: Next scale: ",ISCL,
     >        " ; weight: ",WTAB(NCOUNT)
         J = 0
         DO IRAP=1,NRAPIDITY
            DO IPT=1,NPT(IRAP)
               J = J+1
ckr Markus: Entspricht sigma mit 1/sqrt(n) fuer ungewichtete Daten
c            MEANL(J,K)  = WGTX(J,K) / WGT(J,K)
c            SIGMAL(J,K) = (WGTX2(J,K)-(WGTX(J,K)*WGTX(J,K)/WGT(J,K)))
c     >           /WGT(J,K)
c            SIGMAL(J,K) = SQRT(SIGMAL(J,K)/NCOUNT)
ckr Klaus: Entspricht sigma mit 1/sqrt(n-1) fuer ungew. Daten
ckr        Ein Freiheitsgrad fuer Mittelwertsberechnung!
ckr               write(*,*)"ipt,irap,iscl,iord",ipt,irap,iscl,iord
               MEAN(IPT,IRAP,ISCL,IORD) =
     >              WGTX(IPT,IRAP,ISCL,IORD) /
     >              WGT(IPT,IRAP,ISCL,IORD)
ckr               write(*,*)"wgtx,wgt,mean",
ckr     >              WGTX(IPT,IRAP,ISCL,IORD),
ckr     >              WGT(IPT,IRAP,ISCL,IORD),
ckr     >              MEAN(IPT,IRAP,ISCL,IORD)
               SIGMA(IPT,IRAP,ISCL,IORD) =
     >              (WGTX2(IPT,IRAP,ISCL,IORD) / 
     >              WGT(IPT,IRAP,ISCL,IORD) -
     >              MEAN(IPT,IRAP,ISCL,IORD) *
     >              MEAN(IPT,IRAP,ISCL,IORD)) /
     >              (1D0 - WGT2(IPT,IRAP,ISCL,IORD) /
     >              WGT(IPT,IRAP,ISCL,IORD)/WGT(IPT,IRAP,ISCL,IORD))
               SIGMA(IPT,IRAP,ISCL,IORD) =
     >              SQRT(SIGMA(IPT,IRAP,ISCL,IORD)/NCOUNT)
               WRITE(*,901) J,MEAN(IPT,IRAP,ISCL,IORD),
     >              100D0*SIGMA(IPT,IRAP,ISCL,IORD) /
     >              MEAN(IPT,IRAP,ISCL,IORD),
     >              100D0*(MINE(IPT,IRAP,ISCL,IORD) -
     >              MEAN(IPT,IRAP,ISCL,IORD)) /
     >              MEAN(IPT,IRAP,ISCL,IORD),
     >              100D0*(MAXE(IPT,IRAP,ISCL,IORD) -
     >              MEAN(IPT,IRAP,ISCL,IORD)) / 
     >              MEAN(IPT,IRAP,ISCL,IORD),
     >              NJMIN(IPT,IRAP,ISCL),NJMAX(IPT,IRAP,ISCL)
            ENDDO
         ENDDO
      ENDDO 
      
c - Fill histos
c - Only for sum of subprocesses
      ISUB = 0
      DO ISCL=1,NSCALEVAR
         DO IRAP=1,NRAPIDITY
            DO IPT=1,NPT(IRAP)
               IHIST = IORD*1000000 + ISCL*100000 +
     >              ISUB*10000 + IRAP*100
               PT(IPT,IRAP) = REAL(PTBIN(IRAP,IPT))
ckr If histogram with central result (IHIST+0) filled, derive stat. unc.
ckr rel. to MEAN and then use central result. If not, just use SIGMA. 
               IF (HI(IHIST,IPT).GT.0.) THEN
                  SERR(IPT,IRAP) = REAL(SIGMA(IPT,IRAP,ISCL,IORD) /
     >                 MEAN(IPT,IRAP,ISCL,IORD))*HI(IHIST,IPT)
               ELSE
                  SERR(IPT,IRAP) = REAL(SIGMA(IPT,IRAP,ISCL,IORD))
               ENDIF
Comment:                write(*,*)"ihist+3,pt,input",
Comment:      >              IHIST+3,PT(IPT,IRAP),
Comment:      >              MAX(0.,REAL(100D0*SIGMA(IPT,IRAP,ISCL,IORD) / 
Comment:      >              MEAN(IPT,IRAP,ISCL,IORD)))
Comment:                write(*,*)"ipt,irap,iscl,iord,mean,sigma",
Comment:      >              ipt,irap,iscl,iord,
Comment:      >              MEAN(IPT,IRAP,ISCL,IORD),
Comment:      >              SIGMA(IPT,IRAP,ISCL,IORD)
               CALL HFILL(IHIST+3,PT(IPT,IRAP),0.,
     >              MAX(0.,REAL(100D0*SIGMA(IPT,IRAP,ISCL,IORD) / 
     >              MEAN(IPT,IRAP,ISCL,IORD))))
               CALL HFILL(IHIST+4,PT(IPT,IRAP),0.,
     >              MAX(0.,REAL(50D0*(MAXE(IPT,IRAP,ISCL,IORD) - 
     >              MINE(IPT,IRAP,ISCL,IORD)) / 
     >              MEAN(IPT,IRAP,ISCL,IORD))))
               DO ICOUNT=1,NCOUNT
ckr            CALL HFILL(1000000+J,REAL(CBIN(I,J,K)),0.,REAL(WTAB(I))
                  IHIST = IORD*1000000 + ISCL*100000 +
     >                 ISUB*10000 + IRAP*100
                  CALL HFILL(IHIST + 10,
     >                 REAL( (CBIN(IPT,IRAP,ISCL,IORD,ICOUNT) -
     >                 MEAN(IPT,IRAP,ISCL,IORD)) / 
     >                 SIGMA(IPT,IRAP,ISCL,IORD)),
     >                 0.,REAL(WTAB(I)))
                  write(*,*)"ipt,irap,iscl,iord",
     >                 ipt,irap,iscl,iord
                  write(*,*)"icount,cbin",
     >                 icount,CBIN(IPT,IRAP,ISCL,IORD,ICOUNT)
                  CALL HFILL(IHIST + 2*IPT + 10,
     >                 REAL( (CBIN(IPT,IRAP,ISCL,IORD,ICOUNT) -
     >                 MEAN(IPT,IRAP,ISCL,IORD)) / 
     >                 SIGMA(IPT,IRAP,ISCL,IORD)),
     >                 0.,REAL(WTAB(I)))
                  CALL HFILL(IHIST + 11,
     >                 REAL(2D0*(CBIN(IPT,IRAP,ISCL,IORD,ICOUNT) -
     >                 MEAN(IPT,IRAP,ISCL,IORD)) /
     >                 (MAXE(IPT,IRAP,ISCL,IORD) - 
     >                 MINE(IPT,IRAP,ISCL,IORD))),
     >                 0.,REAL(WTAB(I)))
                  CALL HFILL(IHIST + 2*IPT + 11,
     >                 REAL(2D0*(CBIN(IPT,IRAP,ISCL,IORD,ICOUNT) -
     >                 MEAN(IPT,IRAP,ISCL,IORD)) /
     >                 (MAXE(IPT,IRAP,ISCL,IORD) - 
     >                 MINE(IPT,IRAP,ISCL,IORD))),
     >                 0.,REAL(WTAB(I)))
               ENDDO
            ENDDO
            CALL HPAKE(IHIST,SERR(1,IRAP))
         ENDDO
      ENDDO
      
      ENDDO

 900  FORMAT (I4,E16.5,"  in %:",3F10.3)
 901  FORMAT (I4,E16.5,"  in %:",3F10.3,2X,2I5)
      
      WRITE(*,*)"\n **********************************************"
      WRITE(*,*)"STATERR: Job finished, "//
     >     "all scales used for histo filling!"
      WRITE(*,*)"**********************************************"
      
      END
      


C-------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION NEVTS(I)
      IMPLICIT NONE
      INTEGER I
      INCLUDE 'fnx9999.inc'
      
      IF (I.EQ.0) THEN
         NEVTS=DBLE(NEVT(2))
      ELSE
         NEVTS=DBLE(NEVT(I))
      ENDIF
      
      RETURN
      END
C-------------------------------------------------------------------
      INTEGER FUNCTION NMAX(I)
      IMPLICIT NONE
      INTEGER I,IRAP,IPT,NBINS
      INCLUDE 'fnx9999.inc'
      
      NBINS = 0
      DO IRAP=1,NRAPIDITY
         DO IPT=1,NPT(IRAP)
            NBINS = NBINS + 1
         ENDDO
      ENDDO
      NMAX= NBINS
      
      RETURN
      END
