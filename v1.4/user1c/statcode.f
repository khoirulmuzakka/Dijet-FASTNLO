      SUBROUTINE STATCODE(TABPATH,SCENARIO,BORNN,NLON)
      
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      
      CHARACTER*(*) TABPATH,SCENARIO
      CHARACTER*255 BORNNAME,NLONAME
      CHARACTER*255 FILEBASE,FILENAME,LOFILE,NLOFILE
      CHARACTER*4 NO
      
      INTEGER BORNN,NLON,NTAB,ICOUNT,NCOUNT
      INTEGER ITAB,J,IBIN,NBINS,NMAX,ICYCLE,ISTAT,LENOCC
      INTEGER IORD,ISCL,ISUB,IRAP,IPT,IHIST,IPTMAX
      INTEGER NJMIN(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)
      INTEGER NJMAX(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)
      
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
      DOUBLE PRECISION NEFF(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)
      DOUBLE PRECISION MEANE(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)
ckr
ckr      DOUBLE PRECISION SIGMMW(NPTMAX,NRAPIDITY,NSCALEVAR,0:NORD)

      DOUBLE PRECISION MUR(NSCALEVAR),MUF(NSCALEVAR)
      DOUBLE PRECISION BWGT,NEVTS,VAL


 
c - HBOOK common
      INTEGER NWPAWC
      PARAMETER (NWPAWC=2500000)
      REAL HMEMOR(NWPAWC)
      COMMON /PAWC/ HMEMOR



c - Initialization
      BORNNAME = TABPATH(1:LENOCC(TABPATH))//"/stat/"//
     >     SCENARIO(1:LENOCC(SCENARIO))//"-hhc-born-"
      NLONAME  = TABPATH(1:LENOCC(TABPATH))//"/stat/"//
     >     SCENARIO(1:LENOCC(SCENARIO))//"-hhc-nlo-"
      IPTMAX = 44
      DO ISCL = 1,NSCALEVAR
         MUR(ISCL) = MURSCALE(ISCL)
         MUF(ISCL) = MUFSCALE(ISCL)
      ENDDO
      


c ==================================================================
c === Loop over all orders up to NLO ===============================
c ==================================================================
      DO IORD=0,MIN(2,NORD)

c - Loop initialization
         DO ISCL=1,NSCALEVAR
            DO IRAP=1,NRAPIDITY
               DO IPT=1,NPT(IRAP)
                  NJMIN(IPT,IRAP,ISCL,IORD) = -1
                  NJMAX(IPT,IRAP,ISCL,IORD) = -1
                  WGT(IPT,IRAP,ISCL,IORD)   = 0D0
                  WGT2(IPT,IRAP,ISCL,IORD)  = 0D0
                  WGTX(IPT,IRAP,ISCL,IORD)  = 0D0
                  WGTX2(IPT,IRAP,ISCL,IORD) = 0D0
                  MINE(IPT,IRAP,ISCL,IORD)  =  1D99
                  MAXE(IPT,IRAP,ISCL,IORD)  = -1D99
                  MEAN(IPT,IRAP,ISCL,IORD)  =  0D0
                  SIGMA(IPT,IRAP,ISCL,IORD) =  0D0
                  NEFF(IPT,IRAP,ISCL,IORD)  =  0D0
                  MEANE(IPT,IRAP,ISCL,IORD) =  0D0
ckr                  SIGMMW(IPT,IRAP,ISCL,IORD)=  0D0
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

      DO ITAB=0,NTAB
         WRITE(*,*) " ###############################################"
         IF (IORD.EQ.1) THEN
            WRITE(*,*)
     >           " ############  NEXT LO TABLE  ##################"
         ELSE
            WRITE(*,*)
     >           " ############  NEXT NLO TABLE ##################"
         ENDIF
         WRITE(*,*) " ###############################################"
         WRITE(NO,'(I4.4)'),ITAB
         FILENAME = FILEBASE(1:LENOCC(FILEBASE))//"2jet_"//NO//".tab"
         OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=ISTAT)
         IF (ISTAT.NE.0) THEN
            FILENAME = FILEBASE(1:LENOCC(FILEBASE))//"3jet_"//NO//".tab"
            OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=ISTAT)
            IF (ISTAT.NE.0) THEN
               WRITE(*,*)"STATERR: WARNING! Table file not found, "//
     >              "skipped! IOSTAT = ",ISTAT
ckr While using f90 DO-ENDDO ... one could also use the EXIT statement
ckr instead of GOTO. However, EXIT leaves the loop!
ckr To continue the loop use CYCLE instead! 
               GOTO 10
ckr            CYCLE
            ENDIF
         ENDIF
         CLOSE(2)
         WRITE(*,*)"Filename for order",IORD,":",FILENAME
         NCOUNT = NCOUNT + 1
         
c - Loop over scales
         DO ISCL=1,NSCALEVAR
            CALL FX9999CC(FILENAME,MUR(ISCL),MUF(ISCL),0,XSECT)
c - Take LO/NLO file with 1D9/1D8 events as weight 1 
            WTAB(NCOUNT) = NEVTS(IORD)*BWGT
            IF (ITAB.EQ.0 .AND. ISCL.EQ.1) THEN 
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
c - Write out error when the total cross section is negative for
c - the primary scale (taken to be 3 here)!
                     IF (VAL.LT.0.D0.AND.ISCL.EQ.3) THEN
                        WRITE(*,*)"STATERR: ERROR! "//
     >                       "Total cross section negative for file:",
     >                       NCOUNT,FILENAME
                        WRITE(*,*)"STATERR: IORD, ISCL, IRAP, IPT, VAL",
     >                       IORD,ISCL,IRAP,IPT,VAL
                        WRITE(*,*)"It is strongly suggested to "//
     >                       "exclude this file from the evaluation!"
                     ENDIF
                  ELSEIF (IORD.EQ.1) THEN
                     VAL = XSECT(J,1)
c - Write warning when the LO cross section is negative!
                     IF (VAL.LT.0.D0.AND.ISCL.EQ.3) THEN
                        WRITE(*,*)"STATERR: ERROR! "//
     >                       "LO cross section negative for file:",
     >                       NCOUNT,FILENAME
                        WRITE(*,*)"STATERR: IORD, ISCL, IRAP, IPT, VAL",
     >                       IORD,ISCL,IRAP,IPT,VAL
                        WRITE(*,*)"It is strongly suggested to "//
     >                       "exclude this file from the evaluation!"
                     ENDIF
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
ckr debug
                     if (iscl.eq.3) then
                        write(*,*)"MINIMUM: ipt,itab,mine,val",ipt,itab,
     >                       MINE(IPT,IRAP,ISCL,IORD),val
                     endif
                     MINE(IPT,IRAP,ISCL,IORD) = VAL
                     NJMIN(IPT,IRAP,ISCL,IORD) = ITAB
                  ENDIF
                  IF (VAL.GT.MAXE(IPT,IRAP,ISCL,IORD)) THEN
ckr debug
                     if (iscl.eq.3) then
                        write(*,*)"MAXIMUM: ipt,itab,maxe,val",ipt,itab,
     >                       MAXE(IPT,IRAP,ISCL,IORD),val
                     endif
                     MAXE(IPT,IRAP,ISCL,IORD) = VAL
                     NJMAX(IPT,IRAP,ISCL,IORD) = ITAB
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
 10      WRITE(*,*)"STATERR: NCOUNT, total weight:",NCOUNT,WTAB(NCOUNT)
      ENDDO
ckr      ENDDO

c - Extract mean values and standard deviations
      WRITE(*,*)"\n *************************************************"
      WRITE(*,*)"STATERR: Looping over scales for order:",IORD
      WRITE(*,*)"*************************************************"
      DO ISCL=1,NSCALEVAR
         WRITE(*,*)"STATERR: Next scale: ",ISCL,
     >        " ; weight: ",WTAB(NCOUNT)
         WRITE(*,*)"========================================"//
     >        "===================================="//
     >        "=========================="//
     >        "==========================="
         WRITE(*,*)"#IBIN #IMIN #IMAX         <s>           "//
     >        "s_min           s_max       ds/<s>/%"//
     >        " ds_min/<s>/% ds_max/<s>/%"//
     >        "     ds_min/ds    ds_max/ds"
         WRITE(*,*)"----------------------------------------"//
     >        "------------------------------------"//
     >        "--------------------------"//
     >        "---------------------------"
         J = 0
         DO IRAP=1,NRAPIDITY
            DO IPT=1,NPT(IRAP)
               J = J+1
ckr Markus: Entspricht mittl. Fehler des Mittelwerts sigma/sqrt(n)
ckr         fuer ungewichtete Daten
ckr            MEANL(J,K)  = WGTX(J,K) / WGT(J,K)
ckr            SIGMAL(J,K) = (WGTX2(J,K)-(WGTX(J,K)*WGTX(J,K)/WGT(J,K)))
ckr     >           /WGT(J,K)
ckr            SIGMAL(J,K) = SQRT(SIGMAL(J,K)/NCOUNT)
ckr Klaus: Entspricht mittl. Fehler des Mittelwerts sigma/sqrt(n-1)
ckr        auch fuer gew. Daten
ckr        Ein Freiheitsgrad weniger fuer Mittelwertsberechnung!
ckr               WRITE(*,*)"STATERR: DEBUG: IPT,IRAP,ISCL,IORD",IPT,IRAP,ISCL,IORD
ckr Sample average
               MEAN(IPT,IRAP,ISCL,IORD) =
     >              WGTX(IPT,IRAP,ISCL,IORD) /
     >              WGT(IPT,IRAP,ISCL,IORD)
ckr Effective number of events
               NEFF(IPT,IRAP,ISCL,IORD) =
     >              WGT(IPT,IRAP,ISCL,IORD)*WGT(IPT,IRAP,ISCL,IORD) /
     >              WGT2(IPT,IRAP,ISCL,IORD)
ckr The right formula for the spread of x section values of single jobs
               SIGMA(IPT,IRAP,ISCL,IORD) =
     >              (WGTX2(IPT,IRAP,ISCL,IORD) / 
     >              WGT(IPT,IRAP,ISCL,IORD) -
     >              MEAN(IPT,IRAP,ISCL,IORD) *
     >              MEAN(IPT,IRAP,ISCL,IORD))
               SIGMA(IPT,IRAP,ISCL,IORD) =
     >              SQRT(SIGMA(IPT,IRAP,ISCL,IORD))
ckr Markus formula: Mittl. Fehler des Mittelwerts
ckr               SIGMMW(IPT,IRAP,ISCL,IORD) =
ckr     >              (WGTX2(IPT,IRAP,ISCL,IORD) - 
ckr     >              (WGTX(IPT,IRAP,ISCL,IORD)*WGTX(IPT,IRAP,ISCL,IORD) /
ckr     >              WGT(IPT,IRAP,ISCL,IORD))) / WGT(IPT,IRAP,ISCL,IORD)
ckr               SIGMMW(IPT,IRAP,ISCL,IORD) =
ckr     >              SQRT(SIGMMW(IPT,IRAP,ISCL,IORD) /
ckr     >              WGT(IPT,IRAP,ISCL,IORD)) 
ckr My formula:
               MEANE(IPT,IRAP,ISCL,IORD) =
     >              SIGMA(IPT,IRAP,ISCL,IORD)/
     >              SQRT(NEFF(IPT,IRAP,ISCL,IORD)-1D0)
c - Write warning when the minimal or maximal cross sections have not
c - the same sign as the average!
Comment:                IF (MINE(IPT,IRAP,ISCL,IORD) *
Comment:      >              MEAN(IPT,IRAP,ISCL,IORD).LT.0D0) THEN
Comment:                   WRITE(*,*)"STATERR: WARNING! "//
Comment:      >                 "Minimal cross section has sign different from"//
Comment:      >                 " average:",
Comment:      >                 MINE(IPT,IRAP,ISCL,IORD),MEAN(IPT,IRAP,ISCL,IORD)
Comment:                   WRITE(*,*)"STATERR: IORD, ISCL, IRAP, IPT",
Comment:      >                 IORD,ISCL,IRAP,IPT
Comment:                ENDIF
Comment:                IF (MAXE(IPT,IRAP,ISCL,IORD) *
Comment:      >              MEAN(IPT,IRAP,ISCL,IORD).LT.0D0) THEN
Comment:                   WRITE(*,*)"STATERR: WARNING! "//
Comment:      >                 "Maximal cross section has sign different from"//
Comment:      >                 " average:",
Comment:      >                 MAXE(IPT,IRAP,ISCL,IORD),MEAN(IPT,IRAP,ISCL,IORD)
Comment:                   WRITE(*,*)"STATERR: IORD, ISCL, IRAP, IPT",
Comment:      >                 IORD,ISCL,IRAP,IPT
Comment:                ENDIF
               WRITE(*,900) J,
     >              NJMIN(IPT,IRAP,ISCL,IORD),
     >              NJMAX(IPT,IRAP,ISCL,IORD),
     >              MEAN(IPT,IRAP,ISCL,IORD),
     >              MINE(IPT,IRAP,ISCL,IORD),
     >              MAXE(IPT,IRAP,ISCL,IORD),
     >              100D0*MEANE(IPT,IRAP,ISCL,IORD) /
     >              DABS(MEAN(IPT,IRAP,ISCL,IORD)),
     >              100D0*(MINE(IPT,IRAP,ISCL,IORD) -
     >              MEAN(IPT,IRAP,ISCL,IORD)) /
     >              DABS(MEAN(IPT,IRAP,ISCL,IORD)),
     >              100D0*(MAXE(IPT,IRAP,ISCL,IORD) -
     >              MEAN(IPT,IRAP,ISCL,IORD)) / 
     >              DABS(MEAN(IPT,IRAP,ISCL,IORD)),
     >              (MINE(IPT,IRAP,ISCL,IORD)-
     >              MEAN(IPT,IRAP,ISCL,IORD))/SIGMA(IPT,IRAP,ISCL,IORD),
     >              (MAXE(IPT,IRAP,ISCL,IORD)-
     >              MEAN(IPT,IRAP,ISCL,IORD))/SIGMA(IPT,IRAP,ISCL,IORD)
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
ckr rel. to MEAN and then use central result. If not, just use MEANE.
               IF (HI(IHIST,IPT).GT.0.) THEN
                  SERR(IPT,IRAP) = REAL(MEANE(IPT,IRAP,ISCL,IORD) /
     >                 MEAN(IPT,IRAP,ISCL,IORD))*HI(IHIST,IPT)
               ELSE
                  SERR(IPT,IRAP) = REAL(MEANE(IPT,IRAP,ISCL,IORD))
               ENDIF
               CALL HFILL(IHIST+3,PT(IPT,IRAP),0.,
     >              MAX(0.,REAL(100D0*MEANE(IPT,IRAP,ISCL,IORD) / 
     >              MEAN(IPT,IRAP,ISCL,IORD))))
               CALL HFILL(IHIST+4,PT(IPT,IRAP),0.,
     >              MAX(0.,REAL(50D0*(MAXE(IPT,IRAP,ISCL,IORD) - 
     >              MINE(IPT,IRAP,ISCL,IORD)) / 
     >              MEAN(IPT,IRAP,ISCL,IORD))))
ckr               write(*,*)"SECOND: iord,iscl,irap,ipt"
ckr     >              ,iord,iscl,irap,ipt
               DO ICOUNT=1,NCOUNT
                  IHIST = IORD*1000000 + ISCL*100000 +
     >                 ISUB*10000 + IRAP*100
ckr                  write(*,*)"SECOND: icount,cbin,mean,sigma",icount,
ckr     >                 CBIN(IPT,IRAP,ISCL,IORD,ICOUNT),
ckr     >                 MEAN(IPT,IRAP,ISCL,IORD),
ckr     >                 SIGMA(IPT,IRAP,ISCL,IORD)
                  CALL HFILL(IHIST + 10,
     >                 REAL(
     >                 ( CBIN(IPT,IRAP,ISCL,IORD,ICOUNT) -
     >                 MEAN(IPT,IRAP,ISCL,IORD) ) / 
     >                 SIGMA(IPT,IRAP,ISCL,IORD))
     >                 ,0.,REAL(WTAB(ICOUNT)))
                  CALL HFILL(IHIST + 11,
     >                 REAL(2D0*(CBIN(IPT,IRAP,ISCL,IORD,ICOUNT) -
     >                 MEAN(IPT,IRAP,ISCL,IORD)) /
     >                 (MAXE(IPT,IRAP,ISCL,IORD) - 
     >                 MINE(IPT,IRAP,ISCL,IORD))),
     >                 0.,REAL(WTAB(ICOUNT)))
ckr IHIST limit before next rapidity bin: IHIST+2*IPT+11 < IHIST+100
ckr => Maximal IPT = IPTMAX < 45
                  IF (IPT.LE.IPTMAX) THEN
                     CALL HFILL(IHIST + 2*IPT + 10,
     >                    REAL(
     >                    ( CBIN(IPT,IRAP,ISCL,IORD,ICOUNT) -
     >                    MEAN(IPT,IRAP,ISCL,IORD) ) / 
     >                    SIGMA(IPT,IRAP,ISCL,IORD))
     >                    ,0.,REAL(WTAB(ICOUNT)))
                     CALL HFILL(IHIST + 2*IPT + 11,
     >                    REAL(2D0*(CBIN(IPT,IRAP,ISCL,IORD,ICOUNT) -
     >                    MEAN(IPT,IRAP,ISCL,IORD)) /
     >                    (MAXE(IPT,IRAP,ISCL,IORD) - 
     >                    MINE(IPT,IRAP,ISCL,IORD))),
     >                    0.,REAL(WTAB(ICOUNT)))
                  ENDIF
               ENDDO
            ENDDO
            CALL HPAKE(IHIST,SERR(1,IRAP))
         ENDDO
      ENDDO
      
      ENDDO

 900  FORMAT (3I6,3E16.5,5(F10.3,3X))
      
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
