      PROGRAM STATERR
* -------------------------------------------------------------------
* K. Rabbertz 16.05.2008 Largely reworked version
* M. Wobisch  12/27/2005 
* 
* STATERR - program to compute statistical errors from a large set
*           of LO/NLO tables
*
* Please specify:
*     - scenario name (e.g. fnl0002)
*     - no. of last LO table to be included (max. 999)
*     - no. of last NLO table to be included (max. 999)
*       (non-existing tables in the range from 0000 - end number are
*        skipped automatically)
* Optionally give:
*     - name of desired PDF set (default is cteq65.LHgrid)
*       (The files are expected to be in $LHAPDF/../PDFsets or ./../PDFsets)
*     - name of the output file with the statistics histograms
*       (default is scenario_stat.hbk) 
*
* -------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*255 SCENARIO,BORNNAME,NLONAME,HISTFILE
      CHARACTER*255 PDFSET,PDFPATH,LHAPDF,ASMODE
      CHARACTER*4 CHBORN,CHNLO
      INTEGER BORNN,NLON,LENOCC

c --- Parse command line
ckr 30.01.2008: Some more checks on input arguments
      WRITE(*,*)"\n ##############################################"
      WRITE(*,*)"# STATERR"
      WRITE(*,*)"##############################################"
      WRITE(*,*)"# Program to compute statistical uncertainties"
      WRITE(*,*)"# from a large number of raw fastNLO tables"
      WRITE(*,*)"##############################################"
      WRITE(*,*)"#"
      IF (IARGC().LT.1) THEN
         WRITE(*,*)
     >        "STATERR: ERROR! No scenario name given, "//
     >        "aborting!"
         WRITE(*,*)"      For an explanation of command line "//
     >        "arguments type:"
         WRITE(*,*)"      ./staterr -h"
         STOP
      ELSE
         CALL GETARG(1,SCENARIO)
         IF (SCENARIO(1:LENOCC(SCENARIO)).EQ."-h") THEN
            WRITE(*,*)" "
            WRITE(*,*)"Usage: ./staterr (arguments)"
            WRITE(*,*)"  scenario name"
            WRITE(*,*)"  last LO table number to use"
            WRITE(*,*)"  last NLO table number to use"
            WRITE(*,*)"  HBOOK output file, def. = (scenario)_stat.hbk"

            WRITE(*,*)"  PDF set, def. = cteq65.LHgrid"
            WRITE(*,*)"  PDF path, def. = $(LHAPDF)/"//
     >           "../share/lhapdf/PDFsets"
            WRITE(*,*)"  alpha_s calc., def. from PDF set"
            WRITE(*,*)" "
            STOP
         ENDIF
         WRITE(*,*)"STATERR: Evaluating scenario: ",
     >        SCENARIO(1:LENOCC(SCENARIO))
      ENDIF
      IF (IARGC().LT.2) THEN
         WRITE(*,*)
     >        "STATERR: ERROR! Last number of LO tables not given, "//
     >        "aborting!"
         STOP
      ELSE
         CALL GETARG(2,CHBORN)
         READ(CHBORN,'(I4)'),BORNN
         WRITE(*,*)"STATERR: Last LO table number: ",BORNN
      ENDIF
      IF (IARGC().LT.3) THEN
         WRITE(*,*)
     >        "STATERR: ERROR! Last number of NLO tables not given, "/
     >        /"aborting!"
         STOP
      ELSE
         CALL GETARG(3,CHNLO)
         READ(CHNLO,'(I4)'),NLON
         WRITE(*,*)"STATERR: Last NLO table number: ",NLON
      ENDIF
      IF (IARGC().LT.4) THEN
         HISTFILE = SCENARIO(1:LENOCC(SCENARIO))//"_stat.hbk"
         WRITE(*,*)
     >        "STATERR: No histo filename given, "//
     >        "using "//HISTFILE(1:LENOCC(HISTFILE))//"!"
      ELSE
         CALL GETARG(4,HISTFILE)
         WRITE(*,*)"STATERR: Using histo filename: ",
     >        HISTFILE(1:LENOCC(HISTFILE))
      ENDIF
      IF (IARGC().LT.5) THEN
         PDFSET = "cteq65.LHgrid"
         WRITE(*,*)
     >        "STATERR: WARNING! No PDF set given, "//
     >        "taking "//PDFSET(1:LENOCC(PDFSET))//" instead!"
      ELSE
         CALL GETARG(5,PDFSET)
         WRITE(*,*)"STATERR: Using PDF set: ",
     >        PDFSET(1:LENOCC(PDFSET))
      ENDIF
      IF (IARGC().LT.6) THEN
         PDFPATH = "/../share/lhapdf/PDFsets"
         WRITE(*,*)
     >        "STATERR: No PDF path given, "//
     >        "assuming: $(LHAPDF)"//PDFPATH
c - Initialize path to LHAPDF libs
         CALL GETENV("LHAPDF",LHAPDF)
         IF (LENOCC(LHAPDF).EQ.0) THEN
            WRITE(*,*)"\nSTATERR: ERROR! $LHAPDF not set, aborting!"
            STOP
         ENDIF
         PDFPATH = LHAPDF(1:LENOCC(LHAPDF))//PDFPATH(1:LENOCC(PDFPATH))
      ELSE
         CALL GETARG(6,PDFPATH)
      ENDIF
      WRITE(*,*)"STATERR: Looking for LHAPDF PDF sets in path: "
     >     //PDFPATH(1:LENOCC(PDFPATH))
      PDFSET = PDFPATH(1:LENOCC(PDFPATH))//"/"//PDFSET
      WRITE(*,*)"STATERR: Taking PDF set "
     >     //PDFSET(1:LENOCC(PDFSET))
      IF (IARGC().GT.6) THEN
         WRITE(*,*)"STATERR: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF

c - Initialize path to LHAPDF libs
      CALL GETENV('LHAPDF',LHAPDF)
      IF (LENOCC(LHAPDF).EQ.0) THEN
         LHAPDF = "."
      ENDIF
      WRITE(*,*)"STATERR: Looking for LHAPDF sets in directory "
     &     //LHAPDF(1:LENOCC(LHAPDF))
      PDFSET = LHAPDF(1:LENOCC(LHAPDF))//"/../PDFsets/"//PDFSET
      WRITE(*,*)"STATERR: Taking PDF set "
     &     //PDFSET(1:LENOCC(PDFSET))

c - Initialize LHAPDF
      CALL INITPDFSET(PDFSET(1:LENOCC(PDFSET)))

c - Initialize one member, 0=best fit member
      CALL INITPDF(0)

c - Call statistical error-code for scenario
      BORNNAME = SCENARIO(1:LENOCC(SCENARIO))//"-hhc-born-2jet_"
      NLONAME  = SCENARIO(1:LENOCC(SCENARIO))//"-hhc-nlo-2jet_"
      CALL STATCODE(BORNN,BORNNAME,NLON,NLONAME,HISTFILE)

      END

C -----------------------------------------------------------------

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
      CALL HLIMIT(NWPAWC)
      CALL HROPEN(11,'fastNLO',HISTFILE,'N',1024,ISTAT)
      IF (ISTAT.NE.0) THEN
         WRITE(*,*)"STATERR: ERROR! Could not open histofile: ",ISTAT
      ENDIF

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
         FILENAME = BORNNAME(1:LENOCC(BORNNAME))//NO//".stc"
         WRITE(*,*)"LO filename:",FILENAME
         OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=ISTAT)
         IF (ISTAT.NE.0) THEN
            WRITE(*,*)"STATERR: WARNING! Table file not found, "//
     >           "skipped! IOSTAT = ",ISTAT
ckr While using f90 DO-ENDDO ... use also the EXIT statement instead of GOTO
ckr            GOTO 10
            EXIT
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
               CALL HBOOK1(1000001,"LO sigma/mu in % vs. bin number",
     >              NBINS,0.5,REAL(NBINS+0.5),0)
               CALL HBOOK1(1000002,"LO dmax/2/mu in % vs. bin number",
     >              NBINS,0.5,REAL(NBINS+0.5),0)
               DO IBIN=0,NBINS
                  WRITE(NO,'(I4)'),IBIN
                  CALL HBOOK1(1001000+IBIN,
     >                 "Distr. of LO x sections norm. to "//
     >                 "mu,sigma for bin no. "//NO,
     >                 21,-5.5,5.5,0)
                  CALL HBOOK1(1002000+IBIN,
     >                 "Distr. of LO x sections norm. to "//
     >                 "mu,dmax half for bin no. "//NO,
     >                 31,-1.5,1.5,0)
               ENDDO
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
      K = KHIST
      DO J=1,NBINS
         CALL HFILL(1000001,REAL(J),0.,
     >        MAX(0.,REAL(100D0*SIGMAL(J,K)/MEANL(J,K))))
         CALL HFILL(1000002,REAL(J),0.,
     >        MAX(0.,REAL(50D0*(MAXL(J,K)-MINL(J,K))/MEANL(J,K))))
         DO I=1,NBORN
CKR            CALL HFILL(1000000+J,REAL(CBIN(I,J,K)),0.,REAL(WTAB(I))
c            write(*,*)"hfill cbin,mean,sig,x",
c     >           CBIN(I,J,K),MEANL(J,K),SIGMAL(J,K),
c     >           (CBIN(I,J,K)-MEANL(J,K))/SIGMAL(J,K)
            CALL HFILL(1001000,
     >           REAL((CBIN(I,J,K)-MEANL(J,K))/SIGMAL(J,K)),
     >           0.,REAL(WTAB(I)))
            CALL HFILL(1001000+J,
     >           REAL((CBIN(I,J,K)-MEANL(J,K))/SIGMAL(J,K)),
     >           0.,REAL(WTAB(I)))
            CALL HFILL(1002000,
     >           REAL(2D0*(CBIN(I,J,K)-MEANL(J,K))/(MAXL(J,K)-MINL(J,K))
     >           ),0.,REAL(WTAB(I)))
            CALL HFILL(1002000+J,
     >           REAL(2D0*(CBIN(I,J,K)-MEANL(J,K))/(MAXL(J,K)-MINL(J,K))
     >           ),0.,REAL(WTAB(I)))
         ENDDO
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
         FILENAME = NLONAME(1:LENOCC(NLONAME))//NO//".stc"
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
               CALL HBOOK1(2000001,"NLO sigma/mu in % vs. bin number",
     >              NBINS,0.5,REAL(NBINS+0.5),0)
               CALL HBOOK1(2000002,"NLO dmax/2/mu in % vs. bin number",
     >              NBINS,0.5,REAL(NBINS+0.5),0)
               DO IBIN=0,NBINS
                  WRITE(NO,'(I4)'),IBIN
                  CALL HBOOK1(2001000+IBIN,
     >                 "Distr. of NLO x sections norm. to "//
     >                 "mu,sigma for bin no. "//NO,
     >                 21,-5.5,5.5,0)
                  CALL HBOOK1(2002000+IBIN,
     >                 "Distr. of NLO x sections norm. to "//
     >                 "mu,dmax half for bin no. "//NO,
     >                 31,-1.5,1.5,0)
               ENDDO
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

c Fill histos for scale no. 3, mu = 1.0      
c Fill histos for scale KHIST, normally no. 3, mu = 1.0      
      K = KHIST
      DO J=1,NBINS
         CALL HFILL(2000001,REAL(J),0.,
     >        MAX(0.,REAL(100D0*SIGMAN(J,K)/MEANN(J,K))))
         CALL HFILL(2000002,REAL(J),0.,
     >        MAX(0.,REAL(50D0*(MAXN(J,K)-MINN(J,K))/MEANN(J,K))))
         DO I=1,NNLO
            CALL HFILL(2001000,
     >           REAL((CBIN(I,J,K)-MEANN(J,K))/SIGMAN(J,K)),
     >           0.,REAL(WTAB(I)))
            CALL HFILL(2001000+J,
     >           REAL((CBIN(I,J,K)-MEANN(J,K))/SIGMAN(J,K)),
     >           0.,REAL(WTAB(I)))
            CALL HFILL(2002000,
     >           REAL(2D0*(CBIN(I,J,K)-MEANN(J,K))/(MAXN(J,K)-MINN(J,K))
     >           ),0.,REAL(WTAB(I)))
            CALL HFILL(2002000+J,
     >           REAL(2D0*(CBIN(I,J,K)-MEANN(J,K))/(MAXN(J,K)-MINN(J,K))
     >           ),0.,REAL(WTAB(I)))
         ENDDO
      ENDDO
      
 999  CONTINUE
      
      WRITE(*,*)"\n **********************************************"
      WRITE(*,*)"STATERR: Job finished, storing histos in file: ",
     >     HISTFILE(1:LENOCC(HISTFILE))
      WRITE(*,*)"STATERR: Scale used for histo filling: ",KHIST
      WRITE(*,*)"**********************************************"
      CALL HROUT (0,ICYCLE,' ')
      CALL HREND ('fastNLO')
      
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
