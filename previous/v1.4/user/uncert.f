      SUBROUTINE UNCERT(IPHASE,IMODE,IWEIGHT,IVAR,LRAT,LNRM)

      IMPLICIT NONE
      INTEGER IPHASE,IMODE,IWEIGHT,IVAR
      LOGICAL LRAT,LNRM
      INCLUDE "fnx9999.inc"
      INCLUDE "uncert.inc"
      DOUBLE PRECISION DIFF,SUMM,RATIO,SUMM2,DIFF2
      DOUBLE PRECISION DIFFPROC,SUMMPROC,DIFFORD,SUMMORD
      DOUBLE PRECISION TWGT(3),TEVTS,WEIGHT,NEFF
      INTEGER IBIN,IRAP,IPT,IORD,ISUB,NBIN,NRAP,NRAP2
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST/.TRUE./

*---First run
c - Implement sample averages with different weights
c - Not active/well tested yet
c - Required in case of e.g. tables with changing nos. of events
c - Define weights to be one for 1D9 events for LO and 1D8 for NLO/total
c - If all weights are identical this has to cancel out completely
      IF (FIRST) THEN
         FIRST = .FALSE.
c - This is the order counting in fnx9999 common ...
         IF (NORD.GT.3) THEN
            WRITE(*,*)"UNCERT: ERROR! Orders higher than ",NORD,
     >           " are not supported, stopped!"
            STOP
         ENDIF
         DO IORD=0,NORD
            IF (IORD.EQ.0) THEN
               TWGT(NORD+1) = 1D0/1D8
            ELSEIF (IORD.EQ.1) THEN
               TWGT(IORD)   = 1D0/1D9
            ELSEIF (IORD.EQ.2) THEN
               TWGT(IORD)   = 1D0/1D8
            ELSE
               TWGT(IORD)   = 1D0/1D8
            ENDIF
         ENDDO
      ENDIF



*---Initialization
      IF (IPHASE.EQ.1) THEN
         IBIN = 0
         NRAP = NRAPIDITY
         IF (LRAT.OR.LNRM) NRAP=2*NRAPIDITY
         DO IRAP=1,NRAP
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               IJMIN(IBIN) = -1
               IJMAX(IBIN) = -1
               DO IORD=1,NORD+1
                  DO ISUB=1,NSUBPROC+1
                     WRES(IBIN,ISUB,IORD)   = 0D0
                     NTN(IBIN,ISUB,IORD)    = 0
                     WT(IBIN,ISUB,IORD)     = 0D0
                     WT2(IBIN,ISUB,IORD)    = 0D0
                     WTX(IBIN,ISUB,IORD)    = 0D0
                     WTX2(IBIN,ISUB,IORD)   = 0D0
                     WTXMIN(IBIN,ISUB,IORD) = +1D99
                     WTXMAX(IBIN,ISUB,IORD) = -1D99
                     WTDXL2(IBIN,ISUB,IORD) = 0D0
                     WTDXU2(IBIN,ISUB,IORD) = 0D0
                     WTDXLM(IBIN,ISUB,IORD) = 0D0
                     WTDXUM(IBIN,ISUB,IORD) = 0D0
                     WTDXMN(IBIN,ISUB,IORD) = 0D0
                     WTDXUL(IBIN,ISUB,IORD) = 0D0
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         NBIN = IBIN
         IF (LRAT.OR.LNRM) NBIN = NBIN/2
         RETURN
      ENDIF



*---Loop
      IF (IPHASE.EQ.2) THEN
         IBIN = 0
         NRAP2 = NRAPIDITY
         IF (LNRM) NRAP2 = NRAP
         DO IRAP=1,NRAP2
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               SUMMORD = 0D0
               DIFFORD = 0D0
               DO IORD=1,NORD
                  IF (IORD.LT.NORD) THEN
                     TEVTS = DBLE(NEVT(IORD))
                  ELSE
                     TEVTS = DBLE(NEVT(2))
                  ENDIF
                  WEIGHT = 1D0
                  IF (IWEIGHT.EQ.1) WEIGHT = TWGT(IORD)*TEVTS
                  SUMMPROC = 0D0
                  DIFFPROC = 0D0
                  DO ISUB=1,NSUBPROC
                     SUMM = RESULT(IBIN,ISUB,IORD)
                     SUMMPROC = SUMMPROC + SUMM
                     DIFF = SUMM - MYRESN(IBIN,ISUB,IORD)
                     DIFFPROC = DIFFPROC + DIFF
cdebug
Comment:                      WRITE(*,*)"ibin,isub,iord,result,myresn",
Comment:      >                    IBIN,ISUB,IORD,
Comment:      >                    result(ibin,isub,iord),
Comment:      >                    myresn(ibin,isub,iord)
Comment:                      WRITE(*,*)"DEF : summ,diff",SUMM,DIFF
cdebug
                     IF (LNRM) THEN
                        SUMM = MYRES(IBIN,ISUB,IORD)
                        DIFF = MYRES(IBIN,ISUB,IORD) -
     >                       MYRESN(IBIN,ISUB,IORD)
cdebug
Comment:                         WRITE(*,*)"ibin,isub,iord,myres,myresn",
Comment:      >                       IBIN,ISUB,IORD,
Comment:      >                       MYRES(IBIN,ISUB,IORD),
Comment:      >                       MYRESN(IBIN,ISUB,IORD)
Comment:                         WRITE(*,*)"LNRM: summ,diff",SUMM,DIFF
cdebug
                     ENDIF
                     CALL SUMMUP(IVAR,IBIN,ISUB,IORD,
     >                    SUMM,DIFF,WEIGHT,
     >                    LRAT,IPT,NBIN,LNRM)
                  ENDDO
                  SUMM = SUMMPROC
                  SUMMORD = SUMMORD + SUMMPROC
                  DIFF = SUMM - MYRESN(IBIN,NSUBPROC+1,IORD)
                  DIFFORD = DIFFORD + DIFF
                  IF (LNRM) THEN
                     SUMMPROC = MYRES(IBIN,ISUB,IORD)
                     DIFFPROC = MYRES(IBIN,ISUB,IORD) -
     >                    MYRESN(IBIN,ISUB,IORD)
                  ENDIF
                  CALL SUMMUP(IVAR,IBIN,NSUBPROC+1,IORD,
     >                 SUMMPROC,DIFFPROC,WEIGHT,
     >                 LRAT,IPT,NBIN,LNRM)
               ENDDO
               IF (LNRM) THEN
                  SUMMORD = MYRES(IBIN,ISUB,IORD)
                  DIFFORD = MYRES(IBIN,ISUB,IORD) -
     >                 MYRESN(IBIN,ISUB,IORD)
               ENDIF
               CALL SUMMUP(IVAR,IBIN,NSUBPROC+1,NORD+1,
     >              SUMMORD,DIFFORD,WEIGHT,
     >              LRAT,IPT,NBIN,LNRM)
               DO ISUB=1,NSUBPROC
                  SUMMORD = 0D0
                  DIFFORD = 0D0
                  DO IORD=1,NORD
                     SUMM = RESULT(IBIN,ISUB,IORD)
                     SUMMORD = SUMMORD + SUMM
                     DIFF = SUMM - MYRESN(IBIN,ISUB,IORD)
                     DIFFORD = DIFFORD + DIFF
                  ENDDO
                  IF (LNRM) THEN
                     SUMMORD = MYRES(IBIN,ISUB,IORD)
                     DIFFORD = MYRES(IBIN,ISUB,IORD) -
     >                    MYRESN(IBIN,ISUB,IORD)
                  ENDIF
                  CALL SUMMUP(IVAR,IBIN,ISUB,NORD+1,
     >                 SUMMORD,DIFFORD,WEIGHT,
     >                 LRAT,IPT,NBIN,LNRM)
               ENDDO
            ENDDO
         ENDDO
         RETURN
      ENDIF



*---Final evaluation
      IF (IPHASE.EQ.3) THEN
         IBIN = 0
         DO IRAP=1,NRAP
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               DO IORD=1,NORD+1
                  DO ISUB=1,NSUBPROC+1
ckr Eigen vector method (as needed for CTEQ or MSTW PDF uncertainty)
                     IF (IMODE.EQ.1) THEN
                        WTDXU2(IBIN,ISUB,IORD) =  SQRT(WTDXU2(IBIN
     >                       ,ISUB,IORD))
                        WTDXL2(IBIN,ISUB,IORD) = -SQRT(WTDXL2(IBIN
     >                       ,ISUB,IORD))
cdebug
Comment:                         IF (DABS(MYRES(IBIN,ISUB,IORD)).GT.1D-99) THEN
Comment:                            WTDXU2(IBIN,ISUB,IORD) =
Comment:      >                          WTDXU2(IBIN,ISUB,IORD) /
Comment:      >                          DABS(MYRES(IBIN,ISUB,IORD))
Comment:                            WTDXL2(IBIN,ISUB,IORD) =
Comment:      >                          WTDXL2(IBIN,ISUB,IORD) /
Comment:      >                          DABS(MYRES(IBIN,ISUB,IORD))
cdebug
                        IF (DABS(MYRESN(IBIN,ISUB,IORD)).GT.1D-99) THEN
                           WTDXU2(IBIN,ISUB,IORD) =
     >                          WTDXU2(IBIN,ISUB,IORD) /
     >                          DABS(MYRESN(IBIN,ISUB,IORD))
                           WTDXL2(IBIN,ISUB,IORD) =
     >                          WTDXL2(IBIN,ISUB,IORD) /
     >                          DABS(MYRESN(IBIN,ISUB,IORD))
                        ELSE
                           WTDXU2(IBIN,ISUB,IORD) = -1D0
                           WTDXL2(IBIN,ISUB,IORD) = -1D0
                        ENDIF
ckr Sampling method (as needed for statistical uncertainties or NNPDF)
                     ELSEIF (IMODE.EQ.2) THEN
                        WTDXU2(IBIN,ISUB,IORD) = -1D0
                        WTDXL2(IBIN,ISUB,IORD) = -1D0
                        WTDXMN(IBIN,ISUB,IORD) = -1D0
                        WTDXUL(IBIN,ISUB,IORD) = -1D0
                        IF (DABS(WT(IBIN,ISUB,IORD)).GT.1D-99) THEN
                           IF (IWEIGHT.EQ.0) THEN
                              NEFF = WT(IBIN,ISUB,IORD)
                           ELSE
                              NEFF = WT(IBIN,ISUB,IORD)*WT(IBIN,ISUB
     >                             ,IORD)/WT2(IBIN,ISUB,IORD)
                           ENDIF
ckr Attention: This overwrites the precalculated result by averages
ckr            from the sampling method as this should be considered the
ckr            central result according to NNPDF.
                           MYRES(IBIN,ISUB,IORD) =
     >                          WTX(IBIN,ISUB,IORD)/WT(IBIN,ISUB,IORD)
                           WTDXU2(IBIN,ISUB,IORD) =
     >                          (WTX2(IBIN,ISUB,IORD)/WT(IBIN,ISUB,IORD)
     >                          -MYRES(IBIN,ISUB,IORD)
     >                          *MYRES(IBIN,ISUB,IORD))
                           IF (NEFF.GE.2D0) THEN
                              WTDXMN(IBIN,ISUB,IORD) =
     >                             WTDXU2(IBIN,ISUB,IORD) /
     >                             (NEFF-1D0)
                           ELSE
                              WTDXMN(IBIN,ISUB,IORD) = -1D0
                           ENDIF
                           IF (WTDXU2(IBIN,ISUB,IORD).GE.0D0) THEN
                              WTDXU2(IBIN,ISUB,IORD) =
     >                             SQRT(WTDXU2(IBIN,ISUB,IORD))
                              WTDXL2(IBIN,ISUB,IORD) =
     >                             -WTDXU2(IBIN,ISUB,IORD)
                              IF (WTDXMN(IBIN,ISUB,IORD).GE.0D0) THEN
                                 WTDXMN(IBIN,ISUB,IORD) =
     >                                SQRT(WTDXMN(IBIN,ISUB,IORD))
                              ENDIF
                              IF (DABS(MYRES(IBIN,ISUB,IORD)).GT.1D-99)
     >                             THEN
                                 WTXMIN(IBIN,ISUB,IORD) =
     >                                WTXMIN(IBIN,ISUB,IORD) /
     >                                DABS(MYRES(IBIN,ISUB,IORD))
                                 WTXMAX(IBIN,ISUB,IORD) =
     >                                WTXMAX(IBIN,ISUB,IORD) /
     >                                DABS(MYRES(IBIN,ISUB,IORD))
                                 WTDXU2(IBIN,ISUB,IORD) =
     >                                WTDXU2(IBIN,ISUB,IORD) /
     >                                DABS(MYRES(IBIN,ISUB,IORD))
                                 WTDXL2(IBIN,ISUB,IORD) =
     >                                WTDXL2(IBIN,ISUB,IORD) /
     >                                DABS(MYRES(IBIN,ISUB,IORD))
                                 WTDXMN(IBIN,ISUB,IORD) =
     >                                WTDXMN(IBIN,ISUB,IORD) /
     >                                DABS(MYRES(IBIN,ISUB,IORD))
                                 WTDXUL(IBIN,ISUB,IORD) =
     >                                (WTXMAX(IBIN,ISUB,IORD) -
     >                                WTXMIN(IBIN,ISUB,IORD)) / 2D0
                              ELSE
                                 WTXMIN(IBIN,ISUB,IORD) = +1D99
                                 WTXMAX(IBIN,ISUB,IORD) = -1D99
                                 WTDXU2(IBIN,ISUB,IORD) = -1D0
                                 WTDXL2(IBIN,ISUB,IORD) = -1D0
                                 WTDXMN(IBIN,ISUB,IORD) = -1D0
                                 WTDXUL(IBIN,ISUB,IORD) = -1D0
                              ENDIF
                           ELSE
                              WTXMIN(IBIN,ISUB,IORD) = +1D99
                              WTXMAX(IBIN,ISUB,IORD) = -1D99
                              WTDXU2(IBIN,ISUB,IORD) = -1D0
                              WTDXL2(IBIN,ISUB,IORD) = -1D0
                              WTDXMN(IBIN,ISUB,IORD) = -1D0
                              WTDXUL(IBIN,ISUB,IORD) = -1D0
                           ENDIF
                        ENDIF
ckr Minimax method (as needed for scale uncertainties)
                     ELSEIF (IMODE.EQ.3) THEN
cdebug
Comment:                         IF (DABS(MYRES(IBIN,ISUB,IORD)).GT.1D-99) THEN
Comment:                            WTDXUM(IBIN,ISUB,IORD) =
Comment:      >                          WTDXUM(IBIN,ISUB,IORD) /
Comment:      >                          DABS(MYRES(IBIN,ISUB,IORD))
Comment:                            WTDXLM(IBIN,ISUB,IORD) =
Comment:      >                          WTDXLM(IBIN,ISUB,IORD) /
Comment:      >                          DABS(MYRES(IBIN,ISUB,IORD))
cdebug
                        IF (DABS(MYRESN(IBIN,ISUB,IORD)).GT.1D-99) THEN
                           WTDXUM(IBIN,ISUB,IORD) =
     >                          WTDXUM(IBIN,ISUB,IORD) /
     >                          DABS(MYRESN(IBIN,ISUB,IORD))
                           WTDXLM(IBIN,ISUB,IORD) =
     >                          WTDXLM(IBIN,ISUB,IORD) /
     >                          DABS(MYRESN(IBIN,ISUB,IORD))
                        ELSE
                           WTDXUM(IBIN,ISUB,IORD) = -1D0
                           WTDXLM(IBIN,ISUB,IORD) = -1D0
                        ENDIF
ckr Deviation from reference (as needed for algorithmic uncertainties)
                     ELSEIF (IMODE.EQ.4) THEN
                        IF (IRAP.LE.INT(NRAPIDITY/2)) THEN
                           IF (MYRES(IBIN+NBIN/2,ISUB,IORD).GT.
     >                          1D-99) THEN
                              WTDXUM(IBIN,ISUB,IORD) =
     >                             (WTX(IBIN,ISUB,IORD) /
     >                             DABS(MYRES(IBIN+NBIN/2,ISUB,IORD))) -
     >                             1D0
                           ELSEIF (MYRES(IBIN+NBIN/2,ISUB,IORD).LT.
     >                             -1D-99) THEN
                              WTDXUM(IBIN,ISUB,IORD) =
     >                             (WTX(IBIN,ISUB,IORD) /
     >                             DABS(MYRES(IBIN+NBIN/2,ISUB,IORD))) +
     >                             1D0
                           ELSE
                              WTDXUM(IBIN,ISUB,IORD) = -1D0
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
cdebug
Comment:          IBIN = 0
Comment:          DO IRAP=1,NRAPIDITY
Comment:             DO IPT=1,NPT(IRAP)
Comment:                IBIN = IBIN+1
Comment:                DO IORD=1,NORD+1
Comment:                   DO ISUB=1,NSUBPROC+1
Comment:                      if (iord.eq.3.and.isub.eq.8.and.ibin.eq.151)then
Comment:                      write(*,*)"ZZZZ: ibin,irap,ipt,iord,isub,i+n",
Comment:      >                    ibin,irap,ipt,iord,isub,ibin+nbin
Comment:                      write(*,*)"ENDW: wres,wtx,myres",
Comment:      >                    WRES(IBIN,ISUB,IORD),
Comment:      >                    WTX(IBIN,ISUB,IORD),
Comment:      >                    MYRES(IBIN,ISUB,IORD)
Comment:                      write(*,*)"ENDW: wres,wtx,myresnbin",
Comment:      >                    WRES(IBIN+NBIN,ISUB,IORD),
Comment:      >                    WTX(IBIN+NBIN,ISUB,IORD),
Comment:      >                    MYRES(IBIN+NBIN,ISUB,IORD)
Comment:                      write(*,*)"ZZZZ: myres, wres, wtdxmn",
Comment:      >                    myres(ibin,isub,iord),
Comment:      >                    wres(ibin,isub,iord),
Comment:      >                    wtdxmn(ibin,isub,iord)
Comment:                   endif
Comment:                   ENDDO
Comment:                ENDDO
Comment:             ENDDO
Comment:          ENDDO
cdebug

         RETURN
      ENDIF

      RETURN
      END



      SUBROUTINE CENRES(ISTEP,LRAT,LNRM,SCENARIO)

      IMPLICIT NONE
      INTEGER ISTEP
      LOGICAL LRAT,LNRM
      CHARACTER*255 SCENARIO
      INCLUDE "fnx9999.inc"
      INCLUDE "uncert.inc"
      INTEGER IBIN,IRAP,IPT,ISUB,IORD,NBIN,IBINN,IBNRM
      INTEGER NRAPIDITYN,NPTN(NRAPMAX),NORDN,NSUBPROCN



*---Initialization
ckr Table to normalize loaded
ckr      WRITE(*,*)"AAAAA: CENRES STEP = ",ISTEP
      IF (ISTEP.EQ.0.OR.ISTEP.EQ.3) THEN
         IBIN = 0
         DO IRAP=1,NRAPIDITY
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               DO IORD=1,NORD+1
                  DO ISUB=1,NSUBPROC+1
                     MYRES(IBIN,ISUB,IORD)  = 0D0
                     IF (ISTEP.EQ.0) THEN
                        MYRESN(IBIN,ISUB,IORD) = 0D0
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         NBIN = IBIN
         IF (LNRM) THEN
            DO IRAP=1,NRAPIDITY
               DO IPT=1,NPT(IRAP)
                  IBIN = IBIN+1
                  DO IORD=1,NORD+1
                     DO ISUB=1,NSUBPROC+1
                        MYRES(IBIN,ISUB,IORD)  = 0D0
                        IF (ISTEP.EQ.0) THEN
                           MYRESN(IBIN,ISUB,IORD) = 0D0
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
ckr Save basic structure of table to normalize
            IF (ISTEP.EQ.0) THEN
               NRAPIDITYN = NRAPIDITY
               DO IRAP=1,NRAPIDITY
                  NPTN(IRAP) = NPT(IRAP)
               ENDDO
               NORDN      = NORD
               NSUBPROCN  = NSUBPROC
            ENDIF
         ENDIF
      ENDIF

ckr Normalization table loaded, if necessary
      IF ((ISTEP.EQ.1.OR.ISTEP.EQ.4).AND.LNRM) THEN
         IF (SCENARIO(1:8).EQ."fnl2332c") THEN
            IBIN = 0
            DO IRAP=1,NRAPIDITYN
               DO IPT=1,NPTN(IRAP)
                  IBIN = IBIN+1
                  IF (IPT.EQ.1) IBNRM=IBIN
                  DO IORD=1,NORDN+1
                     DO ISUB=1,NSUBPROCN+1
ckr Only first pT bin needed to normalize in each rapidity region
ckr Might be required to sum up multiple bins in pT and y for CMS analysis
ckr                        IF (IPT.EQ.1) THEN
                           IF (IORD.LE.NORDN) THEN
                              MYTMP(IBIN,ISUB,IORD) =
     >                             RESULT(IBIN,ISUB,IORD)
                           ELSE
                              MYTMP(IBIN,ISUB,IORD) =
     >                             (RESULT(IBIN,ISUB,1) +
     >                             RESULT(IBIN,ISUB,2))
                           ENDIF
ckr                        ENDIF
                           IF (IPT.EQ.13) THEN
ckr                           IF (IPT.EQ.1) THEN
                              MYTMP(IBNRM,ISUB,IORD) =
     >                             MYTMP(IBIN,ISUB,IORD)
                           ENDIF
                           IF (13.LT.IPT.AND.IPT.LE.15) THEN
                              MYTMP(IBNRM,ISUB,IORD) =
     >                             MYTMP(IBNRM,ISUB,IORD) +
     >                             MYTMP(IBIN,ISUB,IORD)
                           ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            IBIN = 0
ckr Counter for bin containing normalization factor
            IBNRM = 1
            DO IRAP=1,NRAPIDITYN
               DO IPT=1,NPTN(IRAP)
                  IBIN = IBIN+1
                  DO IORD=1,NORDN+1
                     DO ISUB=1,NSUBPROCN+1
                        MYRES(IBIN+NBIN,ISUB,IORD) =
     >                       1D0/MYTMP(IBNRM,ISUB,IORD)
                        IF (ISTEP.EQ.1) THEN
                           MYRESN(IBIN+NBIN,ISUB,IORD) =
     >                          MYRES(IBIN+NBIN,ISUB,IORD)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
               IBNRM = IBNRM + NPTN(IRAP)
            ENDDO
         ELSEIF (SCENARIO(1:7).EQ."fnl2442") THEN
            IBIN = 0
            DO IRAP=1,NRAPIDITYN
               DO IPT=1,NPTN(IRAP)
                  IBIN = IBIN+1
                  DO IORD=1,NORDN+1
                     DO ISUB=1,NSUBPROCN+1
                        IF (IRAP.EQ.2) THEN
                           IF (IORD.LE.NORDN) THEN
                              MYTMP(IBIN,ISUB,IORD) =
     >                             RESULT(IBIN,ISUB,IORD)
                           ELSE
                              MYTMP(IBIN,ISUB,IORD) =
     >                             (RESULT(IBIN,ISUB,1) +
     >                             RESULT(IBIN,ISUB,2))
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            IBIN = 0
            DO IRAP=1,NRAPIDITYN
               DO IPT=1,NPTN(IRAP)
                  IBIN = IBIN+1
                  DO IORD=1,NORDN+1
                     DO ISUB=1,NSUBPROCN+1
                        MYRES(IBIN+NBIN,ISUB,IORD) =
     >                       1D0/MYTMP(NPTN(1)+IPT,ISUB,IORD)
                        IF (ISTEP.EQ.1) THEN
                           MYRESN(IBIN+NBIN,ISUB,IORD) =
     >                          MYRES(IBIN+NBIN,ISUB,IORD)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ELSEIF (SCENARIO(1:7).EQ."fnl2522") THEN
            IBINN = 1
            IF (SCENARIO(8:14).EQ."diffpt1") THEN
               IBINN = 1
            ELSEIF (SCENARIO(8:14).EQ."diffpt2") THEN
               IBINN = 2
            ELSEIF (SCENARIO(8:14).EQ."diffpt3") THEN
               IBINN = 3
            ELSEIF (SCENARIO(8:14).EQ."diffpt4") THEN
               IBINN = 4
            ELSEIF (SCENARIO(8:14).EQ."diffpt5") THEN
               IBINN = 5
            ENDIF
            IBIN = 0
            DO IRAP=1,NRAPIDITYN
               DO IPT=1,NPTN(IRAP)
                  IBIN = IBIN+1
                  DO IORD=1,NORDN+1
                     DO ISUB=1,NSUBPROCN+1
                        IF (IORD.LE.NORDN) THEN
                           MYRES(IBIN+NBIN,ISUB,IORD) =
     >                          1D0/RESULT(IBINN,ISUB,IORD)
                           IF (ISTEP.EQ.1) THEN
                              MYRESN(IBIN+NBIN,ISUB,IORD) =
     >                             MYRES(IBIN+NBIN,ISUB,IORD)
                           ENDIF
                        ELSE
                           MYRES(IBIN+NBIN,ISUB,IORD) =
     >                          1D0/(RESULT(IBINN,ISUB,1) +
     >                          RESULT(IBINN,ISUB,2))
                           IF (ISTEP.EQ.1) THEN
                              MYRESN(IBIN+NBIN,ISUB,IORD) =
     >                             MYRES(IBIN+NBIN,ISUB,IORD)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ELSEIF (SCENARIO(1:7).EQ."fnl2622".OR.
     >           SCENARIO(1:7).EQ."fnl2652") THEN
            IBIN = 0
            DO IRAP=1,NRAPIDITYN
               DO IPT=1,NPTN(IRAP)
                  IBIN = IBIN+1
                  DO IORD=1,NORDN+1
                     DO ISUB=1,NSUBPROCN+1
                        IF (IPT.EQ.1) THEN
                           MYTMP(IRAP,ISUB,IORD) = 0.D0
Comment:                            write(*,*)"AA ibin,iord,isub,irap,ipt,"//
Comment:      >                          "mytmp,drap,dpt",
Comment:      >                          ibin,iord,isub,irap,ipt,
Comment:      >                          mytmp(irap,isub,iord),
Comment:      >                          rapbin(irap+1)-rapbin(irap),
Comment:      >                          ptbin(irap,ipt+1)-ptbin(irap,ipt)
                        ENDIF
ckr Attention: For CMS Publ. normalize only up to 16 in Chi like in exp. analysis!
                        IF ((SCENARIO(1:7).EQ."fnl2622".AND.IPT.LT.13)
     >                       .OR.(SCENARIO(1:7).EQ."fnl2652"))THEN
                           IF (IORD.LE.NORDN) THEN
                              MYTMP(IRAP,ISUB,IORD) =
     >                             MYTMP(IRAP,ISUB,IORD) +
     >                             RESULT(IBIN,ISUB,IORD) *
     >                             (PTBIN(IRAP,IPT+1)-PTBIN(IRAP,IPT))
                           ELSE
                              MYTMP(IRAP,ISUB,IORD) =
     >                             MYTMP(IRAP,ISUB,IORD) +
     >                             (RESULT(IBIN,ISUB,1) +
     >                             RESULT(IBIN,ISUB,2)) *
     >                             (PTBIN(IRAP,IPT+1)-PTBIN(IRAP,IPT))
                           ENDIF
                        ENDIF
Comment:                         write(*,*)"BB ibin,iord,isub,irap,ipt,mytmp",
Comment:      >                       ibin,iord,isub,irap,ipt,
Comment:      >                       mytmp(irap,isub,iord)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            IBIN = 0
            DO IRAP=1,NRAPIDITYN
               DO IPT=1,NPTN(IRAP)
                  IBIN = IBIN+1
                  DO IORD=1,NORDN+1
                     DO ISUB=1,NSUBPROCN+1
                        MYRES(IBIN+NBIN,ISUB,IORD) =
     >                       1D0/MYTMP(IRAP,ISUB,IORD)
                        IF (ISTEP.EQ.1) THEN
                           MYRESN(IBIN+NBIN,ISUB,IORD) =
     >                          MYRES(IBIN+NBIN,ISUB,IORD)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ELSEIF (SCENARIO(1:7).EQ."fnl2722".OR.
     >           SCENARIO(1:7).EQ."fnl2732".OR.
     >           SCENARIO(1:7).EQ."fnl2742") THEN
            IBIN = 0
            DO IRAP=1,NRAPIDITYN
               DO IPT=1,NPTN(IRAP)
                  IBIN = IBIN+1
                  DO IORD=1,NORDN+1
                     DO ISUB=1,NSUBPROCN+1
                        IF (IORD.LE.NORDN) THEN
                           MYRES(IBIN+NBIN,ISUB,IORD) =
     >                          1D0/RESULT(IBIN,ISUB,IORD)
                           IF (ISTEP.EQ.1) THEN
                              MYRESN(IBIN+NBIN,ISUB,IORD) =
     >                             MYRES(IBIN+NBIN,ISUB,IORD)
                           ENDIF
                        ELSE
                           MYRES(IBIN+NBIN,ISUB,IORD) =
     >                          1D0/(RESULT(IBIN,ISUB,1) +
     >                          RESULT(IBIN,ISUB,2))
                           IF (ISTEP.EQ.1) THEN
                              MYRESN(IBIN+NBIN,ISUB,IORD) =
     >                             MYRES(IBIN+NBIN,ISUB,IORD)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF



ckr Reload table to normalize
      IF (ISTEP.EQ.2.OR.ISTEP.EQ.5) THEN
*---Fill array
      IBIN = 0
      DO IRAP=1,NRAPIDITY
         DO IPT=1,NPT(IRAP)
            IBIN = IBIN+1
            DO IORD=1,NORD
               DO ISUB=1,NSUBPROC
                  MYRES(IBIN,ISUB,IORD) = RESULT(IBIN,ISUB,IORD)
                  MYRES(IBIN,NSUBPROC+1,IORD) =
     >                 MYRES(IBIN,NSUBPROC+1,IORD) +
     >                 MYRES(IBIN,ISUB,IORD)
ckr Centrality ratio: LRAT
ckr Normalization: LNRM
                  IF (LRAT) THEN
                     MYRES(IBIN+NBIN,ISUB,IORD) =
     >                    MYRES(IPT,ISUB,IORD) /
     >                    RESULT(IBIN,ISUB,IORD)
Comment:                      write(*,*)"LRATA: ibin,iord,isub,irap,ipt,a,b,a/b",
Comment:      >                    ibin,iord,isub,irap,ipt,
Comment:      >                    MYRES(IPT,ISUB,IORD),
Comment:      >                    RESULT(IBIN,ISUB,IORD),
Comment:      >                    MYRES(IBIN+NBIN,ISUB,IORD)
                  ELSEIF (LNRM) THEN
                     IF (ISTEP.EQ.2.OR.ISTEP.EQ.5) THEN
Comment:                         write(*,*)"EE: ibin,iord,isub,irap,ipt,mi,mi+n",
Comment:      >                       ibin,iord,isub,irap,ipt,
Comment:      >                       MYRES(IBIN,ISUB,IORD),
Comment:      >                       MYRES(IBIN+NBIN,ISUB,IORD)
                        MYRES(IBIN+NBIN,ISUB,IORD) =
     >                       MYRES(IBIN,ISUB,IORD) *
     >                       MYRES(IBIN+NBIN,ISUB,IORD)
Comment:                         write(*,*)"FF: ibin,iord,isub,irap,ipt,mi,mi+n",
Comment:      >                       ibin,iord,isub,irap,ipt,
Comment:      >                       MYRES(IBIN,ISUB,IORD),
Comment:      >                       MYRES(IBIN+NBIN,ISUB,IORD)
                     ENDIF
                  ENDIF
                  IF (ISTEP.EQ.2) THEN
                     MYRESN(IBIN,ISUB,IORD) =
     >                    MYRES(IBIN,ISUB,IORD)
                     MYRESN(IBIN+NBIN,ISUB,IORD) =
     >                    MYRES(IBIN+NBIN,ISUB,IORD)
                  ENDIF
               ENDDO
               MYRES(IBIN,NSUBPROC+1,NORD+1) =
     >              MYRES(IBIN,NSUBPROC+1,NORD+1) +
     >              MYRES(IBIN,NSUBPROC+1,IORD)
               IF (LRAT) THEN
                  MYRES(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                 MYRES(IPT,NSUBPROC+1,IORD) /
     >                 MYRES(IBIN,NSUBPROC+1,IORD)
Comment:                   write(*,*)"LRATB: ibin,iord,isub,irap,ipt,a,b,a/b",
Comment:      >                 ibin,iord,isub,irap,ipt,
Comment:      >                 MYRES(IPT,NSUBPROC+1,IORD),
Comment:      >                 MYRES(IBIN,NSUBPROC+1,IORD),
Comment:      >                 MYRES(IBIN+NBIN,NSUBPROC+1,IORD)
               ELSEIF (LNRM) THEN
                  IF (ISTEP.EQ.2.OR.ISTEP.EQ.5) THEN
                     MYRES(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                    MYRES(IBIN,NSUBPROC+1,IORD) *
     >                    MYRES(IBIN+NBIN,NSUBPROC+1,IORD)
                  ENDIF
               ENDIF
               IF (ISTEP.EQ.2) THEN
                  MYRESN(IBIN,NSUBPROC+1,IORD) =
     >                 MYRES(IBIN,NSUBPROC+1,IORD)
                  MYRESN(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                 MYRES(IBIN+NBIN,NSUBPROC+1,IORD)
               ENDIF
            ENDDO
            IF (LRAT) THEN
               MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1) =
     >              MYRES(IPT,NSUBPROC+1,NORD+1) /
     >              MYRES(IBIN,NSUBPROC+1,NORD+1)
Comment:                write(*,*)"LRATC: ibin,iord,isub,irap,ipt,a,b,a/b",
Comment:      >              ibin,iord,isub,irap,ipt,
Comment:      >              MYRES(IPT,NSUBPROC+1,NORD+1),
Comment:      >              MYRES(IBIN,NSUBPROC+1,NORD+1),
Comment:      >              MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1)
            ELSEIF (LNRM) THEN
               IF (ISTEP.EQ.2.OR.ISTEP.EQ.5) THEN
                  MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1) =
     >                 MYRES(IBIN,NSUBPROC+1,NORD+1) *
     >                 MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1)
               ENDIF
            ENDIF
            IF (ISTEP.EQ.2) THEN
               MYRESN(IBIN,NSUBPROC+1,NORD+1) =
     >              MYRES(IBIN,NSUBPROC+1,NORD+1)
               MYRESN(IBIN+NBIN,NSUBPROC+1,NORD+1) =
     >              MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1)
            ENDIF
            DO ISUB=1,NSUBPROC
               DO IORD=1,NORD
                  MYRES(IBIN,ISUB,NORD+1) =
     >                 MYRES(IBIN,ISUB,NORD+1) +
     >                 MYRES(IBIN,ISUB,IORD)
               ENDDO
               IF (LRAT) THEN
                  MYRES(IBIN+NBIN,ISUB,NORD+1) =
     >                 MYRES(IPT,ISUB,NORD+1) /
     >                 MYRES(IBIN,ISUB,NORD+1)
Comment:                   write(*,*)"LRATD: ibin,iord,isub,irap,ipt,a,b,a/b",
Comment:      >                 ibin,iord,isub,irap,ipt,
Comment:      >                 MYRES(IPT,ISUB,NORD+1),
Comment:      >                 MYRES(IBIN,ISUB,NORD+1),
Comment:      >                 MYRES(IBIN+NBIN,ISUB,NORD+1)
               ELSEIF (LNRM) THEN
                  IF (ISTEP.EQ.2.OR.ISTEP.EQ.5) THEN
                     MYRES(IBIN+NBIN,ISUB,NORD+1) =
     >                    MYRES(IBIN,ISUB,NORD+1) *
     >                    MYRES(IBIN+NBIN,ISUB,NORD+1)
                  ENDIF
               ENDIF
               IF (ISTEP.EQ.2) THEN
                  MYRESN(IBIN,ISUB,NORD+1) =
     >                 MYRES(IBIN,ISUB,NORD+1)
                  MYRESN(IBIN+NBIN,ISUB,NORD+1) =
     >                 MYRES(IBIN+NBIN,ISUB,NORD+1)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      ENDIF

cdebug
Comment:       WRITE(*,*)"BBBBB: CENRES STEP = ",ISTEP
Comment:       IBIN = 0
Comment:       DO IRAP=1,NRAPIDITY
Comment:          DO IPT=1,NPT(IRAP)
Comment:             IBIN = IBIN+1
Comment:             DO IORD=1,NORD+1
Comment:                DO ISUB=1,NSUBPROC+1
Comment:                   WRITE(*,*)"DEBUG: IBIN,IRAP,IPT,IORD,ISUB,I+N",
Comment:      >                 IBIN,IRAP,IPT,IORD,ISUB,IBIN+NBIN
Comment:                   WRITE(*,*)"DEBUG: "//
Comment:      >                 "MYRES(IBIN), MYRESN(IBIN)",
Comment:      >                 MYRES(IBIN,ISUB,IORD),MYRESN(IBIN,ISUB,IORD)
Comment:                   WRITE(*,*)"DEBUG: "//
Comment:      >                 "MYRES(IBIN+NBIN), MYRESN(IBIN+NBIN)",
Comment:      >                 MYRES(IBIN+NBIN,ISUB,IORD),
Comment:      >                 MYRESN(IBIN+NBIN,ISUB,IORD)
Comment:                ENDDO
Comment:             ENDDO
Comment:          ENDDO
Comment:       ENDDO
cdebug

      RETURN
      END



      SUBROUTINE SUMMUP(IVAR,IBIN,ISUB,IORD,SUMM,DIFF,WEIGHT,
     >     LRAT,IBREF,NBIN,LNRM)

      IMPLICIT NONE
      INTEGER IVAR,IBIN,IORD,ISUB,IBREF,NBIN
      DOUBLE PRECISION SUMM,DIFF,WEIGHT
      LOGICAL LRAT,LNRM
      INCLUDE "fnx9999.inc"
      INCLUDE "uncert.inc"
      DOUBLE PRECISION RATIO,RDIFF

      WRES(IBIN,ISUB,IORD) = SUMM
      NTN(IBIN,ISUB,IORD)  =
     >     NTN(IBIN,ISUB,IORD)  +
     >     1
      WT(IBIN,ISUB,IORD)   =
     >     WT(IBIN,ISUB,IORD)   +
     >     WEIGHT
      WT2(IBIN,ISUB,IORD)  =
     >     WT2(IBIN,ISUB,IORD)  +
     >     WEIGHT*WEIGHT
      WTX(IBIN,ISUB,IORD) =
     >     WTX(IBIN,ISUB,IORD)  +
     >     SUMM*WEIGHT
      WTX2(IBIN,ISUB,IORD) =
     >     WTX2(IBIN,ISUB,IORD) +
     >     SUMM*SUMM*WEIGHT
cdebug
Comment:       WRITE(*,*)"ivar,ibin,isub,iord",ivar,ibin,isub,iord
Comment:       WRITE(*,*)"wres,ntn,wt,wt2,wtx,wtx2",
Comment:      >     wres(ibin,isub,iord),ntn(ibin,isub,iord),
Comment:      >     wt(ibin,isub,iord),wt2(ibin,isub,iord),
Comment:      >     wtx(ibin,isub,iord),wtx2(ibin,isub,iord)
cdebug
      IF (SUMM.LT.WTXMIN(IBIN,ISUB,IORD)) THEN
         WTXMIN(IBIN,ISUB,IORD) = SUMM
         IF (ISUB.EQ.NSUBPROC+1.AND.IORD.EQ.NORD+1) THEN
            IJMIN(IBIN) = IVAR
         ENDIF
      ENDIF
      IF (SUMM.GT.WTXMAX(IBIN,ISUB,IORD)) THEN
         WTXMAX(IBIN,ISUB,IORD) = SUMM
         IF (ISUB.EQ.NSUBPROC+1.AND.IORD.EQ.NORD+1) THEN
            IJMAX(IBIN) = IVAR
         ENDIF
      ENDIF
      IF (DIFF.GT.0D0) THEN
         WTDXU2(IBIN,ISUB,IORD) =
     >        WTDXU2(IBIN,ISUB,IORD) + DIFF*DIFF
         IF (DIFF.GT.WTDXUM(IBIN,ISUB,IORD)) THEN
            WTDXUM(IBIN,ISUB,IORD) = DIFF
         ENDIF
      ELSE
         WTDXL2(IBIN,ISUB,IORD) =
     >        WTDXL2(IBIN,ISUB,IORD) + DIFF*DIFF
         IF (DIFF.LT.WTDXLM(IBIN,ISUB,IORD)) THEN
            WTDXLM(IBIN,ISUB,IORD) = DIFF
         ENDIF
      ENDIF
      IF (LRAT) THEN
         RATIO =
     >        WRES(IBREF,ISUB,IORD) /
     >        WRES(IBIN,ISUB,IORD)
         WRES(IBIN+NBIN,ISUB,IORD) = RATIO
         NTN(IBIN+NBIN,ISUB,IORD) =
     >        NTN(IBIN+NBIN,ISUB,IORD)  +
     >        1
         WT(IBIN+NBIN,ISUB,IORD) =
     >        WT(IBIN+NBIN,ISUB,IORD)   +
     >        WEIGHT
         WT2(IBIN+NBIN,ISUB,IORD) =
     >        WT2(IBIN+NBIN,ISUB,IORD)  +
     >        WEIGHT*WEIGHT
         WTX(IBIN+NBIN,ISUB,IORD) =
     >        WTX(IBIN+NBIN,ISUB,IORD)  +
     >        RATIO*WEIGHT
         WTX2(IBIN+NBIN,ISUB,IORD) =
     >        WTX2(IBIN+NBIN,ISUB,IORD) +
     >        RATIO*RATIO*WEIGHT
         IF (RATIO.LT.WTXMIN(IBIN+NBIN,ISUB,IORD)) THEN
            WTXMIN(IBIN+NBIN,ISUB,IORD) = RATIO
            IF (ISUB.EQ.NSUBPROC+1.AND.IORD.EQ.NORD+1) THEN
               IJMIN(IBIN+NBIN) = IVAR
            ENDIF
         ENDIF
         IF (RATIO.GT.WTXMAX(IBIN+NBIN,ISUB,IORD)) THEN
            WTXMAX(IBIN+NBIN,ISUB,IORD) = RATIO
            IF (ISUB.EQ.NSUBPROC+1.AND.IORD.EQ.NORD+1) THEN
               IJMAX(IBIN+NBIN) = IVAR
            ENDIF
         ENDIF
         RDIFF = RATIO -
     >        MYRES(IBIN+NBIN,ISUB,IORD)
         IF (RDIFF.GT.0D0) THEN
            WTDXU2(IBIN+NBIN,ISUB,IORD) =
     >           WTDXU2(IBIN+NBIN,ISUB,IORD) + RDIFF*RDIFF
            IF (RDIFF.GT.WTDXUM(IBIN+NBIN,ISUB,IORD)) THEN
               WTDXUM(IBIN+NBIN,ISUB,IORD) = RDIFF
            ENDIF
         ELSE
            WTDXL2(IBIN+NBIN,ISUB,IORD) =
     >           WTDXL2(IBIN+NBIN,ISUB,IORD) + RDIFF*RDIFF
            IF (RDIFF.LT.WTDXLM(IBIN+NBIN,ISUB,IORD)) THEN
               WTDXLM(IBIN+NBIN,ISUB,IORD) = RDIFF
            ENDIF
         ENDIF
      ENDIF
Comment:       ELSEIF (LNRM) THEN
Comment:          RATIO = WRES(IBIN+NBIN,ISUB,IORD)
Comment: ckr         WRES(IBIN+NBIN,ISUB,IORD) = RATIO
Comment: ckr         RATIO =
Comment: ckr     >        WRES(IBIN,ISUB,IORD) /
Comment: ckr     >        WRES(IBIN+NBIN,ISUB,IORD)
Comment: ckr         WRES(IBIN+NBIN,ISUB,IORD) = RATIO
Comment:          NTN(IBIN+NBIN,ISUB,IORD) =
Comment:      >        NTN(IBIN+NBIN,ISUB,IORD)  +
Comment:      >        1
Comment:          WT(IBIN+NBIN,ISUB,IORD) =
Comment:      >        WT(IBIN+NBIN,ISUB,IORD)   +
Comment:      >        WEIGHT
Comment:          WT2(IBIN+NBIN,ISUB,IORD) =
Comment:      >        WT2(IBIN+NBIN,ISUB,IORD)  +
Comment:      >        WEIGHT*WEIGHT
Comment:          WTX(IBIN+NBIN,ISUB,IORD) =
Comment:      >        WTX(IBIN+NBIN,ISUB,IORD)  +
Comment:      >        RATIO*WEIGHT
Comment:          WTX2(IBIN+NBIN,ISUB,IORD) =
Comment:      >        WTX2(IBIN+NBIN,ISUB,IORD) +
Comment:      >        RATIO*RATIO*WEIGHT
Comment:          IF (RATIO.LT.WTXMIN(IBIN+NBIN,ISUB,IORD)) THEN
Comment:             WTXMIN(IBIN+NBIN,ISUB,IORD) = RATIO
Comment:             IF (ISUB.EQ.NSUBPROC+1.AND.IORD.EQ.NORD+1) THEN
Comment:                IJMIN(IBIN+NBIN) = IVAR
Comment:             ENDIF
Comment:          ENDIF
Comment:          IF (RATIO.GT.WTXMAX(IBIN+NBIN,ISUB,IORD)) THEN
Comment:             WTXMAX(IBIN+NBIN,ISUB,IORD) = RATIO
Comment:             IF (ISUB.EQ.NSUBPROC+1.AND.IORD.EQ.NORD+1) THEN
Comment:                IJMAX(IBIN+NBIN) = IVAR
Comment:             ENDIF
Comment:          ENDIF
Comment:          RDIFF = RATIO -
Comment:      >        MYRES(IBIN+NBIN,ISUB,IORD)
Comment:          IF (RDIFF.GT.0D0) THEN
Comment:             WTDXU2(IBIN+NBIN,ISUB,IORD) =
Comment:      >           WTDXU2(IBIN+NBIN,ISUB,IORD) + RDIFF*RDIFF
Comment:             IF (RDIFF.GT.WTDXUM(IBIN+NBIN,ISUB,IORD)) THEN
Comment:                WTDXUM(IBIN+NBIN,ISUB,IORD) = RDIFF
Comment:             ENDIF
Comment:          ELSE
Comment:             WTDXL2(IBIN+NBIN,ISUB,IORD) =
Comment:      >           WTDXL2(IBIN+NBIN,ISUB,IORD) + RDIFF*RDIFF
Comment:             IF (RDIFF.LT.WTDXLM(IBIN+NBIN,ISUB,IORD)) THEN
Comment:                WTDXLM(IBIN+NBIN,ISUB,IORD) = RDIFF
Comment:             ENDIF
Comment:          ENDIF
Comment:       ENDIF

      RETURN
      END
