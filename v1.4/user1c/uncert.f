      SUBROUTINE UNCERT(ISTAT,IMODE)
      
      IMPLICIT NONE
      INTEGER ISTAT,IMODE
      INCLUDE "fnx9999.inc"
      INCLUDE "uncert.inc"
      DOUBLE PRECISION DIFF,SUMM,DIFFPROC,SUMMPROC,DIFFORD,SUMMORD
      INTEGER IBIN,IRAP,IPT,IORD,ISUB,NBIN,NRAP

      

*---Initialization	
      IF (ISTAT.EQ.1) THEN
         IBIN = 0
         NRAP = NRAPIDITY
         IF (IMODE.EQ.5) NRAP=2*NRAPIDITY
         DO IRAP=1,NRAP
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               DO IORD=1,NORD+1
                  DO ISUB=1,NSUBPROC+1
                     WT(IBIN,ISUB,IORD)     = 0.D0
                     WT2(IBIN,ISUB,IORD)    = 0.D0
                     WTX(IBIN,ISUB,IORD)    = 0.D0
                     WTX2(IBIN,ISUB,IORD)   = 0.D0
                     WTDXL2(IBIN,ISUB,IORD) = 0.D0
                     WTDXU2(IBIN,ISUB,IORD) = 0.D0
                     WTDXLM(IBIN,ISUB,IORD) = 0.D0
                     WTDXUM(IBIN,ISUB,IORD) = 0.D0
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         NBIN = IBIN
         RETURN
      ENDIF
      
        
        
*---Loop
      IF (ISTAT.EQ.2) THEN
         IBIN = 0
         DO IRAP=1,NRAPIDITY
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               SUMMORD = 0.D0
               DIFFORD = 0.D0
               DO IORD=1,NORD
                  SUMMPROC = 0.D0
                  DIFFPROC = 0.D0
                  DO ISUB=1,NSUBPROC
                     SUMM = RESULT(IBIN,ISUB,IORD)
                     SUMMPROC = SUMMPROC + SUMM
                     DIFF = SUMM - MYRES(IBIN,ISUB,IORD)
                     DIFFPROC = DIFFPROC + DIFF
                     WT(IBIN,ISUB,IORD)   =
     >                    WT(IBIN,ISUB,IORD)   + 1.D0 
                     WT2(IBIN,ISUB,IORD)  =
     >                    WT2(IBIN,ISUB,IORD)  + 1.D0 
                     WTX(IBIN,ISUB,IORD) =
     >                    WTX(IBIN,ISUB,IORD)  + RESULT(IBIN,ISUB,IORD)
                     WTX2(IBIN,ISUB,IORD) =
     >                    WTX2(IBIN,ISUB,IORD) +
     >                    RESULT(IBIN,ISUB,IORD)*RESULT(IBIN,ISUB,IORD)
                     IF (DIFF.GT.0D0) THEN
                        WTDXU2(IBIN,ISUB,IORD) =
     >                       WTDXU2(IBIN,ISUB,IORD) + DIFF*DIFF
                        IF (DIFF.GT.WTDXUM(IBIN,ISUB,IORD)) THEN
                           WTDXUM(IBIN,ISUB,IORD) = DIFF
                        ENDIF
                     ELSE
                        WTDXL2(IBIN,ISUB,IORD) =
     >                       WTDXL2(IBIN,ISUB,IORD) + DIFF*DIFF
                        IF (DIFF.LT.WTDXLM(IBIN,ISUB,IORD)) THEN
                           WTDXLM(IBIN,ISUB,IORD) = DIFF
                        ENDIF
                     ENDIF
                     IF (IMODE.EQ.5) THEN
                        WTX(IBIN+NBIN,ISUB,IORD) = 
     >                       WTX(IPT,ISUB,IORD) /
     >                       RESULT(IBIN,ISUB,IORD)
                        DIFF = WTX(IBIN+NBIN,ISUB,IORD) -
     >                       MYRES(IBIN+NBIN,ISUB,IORD)
                        IF (DIFF.GT.0D0) THEN
                           WTDXU2(IBIN+NBIN,ISUB,IORD) =
     >                          WTDXU2(IBIN+NBIN,ISUB,IORD) + DIFF*DIFF
                           IF (DIFF.GT.WTDXUM(IBIN+NBIN,ISUB,IORD)) THEN
                              WTDXUM(IBIN+NBIN,ISUB,IORD) = DIFF
                           ENDIF
                        ELSE
                           WTDXL2(IBIN+NBIN,ISUB,IORD) =
     >                          WTDXL2(IBIN+NBIN,ISUB,IORD) + DIFF*DIFF
                           IF (DIFF.LT.WTDXLM(IBIN+NBIN,ISUB,IORD)) THEN
                              WTDXLM(IBIN+NBIN,ISUB,IORD) = DIFF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
                  SUMM = SUMMPROC
                  SUMMORD = SUMMORD + SUMMPROC
                  DIFF = SUMM - MYRES(IBIN,NSUBPROC+1,IORD)
                  DIFFORD = DIFFORD + DIFF
                  WT(IBIN,NSUBPROC+1,IORD)    =
     >                 WT(IBIN,NSUBPROC+1,IORD)   + 1.D0 
                  WT2(IBIN,NSUBPROC+1,IORD)   =
     >                 WT2(IBIN,NSUBPROC+1,IORD)  + 1.D0 
                  WTX(IBIN,NSUBPROC+1,IORD)   =
     >                 WTX(IBIN,NSUBPROC+1,IORD)  + SUMMPROC
                  WTX2(IBIN,NSUBPROC+1,IORD)  =
     >                 WTX2(IBIN,NSUBPROC+1,IORD) + SUMMPROC*SUMMPROC
                  IF (DIFFPROC.GT.0D0) THEN
                     WTDXU2(IBIN,NSUBPROC+1,IORD) =
     >                    WTDXU2(IBIN,NSUBPROC+1,IORD) +
     >                    DIFFPROC*DIFFPROC
                     IF (DIFFPROC.GT.WTDXUM(IBIN,NSUBPROC+1,IORD)) THEN
                        WTDXUM(IBIN,NSUBPROC+1,IORD) = DIFFPROC
                     ENDIF
                  ELSE
                     WTDXL2(IBIN,NSUBPROC+1,IORD) =
     >                    WTDXL2(IBIN,NSUBPROC+1,IORD) +
     >                    DIFFPROC*DIFFPROC
                     IF (DIFFPROC.LT.WTDXLM(IBIN,NSUBPROC+1,IORD)) THEN
                        WTDXLM(IBIN,NSUBPROC+1,IORD) = DIFFPROC
                     ENDIF
                  ENDIF
                  IF (IMODE.EQ.5) THEN
                     WTX(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                    WTX(IPT,NSUBPROC+1,IORD) / 
     >                    WTX(IBIN,NSUBPROC+1,IORD)
                     DIFF = WTX(IBIN+NBIN,NSUBPROC+1,IORD) -
     >                    MYRES(IBIN+NBIN,NSUBPROC+1,IORD)
                     IF (DIFF.GT.0D0) THEN
                        WTDXU2(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                       WTDXU2(IBIN+NBIN,NSUBPROC+1,IORD) +
     >                       DIFF*DIFF
                        IF (DIFF.GT.WTDXUM(IBIN+NBIN,NSUBPROC+1,IORD))
     >                       THEN
                           WTDXUM(IBIN+NBIN,NSUBPROC+1,IORD) = DIFF
                        ENDIF
                     ELSE
                        WTDXL2(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                       WTDXL2(IBIN+NBIN,NSUBPROC+1,IORD) +
     >                       DIFF*DIFF
                        IF (DIFF.LT.WTDXLM(IBIN+NBIN,NSUBPROC+1,IORD))
     >                       THEN
                           WTDXLM(IBIN+NBIN,NSUBPROC+1,IORD) = DIFF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
               WT(IBIN,NSUBPROC+1,NORD+1)   =
     >              WT(IBIN,NSUBPROC+1,NORD+1)   + 1.D0 
               WT2(IBIN,NSUBPROC+1,NORD+1)  =
     >              WT2(IBIN,NSUBPROC+1,NORD+1)  + 1.D0 
               WTX(IBIN,NSUBPROC+1,NORD+1)  = 
     >              WTX(IBIN,NSUBPROC+1,NORD+1)  + SUMMORD
               WTX2(IBIN,NSUBPROC+1,NORD+1) =
     >              WTX2(IBIN,NSUBPROC+1,NORD+1) + SUMMORD*SUMMORD
               IF (DIFFORD.GT.0D0) THEN
                  WTDXU2(IBIN,NSUBPROC+1,NORD+1) =
     >                 WTDXU2(IBIN,NSUBPROC+1,NORD+1) +
     >                 DIFFORD*DIFFORD
                  IF (DIFFORD.GT.WTDXUM(IBIN,NSUBPROC+1,NORD+1)) THEN
                     WTDXUM(IBIN,NSUBPROC+1,NORD+1) = DIFFORD
                  ENDIF
               ELSE
                  WTDXL2(IBIN,NSUBPROC+1,NORD+1) =
     >                 WTDXL2(IBIN,NSUBPROC+1,NORD+1) +
     >                 DIFFORD*DIFFORD
                  IF (DIFFORD.LT.WTDXLM(IBIN,NSUBPROC+1,NORD+1)) THEN
                     WTDXLM(IBIN,NSUBPROC+1,NORD+1) = DIFFORD
                  ENDIF
               ENDIF
               IF (IMODE.EQ.5) THEN
                  WTX(IBIN+NBIN,NSUBPROC+1,NORD+1) = 
     >                 WTX(IPT,NSUBPROC+1,NORD+1) /
     >                 WTX(IBIN,NSUBPROC+1,NORD+1)
                  DIFF = WTX(IBIN+NBIN,NSUBPROC+1,NORD+1) -
     >                 MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1)
                  write(*,*)"irap,ipt,iord,isub,ibin",irap,ipt,iord,isub
     >                 ,ibin
                  write(*,*)"wtx1,2,2",WTX(IBIN+NBIN,NSUBPROC+1,NORD+1),
     >                 WTX(IPT,NSUBPROC+1,NORD+1),
     >                 WTX(IBIN,NSUBPROC+1,NORD+1)
                  write(*,*)"wtx,myres,diff",
     >                 WTX(IBIN+NBIN,NSUBPROC+1,NORD+1),
     >                 MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1),diff
                  IF (DIFF.GT.0D0) THEN
                     WTDXU2(IBIN+NBIN,NSUBPROC+1,NORD+1) =
     >                    WTDXU2(IBIN+NBIN,NSUBPROC+1,NORD+1) +
     >                    DIFF*DIFF
                     IF (DIFF.GT.WTDXUM(IBIN+NBIN,NSUBPROC+1,NORD+1))
     >                    THEN
                        WTDXUM(IBIN+NBIN,NSUBPROC+1,NORD+1) = DIFF
                     ENDIF
                  ELSE
                     WTDXL2(IBIN+NBIN,NSUBPROC+1,NORD+1) =
     >                    WTDXL2(IBIN+NBIN,NSUBPROC+1,NORD+1) +
     >                    DIFF*DIFF
                     IF (DIFF.LT.WTDXLM(IBIN+NBIN,NSUBPROC+1,NORD+1))
     >                    THEN
                        WTDXLM(IBIN+NBIN,NSUBPROC+1,NORD+1) = DIFF
                     ENDIF
                  ENDIF
               ENDIF
               DO ISUB=1,NSUBPROC
                  SUMMORD = 0.D0
                  DIFFORD = 0.D0
                  DO IORD=1,NORD
                     SUMM = RESULT(IBIN,ISUB,IORD)
                     SUMMORD = SUMMORD + SUMM
                     DIFF = SUMM - MYRES(IBIN,ISUB,IORD)
                     DIFFORD = DIFFORD + DIFF
                  ENDDO
                  WT(IBIN,ISUB,NORD+1)   =
     >                 WT(IBIN,ISUB,NORD+1)   + 1.D0 
                  WT2(IBIN,ISUB,NORD+1)  =
     >                 WT2(IBIN,ISUB,NORD+1)  + 1.D0 
                  WTX(IBIN,ISUB,NORD+1)  = 
     >                 WTX(IBIN,ISUB,NORD+1)  + SUMMORD
                  WTX2(IBIN,ISUB,NORD+1) =
     >                 WTX2(IBIN,ISUB,NORD+1) + SUMMORD*SUMMORD
                  IF (DIFFORD.GT.0D0) THEN
                     WTDXU2(IBIN,ISUB,NORD+1) =
     >                    WTDXU2(IBIN,ISUB,NORD+1) +
     >                    DIFFORD*DIFFORD
                     IF (DIFFORD.GT.WTDXUM(IBIN,ISUB,NORD+1)) THEN
                        WTDXUM(IBIN,ISUB,NORD+1) = DIFFORD
                     ENDIF
                  ELSE
                     WTDXL2(IBIN,ISUB,NORD+1) =
     >                    WTDXL2(IBIN,ISUB,NORD+1) +
     >                    DIFFORD*DIFFORD
                     IF (DIFFORD.LT.WTDXLM(IBIN,ISUB,NORD+1)) THEN
                        WTDXLM(IBIN,ISUB,NORD+1) = DIFFORD
                     ENDIF
                  ENDIF
                  IF (IMODE.EQ.5) THEN
                     WTX(IBIN+NBIN,ISUB,NORD+1) = 
     >                    WTX(IPT,ISUB,NORD+1) /
     >                    WTX(IBIN,ISUB,NORD+1)
                     DIFF = WTX(IBIN+NBIN,ISUB,NORD+1) -
     >                    MYRES(IBIN+NBIN,ISUB,NORD+1)
                     IF (DIFF.GT.0D0) THEN
                        WTDXU2(IBIN+NBIN,ISUB,NORD+1) =
     >                       WTDXU2(IBIN+NBIN,ISUB,NORD+1) +
     >                       DIFF*DIFF
                        IF (DIFF.GT.WTDXUM(IBIN+NBIN,ISUB,NORD+1))
     >                       THEN
                           WTDXUM(IBIN+NBIN,ISUB,NORD+1) = DIFF
                        ENDIF
                     ELSE
                        WTDXL2(IBIN+NBIN,ISUB,NORD+1) =
     >                       WTDXL2(IBIN+NBIN,ISUB,NORD+1) +
     >                       DIFF*DIFF
                        IF (DIFF.LT.WTDXLM(IBIN+NBIN,ISUB,NORD+1))
     >                       THEN
                           WTDXLM(IBIN+NBIN,ISUB,NORD+1) = DIFF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         RETURN
      ENDIF 
      
      
      
*---Final evaluation
      IF (ISTAT.EQ.3) THEN
         IBIN = 0
         DO IRAP=1,NRAP
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               DO IORD=1,NORD+1
                  DO ISUB=1,NSUBPROC+1
ckr Eigen vector method (as needed for CTEQ or MSTW PDF uncertainty)
                     IF (IMODE.EQ.1.OR.IMODE.EQ.5) THEN
                        WTDXU2(IBIN,ISUB,IORD) =  SQRT(WTDXU2(IBIN
     >                       ,ISUB,IORD))
                        WTDXL2(IBIN,ISUB,IORD) = -SQRT(WTDXL2(IBIN
     >                       ,ISUB,IORD))
                        IF (DABS(MYRES(IBIN,ISUB,IORD)).GT.1.D-99) THEN
                           WTDXU2(IBIN,ISUB,IORD) = 1D0 +
     >                          WTDXU2(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
                           WTDXL2(IBIN,ISUB,IORD) = 1D0 +
     >                          WTDXL2(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
                        ELSE
                           WTDXU2(IBIN,ISUB,IORD) = -1D0
                           WTDXL2(IBIN,ISUB,IORD) = -1D0
                        ENDIF
ckr Sampling method (as needed for NNPDF PDF or stat. uncertainty) 
                     ELSEIF (IMODE.EQ.2) THEN
ckr Attention: This overwrites the precalculated result by averages
ckr            from the sampling method as this should be considered the
ckr            central result according to NNPDF.
                        MYRES(IBIN,ISUB,IORD) = WTX(IBIN,ISUB,IORD)
     >                       /WT(IBIN,ISUB,IORD)
                        WTDXU2(IBIN,ISUB,IORD) = 
     >                       (WTX2(IBIN,ISUB,IORD)/WT(IBIN,ISUB,IORD) -
     >                       MYRES(IBIN,ISUB,IORD)*MYRES(IBIN,ISUB,IORD)
     >                       )
                        WTDXU2(IBIN,ISUB,IORD) =
     >                       SQRT(WTDXU2(IBIN,ISUB,IORD))
                        WTDXL2(IBIN,ISUB,IORD) = -WTDXU2(IBIN,ISUB,IORD)
                        IF (DABS(MYRES(IBIN,ISUB,IORD)).GT.1.D-99) THEN
                           WTDXU2(IBIN,ISUB,IORD) = 1D0 +
     >                          WTDXU2(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
                           WTDXL2(IBIN,ISUB,IORD) = 1D0 +
     >                          WTDXL2(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
                        ELSE
                           WTDXU2(IBIN,ISUB,IORD) = -1D0
                           WTDXL2(IBIN,ISUB,IORD) = -1D0
                        ENDIF
ckr Minimax method (as needed for scale uncertainties) 
                     ELSEIF (IMODE.EQ.3) THEN
                        IF (DABS(MYRES(IBIN,ISUB,IORD)).GT.1.D-99) THEN
                           WTDXUM(IBIN,ISUB,IORD) = 1D0 +
     >                          WTDXUM(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
                           WTDXLM(IBIN,ISUB,IORD) = 1D0 +
     >                          WTDXLM(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
                        ELSE
                           WTDXUM(IBIN,ISUB,IORD) = -1D0
                           WTDXLM(IBIN,ISUB,IORD) = -1D0
                        ENDIF
ckr Deviation from reference (as needed for algorithmic uncertainties) 
                     ELSEIF (IMODE.EQ.4) THEN
                        IF (IRAP.LE.INT(NRAPIDITY/2)) THEN
                           IF (DABS(MYRES(IBIN+NBIN/2,ISUB,IORD)).GT.
     >                          1D-99) THEN
                              WTDXUM(IBIN,ISUB,IORD) =
     >                             WTX(IBIN,ISUB,IORD) /
     >                             MYRES(IBIN+NBIN/2,ISUB,IORD) 
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
         IBIN = 0
         DO IRAP=1,2*NRAPIDITY
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               DO IORD=1,NORD+1
                  DO ISUB=1,NSUBPROC+1
                     write(*,*)"ibin,iord,isub,wtxxx",
     >                    ibin,iord,isub,wtx(IBIN,ISUB,IORD)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         NBIN = IBIN
cdebug

         RETURN
      ENDIF    
      
      RETURN
      END



      SUBROUTINE CENRES(IMODE)
      
      IMPLICIT NONE
      INTEGER IMODE
      INCLUDE "fnx9999.inc"
      INCLUDE "uncert.inc"
      INTEGER IBIN,IRAP,IPT,ISUB,IORD,NBIN

      

*---Initialization	
      IBIN = 0
      DO IRAP=1,NRAPIDITY
         DO IPT=1,NPT(IRAP)
            IBIN = IBIN+1
            DO IORD=1,NORD+1
               DO ISUB=1,NSUBPROC+1
                  MYRES(IBIN,ISUB,IORD) = 0.D0
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      NBIN = IBIN


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
ckr Centrality ratio: IMODE = 5
                  IF (IMODE.EQ.5) THEN
                     MYRES(IBIN+NBIN,ISUB,IORD) = 
     >                    MYRES(IPT,ISUB,IORD) /
     >                    RESULT(IBIN,ISUB,IORD)
                  ENDIF
                  write(*,*)"ibin,irap,ipt,iord,isub",
     >                 ibin,irap,ipt,iord,isub
                  write(*,*)"r,a,b",MYRES(IBIN+NBIN,ISUB,IORD),
     >                 MYRES(IPT,ISUB,IORD),RESULT(IBIN,ISUB,IORD)
               ENDDO
               IF (IMODE.EQ.5) THEN
                  MYRES(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                 MYRES(IPT,NSUBPROC+1,IORD) / 
     >                 MYRES(IBIN,NSUBPROC+1,IORD)
               ENDIF
               write(*,*)"sub r,a,b",MYRES(IBIN+NBIN,NSUBPROC+1,IORD),
     >              MYRES(IPT,NSUBPROC+1,IORD),MYRES(IBIN,NSUBPROC+1
     >              ,IORD)
               MYRES(IBIN,NSUBPROC+1,NORD+1) = 
     >              MYRES(IBIN,NSUBPROC+1,NORD+1) +
     >              MYRES(IBIN,NSUBPROC+1,IORD)
            ENDDO
            IF (IMODE.EQ.5) THEN
               MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1) = 
     >              MYRES(IPT,NSUBPROC+1,NORD+1) /
     >              MYRES(IBIN,NSUBPROC+1,NORD+1)
            ENDIF
            write(*,*)"tot r,a,b",MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1),
     >           MYRES(IPT,NSUBPROC+1,NORD+1),MYRES(IBIN,NSUBPROC+1
     >           ,NORD+1)
            DO ISUB=1,NSUBPROC
               DO IORD=1,NORD
                  MYRES(IBIN,ISUB,NORD+1) =
     >                 MYRES(IBIN,ISUB,NORD+1) + 
     >                 MYRES(IBIN,ISUB,IORD)
               ENDDO
               IF (IMODE.EQ.5) THEN
                  MYRES(IBIN+NBIN,ISUB,NORD+1) = 
     >                 MYRES(IPT,ISUB,NORD+1) /
     >                 MYRES(IBIN,ISUB,NORD+1)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

cdebug
      IBIN = 0
      DO IRAP=1,2*NRAPIDITY
         DO IPT=1,NPT(IRAP)
            IBIN = IBIN+1
            DO IORD=1,NORD+1
               DO ISUB=1,NSUBPROC+1
                  write(*,*)"ibin,iord,isub,myres",
     >                 ibin,iord,isub,MYRES(IBIN,ISUB,IORD)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      NBIN = IBIN
cdebug


      RETURN
      END
