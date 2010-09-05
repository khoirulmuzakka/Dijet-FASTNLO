      SUBROUTINE UNCERT(ISTAT,IMODE,LRAT)
      
      IMPLICIT NONE
      INTEGER ISTAT,IMODE
      LOGICAL LRAT
      INCLUDE "fnx9999.inc"
      INCLUDE "uncert.inc"
      DOUBLE PRECISION DIFF,SUMM,DIFFPROC,SUMMPROC,DIFFORD,SUMMORD
      INTEGER IBIN,IRAP,IPT,IORD,ISUB,NBIN,NRAP

      

*---Initialization	
      IF (ISTAT.EQ.1) THEN
         IBIN = 0
         NRAP = NRAPIDITY
         IF (LRAT) NRAP=2*NRAPIDITY
         DO IRAP=1,NRAP
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               DO IORD=1,NORD+1
                  DO ISUB=1,NSUBPROC+1
                     WT(IBIN,ISUB,IORD)     = 0.D0
                     WT2(IBIN,ISUB,IORD)    = 0.D0
                     WRES(IBIN,ISUB,IORD)   = 0.D0
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
         IF (LRAT) NBIN = NBIN/2
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
                     WRES(IBIN,ISUB,IORD) = SUMM
                     WTX(IBIN,ISUB,IORD) =
     >                    WTX(IBIN,ISUB,IORD)  + SUMM
                     WTX2(IBIN,ISUB,IORD) =
     >                    WTX2(IBIN,ISUB,IORD) + SUMM*SUMM
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
                     IF (LRAT) THEN
                        WTX(IBIN+NBIN,ISUB,IORD) = 
     >                       WRES(IPT,ISUB,IORD) /
     >                       WRES(IBIN,ISUB,IORD)
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
                  WRES(IBIN,NSUBPROC+1,IORD)  = SUMMPROC
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
                  IF (LRAT) THEN
                     WTX(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                    WRES(IPT,NSUBPROC+1,IORD) / 
     >                    WRES(IBIN,NSUBPROC+1,IORD)
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
               WRES(IBIN,NSUBPROC+1,NORD+1) = SUMMORD
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
               IF (LRAT) THEN
                  WTX(IBIN+NBIN,NSUBPROC+1,NORD+1) = 
     >                 WRES(IPT,NSUBPROC+1,NORD+1) /
     >                 WRES(IBIN,NSUBPROC+1,NORD+1)
                  DIFF = WTX(IBIN+NBIN,NSUBPROC+1,NORD+1) -
     >                 MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1)
Comment:                   write(*,*)"COMP: ibin,irap,ipt,iord,isub,i+n",
Comment:      >                 ibin,irap,ipt,iord,isub,ibin+nbin
Comment:                   write(*,*)"COMPA: num,den,wtx,myres,"//
Comment:      >                 "diff,sum^2ul,maxul",
Comment:      >                 WRES(IPT,NSUBPROC+1,NORD+1),
Comment:      >                 WRES(IBIN,NSUBPROC+1,NORD+1),
Comment:      >                 WTX(IBIN+NBIN,NSUBPROC+1,NORD+1),
Comment:      >                 MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1),
Comment:      >                 diff,wtdxu2(ibin+nbin,nsubproc+1,nord+1),
Comment:      >                 wtdxl2(ibin+nbin,nsubproc+1,nord+1),
Comment:      >                 wtdxum(ibin+nbin,nsubproc+1,nord+1),
Comment:      >                 wtdxlm(ibin+nbin,nsubproc+1,nord+1)
                  IF (DIFF.GT.0D0) THEN
                     WTDXU2(IBIN+NBIN,NSUBPROC+1,NORD+1) =
     >                    WTDXU2(IBIN+NBIN,NSUBPROC+1,NORD+1) +
     >                    DIFF*DIFF
Comment:                      write(*,*)"DIFFA diff, wtdxum old",diff, 
Comment:      >                    WTDXUM(IBIN+NBIN,NSUBPROC+1,NORD+1)
                     IF (DIFF.GT.WTDXUM(IBIN+NBIN,NSUBPROC+1,NORD+1))
     >                    THEN
                        WTDXUM(IBIN+NBIN,NSUBPROC+1,NORD+1) = DIFF
Comment:                         write(*,*)"DIFFB diff, wtdxum old",diff, 
Comment:      >                       WTDXUM(IBIN+NBIN,NSUBPROC+1,NORD+1)
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
Comment:                   write(*,*)"COMPB: num,den,wtx,myres,"//
Comment:      >                 "diff,sum^2ul,maxul",
Comment:      >                 WRES(IPT,NSUBPROC+1,NORD+1),
Comment:      >                 WRES(IBIN,NSUBPROC+1,NORD+1),
Comment:      >                 WTX(IBIN+NBIN,NSUBPROC+1,NORD+1),
Comment:      >                 MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1),
Comment:      >                 diff,wtdxu2(ibin+nbin,nsubproc+1,nord+1),
Comment:      >                 wtdxl2(ibin+nbin,nsubproc+1,nord+1),
Comment:      >                 wtdxum(ibin+nbin,nsubproc+1,nord+1),
Comment:      >                 wtdxlm(ibin+nbin,nsubproc+1,nord+1)
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
                  WRES(IBIN,ISUB,NORD+1) = SUMMORD
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
                  IF (LRAT) THEN
                     WTX(IBIN+NBIN,ISUB,NORD+1) = 
     >                    WRES(IPT,ISUB,NORD+1) /
     >                    WRES(IBIN,ISUB,NORD+1)
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
                     IF (IMODE.EQ.1) THEN
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
Comment:                            write(*,*)"FINA: ibin,irap,ipt,iord,isub,i+n"
Comment:      >                          ,ibin,irap,ipt,iord,isub,ibin+nbin
Comment:                            write(*,*)"FINA: dxu, dxl, myres",
Comment:      >                          WTDXUM(IBIN,ISUB,IORD),
Comment:      >                          WTDXLM(IBIN,ISUB,IORD),
Comment:      >                          MYRES(IBIN,ISUB,IORD) 
                           WTDXUM(IBIN,ISUB,IORD) = 1D0 +
     >                          WTDXUM(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
                           WTDXLM(IBIN,ISUB,IORD) = 1D0 +
     >                          WTDXLM(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
Comment:                            write(*,*)"FINB: dxu, dxl, myres",
Comment:      >                          WTDXUM(IBIN,ISUB,IORD),
Comment:      >                          WTDXLM(IBIN,ISUB,IORD),
Comment:      >                          MYRES(IBIN,ISUB,IORD) 
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
Comment:          IBIN = 0
Comment:          DO IRAP=1,NRAPIDITY
Comment:             DO IPT=1,NPT(IRAP)
Comment:                IBIN = IBIN+1
Comment:                DO IORD=1,NORD+1
Comment:                   DO ISUB=1,NSUBPROC+1
Comment:                      write(*,*)"ENDW: ibin,irap,ipt,iord,isub,i+n",
Comment:      >                    ibin,irap,ipt,iord,isub,ibin+nbin
Comment:                      write(*,*)"ENDW: wres,wtx,myres",
Comment:      >                    WRES(IBIN,ISUB,IORD),
Comment:      >                    WTX(IBIN,ISUB,IORD),
Comment:      >                    MYRES(IBIN,ISUB,IORD)
Comment:                      write(*,*)"ENDW: wres,wtx,myresnbin",
Comment:      >                    WRES(IBIN+NBIN,ISUB,IORD),
Comment:      >                    WTX(IBIN+NBIN,ISUB,IORD),
Comment:      >                    MYRES(IBIN+NBIN,ISUB,IORD)
Comment:                   ENDDO
Comment:                ENDDO
Comment:             ENDDO
Comment:          ENDDO
cdebug

         RETURN
      ENDIF    
      
      RETURN
      END



      SUBROUTINE CENRES(LRAT)
      
      IMPLICIT NONE
      LOGICAL LRAT
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
ckr Centrality ratio: LRAT
                  IF (LRAT) THEN
                     MYRES(IBIN+NBIN,ISUB,IORD) = 
     >                    MYRES(IPT,ISUB,IORD) /
     >                    RESULT(IBIN,ISUB,IORD)
                  ENDIF
c                  write(*,*)"ibin,irap,ipt,iord,isub",
c     >                 ibin,irap,ipt,iord,isub
c                  write(*,*)"r,a,b",MYRES(IBIN+NBIN,ISUB,IORD),
c     >                 MYRES(IPT,ISUB,IORD),RESULT(IBIN,ISUB,IORD)
               ENDDO
               IF (LRAT) THEN
                  MYRES(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                 MYRES(IPT,NSUBPROC+1,IORD) / 
     >                 MYRES(IBIN,NSUBPROC+1,IORD)
               ENDIF
c               write(*,*)"sub r,a,b",MYRES(IBIN+NBIN,NSUBPROC+1,IORD),
c     >              MYRES(IPT,NSUBPROC+1,IORD),MYRES(IBIN,NSUBPROC+1
c     >              ,IORD)
               MYRES(IBIN,NSUBPROC+1,NORD+1) = 
     >              MYRES(IBIN,NSUBPROC+1,NORD+1) +
     >              MYRES(IBIN,NSUBPROC+1,IORD)
            ENDDO
            IF (LRAT) THEN
               MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1) = 
     >              MYRES(IPT,NSUBPROC+1,NORD+1) /
     >              MYRES(IBIN,NSUBPROC+1,NORD+1)
            ENDIF
c            write(*,*)"tot r,a,b",MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1),
c     >           MYRES(IPT,NSUBPROC+1,NORD+1),MYRES(IBIN,NSUBPROC+1
c     >           ,NORD+1)
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
               ENDIF
            ENDDO
         ENDDO
      ENDDO

cdebug
Comment:       IBIN = 0
Comment:       DO IRAP=1,NRAPIDITY
Comment:          DO IPT=1,NPT(IRAP)
Comment:             IBIN = IBIN+1
Comment:             DO IORD=1,NORD+1
Comment:                DO ISUB=1,NSUBPROC+1
Comment:                   write(*,*)"ENDM: ibin,irap,ipt,iord,isub,i+n",
Comment:      >                 ibin,irap,ipt,iord,isub,ibin+nbin
Comment:                   write(*,*)"ENDM: myres, myresnbin",
Comment:      >                 MYRES(IBIN,ISUB,IORD),MYRES(IBIN+NBIN,ISUB,IORD)
Comment:                ENDDO
Comment:             ENDDO
Comment:          ENDDO
Comment:       ENDDO
cdebug


      RETURN
      END
