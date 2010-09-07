      SUBROUTINE UNCERT(IPHASE,IMODE,IWEIGHT,LRAT)
      
      IMPLICIT NONE
      INTEGER IPHASE,IMODE,IWEIGHT
      LOGICAL LRAT
      INCLUDE "fnx9999.inc"
      INCLUDE "uncert.inc"
      DOUBLE PRECISION DIFF,SUMM,RATIO
      DOUBLE PRECISION DIFFPROC,SUMMPROC,DIFFORD,SUMMORD
      DOUBLE PRECISION TWGT(3),TEVTS,WEIGHT,NEFF
      INTEGER IBIN,IRAP,IPT,IORD,ISUB,NBIN,NRAP
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST/.TRUE./

*---First run
c - Implement sample averages with different weights
c - Required in case of e.g. tables with changing nos. of events
c - Define weights to be one for 1D9 events for LO and 1D8 for NLO/total 
c - If all weights are identical this has to cancel out completely
      IF (FIRST) THEN
         FIRST = .FALSE.
c - This is the order counting in fnx9999 common ...
         IF (NORD.GT.2) THEN
            WRITE(*,*)"UNCERT: ERROR! Orders higher than ",NORD,
     >           " are currently not supported with LSTAT, stopped!"
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
         IF (LRAT) NRAP=2*NRAPIDITY
         DO IRAP=1,NRAP
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               DO IORD=1,NORD+1
                  DO ISUB=1,NSUBPROC+1
                     WRES(IBIN,ISUB,IORD)   = 0.D0
                     WT(IBIN,ISUB,IORD)     = 0.D0
                     WT2(IBIN,ISUB,IORD)    = 0.D0
                     WTN(IBIN,ISUB,IORD)    = 0.D0
                     WTX(IBIN,ISUB,IORD)    = 0.D0
 1                   WTX2(IBIN,ISUB,IORD)   = 0.D0
                     WTDXL2(IBIN,ISUB,IORD) = 0.D0
                     WTDXU2(IBIN,ISUB,IORD) = 0.D0
                     WTDXLM(IBIN,ISUB,IORD) = 0.D0
                     WTDXUM(IBIN,ISUB,IORD) = 0.D0
                     WTDXMN(IBIN,ISUB,IORD) = 0.D0
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         NBIN = IBIN
         IF (LRAT) NBIN = NBIN/2
         RETURN
      ENDIF
      
        
        
*---Loop
      IF (IPHASE.EQ.2) THEN
         IBIN = 0
         DO IRAP=1,NRAPIDITY
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               SUMMORD = 0.D0
               DIFFORD = 0.D0
               DO IORD=1,NORD
                  IF (IORD.LT.NORD) THEN
                     TEVTS = DBLE(NEVT(IORD))
                  ELSE
                     TEVTS = DBLE(NEVT(2))
                  ENDIF
                  WEIGHT = 1D0
                  IF (IWEIGHT.EQ.1) WEIGHT = TWGT(IORD)*TEVTS
Comment:                   write(*,*)"iord,tevts,twgt,weight",
Comment:      >                 iord,tevts,twgt(iord),weight
                  SUMMPROC = 0.D0
                  DIFFPROC = 0.D0
                  DO ISUB=1,NSUBPROC
                     SUMM = RESULT(IBIN,ISUB,IORD)
                     SUMMPROC = SUMMPROC + SUMM
                     DIFF = SUMM - MYRES(IBIN,ISUB,IORD)
                     DIFFPROC = DIFFPROC + DIFF
                     WT(IBIN,ISUB,IORD)   =
     >                    WT(IBIN,ISUB,IORD)  +
     >                    WEIGHT
                     WT2(IBIN,ISUB,IORD)  =
     >                    WT2(IBIN,ISUB,IORD) +
     >                    WEIGHT*WEIGHT 
                     WTN(IBIN,ISUB,IORD)  =
     >                    WTN(IBIN,ISUB,IORD) +
     >                    1.D0 
                     WRES(IBIN,ISUB,IORD) = SUMM
                     WTX(IBIN,ISUB,IORD) =
     >                    WTX(IBIN,ISUB,IORD)  +
     >                    SUMM*WEIGHT
                     WTX2(IBIN,ISUB,IORD) =
     >                    WTX2(IBIN,ISUB,IORD) +
     >                    SUMM*WEIGHT*SUMM*WEIGHT
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
                        RATIO =
     >                       WRES(IPT,ISUB,IORD) /
     >                       WRES(IBIN,ISUB,IORD)
                        WRES(IBIN+NBIN,ISUB,IORD) = RATIO
                        WT(IBIN+NBIN,ISUB,IORD) =
     >                       WT(IBIN+NBIN,ISUB,IORD) + WEIGHT 
                        WT2(IBIN+NBIN,ISUB,IORD) =
     >                       WT2(IBIN+NBIN,ISUB,IORD) + WEIGHT*WEIGHT 
                        WTX(IBIN+NBIN,ISUB,IORD) = 
     >                       WTX(IBIN+NBIN,ISUB,IORD) +
     >                       RATIO*WEIGHT
                        WTX2(IBIN+NBIN,ISUB,IORD) = 
     >                       WTX2(IBIN+NBIN,ISUB,IORD) +
     >                       RATIO*WEIGHT*RATIO*WEIGHT
                        DIFF = RATIO -
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
     >                 WT(IBIN,NSUBPROC+1,IORD)  +
     >                 WEIGHT
                  WT2(IBIN,NSUBPROC+1,IORD)   =
     >                 WT2(IBIN,NSUBPROC+1,IORD) +
     >                 WEIGHT*WEIGHT 
                  WTN(IBIN,NSUBPROC+1,IORD)   =
     >                 WTN(IBIN,NSUBPROC+1,IORD) +
     >                 1D0
                  WRES(IBIN,NSUBPROC+1,IORD)  = SUMMPROC
                  WTX(IBIN,NSUBPROC+1,IORD)   =
     >                 WTX(IBIN,NSUBPROC+1,IORD)  +
     >                 SUMMPROC*WEIGHT
                  WTX2(IBIN,NSUBPROC+1,IORD)  =
     >                 WTX2(IBIN,NSUBPROC+1,IORD) +
     >                 SUMMPROC*WEIGHT*SUMMPROC*WEIGHT
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
                     RATIO =
     >                    WRES(IPT,NSUBPROC+1,IORD) / 
     >                    WRES(IBIN,NSUBPROC+1,IORD)
                     WRES(IBIN+NBIN,NSUBPROC+1,IORD) = RATIO
                     WT(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                    WT(IBIN+NBIN,NSUBPROC+1,IORD) + WEIGHT 
                     WT2(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                    WT2(IBIN+NBIN,NSUBPROC+1,IORD) + WEIGHT*WEIGHT
                     WTX(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                    WTX(IBIN+NBIN,NSUBPROC+1,IORD) +
     >                    RATIO*WEIGHT
                     WTX2(IBIN+NBIN,NSUBPROC+1,IORD) = 
     >                    WTX2(IBIN+NBIN,NSUBPROC+1,IORD) +
     >                    RATIO*WEIGHT*RATIO*WEIGHT
                     DIFF = RATIO -
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
     >              WT(IBIN,NSUBPROC+1,NORD+1)  +
     >              WEIGHT
               WT2(IBIN,NSUBPROC+1,NORD+1)  =
     >              WT2(IBIN,NSUBPROC+1,NORD+1) +
     >              WEIGHT*WEIGHT 
               WTN(IBIN,NSUBPROC+1,NORD+1)  =
     >              WTN(IBIN,NSUBPROC+1,NORD+1) +
     >              1D0
               WRES(IBIN,NSUBPROC+1,NORD+1) = SUMMORD
               WTX(IBIN,NSUBPROC+1,NORD+1)  = 
     >              WTX(IBIN,NSUBPROC+1,NORD+1)  +
     >              SUMMORD*WEIGHT
               WTX2(IBIN,NSUBPROC+1,NORD+1) =
     >              WTX2(IBIN,NSUBPROC+1,NORD+1) +
     >              SUMMORD*WEIGHT*SUMMORD*WEIGHT
Comment:                write(*,*)"MOD1: ibin,irap,ipt,iord,isub,i+n"
Comment:      >              ,ibin,irap,ipt,iord,isub,ibin+nbin
Comment:                write(*,*)"wt,summord,wtx,wtx2",
Comment:      >              WT(IBIN,NSUBPROC+1,NORD+1),
Comment:      >              summord,
Comment:      >              WTX(IBIN,NSUBPROC+1,NORD+1),
Comment:      >              WTX2(IBIN,NSUBPROC+1,NORD+1)
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
                  RATIO =
     >                 WRES(IPT,NSUBPROC+1,NORD+1) / 
     >                 WRES(IBIN,NSUBPROC+1,NORD+1)
                  WRES(IBIN+NBIN,NSUBPROC+1,NORD+1) = RATIO
                  WT(IBIN+NBIN,NSUBPROC+1,NORD+1) =
     >                 WT(IBIN+NBIN,NSUBPROC+1,NORD+1) + WEIGHT 
                  WT2(IBIN+NBIN,NSUBPROC+1,NORD+1) =
     >                 WT2(IBIN+NBIN,NSUBPROC+1,NORD+1) + WEIGHT*WEIGHT 
                  WTX(IBIN+NBIN,NSUBPROC+1,NORD+1) =
     >                 WTX(IBIN+NBIN,NSUBPROC+1,NORD+1) +
     >                 RATIO*WEIGHT
                  WTX2(IBIN+NBIN,NSUBPROC+1,NORD+1) = 
     >                 WTX2(IBIN+NBIN,NSUBPROC+1,NORD+1) +
     >                 RATIO*WEIGHT*RATIO*WEIGHT
                  DIFF = RATIO -
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
     >                 WT(IBIN,ISUB,NORD+1)  +
     >                 WEIGHT
                  WT2(IBIN,ISUB,NORD+1)  =
     >                 WT2(IBIN,ISUB,NORD+1) +
     >                 WEIGHT*WEIGHT 
                  WTN(IBIN,ISUB,NORD+1)  =
     >                 WTN(IBIN,ISUB,NORD+1) +
     >                 1D0
                  WRES(IBIN,ISUB,NORD+1) = SUMMORD
                  WTX(IBIN,ISUB,NORD+1)  = 
     >                 WTX(IBIN,ISUB,NORD+1)  +
     >                 SUMMORD*WEIGHT
                  WTX2(IBIN,ISUB,NORD+1) =
     >                 WTX2(IBIN,ISUB,NORD+1) +
     >                 SUMMORD*WEIGHT*SUMMORD*WEIGHT
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
                     RATIO =
     >                    WRES(IPT,ISUB,NORD+1) / 
     >                    WRES(IBIN,ISUB,NORD+1)
                     WRES(IBIN+NBIN,ISUB,NORD+1) = RATIO
                     WT(IBIN+NBIN,ISUB,NORD+1) =
     >                    WT(IBIN+NBIN,ISUB,NORD+1) + WEIGHT 
                     WT2(IBIN+NBIN,ISUB,NORD+1) =
     >                    WT2(IBIN+NBIN,ISUB,NORD+1) + WEIGHT*WEIGHT 
                     WTX(IBIN+NBIN,ISUB,NORD+1) =
     >                    WTX(IBIN+NBIN,ISUB,NORD+1) +
     >                    RATIO*WEIGHT
                     WTX2(IBIN+NBIN,ISUB,NORD+1) = 
     >                    WTX2(IBIN+NBIN,ISUB,NORD+1) +
     >                    RATIO*WEIGHT*RATIO*WEIGHT
                     DIFF = RATIO -
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
                        IF (DABS(MYRES(IBIN,ISUB,IORD)).GT.1.D-99) THEN
                           WTDXU2(IBIN,ISUB,IORD) =
     >                          WTDXU2(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
                           WTDXL2(IBIN,ISUB,IORD) =
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
Comment:                         write(*,*)"MOD2: ibin,irap,ipt,iord,isub,i+n"
Comment:      >                       ,ibin,irap,ipt,iord,isub,ibin+nbin
                        IF (WT(IBIN,ISUB,IORD).GT.1D-99) THEN
                           MYRES(IBIN,ISUB,IORD) = 
     >                          WTX(IBIN,ISUB,IORD)/WT(IBIN,ISUB,IORD)
                           WTDXU2(IBIN,ISUB,IORD) = 
     >                          (WTX2(IBIN,ISUB,IORD)/WT(IBIN,ISUB,IORD)
     >                          -MYRES(IBIN,ISUB,IORD)
     >                          *MYRES(IBIN,ISUB,IORD))
                        ELSE
                           WTDXU2(IBIN,ISUB,IORD) = -1D0 
                        ENDIF
Comment:                         if (iord.eq.3.and.isub.eq.8.and.ibin.eq.151)then
Comment:                            write(*,*)"MMMM: wt, wtx, wtx2, wtxu2 A",
Comment:      >                          WT(IBIN,ISUB,IORD), WTX(IBIN,ISUB,IORD),
Comment:      >                          WTX2(IBIN,ISUB,IORD),WTDXU2(IBIN,ISUB
Comment:      >                          ,IORD)
Comment:                            write(*,*)"MMMM: wtx2/wt",
Comment:      >                          WTX2(IBIN,ISUB,IORD)/WT(IBIN,ISUB,IORD)
Comment:                            write(*,*)"MMMM: myres^2",
Comment:      >                          MYRES(IBIN,ISUB,IORD)*MYRES(IBIN,ISUB
Comment:      >                          ,IORD)
Comment:                         endif
                        IF (WTDXU2(IBIN,ISUB,IORD).GE.0D0) THEN
                           WTDXU2(IBIN,ISUB,IORD) =
     >                          SQRT(WTDXU2(IBIN,ISUB,IORD))
                           WTDXL2(IBIN,ISUB,IORD) =
     >                          -WTDXU2(IBIN,ISUB,IORD)
                           NEFF = WT(IBIN,ISUB,IORD)*WT(IBIN,ISUB
     >                          ,IORD)/WT2(IBIN,ISUB,IORD)
Comment:                         if (iord.eq.3.and.isub.eq.8.and.ibin.eq.151)then
Comment:                            write(*,*)"NNNN: wt, wt2, neff",
Comment:      >                          WT(IBIN,ISUB,IORD),
Comment:      >                          WT2(IBIN,ISUB,IORD),neff
Comment:                         endif
                           IF (NEFF.GE.2D0) THEN
                              WTDXMN(IBIN,ISUB,IORD) =
     >                             WTDXU2(IBIN,ISUB,IORD) /
     >                             SQRT(NEFF-1.D0)
                           ELSE
                              WTDXMN(IBIN,ISUB,IORD) = -1D0
                           ENDIF
                        ELSE
                           WTDXU2(IBIN,ISUB,IORD) = -1D0
                           WTDXL2(IBIN,ISUB,IORD) = -1D0
                           WTDXMN(IBIN,ISUB,IORD) = -1D0
                        ENDIF
Comment:                         if (iord.eq.3.and.isub.eq.8.and.ibin.eq.151)then
Comment:                            write(*,*)"TTTT: myres, wres, wtdxmn",
Comment:      >                          myres(ibin,isub,iord),
Comment:      >                          wres(ibin,isub,iord),
Comment:      >                          wtdxmn(ibin,isub,iord)
Comment:                         endif
                        IF (DABS(MYRES(IBIN,ISUB,IORD)).GT.1.D-99)
     >                       THEN
                           WTDXU2(IBIN,ISUB,IORD) =
     >                          WTDXU2(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
                           WTDXL2(IBIN,ISUB,IORD) =
     >                          WTDXL2(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
                           WTDXMN(IBIN,ISUB,IORD) =
     >                          WTDXMN(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
                        ELSE
                           WTDXU2(IBIN,ISUB,IORD) = -1D0
                           WTDXL2(IBIN,ISUB,IORD) = -1D0
                           WTDXMN(IBIN,ISUB,IORD) = -1D0
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
                           WTDXUM(IBIN,ISUB,IORD) =
     >                          WTDXUM(IBIN,ISUB,IORD) /
     >                          MYRES(IBIN,ISUB,IORD) 
                           WTDXLM(IBIN,ISUB,IORD) =
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
     >                             (WTX(IBIN,ISUB,IORD) /
     >                             MYRES(IBIN+NBIN/2,ISUB,IORD)) - 1D0
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
         DO IRAP=1,NRAPIDITY
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               DO IORD=1,NORD+1
                  DO ISUB=1,NSUBPROC+1
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
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
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
