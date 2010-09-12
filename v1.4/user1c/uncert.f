      SUBROUTINE UNCERT(IPHASE,IMODE,IWEIGHT,IVAR,LRAT)
      
      IMPLICIT NONE
      INTEGER IPHASE,IMODE,IWEIGHT,IVAR
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
c - Not active/well tested yet
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
                     WTXMIN(IBIN,ISUB,IORD) = 1D99
                     WTXMAX(IBIN,ISUB,IORD) = 1D-99
                     WTDXL2(IBIN,ISUB,IORD) = 0D0
                     WTDXU2(IBIN,ISUB,IORD) = 0D0
                     WTDXLM(IBIN,ISUB,IORD) = 0D0
                     WTDXUM(IBIN,ISUB,IORD) = 0D0
                     WTDXMN(IBIN,ISUB,IORD) = 0D0
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
               SUMMORD = 0D0
               DIFFORD = 0D0
Comment:                write(*,*)"AAA: wtdxum(1,8,3),so,do,sp,dp,s,d",
Comment:      >              wtdxum(1,8,3),summord,difford,
Comment:      >              summproc,diffproc,summ,diff
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
Comment:                   write(*,*)"BBB: wtdxum(1,8,3),so,do,sp,dp,s,d",
Comment:      >                 wtdxum(1,8,3),summord,difford,
Comment:      >                 summproc,diffproc,summ,diff
                  DO ISUB=1,NSUBPROC
                     SUMM = RESULT(IBIN,ISUB,IORD)
                     SUMMPROC = SUMMPROC + SUMM
                     DIFF = SUMM - MYRES(IBIN,ISUB,IORD)
                     DIFFPROC = DIFFPROC + DIFF
                     CALL SUMMUP(IVAR,IBIN,ISUB,IORD,
     >                    SUMM,DIFF,WEIGHT,
     >                    LRAT,IPT,NBIN)
                  ENDDO
Comment:                   write(*,*)"CCC: wtdxum(1,8,3),so,do,sp,dp,s,d",
Comment:      >                 wtdxum(1,8,3),summord,difford,
Comment:      >                 summproc,diffproc,summ,diff
                  SUMM = SUMMPROC
                  SUMMORD = SUMMORD + SUMMPROC
                  DIFF = SUMM - MYRES(IBIN,NSUBPROC+1,IORD)
                  DIFFORD = DIFFORD + DIFF
c                  write(*,*)"XXX: DIFFPROC1,2",DIFFPROC,
c     >                 SUMM-myres(ibin,nsubproc+1,iord)
                  CALL SUMMUP(IVAR,IBIN,NSUBPROC+1,IORD,
     >                 SUMMPROC,DIFFPROC,WEIGHT,
     >                 LRAT,IPT,NBIN)
               ENDDO
Comment:                write(*,*)"DDD: wtdxum(1,8,3),so,do,sp,dp,s,d",
Comment:      >              wtdxum(1,8,3),summord,difford,
Comment:      >              summproc,diffproc,summ,diff
c               write(*,*)"YYY: DIFFORD1,2",DIFFORD,
c     >              SUMMORD-myres(ibin,nsubproc+1,nord+1)
               CALL SUMMUP(IVAR,IBIN,NSUBPROC+1,NORD+1,
     >              SUMMORD,DIFFORD,WEIGHT,
     >              LRAT,IPT,NBIN)
Comment:                if (ibin.eq.1) then
Comment:                   write(*,*)"UNCERT: ibin,itab,difford",
Comment:      >                 ibin,ivar,difford
Comment:                   write(*,*)"UNCERT:wt,wt2,wtx,wtx2,s+wtdxlm,s+wtdxum"//
Comment:      >                 "wtxmin,wtxmax",
Comment:      >                 WT(IBIN,NSUBPROC+1,NORD+1),
Comment:      >                 Wt2(IBIN,NSUBPROC+1,NORD+1),
Comment:      >                 Wtx(IBIN,NSUBPROC+1,nORD+1),
Comment:      >                 WTx2(IBIN,NSUBPROC+1,nORD+1),
Comment:      >                 Wtx(IBIN,NSUBPROC+1,nORD+1) /
Comment:      >                 WT(IBIN,NSUBPROC+1,NORD+1) +
Comment:      >                 WTdxlm(IBIN,NSUBPROC+1,nORD+1),
Comment:      >                 Wtx(IBIN,NSUBPROC+1,nORD+1) /
Comment:      >                 WT(IBIN,NSUBPROC+1,NORD+1) +
Comment:      >                 WTdxum(IBIN,NSUBPROC+1,nORD+1),
Comment:      >                 WTxmin(IBIN,NSUBPROC+1,nORD+1),
Comment:      >                 WTxmax(IBIN,NSUBPROC+1,nORD+1)
Comment:                endif
               DO ISUB=1,NSUBPROC
                  SUMMORD = 0D0
                  DIFFORD = 0D0
                  DO IORD=1,NORD
                     SUMM = RESULT(IBIN,ISUB,IORD)
                     SUMMORD = SUMMORD + SUMM
                     DIFF = SUMM - MYRES(IBIN,ISUB,IORD)
                     DIFFORD = DIFFORD + DIFF
                  ENDDO
                  CALL SUMMUP(IVAR,IBIN,ISUB,NORD+1,
     >                 SUMMORD,DIFFORD,WEIGHT,
     >                 LRAT,IPT,NBIN)
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
                        IF (DABS(MYRES(IBIN,ISUB,IORD)).GT.1D-99) THEN
                           WTDXU2(IBIN,ISUB,IORD) =
     >                          WTDXU2(IBIN,ISUB,IORD) /
     >                          DABS(MYRES(IBIN,ISUB,IORD))
                           WTDXL2(IBIN,ISUB,IORD) =
     >                          WTDXL2(IBIN,ISUB,IORD) /
     >                          DABS(MYRES(IBIN,ISUB,IORD))
                        ELSE
                           WTDXU2(IBIN,ISUB,IORD) = -1D0
                           WTDXL2(IBIN,ISUB,IORD) = -1D0
                        ENDIF
ckr Sampling method (as needed for statistical uncertainties or NNPDF) 
                     ELSEIF (IMODE.EQ.2) THEN
                        WTDXU2(IBIN,ISUB,IORD) = -1D0
                        WTDXL2(IBIN,ISUB,IORD) = -1D0
                        WTDXMN(IBIN,ISUB,IORD) = -1D0
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
                              ELSE
                                 WTXMIN(IBIN,ISUB,IORD) = -1D0
                                 WTXMAX(IBIN,ISUB,IORD) = -1D0
                                 WTDXU2(IBIN,ISUB,IORD) = -1D0
                                 WTDXL2(IBIN,ISUB,IORD) = -1D0
                                 WTDXMN(IBIN,ISUB,IORD) = -1D0
                              ENDIF
                           ELSE
                              WTDXU2(IBIN,ISUB,IORD) = -1D0
                              WTDXL2(IBIN,ISUB,IORD) = -1D0
                              WTDXMN(IBIN,ISUB,IORD) = -1D0
                           ENDIF
                        ENDIF
ckr Minimax method (as needed for scale uncertainties) 
                     ELSEIF (IMODE.EQ.3) THEN
                        IF (DABS(MYRES(IBIN,ISUB,IORD)).GT.1D-99) THEN
                           WTDXUM(IBIN,ISUB,IORD) =
     >                          WTDXUM(IBIN,ISUB,IORD) /
     >                          DABS(MYRES(IBIN,ISUB,IORD))
                           WTDXLM(IBIN,ISUB,IORD) =
     >                          WTDXLM(IBIN,ISUB,IORD) /
     >                          DABS(MYRES(IBIN,ISUB,IORD))
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
     >                             DABS(MYRES(IBIN+NBIN/2,ISUB,IORD))) -
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
                  MYRES(IBIN,ISUB,IORD) = 0D0
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
               ENDDO
               IF (LRAT) THEN
                  MYRES(IBIN+NBIN,NSUBPROC+1,IORD) =
     >                 MYRES(IPT,NSUBPROC+1,IORD) / 
     >                 MYRES(IBIN,NSUBPROC+1,IORD)
               ENDIF
               MYRES(IBIN,NSUBPROC+1,NORD+1) = 
     >              MYRES(IBIN,NSUBPROC+1,NORD+1) +
     >              MYRES(IBIN,NSUBPROC+1,IORD)
            ENDDO
            IF (LRAT) THEN
               MYRES(IBIN+NBIN,NSUBPROC+1,NORD+1) = 
     >              MYRES(IPT,NSUBPROC+1,NORD+1) /
     >              MYRES(IBIN,NSUBPROC+1,NORD+1)
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
      
      
      
      SUBROUTINE SUMMUP(IVAR,IBIN,ISUB,IORD,SUMM,DIFF,WEIGHT,
     >     LRAT,IBREF,NBIN)
      
      IMPLICIT NONE
      INTEGER IVAR,IBIN,IORD,ISUB,IBREF,NBIN
      DOUBLE PRECISION SUMM,DIFF,WEIGHT
      LOGICAL LRAT
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
      
      RETURN
      END
