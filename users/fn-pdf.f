*******************************************************************
*  M. Wobisch - July 26, 2005
*
*  code to compute PDF linear combinations for different subprocesses
*
* included:
*      SUBROUTINE FNPAPDF(X1,X2,MUF,PDF)   PA for proton-antiproton
*      SUBROUTINE FNPPPDF(X1,X2,MUF,PDF)   PP for proton-proton
*      SUBROUTINE FNEPPDF(X,MUF,PDF)       EP for Electron-Proton
*
*  all routines access FNPDF(X,MUF,XPDF) to provide PDF information
*
*******************************************************************
      SUBROUTINE FNPAPDF(X1,X2,MUF,PDF)
*-----------------------------------------------------------------
* MW 05/12/2005  
*
* interface to PDFs from LHAPDF - for proton-antiproton collisions
* input   X1,X2   values for both hadrons
*         MUF     factorization scale in GeV
* output  PDF(7)  array with seven PDF combinations H1-H7
*                     note: PDFs are *not* multiplied with x
*-----------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION X1,X2,MUF,PDF(7)
      INTEGER i
      DOUBLE PRECISION 
     +     G1, G2,              ! gluon densities from both hadrons
     +     SumQ1, SumQ2,        ! sum of quark densities
     +     SumQB1, SumQB2,      ! sum of anti-quark densities
     +     Q1(6),Q2(6), QB1(6),QB2(6), ! arrays of 6 (anti-)quark densities
     +     S,A                  ! products S,A

c - define variables for PDF cache
      DOUBLE PRECISION XCACHE,MUFCACHE,DIFF
      PARAMETER (DIFF=1D-10)
      SAVE XCACHE,MUFCACHE,G1,SumQ1,SumQB1,Q1,QB1
      DATA XCACHE, MUFCACHE/2*0d0/
      DATA SumQ1,SumQB1,Q1,QB1,G1 /15*0d0/

c - define pdf arrays for LHAPDF
      DOUBLE PRECISION XPDF1(-6:6), XPDF2(-6:6) 

c - get PDFs for first Hadron (always Proton)
c         apply the strategy as used in DISENT:
c          - make a cache for X1,MUF (the PDFs are SAVE-ed)
c          - first: check if X1 and MUF are the same as in prev. call
      if (abs(X1-XCACHE).gt.DIFF .or. abs(MUF-MUFCACHE).gt.DIFF) then

c  - note: the momentum densities are returned, i.e. x X pdf! 
c  - tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t
c  -  -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 , 0, 1, 2, 3, 4, 5, 6

         call FNPDF(X1,MUF,XPDF1)
         SumQ1  = 0d0
         SumQB1 = 0d0
         do i=1,6
            Q1(i)  = XPDF1(i)
            QB1(i) = XPDF1(-i)
            SumQ1  = SumQ1  + Q1(i)
            SumQB1 = SumQB1 + QB1(i)
         enddo
         G1     = XPDF1(0)
      endif

c - get PDFs for second Hadron    (no cache needed here)
      call FNPDF(X2,MUF,XPDF2)
      SumQ2  = 0d0
      SumQB2 = 0d0
      do i=1,6
         Q2(i)  = XPDF2(i)
         QB2(i) = XPDF2(-i)
         SumQ2  = SumQ2  + Q2(i)
         SumQB2 = SumQB2 + QB2(i)
      enddo
      G2     = XPDF2(0)
        
c - compute S,A
      S = 0d0
      A = 0d0
      do i=1,6
         S = S + (Q1(i)*Q2(i)) + (QB1(i)*QB2(i)) 
         A = A + (Q1(i)*QB2(i)) + (QB1(i)*Q2(i)) 
      enddo

c - compute seven combinatins  (for pp)
      PDF(1) = G1*G2
      PDF(2) = (SumQ1+SumQB1)*G2
      PDF(3) = G1*(SumQ2+SumQB2)
c      PDF(4) = SumQ1*SumQ2 + SumQB1*SumQB2 - S
c      PDF(5) = S
c      PDF(6) = A
c      PDF(7) = SumQ1*SumQB2 + SumQB1*SumQ2 - A

c   - for p-pbar: swap 4<->7 and 5<->6
      PDF(7) = SumQ1*SumQ2 + SumQB1*SumQB2 - S
      PDF(6) = S
      PDF(5) = A
      PDF(4) = SumQ1*SumQB2 + SumQB1*SumQ2 - A

c - test phase: for tables that already include PDF info 
c      do i=1,7
c         PDF(i)=1d0
c      enddo

      RETURN
      END
*******************************************************************
      SUBROUTINE FNPPPDF(X1,X2,MUF,PDF)
*-----------------------------------------------------------------
* MW 07/25/2005  
*
* interface to PDFs from LHAPDF - for proton-proton collisions
* input   X1,X2   values for both hadrons
*         MUF     factorization scale in GeV
* output  PDF(7)  array with seven PDF combinations H1-H7
*                     note: PDFs are *not* multiplied with x
*-----------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION X1,X2,MUF,PDF(7)
      INTEGER i
      DOUBLE PRECISION 
     +     G1, G2,              ! gluon densities from both hadrons
     +     SumQ1, SumQ2,        ! sum of quark densities
     +     SumQB1, SumQB2,      ! sum of anti-quark densities
     +     Q1(6),Q2(6), QB1(6),QB2(6), ! arrays of 6 (anti-)quark densities
     +     S,A                  ! products S,A

c - define variables for PDF cache
      DOUBLE PRECISION XCACHE,MUFCACHE,DIFF
      PARAMETER (DIFF=1D-10)
      SAVE XCACHE,MUFCACHE,G1,SumQ1,SumQB1,Q1,QB1
      DATA XCACHE, MUFCACHE/2*0d0/
      DATA SumQ1,SumQB1,Q1,QB1,G1 /15*0d0/

c - define pdf arrays for LHAPDF
      DOUBLE PRECISION XPDF1(-6:6), XPDF2(-6:6) 

c - get PDFs for first Hadron (always Proton)
c         apply the strategy as used in DISENT:
c          - make a cache for X1,MUF (the PDFs are SAVE-ed)
c          - first: check if X1 and MUF are the same as in prev. call
      if (abs(X1-XCACHE).gt.DIFF .or. abs(MUF-MUFCACHE).gt.DIFF) then

c  - note: the momentum densities are returned, i.e. x X pdf! 
c  - tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t
c  -  -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 , 0, 1, 2, 3, 4, 5, 6

         call FNPDF(X1,MUF,XPDF1)
         SumQ1  = 0d0
         SumQB1 = 0d0
         do i=1,6
            Q1(i)  = XPDF1(i)
            QB1(i) = XPDF1(-i)
            SumQ1  = SumQ1  + Q1(i)
            SumQB1 = SumQB1 + QB1(i)
         enddo
         G1     = XPDF1(0)
      endif

c - get PDFs for second Hadron    (no cache needed here)
      call FNPDF(X2,MUF,XPDF2)
      SumQ2  = 0d0
      SumQB2 = 0d0
      do i=1,6
         Q2(i)  = XPDF2(i)
         QB2(i) = XPDF2(-i)
         SumQ2  = SumQ2  + Q2(i)
         SumQB2 = SumQB2 + QB2(i)
      enddo
      G2     = XPDF2(0)
        
c - compute S,A
      S = 0d0
      A = 0d0
      do i=1,6
         S = S + (Q1(i)*Q2(i)) + (QB1(i)*QB2(i)) 
         A = A + (Q1(i)*QB2(i)) + (QB1(i)*Q2(i)) 
      enddo

c - compute seven combinatins  (for pp)
      PDF(1) = G1*G2
      PDF(2) = (SumQ1+SumQB1)*G2
      PDF(3) = G1*(SumQ2+SumQB2)
      PDF(4) = SumQ1*SumQ2 + SumQB1*SumQB2 - S
      PDF(5) = S
      PDF(6) = A
      PDF(7) = SumQ1*SumQB2 + SumQB1*SumQ2 - A

      RETURN
      END
*******************************************************************
      SUBROUTINE FNEPPDF(X,MUF,PDF)
*-----------------------------------------------------------------
* MW 07/25/2005  
*
* interface to PDFs from LHAPDF - for electron-proton collisions
* input   X       momentum fraction of hadron
*         MUF     factorization scale in GeV
* output  PDF(3)  array with three PDF combinations
*                          D1 gluon
*                          D2 sigma
*                          D3 delta
*                     note: PDFs are *not* multiplied with x
*-----------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION X,MUF,PDF(3)
      INTEGER i
c - define pdf array for LHAPDF
      DOUBLE PRECISION XPDF(-6:6)

c  - note: the momentum densities are returned, i.e. x X pdf! 
c  - tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t
c  -  -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 , 0, 1, 2, 3, 4, 5, 6
     
      call FNPDF(X,MUF,XPDF)
      PDF(1) = XPDF(0)          ! gluon
      PDF(2) = 0d0              ! Sigma 
      PDF(3) = 0d0              ! Delta
      do i=1,6
         PDF(2) = PDF(2)+ XPDF(i)+XPDF(-i)
      enddo
      do i=1,5,2
         PDF(3) = PDF(3)+ (XPDF(i)+XPDF(-i)+
     +        4d0*(XPDF(i+1)+XPDF(-i-1)))/9d0
      enddo

      RETURN
      END
