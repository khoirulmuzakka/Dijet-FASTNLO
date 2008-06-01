***********************************************************************
***********************************************************************
* fastNLO user code                 T. Kluge, M. Wobisch 02/01/2006
* Cleaned up and improved version   K. Rabbertz          01/06/2008
*                            
*   >> this code should not be edited
*   >> if you have questions, please contact 
*      the authors at:  fastnlo@cedar.ac.uk
*
* fnxNNNN.f     
*     -> where "Xnnnn" is the respective fastNLO scenario name
*            "X" is an abbreviation for the collider
*                   "L" LHC
*                   "T" Tevatron
*                   "H" HERA 
*                   "R" RHIC
*            "nnnn" is a scenario number
*
* contains the following routines
*     fXnnnnCC        computes cross sections / fills arrays
*     fXnnnnGP        gets the PDF values for all x-bins
*     fXnnnnPL        computes PDF linear combinations
*     fXnnnnMT        multiply PDF array and coefficient array
*     fXnnnnPR        prints x-section results 
*     fXnnnnRD        reads coefficient tables
*
* uses commonblock definitions in: fnXNNNN.inc
*
* needs the following routines from "fn-interface.f":
*     FNALPHAS        alpha_s interface (double precision function)
*     FNPDF           PDF interface
*
***********************************************************************
***********************************************************************



***********************************************************************
      SUBROUTINE FX9999CC(FILENAME, XMUR, XMUF, IPRINTFLAG, XSECT)
*----------------------------------------------------------------------
* fastNLO user code - main routine 
*
* input:
*   FILENAME    name of input table
*   XMUR        prefactor for nominal renormalization scale     
*                    any choice is possible, but please note 
*                    that 2-loop threshold corrections work
*                    only for xmur=xmuf
*   XMUF        prefactor for nominal fact-scale
*                     only a few choices are possible
*                     (see output or table documentation)
*   IPRINTFALG  =1 print results / =0 don't print
*
* output:
*   XSECT(nbin,3)    array of cross sections 
*                      first dimension: Bin number
*                      second dimension: 1 LO   2 NLO   
*                        3 NNLO / or other correction (if available)
*
* MW 04/15/2005 initial version
* MW 09/02/2005 implement flexible scale variations
* TK 12/07/2005 table format now LO + NLO with 5 scale variations
* MW 2006/01/17 implement tableformat version 1c
* MW 2006/02/01 implement tableformat version 1.4 
* TK 2006/08/06 implement scalebins for N>2
* MW 2006/08/09 add normalization feature: 1/sigma dsigma/d[s.th.] 
*-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
ckr 30.01.2008: Use LENOCC string length function 
      INTEGER IFIRST,IORD,ISUB,LENOCC,I,J,K,L,M, 
     >     IPRINTFLAG,INORMFLAG, 
     >     NBIN,IXMUR,IXMUF
      CHARACTER*(*) FILENAME
ckr 01.06.2008 Allow longer paths/filenames
      CHARACTER*255 OLDFILENAME
ckr 30.01.2008: Define temp. CSTRNG for output
      CHARACTER*255 CSTRNG
      DOUBLE PRECISION XMUR,XMUF,SUM(3)
      SAVE IFIRST, OLDFILENAME
      INCLUDE 'strings.inc'
      DATA IFIRST/0/, OLDFILENAME/'xxxx'/

c - Output in first call
      IF (IFIRST.EQ.0) THEN
         DO I=1,12
ckr 30.01.2008: Print strings only until end, not complete
ckr             CHARACTER dimension
ckr             Similar usage below with CSTRNG at many places
            CSTRNG = CHEADER(I)
            WRITE(*,5000) " #",CSTRNG(1:LENOCC(CSTRNG))
         ENDDO
      ENDIF

c - Read the fastNLO coefficient table
      IF (FILENAME.NE.OLDFILENAME) THEN
         CALL FX9999RD(FILENAME)
         OLDFILENAME = FILENAME
         
c - Check consistency of array dimensions / commonblock parameters
         IF (NXSUM.GT.NXMAX) THEN
            WRITE(*,*)"fastNLO: ERROR! Inconsistent x binning,"//
     >           " NXSUM =",NXSUM," > NXMAX =",NXMAX
            STOP
         ENDIF
         IF (NSCALEVAR.GT.NSCALEMAX) THEN
            WRITE(*,*)"fastNLO: ERROR! Inconsistent scale variations,"
     >           //", NSCALEVAR =",NSCALEVAR," > NSCALEMAX =",NSCALEMAX
            STOP
         ENDIF
         IF (NRAPIDITY.GT.NRAPMAX) THEN
            WRITE(*,*)"fastNLO: ERROR! Inconsistent rap. binning,"//
     >           " NRAPIDITY =",NRAPIDITY," > NRAPMAX =",NRAPMAX
            STOP
         ENDIF
         DO I=1,NRAPIDITY
            IF (NPT(I).GT.NPTMAX) THEN 
               WRITE(*,*)"fastNLO: ERROR! Inconsistent pt binning,"//
     >              " NPT(",I,") =",NPT(I)," > NPTMAX =",NPTMAX
               STOP
            ENDIF
         ENDDO
         K = 0
         DO I=1,NRAPIDITY
            DO J=1,NPT(I)
               K=K+1
            ENDDO
         ENDDO
         IF (K.GT.NBINTOTMAX) THEN
            WRITE(*,*)"fastNLO: ERROR! Too many bins =",K,
     >           " > NBINTOTMAX=",NBINTOTMAX
            STOP
         ENDIF
         IF (NBINTOT.GT.NBINTOTMAX) THEN
            WRITE(*,*)"fastNLO: ERROR! NBINTOT =",NBINTOT,
     >           " > NBINTOTMAX =",NBINTOTMAX
            STOP
         ENDIF

c - Print further info
         WRITE(*,*)"#      "
         CSTRNG = NAMELABEL(1)
         WRITE(*,5000)" #      this table contains: ",
     >        CSTRNG(1:LENOCC(CSTRNG))
         CSTRNG = NAMELABEL(2)
         WRITE(*,5000)" #      as published in:     ",
     >        CSTRNG(1:LENOCC(CSTRNG))
         CSTRNG = NAMELABEL(3)
         WRITE(*,5000)" #      by:                  ",
     >        CSTRNG(1:LENOCC(CSTRNG))
         WRITE(*,*)"#      "
         CSTRNG = CIREACTION(IREACTION)
         WRITE(*,*)"#      reaction: ",CSTRNG(1:LENOCC(CSTRNG))
         CSTRNG = CIPROC(IPROC)
         WRITE(*,*)"#      process: ",CSTRNG(1:LENOCC(CSTRNG))
         WRITE(*,*)"#      total No. of observable bins:",Nbintot
         CSTRNG = CIALGO(IALGO)
         WRITE(*,*)"#      jet algo: ",CSTRNG(1:LENOCC(CSTRNG))
         CSTRNG = CJETRES1(IALGO)
ckr 30.01.2008: Use explicit format for comparisons
         WRITE(*,FMT='(A,A,A,F6.4)')
     >        " #         parameter 1: ",CSTRNG(1:LENOCC(CSTRNG)),
     >        " = ",JETRES1
         CSTRNG = CJETRES2(IALGO)
         WRITE(*,FMT='(A,A,A,F6.4)')
     >        " #         parameter 2: ",CSTRNG(1:LENOCC(CSTRNG)),
     >        " = ",JETRES2
         WRITE(*,*)"#"
         WRITE(*,*)"#"
         WRITE(*,*)"#      the single contributions have been computed"
         WRITE(*,*)"#      using the following codes:"
         DO I=1,NORD
            CSTRNG = POWLABEL(I)
            WRITE(*,5000)" #      ",CSTRNG(1:LENOCC(CSTRNG))
            CSTRNG = CODELABEL(I)
            WRITE(*,5000)" #      by: ",CSTRNG(1:LENOCC(CSTRNG))
         ENDDO
         WRITE(*,*)"#"
         IF (NORD.GT.0) THEN
            DO I=1,4
               CSTRNG = CNLOJET(I)
               WRITE(*,5000)" # ",CSTRNG(1:LENOCC(CSTRNG))
            ENDDO
         ENDIF 
         IF (NORD.EQ.3) THEN
            DO I=1,5
               CSTRNG = CTHRCOR(I)
               WRITE(*,5000)" # ",CSTRNG(1:LENOCC(CSTRNG))
            ENDDO
         ENDIF 

c - Print scale-variations available in the table
         WRITE(*,*)"#"
         WRITE(*,*)"#   --- the renormalization and factorization "//
     >        "scales mur, muf"
         WRITE (*,*)"#       are proportional to"
         CSTRNG = SCALELABEL
         WRITE (*,5000)" #           mu0 = ",CSTRNG(1:LENOCC(CSTRNG))
         WRITE (*,*)"#"
         WRITE (*,*)"#   --- available No. of scale",
     >        " variations:",NSCALEVAR
         WRITE (*,*)"#         available factorization scale settings:"
         DO I=1,NSCALEVAR
            WRITE (*,FMT='(A,I1,A,F6.4)')
     >           " #           ",I,"  (muf/mu0) = ",MUFSCALE(I)
         ENDDO
         WRITE (*,*)"#"
         WRITE (*,*)"#         available renormalization scale "//
     >        "settings:"
         DO I=1,NSCALEVAR
            WRITE (*,FMT='(A,I1,A,F6.4)')
     >           " #           ",I,"  (mur/mu0) = ",MURSCALE(I)
         ENDDO
         WRITE (*,*)"#   (In LO and NLO, the renormalization scale"
         WRITE (*,*)"#    can be varied arbitrarily afterwards."
         WRITE (*,*)"#    This is, however, not possible for the"
         WRITE (*,*)"#    2-loop threshold corrections.)"
         WRITE (*,*)"# "
      ENDIF

c - Identify the scales chosen in this call
      IXMUR = 0
      IXMUF = 0
      DO I=1,NSCALEVAR
         IF (ABS(XMUR/MURSCALE(I)-1D0).LT.0.000001) IXMUR=I
         IF (ABS(XMUF/MUFSCALE(I)-1D0).LT.0.000001) IXMUF=I
      ENDDO

      IF (IFIRST.EQ.0) THEN
         WRITE (*,*)"#    --> in the first call the scales are "//
     >        "chosen to be:"
         WRITE (*,5002)" #     (mur/mu0) =",XMUR,
     >        "(muf/mu0) =",MUFSCALE(IXMUF)
      ENDIF
      IF (IXMUF.EQ.0) THEN
         WRITE(*,*)"# factorization scale ",XMUF,
     >        " not available in table"
         GOTO 998
      ENDIF
      IF (NORD.EQ.3 .AND. IXMUR.NE.IXMUF) THEN
         WRITE(*,*)"# renormalization scale  mur<>muf  not available ",
     >        "for threshold corrections"
         WRITE(*,*)"#                (only mur=muf)"
      ENDIF

      IF (IFIRST.EQ.0) THEN
         IFIRST = 1
         WRITE(*,*)"#"
         WRITE(*,*)"#################################################",
     >        "################"
         WRITE(*,*)"#"
      ENDIF


c - Reset result arrays
      DO NBIN=1,NBINTOT         ! Cont. numbering of the final array
         DO L=1,NORD            ! All orders LO, NLO, NNLO, ...
            DO M=1,(NSUBPROC+1) ! No. of Subproc + 1 for total sum
               RESULT(NBIN,M,L) = 0D0 ! Reset 
            ENDDO
            XSECT(NBIN,L) = 0D0 ! Reset 
         ENDDO
      ENDDO

c - Now get the PDFs ...
      CALL FX9999GP(MUFSCALE(IXMUF))
c - ... and multiply with perturbative coefficients 
      CALL FX9999MT(XMUR,IXMUF)

c ----------------- final touches ------------------------------------
c - Sum subprocesses / fill result array / fill 'XSECT' array
      DO NBIN=1,NBINTOT         ! Cont. numbering of the final array
         DO IORD=1,NORD         ! Loop over all orders
            RESULT(NBIN,(NSUBPROC+1),IORD) = 0D0
            XSECT(NBIN,IORD) = 0D0
            DO M=1,NSUBPROC     ! No. of Subprocesses
               RESULT(NBIN,(NSUBPROC+1),IORD) =
     >              RESULT(NBIN,(NSUBPROC+1),IORD)+RESULT(NBIN,M,IORD)
            ENDDO
            XSECT(NBIN,IORD) = RESULT(NBIN,(NSUBPROC+1),IORD)
         ENDDO
      ENDDO

c -----------------------------------------------------------------
c --- Special feature:    compute   1/sigma * dsigma/d[whatever]
c --------- normalization - if required for observable ------------
c
      INORMFLAG = 0             ! switch normnalization Off(=0)/On(=1)
c
      IF (INORMFLAG.EQ.1) THEN
c --- Assume: observable has been divided by binwidth in 2nd dimension
c     here: normalize LO & NLO results (the NLO numbers are normalized
c                    such that the sum (LO+NLO) is properly normalized)
c   
         DO ISUB=1,(NSUBPROC+1) ! subprocess: Nsubproc + 1 tot
            NBIN=0
            DO I=1,NRAPIDITY                   
               DO IORD=1,NORD
                  SUM(IORD) = 0D0
               ENDDO
               DO J=1,NPT(I)    ! *** step 1: compute norm. factor
                  NBIN = NBIN + 1
                  DO IORD=1,NORD ! Order: tot, LO, NLO-corr, 3 NNLOcorr
                     DO K=1,IORD ! Sum over all lower orders
                        SUM(IORD) = SUM(IORD) +  
     >                       RESULT(NBIN,ISUB,K)*(PTBIN(I,J+1) -
     >                       PTBIN(I,J))
                     ENDDO
                  ENDDO
               ENDDO
               NBIN = NBIN - NPT(I)
               DO J=1,NPT(I)    ! *** step 2: Apply norm. factor
                  NBIN = NBIN + 1
                  IF (NORD.GE.1) RESULT(NBIN,ISUB,1) =
     >                 RESULT(NBIN,ISUB,1) / SUM(1)
                  IF (NORD.GE.2) RESULT(NBIN,ISUB,2) =
     >                 (RESULT(NBIN,ISUB,1)*(SUM(1)-SUM(2))
     >                 + RESULT(NBIN,ISUB,2))  / SUM(2)
               ENDDO 
            ENDDO               ! Rapidity-Loop (or: 1st Dimension)
         ENDDO                  ! Subprocess-Loop
c - Fill "xsect" array
         DO NBIN=1,NBINTOT      ! Cont. numbering of the final array
            DO IORD=1,NORD      ! Loop over all orders
               XSECT(NBIN,IORD) = RESULT(NBIN,(NSUBPROC+1),IORD)
            ENDDO
         ENDDO
      ENDIF                     ! End: norm. 1/sigma dsigma/d[s.th.]

c - Print results - if requested
      IF (IPRINTFLAG.EQ.1) CALL FX9999PR(XMUR,IXMUF)

      RETURN

ckr 30.01.2008: Use simple A format, string length via LENOCC
 5000 FORMAT (A,A)
 5001 FORMAT (A,A,A)
 5002 FORMAT (A,F8.4,4X,A,F8.4)
 998  CONTINUE
      END



*******************************************************************
      SUBROUTINE FX9999MT(XMUR,IXMUF)
*-----------------------------------------------------------------
* MW 08/26/2005
*
* multiply the PDFs and the coefficients for a single scale setting
* 
* input:   XMUR  prefactor for nominal renormalization scale
*          IXMUF No.of factorization scale setting (as stored in table)
*
*-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INTEGER IXMUF,I,J,K,L,M,IORD,JORD,NBIN
      INTEGER IPOSITION(5)
      DOUBLE PRECISION  MUR,AS(NSCALEBINMAX),ASPOW(5),FNALPHAS,COEFF 
      DOUBLE PRECISION BETA0,BETA1,PI,XMUR,LOGMU,SCFAC
ckr      DOUBLE PRECISION SCFAC2A,SCFAC2B ! For ren-scale variation
      DOUBLE PRECISION MU0SCALE,LLPTLO,LLPTHI,T1,T2,BWEIGHT
      INTEGER NF,CA
      DOUBLE PRECISION CF
      PARAMETER (PI=3.14159265358979323846, NF=5, CA=3, CF=4D0/3D0)
      PARAMETER (BETA0=(11D0*CA-2D0*NF)/3D0) 
      PARAMETER (BETA1=34*CA*CA/3D0-2D0*NF*(CF+5D0*CA/3D0))
      PARAMETER (MU0SCALE=0.25)

c - Get the absolute order in alpha_s of the LO contribution
      JORD = NPOW(1)

c - Vary renormalization scale around the value used in orig. calc.
      LOGMU = LOG(XMUR/MURSCALE(IXMUF)) ! change w.r.t. orig. calc.
      SCFAC = DBLE(JORD)*BETA0*LOGMU ! NLO contrib.
ckr 30.01.2008: Comment unused defs
c      scfac2a= dble(jord+1)*beta0 *logmu          ! NNLO contrib.
c      scfac2b= dble(jord*(jord+1))/2d0*beta0*beta0*logmu*logmu  
c     +     + dble(jord)*beta1/2d0*logmu           ! NNLO contrib.cont.

c - MW:  we may save time if we make the mur-variation later
c        for the whole contribution - instead of doing it for
c        each array element.

c - Position of scale/order in array
      IPOSITION(1) = 1
      IPOSITION(2) = 1+IXMUF+(2-2)*NSCALEVAR
      IPOSITION(3) = 1+IXMUF+(3-2)*NSCALEVAR

c - Loop over coefficient array
      NBIN = 0                  ! Cont. numbering for the final array
      DO I=1,NRAPIDITY          ! (Pseudo-)Rapidity Bins
         DO J=1,NPT(I)          ! ET/pT Bins
            NBIN=NBIN+1         ! Continuous bin No.
            LLPTLO = LOG(LOG(XMUR * MURVAL(I,J,1)/MU0SCALE))
            LLPTHI = LOG(LOG(XMUR * MURVAL(I,J,NSCALEBIN)/MU0SCALE))
            IF (NSCALEBIN.EQ.1) THEN
               AS(1) = FNALPHAS(XMUR * MURVAL(I,J,1))
            ELSE
               DO L=1,NSCALEBIN ! Loop over all scale bins, calc alphas
                  MUR = MU0SCALE * EXP(EXP((LLPTLO + (DBLE(L)-1.)
     >                 /(NSCALEBIN-1)*(LLPTHI-LLPTLO) )))
                  AS(L) = FNALPHAS(MUR) !    ... and get alpha_s
               ENDDO
            ENDIF
            DO L=1,NSCALEBIN
               IF (NSCALEBIN.EQ.1) THEN
                  BWEIGHT = AS(1)
               ELSEIF (NSCALEBIN.EQ.2) THEN
                  IF (L.EQ.1) BWEIGHT = AS(1)
                  IF (L.EQ.2) BWEIGHT = AS(2)
               ELSEIF (NSCALEBIN.EQ.3) THEN 
                  T1 = 1./2.
                  IF (L.EQ.1) BWEIGHT = AS(1)
                  IF (L.EQ.2) BWEIGHT = 1./(2.*(1.-T1)*T1)*
     >                 (AS(2)-AS(1)*(1.-T1)**2-AS(3)*(T1)**2)
                  IF (L.EQ.3) BWEIGHT = AS(3)
               ELSEIF (NSCALEBIN.EQ.4) THEN 
                  T1 = 1./3.
                  T2 = 2./3.
                  IF (L.EQ.1) BWEIGHT = AS(1)
                  IF (L.EQ.2) BWEIGHT = 1./(3.*(1.-T1)**2*T1 * 
     >                 3.*(1.-T2)*(T2)**2 - 3.*(1.-T2)**2*T2 * 
     >                 3.*(1.-T1)*(T1)**2)
     >                 *(3.*(1.-T2)*(T2)**2*(AS(2)-AS(1)*(1.-T1)**3
     >                 -AS(4)*(T1)**3) 
     >                 - 3.*(1.-T1)*(T1)**2*(AS(3)-AS(1)*(1.-T2)**3
     >                 -AS(4)*(T2)**3))
                  IF (L.EQ.3) BWEIGHT = 1./(3.*(1.-T1)**2*T1 *
     >                 3.*(1.-T2)*(T2)**2
     >                 - 3.*(1.-T2)**2*T2 * 3.*(1.-T1)*(T1)**2)
     >                 *(3.*(1.-T1)**2*(T1)*(AS(3)-AS(1)*(1.-T2)**3
     >                 -AS(4)*(T2)**3) 
     >                 - 3.*(1.-T2)**2*(T2)*(AS(2)-AS(1)*(1.-T1)**3
     >                 -AS(4)*(T1)**3))
                  IF (L.EQ.4) BWEIGHT = AS(4)
               ELSE
                  WRITE(*,*)"fastNLO: ERROR! NSCALEBIN = ",NSCALEBIN,
     >                 " not yet supported."
                  STOP
               ENDIF
               DO IORD=1,NORD
                  ASPOW(IORD) = BWEIGHT**NPOW(IORD)
               ENDDO
               DO K=1,NXSUM     ! Loop over all x bins
                  DO M=1,NSUBPROC ! Loop over subprocesses
                     DO IORD=1,NORD ! Rel. order: 1 LO 2 NLO 3 NNLO ...
                        IF (IORD.EQ.1) THEN ! LO contribution
                           COEFF = ARRAY(NBIN,K,M,IPOSITION(1),L)
                        ELSEIF (IORD.EQ.2) THEN ! NLO CONTRIBUTIONS
                           COEFF = 
     >                          ARRAY(NBIN,K,M,IPOSITION(2),L)
     >                          + SCFAC*ARRAY(NBIN,K,M,1,L) 
                        ELSEIF (IORD.EQ.3) THEN !2-LOOP THRESHOLD CORR.
c
c     - the following works only for "true" higher orders (NNLO)
c     -> not for 2-loop threshold corrections (N. Kidonakis,10.01.2006)
c     coeff = 
c     +                       array(nbin,k,m,(1+ixmuf+(iord-2)*nscalevar),l)
c     +                       array(nbin,k,m,iposition(3),l)
c     +                       + scfac2a*array(nbin,k,m,(1+ixmuf))
c     +                       + scfac2b*array(nbin,k,m,1) 
c     
c     ... therefore the NLLO-NLL contributions are only available for mu_r=mu_f
c     -           in other words: for  log(mur/muf)=0
                           IF (LOGMU.EQ.0D0) THEN
                              COEFF = ARRAY(NBIN,K,M,IPOSITION(3),L)
                           ELSE
                              COEFF = 0D0
                           ENDIF
                        ENDIF

c - For 'standard' fastNLO tables 
                        IF (IREF.EQ.0 .OR. I.LE.(NRAPIDITY/2)) THEN 
                           RESULT(NBIN,M,IORD) = RESULT(NBIN,M,IORD)
     >                          + COEFF 
     >                          * ASPOW(IORD) ! Mult.w.(alpha_s/2Pi)**N
     >                          * PDF(NBIN,K,M,L) ! Mult.with PDFs
c - For 'reference' fastNLO tables including PDF/alphas
c - Only relevant for fastNLO authors -> for precision studies
                        ELSE
                           IF (L.EQ.1) THEN !Ref.is stored in sc. bin 1
                              RESULT(NBIN,M,IORD) =
     >                             RESULT(NBIN,M,IORD) + COEFF
                           ENDIF
                        ENDIF
                     ENDDO      ! iord perturbative order
                  ENDDO         ! l scale-bins
               ENDDO            ! m subprocess
            ENDDO               ! k x-bin
         ENDDO                  ! j pt
      ENDDO                     ! i rapidity
      
      RETURN
      END



*******************************************************************
      SUBROUTINE FX9999GP(muffactor)
*-----------------------------------------------------------------
* MW 08/26/2005
*
* read the PDFs in all (rap,pT) bins - at all xmax, xmin bins
* the default factorization scale (in GeV) is multiplied by muffactor
*
* input: muffactor - prefactor for nominal factorization scale
*
*-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      DOUBLE PRECISION MUFFACTOR
      INTEGER NBIN,NX, i,j,k,l,m,p,nx2limit
      DOUBLE PRECISION x1, xlim, hx, hxlim, muf, H(7),
     +   reweight, NEWPDF(-6:6), XPDF(nxmax,-6:6)
c      DOUBLE PRECISION hxinv3   ! invert h(x) for ixscheme=3: log(1/x)+x-1
      DOUBLE PRECISION mu0scale,llptlo,llpthi
      PARAMETER (mu0scale=0.25)

      nbin=0
      do i=1,nrapidity          ! Rapidity Bins
         do j=1,NPT(i)          ! pT Bins
            nbin=nbin+1
            xlim = xlimit(i,j)
            if (ixscheme.eq.2) then
               hxlim = - sqrt(-log10(xlim))
            elseif (ixscheme.eq.1) then
               hxlim = log10(xlim)
            elseif (ixscheme.eq.3) then
               hxlim = log10(xlim)+xlim-1d0
            else
               WRITE(*,*)" fastNLO - IXSCHEME ",ixscheme
     >              ," not available"
               Stop
            endif

            llptlo = log(log( muffactor * mufval(i,j,1)/mu0scale))
            llpthi = log(log( muffactor * mufval(i,j,NSCALEBIN)/
     >           mu0scale))
            do p=1,NSCALEBIN
               if(NSCALEBIN.eq.1) then
                  muf =  muffactor * mufval(i,j,1)
               else
                  muf = mu0scale * exp(exp((llptlo + (dble(p)-1.)
     >                 /(NSCALEBIN-1)*(llpthi-llptlo) )))
               endif
c - get PDFs(-6:6) directly from interface, reweight, copy into linear array
               do k=1,NXTOT     ! loop over all x-values
                  hx = hxlim *(1d0 - dble(k-1)/dble(nxtot)) ! compute x1-value
                  if (ixscheme.eq.2) then
                     x1 = 10**(-(hx*hx))  !best scheme: sqrt(log10(1/x)
                  elseif (ixscheme.eq.1) then
                     x1 = 10**(hx)      ! simple log10(1/x)
c                  elseif (ixscheme.eq.3) then
c                     x1 = hxinv3(hx)    ! inefficient: log10(1/x)+x-1
                  else
                     WRITE(*,*)" fastNLO - IXSCHEME ",ixscheme,
     +                    " not available"
                     Stop
                  endif

                  call FNPDF(x1,muf,newpdf)
                  reweight = 1d0
                  IF (IPDFWGT.eq.1) then ! standard fastNLO reweighting
                     reweight = sqrt(x1)/(1d0-0.99d0*x1)**3
                  elseif (IPDFWGT.eq.0) then ! no reweighting
                     reweight = 1d0
                  else
                     WRITE(*,*)" fastNLO - reweighting scheme not "//
     >                    "available: ",IPDFWGT
                  ENDIF
                  do l=-6,6
                     XPDF(k,l) = newpdf(l) * reweight
                  enddo
               enddo  

c - now fill main PDF array - compute different lin. comb for diff sub-proc
               nx = 0
               do k=1,NXTOT
                  nx2limit = k
                  if (ireaction.eq.1) nx2limit=1
                  do l=1,nx2limit
                     nx=nx+1
                     call fx9999pl(ireaction,k,l,XPDF,H)
                     do m=1,Nsubproc
                        pdf(nbin,nx,m,p) = H(m)
                     enddo
                  enddo
               enddo 

            enddo
         enddo
      enddo

      RETURN 
      END

*******************************************************************
      SUBROUTINE FX9999PL(ireact,i,j,XPDF,H)
* ---------------------------------------------------------------
* MW 03/27/06
* compute PDF linear combinations - for different reactions
*
* depending on ireaction, the product of the i-th and j-th entries
* of the PDF array XPDF are multiplied into the relevant linear 
* combinations in the array H
*
* input:
*    ireact             flag for reaction (1:DIS, 2:pp-jets, 3:ppbar-jets)
*    i                  x-index of first hadron        
*    j                  x-index for second hadron (if two-hadron process)
*    XPDF(nxmax,-6:6)   PDF array for all x-bins
*
* output:
*    H(10)              PDF linear combinations
* ---------------------------------------------------------------
      Implicit None
      INCLUDE 'fnx9999.inc'
      Integer ireact, i,j,k
      Double Precision XPDF(nxmax,-6:6), H(10),
     +     G1, G2,              ! gluon densities from both hadrons
     +     SumQ1, SumQ2,        ! sum of quark densities
     +     SumQB1, SumQB2,      ! sum of anti-quark densities
     +     Q1(6),Q2(6), QB1(6),QB2(6), ! arrays of 6 (anti-)quark densities
     +     S,A                  ! products S,A

c --- for DIS ---
      if (Ireaction .eq. 1) then 
         H(1) = XPDF(i,0)       ! Gluon
         H(2) = 0d0             ! Sigma
         H(3) = 0d0             ! Delta
         do k=1,6
            H(2) = H(2)+XPDF(i,k)+XPDF(i,-k)
         enddo
         do k=1,5,2
            H(3) = H(3) + (XPDF(i,k)+XPDF(i,-k)+
     +           4d0*(XPDF(i,k+1)+XPDF(i,-k-1)))/9d0
         enddo

c --- for hadron-hadron ---
      elseif (Ireact .eq. 2 .or. Ireact .eq. 3) then ! pp/ppbar
         SumQ1  = 0d0
         SumQB1 = 0d0
         SumQ2  = 0d0
         SumQB2 = 0d0
         do k=1,6
            Q1(k)  = XPDF(i,k)  ! read x1
            QB1(k) = XPDF(i,-k)
            SumQ1  = SumQ1  + Q1(k)
            SumQB1 = SumQB1 + QB1(k)
            Q2(k)  = XPDF(j,k)  ! read x2
            QB2(k) = XPDF(j,-k)
            SumQ2  = SumQ2  + Q2(k)
            SumQB2 = SumQB2 + QB2(k)
         enddo
         G1     = XPDF(i,0)
         G2     = XPDF(j,0)
c  -  compute S,A
         S = 0d0
         A = 0d0
         do k=1,6
            S = S + (Q1(k)*Q2(k)) + (QB1(k)*QB2(k)) 
            A = A + (Q1(k)*QB2(k)) + (QB1(k)*Q2(k)) 
         enddo
c  - compute seven combinations
         H(1) = G1*G2
         H(2) = (SumQ1+SumQB1)*G2
         H(3) = G1*(SumQ2+SumQB2)
c  - for pp
         if (Ireact .eq. 2) then
            H(4) = SumQ1*SumQ2 + SumQB1*SumQB2 - S
            H(5) = S
            H(6) = A
            H(7) = SumQ1*SumQB2 + SumQB1*SumQ2 - A
c  - for p-pbar: swap combinations 4<->7 and 5<->6
         elseif (Ireact .eq. 3) then
            H(7) = SumQ1*SumQ2 + SumQB1*SumQB2 - S
            H(6) = S
            H(5) = A
            H(4) = SumQ1*SumQB2 + SumQB1*SumQ2 - A
         endif
      else
         WRITE(*,*)"    ireaction =",ireact
         WRITE(*,*)" this reaction is not yet defined"
         stop
      endif

      RETURN 
      END
*******************************************************************

      SUBROUTINE FX9999PR(XMUR,IMUFFLAG)
*-----------------------------------------------------------------
* m. Wobisch  - print fastNLO cross section results
*
* input:
*   XMUR        prefactor of nominal renormalization scale
*   IMUFFLAG    print results for factorization scale No. IMUFFLAG
*
* MW 04/18/2005
*-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IMUFFLAG, nbin, i,j, nproc
      DOUBLE PRECISION XMUR
      INCLUDE 'fnx9999.inc'

      WRITE(*,5000)" --  fastNLO - results for  ",NAMELABEL(1)

      nproc = Nsubproc+1        ! print sum of all subprocesses
c      nproc = 1                 ! print only gg->jets subprocess

      if (nord.le.2) then ! - LO and NLO output only
         WRITE(*,*)"--    cross sections:              ",
     +        "LO          NLOcorr      total"
      elseif (nord.eq.3) then   ! incl 2-loop threshold correction term
         WRITE(*,*)"--    cross sections:              ",
     +        "LO         NLOcorr      2-loop       total"
      else
         WRITE(*,*)" fastNLO error: Nord=",nord," is not possible"
         stop
      endif

      WRITE(*,*)"--------  muf/mu0=",mufscale(imufflag),
     +     "     mur/mu0=",xmur
   
      nbin = 0                  ! linear bin for the final observable
      do i=1,nrapidity          ! Rapidity Bins
         if (RAPBIN(i).lt.RAPBIN(i+1)) then
            WRITE(*,*)"       from ",RAPBIN(i)," - ",RAPBIN(i+1),
     +           "  in:  ",dimlabel(1)
         else                   !  happens only for reference tables
            WRITE(*,*)"       from ",RAPBIN(1)," - ",RAPBIN(i+1),
     +           "  in:  ",dimlabel(1)
         endif
         do j=1,NPT(i)          ! pT Bins
            nbin=nbin+1
            if (nord.le.2) then ! - LO and NLO output only
               WRITE(*,900) dimlabel(2),PTBIN(i,j),PTBIN(i,j+1),
     +              real(result(nbin,(nproc),1)),
     +              real(result(nbin,(nproc),2)),
     +              real(result(nbin,(nproc),1)
     +              + result(nbin,(nproc),2))
            elseif (nord.eq.3) then  ! incl NNLO or 2-loop term
               WRITE(*,901) dimlabel(2),PTBIN(i,j),PTBIN(i,j+1),
     +              real(result(nbin,(nproc),1)),
     +              real(result(nbin,(nproc),2)),
     +              real(result(nbin,(nproc),3)),
     +              real(result(nbin,(nproc),1)
     +              + result(nbin,(nproc),2)
     +              + result(nbin,(nproc),3))
            endif
         enddo
      enddo
      WRITE(*,*)" "

      RETURN
ckr 30.01.2008: Change format for better comparison
c 900  Format (A12,F8.2,"-",F8.2,":",3E13.4)
 900  Format (A12,F8.2,"-",F8.2,":",3E17.8)
 901  Format (A12,F8.2,"-",F8.2,":",4E13.4)
 5000 FORMAT (A,A64)
      END

*******************************************************************
      SUBROUTINE FX9999RD(FILENAME)
*-----------------------------------------------------------------
* M. Wobisch   read ASCII table of perturbative coefficients
*
* input: FILENAME  name of table
*
* MW 04/15/2005
* MW 01/30/2006 -> now reading tableformat v1.4
* MW 06/02/2006 -> check table format
*-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) FILENAME
      CHARACTER*255 CSTRNG
      INTEGER IFILE, LENOCC, I,J,K,L,M,N,   NBIN
      INCLUDE 'fnx9999.inc'

      OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=IFILE)
      IF (IFILE .ne. 0) THEN
         WRITE(*,*)"          fastNLO:  table file not found ",
     +        "  -  IOSTAT = ",IFILE
         STOP
      ENDIF

      READ(2,*) i
      if (i.ne.iseparator) goto 999
      READ(2,*) ITABVERSION
      WRITE(*,*)"#       tableformat is version",
     >     real(itabversion)/10000d0
      if (ITABVERSION.ne.14000) then
         WRITE(*,*)"#     ==> this usercode works only for version 1.4"
         WRITE(*,*)"#     ==> please get updated usercode from"
         WRITE(*,*)"#         http://hepforge.cedar.ac.uk/fastnlo "
         stop
      endif
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c  -----------------------------------
      READ(2,*) IREACTION
      READ(2,*) ECMS
      READ(2,*) IXSECTUNITS
      do i=1,5
         READ(2,FMT='(A)') NAMELABEL(i)
      enddo
      READ(2,*) IPROC
      READ(2,*) IALGO
      READ(2,*) JETRES1
      READ(2,*) JETRES2
      READ(2,*) NORD
      do i=1,nord
         READ(2,*) NPOW(i)
      enddo
      do i=1,nord
         READ(2,*) POWLABEL(i)
         READ(2,*) CODELABEL(i)
      enddo
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      do i=1,nord
         READ(2,*) NEVT(i)
      enddo
      do i=1,nord
         CSTRNG = POWLABEL(I)
         WRITE(*,5000)
     &        " #      ",NEVT(i)," events in ",CSTRNG(1:LENOCC(CSTRNG))
      enddo
      READ(2,*) NXTOT
      WRITE(*,*)"#           No. of x bins:",NXTOT
      READ(2,*) IXSCHEME
      READ(2,*) IPDFWGT
      READ(2,*) IREF
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      READ(2,*) NBINTOT
      READ(2,*) NDIMENSION
      do i=1,ndimension
         READ(2,*) DIMLABEL(i)
      enddo
      
      READ(2,*) NRAPIDITY
      do i=1,nrapidity+1
         READ(2,*) RAPBIN(i)
      enddo
      do i=1,nrapidity
         READ(2,*) NPT(i)
      enddo
      do i=1,nrapidity
         do j=1,NPT(i)+1
            READ(2,*) PTBIN(i,j)
         enddo
      enddo
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      do i=1,nrapidity
         do j=1,NPT(i)
            READ(2,*) XLIMIT(i,j)
         enddo
      enddo
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      READ(2,*) SCALELABEL
      READ(2,*) NSCALEBIN
      do i=1,nrapidity
         do j=1,NPT(i)
            do k=1,NSCALEBIN
               READ(2,*) MURVAL(i,j,k)
            enddo
         enddo
      enddo
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      do i=1,nrapidity
         do j=1,NPT(i)
            do k=1,NSCALEBIN
               READ(2,*) MUFVAL(i,j,k)
            enddo
         enddo
      enddo
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      READ(2,*) NSCALEVAR
      do i=1,NSCALEVAR
         READ(2,*) MURSCALE(i)
      enddo
      do i=1,NSCALEVAR
         READ(2,*) MUFSCALE(i)
      enddo
      
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      if((IREACTION.eq.2).OR.(IREACTION.eq.3)) then ! pp or ppbar
         NXSUM = (NXTOT*NXTOT+NXTOT)/2
         Nsubproc = 7
      elseif(IREACTION.eq.1) then ! DIS
         NXSUM = NXTOT
         Nsubproc = 3
      endif

      nbin=0
      do i=1,nrapidity
         do j=1,NPT(i)
            nbin=nbin+1         ! linear numbering of (rap,pT) bins
            do k=1,NXSUM        ! tot. No. x-Bins
               do m=1,Nsubproc  ! No. of Subproc
                  do n=1,1+NSCALEVAR*(NORD-1) ! LO & NLO & w/ scale var
                     do l=1,NSCALEBIN ! No. of Bins in Scale
                        READ(2,*) array(nbin,k,m,n,l)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      READ(2,*) i
      if (i.ne.iseparator) goto 999
      CLOSE(2)

      RETURN
 999  continue

      close (2) 
      WRITE(*,*)" >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      WRITE(*,*)" >>>>>   fastNLO error in table format "
      WRITE(*,*)" >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      stop
      RETURN

ckr 30.01.2008: Use free format A, string length det. by LENOCC
 5000 FORMAT (A,I12,A,A) 
      END

*******************************************************************
c      Double Precision Function hxinv3(hx)
*----------------------------------------------------------------
*    compute inverse of: h(x)=log10(1/x)+x-1
*      -> return x-value for given h(x) 
*
* MW 04/26/2006  early poor version - extremely slow & silly
*                only for test purposes - not meant for real work
*                too slow - and maybe poor precision - 
*                and limited to x>10**-4
*----------------------------------------------------------------
c      implicit none
c      Double Precision hx,x, hxtest
c      Integer i,j,nmax
c
c      nmax=10000 ! defines the precision
c
c      do i=0,nmax
c         x = 10d0**(-4d0*dble(i)/dble(nmax))
c         hxtest = log10(x)+x-1d0
c         if (hxtest.lt.hx) goto 100
c      enddo
c
c 100  Continue
c      hxinv3=x
c
c      Return
c      End
