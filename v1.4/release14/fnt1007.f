*******************************************************************
*******************************************************************
* fastNLO user code                 T. Kluge, M. Wobisch 02/01/2006      
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
*******************************************************************
*******************************************************************



*******************************************************************
      SUBROUTINE FT1007CC(FILENAME, XMUR, XMUF, IPRINTFLAG, XSECT)
*-----------------------------------------------------------------
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
*                                3 NNLO / or other correction (if available)
*
* MW 04/15/2005 initial version
* MW 09/02/2005 implement flexible scale variations
* TK 12/07/2005 table format contains now LO + NLO with 5 scale variations
* MW 2006/01/17 implement tableformat version 1c
* MW 2006/02/01 implement tableformat version 1.4 
*-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fnt1007.inc'
      INTEGER IFIRST, IFILE, iord, I,J,K,L,M, IPRINTFLAG, 
     +     maxscale, nbin,nx, ixmur,ixmuf
      CHARACTER*(*) FILENAME
      CHARACTER*50 OLDFILENAME
      DOUBLE PRECISION xmur, xmuf
      DATA IFIRST/0/, OLDFILENAME/'xxxx'/
      SAVE IFIRST, OLDFILENAME

c - output in first call
      if (IFIRST.EQ.0) then
         Do i=1,12
            write(*,5000) " # ",CHEADER(i)
         enddo
      endif

c --- read the fastNLO coefficient table
      IF (FILENAME.ne.OLDFILENAME) THEN
         call FT1007RD(FILENAME)
         OLDFILENAME=FILENAME

c - check consistency of array dimensions / commonblock parameters
         if (NXSUM.gt.NXMAX) goto 999
         if (NSCALEVAR.gt.NSCALEMAX) goto 999
         if (NRAPIDITY.gt.NRAPMAX) goto 999
         do i=1,nrapidity
            if (npt(i).gt.NPTMAX) goto 999
         enddo
         k = 0
         do i=1,nrapidity
            do j=1,npt(i)
               k=k+1
            enddo
         enddo
         if (k.gt.NBINTOTMAX) goto 999
         if (NBINTOT.gt.NBINTOTMAX) goto 999

c - print further info
         write(*,*) "#      "
         write(*,5000) " #      this table contains: ",NAMELABEL(1) 
         write(*,5000) " #      as published in:     ",NAMELABEL(2)
         write(*,5000) " #      by:                  ",NAMELABEL(3)
         write(*,*) "#      "
         write(*,*) "#      reaction: ",cireaction(ireaction)
         write(*,*) "#      process:  ",ciproc(iproc)
         write(*,*) "#      total No. of observable bins: ",Nbintot        
         write(*,*) "#      jet algo: ",cialgo(ialgo)
         write(*,*) "#         parameter 1:  ",cjetres1(ialgo),'=',jetres1
         write(*,*) "#         parameter 2:  ",cjetres2(ialgo),'=',jetres2
         write(*,*) "#"
         write(*,*) "#"
         write(*,*) "#      the single contributions have been computed"
         write(*,*) "#      using the following codes:"
         do i=1,nord
            write(*,5000) " #      ",powlabel(i)
            write(*,5000) " #         by: ",codelabel(i)
         enddo
         write(*,*)"#"
         IF (nord.gt.0) then
            do i=1,4
               write(*,5000) " # ",CNLOJET(i)
            enddo
         endif 
         IF (nord.eq.3) then
            do i=1,5
               write(*,5000) " # ",CTHRCOR(i)
            enddo
         endif 

c - print scale-variations available in the table
         write(*,*) "#"
         write (*,*)
     +        "#   --- the renormalization and factorization scales mur, muf"
         write (*,*) "#       are proportional to"
         write (*,5000) " #           mu0 = ",scalelabel
         write (*,*) "#  "

         write (*,*) "#   --- available No. of scales",
     +        " variations: ",nscalevar
         write (*,*) "#         available factorization scale settings:"
         do i=1,nscalevar
            write (*,*) "#           ",i,"  (muf/mu0)",mufscale(i)
         enddo
         write (*,*) "#"
         write (*,*) "#         available renormalization scale settings:"
         do i=1,nscalevar
            write (*,*) "#           ",i,"  (mur/mu0)=",murscale(i)
         enddo
         write (*,*) "#         (in LO and NLO, the renormalization scale"
         write (*,*) "#          can be varied arbitrarily afterwards."
         write (*,*) "#          This is, however, not possible for the"
         write (*,*) "#          2-loop threshold corrections.)"
         write (*,*) "# "
      ENDIF


c - identify the scales chosen in this call
      ixmur = 0
      ixmuf = 0
      do i=1,nscalevar
         if ( abs(xmur/murscale(i)-1d0).lt.0.000001 ) ixmur=i
         if ( abs(xmuf/mufscale(i)-1d0).lt.0.000001 ) ixmuf=i
      enddo

      if (IFIRST.EQ.0) then
         write (*,*) "#    --> in the first call the scales are chosen to be:"
         write (*,5002) " #     (mur/mu0) =",xmur,"(muf/mu0) =",mufscale(ixmuf)
      endif
      if (ixmuf.eq.0) then
         write(*,*) "# factorization scale ",xmuf," not available in table"
         goto 998
      endif
      if (nord.eq.3 .and. ixmur.ne.ixmuf) then
         write(*,*) "# renormalization scale  mur<>muf  not available ",
     +        "for threshold corrections"
         write(*,*) "#                (only mur=muf)"
      endif

      if (IFIRST.EQ.0) then
         IFIRST = 1
         write(*,*) "# "
         write(*,*) "#################################################",
     +        "#############"
         write(*,*) " "
      endif


c - reset result arrays
      do nbin=1,nbintot         ! continuous numbering of the final array
         do l=1,Nord            ! all orders LO, NLO, NNLO, ...
            do m=1,(Nsubproc+1) ! No.of Subproc + 1 for total sum
               result(nbin,m,l) = 0d0 ! reset 
            enddo
            xsect(nbin,l) = 0d0   ! reset 
         enddo
      enddo


c - now get the PDFs ...
      call FT1007GP(mufscale(ixmuf))
c - ... and multiply with perturbative coefficients 
      call FT1007MT(xmur,ixmuf)


c ----------------- final touches ------------------------------------
c - sum subprocesses / fill result array / fill 'XSECT' array
      do nbin=1,nbintot         ! continuos numbering of the final array
         do iord=1,Nord         ! loop over all orders
            result(nbin,(Nsubproc+1),iord) = 0d0
            xsect(nbin,iord) = 0d0
            do m=1,Nsubproc     ! No. of Subprocesses
               result(nbin,m,iord) = result(nbin,m,iord)
               result(nbin,(Nsubproc+1),iord)=
     +              result(nbin,(Nsubproc+1),iord) + result(nbin,m,iord)
            enddo
            xsect(nbin,iord) = result(nbin,(Nsubproc+1),iord)
         enddo
      enddo

c - print results - if requested
      if (IPRINTFLAG.eq.1) call FT1007PR(xmur,IXMUF)

      RETURN

 999  continue
      WRITE(*,*) 'fastNLO: error occured - parameter outside range -',
     +     ' check commonblock'
      STOP
 5000 FORMAT (A,A64)
 5001 FORMAT (A,A,A)
 5002 FORMAT (A,F9.4,4X,A,F9.4)
 998  continue
      END

*******************************************************************
      SUBROUTINE FT1007MT(xmur,ixmuf)
*-----------------------------------------------------------------
* MW 08/26/2005
*
* multiply the PDFs and the coefficients for a single scale setting
* 
* input:    XMUR  prefactor for nominal renormalization scale
*           IXMUF No.of factorization scale setting (as stored in table) 
*
*-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fnt1007.inc'
      INTEGER IXMUF, i,j,k,l,m,iord, jord,    nbin,nx
      INTEGER iposition(5)
      DOUBLE PRECISION  mur, as, aspow(5), FNALPHAS, coeff 
      DOUBLE PRECISION beta0,beta1,PI,xmur, logmu, scfac,
     +     scfac2a,scfac2b      ! for ren-scale variation
      INTEGER NF, CA
      DOUBLE PRECISION CF
      PARAMETER (PI=3.14159265358979323846, NF=5, CA=3, CF=4d0/3d0)
      PARAMETER (beta0=(11d0*CA-2d0*NF)/3d0) 
      PARAMETER (beta1=34*CA*CA/3d0-2d0*NF*(CF+5d0*CA/3d0))

c - get the absolute order in alpha_s of the LO contribution
      jord = NPOW(1)

c - vary renormalization scale around the value used in orig. calculation
      logmu = log(xmur/murscale(ixmuf)) ! change w.r.t. orig. calculation
      scfac  = dble(jord)  *beta0 *logmu          ! NLO contrib.
      scfac2a= dble(jord+1)*beta0 *logmu          ! NNLO contrib.
      scfac2b= dble(jord*(jord+1))/2d0*beta0*beta0*logmu*logmu  
     +     + dble(jord)*beta1/2d0*logmu           ! NNLO contrib. continued


c - MW:  we may save time if we make the mur-variation later
c        for the whole contribution - instead of doing it for
c        each array element.

c - position of scale/order in array
      iposition(1) = 1
      iposition(2) = 1+ixmuf+(2-2)*nscalevar
      iposition(3) = 1+ixmuf+(3-2)*nscalevar

c - loop over coefficient array
      nbin = 0                  ! continuos numbering for the final array
      do i=1,nrapidity          ! (Pseudo-)Rapidity Bins
         do j=1,NPT(i)          ! ET/pT Bins
            nbin=nbin+1         ! continuous bin No.
            do l=1,NSCALEBIN    ! loop over all scale bins
               mur = xmur * murval(i,j,l) ! set the ren. scale ...
               as = FNALPHAS(mur) !    ... and get alpha_s
               do iord=1,Nord
                  aspow(iord) = as**npow(iord)
               enddo
               do k=1,NXSUM     ! loop over all x bins
                  do m=1,Nsubproc ! loop over subprocesses
                     do iord=1,Nord ! relative order: 1 LO  2 NLO  3 NNLO ...
                        if (iord.eq.1) then ! LO contribution
                           coeff = array(nbin,k,m,iposition(1),l)
                        elseif (iord.eq.2) then ! NLO contributions
                           coeff = 
     +                          array(nbin,k,m,iposition(2),l)
     +                          + scfac*array(nbin,k,m,1,l) 
                        elseif (iord.eq.3) then !2-loop threshold corr.
c
c     - the following works only for "true" higher orders (NNLO)
c     -> not for 2-loop threshold corrections (N. Kidonakis, Jan 10, 2006)
c     coeff = 
cc     +                       array(nbin,k,m,(1+ixmuf+(iord-2)*nscalevar),l)
c     +                       array(nbin,k,m,iposition(3),l)
c     +                       + scfac2a*array(nbin,k,m,(1+ixmuf))
c     +                       + scfac2b*array(nbin,k,m,1) 
c     
c     ... therefore the NLLO-NLL contributions are only available for mu_r=mu_f
c     -           in other words: for  log(mur/muf)=0
                           if (logmu.eq.0d0) then
                              coeff = array(nbin,k,m,iposition(3),l)
                           else
                              coeff = 0d0
                           endif
                        endif

c - for 'standard' fastNLO tables 
                        if (iref.eq.0 .or. i.le.(nrapidity/2)) then 
                           result(nbin,m,iord) = result(nbin,m,iord)
     +                          + coeff
     +                          * aspow(iord) ! multiply w/ (alpha-s/2pi)**n
     +                          * pdf(nbin,k,m,l) ! multiply with PDFs
c - for 'reference' fastNLO tables including PDF/alphas
c - only relevant for fastNLO authors -> for precision studies
                        else
                           if(l.eq.1) then ! reference is stored in scale bin#1
                              result(nbin,m,iord) = result(nbin,m,iord) + coeff
                           endif
                        endif
                     enddo      ! iord perturbative order
                  enddo         ! l scale-bins
               enddo            ! m subprocess
            enddo               ! k x-bin
         enddo                  ! j pt
      enddo                     ! i rapidity

      RETURN
      END

*******************************************************************
      SUBROUTINE FT1007GP(muffactor)
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
      INCLUDE 'fnt1007.inc'
      DOUBLE PRECISION MUFFACTOR
      INTEGER NBIN,NX, i,j,k,l,m,n,p,nx2limit
      DOUBLE PRECISION x1, x2, xlim, hx, hxlim, muf, H(7), D(3), sum,
     +   reweight, NEWPDF(-6:6), XPDF(nxmax,-6:6)
c      DOUBLE PRECISION hxinv3   ! invert h(x) for ixscheme=3: log(1/x)+x-1

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
               write(*,*) ' fastNLO - IXSCHEME ',ixscheme,' not available'
               Stop
            endif

            do p=1,NSCALEBIN
               muf = mufval(i,j,p) * muffactor

c - get PDFs(-6:6) directly from interface, reweight, copy into linear array
               do k=1,NXTOT     ! loop over all x-values
                  hx = hxlim *(1d0 - dble(k-1)/dble(nxtot)) ! compute x1-value
                  if (ixscheme.eq.2) then
                     x1 = 10**-(hx*hx)  ! best scheme: sqrt(log10(1/x)
                  elseif (ixscheme.eq.1) then
                     x1 = 10**(hx)      ! simple log10(1/x)
c                  elseif (ixscheme.eq.3) then
c                     x1 = hxinv3(hx)    ! inefficient: log10(1/x)+x-1
                  else
                     write(*,*) ' fastNLO - IXSCHEME ',ixscheme,
     +                    ' not available'
                     Stop
                  endif

                  call FNPDF(x1,muf,newpdf)
                  reweight = 1d0
                  IF (IPDFWGT.eq.1) then ! standard fastNLO reweighting
                     reweight = sqrt(x1)/(1d0-0.99d0*x1)**3
                  elseif (IPDFWGT.eq.0) then ! no reweighting
                     reweight = 1d0
                  else
                     write(*,*) ' fastNLO - reweighting scheme not available: '
     +                    ,IPDFWGT
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
                     call ft1007pl(ireaction,k,l,XPDF,H)
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
      SUBROUTINE FT1007PL(ireact,i,j,XPDF,H)
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
      INCLUDE 'fnt1007.inc'
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
         write(*,*) '    ireaction =',ireact
         write(*,*) ' this reaction is not yet defined'
         stop
      endif

      RETURN 
      END
*******************************************************************

      SUBROUTINE FT1007PR(XMUR,IMUFFLAG)
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
      INTEGER IMUFFLAG, nbin, i,j, maxscale
      DOUBLE PRECISION XMUR
      INCLUDE 'fnt1007.inc'

      write(*,5000) " --  fastNLO - results for  ",NAMELABEL(1)

      if (nord.le.2) then ! - LO and NLO output only
         write(*,*) "--    cross sections:              ",
     +        "LO          NLOcorr      total"
      elseif (nord.eq.3) then   ! incl 2-loop threshold correction term
         write(*,*) "--    cross sections:              ",
     +        "LO         NLOcorr      2-loop       total"
      else
         write(*,*) " fastNLO error: Nord=",nord," is not possible"
         stop
      endif

      write(*,*) '--------  muf/mu0=',mufscale(imufflag),
     +     '     mur/mu0=',xmur
   
      nbin = 0                  ! linear bin for the final observable
      do i=1,nrapidity          ! Rapidity Bins
         if (RAPBIN(i).lt.RAPBIN(i+1)) then
            WRITE(*,*) "       from ",RAPBIN(i)," - ",RAPBIN(i+1),
     +           "  in:  ",dimlabel(1)
         else                   !  happens only for reference tables
            WRITE(*,*) "       from ",RAPBIN(1)," - ",RAPBIN(i+1),
     +           "  in:  ",dimlabel(1)
         endif
         do j=1,NPT(i)          ! pT Bins
            nbin=nbin+1
            if (nord.le.2) then ! - LO and NLO output only
               WRITE(*,900) dimlabel(2),PTBIN(i,j),PTBIN(i,j+1),
     +              real(result(nbin,(Nsubproc+1),1)),
     +              real(result(nbin,(Nsubproc+1),2)),
     +              real(result(nbin,(Nsubproc+1),1)
     +              + result(nbin,(Nsubproc+1),2))
            elseif (nord.eq.3) then  ! incl NNLO or 2-loop term
               WRITE(*,901) dimlabel(2),PTBIN(i,j),PTBIN(i,j+1),
     +              real(result(nbin,(Nsubproc+1),1)),
     +              real(result(nbin,(Nsubproc+1),2)),
     +              real(result(nbin,(Nsubproc+1),3)),
     +              real(result(nbin,(Nsubproc+1),1)
     +              + result(nbin,(Nsubproc+1),2)
     +              + result(nbin,(Nsubproc+1),3))
            endif
         enddo
      enddo
      write(*,*) " "

      RETURN
 900  Format (A12,F8.2,"-",F8.2,":",3E13.4)
 901  Format (A12,F8.2,"-",F8.2,":",4E13.4)
 5000 FORMAT (A,A64)
      END

*******************************************************************
      SUBROUTINE FT1007RD(FILENAME)
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
      CHARACTER*255 BUFFER
      INTEGER IFIRST, IFILE, I,J,K,L,M,N,   NBIN,NX
      INCLUDE 'fnt1007.inc'

      DATA IFIRST/0/
      SAVE IFIRST

      OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=IFILE)
      IF (IFILE .ne. 0) THEN
         WRITE(*,*) '          fastNLO:  table file not found ',
     +        '  -  IOSTAT = ',IFILE
         STOP
      ENDIF

      READ(2,*) i
      if (i.ne.iseparator) goto 999
      READ(2,*) ITABVERSION
      WRITE(*,*) "#       tableformat is version",real(itabversion)/10000d0
      if (ITABVERSION.ne.14000) then
         write(*,*)"#     ==> this usercode works only for version 1.4"
         write(*,*)"#     ==> please get updated usercode from"
         write(*,*)"#         http://hepforge.cedar.ac.uk/fastnlo "
         stop
      endif
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c  -----------------------------------
      READ(2,*) IREACTION
      READ(2,*) ECMS
      READ(2,*) IXSECTUNITS
      do i=1,5
         READ(2,*) NAMELABEL(i)
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
         write(*,5000) ' #      ',NEVT(i),
     +        ' events in ',POWLABEL(i)
      enddo
      READ(2,*) NXTOT
      WRITE(*,*) "#           No. of x bins: ",NXTOT
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
      write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      write(*,*) " >>>>>   fastNLO error in table format "
      write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      stop
      RETURN

 5000 FORMAT (A,I12,A,A64) 
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
