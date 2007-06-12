
* -------------------------------------------------------------------
* fastNLO user package                     M. Wobisch 06/29/2005
*
* fnx9999.f
*
* contains the following routines
*     FX9999CC        computes cross sections / fills arrays
*     FX9999RD        reads coefficient tables
*     FX9999PR        prints x-section results 
*     FX9999GP        gets the PDF values for all x-bins
*     FX9999MT        multiply PDF array and coefficient array
*
* uses commonblock definitions in: fnx9999.inc
*
* needs the following routines from fn-pdf.f 
*     FNPAPDF         computes the 7 PDF linear combinations
*                         for "PA": proton-antiproton
*     FNPPPDF         computes the 7 PDF linear combinations
*                         for "PP": proton-proton
*     FNEPPDF         computes the 3 PDF linear combinations
*                         for "EP": electron-proton
*
* needs the following routines from "fnpdfalpha.f":
*     FNALPHAS        alpha_s interface (double precision function)
*     FNPDF           PDF interface
*
* question: change 2-d x-Array to 1-dim (save factor 2 in memory)?
*          - also store PDFs in 'parallel 1d array
* -------------------------------------------------------------------


*******************************************************************
      SUBROUTINE FX9999CC(FILENAME, XMUR, IMUFFLAG, IPRINTFLAG, XSECT)
*-----------------------------------------------------------------
* MW 09/02/2005 implement flexible scale variations
* TK 12/07/2005 table format contains now LO + NLO with 5 scale variations, 
* support for zipped files
* MW 04/15/2005
*
* compute cross section
*
* input:
*   FILENAME     name of input table
*   XMUR        negative values: multiply fact. scale:        mur=xmur*muf
*               positive values: multiply default ren. scale: mur=xmur*mur0
*   IMUFFLAG     flag to enable/disable scale variations
*                 0: compute x-sect only for central scale
*                 1: compute x-sect for all different scales
*                    (takes more time for additional PDF access)
* output:
*   XSECT(nbin,2,nscalevar)    array of cross sections 
*                              first dimension: Bin number
*                              second dim: 1 LO   2 NLO
*                              third dim: 1 default   2....nscalevar varaiations
*-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INTEGER IFIRST, IFILE, iord,iscale, I,J,K,L,M, IMUFFLAG, IPRINTFLAG, 
     +     maxscale, nbin,nx
      CHARACTER*(*) FILENAME
      CHARACTER*50 OLDFILENAME
      DOUBLE PRECISION xmur, xmur1, binwidth, renscale
      DATA IFIRST/0/, OLDFILENAME/'xxxx'/
      SAVE IFIRST, OLDFILENAME

c - output in first call
      if (IFIRST.EQ.0) then
         write(*,*) " "
         write(*,*) "  ##################################################",
     +        "#############"
         write(*,*) "  #" 
         write(*,*) "  #   fastNLO - beta release            September 14, 2005"
         write(*,*) "  #" 
         write(*,*) "  #   Thomas Kluge, Klaus Rabbertz, Markus Wobisch" 
         write(*,*) "  #" 
         write(*,*) "  #   user code to compute the inclusive jet ",
     +        "cross section"
         write(*,*) "  #   for the Run I measurements at the Tevatron by the"
         write(*,*) "  #   D0 and CDF collaborations:" 
         write(*,*) "  #      ",
     +        "B. Abbott et al. (D0 Collaboration), "
         write(*,*) "  #         ",
     +        "          Phys. Rev. Lett. {86} 1707 (2001)."
         write(*,*) "  #      ",
     +        "T. Affolder et al. (CDF Collaboration), "
         write(*,*) "  #         ",
     +        "          Phys. Rev. D64, 032001 (2001)."
         write(*,*) "  #  " 
         write(*,*) "  #  " 
c         write(*,*) "  #  if you use this code, please cite the ",
c     +        "following reference:" 
         write(*,*) "  #  !!!!  The purpose of this beta-release is to debug" 
         write(*,*) "  #  !!!!  the code and to collect suggestions for the" 
         write(*,*) "  #  !!!!  user-interface" 
         write(*,*) "  #  !!!!  Results obtained using this code should not" 
         write(*,*) "  #  !!!!  be used in publications!!!" 
         write(*,*) "  #        " 
         write(*,*) "  #  " 
         write(*,*) "  #  the perturbative coefficients have been computed" 
         write(*,*) "  #  for LO and NLO using NLOJET++ by Zoltan Nagy ",
     +        "- please cite:" 
         write(*,*) "  #      ",
     +        "Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),"
         write(*,*) "  #      ",
     +        "Z. Nagy, Phys. Rev. D68, 094002 (2003)."
         write(*,*) "  #  " 
         write(*,*) "  #  the NNLO-NLL threshold corrections have been computed"
         write(*,*) "  #  by N. Kidonakis and J.F. Owens - please cite:" 
         write(*,*) "  #      N. Kidonakis, J.F. Owens, ",
     +        "Phys. Rev. D63, 054019 (2001)." ! hep-ph/0007268
c   EFFECTS OF HIGHER ORDER THRESHOLD CORRECTIONS IN HIGH E(T) JET PRODUCTION.
         write(*,*) "  #  " 
      endif

c --- read the fastNLO coefficient tables
      IF (FILENAME.ne.OLDFILENAME) THEN
c         write (*,*) ' old/new table ',OLDFILENAME,'   ',FILENAME
         call FX9999RD(FILENAME)  ! read combined LO & NLO table
         OLDFILENAME=FILENAME

c - make consistency checks of array dimensions / commonblock parameters
         if (((NXTOT**2+NXTOT)/2).gt.NXMAX) goto 999
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
         if (k.gt.NRAPPT) goto 999

c - print further info
         write(*,*) "  #"
         write(*,*) "  #      reaction: ",cireaction(ireaction)
         write(*,*) "  #       process: ",ciproc(iproc)
         write(*,*) "  #      jet algo: ",cialgo(ialgo)
         write(*,*) "  #         parameter 1:  ",cjetres1(ialgo),'=',jetres1
         write(*,*) "  #         parameter 2:  ",cjetres2(ialgo),'=',jetres2
         write(*,*) "  #"
c - print scale-variations used in the table


         write (*,*) "  #   --- available scale variations:",
     +        nscalevar
         write (*,*) "  #         factorization scale variations:"
         do i=1,nscalevar
            write (*,*) "  #           ",i,"  (muf/ET)",mufscale(i)
         enddo
         write (*,*) "  #         default renormalization scale settings:"
         do i=1,nscalevar
            write (*,*) "  #           ",i,"  (mur/ET)=",murscale(i)
         enddo
         write (*,*) "  #     ->  your chosen renormalization scales in 1st call:"
         do i=1,nscalevar
            if (xmur.lt.0d0) then
               renscale = -xmur * mufscale(i)
            elseif (xmur.gt.0d0) then
               renscale = xmur * murscale(i)
            endif
c            write (*,*) " xmur,renscale ",xmur,renscale, mufscale(i),murscale(i)
            write (*,*) "  #           ",i," (mur/ET)=",renscale,
     +           "   (muf/ET)",mufscale(i)
         enddo
      ENDIF


      if (IFIRST.EQ.0) then
         IFIRST = 1
         write(*,*) "  #     -- coefficient tables are now in memory --"
         write(*,*) "  #################################################",
     +        "#############"
         write(*,*) " "
      endif


c - reset result arrays
      nbin = 0                  ! continuos numbering for the final array
      do i=1,nrapidity          ! (Pseudo-)Rapidity Bins
         do j=1,NPT(i)          ! ET/pT Bins
            nbin=nbin+1
            do k=1,NSCALEVAR    ! central + scale variations
               do l=1,Nord      ! all orders LO, NLO, NNLO, ...
                  do m=1,8      ! 7 Subproc + total
                     result(nbin,m,l,k)=0d0 ! reset 
                  enddo
                  xsect(nbin,l,k)=0d0 ! reset 
               enddo
            enddo
         enddo
      enddo


c - compute xsect - for different scales / if selected
      if (imufflag.eq.0) then
         maxscale = 1
      elseif (imufflag.eq.1) then
         maxscale = nscalevar
      else
         write(*,*) '     you can''t select   imufflag=',imufflag
         stop
      endif 
      do i=1,maxscale
         call FX9999GP(mufscale(i))
         call FX9999MT(xmur,i)
      enddo


c ----------------- final touches ------------------------------------
c - sum up subprocesses / divide by pT Bins width / fill result array
c -                   multiply by 1000 to get pb / fill array xsect
      nbin=0
      do i=1,nrapidity          ! Rapidity Bins
         do j=1,NPT(i)          ! pT Bins
            nbin=nbin+1
c - divide by binwidth / factor "2" because we use: abs(eta)
            IF (RAPBIN(i+1).gt.RAPBIN(i)) then
               binwidth=(PTBIN(i,j+1)-PTBIN(i,j))*(RAPBIN(i+1)-RAPBIN(i))*2d0
            else                ! needed for CDF data in D0 array
               binwidth=(PTBIN(i,j+1)-PTBIN(i,j))*(0.6d0)*2d0
            endif
            do iscale=1,NSCALEVAR
               do iord=1,Nord
                  result(nbin,8,iord,iscale) = 0d0
                  xsect(nbin,iord,iscale) = 0d0
                  do m=1,7      ! Subprocess
                     result(nbin,m,iord,iscale) = result(nbin,m,iord,iscale)
     +                    *1000d0/binwidth
                     result(nbin,8,iord,iscale) = result(nbin,8,iord,iscale)
     +                    + result(nbin,m,iord,iscale)
                  enddo
                  xsect(nbin,iord,iscale)=result(nbin,8,iord,iscale)
               enddo
            enddo
         enddo
      enddo

c - compute xsect for different scales
      if (IPRINTFLAG.eq.1) call FX9999PR(IMUFFLAG)

      RETURN

 999  continue
      WRITE(*,*) 'fastNLO: error occured - parameter outside range - check commonblock'
      STOP
      END

*******************************************************************
      SUBROUTINE FX9999MT(xmur,iscale)
*-----------------------------------------------------------------
* MW 08/26/2005
*
* multiply the PDFs and the coefficients for a single scale setting
* 
* input:    XMUR variation of renormalization scale
*         ISCALE = 1,..,Nscalevar: compute for different fact. scale
*                                (as stored in table) 
*
*   XMUR     negative values: multiply fact. scale:    mur=abs(xmur)*muf
*            positive values: multiply default ren. scale: mur=xmur*mur0
*-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INTEGER ISCALE, i,j,k,l,m,iord, jord,    nbin,nx
      DOUBLE PRECISION  mur, as, FNALPHAS, coeff , xmur1
      DOUBLE PRECISION beta02pi, PI, xmur, scfac ! for ren-scale variation
      INTEGER NF, NC
      PARAMETER (PI=3.14159265358979323846, NF=5, NC=3)
      PARAMETER (beta02pi=(11*NC-2*NF)/6) ! = (11*NC-2*NF)/12/PI*2*PI

c - set absolute order in alpha_s of LO process
      jord = NPOW(1)

c - evaluate real xmur (=> xmur1)
      if (xmur.lt.0d0) then
         xmur1 = -xmur * mufscale(iscale)/murscale(iscale)
      elseif (xmur.gt.0d0) then
         xmur1 = xmur
      endif
               
c - vary renormalization scale around original vaule
      scfac = dble(jord) *beta02pi *2d0*log(xmur1)

c - loop over array
      nbin = 0                  ! continuos numbering for the final array
      do i=1,nrapidity          ! (Pseudo-)Rapidity Bins
         do j=1,NPT(i)          ! ET/pT Bins
            nbin=nbin+1            

            mur = xmur1 * murval(i,j) * murscale(iscale) ! set scale
            as = FNALPHAS(mur)  ! get alpha_s

            nx=0
            do k=1,NXTOT         ! Nxmax
               do l=1,k         ! Nxmin 
                  nx=nx+1
                  do m=1,7      ! Subprocess
                     do iord=1,Nord ! relative order: 1 LO  2 NLO  3 NNLO ...
                        if (iord.eq.1) then ! LO contribution
                           coeff = array(nbin,nx,m,1)
                        elseif (iord.eq.2) then ! NLO contributions
                           coeff = array(nbin,nx,m,
     +                          (1+iscale+(iord-2)*nscalevar)) +
     +                          scfac*array(nbin,nx,m,1) 
                        elseif (iord.eq.3) then ! NNLO-NLL threshold corr.
c - ren. scale variation works not (yet?) for NNLO-NLL(threshold-correction)
c   -> the NLLO-NLL contributions are only available for  mu_r = mu_f
c -           in other words: for  xmur=-1
                           if (XMUR.eq.-1d0) then
                              coeff = array(nbin,nx,m,
     +                          (1+iscale+(iord-2)*nscalevar))
                           else
                              coeff = 0d0
                           endif
c                           write(*,*) i,j,k,l,m,iord,coeff
                        endif
                        result(nbin,m,iord,iscale)=result(nbin,m,iord,iscale)
     +                       + coeff
     +                       * as**(iord+jord-1) ! multiply with (alpha-s/2pi)
     +                       * pdf(nbin,nx,m) ! multiply with PDFs
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      RETURN
      END

*******************************************************************
      SUBROUTINE FX9999GP(muffactor)
*-----------------------------------------------------------------
* MW 08/26/2005
*
* read the PDFs in all (rap,pT) bins - at all xmax, xmin bins
* the default factorization scale (in GeV) is multiplied by muffactor
*-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      DOUBLE PRECISION MUFFACTOR
      INTEGER NBIN,NX, i,j,k,l,m
      DOUBLE PRECISION x1, x2, xlim, hx, hxlim, muf, H(7), sum

      nbin=0
      do i=1,nrapidity          ! Rapidity Bins
         do j=1,NPT(i)          ! pT Bins
            nbin=nbin+1

            xlim = xlimit(i,j)
            hxlim = - sqrt(-log10(xlim))
            muf = mufval(i,j) * muffactor
c           write(*,*) '>>> ',i,j, muffactor,mufval(i,j)
            nx=0 
            do k=1,NXTOT         ! Nxmax
               hx = hxlim *(1d0 - dble(k-1)/dble(nxtot)) ! compute x1-value  
               x1 = 10**-(hx*hx)
               do l=1,k         ! Nxmin 
                  hx = hxlim *(1d0 - dble(l-1)/dble(nxtot)) ! compute x2
                  x2 = 10**-(hx*hx)

                  nx=nx+1

c * save time by reducing LHAPDF access - no access if coeff's are zero
c - need to check both LO and NLO array elements (and NNLO if available)
c - assume: if zero at one scale - then zero at all scales (trivial!)
c - currently we check the contributions from all subprocesses
c -    -> probably not needed - but does not hurt CPU usage
                  sum = 0d0
                  do m=1,7
                     sum=sum+abs(array(nbin,nx,m,1))+abs(array(nbin,nx,m,2))
                     if (nord.gt.2) sum=sum+abs(array(nbin,nx,m,(2+nscalevar)))
                  enddo
                  if (sum.eq.0d0) then
                     do m=1,7
                        pdf(nbin,nx,m) = 0d0
                     enddo
                  else

c -     get PDFs at x1,x2 at scale muf
                     if (ireaction.eq.2) then
                        call FNPPPDF(x1,x2,muf,H) ! for proton-proton
                     elseif (ireaction.eq.3) then
                        call FNPAPDF(x1,x2,muf,H) ! for proton-antiproton
                     else
                        write(*,*) '    ireaction =',ireaction
                        write(*,*) ' this is not defined in hadron collisions'
                        stop
                     endif
                     do m=1,7   ! Subprocess
                        pdf(nbin,nx,m) = H(m)
     +                       / (1d0-0.99d0*x1)**3/(1d0-0.99d0*x2)**3
     +                       * sqrt(x1*x2) ! 2 lines PDF reweighting
                     enddo
                  endif         ! end: check if coeff's are zero
               enddo
            enddo
         enddo
      enddo

      RETURN 
      END

*******************************************************************
      SUBROUTINE FX9999PR(IMUFFLAG)
*-----------------------------------------------------------------
* MW 04/18/2005
*
* prints cross section results
*  IMUFFLAG    0: print only results for central scale
*              1: print also results for scale variations 
*-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IMUFFLAG, iscale, nbin, i,j, maxscale
      INCLUDE 'fnx9999.inc'

      write(*,*) "      fastNLO - results (in pb):"
      write(*,*) "         cross sections:    LO       NLOcorr     NNLOcorr",
     +     "      tot"
      write(*,*) " "

      if (IMUFFLAG.eq.1) then
         maxscale = NSCALEVAR
      else
         maxscale = 1
      endif

      do iscale=1,maxscale
         write(*,*) "  ----------- scale variation No.",iscale
         
         nbin = 0               ! linear bin for the final observable
         do i=1,nrapidity       ! Rapidity Bins
            if (RAPBIN(i).lt.RAPBIN(i+1)) then
               WRITE(*,*) "   rapidities from ",RAPBIN(i)," - ",RAPBIN(i+1),
     +              "         No pT bins: ",NPT(i)
            else                !  happens only for CDF bins
               WRITE(*,*) "   rapidities from ",0.1d0," - ",RAPBIN(i+1),
     +              "         No pT bins: ",NPT(i)
            endif
            do j=1,NPT(i)       ! pT Bins
               nbin=nbin+1
               WRITE(*,900) PTBIN(i,j),PTBIN(i,j+1),
     +              real(result(nbin,8,1,iscale)),
     +              real(result(nbin,8,2,iscale)),
     +              real(result(nbin,8,3,iscale)),
     +              real(result(nbin,8,1,iscale)
     +              + result(nbin,8,2,iscale)
     +              + result(nbin,8,3,iscale))
            enddo
         enddo
         write(*,*) " "
      enddo
      RETURN
 900  Format ("pT:",F8.2,"-",F8.2,":",4E13.4)
      END

*******************************************************************
      SUBROUTINE FX9999RD(FILENAME)
*-----------------------------------------------------------------
* MW 04/15/2005
*
* read ASCII table
* iord: order    1:LO  2:NLO table
*-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) FILENAME
      CHARACTER*255 BUFFER
      INTEGER IFIRST, IFILE, I,J,K,L,M,N,   NBIN,NX
      INTEGER FLNLEN,LENOCC
      LOGICAL ZIPPED
      INCLUDE 'fnx9999.inc'

      DATA IFIRST/0/
      SAVE IFIRST

      FLNLEN = LENOCC (FILENAME)
      IF (FILENAME(FLNLEN-2:FLNLEN) .EQ. '.gz')  THEN
         ZIPPED = .TRUE.
         write(*,*) '  #      reading the table (in zipped format)' 
      ELSE
         ZIPPED = .FALSE.
         write(*,*) '  #      reading the table (in ASCII format)' 
      ENDIF

      IF(ZIPPED) THEN
         CALL GZOPEN(2,'R' , FILENAME, IFILE)
      ELSE
      OPEN(2,STATUS='OLD',
     +        FILE=FILENAME,IOSTAT=IFILE)
      ENDIF

      IF (IFILE .ne. 0) THEN
         WRITE(*,*) '          fastNLO:  table file not found ',
     +        '  -  IOSTAT = ',IFILE
         STOP
      ENDIF

c      WRITE(*,*) "  #      No. max bins (eta,ET,x): ",NRAPMAX,NPTMAX,NXMAX
c     -----------------------------------
      IF(ZIPPED) THEN
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) IREACTION
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) ECMS
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) IPROC
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) IALGO
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) JETRES1
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) JETRES2
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) NORD
         do i=1,NORD
            CALL GZREAD(2,BUFFER,255,IFILE)
            READ(BUFFER,*) NPOW(i)
         enddo
         do i=1,NORD
            CALL GZREAD(2,BUFFER,255,IFILE)
            READ(BUFFER,*) POWLABEL(i)
         enddo
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,NORD
            CALL GZREAD(2,BUFFER,255,IFILE)
            READ(BUFFER,*) NEVT(i)
         enddo
         do i=1,nord
            write(*,5000) '   #      ',NEVT(i),
     +           ' events in ',POWLABEL(i)
         enddo
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) NXTOT
         WRITE(*,*) "  #       No. of x bins: ",NXTOT
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) IXSCHEME
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) IPDFWGT
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) NRAPIDITY
         do i=1,nrapidity+1
            CALL GZREAD(2,BUFFER,255, IFILE)
            READ(BUFFER,*) RAPBIN(i)
         enddo
         do i=1,nrapidity
            CALL GZREAD(2,BUFFER,255, IFILE)
            READ(BUFFER,*) NPT(i)
         enddo
         do i=1,nrapidity
            do j=1,NPT(i)+1
               CALL GZREAD(2,BUFFER,255, IFILE)
               READ(BUFFER,*) PTBIN(i,j)
            enddo
         enddo
         CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,nrapidity
            do j=1,NPT(i)
               CALL GZREAD(2,BUFFER,255, IFILE)
               READ(BUFFER,"(F20.17)") XLIMIT(i,j)
            enddo
         enddo
         CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,nrapidity
            do j=1,NPT(i)
               CALL GZREAD(2,BUFFER,255, IFILE)
               READ(BUFFER,*) MURVAL(i,j)
            enddo
         enddo
               CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,nrapidity
            do j=1,NPT(i)
               CALL GZREAD(2,BUFFER,255, IFILE)
               READ(BUFFER,*) MUFVAL(i,j)
            enddo
         enddo
         CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) NSCALEVAR
         do i=1,NSCALEVAR
            CALL GZREAD(2,BUFFER,255, IFILE)
            READ(BUFFER,*) MURSCALE(i)
         enddo
         do i=1,NSCALEVAR
            CALL GZREAD(2,BUFFER,255, IFILE)
            READ(BUFFER,*) MUFSCALE(i)
         enddo

         CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         nbin=0
         do i=1,nrapidity
            do j=1,NPT(i)
               nbin=nbin+1      ! linear numbering of rap,pT bins
               nx=0
               do k=1,NXTOT      ! Nxmax
                  do l=1,k      ! Nxmin 
                     nx=nx+1    ! linear numbering of xmax,xmin bins
                     do m=1,7   ! Subproc
                        do n=1,1+NSCALEVAR*(NORD-1) ! LO & NLO & w/ scale var
                           CALL GZREAD(2,BUFFER,255, IFILE)
                           READ(BUFFER,*) array(nbin,nx,m,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
         CALL GZCLOSE(2,IFILE)
      ELSE
         READ(2,*) IREACTION
         READ(2,*) ECMS
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
         enddo
         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,nord
            READ(2,*) NEVT(i)
         enddo
         do i=1,nord
            write(*,5000) '   #      ',NEVT(i),
     +           ' events in ',POWLABEL(i)
         enddo
         READ(2,*) NXTOT
         WRITE(*,*) "  #      No. of x bins: ",NXTOT
         READ(2,*) IXSCHEME
         READ(2,*) IPDFWGT
         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
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
c     -----------------------------------
         do i=1,nrapidity
            do j=1,NPT(i)
               READ(2,*) XLIMIT(i,j)
            enddo
         enddo
         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,nrapidity
            do j=1,NPT(i)
               READ(2,*) MURVAL(i,j)
            enddo
         enddo
         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,nrapidity
            do j=1,NPT(i)
               READ(2,*) MUFVAL(i,j)
            enddo
         enddo
         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         READ(2,*) NSCALEVAR
         do i=1,NSCALEVAR
            READ(2,*) MURSCALE(i)
         enddo
         do i=1,NSCALEVAR
            READ(2,*) MUFSCALE(i)
         enddo

         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         nbin=0
         do i=1,nrapidity
            do j=1,NPT(i)
               nbin=nbin+1      ! linear numbering of rap,pT bins
               nx=0
               do k=1,NXTOT      ! Nxmax
                  do l=1,k      ! Nxmin
                     nx=nx+1    ! linear numbering of xmax,xmin bins 
                     do m=1,7   ! Subproc
c                        do n=1,1+NSCALEVAR ! LO &  NLO with  scale variations
                        do n=1,1+NSCALEVAR*(NORD-1) ! LO & NLO & w/ scale var
                           READ(2,*) array(nbin,nx,m,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         READ(2,*) i
         if (i.ne.iseparator) goto 999
         CLOSE(2)
      ENDIF

      goto 998 

      WRITE(*,*) "nxtot  ",NXTOT
      do i=1,nrapidity
         WRITE(*,*) i,"  ",NPT(i)
      enddo
      goto 998
      WRITE(*,*) IREACTION
      WRITE(*,*) ECMS
      WRITE(*,*) IPROC
      WRITE(*,*) IALGO
      WRITE(*,*) JETRES1, JETRES2
      WRITE(*,*) NORD
      do i=1,nord
         WRITE(*,*) NPOW(i)
      enddo
      do i=1,nord
         WRITE(*,*) NEVT(i)
      enddo
      WRITE(*,*) NXTOT
      WRITE(*,*) IXSCHEME
      WRITE(*,*) IPDFWGT
      WRITE(*,*) NRAPIDITY
      do i=1,nrapidity+1
         WRITE(*,*) RAPBIN(i)
      enddo
      do i=1,nrapidity
         WRITE(*,*) i,"  ",NPT(i)
      enddo
      do i=1,nrapidity
         do j=1,NPT(i)+1
            WRITE(*,*) PTBIN(i,j)
         enddo
      enddo
      do i=1,nrapidity
         do j=1,NPT(i)
            WRITE(*,*) XLIMIT(i,j)
         enddo
      enddo



      WRITE(*,*) 
      WRITE(*,*) 

 998  continue
c      write(*,*) '  *   table is in memory'
      IF(ZIPPED) THEN
         CALL GZCLOSE(2,IFILE)
      ELSE
         close (2) 
      ENDIF
      RETURN

 999  continue
      IF(ZIPPED) THEN
         CALL GZCLOSE(2,IFILE)
      ELSE
         close (2) 
      ENDIF
      write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      write(*,*) " >>>>>   fastNLO error in table format "
c      write(*,*) " >>>>>       separator at wrong place!"
      write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      RETURN

 5000    FORMAT (A,I12,A,A64)

      END

*******************************************************************
 




