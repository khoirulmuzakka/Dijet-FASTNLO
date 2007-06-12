*******************************************************************
* M. Wobisch 04/18/05          fastnlo00-hbook.f
*
* this package contains the tools to create and fill HBOOK histogram
* files after the fastNLO user routine has been called at least once
* to
* to do: book x-histos not in fixed numbers (e.g.:1-25)
*        but in h(x) using info of xlimit
*        -> easier to print later in PAW ith axis labels
*
*
* included::
*      SUBROUTINE FNHIST(filename)
*      SUBROUTINE FNHFILL(nscale)
*      SUBROUTINE FNHBOOK(n,filename,m)
*
*******************************************************************
      SUBROUTINE FNHIST(tablefile,histofile,ix)
*-----------------------------------------------------------------
* MW 04/18/2005
*
* book & fill HBOOK histograms
* requires that the corresponding tableinformation is already in
* the commonblock - from a previous call of FX9999CC
*
* Input:
*    tablefile  name of table-file
*    histofile  name of HBOOK-file
*    ix         flag:  1: only x-sections no x-histos
*                      2: also x-histos
*-----------------------------------------------------------------
      IMPLICIT NONE
      Character*(*) tablefile, histofile
      integer ix, nbin, i,j,l,m
      double precision mur,muf
      INCLUDE 'fnx9999.inc'

      write(*,*) '   #     store HBOOK histograms in file >>',histofile,'<<'
 
c - reset result arrays
      nbin = 0                  ! continuos numbering for the final array
      do i=1,nrapidity          ! (Pseudo-)Rapidity Bins
         do j=1,NPT(i)          ! ET/pT Bins
            nbin=nbin+1
            do l=1,Nord         ! all orders LO, NLO, NNLO, ...
               do m=1,(Nsubproc+1) ! No.of Subproc + 1 for sigma-tot
                  result(nbin,m,l)=0d0 ! reset 
               enddo
               xsect(nbin,l)=0d0 ! reset 
            enddo
         enddo
      enddo

c - open Hbook file & book histograms
      write(*,*)'FNHIST: booking histograms'
      call FNHBOOK(1,histofile,ix) 
c - loop over all scale variations
      do i=1,nscalevar
         mur = murscale(i)
         muf = mufscale(i)
         call FX9999CC(TABLEFILE, mur, muf, 0, XSECT)
c --- fill histograms
         write(*,*)'FNHIST: filling histograms for scale No. ',i,' mur/muf=',mur,muf
         call FNHFILL(i)
      enddo

c - fill x-histos for all scales
      if (ix.eq.2) call FNHFILLX

c - and close Hbook file
      call FNHBOOK(2,histofile,ix)  

      RETURN
      END
C -------------------------------------------------------------

      SUBROUTINE FNHFILL(nscale)
*-----------------------------------------------------------------
* MW 02/08/2006 make consistent with commonblock in v1.4
* MW 12/15/2005 add histograms for higher orders
* MW 08/26/2005 change to new array structure
* MW 04/18/2005
*   
*  input:  nscale   No of scale choice where results are to be stored
*
* fill HBOOK histograms
*-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER nscale,   I,J,K,L, nbin,nx
      Integer iord,isub,isub2, iscale,ihist
      INCLUDE 'fnx9999.inc'
      real value,v

c - identify the scale bin
      if (nscale.lt.1 .or. nscale.gt.Nscalevar) then
         write(*,*) " FNHFILL:   nscale=",nscale," is out of range"
         write(*,*) "            max: ",nscalevar
         stop
      endif
      iscale = nscale
         
c - fill all histograms for the given scale
      do iord=0,Nord            ! order: tot, LO, NLO-corr, 3 NNLOcorr
         do isub2=1,(Nsubproc+1) ! subprocess: Nsubproc + 1 tot
            isub=isub2
            if (isub.eq.(Nsubproc+1)) isub=0
            nbin=0
            do i=1,nrapidity                   
               do j=1,npt(i)
                  nbin = nbin + 1
                  if (iord.gt.0) then
                     value = 0d0
                     do k=1,iord
                        value = value+real(result(nbin,isub2,k))
                     enddo
                  else
                     value = 0d0
                     do k=1,Nord
                        value = value+real(result(nbin,isub2,k))
                     enddo
                  endif
c - avoid huge numbers - outside real precision range
                  if (abs(value).lt.1E38) then
                     continue
                  else          ! large number / or NAN / or INF
                     write(*,*) 'fn-hbook: found extreme value ',value,
     +                    ' in bin ',i,j
                     value = -1E20
                  endif
                  ihist = iord*1000000+iscale*100000+isub*10000+i*100
                  call hfill(ihist, real(PTBIN(i,j)+0.01) ,0.0, value)
               enddo            ! pT-loop
            enddo               ! rap-loop
         enddo                  ! isub loop
      enddo                     ! iord-loop

      RETURN
      END

C -------------------------------------------------------------

      SUBROUTINE FNHFILLX
*-----------------------------------------------------------------
* MW 02/08/2006 -> x-histos are for all orders / for a scales / for v1.4
* MW 12/15/2005 -> x-histos are only for NLO (no higher orders yet)
* MW 08/26/2005 change to new array structure
* MW 04/18/2005
*   
* fill HBOOK histograms
*-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER I,J,K,L,P,iord,isub, iscale,  nx2limit    
      Integer  nbin,nx , index ,ihist,ihi2,ihi3,ihi4
      INCLUDE 'fnx9999.inc'
      double precision FNALPHAS, as(nscalebinmax), mur, contr, sum, hxlimit,hx1,hx2
      DOUBLE PRECISION mu0scale,llptlo,llpthi,t1,t2,bweight
      real v
      PARAMETER (mu0scale=0.25)

      WRITE(*,*) '          fastNLO: fill x-histograms'

c - if the x-histos should be extended to different scales, then first
c   the booking needs to be extended .... in FNHBOOK
      do iscale=1,Nscalevar
c - first get PDFs at right scale
         call FX9999GP(mufscale(iscale))
         nbin = 0               ! continuos numbering for the final array
         do i=1,nrapidity       ! (Pseudo-)Rapidity Bins
            do j=1,NPT(i)       ! ET/pT Bins
               nbin=nbin+1                           
               llptlo = log(log(murscale(iscale) * murval(i,j,1)/mu0scale))
               llpthi = log(log(murscale(iscale) * murval(i,j,NSCALEBIN)/mu0scale))
            
               sum = 0d0
               
               nx=0
               do k=1,NXTOT     ! Nxmax
                  nx2limit = k
                  if (ireaction.eq.1) nx2limit=1 ! for DIS integrate only one x
                  do l=1,nx2limit ! Nxmin (for DIS: no loop)
                     nx=nx+1
                     do isub=1,Nsubproc ! Subprocess

                        if(NSCALEBIN.eq.1) then
                           as(1) =  FNALPHAS(murval(i,j,1) * murscale(iscale))
                        else
                           do p=1,NSCALEBIN ! loop over all scale bins, calc alphas
                              mur = mu0scale * exp(exp((llptlo + (dble(l)-1.)/(NSCALEBIN-1)*(llpthi-llptlo) )))
                              as(p) = FNALPHAS(mur) !    ... and get alpha_s
                           enddo
                        endif
                        do p=1,NSCALEBIN
                           if(nscalebin.eq.1) then
                              bweight =  as(1)
                           elseif(nscalebin.eq.2) then
                              if(p.eq.1) bweight =  as(1)
                              if(p.eq.2) bweight =  as(2)
                           elseif(nscalebin.eq.3) then 
                              t1 = 1./2.
                              if(p.eq.1) bweight =  as(1)
                              if(p.eq.2) bweight = 1./(2.*(1.-t1)*t1)*
     &                             (as(2)-as(1)*(1.-t1)**2-as(3)*(t1)**2)
                              if(p.eq.3) bweight = as(3)
                           elseif(nscalebin.eq.4) then 
                              t1 = 1./3.
                              t2 = 2./3.
                              if(p.eq.1) bweight =  as(1)
                              if(p.eq.2) bweight = 1./(3.*(1.-t1)**2*t1*3.*(1.-t2)*(t2)**2 - 3.*(1.-t2)**2*t2 * 3.*(1.-t1)*(t1)**2)
     &                             *(3.*(1.-t2)*(t2)**2*(as(2)-as(1)*(1.-t1)**3-as(4)*(t1)**3) 
     &                             - 3.*(1.-t1)*(t1)**2*(as(3)-as(1)*(1.-t2)**3-as(4)*(t2)**3  ) )
                              if(p.eq.3) bweight =1./(3.*(1.-t1)**2*t1*3.*(1.-t2)*(t2)**2  - 3.*(1.-t2)**2*t2 * 3.*(1.-t1)*(t1)**2)
     &                             *(3.*(1.-t1)**2*(t1)*(as(3)-as(1)*(1.-t2)**3-as(4)*(t2)**3) 
     &                             - 3.*(1.-t2)**2*(t2)*(as(2)-as(1)*(1.-t1)**3-as(4)*(t1)**3  ) )
                              if(p.eq.4) bweight = as(4)
                           else
                              write(*,*)'NSCALEBIN = ',nscalebin,' not yet supported.'
                              stop
                           endif

                           do iord=1,Nord ! order: 1 LO  2 NLO 3 NNLO
                              ihist= iord*1000000+iscale*100000+isub*10000+i*100+j
                              ihi2 = iord*1000000+iscale*100000 +i*100+j ! sum subproc
                              ihi3 = iscale*100000+isub*10000+i*100+j ! sum orders
                              ihi4 = iscale*100000+i*100+j ! sum orders and subpr.

                              if (iord.eq.1) index=1
                              if (iord.eq.2) index=iscale+1
                              if (iord.eq.3) index=Nscalevar+iscale+1

                              contr = array(nbin,nx,isub,index,p)
     +                             * bweight**(iord+1) ! multiply alpha-s
     +                             * pdf(nbin,nx,isub,p) ! multiply PDFs
                              sum=sum+contr
c     - - - - - - - - - - see below for test-output for: sum

                              v = real(contr)
                              if (ixscheme.eq.2) then ! default
                                 hxlimit = - sqrt(log10(1D0/xlimit(i,j)))
                              elseif (ixscheme.eq.1) then ! simple log
                                 hxlimit = log10(xlimit(i,j))
                              elseif (ixscheme.eq.3) then ! proposed log(x) + x -1
                                 hxlimit = log10(xlimit(i,j))+xlimit(i,j)-1d0
                              else
                                 write(*,*) 'FNHBOOK: ixscheme not implemented: ',ixscheme
                              endif
                              hx1 = hxlimit*(1d0 - dble(k-0.5d0)/dble(nxtot))
                              hx2 = hxlimit*(1d0 - dble(l-0.5d0)/dble(nxtot))
                              call hfill( ihist, real(hx1) ,0.0, v)
                              call hfill(-ihist, real(hx2) ,0.0, v)
                              call hfill( ihi2 , real(hx1) ,0.0, v)
                              call hfill(-ihi2 , real(hx2) ,0.0, v)
                              call hfill( ihi3 , real(hx1) ,0.0, v)
                              call hfill(-ihi3 , real(hx2) ,0.0, v)
                              call hfill( ihi4 , real(hx1) ,0.0, v)
                              call hfill(-ihi4 , real(hx2) ,0.0, v)

                           enddo ! iord
                        enddo   ! scalebin
                     enddo      ! subproc
                  enddo         ! xmin
               enddo            ! xmax
c - test that we get the same integrated cross section as in original calculation
c               if (iscale.eq.3) write(*,*) " rapbin ",i,"  ptbin ",j," -> ",sum
            enddo               ! npt
         enddo                  ! nrap
      enddo                     ! iscale

      RETURN
      END

*********************************************************************

      SUBROUTINE FNHBOOK(n,filename,ix)
*-----------------------------------------------------------------
* MW 04/18/2005
*
* n=1  create HBOOK file - book HBOOK histograms
* n=2  close HBOOK file
*
* filename: name of histogram file
*
* for n=1 only:
* ix=1 create only pT histograms
* ix=2 create also x-histograms
*-----------------------------------------------------------------
      IMPLICIT NONE
      Character*(*) filename
      INTEGER IFIRST, IFILE, J,N,istat2,icycle
      Integer ix,iord,isub, iscale, irap, ihist
      Double Precision hxlimit
      INCLUDE 'fnx9999.inc'
      real pt(nptmax)

c - HBOOK common
      INTEGER NWPAWC
c      PARAMETER (NWPAWC=2500000)
      PARAMETER (NWPAWC=15000000) ! for very large scenarios w/ x-bins
      REAL HMEMOR(NWPAWC)
      COMMON /PAWC/ HMEMOR

      write(*,*) "       FNHBOOK  n=",n,"   ord=",nord

c - open & book
      if (n.eq.1) then
         CALL HLIMIT(NWPAWC)
         CALL HROPEN(11,'fastNLO',filename,'N',1024,ISTAT2)
         if (ISTAT2.NE.0) then
            WRITE(*,*) ' FNHBOOK: could not open histofile ',istat2
         endif

         do iord=0,Nord         ! order: tot, LO, NLO-corr, NNLO-corr -> Nord
            do iscale=1,NSCALEVAR ! scale variations
               do isub=0,Nsubproc      ! subprocess: 0 tot + 7 subproc
               
                  do irap=1, nrapidity
                     do j=1,(npt(irap)+1)

                        pt(j) = real(PTBIN(irap,j))
c                        write(*,*) "  irap,j=",irap,j,"  pt=",pt(j)
                        ihist = iord*1000000+iscale*100000+isub*10000+irap*100

c --------------------- x-distributions in pT-bins 
                        if (ix.eq.2 .and.j.ne.(npt(irap)+1)) then
                           if (ixscheme.eq.2) then    ! default
                              hxlimit = - sqrt(log10(1D0/xlimit(irap,j)))
                           elseif (ixscheme.eq.1) then   ! simple log
                              hxlimit = log10(xlimit(irap,j))
                           elseif (ixscheme.eq.3) then   ! proposed log(x) + x -1
                              hxlimit = log10(xlimit(irap,j))+xlimit(irap,j)-1d0
                           else
                              write(*,*) 'FNHBOOK: ixscheme not implemented: ',ixscheme
                           endif
c                           write(*,*) "xlim ",irap,j,xlimit(irap,j),hxlimit
                           call hbook1(ihist+j,'x-Bin',nxtot,real(hxlimit),0.0,0)
                           call hbook1(-(ihist+j),'x-Bin',nxtot,real(hxlimit),0.0,0)
                        endif

                     enddo
                     call hbookb(ihist,'p?T! (GeV)' , NPT(irap) ,PT ,0)
                  enddo
               enddo
            enddo
         enddo

c - close
      elseif (n.eq.2) then
         CALL HROUT (0,ICYCLE,' ')
         CALL HREND ('fastNLO')
      endif

      RETURN
      END


