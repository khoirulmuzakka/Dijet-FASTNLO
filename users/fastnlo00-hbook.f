*******************************************************************
* M. Wobisch 04/18/05          fastnlo00-hbook.f
*
* this package contains the tools to create HBOOK histogram files
* from the fastNLO user routine
*
* included::
*      SUBROUTINE FNHIST(filename)
*      SUBROUTINE FNHFILL
*      SUBROUTINE FNHBOOK(n,filename,m)
*
*******************************************************************
      SUBROUTINE FNHIST(filename,ix)
*-----------------------------------------------------------------
* MW 04/18/2005
*
* book & fill HBOOK histograms
*
* Input
*    filename   name of HBOOK-file
*    ix         flag:  0: only x-sections no x-histos
*                      1: also x-histos
*-----------------------------------------------------------------
      IMPLICIT NONE
      Character*(*) filename
      integer ix

      write(*,*) '   #     store HBOOK histograms in file >>',filename,'<<'
 
c - open Hbook file & book histograms
      call FNHBOOK(1,filename,ix) 
c - fill histograms
      call FNHFILL
      if (ix.eq.2) call FNHFILLX
c - and close Hbook file
      call FNHBOOK(2,filename,ix)  

      RETURN
      END
C -------------------------------------------------------------

      SUBROUTINE FNHFILL
*-----------------------------------------------------------------
* MW 08/26/2005 change to new array structure
* MW 04/18/2005
*   
* fill HBOOK histograms
*-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER I,J,K,L, nbin,nx
      Integer iord,isub,isub2, iscale,ihist
      INCLUDE 'fnx9999.inc'
      real value,v

      do iord=0,2               ! ,2            ! order: tot, LO, NLO-corr
c         write(*,*) '  fill histograms of order ',iord,'/3 '
         do iscale=1,NSCALEMAX  ! scale variations
            do isub2=1,8         ! subprocess: 8 tot + 7 subproc
               isub=isub2
               if (isub.eq.8) isub=0
               nbin=0
               do i=1,nrapidity                   
                  do j=1,npt(i)
                     nbin = nbin + 1
                     if (iord.gt.0) then
                        value = real(result(nbin,isub2,iord,iscale))
                     else
                        value = real(result(nbin,isub2,1,iscale)
     +                       +result(nbin,isub2,2,iscale))
                     endif
                     ihist = iord*1000000+iscale*100000+isub*10000+i*100
                     call hfill(ihist, real(PTBIN(i,j)+0.01) ,0.0, value)
                  enddo         ! pT-loop
               enddo            ! rap-loop
            enddo               ! isub loop
         enddo                  ! iscale loop
      enddo                     ! iord-loop

      RETURN
      END

C -------------------------------------------------------------

      SUBROUTINE FNHFILLX
*-----------------------------------------------------------------
* MW 08/26/2005 change to new array structure
* MW 04/18/2005
*   
* fill HBOOK histograms
*-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER I,J,K,L,iord,isub, iscale      
      Integer  nbin,nx , index ,ihist,ihi2,ihi3,ihi4
      double precision FNALPHAS, as, mur, contr, sum
      INCLUDE 'fnx9999.inc'
      real v


      WRITE(*,*) '          fastNLO: fill x-histograms - test-output:'
      WRITE(*,*) '              test-output: - all non-normalized cross sections'

c - if the x-histos should be extended to different scales, then first
c   the booking needs to be extended .... in FNHBOOK
      do iscale=1,1             ! nscalemax 
c - first get PDFs at right scale
         call FX9999GP(mufscale(1))
         nbin = 0               ! continuos numbering for the final array
         do i=1,nrapidity       ! (Pseudo-)Rapidity Bins
            do j=1,NPT(i)       ! ET/pT Bins
               nbin=nbin+1                           
            
               sum = 0d0
               mur = murval(i,j) * murscale(iscale)
               as = FNALPHAS(mur) ! get alpha_s
               
               nx=0
               do k=1,NTOT      ! Nxmax
                  do l=1,k      ! Nxmin 
                     nx=nx+1
                     do isub=1,7 ! Subprocess
                        do iord=1,2 ! order: 1 LO  2 NLO
                           ihist = iord*1000000+iscale*100000+isub*10000+i*100+j
                           ihi2  = iord*1000000+iscale*100000 +i*100+j ! sum subproc
                           ihi3  = iscale*100000+isub*10000+i*100+j ! sum orders
                           ihi4  = iscale*100000+i*100+j ! sum orders and subproc

                           if (iord.eq.1) index=1
                           if (iord.eq.2) index=iscale+1

                           contr = array(nbin,nx,isub,index)
     +                          * as**(iord+1) ! multiply alpha-s
     +                          * pdf(nbin,nx,isub) ! multiply PDFs
                           sum=sum+contr
                           v = real(contr)
                           call hfill( ihist, real(k-0.5d0) ,0.0, v)
                           call hfill(-ihist, real(l-0.5d0) ,0.0, v)
                           call hfill( ihi2 , real(k-0.5d0) ,0.0, v)
                           call hfill(-ihi2 , real(l-0.5d0) ,0.0, v)
                           call hfill( ihi3 , real(k-0.5d0) ,0.0, v)
                           call hfill(-ihi3 , real(l-0.5d0) ,0.0, v)
                           call hfill( ihi4 , real(k-0.5d0) ,0.0, v)
                           call hfill(-ihi4 , real(l-0.5d0) ,0.0, v)

                        enddo   ! iord
                     enddo      ! subproc
                  enddo         ! xmin
               enddo            ! xmax
               if (iscale.eq.1) write(*,*) " rapbin ",i,"  ptbin ",j," -> ",sum
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
      INCLUDE 'fnx9999.inc'
      real pt(nptmax)

c - HBOOK common
      INTEGER NWPAWC
      PARAMETER (NWPAWC=2500000)
      REAL HMEMOR(NWPAWC)
      COMMON /PAWC/ HMEMOR

      write(*,*) "       FNHBOOK  n=",n

c - open & book
      if (n.eq.1) then
         CALL HLIMIT(NWPAWC)
         CALL HROPEN(11,'fastNLO',filename,'N',1024,ISTAT2)
         if (ISTAT2.NE.0) then
            WRITE(*,*) ' FNHBOOK: could not open histofile ',istat2
         endif

         do iord=0,2            ! order: tot, LO, NLO-corr
            do iscale=1,NSCALEMAX ! scale variations
               do isub=0,7      ! subprocess: 0 tot + 7 subproc
               
                  do irap=1, nrapidity
                     do j=1,(npt(irap)+1)

                        pt(j) = real(PTBIN(irap,j))
c                        write(*,*) "  irap,j=",irap,j,"  pt=",pt(j)
                        ihist = iord*1000000+iscale*100000+isub*10000+irap*100

c --------------------- x-distributions in pT-bins - only for default scale
                        if (ix.eq.2 .and. iscale.eq.1 .and.j.ne.(npt(irap)+1)) then
                           call hbook1(ihist+j,'x-Bin',ntot,0.0,real(ntot),0)
                           call hbook1(-(ihist+j),'x-Bin',ntot,0.0,real(ntot),0)
                        endif

                     enddo
                     call hbookb(ihist,'E?T! (GeV)' , NPT(irap) ,PT ,0)
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


