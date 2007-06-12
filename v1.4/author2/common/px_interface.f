      SUBROUTINE PXINTERFACE(np,nj, vpp,vjj)
************************************************************************
* M Wobisch 03/18/04   interface to px_cone
*
* right now: max of 4 particles
*
* input   
*    integer  np       No Particles
*    double   vp(.,.)  particle array  - first index: particle
* output
*    integer  nj       No jets
*    double   vj(.,.)  jet array
*
**************************************************************************
      implicit none
      integer  np, nj
      double precision vpp(16),vjj(16), vp(4,4), vj(4,4)
      integer i,ii,j, ierr
      double precision PJET(5,4), ipass(4), ijmul(4), pt1,pt2

      nj = 0
      do i=1,np
         do j=1,4
            vp(j,i) = vpp((i-1)*4+j)
         enddo
      enddo

      call PXCONE(3, np, 4, vp, 0.7d0, 8d0, 0.5d0,
     +                   4, nj, PJET, IPASS,IJMUL,IERR)

      if (ierr.ne.0) then
         write(*,*) "PXinterface:  IERR = ",ierr,"  Njets=",nj
         do i=1,nj
            write(*,*)'  jets:',vj(1,i),vj(2,i),vj(3,i),vj(4,i)
         enddo 
         write(*,*)'    -----------------------------'
         nj = 0
      endif

c      write(*,*) "..",nj
c      do i=1,np
c         write(*,*) "..",vp
c      enddo

      do i=1,nj
         do j=1,4
            vj(j,i) = PJET(j,i)
            vjj((i-1)*4+j) = vj(j,i)
         enddo
      enddo

c      if (nj.lt.np) then
c         do i=1,np
cc            write(*,*)'in: ',vp(1,i),vp(2,i),vp(3,i),vp(4,i)
c         enddo
c         do i=1,nj
cc            write(*,*)'out:',vj(1,i),vj(2,i),vj(3,i),vj(4,i)
c         enddo
cc         write(*,*)'   -------------------------------'
c      endif


c - check: are outgoing jets ordered in energy or pT??
c      do i=1,(nj-1)
c         pt1 = sqrt(vj(1,i)**2   + vj(2,i)**2)
c         pt2 = sqrt(vj(1,i+1)**2 + vj(2,i+1)**2)
c         if ((pt1+0.00000001d0).lt.pt2) then
c            write(*,*) " MW  no pT ordering ",i,nj," ",pt1,"  ",pt2
c            do ii=1,np
c               write(*,*)'in: ',vp(1,ii),vp(2,ii),vp(3,ii),vp(4,ii)
c            enddo
c            do ii=1,nj
c               write(*,*)'out:',vj(1,ii),vj(2,ii),vj(3,ii),vj(4,ii)
c            enddo
c            write(*,*)'   -------------------------------'
c         endif
c      enddo
c
c   --> they  *are* ordered in pT !!!!!!!!!!!!!!! (checked in E-scheme)
c
      RETURN
      END
