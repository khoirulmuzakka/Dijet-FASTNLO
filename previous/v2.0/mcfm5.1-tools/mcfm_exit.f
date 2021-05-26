      subroutine mcfm_exit(itmx,xinteg,xinteg_err)
************************************************************************
*                                                                      *
*  This routine should perform the final processing and print-outs     *
*                                                                      *
************************************************************************
      implicit none
      include 'efficiency.f'
      include 'PDFerrors.f'
      integer j,k,itno,itmx
      double precision xinteg,xinteg_err,minPDFxsec,maxPDFxsec
      double precision PDFerror,PDFperror,PDFnerror
      double precision lord_bypart(-1:1,-1:1),lordnorm
      double precision ggpart,gqpart,qgpart,qqpart,qqbpart
      character*4 part
      logical creatent,dswhisto
      common/outputflags/creatent,dswhisto
      common/part/part
      common/bypart/lord_bypart
      common/finalpart/ggpart,gqpart,qgpart,qqpart,qqbpart
      
c--- Print-out the value of the integral and its error
      write(6,*) 
      write(6,53)'Value of final ',part,' integral is',
     . xinteg,' +/-',xinteg_err, ' fb'
     
   53 format(a15,a4,a12,f13.3,a4,f10.3,a3)

c--- Print-out a summary of the effects of jets and cuts
      write(6,*) 
      write(6,*) 'Total number of shots       : ',ntotshot
      write(6,*) 'Total no. failing cuts      : ',ntotzero
      write(6,*) 'Number failing jet cuts     : ',njetzero
      write(6,*) 'Number failing process cuts : ',ncutzero
      write(6,*) 
      call flush(6)

c--- Calculate the actual number of shots that were passed
c--- through the jet and cut routines
      ntotshot=ntotshot-(ntotzero-njetzero-ncutzero)
      write(6,54) 'Jet efficiency : ',
     .  100d0-100d0*dfloat(njetzero)/dfloat(ntotshot)
      write(6,54) 'Cut efficiency : ',
     .  100d0-100d0*dfloat(ncutzero)/dfloat((ntotshot-njetzero))
      write(6,54) 'Total efficiency : ',
     .  100d0-100d0*dfloat((njetzero+ncutzero))/dfloat(ntotshot)
      write(6,*) 
      
      lordnorm=0d0
      do j=-1,1
      do k=-1,1
        lordnorm=lordnorm+lord_bypart(j,k)
      enddo
      enddo
      ggpart=lord_bypart( 0, 0)/lordnorm
      gqpart=(lord_bypart( 0,+1)+lord_bypart( 0,-1))/lordnorm
      qgpart=(lord_bypart(+1, 0)+lord_bypart(-1, 0))/lordnorm
      qqpart=(lord_bypart(+1,+1)+lord_bypart(-1,-1))/lordnorm
      qqbpart=(lord_bypart(+1,-1)+lord_bypart(-1,+1))/lordnorm
      write(6,*) 'Contribution from parton sub-processes:'
      write(6,*) '---------------------------------------'      
      write(6,55) '   GG    ',ggpart*xinteg,ggpart*100d0
      write(6,55) 'GQ + GQB ',gqpart*xinteg,gqpart*100d0
      write(6,55) 'QG + QBG ',qgpart*xinteg,qgpart*100d0
      write(6,55) 'QQ + QBQB',qqpart*xinteg,qqpart*100d0
      write(6,55) '   QQB   ',qqbpart*xinteg,qqbpart*100d0
      write(6,*) '---------------------------------------'
      call flush(6)

   54 format(a20,f6.2,'%')
   55 format(4x,a9,' |',f15.5,f8.2,'%')

c--- if we've calculated PDF errors, present results      
      if (PDFerrors) then
        write(6,*)
        write(6,58) '************ PDF error analysis ************'
        write(6,58) '*                                          *'
        minPDFxsec=PDFxsec(0)
        maxPDFxsec=PDFxsec(0)
        PDFerror=0d0
        do j=0,maxPDFsets
          write(6,56) j,PDFxsec(j)
          if     (PDFxsec(j) .lt. minPDFxsec) then
            minPDFxsec=PDFxsec(j)
          elseif (PDFxsec(j) .gt. maxPDFxsec) then
            maxPDFxsec=PDFxsec(j)
         endif
         if ( (j .gt. 0) .and. (j/2 .eq. (j-1)/2) ) then
           PDFerror=PDFerror+(PDFxsec(j)-PDFxsec(j+1))**2
         endif
         if ( (j .gt. 0) .and. (PDFxsec(j) .gt. PDFxsec(0)) ) then
           PDFperror=PDFperror+(PDFxsec(j)-PDFxsec(0))**2
         endif
         if ( (j .gt. 0) .and. (PDFxsec(j) .lt. PDFxsec(0)) ) then
           PDFnerror=PDFnerror+(PDFxsec(j)-PDFxsec(0))**2
         endif
        enddo
        PDFerror=0.5d0*dsqrt(PDFerror)
        PDFperror=dsqrt(PDFperror)
        PDFnerror=dsqrt(PDFnerror)
        write(6,58) '*                                          *'
        write(6,58) '* --------------- SUMMARY ---------------- *'
        write(6,58) '*                                          *'
        write(6,57) 'Minimum value',minPDFxsec
        write(6,57) 'Central value',PDFxsec(0)
        write(6,57) 'Maximum value',maxPDFxsec
        write(6,58) '*                                          *'
        write(6,57) 'Err estimate +/-',PDFerror
        write(6,57) '   +ve direction',PDFperror
        write(6,57) '   -ve direction',PDFnerror
        write(6,58) '********************************************'
      endif
      
   56 format('* PDF error set ',i3,'  --->',f13.3,' fb  *')
   57 format('*   ',a16,f14.3,' fb      *')
   58 format(a44)
 
c--- Finalize the histograms, if we're not filling ntuples instead
      if (creatent .eqv. .false.) then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
          call histofin(xinteg,xinteg_err,0,itmx)
        else
c--- DSW histograms - store the information
          call dswhbook(200,'Sigma   ',1.0d0,0.0d0,10.0d0)
          call dswhfill(200,0.5d0,xinteg)
          call dswhfill(200,1.5d0,xinteg_err)
c--- DSW histograms - output and close file
          call dswhrout
          call dswclose
        endif
      else
c--- ADDED - to produce normal histograms as well
        call histofin(xinteg,xinteg_err,0,itmx)
        call dswhrout
        call dswclose
      endif

c === fastNLO
      call FnTerm
c ==============================================
      return
      
      end
