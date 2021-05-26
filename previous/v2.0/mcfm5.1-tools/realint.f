      double precision function realint(vector,wgt)
      implicit none
      include 'constants.f'
      include 'debug.f'
      include 'realonly.f'
      include 'virtonly.f'
      include 'noglue.f'
      include 'nflav.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'new_pspace.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'realwt.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'process.f'
      include 'PDFerrors.f'
      include 'wts_bypart.f'
c === fastNLO
      include 'fastnlo.f'
c ============================================
cz
      integer nproc
      common/nproc/nproc
cz //
      integer ih1,ih2,j,k,nd,nmax,nmin,nvec
      double precision vector(mxdim),W,val,val2,xint,dot,tmp
      double precision sqrts,fx1(-nf:nf),fx2(-nf:nf)
      double precision p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      double precision pswt,rscalestart,fscalestart
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd)
      double precision flux,BrnRat,xreal,xreal2
      double precision xx1,xx2,q(mxpart,4),msqa(-nf:nf,-nf:nf)
      double precision m3,m4,m5,R,Rbbmin
      double precision xmsq_bypart(0:maxd,-1:1,-1:1),xmsqjk,
     . lord_bypart(-1:1,-1:1)
      integer n2,n3,sgnj,sgnk
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      common/xreal/xreal,xreal2
      common/Rbbmin/Rbbmin
      logical bin,first,failed
      logical incldip(0:maxd),includedipole,includereal
      logical creatent,dswhisto
      external qqb_w2jet_g,qqb_w2jet_gs,qqb_z2jet_g,qqb_z2jet_gs,
     . qqb_w2jet,qqb_w1jet_gs,qqb_z2jet,qqb_z1jet_gs,qqb_Hg_g,qqb_Hg_gs,
     . qqb_hww_g,qqb_hww_gs,qqb_zbb_g,qqb_zbb_gs,
     . qqb_wbb_g,qqb_wbb_gs,
     . qqb_w_g,qqb_w_gs,qqb_z1jet,qqb_z_gs,qqb_ww_g,qqb_ww_gs,
     . qqb_wz_g,qqb_wz_gs,qqb_zz_g,qqb_zz_gs,qqb_wgam_g,qqb_wgam_gs,
     . qqb_QQb_g,qqb_QQb_gs,
     . VV_Hqq_g,VV_Hqq_gs,gg_Hg,gg_H_gs,gg_Hgg,gg_Hg_gs,
     . gQ_zQ_g,gQ_zQ_gs,qqb_tbb_g,qqb_tbb_gs,
     . qqb_w_tndk_g,qqb_w_tndk_gs,
     . qqb_w_twdk_g,qqb_w_twdk_gs,qqb_w_twdk_gdk,qqb_w_twdk_gsdk,
     . qqb_zbjet_g,qqb_zbjet_gs
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/Pext/p1ext,p2ext
      common/nmax/nmax
      common/BrnRat/BrnRat
      common/nmin/nmin
      common/bypart/lord_bypart
      common/incldip/incldip
      common/outputflags/creatent,dswhisto
cz Add b fraction
      double precision bwgt
      common/btagging/ bwgt
      double precision msqtmp(0:maxd),bwgttmp(0:maxd)
      double precision realeventp(mxpart,4)
      common/realeventp/realeventp
      
      data bwgt / 0d0 /  ! in common block
cz // Add b fraction   Note: only msqtmp(0), bwgttmp(0) are used in nplotter.f
      data p/48*0d0/
      data first/.true./
      save first,rscalestart,fscalestart
      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif
      ntotshot=ntotshot+1
      pswt=0d0
      realint=0d0      

      W=sqrts**2
      
      if (first) then
         write(6,*)
         write(6,*) 'nmin=',nmin,',nmax=',nmax
         write(6,*)
         first=.false.
      endif
      
c--- processes that use "gen3"     
      if     ( (case .eq. 'W_only')
     .    .or. (case .eq. 'Z_only')
     .    .or. (case .eq. 'ggfus0') ) then
        npart=3
        if (new_pspace) then
          call gen3a(vector,p,pswt,*999)
        else
          call gen3(vector,p,pswt,*999)
        endif

c--- processes that use "gen3m"     
      elseif ( (case .eq. 'tt_tot')
     .    .or. (case .eq. 'bb_tot')
     .    .or. (case .eq. 'cc_tot') ) then
        m3=mass2
        m4=mass2
        m5=0d0
        npart=3
        call gen3m(vector,p,m3,m4,m5,pswt,*999)
          
c--- processes that use "gen4"     
      elseif ( (case .eq. 'W_cjet') 
     .   .or.  (case .eq. 'Wgamma')
     .   .or.  (case .eq. 'Zgamma')
     .   .or.  (case .eq. 'W_tndk') ) then
        npart=4
        call gen4(vector,p,pswt,*999)
                  
c--- processes that use "gen6"     
      elseif ( (case .eq. 'W_twdk') .or. (case .eq. 'Wtdkay') ) then
        npart=6
        call gen6(vector,p,pswt,*999)
                  
c--- processes that use "gen_njets" with an argument of "2"     
      elseif ( (case .eq. 'W_1jet')
     .    .or. (case .eq. 'Wcjet0')
     .    .or. (case .eq. 'Z_1jet')
     .    .or. (case .eq. 'H_1jet')
     .    .or. (case .eq. 'ggfus1')
     .    .or. (case .eq. 'gQ__ZQ') ) then
        npart=4
        if (new_pspace) then
          call gen4a(vector,p,pswt,*999)      
        else
          call gen_njets(vector,2,p,pswt,*999)
        endif 
        
c--- processes that use "gen_njets" with an argument of "3"
      elseif ( (case .eq. 'Wbbbar')
     .    .or. (case .eq. 'W_2jet')
     .    .or. (case .eq. 'Z_2jet')
     .    .or. (case .eq. 'Zbbbar')
     .    .or. (case .eq. 'W_bjet')
     .    .or. (case .eq. 'Z_bjet')
     .    .or. (case .eq. 'qq_Hqq')
     .    .or. (case .eq. 'ggfus2') ) then
        npart=5
        call gen_njets(vector,3,p,pswt,*999) 
        
c--- processes that use "gen_stop" with an argument of "1"
      elseif ( (case .eq. 'ttdkay')
     .    .or. (case .eq. 'tdecay') ) then
        npart=5
        call gen_stop(vector,1,p,pswt,*999)

c--- processes that use "gen_stop" with an argument of "2"
      elseif ( (case .eq. 'bq_tpq')
     .    .or. (case .eq. 't_bbar') ) then
        npart=5
        call gen_stop(vector,2,p,pswt,*999)
        
c--- DEFAULT: processes that use "gen5"
      else
        npart=5
        if (new_pspace) then
          call gen5a(vector,p,pswt,*999)
        else
          call gen5(vector,p,pswt,*999)    
        endif
      endif
            
      nvec=npart+2
      call dotem(nvec,p,s)
      
c---impose cuts on final state
      call masscuts(s,*999)
c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
     
c--- extra cut to divide WQj/ZQj regions
      if ( (nproc .eq. 312) .or. (nproc .eq. 317)
     . .or.(nproc .eq. 322) .or. (nproc .eq. 327)
     . .or.(nproc .eq. 342) .or. (nproc .eq. 352)) then
        if (R(p,5,6) .lt. Rbbmin) goto 999
      endif

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead set to zero
      includereal=includedipole(0,p)
      incldip(0)=includereal 
 
      if (includereal .eqv. .false.) then
        do j=-nf,nf
        do k=-nf,nf
          msq(j,k)=0d0
        enddo
        enddo
      endif
      
      do j=1,mxpart
      do k=1,4
        realeventp(j,k)=p(j,k)
      enddo
      enddo
      
      if (dynamicscale) call scaleset(rscalestart,fscalestart,p)
      
c---- generate collinear points that satisfy the jet cuts (for checking)
c      call singgen(p,s,*998)
            
c----calculate the x's for the incoming partons from generated momenta

      xx1=two*(p(1,4)*p2ext(4)-p(1,3)*p2ext(3))/W
      xx2=two*(p(2,4)*p1ext(4)-p(2,3)*p1ext(3))/W

      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx1,xx2
      
      if ((xx1 .gt. 1d0) .or. (xx2 .gt. 1d0)) then
         realint=0d0
         return
      endif

c--- Calculate the required matrix elements      
      if     (case .eq. 'W_only') then
c        call singcheck(qqb_w_g,qqb_w_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_w_g(p,msq)      
        call qqb_w_gs(p,msqc)     
      elseif (case .eq. 'W_1jet') then
c        call singcheck(qqb_w2jet,qqb_w1jet_gs,p)   ! Checked 11/16/01
        if (includereal) call qqb_w2jet(p,msq)      
        call qqb_w1jet_gs(p,msqc)  
      elseif (case .eq. 'Wgamma') then
c        call singcheck(qqb_wgam_g,qqb_wgam_gs,p)   ! Checked 08/27/02
        if (includereal) call qqb_wgam_g(p,msq)      
        call qqb_wgam_gs(p,msqc)  
      elseif (case .eq. 'Zgamma') then
c        call singcheck(qqb_zgam_g,qqb_zgam_gs,p)
        if (includereal) call qqb_zgam_g(p,msq)      
        call qqb_zgam_gs(p,msqc)  
      elseif (case .eq. 'Wbbbar') then
c        call singcheck(qqb_wbb_g,qqb_wbb_gs,p)     ! Checked 11/30/01
        if (includereal) call qqb_wbb_g(p,msq)      
        call qqb_wbb_gs(p,msqc)      
      elseif (case .eq. 'W_2jet') then      
c        call singcheck(qqb_w2jet_g,qqb_w2jet_gs,p) ! Re-checked 4/26/05
        if (includereal)  call qqb_w2jet_g(p,msq)  
        call qqb_w2jet_gs(p,msqc)
      elseif (case .eq. 'Z_only') then
c        call singcheck(qqb_z1jet,qqb_z_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_z1jet(p,msq)      
        call qqb_z_gs(p,msqc)     
      elseif (case .eq. 'Z_1jet') then
c        call singcheck(qqb_z2jet,qqb_z1jet_gs,p)   ! Checked 11/16/01
        if (includereal) call qqb_z2jet(p,msq)      
        call qqb_z1jet_gs(p,msqc)  
      elseif (case .eq. 'Z_2jet') then
c        call singcheck(qqb_z2jet_g,qqb_z2jet_gs,p) ! Checked 11/16/01
        if (includereal) call qqb_z2jet_g(p,msq)  
        call qqb_z2jet_gs(p,msqc) 
      elseif (case .eq. 'Zbbbar') then
c        call singcheck(qqb_zbb_g,qqb_zbb_gs,p)     ! Checked 11/30/01
        if (includereal) call qqb_zbb_g(p,msq)
        call qqb_zbb_gs(p,msqc) 
      elseif (case .eq. 'WWqqbr') then
c        call singcheck(qqb_ww_g,qqb_ww_gs,p)       ! Checked 11/30/01
        if (includereal) call qqb_ww_g(p,msq)      
        call qqb_ww_gs(p,msqc)      
      elseif (case .eq. 'WZbbar') then
c        call singcheck(qqb_wz_g,qqb_wz_gs,p)       ! Checked 12/05/01
        if (includereal) call qqb_wz_g(p,msq)      
        call qqb_wz_gs(p,msqc)      
      elseif (case .eq. 'ZZlept') then
c        call singcheck(qqb_zz_g,qqb_zz_gs,p)       ! Checked 12/05/01
        if (includereal) call qqb_zz_g(p,msq)      
        call qqb_zz_gs(p,msqc)      
      elseif (case .eq. 'WHbbar') then
c        call singcheck(qqb_wh_g,qqb_wh_gs,p)
        if (includereal) call qqb_wh_g(p,msq)      
        call qqb_wh_gs(p,msqc)     
      elseif (case .eq. 'ZHbbar') then
c        call singcheck(qqb_zh_g,qqb_zh_gs,p)
        if (includereal) call qqb_zh_g(p,msq)      
        call qqb_zh_gs(p,msqc)     
       elseif (case .eq. 'ggfus0') then
c         call singcheck(gg_hg,gg_h_gs,p)       ! Checked 28/02/03
         if (includereal) call gg_hg(p,msq)
         call gg_h_gs(p,msqc)
      elseif (case .eq. 'HWW_4l') then
c        call singcheck(qqb_hww_g,qqb_hww_gs,p)
        if (includereal) call qqb_hww_g(p,msq)      
        call qqb_hww_gs(p,msqc)     
      elseif (case .eq. 'HZZ_4l') then
c        call singcheck(qqb_hzz_g,qqb_hzz_gs,p)
        if (includereal) call qqb_hzz_g(p,msq)      
        call qqb_hzz_gs(p,msqc)  
      elseif (case .eq. 'H_1jet') then
c        call singcheck(qqb_Hg_g,qqb_Hg_gs,p)       ! Checked 19/02/02
        if (includereal) call qqb_Hg_g(p,msq)  
        call qqb_Hg_gs(p,msqc) 
      elseif ((case .eq. 'tt_tot') .or. (case .eq. 'cc_tot')
     .   .or. (case .eq. 'bb_tot')) then
c        call singcheck(qqb_QQb_g,qqb_QQb_gs,p)
        if (includereal) call qqb_QQb_g(p,msq)
        call qqb_QQb_gs(p,msqc)
      elseif ((case .eq. 'bq_tpq') .or. (case .eq. 't_bbar')) then
c        call singcheck(qqb_tbb_g,qqb_tbb_gs,p)
        if (includereal) call qqb_tbb_g(p,msq)
        call qqb_tbb_gs(p,msqc)
      elseif (case .eq. 'ttdkay') then
c        call singcheck(qqb_tbb_g_new,qqb_tbb_gs,p)
        if (includereal) call bq_tpq_gdk(p,msq)
        call bq_tpq_gsdk(p,msqc)
      elseif (case .eq. 'tdecay') then
c        call singcheck(qqb_tbb_g_new,qqb_tbb_gs,p)
        if (includereal)  call qqb_tbb_gdk(p,msq)
        call qqb_tbb_gsdk(p,msqc)       
       elseif (case .eq. 'W_tndk') then
c        call singcheck(qqb_w_tndk_g,qqb_w_tndk_gs,p)	! Checked 12/3/04
        if (includereal) call qqb_w_tndk_g(p,msq)
        call qqb_w_tndk_gs(p,msqc)
      elseif (case .eq. 'W_twdk') then
c        call singcheck(qqb_w_twdk_g,qqb_w_twdk_gs,p)		! Checked 2/4/05
        if (includereal) call qqb_w_twdk_g(p,msq)
        call qqb_w_twdk_gs(p,msqc)
      elseif (case .eq. 'Wtdkay') then
c        call singcheck(qqb_w_twdk_gdk,qqb_w_twdk_gsdk,p)	! Checked 2/2/05
        if (includereal) call qqb_w_twdk_gdk(p,msq)
        call qqb_w_twdk_gsdk(p,msqc)
      elseif (case .eq. 'ggfus1') then
c        call singcheck(gg_hgg,gg_hg_gs,p)
        if (includereal) call gg_hgg(p,msq)
        call gg_hg_gs(p,msqc)
      elseif (case .eq. 'qq_Hqq') then
c        call singcheck(VV_Hqq_g,VV_Hqq_gs,p)
        if (includereal) call VV_Hqq_g(p,msq)
        call VV_Hqq_gs(p,msqc)
      elseif (case .eq. 'gQ__ZQ') then
c        call singcheck(gQ_zQ_g,gQ_zQ_gs,p)
        if (includereal) call gQ_zQ_g(p,msq)
        call gQ_zQ_gs(p,msqc)
      elseif (case .eq. 'Z_bjet') then
c        call singcheck(qqb_zbjet_g,qqb_zbjet_gs,p)	! Checked 07/18/05
        if (includereal) call qqb_zbjet_g(p,msq)
        call qqb_zbjet_gs(p,msqc)
      endif
      
      do nd=0,ndmax
      xmsq(nd)=0d0
cz
      msqtmp(nd)=0d0
      bwgttmp(nd)=0d0
cz //
      do j=-1,1
      do k=-1,1
      xmsq_bypart(nd,j,k)=0d0
      enddo
      enddo
      enddo
      
      currentPDF=0
            
      flux=fbGeV2/(two*xx1*xx2*W)

c--- initialize a PDF set here, if calculating errors
  777 continue    
      do nd=0,ndmax
      xmsq(nd)=0d0
cz
      msqtmp(nd)=0d0
      bwgttmp(nd)=0d0
cz //
      enddo
      if (PDFerrors) then
        call InitPDF(currentPDF)
      endif
            
c--- calculate PDF's  
      call fdist(ih1,xx1,facscale,fx1)
      call fdist(ih2,xx2,facscale,fx2)
            
      do j=-nflav,nflav
      do k=-nflav,nflav
c      do j=2,2
c      do k=0,0

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
      endif      
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif

      if (realonly) then 
        xmsq(0)=xmsq(0)+fx1(j)*fx2(k)*msq(j,k)
        do nd=1,ndmax
        xmsq(nd)=0d0
        enddo
      elseif (virtonly) then
         xmsq(0)=0d0
         do nd=1,ndmax
           xmsq(nd)=xmsq(nd)+fx1(j)*fx2(k)*(-msqc(nd,j,k))
         enddo
      else

         if     (j .gt. 0) then
           sgnj=+1
         elseif (j .lt. 0) then
           sgnj=-1
         else
           sgnj=0
         endif
         if     (k .gt. 0) then
           sgnk=+1
         elseif (k .lt. 0) then
           sgnk=-1
         else
           sgnk=0
         endif

         xmsqjk=fx1(j)*fx2(k)*msq(j,k)
         xmsq(0)=xmsq(0)+xmsqjk
cz
cz Extract fraction with b in final state, store in common
cz
c     nproc=161 then t-channel: t
c     isub = 1, nwz = 1
         if (nproc.eq.161) then ! t
            msqtmp(0)=msqtmp(0)+xmsqjk
            if ((j.eq.0).and.((k.lt.0).or.((k.gt.0).and.
     &           (k.ne.5)))) then
               bwgttmp(0)=bwgttmp(0)+xmsqjk
            endif
            if ((k.eq.0).and.((j.lt.0).or.((j.gt.0).and.
     &           (j.ne.5)))) then
               bwgttmp(0)=bwgttmp(0)+xmsqjk
            endif
         endif
c     nproc=166 then t-channel: t~
         if (nproc.eq.166) then ! t
            msqtmp(0)=msqtmp(0)+xmsqjk
            if ((j.eq.0).and.((k.gt.0).or.((k.lt.0).and.
     &           (k.ne.-5)))) then
               bwgttmp(0)=bwgttmp(0)+xmsqjk
            endif
            if ((k.eq.0).and.((j.gt.0).or.((j.lt.0).and.
     &           (j.ne.-5)))) then
               bwgttmp(0)=bwgttmp(0)+xmsqjk
            endif
         endif
cz // end fill index 0

         if (currentPDF .eq. 0) then
           xmsq_bypart(0,sgnj,sgnk)=xmsq_bypart(0,sgnj,sgnk)+xmsqjk
         endif

c === fastNLO ---------
      fnpdf(j,k) = fx1(j)*fx2(k)
      fnx(1) = xx1
      fnx(2) = xx2
      fncftmp(0,j,k) = msq(j,k)
      Do nd=1,ndmax
         fncftmp(nd,j,k) = -msqc(nd,j,k)
      Enddo
c ==============================

         do nd=1,ndmax
           xmsqjk=fx1(j)*fx2(k)*(-msqc(nd,j,k))
           xmsq(nd)=xmsq(nd)+xmsqjk
           if (currentPDF .eq. 0) then
             xmsq_bypart(nd,sgnj,sgnk)=xmsq_bypart(nd,sgnj,sgnk)+xmsqjk
           endif
         enddo
         
      endif
 20   continue


      enddo
      enddo

c--- loop over all PDF error sets, if necessary
      if (PDFerrors) then
        do nd=0,ndmax 
          PDFxsec_nd(currentPDF,nd)=xmsq(nd)
        enddo
        currentPDF=currentPDF+1
        if (currentPDF .le. maxPDFsets) goto 777
c--- reset xmsq to the central PDF values
        do nd=0,ndmax 
          xmsq(nd)=PDFxsec_nd(0,nd)
        enddo
      endif    

      realint=0d0
      xint=0d0
cz
      bwgt=0d0
cz //

c---trial with weight of real alone
c---first set up all dipole contributions
c---this is the value of integral including subtractions
      do nd=0,ndmax
        xmsq(nd)=xmsq(nd)*flux*pswt/BrnRat
        if (creatent) then
          wt_gg=xmsq_bypart(nd,0,0)*wgt*flux*pswt/BrnRat/dfloat(itmx)
          wt_gq=(xmsq_bypart(nd,+1,0)+xmsq_bypart(nd,-1,0)
     .          +xmsq_bypart(nd,0,+1)+xmsq_bypart(nd,0,-1)
     .          )*wgt*flux*pswt/BrnRat/dfloat(itmx)
          wt_qq=(xmsq_bypart(nd,+1,+1)+xmsq_bypart(nd,-1,-1)
     .          )*wgt*flux*pswt/BrnRat/dfloat(itmx)
          wt_qqb=(xmsq_bypart(nd,+1,-1)+xmsq_bypart(nd,-1,+1)
     .          )*wgt*flux*pswt/BrnRat/dfloat(itmx)
        endif
        failed=.false.
        
c--- if this dipole has no contribution, go to end of loop
c        if (xmsq(nd) .eq. 0d0) goto 997         
         
        if (nd .eq. 0) then
c---if there's no real contribution, record the event as failing to pass cuts
          if (xmsq(nd) .eq. 0d0) then
             failed=.true.
             goto 996
          endif
        else
c--- if this dipole has no contribution, go to end of loop
          if (xmsq(nd) .eq. 0d0) goto 997         
c---check whether each counter-event passes the cuts
          do j=1,mxpart
          do k=1,4
          q(j,k)=ptilde(nd,j,k)
          enddo
          enddo
c          write(6,*)
c          write(6,*) 'Dipole nd=',nd
          if (incldip(nd)) incldip(nd)=includedipole(nd,q)
          if (incldip(nd) .eqv. .false.) failed=.true.
        endif

 996    if (failed) then
          if (nd .eq. 0) then
            ncutzero=ncutzero+1
            ntotzero=ntotzero+1
          endif
          call dotem(nvec,p,s)
          xmsq(nd)=0d0
          goto 997         
        endif

c---if it does, add to total
        xint=xint+xmsq(nd)
        do j=-1,1
        do k=-1,1
          lord_bypart(j,k)=lord_bypart(j,k)+
     .         wgt*flux*pswt*xmsq_bypart(nd,j,k)/BrnRat
        enddo
        enddo

        val=xmsq(nd)*wgt
        val2=xmsq(nd)**2*wgt
cz Fill bwgt if needed
        if(dabs(msqtmp(nd)).gt.0d0) bwgt=bwgttmp(nd)/msqtmp(nd)
cz //
        
c--- update PDF errors
        if (PDFerrors) then
          do currentPDF=0,maxPDFsets        
          PDFwgt(currentPDF)=
     .       flux*pswt*PDFxsec_nd(currentPDF,nd)/BrnRat*wgt/itmx
          PDFxsec(currentPDF)=PDFxsec(currentPDF)
     .       +PDFwgt(currentPDF)
          enddo           
        endif
                
c--- update the maximum weight so far, if necessary
        if (dabs(val) .gt. wtmax) then
          wtmax=dabs(val)
        endif

c === fastNLO
        fnsigma = flux*pswt/BrnRat*wgt/dfloat(itmx)
        Do j=-nflav,nflav
           Do k=-nflav,nflav
              fncoeff(j,k) = fncftmp(nd,j,k)
           Enddo
        Enddo
        ifnctrb = 3
c ===================================================

c---if we're binning, add to histo too
        if (bin) then
          call getptildejet(nd,pjet)
          call dotem(nvec,pjet,s)
          val=val/dfloat(itmx)
          val2=val2/dfloat(itmx)**2
          if (nd .eq. 0) then
            call nplotter(pjet,val,val2,0)
          else
            call nplotter(pjet,val,val2,1)
          endif
        endif
c---otherwise, skip contribution
 997    continue
      enddo

c 998  continue

      if (realwt) then
      realint=xmsq(0)
      else
      realint=xint
      endif
      xreal=xreal+xint*wgt/dfloat(itmx)
      xreal2=xreal2+(xint*wgt)**2/dfloat(itmx)
      

      return

 999  realint=0d0
      ntotzero=ntotzero+1
 
      return
      
      end
















