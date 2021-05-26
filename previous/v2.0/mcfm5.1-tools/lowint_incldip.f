      double precision function lowint(r,wgt)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'limits.f'
      include 'npart.f'
      include 'debug.f'
      include 'new_pspace.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'noglue.f'
      include 'process.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'phasemin.f'
      include 'PDFerrors.f'
      include 'wts_bypart.f'
c --- DSW. To store flavour information :
      include 'nflav.f'
c      include 'b0.f'
c --- DSW.
c === fastNLO
      include 'fastnlo.f'
c ============================================

      integer pflav,pbarflav
c --- To use VEGAS random number sequence :
      double precision ran2
      integer ih1,ih2,j,k,nvec,sgnj,sgnk
      double precision r(mxdim),W,sqrts,xmsq,val,val2,
     . fx1(-nf:nf),fx2(-nf:nf),p(mxpart,4),pjet(mxpart,4),
     . pswt,rscalestart,fscalestart
      double precision wgt,msq(-nf:nf,-nf:nf),m3,m4,m5,xmsqjk
      double precision xx(2),flux,vol,vol_mass,vol3_mass,vol_wt,BrnRat
      double precision xmsq_bypart(-1:1,-1:1),lord_bypart(-1:1,-1:1)
      logical bin,first,includedipole
      logical creatent,dswhisto
c      double precision gx1(-nf:nf),gx2(-nf:nf)
c      integer idum
c      COMMON/ranno/idum
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/x1x2/xx
      common/BrnRat/BrnRat
      common/bypart/lord_bypart
      common/outputflags/creatent,dswhisto
      data p/48*0d0/
      data first/.true./
      save first,rscalestart,fscalestart
      
      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif

      ntotshot=ntotshot+1
      lowint=0d0

      W=sqrts**2

c--- processes that use "gen2"     
      if     ( (case .eq. 'W_only')
     .    .or. (case .eq. 'Z_only')
     .    .or. (case .eq. 'ggfus0')
     .    .or. (case .eq. 'vlchk2') ) then
        if (case .eq. 'vlchk2') then
          wsqmin=0d0
          wsqmax=sqrts**2
        endif
        npart=2
        if (new_pspace) then
          call gen2a(r,p,pswt,*999)
        else
          call gen2(r,p,pswt,*999)
        endif

c--- processes that use "gen2m"     
      elseif ( (case .eq. 'tt_tot')
     .    .or. (case .eq. 'bb_tot')
     .    .or. (case .eq. 'cc_tot') ) then
        npart=2
        call gen2m(r,p,pswt,*999)
          
c--- processes that use "gen3"     
      elseif ( (case .eq. 'W_cjet') 
     .   .or.  (case .eq. 'W_tndk')
     .   .or.  (case .eq. 'vlchwn') ) then
        npart=3
        call gen3(r,p,pswt,*999)

c--- processes that use "gen3jet"     
      elseif ( (case .eq. 'Wgamma') 
     .   .or.  (case .eq. 'Zgamma')
     .   .or.  (case .eq. 'vlchk3') ) then
        if (case .eq. 'vlchk3') then
          wsqmin=0d0
          wsqmax=sqrts**2
        endif
        npart=3
        call gen3jet(r,p,pswt,*999)

c--- processes that use "gen3m"     
      elseif ( (case .eq. 'tottth') ) then
        m3=mt
        m4=mt
        m5=hmass
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)
          
c--- processes that use "gen3m"     
      elseif (case .eq. 'tt_glu') then
        m3=mt
        m4=mt
        m5=0d0
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)
c--- processes that use "gen3m_rap"     
      elseif ( (case .eq. 'vlchm3') ) then
        taumin=(2d0*mt/sqrts)**2
        m3=mt
        m4=mt
        npart=3
        call gen3m_rap(r,p,m3,m4,pswt,*999)
          
c--- processes that use "gen4h"     
      elseif ( (case .eq. 'HWW_4l')
     .    .or. (case .eq. 'HZZ_4l') ) then
        npart=4
        call gen4h(r,p,pswt,*999)
          
c--- processes that use "gen5" 
      elseif ( (case .eq. 'W_twdk') .or. (case .eq. 'vlchwt') ) then 
        npart=5 
        call gen5(r,p,pswt,*999)
    
c--- processes that use "gen6"     
      elseif ( (case .eq. 'tt_bbl')
     .    .or. (case .eq. 'tt_bbh')
     .    .or. (case .eq. 'tautau')
     .    .or. (case .eq. 'Wtbwdk')
     .    .or. (case .eq. 'vlchk6')
     .    .or. (case .eq. 'vlchwg')
     .    .or. (case .eq. 'vlchwh') ) then
        npart=6
        call gen6(r,p,pswt,*999)
        
c--- processes that use "gen7"     
      elseif ( (case .eq. 'qq_ttg') ) then
        npart=7
        call gen7(r,p,pswt,*999)

c--- processes that use "gen8"     
      elseif ( (case .eq. 'qq_tth') 
     .    .or. (case .eq. 'qq_ttz') 
     .    .or. (case .eq. 'vlchk8') ) then
        npart=8
        call gen8(r,p,pswt,*999)
          
c--- processes that use "gen_njets" with an argument of "1"     
      elseif ( (case .eq. 'W_1jet')
     .    .or. (case .eq. 'Wcjet0')
     .    .or. (case .eq. 'Z_1jet')
     .    .or. (case .eq. 'H_1jet')
     .    .or. (case .eq. 'httjet')
     .    .or. (case .eq. 'ggfus1')
     .    .or. (case .eq. 'attjet')
     .    .or. (case .eq. 'gQ__ZQ') ) then
        npart=3
        call gen_njets(r,1,p,pswt,*999)
         
c--- processes that use "gen_njets" with an argument of "2"
      elseif ( (case .eq. 'Wbbbar')
     .    .or. (case .eq. 'W_2jet')
     .    .or. (case .eq. 'Z_2jet')
     .    .or. (case .eq. 'Zbbbar')
     .    .or. (case .eq. 'qq_Hqq')
     .    .or. (case .eq. 'ggfus2')
     .    .or. (case .eq. 'W_bjet')
     .    .or. (case .eq. 'Wcjetg')
     .    .or. (case .eq. 'Z_bjet') ) then
        npart=4
        call gen_njets(r,2,p,pswt,*999)         
        
c--- processes that use "gen_njets" with an argument of "3"
      elseif ( (case .eq. 'W_3jet') 
     .    .or. (case .eq. 'Wbbjet') 
     .    .or. (case .eq. 'Z_3jet') 
     .    .or. (case .eq. 'Zbbjet') 
     .    .or. (case .eq. 'Wb2jet') 
     .    .or. (case .eq. 'qqHqqg')
     .    .or. (case .eq. 'ggfus3')
     .    .or. (case .eq. 'Zbjetg')
     .    .or. (case .eq. 'vlchk5') ) then
        npart=5
        call gen_njets(r,3,p,pswt,*999)      
 
c--- processes that use "gen_stop" with an argument of "1" (number of extra jets)
      elseif ( (case .eq. 'bq_tpq')
     .    .or. (case .eq. 'ttdkay')
     .    .or. (case .eq. 't_bbar')
     .    .or. (case .eq. 'tdecay') ) then
        npart=4
        call gen_stop(r,1,p,pswt,*999)

c--- DEFAULT: processes that use "gen4"
      else
        if ((case .eq. 'vlchk4') .or. (case .eq. 'vlchkm')) then
          wsqmin=0d0
          wsqmax=sqrts**2
        endif
        npart=4
        call gen4(r,p,pswt,*999)      
      endif
     
      nvec=npart+2
      call dotem(nvec,p,s)
      
c---impose cuts on final state
      if (case .ne. 'vlchk6' .and. case .ne. 'tautau') then
        call masscuts(s,*999)
      endif
c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)                                             

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif
      
      if (dynamicscale) call scaleset(rscalestart,fscalestart,p)
      
      xx(1)=-2d0*p(1,4)/sqrts
      xx(2)=-2d0*p(2,4)/sqrts

      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx(1),xx(2)      

c--- Calculate the required matrix elements      
      if     (case .eq. 'W_only') then
        call qqb_w(p,msq)
      elseif (case .eq. 'W_1jet') then
        call qqb_w_g(p,msq)
      elseif (case .eq. 'Wgamma') then
        call qqb_wgam(p,msq)
      elseif (case .eq. 'W_cjet') then
        call qqb_w_cjet(p,msq)
      elseif (case .eq. 'Wcjet0') then
        call qqb_w_cjet_massless(p,msq)
      elseif (case .eq. 'Wbbmas') then
        call qqb_wbbm(p,msq)
      elseif (case .eq. 'Wbbbar') then
        call qqb_wbb(p,msq)
      elseif (case .eq. 'W_2jet') then
        call qqb_w2jet(p,msq)
      elseif (case .eq. 'W_3jet') then
        call qqb_w2jet_g(p,msq)
      elseif (case .eq. 'Wbbjet') then
        call qqb_wbb_g(p,msq)
      elseif (case .eq. 'Z_only') then
        call qqb_z(p,msq)
      elseif (case .eq. 'Z_1jet') then
        call qqb_z1jet(p,msq)
      elseif (case .eq. 'Z_2jet') then
        call qqb_z2jet(p,msq)
      elseif (case .eq. 'Z_3jet') then
        call qqb_z2jet_g(p,msq)
      elseif (case .eq. 'Zgamma') then
        call qqb_zgam(p,msq)
      elseif (case .eq. 'Zbbmas') then
        call qqb_zbbm(p,msq)
      elseif (case .eq. 'Zbbbar') then
        call qqb_zbb(p,msq)
      elseif (case .eq. 'Zbbjet') then
        call qqb_zbb_g(p,msq)
      elseif (case .eq. 'WWqqbr') then
        call qqb_ww(p,msq)
      elseif (case .eq. 'WWnpol') then
        call qqb_ww_unpol(p,msq)
      elseif (case .eq. 'WZbbar') then
        call qqb_wz(p,msq)
      elseif (case .eq. 'ZZlept') then
        call qqb_zz(p,msq)
      elseif (case .eq. 'WHbbar') then
        call qqb_wh(p,msq)
      elseif (case .eq. 'ZHbbar') then
        call qqb_zh(p,msq)
      elseif (case .eq. 'ggfus0') then
        call gg_h(p,msq)
      elseif (case .eq. 'HWW_4l') then
        call qqb_hww(p,msq)
      elseif (case .eq. 'HZZ_4l') then
        call qqb_hzz(p,msq)
      elseif (case .eq. 'H_1jet') then
        call qqb_hg(p,msq)
      elseif ((case .eq. 'tt_bbl') .or. (case .eq. 'tt_bbh')) then
        call qqb_ttb(p,msq)
      elseif (case .eq. 'qq_ttg') then
        call qqb_ttb_g(p,msq)
      elseif (case .eq. 'tt_tot') then
        call qqb_QQb(p,msq)
      elseif (case .eq. 'bb_tot') then
        call qqb_QQb(p,msq)
      elseif (case .eq. 'cc_tot') then
        call qqb_QQb(p,msq)
      elseif (case .eq. 'tt_glu') then
        call qqb_QQb_g(p,msq)
      elseif (case .eq. 'bq_tpq') then
        call bq_tpq(p,msq)
      elseif (case .eq. 'ttdkay') then
        write(6,*) 'This process is not a leading order contribution'
        stop
      elseif (case .eq. 't_bbar') then
        call qqb_tbb(p,msq)
      elseif (case .eq. 'tdecay') then
        write(6,*) 'This process is not a leading order contribution'
        stop
      elseif (case .eq. 'W_tndk') then
        call qqb_w_tndk(p,msq)
      elseif (case .eq. 'W_twdk') then
        call qqb_w_twdk(p,msq)
      elseif (case .eq. 'Wtbwdk') then
        call qqb_wtbwdk(p,msq)
      elseif (case .eq. 'Wtbndk') then
        call qqb_wtbndk(p,msq)
      elseif (case .eq. 'tottth') then
        call qqb_tottth(p,msq)
      elseif (case .eq. 'qq_tth') then
        call qqb_tth(p,msq)
      elseif (case .eq. 'qq_ttz') then
        call qqb_ttz(p,msq)
      elseif (case .eq. 'httjet') then
        call qqb_higgs(p,msq)
      elseif (case .eq. 'ggfus1') then
        call gg_hg(p,msq)
      elseif (case .eq. 'attjet') then
        call qqb_higgs_odd(p,msq)
      elseif (case .eq. 'qq_Hqq') then
        call VV_hqq(p,msq)
      elseif (case .eq. 'qqHqqg') then
        call VV_hqq_g(p,msq)
      elseif (case .eq. 'tautau') then
        call qqb_tautau(p,msq)
      elseif (case .eq. 'gQ__ZQ') then
        call gQ_zQ(p,msq)
      elseif (case .eq. 'Zccmas') then
        call qqb_zccm(p,msq)
      elseif (case .eq. 'ggfus2') then
        call gg_hgg(p,msq)
      elseif (case .eq. 'ggfus3') then
        call gg_hggg(p,msq)
      elseif (case .eq. 'W_bjet') then
        call qqb_wbjet(p,msq)
      elseif (case .eq. 'Wcjetg') then
        call qqb_w_cjet_massless_g(p,msq)
      elseif (case .eq. 'Z_bjet') then
        call qqb_zbjet(p,msq)
      elseif (case .eq. 'Zbjetg') then
        call qqb_zbjet_g(p,msq)
      elseif (case .eq. 'vlchk2') then
        call qqb_vol(p,msq)
        flux=one/vol(W,2)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=1d0
        fx2(-1)=1d0
      elseif (case .eq. 'vlchk3') then
        call qqb_vol(p,msq)
        flux=one/vol(W,3)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=2d0
        fx2(-1)=2d0
      elseif (case .eq. 'vlchk4') then
        taumin=0.0001d0
        bbsqmax=W
        bbsqmin=0d0
        call qqb_vol(p,msq)
        flux=one/vol(W,4)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=4d0*xx(1)
        fx2(-1)=2d0/xx(2)
      elseif (case .eq. 'vlchk5') then
        taumin=0.0001d0
        bbsqmax=W
        bbsqmin=0d0
        call qqb_vol(p,msq)
        flux=one/vol(W,5)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=4d0
        fx2(-1)=4d0
      elseif (case .eq. 'vlchk6') then
        call qqb_vol(p,msq)
        flux=one/vol(W,6)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=6d0*xx(1)
        fx2(-1)=4d0/xx(2)
      elseif (case .eq. 'vlchk8') then
        call qqb_vol(p,msq)
        flux=one/vol(W,8)
        bbsqmax=W
        bbsqmin=0d0
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=6d0/xx(1)
        fx2(-1)=6d0/xx(2)
      elseif (case .eq. 'vlchkm') then
        taumin=0.0001d0
        bbsqmax=W
        bbsqmin=0d0
        call qqb_vol(p,msq)
        flux=one/vol_mass(mb,W)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=4d0*xx(1)
        fx2(-1)=2d0/xx(2)      
      elseif (case .eq. 'vlchm3') then
        taumin=(2d0*mt/sqrts)**2
        call qqb_vol(p,msq)
        flux=one/vol3_mass(mt,W)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=2d0
        fx2(-1)=2d0
      elseif ((case .eq. 'vlchwt') .or. (case .eq. 'vlchwn')
     .   .or. (case .eq. 'vlchwg') .or. (case .eq. 'vlchwh')) then
        taumin=0.0001d0
        call qqb_vol(p,msq)
        flux=one/vol_wt(W)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=1d0
        fx2(-1)=1d0
      else
        write(6,*) 'Unimplemented process in lowint : case=',case
        stop 
      endif
      

      do j=-1,1
      do k=-1,1
      xmsq_bypart(j,k)=0d0
      enddo
      enddo 


      currentPDF=0

c--- do not calculate the flux if we're only checking the volume      
      if (case(1:4) .ne. 'vlch') then      
        flux=fbGeV2/(2d0*xx(1)*xx(2)*W)
      endif
      
c--- initialize a PDF set here, if calculating errors
  777 continue    
      xmsq=0d0
      if (PDFerrors) then
        call InitPDF(currentPDF)
      endif
      
c--- calculate PDF's  
      if (case(1:4) .ne. 'vlch') then      
        call fdist(ih1,xx(1),facscale,fx1)
        call fdist(ih2,xx(2),facscale,fx2)
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav    
      
      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
      endif
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif

c === fastNLO
      fncoeff(j,k) = msq(j,k)
      fnpdf(j,k) = fx1(j)*fx2(k)
      fnx(1) = xx(1)
      fnx(2) = xx(2)
c ==============================
      
      xmsqjk=fx1(j)*fx2(k)*msq(j,k)
      xmsq=xmsq+xmsqjk
      
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

      if (currentPDF .eq. 0) then
        xmsq_bypart(sgnj,sgnk)=xmsq_bypart(sgnj,sgnk)+xmsqjk
      endif
      
 20   continue
      enddo
      enddo

      if (currentPDF .eq. 0) then
        lowint=flux*pswt*xmsq/BrnRat
      endif
      
c--- loop over all PDF error sets, if necessary
      if (PDFerrors) then
        PDFwgt(currentPDF)=flux*pswt*xmsq/BrnRat*wgt/itmx
        PDFxsec(currentPDF)=PDFxsec(currentPDF)
     .     +PDFwgt(currentPDF)
        currentPDF=currentPDF+1
        if (currentPDF .le. maxPDFsets) goto 777
      endif    

      if (creatent) then
        wt_gg=xmsq_bypart(0,0)*wgt*flux*pswt/BrnRat/dfloat(itmx)
        wt_gq=(xmsq_bypart(+1,0)+xmsq_bypart(-1,0)
     .        +xmsq_bypart(0,+1)+xmsq_bypart(0,-1)
     .        )*wgt*flux*pswt/BrnRat/dfloat(itmx)
        wt_qq=(xmsq_bypart(+1,+1)+xmsq_bypart(-1,-1)
     .        )*wgt*flux*pswt/BrnRat/dfloat(itmx)
        wt_qqb=(xmsq_bypart(+1,-1)+xmsq_bypart(-1,+1)
     .        )*wgt*flux*pswt/BrnRat/dfloat(itmx)
      endif

      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)

      do j=-1,1
      do k=-1,1
        lord_bypart(j,k)=lord_bypart(j,k)+
     .       wgt*flux*pswt*xmsq_bypart(j,k)/BrnRat
      enddo
      enddo


      val=lowint*wgt
      val2=lowint**2*wgt
c--- update the maximum weight so far, if necessary
c---  but not if we are already unweighting ...
      if ((.not.unweight) .and. (dabs(val) .gt. wtmax)) then
        wtmax=dabs(val)
      endif

      if (bin) then
c === fastNLO
        fnsigma = flux*pswt/BrnRat*wgt/dfloat(itmx)
        ifnctrb = 1
c ===================================================
        val=val/dfloat(itmx)
        val2=val2/dfloat(itmx)**2
c ---   DSW. If the user has not selected to generate
c ---   events, still call nplotter here in order to
c ---   fill histograms/ntuples with weighted events :
        if (.not.evtgen) then
          call nplotter(pjet,val,val2,0)
        endif
      endif

c --- Check weights :
      if (unweight) then
c       write(6,*) 'Test weight ',val,' against max ',wtmax
        wtabs = dabs(val)
        if (ran2() .lt. (wtabs/wtmax)) then
c         write(6,*) 'Keep event with weight',val
          if (wtabs.lt.wtmax) then
            newwt = 1d0
          else
            newwt = wtabs/wtmax
          endif
          if (newwt .gt. 1.0d0) then
            write(6,*) 'WARNING : lowint : event with |weight| > 1.',
     +            ' |weight| = ',newwt
          endif
c ---     just in case the weight was negative :
          newwt = newwt*dsign(1d0,val)
          call nplotter(pjet,newwt,newwt,0)
c ---     DSW. If I'm storing the event, I need to make a decision
c ---     about the flavours :
          call decide_flavour(pflav,pbarflav)
          call storeevent(pjet,newwt,pflav,pbarflav)
        endif
      endif

      return

 999  continue
      ntotzero=ntotzero+1
      
      return
      end


