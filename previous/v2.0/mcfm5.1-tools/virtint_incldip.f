      double precision function virtint(r,wgt)
      implicit none
      include 'constants.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'npart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'agq.f'
      include 'PR_new.f'
      include 'PR_cs_new.f'
      include 'PR_h2j.f'
      include 'PR_twojet.f'
      include 'PR_stop.f'
      include 'msq_cs.f'
      include 'msq_struc.f'
      include 'new_pspace.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'clustering.f'
      include 'efficiency.f'
      include 'lc.f'
      include 'process.f'
      include 'maxwt.f'
      include 'limits.f'
      include 'heavyflav.f'
      include 'nflav.f'
      include 'b0.f'
      include 'PDFerrors.f'
      include 'masses.f'
      include 'wts_bypart.f'
      include 'nores.f'
c     === fastNLO
      include 'fastnlo.f'
c     ============================================
      double precision mqq(0:2,fn:nf,fn:nf)
      double precision msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      double precision msqx_cs(0:2,-nf:nf,-nf:nf)
      double precision AP(-1:1,-1:1,3)

      integer ih1,ih2,j,k,m,n,cs,ics,csmax,nvec,is,iq,ia,ib,ic
      double precision p(mxpart,4),pjet(mxpart,4),r(mxdim),W,sqrts,xmsq,
     .     val,val2,fx1(-nf:nf),fx2(-nf:nf),fx1z(-nf:nf),fx2z(-nf:nf),xmsqt
      double precision pswt,xjac,rscalestart,fscalestart,
     .     wgt,msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),msqvdk(-nf:nf,-nf:nf),
     .     msq_qq,msq_aa,msq_aq,msq_qa,msq_qg,msq_gq,epcorr
      double precision xx(2),z,x1onz,x2onz,flux,omz,
     .     BrnRat,xmsq_old,tmp
      double precision xmsq_bypart(-1:1,-1:1),lord_bypart(-1:1,-1:1)
      integer nshot,rvcolourchoice,sgnj,sgnk
      logical bin,first,includedipole
      logical creatent,dswhisto
      character*4 mypart
      character*30 runstring
c     double complex ee(mxpart,mxpart),ek(mxpart,mxpart)
      common/runstring/runstring
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/x1x2/xx
      common/BrnRat/BrnRat
      common/rvcolourchoice/rvcolourchoice
      common/bypart/lord_bypart
      common/mypart/mypart
      common/outputflags/creatent,dswhisto
      data p/48*0d0/
      data nshot/1/
      data first/.true./
      save first,rscalestart,fscalestart
      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif

      ntotshot=ntotshot+1
      virtint=0d0

      W=sqrts**2

c---  processes that use "gen2"     
      if     ( (case .eq. 'W_only')
     .     .or. (case .eq. 'Z_only')
     .     .or. (case .eq. 'ggfus0')
     .     .or. (case .eq. 'vlchk2') ) then
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

c---  processes that use "gen2m"     
      elseif ( (case .eq. 'tt_tot')
     .        .or. (case .eq. 'bb_tot')
     .        .or. (case .eq. 'cc_tot') ) then
         npart=2
         call gen2m(r,p,pswt,*999)
         
c---  processes that use "gen3"     
      elseif ( (case .eq. 'W_cjet') 
     .        .or.  (case .eq. 'W_tndk') ) then
         npart=3
         call gen3(r,p,pswt,*999)

c---  processes that use "gen3jet"     
      elseif ( (case .eq. 'Wgamma') 
     .        .or.  (case .eq. 'Zgamma')
     .        .or.  (case .eq. 'vlchk3') ) then
         if (case .eq. 'vlchk3') then
            wsqmin=0d0
            wsqmax=sqrts**2
         endif
         npart=3
         call gen3jet(r,p,pswt,*999)

c---  processes that use "gen4h"     
      elseif ( (case .eq. 'HWW_4l')
     .        .or. (case .eq. 'HZZ_4l') ) then
         npart=4
         call gen4h(r,p,pswt,*999)
         
c---  processes that use "gen5" 
      elseif ( (case .eq. 'W_twdk') .or. (case .eq. 'Wtdkay') ) then 
         npart=5 
         call gen5(r,p,pswt,*999)
         
c---  processes that use "gen6"     
      elseif ( (case .eq. 'tt_bbl')
     .        .or. (case .eq. 'tt_bbh')
     .        .or. (case .eq. 'tautau')
     .        .or. (case .eq. 'vlchk6') ) then
         npart=6
         call gen6(r,p,pswt,*999)

c---  processes that use "gen7"     
      elseif ( (case .eq. 'qq_ttg') ) then
         npart=7
         call gen7(r,p,pswt,*999)

c---  processes that use "gen8"     
      elseif ( (case .eq. 'qq_tth') 
     .        .or. (case .eq. 'qq_ttz') 
     .        .or. (case .eq. 'vlchk8') ) then
         npart=8
         call gen8(r,p,pswt,*999)
         
c---  processes that use "gen_njets" with an argument of "1"     
      elseif ( (case .eq. 'W_1jet')
     .        .or. (case .eq. 'Wcjet0')
     .        .or. (case .eq. 'Z_1jet')
     .        .or. (case .eq. 'H_1jet')
     .        .or. (case .eq. 'httjet')
     .        .or. (case .eq. 'ggfus1')
     .        .or. (case .eq. 'attjet')
     .        .or. (case .eq. 'gQ__ZQ') ) then
         npart=3
         call gen_njets(r,1,p,pswt,*999)
         
c---  processes that use "gen_njets" with an argument of "2"
      elseif ( (case .eq. 'Wbbbar')
     .        .or. (case .eq. 'W_2jet')
     .        .or. (case .eq. 'Z_2jet')
     .        .or. (case .eq. 'Zbbbar')
     .        .or. (case .eq. 'qq_Hqq')
     .        .or. (case .eq. 'ggfus2')
     .        .or. (case .eq. 'W_bjet')
     .        .or. (case .eq. 'Wcjetg')
     .        .or. (case .eq. 'Z_bjet') ) then
         npart=4
         call gen_njets(r,2,p,pswt,*999)         
         
c---  processes that use "gen_njets" with an argument of "3"
      elseif ( (case .eq. 'W_3jet') 
     .        .or. (case .eq. 'Wbbjet') 
     .        .or. (case .eq. 'Z_3jet') 
     .        .or. (case .eq. 'Zbbjet') 
     .        .or. (case .eq. 'Wb2jet') 
     .        .or. (case .eq. 'qqHqqg')
     .        .or. (case .eq. 'vlchk5') ) then
         npart=5
         call gen_njets(r,3,p,pswt,*999)      
         
c---  processes that use "gen_stop" with an argument of "1"
      elseif ( (case .eq. 'bq_tpq')
     .        .or. (case .eq. 'ttdkay')
     .        .or. (case .eq. 't_bbar')
     .        .or. (case .eq. 'tdecay') ) then
         npart=4
         call gen_stop(r,1,p,pswt,*999)

c---  DEFAULT: processes that use "gen4"
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

c---  impose mass cuts on final state
      call masscuts(s,*999)
c---- reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
      
c---  see whether this point will pass cuts - if it will not, do not
c---  bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
         goto 999
      endif
      
      if (dynamicscale) call scaleset(rscalestart,fscalestart,p)
      
      z=r(ndim)**2
      if (nshot .eq. 1) z=0.95d0
      xjac=two*dsqrt(z)
      
      omz=1d0-z

      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)

c---  to test poles, we need colourchoice=0, but save real value
      if (nshot .eq. 1) then
         rvcolourchoice=colourchoice
         colourchoice=0
      endif
      
 12   continue
c---  point to restart from when checking epsilon poles

c---  correction to epinv from AP subtraction when mu_FAC != mu_REN,
c---  corresponding to subtracting -1/epinv*Pab*log(musq_REN/musq_FAC)
      epcorr=epinv+2d0*dlog(scale/facscale)
      
c---  for the case of virtual correction in the top quark decayy,
c---  ('tdecay','ttdkay','Wtdkay') there are no extra initial-state
c---  contributions, so all these should be set to zero
      if  ( (case .eq. 'tdecay') .or. (case .eq. 'ttdkay') 
     .     .or. (case .eq. 'Wtdkay') ) epcorr=0d0      

      AP(q,q,1)=+ason2pi*Cf*1.5d0*epcorr
      AP(q,q,2)=+ason2pi*Cf*(-1d0-z)*epcorr
      AP(q,q,3)=+ason2pi*Cf*2d0/omz*epcorr
      AP(a,a,1)=+ason2pi*Cf*1.5d0*epcorr
      AP(a,a,2)=+ason2pi*Cf*(-1d0-z)*epcorr
      AP(a,a,3)=+ason2pi*Cf*2d0/omz*epcorr

      AP(q,g,1)=0d0
      AP(q,g,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
      AP(q,g,3)=0d0
      AP(a,g,1)=0d0
      AP(a,g,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
      AP(a,g,3)=0d0

c---  modifications for running with mb>0
      if ( ((case .eq. 'W_twdk') .or. (case .eq. 'W_tndk'))
     .     .and. (runstring(1:4) .eq. 'mass')) then
         AP(q,g,2)=-ason2pi*Tr*(z**2+omz**2)*dlog(facscale**2/mbsq)
         AP(a,g,2)=AP(q,g,2)
      endif
      
      if ( ((case .eq. 'W_twdk') .or. (case .eq. 'W_tndk'))
     .     .and. (nores) ) then
         AP(q,g,2)=0d0
         AP(a,g,2)=AP(q,g,2)      
      endif
      
      AP(g,q,1)=0d0
      AP(g,q,2)=ason2pi*Cf*(1d0+omz**2)/z*epcorr
      AP(g,q,3)=0d0
      AP(g,a,1)=0d0
      AP(g,a,2)=ason2pi*Cf*(1d0+omz**2)/z*epcorr
      AP(g,a,3)=0d0

      AP(g,g,1)=+ason2pi*b0*epcorr
      AP(g,g,2)=+ason2pi*xn*2d0*(1d0/z+z*omz-2d0)*epcorr
      AP(g,g,3)=+ason2pi*xn*2d0/omz*epcorr

      if (case .eq. 'bq_tpq') then
         do ia=-1,+2
            do ib=-1,+2
               do ic=-1,+2
                  do is=1,3
                     B1(ia,ib,ic,is)=0d0
                     B2(ia,ib,ic,is)=0d0
                  enddo
               enddo
            enddo
         enddo
      endif
      
      do ia=-1,+1
         do ib=-1,+1
            do ic=-1,+1
               do is=1,3
                  Q1(ia,ib,ic,is)=0d0
                  Q2(ia,ib,ic,is)=0d0
                  do cs=1,8
                     H1(ia,ib,ic,cs,is)=0d0
                     H2(ia,ib,ic,cs,is)=0d0
                  enddo
                  do cs=0,2
                     R1(ia,ib,ic,cs,is)=0d0
                     R2(ia,ib,ic,cs,is)=0d0
                     do j=1,8
                        S1(ia,ib,ic,j,cs,is)=0d0
                        S2(ia,ib,ic,j,cs,is)=0d0
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      
c---  Calculate the required matrix elements      
      if     (case .eq. 'W_only') then
         call qqb_w(p,msq)
         call qqb_w_v(p,msqv)
         call qqb_w_z(p,z)
      elseif (case .eq. 'W_1jet') then
         call qqb_w_g(p,msq)
         call qqb_w1jet_v(p,msqv)
         call qqb_w1jet_z(p,z)
      elseif (case .eq. 'Wgamma') then
         call qqb_wgam(p,msq)
         call qqb_wgam_v(p,msqv)
         call qqb_wgam_z(p,z)
      elseif (case .eq. 'W_cjet') then
         call qqb_w_cjet(p,msq)
         write(6,*) 'Process not yet fully implemented'
         stop
c     call qqb_w_cjet_v(p,msqv)
c     call qqb_w_cjet_z(p,z)
      elseif (case .eq. 'Wbbbar') then
         call qqb_wbb(p,msq)
         call qqb_wbb_v(p,msqv)
         call qqb_wbb_z(p,z)
      elseif (case .eq. 'W_2jet') then         
         call qqb_w2jetx(p,msq,mqq,msqx,msqx_cs)
         call qqb_w2jet_v(p,msqv)
         call qqb_w2jet_z(p,z)
      elseif (case .eq. 'Z_only') then
         call qqb_z(p,msq)
         call qqb_z_v(p,msqv)
         call qqb_z_z(p,z)
      elseif (case .eq. 'Z_1jet') then
         call qqb_z1jet(p,msq)
         call qqb_z1jet_v(p,msqv)
         call qqb_z1jet_z(p,z)
      elseif (case .eq. 'Z_2jet') then
         call qqb_z2jetx(p,msq,mqq,msqx,msqx_cs)
         call qqb_z2jet_v(p,msqv)
         call qqb_z2jet_z(p,z)
      elseif (case .eq. 'Zgamma') then
         call qqb_zgam(p,msq)
         call qqb_zgam_v(p,msqv)
         call qqb_zgam_z(p,z)
      elseif (case .eq. 'Zbbbar') then
         call qqb_zbb(p,msq)
         call qqb_zbb_v(p,msqv)
         call qqb_zbb_z(p,z)
      elseif (case .eq. 'WWqqbr') then
         call qqb_ww(p,msq)
         call qqb_ww_v(p,msqv)
         call qqb_ww_z(p,z)
      elseif (case .eq. 'WZbbar') then
         call qqb_wz(p,msq)
         call qqb_wz_v(p,msqv)
         call qqb_wz_z(p,z)
      elseif (case .eq. 'ZZlept') then
         call qqb_zz(p,msq)
         call qqb_zz_v(p,msqv)
         call qqb_zz_z(p,z)
      elseif (case .eq. 'WHbbar') then
         call qqb_wh(p,msq)
         call qqb_wh_v(p,msqv)
         call qqb_wh_z(p,z)
      elseif (case .eq. 'ZHbbar') then
         call qqb_zh(p,msq)
         call qqb_zh_v(p,msqv)
         call qqb_zh_z(p,z)
      elseif (case .eq. 'ggfus0') then
         call gg_h(p,msq)
         call gg_h_v(p,msqv)
         call gg_h_z(p,z)
      elseif (case .eq. 'HWW_4l') then
         call qqb_hww(p,msq)
         call qqb_hww_v(p,msqv)
         call qqb_hww_z(p,z)
      elseif (case .eq. 'HZZ_4l') then
         call qqb_hzz(p,msq)
         call qqb_hzz_v(p,msqv)
         call qqb_hzz_z(p,z)
      elseif (case .eq. 'H_1jet') then
         call qqb_hg(p,msq)
         call qqb_hg_v(p,msqv)
         call qqb_hg_z(p,z)
      elseif ((case .eq. 'tt_tot')
     .        .or. (case .eq. 'bb_tot')
     .        .or. (case .eq. 'cc_tot')) then
         call qqb_QQb(p,msq)
         call qqb_QQb_v(p,msqv)
         call qqb_QQb_z(p,z)
      elseif (case .eq. 'bq_tpq') then
         call bq_tpq(p,msq)
         call bq_tpq_v(p,msqv)
         call bq_tpq_z(p,z)
         if (mypart .eq. 'todk') then
            call bq_tpq_vdk(p,msqvdk)
            do j=-nf,nf
               do k=-nf,nf
                  msqv(j,k)=msqv(j,k)+msqvdk(j,k)
               enddo
            enddo
         endif
      elseif (case .eq. 'ttdkay') then
         do j=-nf,nf
            do k=-nf,nf
               msq(j,k)=0d0
            enddo
         enddo
         call bq_tpq_vdk(p,msqv)
      elseif (case .eq. 't_bbar') then
         call qqb_tbb(p,msq)
         call qqb_tbb_v(p,msqv)
         call qqb_tbb_z(p,z)
         if (mypart .eq. 'todk') then
            call qqb_tbb_vdk(p,msqvdk)
            do j=-nf,nf
               do k=-nf,nf
                  msqv(j,k)=msqv(j,k)+msqvdk(j,k)
               enddo
            enddo
         endif
      elseif (case .eq. 'tdecay') then
         do j=-nf,nf
            do k=-nf,nf
               msq(j,k)=0d0
            enddo
         enddo
         call qqb_tbb_vdk(p,msqv)    
      elseif (case .eq. 'W_tndk') then
         call qqb_w_tndk(p,msq)    
         call qqb_w_tndk_v(p,msqv)    
         call qqb_w_tndk_z(p,z)    
      elseif (case .eq. 'W_twdk') then
         call qqb_w_twdk(p,msq)    
         call qqb_w_twdk_v(p,msqv)    
         call qqb_w_twdk_z(p,z)    
         if (mypart .eq. 'todk') then
            call qqb_w_twdk_vdk(p,msqvdk)
            do j=-nf,nf
               do k=-nf,nf
                  msqv(j,k)=msqv(j,k)+msqvdk(j,k)
               enddo
            enddo
         endif
      elseif (case .eq. 'Wtdkay') then
         do j=-nf,nf
            do k=-nf,nf
               msq(j,k)=0d0
            enddo
         enddo
         call qqb_w_twdk_vdk(p,msqv)
      elseif (case .eq. 'ggfus1') then
         call gg_hg(p,msq)
         call gg_hg_v(p,msqv)
         call gg_hg_z(p,z)
      elseif (case .eq. 'qq_Hqq') then
         call VV_hqq(p,msq)
         call VV_hqq_v(p,msqv)
         call VV_hqq_z(p,z)
      elseif (case .eq. 'gQ__ZQ') then
         call gQ_zQ(p,msq)
         call gQ_zQ_v(p,msqv)
         call gQ_zQ_z(p,z)
      elseif (case .eq. 'Z_bjet') then
         call qqb_zbjet(p,msq)
         call qqb_zbjet_v(p,msqv)
         call qqb_zbjet_z(p,z)
      endif
      
C---  initialize to zero
      do j=-1,1
         do k=-1,1
            xmsq_bypart(j,k)=0d0
         enddo
      enddo

      currentPDF=0
      
c---  initialize a PDF set here, if calculating errors
 777  continue    
      xmsq=0d0
      if (PDFerrors) then
         call InitPDF(currentPDF)
      endif
      call fdist(ih1,xx(1),facscale,fx1)
      call fdist(ih2,xx(2),facscale,fx2)

      do j=-nf,nf
         fx1z(j)=0d0
         fx2z(j)=0d0
      enddo
      
      if (z .gt. xx(1)) then
         x1onz=xx(1)/z
         call fdist(ih1,x1onz,facscale,fx1z)
      endif
      if (z .gt. xx(2)) then
         x2onz=xx(2)/z
         call fdist(ih2,x2onz,facscale,fx2z)
      endif         

c === fastNLO - reset coeff. arrays
      do j=-nflav,nflav
         do k=-nflav,nflav
            fncoeff(j,k) = 0d0
            fncfv2(j,k) = 0d0
            fncfv3(j,k) = 0d0
         Enddo
      Enddo
c =========================================

      do j=-nflav,nflav
         do k=-nflav,nflav
            
c     === fastNLO
            fnpdf(j,k) = fx1(j)*fx2(k)
            fnx(1) = xx(1)
            fnx(2) = xx(2)
            fnpdf2(j,k) = fx1z(j)*fx2(k)
            fnpdf3(j,k) = fx1(j)*fx2z(k)
            fnxz(1) = x1onz
            fnxz(2) = x2onz
c     ==============================
            
            if (ggonly) then
               if ((j.ne.0) .or. (k.ne.0)) goto 20
            endif

            if (gqonly) then
               if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
            endif      
            
            if (noglue) then 
               if ((j.eq.0) .or. (k.eq.0)) goto 20
            endif

            tmp=xmsq

c---  The variables R1 and R2 provide the Regular and Plus pieces associated
c---  with radiation from leg 1 (R1(a,b,c,cs,is)) and leg 2 (R2(a,b,c,cs,is))
c---  In each case the parton labelling is using the normal QM notation of 
c---  putting everything backward
c---  emitted line after emission =    a
c---  emitter before emission     =    b
c---  spectator                   =    c
c---  There is no label for he or she who is emitted.
c---  Note that in general each piece will be composed of many different
c---  dipole contributions

c---  SUM BY COLOUR STRUCTURES: H+2jets only
            if     (case .eq. 'ggfus2')then
               write(*,*) "fastNLO case ggfus2"
               xmsq=xmsq+fx1(j)*fx2(k)*(
     .              msqv(j,k)+msq(j,k))
c     write(6,*) j,k,'-> msqv = ',fx1(j)*fx2(k)*(
c     . msqv(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k))
c     tmp=xmsq

c---  quark-quark or antiquark-antiquark
               if (  ((j .gt. 0).and.(k .gt. 0))
     .              .or. ((j .lt. 0).and.(k .lt. 0))) then
c     write(6,*) 'virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
                  if (j. eq. k) then
                     m=+1
                     n=+1
                     csmax=6
                  else
                     m=abs(j)
                     n=abs(k)
                     csmax=6
                  endif
                  xmsqt=0d0
                  do cs=1,csmax
                     xmsqt=xmsqt
     &                    +msq_struc(cs,m,n)*(AP(q,q,1)-AP(q,q,3)
     &                    +H1(q,q,q,cs,1)-H1(q,q,q,cs,3)
     &                    +AP(q,q,1)-AP(q,q,3)
     &                    +H2(q,q,q,cs,1)-H2(q,q,q,cs,3))*fx1(j)*fx2(k)
     &                    +(msq_struc(cs,m,n)*(AP(q,q,2)+AP(q,q,3)
     &                    +H1(q,q,q,cs,2)+H1(q,q,q,cs,3))
     &                    +msq_struc(cs,g,k)*(AP(g,q,2)+H1(g,q,q,cs,2)))*fx1z(j)/z*fx2(k)
     &                    +(msq_struc(cs,m,n)*(AP(q,q,2)+AP(q,q,3)
     &                    +H2(q,q,q,cs,2)+H2(q,q,q,cs,3))
     &                    +msq_struc(cs,j,g)*(AP(g,q,2)+H2(g,q,q,cs,2)))*fx1(j)*fx2z(k)/z
                  enddo      
                  xmsq=xmsq+xmsqt
c     write(6,*) 'virtint: j,k,  ct=',j,k,xmsqt

c---  quark-antiquark or antiquark-quark
               elseif (  ((j .gt. 0).and.(k .lt. 0))
     .                 .or. ((j .lt. 0).and.(k .gt. 0))) then
c     write(6,*) 'virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
                  if (j. eq. -k) then
                     m=+1
                     n=-1
                     csmax=7
                  else
                     m=abs(j)
                     n=-abs(k)
                     csmax=6
                  endif
                  xmsqt=0d0
                  do cs=1,csmax
c     if ((cs .gt. 3) .and. (cs .lt. 7)) goto 67
c     do cs=4,6
                     xmsqt=xmsqt
     &                    +msq_struc(cs,m,n)*(AP(q,q,1)-AP(q,q,3)
     &                    +H1(q,q,a,cs,1)-H1(q,q,a,cs,3)
     &                    +AP(a,a,1)-AP(a,a,3)
     &                    +H2(a,a,q,cs,1)-H2(a,a,q,cs,3))*fx1(j)*fx2(k)
     &                    +(msq_struc(cs,m,n)*(AP(q,q,2)+AP(q,q,3)
     &                    +H1(q,q,a,cs,2)+H1(q,q,a,cs,3))
     &                    + msq_struc(cs,g,k)*(AP(g,q,2)+H1(g,q,a,cs,2)))*fx1z(j)/z*fx2(k)
     &                    +(msq_struc(cs,m,n)*(AP(a,a,2)+AP(a,a,3)
     &                    +H2(a,a,q,cs,2)+H2(a,a,q,cs,3))
     &                    + msq_struc(cs,j,g)*(AP(g,a,2)+H2(g,a,q,cs,2)))*fx1(j)*fx2z(k)/z
c     67 continue   
                  enddo
                  xmsq=xmsq+xmsqt
c     write(6,*) 'virtint: j,k,  ct=',j,k,xmsqt

c---  gluon-gluon
               elseif ((j .eq. g) .and. (k .eq. g)) then
c     write(6,*) 'virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
                  xmsqt=0d0
c---  loop up to 3 to cancel poles from gggg
c     do cs=1,3
                  do cs=1,6
                     msq_qg=msq_struc(cs,+5,g)+msq_struc(cs,+4,g)+msq_struc(cs,+3,g)
     &                    +msq_struc(cs,+2,g)+msq_struc(cs,+1,g)
     &                    +msq_struc(cs,-5,g)+msq_struc(cs,-4,g)+msq_struc(cs,-3,g)
     &                    +msq_struc(cs,-2,g)+msq_struc(cs,-1,g)
                     msq_gq=msq_struc(cs,g,+5)+msq_struc(cs,g,+4)+msq_struc(cs,g,+3)
     &                    +msq_struc(cs,g,+2)+msq_struc(cs,g,+1)
     &                    +msq_struc(cs,g,-5)+msq_struc(cs,g,-4)+msq_struc(cs,g,-3)
     &                    +msq_struc(cs,g,-2)+msq_struc(cs,g,-1)
                     xmsqt=xmsqt
     &                    +msq_struc(cs,g,g)*(AP(g,g,1)-AP(g,g,3)
     &                    +H1(g,g,g,cs,1)-H1(g,g,g,cs,3)
     &                    +AP(g,g,1)-AP(g,g,3)
     &                    +H2(g,g,g,cs,1)-H2(g,g,g,cs,3))*fx1(g)*fx2(g)
     &                    +(msq_struc(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                    +H1(g,g,g,cs,2)+H1(g,g,g,cs,3))
     &                    +msq_qg*(AP(q,g,2)+H1(q,g,g,cs,2)))*fx1z(g)/z*fx2(g)
     &                    +(msq_struc(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                    +H2(g,g,g,cs,2)+H2(g,g,g,cs,3))
     &                    +msq_gq*(AP(q,g,2)+H2(q,g,g,cs,2)))*fx1(g)*fx2z(g)/z
                  enddo
                  xmsq=xmsq+xmsqt
c     write(6,*) 'virtint: j,k,  ct=',j,k,xmsqt

c---  quark-gluon and anti-quark gluon
               elseif ((j .ne. 0) .and. (k .eq. g)) then
c     write(6,*) 'virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
                  m=+1
                  n=0
                  xmsqt=0d0
                  do cs=4,6
                     xmsqt=xmsqt    
     &                    + msq_struc(cs,m,g)*(AP(q,q,1)-AP(q,q,3)
     &                    +H1(q,q,g,cs,1)-H1(q,q,g,cs,3)
     &                    +H2(g,g,q,cs,1)-H2(g,g,q,cs,3)
     &                    +AP(g,g,1)-AP(g,g,3))*fx1(j)*fx2(g)
     &                    +(msq_struc(cs,m,g)*(AP(q,q,2)+AP(q,q,3)
     &                    +H1(q,q,g,cs,2)+H1(q,q,g,cs,3))
     &                    +msq_struc(cs,g,g)*(AP(g,q,2)+H1(g,q,g,cs,2)))*fx1z(j)/z*fx2(g)
     &                    +(msq_struc(cs,m,g)*(AP(g,g,2)+AP(g,g,3)
     &                    +H2(g,g,q,cs,2)+H2(g,g,q,cs,3))
     &                    +msq_struc(cs,m,-m)*(AP(a,g,2)+H2(a,g,q,cs,2)))*fx1(j)*fx2z(g)/z
                  enddo
                  do cs=1,3
                     msq_qa=msq_struc(cs,j,-1)+msq_struc(cs,j,-2)+msq_struc(cs,j,-3)
     &                    +msq_struc(cs,j,-4)+msq_struc(cs,j,-5)
                     msq_qq=msq_struc(cs,j,+1)+msq_struc(cs,j,+2)+msq_struc(cs,j,+3)
     &                    +msq_struc(cs,j,+4)+msq_struc(cs,j,+5)
                     xmsqt=xmsqt    
     &                    +(msq_struc(cs,g,g)*(AP(g,q,2)+H1(g,q,g,cs,2)))*fx1z(j)/z*fx2(g)
     &                    +(+msq_qa*(AP(a,g,2)+H2(a,g,q,cs,2))
     &                    +msq_qq*(AP(q,g,2)+H2(q,g,q,cs,2)))*fx1(j)*fx2z(g)/z
                  enddo
                  xmsqt=xmsqt     
     &                 +msq_struc(iqr,m,-m)*(AP(a,g,2)+H2(a,g,q,iqr,2))*fx1(j)*fx2z(g)/z
                  xmsq=xmsq+xmsqt
c     write(6,*) 'virtint: j,k,  ct=',j,k,xmsqt
c     write(6,*) 'virtint: j,k, SUM=',j,k,xmsqt+fx1(j)*fx2(k)*msqv(j,k)
                  
c---  gluon-quark and gluon anti-quark
               elseif ((j .eq. 0) .and. (k .ne. 0)) then
c     write(6,*) 'virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
                  m=0
                  n=+1
                  xmsqt=0d0
                  do cs=4,6
                     xmsqt=xmsqt     
     &                    +msq_struc(cs,g,n)*(AP(g,g,1)-AP(g,g,3)
     &                    +H1(g,g,q,cs,1)-H1(g,g,q,cs,3)
     &                    +AP(q,q,1)-AP(q,q,3)
     &                    +H2(q,q,g,cs,1)-H2(q,q,g,cs,3))*fx1(g)*fx2(k)
     &                    +(msq_struc(cs,g,n)*(AP(g,g,2)+AP(g,g,3)
     &                    +H1(g,g,q,cs,2)+H1(g,g,q,cs,3))
     &                    +msq_struc(cs,-n,n)*(AP(a,g,2)+H1(a,g,q,cs,2)))*fx1z(g)/z*fx2(k)
     &                    +(msq_struc(cs,g,n)*(AP(q,q,2)+AP(q,q,3)
     &                    +H2(q,q,g,cs,2)+H2(q,q,g,cs,3))
     &                    +msq_struc(cs,g,g)*(AP(g,q,2)+H2(g,q,g,cs,2)))*fx1(g)*fx2z(k)/z
                  enddo
                  do cs=1,3
                     msq_aq=msq_struc(cs,-1,k)+msq_struc(cs,-2,k)+msq_struc(cs,-3,k)
     &                    +msq_struc(cs,-4,k)+msq_struc(cs,-5,k)
                     msq_qq=msq_struc(cs,+1,k)+msq_struc(cs,+2,k)+msq_struc(cs,+3,k)
     &                    +msq_struc(cs,+4,k)+msq_struc(cs,+5,k)
                     xmsqt=xmsqt     
     &                    +(+msq_aq*(AP(a,g,2)+H1(a,g,q,cs,2))
     &                    +msq_qq*(AP(q,g,2)+H1(q,g,q,cs,2)))*fx1z(g)/z*fx2(k)
     &                    +(+msq_struc(cs,g,g)*(AP(g,q,2)+H2(g,q,g,cs,2)))*fx1(g)*fx2z(k)/z
                  enddo
                  xmsqt=xmsqt     
     &                 +msq_struc(iqr,-n,n)*(AP(a,g,2)+H1(a,g,q,iqr,2))*fx1z(g)/z*fx2(k)
                  xmsq=xmsq+xmsqt
c     write(6,*) 'virtint: j,k,  ct=',j,k,xmsqt

               endif
               
            elseif (case .eq. 'twojet')then
               write(*,*) "fastNLO case twojet"
               
               xmsq=xmsq+fx1(j)*fx2(k)*(
     .              msqv(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k))

c---  SUM BY COLOUR STRUCTURES AND FINAL STATES: 2 jets only
               if      ((j .eq. 0) .and. (k .eq. 0)) then
                  tmp=xmsq
                  do cs=0,2
                     msq_qg=dfloat(nf)*msqx(cs,+1,g,+1,g)
                     msq_gq=dfloat(nf)*msqx(cs,g,+1,+1,g)
                     xmsq=xmsq
c     & +msqx(cs,g,g,g,g)*(AP(g,g,1)-AP(g,g,3)
c     &     +S1(g,g,g,gf_gf,cs,1)-S1(g,g,g,gf_gf,cs,3)
c     &     +AP(g,g,1)-AP(g,g,3)
c     &     +S2(g,g,g,gf_gf,cs,1)-S2(g,g,g,gf_gf,cs,3))*fx1(g)*fx2(g)
c     & +msqx(cs,g,g,g,g)*(AP(g,g,2)+AP(g,g,3)
c     &     +S1(g,g,g,gf_gf,cs,3)+S1(g,g,g,gf_gf,cs,2))*fx1z(g)/z*fx2(g)
c     & +msqx(cs,g,g,g,g)*(AP(g,g,2)+AP(g,g,3)
c     &     +S2(g,g,g,gf_gf,cs,3)+S2(g,g,g,gf_gf,cs,2))*fx1(g)*fx2z(g)/z
     &                    +msqx(cs,g,g,q,a)*dfloat(nf)*(AP(g,g,1)-AP(g,g,3)
     &                    +S1(g,g,g,qf_af,cs,1)-S1(g,g,g,qf_af,cs,3)
     &                    +AP(g,g,1)-AP(g,g,3)
     &                    +S2(g,g,g,qf_af,cs,1)-S2(g,g,g,qf_af,cs,3))*fx1(g)*fx2(g)
     &                    +msqx(cs,g,g,q,a)*dfloat(nf)*(AP(g,g,2)+AP(g,g,3)
     &                    +S1(g,g,g,qf_af,cs,3)+S1(g,g,g,qf_af,cs,2))*fx1z(g)/z*fx2(g)
     &                    +msqx(cs,g,g,q,a)*dfloat(nf)*(AP(g,g,2)+AP(g,g,3)
     &                    +S2(g,g,g,qf_af,cs,3)+S2(g,g,g,qf_af,cs,2))*fx1(g)*fx2z(g)/z
     &                    +msq_qg*(AP(q,g,2)+S1(q,g,g,qf_gf,cs,2)
     &                    +AP(a,g,2)+S1(a,g,g,af_gf,cs,2))*fx1z(g)/z*fx2(g)
     &                    +msq_gq*(AP(q,g,2)+S2(q,g,g,qf_gf,cs,2)
     &                    +AP(a,g,2)+S2(a,g,g,af_gf,cs,2))*fx1(g)*fx2z(g)/z
                  enddo
                  write(6,*) '_v: ',tmp
                  write(6,*) '_z: ',xmsq-tmp
               endif

            elseif ((case .eq. 'W_2jet') .or. (case .eq. 'Z_2jet')
     .              .or. (case .eq. 'W_bjet') .or. (case .eq. 'Z_bjet')
     .              .or. (case .eq. 'tt_tot') .or. (case .eq. 'bb_tot')
     .              .or. (case .eq. 'cc_tot')) then
c---  SUM BY COLOUR STRUCTURES: W/Z + 2 jet only
               write(*,*) "fastNLO case (8 cases)"

               xmsq=xmsq+fx1(j)*fx2(k)*(
     .              msqv(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k))
c     write(6,*) j,k,'-> msqv = ',fx1(j)*fx2(k)*(
c     . msqv(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k))
c     tmp=xmsq

               if ((j .gt. 0) .and. (k.gt.0)) then
                  do cs=0,2
c---  handle the Z+b jet case where identity of 5 and 6 are switched
c---  (interchange colour structures 1 and 2)
                     ics=cs
                     if ( (case .eq. 'Z_bjet') .and. (k .eq. +flav)
     .                    .and. (cs .gt. 0) ) ics=3-cs
                     xmsq=xmsq
     &                    +msq_cs(ics,j,k)*(AP(q,q,1)-AP(q,q,3)
     &                    +R1(q,q,q,cs,1)-R1(q,q,q,cs,3)
     &                    +AP(q,q,1)-AP(q,q,3)
     &                    +R2(q,q,q,cs,1)-R2(q,q,q,cs,3))*fx1(j)*fx2(k)
     &                    +(msq_cs(ics,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                    +R1(q,q,q,cs,2)+R1(q,q,q,cs,3))
     &                    + msq_cs(ics,g,k)*(AP(g,q,2)+R1(g,q,q,cs,2)))*fx1z(j)/z*fx2(k)
     &                    +(msq_cs(ics,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                    +R2(q,q,q,cs,2)+R2(q,q,q,cs,3))
     &                    + msq_cs(ics,j,g)*(AP(g,q,2)+R2(g,q,q,cs,2)))*fx1(j)*fx2z(k)/z
                  enddo
               elseif ((j .lt. 0) .and. (k.lt.0)) then
                  do cs=0,2
c---  handle the Z+b jet case where identity of 5 and 6 are switched
c---  (interchange colour structures 1 and 2)
                     ics=cs
                     if ( (case .eq. 'Z_bjet') .and. (k .eq. -flav)
     .                    .and. (cs .gt. 0) ) ics=3-cs
                     xmsq=xmsq
     &                    +msq_cs(ics,j,k)*(AP(a,a,1)-AP(a,a,3)
     &                    +R1(a,a,a,cs,1)-R1(a,a,a,cs,3)
     &                    +AP(a,a,1)-AP(a,a,3)
     &                    +R2(a,a,a,cs,1)-R2(a,a,a,cs,3))*fx1(j)*fx2(k)
     &                    +(msq_cs(ics,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                    +R1(a,a,a,cs,2)+R1(a,a,a,cs,3))
     &                    + msq_cs(ics,g,k)*(AP(g,a,2)+R1(g,a,a,cs,2)))*fx1z(j)/z*fx2(k)
     &                    +(msq_cs(ics,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                    +R2(a,a,a,cs,2)+R2(a,a,a,cs,3))
     &                    + msq_cs(ics,j,g)*(AP(g,a,2)+R2(g,a,a,cs,2)))*fx1(j)*fx2z(k)/z
                  enddo
               elseif ((j .gt. 0) .and. (k.lt.0)) then
c---  handle the Z+b jet case where identity of 5 and 6 are switched
c---  (interchange integrated CT's for q-qbar and qbar-q)
                  iq=q
                  ia=a
                  if ((case .eq. 'Z_bjet') .and. (k .eq. -flav)) then
                     iq=a
                     ia=q
                  endif
                  do cs=0,2
                     xmsq=xmsq
     &                    +msq_cs(cs,j,k)*(AP(q,q,1)-AP(q,q,3)
     &                    +R1(iq,iq,ia,cs,1)-R1(iq,iq,ia,cs,3)
     &                    +AP(a,a,1)-AP(a,a,3)
     &                    +R2(ia,ia,iq,cs,1)-R2(ia,ia,iq,cs,3)
     &                    )*fx1(j)*fx2(k)
     &                    +(msq_cs(cs,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                    +R1(iq,iq,ia,cs,2)+R1(iq,iq,ia,cs,3))
     &                    + msq_cs(cs,g,k)*(AP(g,q,2)+R1(g,q,a,cs,2)))*fx1z(j)/z*fx2(k)
     &                    +(msq_cs(cs,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                    +R2(ia,ia,iq,cs,2)+R2(ia,ia,iq,cs,3))
     &                    + msq_cs(cs,j,g)*(AP(g,a,2)+R2(g,a,q,cs,2)))*fx1(j)*fx2z(k)/z
                  enddo
               elseif ((j .lt. 0) .and. (k.gt.0)) then
c---  handle the Z+b jet case where identity of 5 and 6 are switched
c---  (interchange integrated CT's for q-qbar and qbar-q)
                  iq=q
                  ia=a
                  if ((case .eq. 'Z_bjet') .and. (j .eq. -flav)) then
                     iq=a
                     ia=q
                  endif
                  do cs=0,2
                     xmsq=xmsq
     &                    +msq_cs(cs,j,k)*(AP(a,a,1)-AP(a,a,3)
     &                    +R1(ia,ia,iq,cs,1)-R1(ia,ia,iq,cs,3)
     &                    +AP(q,q,1)-AP(q,q,3)
     &                    +R2(iq,iq,ia,cs,1)-R2(iq,iq,ia,cs,3)
     &                    )*fx1(j)*fx2(k)
     &                    +(msq_cs(cs,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                    +R1(ia,ia,iq,cs,2)+R1(ia,ia,iq,cs,3))
     &                    + msq_cs(cs,g,k)*(AP(g,a,2)+R1(g,a,q,cs,2)))*fx1z(j)/z*fx2(k)
     &                    +(msq_cs(cs,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                    +R2(iq,iq,ia,cs,2)+R2(iq,iq,ia,cs,3))
     &                    + msq_cs(cs,j,g)*(AP(g,q,2)+R2(g,q,a,cs,2)))*fx1(j)*fx2z(k)/z
                  enddo
               elseif ((j .eq. g) .and. (k .eq. g)) then
                  do cs=0,2
                     msq_qg=msq_cs(cs,+5,g)+msq_cs(cs,+4,g)+msq_cs(cs,+3,g)
     &                    +msq_cs(cs,+2,g)+msq_cs(cs,+1,g)
     &                    +msq_cs(cs,-5,g)+msq_cs(cs,-4,g)+msq_cs(cs,-3,g)
     &                    +msq_cs(cs,-2,g)+msq_cs(cs,-1,g)
                     msq_gq=msq_cs(cs,g,+5)+msq_cs(cs,g,+4)+msq_cs(cs,g,+3)
     &                    +msq_cs(cs,g,+2)+msq_cs(cs,g,+1)
     &                    +msq_cs(cs,g,-5)+msq_cs(cs,g,-4)+msq_cs(cs,g,-3)
     &                    +msq_cs(cs,g,-2)+msq_cs(cs,g,-1)
                     xmsq=xmsq
     &                    +msq_cs(cs,g,g)*(AP(g,g,1)-AP(g,g,3)
     &                    +R1(g,g,g,cs,1)-R1(g,g,g,cs,3)
     &                    +AP(g,g,1)-AP(g,g,3)
     &                    +R2(g,g,g,cs,1)-R2(g,g,g,cs,3))*fx1(g)*fx2(g)
     &                    +(msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                    +R1(g,g,g,cs,3)+R1(g,g,g,cs,2))
     &                    + msq_qg*(AP(q,g,2)+R1(q,g,g,cs,2)))*fx1z(g)/z*fx2(g)
     &                    +(msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                    +R2(g,g,g,cs,3)+R2(g,g,g,cs,2))
     &                    + msq_gq*(AP(q,g,2)+R2(q,g,g,cs,2)))*fx1(g)*fx2z(g)/z
                  enddo
               elseif ((j .eq. g) .and. (k .gt. 0)) then
                  do cs=0,2
                     msq_aq=msq_cs(cs,-1,k)+msq_cs(cs,-2,k)+msq_cs(cs,-3,k)
     &                    +msq_cs(cs,-4,k)+msq_cs(cs,-5,k)
                     msq_qq=msq_cs(cs,+1,k)+msq_cs(cs,+2,k)+msq_cs(cs,+3,k)
     &                    +msq_cs(cs,+4,k)+msq_cs(cs,+5,k)
                     xmsq=xmsq     
     &                    +msq_cs(cs,g,k)*(AP(g,g,1)-AP(g,g,3)
     &                    +R1(g,g,q,cs,1)-R1(g,g,q,cs,3)
     &                    +AP(q,q,1)-AP(q,q,3)
     &                    +R2(q,q,g,cs,1)-R2(q,q,g,cs,3))*fx1(g)*fx2(k)
     &                    +(msq_cs(cs,g,k)*(AP(g,g,2)+AP(g,g,3)
     &                    +R1(g,g,q,cs,2)+R1(g,g,q,cs,3))
     &                    + msq_aq*(AP(a,g,2)+R1(a,g,q,cs,2))
     &                    + msq_qq*(AP(q,g,2)+R1(q,g,q,cs,2)))*fx1z(g)/z*fx2(k)
     &                    +(msq_cs(cs,g,k)*(AP(q,q,2)+AP(q,q,3)
     &                    +R2(q,q,g,cs,2)+R2(q,q,g,cs,3))
     &                    + msq_cs(cs,g,g)*(AP(g,q,2)+R2(g,q,g,cs,2)))*fx1(g)*fx2z(k)/z
                  enddo

               elseif ((j .eq. g) .and. (k .lt. 0)) then
                  do cs=0,2
                     msq_qa=msq_cs(cs,+1,k)+msq_cs(cs,+2,k)+msq_cs(cs,+3,k)
     &                    +msq_cs(cs,+4,k)+msq_cs(cs,+5,k)
                     msq_aa=msq_cs(cs,-1,k)+msq_cs(cs,-2,k)+msq_cs(cs,-3,k)
     &                    +msq_cs(cs,-4,k)+msq_cs(cs,-5,k)
                     xmsq=xmsq     
     &                    +msq_cs(cs,g,k)*(AP(g,g,1)-AP(g,g,3)
     &                    +R1(g,g,a,cs,1)-R1(g,g,a,cs,3)
     &                    +AP(a,a,1)-AP(a,a,3)
     &                    +R2(a,a,g,cs,1)-R2(a,a,g,cs,3))*fx1(g)*fx2(k)
     &                    +(msq_cs(cs,g,k)*(AP(g,g,2)+AP(g,g,3)
     &                    +R1(g,g,a,cs,2)+R1(g,g,a,cs,3))
     &                    + msq_qa*(AP(q,g,2)+R1(q,g,a,cs,2))
     &                    + msq_aa*(AP(a,g,2)+R1(a,g,a,cs,2)))*fx1z(g)/z*fx2(k)
     &                    +(msq_cs(cs,g,k)*(AP(a,a,2)+AP(a,a,3)
     &                    +R2(a,a,g,cs,2)+R2(a,a,g,cs,3))
     &                    + msq_cs(cs,g,g)*(AP(g,a,2)+R2(g,a,g,cs,2)))*fx1(g)*fx2z(k)/z
                  enddo

               elseif ((j .gt. 0) .and. (k .eq. g)) then
                  do cs=0,2
                     msq_qa=msq_cs(cs,j,-1)+msq_cs(cs,j,-2)+msq_cs(cs,j,-3)
     &                    +msq_cs(cs,j,-4)+msq_cs(cs,j,-5)
                     msq_qq=msq_cs(cs,j,+1)+msq_cs(cs,j,+2)+msq_cs(cs,j,+3)
     &                    +msq_cs(cs,j,+4)+msq_cs(cs,j,+5)
                     xmsq=xmsq     
     &                    + msq_cs(cs,j,g)*(AP(q,q,1)-AP(q,q,3)
     &                    +R1(q,q,g,cs,1)-R1(q,q,g,cs,3)
     &                    +R2(g,g,q,cs,1)-R2(g,g,q,cs,3)
     &                    +AP(g,g,1)-AP(g,g,3))*fx1(j)*fx2(g)
     &                    +(msq_cs(cs,j,g)*(AP(q,q,2)+AP(q,q,3)
     &                    +R1(q,q,g,cs,2)+R1(q,q,g,cs,3))
     &                    + msq_cs(cs,g,g)*(AP(g,q,2)+R1(g,q,g,cs,2)))*fx1z(j)/z*fx2(g)
     &                    +(msq_cs(cs,j,g)*(AP(g,g,2)+AP(g,g,3)
     &                    +R2(g,g,q,cs,2)+R2(g,g,q,cs,3))
     &                    + msq_qa*(AP(a,g,2)+R2(a,g,q,cs,2))
     &                    + msq_qq*(AP(q,g,2)+R2(q,g,q,cs,2)))*fx1(j)*fx2z(g)/z
                  enddo
                  
               elseif ((j .lt. 0) .and. (k .eq. g)) then
                  do cs=0,2
                     msq_aq=msq_cs(cs,j,+1)+msq_cs(cs,j,+2)+msq_cs(cs,j,+3)
     &                    +msq_cs(cs,j,+4)+msq_cs(cs,j,+5)
                     msq_aa=msq_cs(cs,j,-1)+msq_cs(cs,j,-2)+msq_cs(cs,j,-3)
     &                    +msq_cs(cs,j,-4)+msq_cs(cs,j,-5)
                     xmsq=xmsq     
     &                    + msq_cs(cs,j,g)*(AP(a,a,1)-AP(a,a,3)
     &                    +R1(a,a,g,cs,1)-R1(a,a,g,cs,3)
     &                    +AP(g,g,1)-AP(g,g,3)
     &                    +R2(g,g,a,cs,1)-R2(g,g,a,cs,3))*fx1(j)*fx2(g)
     &                    +(msq_cs(cs,j,g)*(AP(a,a,2)+AP(a,a,3)
     &                    +R1(a,a,g,cs,2)+R1(a,a,g,cs,3))
     &                    + msq_cs(cs,g,g)*(AP(g,a,2)+R1(g,a,g,cs,2)))*fx1z(j)/z*fx2(g)
     &                    +(msq_cs(cs,j,g)*(AP(g,g,2)+AP(g,g,3)
     &                    +R2(g,g,a,cs,2)+R2(g,g,a,cs,3))
     &                    + msq_aq*(AP(q,g,2)+R2(q,g,a,cs,2))
     &                    + msq_aa*(AP(a,g,2)+R2(a,g,a,cs,2)))*fx1(j)*fx2z(g)/z
                  enddo
               endif

c     write(6,*) j,k,'-> msqc = ',xmsq-tmp

            else

c---  SUM BY TOTAL MATRIX ELEMENTS: everything else
c---  special code to remove the Q+Qb -> Z+Q+Qb contribution
               if ((j*k .eq. -flav*flav) .and. (case .eq. 'gQ__ZQ')) then
                  xmsq=xmsq-(
     &                 + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2))*fx1z(j)/z*fx2(k)
     &                 + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2))*fx1(j)*fx2z(k)/z)
               endif
c---  special code to remove the b+b(bar) -> W+t+b(bar) contribution
               if ( (abs(j*k) .eq. nf*nf) .and. 
     .              ((case .eq. 'W_tndk') .or. (case .eq. 'W_twdk')) ) then
                  xmsq=xmsq-(
     &                 + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2))*fx1z(j)/z*fx2(k)
     &                 + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2))*fx1(j)*fx2z(k)/z)
               endif
c === fastNLO comment: the following 7 choices are relevant in W
C--   QQ
               if     ((j .gt. 0) .and. (k.gt.0)) then
c                  write(*,*) "fn here 1 QQ"
                  if ((case .eq. 'bq_tpq') .and. (j .eq. 5)) then
                     write(*,*)'fastNLO a'
                     stop
                     xmsq=xmsq+(msqv(j,k)
     &                    + msq(j,k)*(one+AP(q,q,1)-AP(q,q,3)+B1(b,b,q,1)-B1(b,b,q,3)
     &                    +AP(q,q,1)-AP(q,q,3)+B2(q,q,b,1)-B2(q,q,b,3)))
     &                    *fx1(j)*fx2(k)
     &                    +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B1(b,b,q,2)+B1(b,b,q,3))
     &                    + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))*fx1z(j)/z*fx2(k)
     &                    +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B2(q,q,b,2)+B2(q,q,b,3))
     &                    + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))*fx1(j)*fx2z(k)/z
                     
                  elseif ((case .eq. 'bq_tpq') .and. (k .eq. 5)) then
                     write(*,*)'fastNLO a'
                     stop
                     xmsq=xmsq+(msqv(j,k)
     &                    + msq(j,k)*(one+AP(q,q,1)-AP(q,q,3)+B1(q,q,b,1)-B1(q,q,b,3)
     &                    +AP(q,q,1)-AP(q,q,3)+B2(b,b,q,1)-B2(b,b,q,3)))
     &                    *fx1(j)*fx2(k)
     &                    +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B1(q,q,b,2)+B1(q,q,b,3))
     &                    + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))*fx1z(j)/z*fx2(k)
     &                    +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+B2(b,b,q,2)+B2(b,b,q,3))
     &                    + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))*fx1(j)*fx2z(k)/z
                  else
                     xmsq=xmsq+(msqv(j,k)
     &                    + msq(j,k)*(one+AP(q,q,1)-AP(q,q,3)+Q1(q,q,q,1)-Q1(q,q,q,3)
     &                    +AP(q,q,1)-AP(q,q,3)+Q2(q,q,q,1)-Q2(q,q,q,3)))
     &                    *fx1(j)*fx2(k)
     &                    +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,q,2)+Q1(q,q,q,3))
     &                    + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))*fx1z(j)/z*fx2(k)
     &                    +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,q,2)+Q2(q,q,q,3))
     &                    + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))*fx1(j)*fx2z(k)/z
c === fastNLO
                     fncoeff(j,k) = fncoeff(j,k) + (msqv(j,k)
     +                    + msq(j,k)*(one+AP(q,q,1)-AP(q,q,3)+Q1(q,q,q,1)-Q1(q,q,q,3)
     +                    +AP(q,q,1)-AP(q,q,3)+Q2(q,q,q,1)-Q2(q,q,q,3)))
                     fncfv2(j,k) = fncfv2(j,k) +
     +                    (msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,q,2)+Q1(q,q,q,3))
     +                    + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))/z
                     fncfv3(j,k) = fncfv3(j,k) +
     +                    (msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,q,2)+Q2(q,q,q,3))
     +                    + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))/z
c ======================================================================
                  endif
C--   QbarQbar
               elseif ((j .lt. 0) .and. (k.lt.0)) then
c                  write(*,*) "fn here 2 QbQb"
                  if ((case .eq. 'bq_tpq') .and. (j .eq. -5)) then
                     write(*,*)'fastNLO a'
                     stop
                     xmsq=xmsq+(msqv(j,k)
     &                    + msq(j,k)*(one+AP(a,a,1)-AP(a,a,3)+B1(b,b,a,1)-B1(b,b,a,3)
     &                    +AP(a,a,1)-AP(a,a,3)+B2(a,a,b,1)-B2(a,a,b,3)))
     &                    *fx1(j)*fx2(k)
     &                    +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+B1(b,b,a,2)+B1(b,b,a,3))
     &                    + msq(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))*fx1z(j)/z*fx2(k)
     &                    +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+B2(a,a,b,2)+B2(a,a,b,3))
     &                    + msq(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))*fx1(j)*fx2z(k)/z
                  elseif ((case .eq. 'bq_tpq') .and. (k .eq. -5)) then
                     write(*,*)'fastNLO a'
                     stop
                     xmsq=xmsq+(msqv(j,k)
     &                    + msq(j,k)*(one+AP(a,a,1)-AP(a,a,3)+B1(a,a,b,1)-B1(a,a,b,3)
     &                    +AP(a,a,1)-AP(a,a,3)+B2(b,b,a,1)-B2(b,b,a,3)))
     &                    *fx1(j)*fx2(k)
     &                    +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+B1(a,a,b,2)+B1(a,a,b,3))
     &                    + msq(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))*fx1z(j)/z*fx2(k)
     &                    +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+B2(b,b,a,2)+B2(b,b,a,3))
     &                    + msq(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))*fx1(j)*fx2z(k)/z
                  else
                     xmsq=xmsq+(msqv(j,k)
     &                    + msq(j,k)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,a,1)-Q1(a,a,a,3)
     &                    +AP(a,a,1)-AP(a,a,3)+Q2(a,a,a,1)-Q2(a,a,a,3)))
     &                    *fx1(j)*fx2(k)
     &                    +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,a,2)+Q1(a,a,a,3))
     &                    + msq(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))*fx1z(j)/z*fx2(k)
     &                    +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,a,2)+Q2(a,a,a,3))
     &                    + msq(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))*fx1(j)*fx2z(k)/z
c === fastNLO
                     fncoeff(j,k) = fncoeff(j,k) + (msqv(j,k)
     +                    + msq(j,k)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,a,1)-Q1(a,a,a,3)
     +                    +AP(a,a,1)-AP(a,a,3)+Q2(a,a,a,1)-Q2(a,a,a,3)))
                     fncfv2(j,k) = fncfv2(j,k) +
     +                    (msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,a,2)+Q1(a,a,a,3))
     +                    + msq(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))/z
                     fncfv3(j,k) = fncfv3(j,k) +
     +                    (msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,a,2)+Q2(a,a,a,3))
     +                    + msq(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))/z
c ======================================================================
                  endif
C--   QQbar
               elseif ((j .gt. 0) .and. (k.lt.0)) then
c                  write(*,*) "fn here 3 QQb"
                  xmsq=xmsq+(msqv(j,k)
     &                 + msq(j,k)*(one+AP(q,q,1)-AP(q,q,3)+Q1(q,q,a,1)-Q1(q,q,a,3)
     &                 +AP(a,a,1)-AP(a,a,3)+Q2(a,a,q,1)-Q2(a,a,q,3)))
     &                 *fx1(j)*fx2(k)
     &                 +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,a,3)+Q1(q,q,a,2))
     &                 + msq(g,k)*(AP(g,q,2)+Q1(g,q,a,2)))*fx1z(j)/z*fx2(k)
     &                 +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,q,3)+Q2(a,a,q,2))
     &                 + msq(j,g)*(AP(g,a,2)+Q2(g,a,q,2)))*fx1(j)*fx2z(k)/z
c === fastNLO
                  fncoeff(j,k) = fncoeff(j,k) + (msqv(j,k)
     +                 + msq(j,k)*(one+AP(q,q,1)-AP(q,q,3)+Q1(q,q,a,1)-Q1(q,q,a,3)
     +                 +AP(a,a,1)-AP(a,a,3)+Q2(a,a,q,1)-Q2(a,a,q,3)))
                  fncfv2(j,k) = fncfv2(j,k) +
     +                 (msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,a,3)+Q1(q,q,a,2))
     +                 + msq(g,k)*(AP(g,q,2)+Q1(g,q,a,2)))/z
                  fncfv3(j,k) = fncfv3(j,k) +
     +                 (msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,q,3)+Q2(a,a,q,2))
     +                 + msq(j,g)*(AP(g,a,2)+Q2(g,a,q,2)))/z
c ======================================================================

               elseif ((j .lt. 0) .and. (k.gt.0)) then
c                  write(*,*) "fn here 4 QbQ"
C--   QbarQ
                  xmsq=xmsq+(msqv(j,k)
     &                 +msq(j,k)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,q,1)-Q1(a,a,q,3)
     &                 +AP(q,q,1)-AP(q,q,3)+Q2(q,q,a,1)-Q2(q,q,a,3)))
     &                 *fx1(j)*fx2(k)
     &                 +(msq(j,k)*(AP(a,a,3)+AP(a,a,2)+Q1(a,a,q,3)+Q1(a,a,q,2))
     &                 + msq(g,k)*(AP(g,a,2)+Q1(g,a,q,2)))*fx1z(j)/z*fx2(k)
     &                 +(msq(j,k)*(AP(q,q,3)+AP(q,q,2)+Q2(q,q,a,3)+Q2(q,q,a,2))
     &                 + msq(j,g)*(AP(g,q,2)+Q2(g,q,a,2)))*fx1(j)*fx2z(k)/z
c === fastNLO
                  fncoeff(j,k) = fncoeff(j,k) + (msqv(j,k)
     +                 +msq(j,k)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,q,1)-Q1(a,a,q,3)
     +                 +AP(q,q,1)-AP(q,q,3)+Q2(q,q,a,1)-Q2(q,q,a,3)))
                  fncfv2(j,k) = fncfv2(j,k) +
     +                 (msq(j,k)*(AP(a,a,3)+AP(a,a,2)+Q1(a,a,q,3)+Q1(a,a,q,2))
     +                 + msq(g,k)*(AP(g,a,2)+Q1(g,a,q,2)))/z
                  fncfv3(j,k) = fncfv3(j,k) +
     +                 (msq(j,k)*(AP(q,q,3)+AP(q,q,2)+Q2(q,q,a,3)+Q2(q,q,a,2))
     +                 + msq(j,g)*(AP(g,q,2)+Q2(g,q,a,2)))/z
c ======================================================================

               elseif ((j .eq. g) .and. (k.eq.g)) then
c                 write (*,*) "fn here 5 gg"
C--   gg
                  if (case .eq. 'Zbbbar') then
                     write(*,*)'fastNLO a'
                     stop
                     xmsq=xmsq+msqv(g,g)*fx1(g)*fx2(g)
                     do cs=0,2
                        xmsq=xmsq+(
     &                       +msq_cs(cs,g,g)*(one
     &                       +AP(g,g,1)-AP(g,g,3)+R1(g,g,g,cs,1)-R1(g,g,g,cs,3)
     &                       +AP(g,g,1)-AP(g,g,3)+R2(g,g,g,cs,1)-R2(g,g,g,cs,3))
     &                       )*fx1(g)*fx2(g)
     &                       +msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                       +R1(g,g,g,cs,2)+R1(g,g,g,cs,3))*fx1z(g)/z*fx2(g)
     &                       +msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                       +R2(g,g,g,cs,2)+R2(g,g,g,cs,3))*fx1(g)*fx2z(g)/z
                     enddo
                  else
                     msq_qg=msq(+5,g)+msq(+4,g)+msq(+3,g)+msq(+2,g)+msq(+1,g)
     &                    +msq(-5,g)+msq(-4,g)+msq(-3,g)+msq(-2,g)+msq(-1,g)
                     msq_gq=msq(g,+5)+msq(g,+4)+msq(g,+3)+msq(g,+2)+msq(g,+1)
     &                    +msq(g,-5)+msq(g,-4)+msq(g,-3)+msq(g,-2)+msq(g,-1)
                     xmsq=xmsq+(msqv(g,g)
     &                    +msq(g,g)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,g,1)-Q1(g,g,g,3)
     &                    +AP(g,g,1)-AP(g,g,3)+Q2(g,g,g,1)-Q2(g,g,g,3)))
     &                    *fx1(g)*fx2(g)
     &                    +(msq(g,g)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,g,2)+Q1(g,g,g,3))
     &                    +   msq_qg*(AP(q,g,2)+Q1(q,g,g,2)))*fx1z(g)/z*fx2(g)
     &                    +(msq(g,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,g,2)+Q2(g,g,g,3))
     &                    +   msq_gq*(AP(q,g,2)+Q2(q,g,g,2)))*fx1(g)*fx2z(g)/z
c === fastNLO
                     fncoeff(j,k) = fncoeff(j,k) + (msqv(g,g)
     +                    +msq(g,g)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,g,1)-Q1(g,g,g,3)
     +                    +AP(g,g,1)-AP(g,g,3)+Q2(g,g,g,1)-Q2(g,g,g,3)))  
                     fncfv2(j,k) = fncfv2(j,k) +
     +                    (msq(g,g)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,g,2)+Q1(g,g,g,3))
     +                    +   msq_qg*(AP(q,g,2)+Q1(q,g,g,2)))/z
                     fncfv3(j,k) = fncfv3(j,k) +
     +                    (msq(g,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,g,2)+Q2(g,g,g,3))
     +                    +   msq_gq*(AP(q,g,2)+Q2(q,g,g,2)))
c ======================================================================
                  endif

               elseif (j .eq. g) then
c                  write(*,*) "fn here 6 gQ gQb"
C--   gQ
                  if    (k .gt. 0) then
                     msq_aq=msq(-1,k)+msq(-2,k)+msq(-3,k)+msq(-4,k)+msq(-5,k)
                     msq_qq=msq(+1,k)+msq(+2,k)+msq(+3,k)+msq(+4,k)+msq(+5,k)
                     xmsq=xmsq+(msqv(g,k)
     &                    +msq(g,k)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,q,1)-Q1(g,g,q,3)
     &                    +AP(q,q,1)-AP(q,q,3)+Q2(q,q,g,1)-Q2(q,q,g,3)))
     &                    *fx1(g)*fx2(k)
     &                    +(msq(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,q,2)+Q1(g,g,q,3))
     &                    +   msq_aq*(AP(a,g,2)+Q1(a,g,q,2))
     &                    +   msq_qq*(AP(q,g,2)+Q1(q,g,q,2)))*fx1z(g)/z*fx2(k)
     &                    +(msq(g,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,g,2)+Q2(q,q,g,3))
     &                    + msq(g,g)*(AP(g,q,2)+Q2(g,q,g,2)))*fx1(g)*fx2z(k)/z
c === fastNLO
                     fncoeff(j,k) = fncoeff(j,k) + (msqv(g,k)
     +                    +msq(g,k)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,q,1)-Q1(g,g,q,3)
     +                    +AP(q,q,1)-AP(q,q,3)+Q2(q,q,g,1)-Q2(q,q,g,3))) 
                     fncfv2(j,k) = fncfv2(j,k) +
     +                    (msq(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,q,2)+Q1(g,g,q,3))
     +                    +   msq_aq*(AP(a,g,2)+Q1(a,g,q,2))
     +                    +   msq_qq*(AP(q,g,2)+Q1(q,g,q,2)))/z
                     fncfv3(j,k) = fncfv3(j,k) +
     +                    (msq(g,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,g,2)+Q2(q,q,g,3))
     +                    + msq(g,g)*(AP(g,q,2)+Q2(g,q,g,2)))/z
c ======================================================================

C--   gQbar

                  elseif (k.lt.0) then
                     msq_qa=msq(+1,k)+msq(+2,k)+msq(+3,k)+msq(+4,k)+msq(+5,k)
                     msq_aa=msq(-1,k)+msq(-2,k)+msq(-3,k)+msq(-4,k)+msq(-5,k)
                     xmsq=xmsq+(msqv(g,k)
     &                    +msq(g,k)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,a,1)-Q1(g,g,a,3)
     &                    +AP(a,a,1)-AP(a,a,3)+Q2(a,a,g,1)-Q2(a,a,g,3)))
     &                    *fx1(g)*fx2(k)
     &                    +(msq(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,a,2)+Q1(g,g,a,3))
     &                    +   msq_qa*(AP(q,g,2)+Q1(q,g,a,2))
     &                    +   msq_aa*(AP(a,g,2)+Q1(a,g,a,2)))*fx1z(g)/z*fx2(k)
     &                    +(msq(g,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,g,2)+Q2(a,a,g,3))
     &                    + msq(g,g)*(AP(g,a,2)+Q2(g,a,g,2)))*fx1(g)*fx2z(k)/z
c === fastNLO
                     fncoeff(j,k) = fncoeff(j,k) + (msqv(g,k)
     +                    +msq(g,k)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,a,1)-Q1(g,g,a,3)
     +                    +AP(a,a,1)-AP(a,a,3)+Q2(a,a,g,1)-Q2(a,a,g,3)))
                     fncfv2(j,k) = fncfv2(j,k) +
     +                    (msq(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,a,2)+Q1(g,g,a,3))
     +                    +   msq_qa*(AP(q,g,2)+Q1(q,g,a,2))
     +                    +   msq_aa*(AP(a,g,2)+Q1(a,g,a,2)))/z
                     fncfv3(j,k) = fncfv3(j,k) +
     +                    (msq(g,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,g,2)+Q2(a,a,g,3))
     +                    + msq(g,g)*(AP(g,a,2)+Q2(g,a,g,2)))/z
c ======================================================================
                  endif
C--   Qg
               elseif (k .eq. g) then
c                  write(*,*) "fn here 7 Qg Qbg"
                  if     (j.gt.0) then
                     msq_qa=msq(j,-1)+msq(j,-2)+msq(j,-3)+msq(j,-4)+msq(j,-5)
                     msq_qq=msq(j,+1)+msq(j,+2)+msq(j,+3)+msq(j,+4)+msq(j,+5)
                     xmsq=xmsq+(msqv(j,g)
     &                    +msq(j,g)*(one
     &                    +AP(q,q,1)-AP(q,q,3)+Q1(q,q,g,1)-Q1(q,q,g,3)
     &                    +AP(g,g,1)-AP(g,g,3)+Q2(g,g,q,1)-Q2(g,g,q,3)))
     &                    *fx1(j)*fx2(g)
     &                    +(msq(j,g)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,g,2)+Q1(q,q,g,3))
     &                    + msq(g,g)*(AP(g,q,2)+Q1(g,q,g,2)))*fx1z(j)/z*fx2(g)
     &                    +(msq(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,q,2)+Q2(g,g,q,3))
     &                    +   msq_qa*(AP(a,g,2)+Q2(a,g,q,2))
     &                    +   msq_qq*(AP(q,g,2)+Q2(q,g,q,2)))*fx1(j)*fx2z(g)/z
c === fastNLO
                     fncoeff(j,k) = fncoeff(j,k) + (msqv(j,g)
     +                    +msq(j,g)*(one
     +                    +AP(q,q,1)-AP(q,q,3)+Q1(q,q,g,1)-Q1(q,q,g,3)
     +                    +AP(g,g,1)-AP(g,g,3)+Q2(g,g,q,1)-Q2(g,g,q,3)))
                     fncfv2(j,k) = fncfv2(j,k) +
     +                    (msq(j,g)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,g,2)+Q1(q,q,g,3))
     +                    + msq(g,g)*(AP(g,q,2)+Q1(g,q,g,2)))/z
                     fncfv3(j,k) = fncfv3(j,k) +
     +                    (msq(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,q,2)+Q2(g,g,q,3))
     +                    +   msq_qa*(AP(a,g,2)+Q2(a,g,q,2))
     +                    +   msq_qq*(AP(q,g,2)+Q2(q,g,q,2)))/z
c ======================================================================
C--   Qbarg
                  elseif (j.lt.0) then
                     msq_aq=msq(j,+1)+msq(j,+2)+msq(j,+3)+msq(j,+4)+msq(j,+5)
                     msq_aa=msq(j,-1)+msq(j,-2)+msq(j,-3)+msq(j,-4)+msq(j,-5)
                     xmsq=xmsq+(msqv(j,g)
     &                    +msq(j,g)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,g,1)-Q1(a,a,g,3)
     &                    +AP(g,g,1)-AP(g,g,3)+Q2(g,g,a,1)-Q2(g,g,a,3)))
     &                    *fx1(j)*fx2(g)
     &                    +(msq(j,g)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,g,2)+Q1(a,a,g,3))
     &                    + msq(g,g)*(AP(g,a,2)+Q1(g,a,g,2)))*fx1z(j)/z*fx2(g)
     &                    +(msq(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,a,3)+Q2(g,g,a,2))
     &                    + msq_aq*(AP(q,g,2)+Q2(q,g,a,2))
     &                    + msq_aa*(AP(a,g,2)+Q2(a,g,a,2)))*fx1(j)*fx2z(g)/z
c === fastNLO
                     fncoeff(j,k) = fncoeff(j,k) + (msqv(j,g)
     +                    +msq(j,g)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,g,1)-Q1(a,a,g,3)
     +                    +AP(g,g,1)-AP(g,g,3)+Q2(g,g,a,1)-Q2(g,g,a,3)))
                     fncfv2(j,k) = fncfv2(j,k) +
     +                    (msq(j,g)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,g,2)+Q1(a,a,g,3))
     +                    + msq(g,g)*(AP(g,a,2)+Q1(g,a,g,2)))/z
                     fncfv3(j,k) = fncfv3(j,k) +
     +                    (msq(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,a,3)+Q2(g,g,a,2))
     +                    + msq_aq*(AP(q,g,2)+Q2(q,g,a,2))
     +                    + msq_aa*(AP(a,g,2)+Q2(a,g,a,2)))/z
c ======================================================================
                  endif

               endif      
 100           continue ! >>>>>>>>>>>>>>>>>>>>>>. fastNLO test
            endif
c            write(*,*) 'fastNLO coeff',fncoeff
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
               xmsq_bypart(sgnj,sgnk)=xmsq_bypart(sgnj,sgnk)+(xmsq-tmp)
            endif
            
 20         continue

         enddo
      enddo

      if (currentPDF .eq. 0) then
         virtint=flux*xjac*pswt*xmsq/BrnRat
      endif
      
c---  loop over all PDF error sets, if necessary
      if (PDFerrors) then
         PDFwgt(currentPDF)=flux*xjac*pswt*xmsq/BrnRat*wgt/itmx
         PDFxsec(currentPDF)=PDFxsec(currentPDF)
     .        +PDFwgt(currentPDF)
         currentPDF=currentPDF+1
         if (currentPDF .le. maxPDFsets) goto 777
      endif    

c---  code to check that epsilon poles cancel      
      if (nshot .eq. 1) then
         if (xmsq .eq. 0d0) goto 999
         xmsq_old=xmsq
         nshot=nshot+1
         epinv=0d0
         epinv2=0d0
c     epinv=1d0
c     epinv2=1d0
         goto 12
      elseif (nshot .eq. 2) then
         nshot=nshot+1

         if (abs(xmsq_old/xmsq-1d0) .gt. 1d-6) then
            write(6,*) 'epsilon fails to cancel'
            write(6,*) 'xmsq (epinv=large) = ',xmsq_old
            write(6,*) 'xmsq (epinv=zero ) = ',xmsq
C     stop
            colourchoice=rvcolourchoice
         else
            write(6,*) 'Poles cancelled!'
c     write(6,*) 'xmsq (epinv=large) = ',xmsq_old
c     write(6,*) 'xmsq (epinv=zero ) = ',xmsq
c     pause
            colourchoice=rvcolourchoice
         endif
      endif
      
      if (creatent) then
         wt_gg=xmsq_bypart(0,0)*wgt*flux*xjac*pswt/BrnRat/dfloat(itmx)
         wt_gq=(xmsq_bypart(+1,0)+xmsq_bypart(-1,0)
     .        +xmsq_bypart(0,+1)+xmsq_bypart(0,-1)
     .        )*wgt*flux*xjac*pswt/BrnRat/dfloat(itmx)
         wt_qq=(xmsq_bypart(+1,+1)+xmsq_bypart(-1,-1)
     .        )*wgt*flux*xjac*pswt/BrnRat/dfloat(itmx)
         wt_qqb=(xmsq_bypart(+1,-1)+xmsq_bypart(-1,+1)
     .        )*wgt*flux*xjac*pswt/BrnRat/dfloat(itmx)
      endif

      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)

      do j=-1,1
         do k=-1,1
            lord_bypart(j,k)=lord_bypart(j,k)+
     .           wgt*flux*xjac*pswt*xmsq_bypart(j,k)/BrnRat
         enddo
      enddo

c     === fastNLO
      fnsigma = flux*xjac*pswt/BrnRat*wgt/dfloat(itmx)
      ifnctrb = 2
c     ===================================================

      val=virtint*wgt 
      val2=virtint**2*wgt 
c---  update the maximum weight so far, if necessary
      if (val .gt. wtmax) then
         wtmax=val
      endif

      if (bin) then
         val=val/dfloat(itmx) 
         val2=val2/dfloat(itmx)**2
         call nplotter(pjet,val,val2,0)
      endif

      return

 999  continue
      ntotzero=ntotzero+1
      
      return
      end


