      PROGRAM THRESHCOR
c ---------------------------------------------------------------------
c --- M. Wobisch 2005/11/17
c --- 
c --- Computation of threshold corrections in NNLO-NLL
c --- for the inclusive jet cross section in hadron collisions
c --- 
c --- based on the code by J.F. Owens
c --- uses an improved sampling of the (xa,s4) phase space
c --- adds integration over pT bin
c --- the definition of the observable is obtained from a fastNLO table
c
c 06/01/04 MW 1st version writes out 2-loop coefficients only
c               (to be included in fixed-order table)
c 06/01/09 MW new version - write out all three tables
c               (Born, 1-loop, 2-loop - can be used in a separate table)
c 06/01/25 MW adopt to tableformat version 1c 
c
c ---------------------------------------------------------------------
      implicit none
      include 'fnthresh-v13.inc'
      Integer  i,j,k,l,m, irap,ipt, icor, ipdf
      Double Precision kd,ld
      Integer  Iset, npoint, nypt, nptpt,nxpt,ns4pt, iscale, nbin
      Double Precision  s,t,u, xamin,xa,xabin, s4max,s4,s4bin,
     +     xb,sh,th,uh, xb0,sh0,th0,uh0, 
     +     mur,mur2, muf,muf2,
     +     rs,y,ymin,ymax,ytest,yup,xt,pt,pt2,ptmin,ptmax
      Double Precision fac,fac1,fac2,  as,alpha_s, pi
      Double Precision temp,A1,A2, als4m, als4
      Double Precision ANS(3), result(nrapmax,nptmax,nscalemax,3)
      Double Precision CF(9) ! coefficients removed from routine COEF
      Double Precision norm, 
     +     xsa(9),xsa1(9),xsa2(9),
     +     xsb(9),xsb1(9),xsb2(9),
     +     xsc(9),xsc1(9),xsc2(9)
      CHARACTER*100 INTABLE, OUTTAB0, OUTTAB1, OUTTAB2, Histofile
      INTEGER IOPTIMIZE

c - variables from original code
      Double Precision 
     +     SIG0(9),SIG0C(9),SIG1(9),SIG1C(9),HAB0(9),HAB1(9),
     +     C10(3,9),C10C(3,9),C20(4,9),C20C(4,9),C1(3,9),C1C(3,9),
     +     C2(4,9),C2C(4,9)

c - variables for PDFs/alpha_s
      Integer NEFF              ! returned from alphas routine
      Integer NFL,IORDER
      Double Precision ALAMBDA
      COMMON/QCDTABLE/ALAMBDA,NFL,IORDER       ! from CTEQ6 code
      DATA PI/3.14159265358979323846/

c - variables from original code
      INTEGER NTT
      COMMON/TARGET/NTT         ! used in PDF routine COEF
      
c - HBOOK common
      INTEGER NWPAWC,ISTAT2,ICYCLE, ihist, ifill
      PARAMETER (NWPAWC=2500000)
      REAL HMEMOR(NWPAWC)
      COMMON /PAWC/ HMEMOR
      real ptbin(100)

c - special variable weight for pt^3/2/pi
      Double Precision ptwgt(2)
      Data ptwgt/1d0, 0d0/      ! ptwgt(2) is defined later in code as pt^3/2pi

c - optimize?
      IOPTIMIZE = 1             ! 0: full calculation (incl PDFs)
                                ! 1: only coefficients for table
                                !     the latter is 4x faster

c - set number of integration points in: y,pT,xa,s4
c -   always make sure that a larger No. of points does
c -   not change the results!
c
c - this gives very good precision - but it may take long for fine bins
c  e.g.                fnt1004midp  (40 bins):   3:30h
c                      fnt1005midp  (40 bins):   3:30h
      nypt  = 16
      nptpt = 32
      nxpt  = 90
      ns4pt = 60
c
c - for faster test - still reasonable results!! (factor 15 faster!!)
c      nypt  = 5
c      nptpt = 12
c      nxpt  = 70
c      ns4pt = 40
c
c - for extremely quick *tests* you can still use much smaller numbers



      write(*,*) ' '
      write(*,*) ' ***********************************************'
      write(*,*) ' *     Resummed Threshold Corrections'
      write(*,*) ' * '
      write(*,*) ' * Markus Wobisch Nov 22, 2005'
      write(*,*) ' * based on code by J.F. Owens and N. Kidonakis'
      write(*,*) ' * '
      write(*,*) ' * Jan 5, 2006 version for tableformat v1c'
      write(*,*) ' * '
      write(*,*) ' *  >>> special version for additional' 
      write(*,*) ' *      weighted x-section bins in 630/1800 ratio'
      write(*,*) ' ***********************************************'
c - three changes in this special version of code
c    1. check if exactly 2 y-ranges
c    2. modified ymin/ymax definition (always use 1st/2nd value)
c    3. add weight pt**3/2pi to SIG0, SIG0C, SIG1, SIG1C

      if (ioptimize.eq.1) then
         write(*,*) ' '
         write(*,*) '   speed optimization is selected (ioptmize=1)'
         write(*,*) '   - this will only print the LO cross section'
         write(*,*) '   - and the coefficients are stored'
         write(*,*) '     in the fastNLO tables'
         write(*,*) ' '         
      endif

      IF ( IARGC().EQ.0)  THEN
         write(*,*) ' '
         write(*,*) ' the code requires 5 arguments:'
         write(*,*) '      [input-fastNLO-table]   (required)'
         write(*,*) '      [output-text-table0]     (optional)'
         write(*,*) '      [output-text-table1]     (optional)'
         write(*,*) '      [output-text-table2]     (optional)'
         write(*,*) '      [output-histogram-file] (optional)'
         write(*,*) ' '
         stop
      ENDIF
      IF ( IARGC().LT.1)  INTABLE = 'table.tab'
      IF ( IARGC().LT.2)  OUTTAB0= 'thrcor-table0.txt'
      IF ( IARGC().LT.3)  OUTTAB1= 'thrcor-table1.txt'
      IF ( IARGC().LT.4)  OUTTAB2= 'thrcor-table2.txt'
      IF ( IARGC().LT.5)  HISTOFILE= 'threshcor.hbk'
      IF ( IARGC().GT.0)  CALL GETARG(1,INTABLE)
      IF ( IARGC().GT.1)  CALL GETARG(2,OUTTAB0)
      IF ( IARGC().GT.2)  CALL GETARG(3,OUTTAB1)
      IF ( IARGC().GT.3)  CALL GETARG(4,OUTTAB2)
      IF ( IARGC().GT.4)  CALL GETARG(5,HISTOFILE)

      write(*,*) ' getting phase space from fastNLO table ',intable
      write(*,*) ' store O(alphas^2) coefficients in table ',OUTTAB0
      write(*,*) ' store O(alphas^3) coefficients in table ',OUTTAB1
      write(*,*) ' store O(alphas^4) coefficients in table ',OUTTAB2
      write(*,*) ' store histograms in ',histofile

c - read fastNLO table to get the phase space definition
      CALL READTABLE(intable)

c ----------------------------------------------------------------------
c - special treatment for D0 ratio 630/1800 where 2nd y-range
c   is the same as 1st, but filled with weight ET^3/2pi
      if (Nrap.ne.2) then
         write(*,*) '  more than 2 y-ranges - not expected for fnt1004'
         stop
      endif

c - initialize PDF
      ISET=200 ! CTEQ6.1M
      CALL SETCTQ6(ISET)

c - output y,pT bins
      do i=1,Nrap
         do j=1,Npt(i)
            write(*,*) i,j,'  ',yarr(i),' -',yarr(i+1),'  ',
     +           ptarr(i,j),' -',ptarr(i,j+1),'  ',xlimit(i,j),murval(i,j)
         enddo
      enddo
c -output scale variations
      do i=1,nscalevar
         write(*,*) ' scale ',i,': ',murscale(i),mufscale(i)
      enddo

c - assign CMS energy, pp vs. pp-bar
      rs = ECMS
      s  = rs*rs
      if (Ireaction.eq.3) then
         NTT=2
         write(*,*) ' >> running in pp-bar mode at sqrt(s)=',rs,' GeV'
      elseif (Ireaction.eq.2) then 
         NTT=1
         write(*,*) ' >> running in pp mode at sqrt(s)=',rs,' GeV'
      else 
         write(*,*) '   reaction No.',Ireaction,' unknown'
         stop
      endif

c - total No. of integration points in y,pT
      npoint = nypt*nptpt
      Nevents = dble(nypt*nptpt*nxpt*ns4pt)
      write(*,*) ' integration points in y,pT,xa,s4: ',nypt,nptpt,nxpt,ns4pt 
      write(*,*) '    -> equivalent No. of events: ',Nevents

c - book histograms 
      CALL HLIMIT(NWPAWC)
      CALL HROPEN(11,'thrcor',histofile,'N',1024,ISTAT2)
      if (ISTAT2.NE.0) then
         WRITE(*,*) ' THRSHCOR: could not open histofile ',istat2
      endif
      do k=1,nscalevar
         do i=1,nrap
            do j=1,(npt(i)+1)
               ptbin(j)=real(ptarr(i,j))
            enddo
            call hbookb(i+1000*k ,'p?T! (GeV)' , NPT(i) ,PTBIN ,0)
            call hbookb(i+100000+1000*k,'p?T! (GeV)',NPT(i),PTBIN,0)
            call hbookb(i+200000+1000*k,'p?T! (GeV)',NPT(i),PTBIN,0)
         enddo
      enddo
      
c - reset result array
      do j=1,nrap
         do i=1,npt(nrap)
            do k=1,nscalemax
               do l=1,3
                  result(j,i,k,l)=0d0
               enddo
            enddo
         enddo
      enddo

c - define coefficients which have been removed from subroutine COEF
c   assume: always NF=5
      cf(1)=1d0
      cf(2)=1d0
      cf(3)=0.5d0
      cf(4)=5d0-1d0   ! =(NEFF-1)
      cf(5)=1d0
      cf(6)=0.5d0
      cf(7)=5d0       ! =NEFF
      cf(8)=1d0
      cf(9)=0.5d0

c - universal factor - for results in units of nb
      FAC1 = 389385.730d0 *2d0*PI*(4d0*PI*PI)

c - ---- adopt units -> if needed modify to pb or fb
c - recognize the fastNLO scenario 
      if (nrap.eq.5 .and. npt(1).eq.24) then ! fnt1002 -> fb
         FAC1 = FAC1 * 1000000d0
         write(*,*) ' ... switching to units of: fb  (for fnt1002)'
      elseif (nrap.eq.6 .and. npt(1).eq.115) then ! fnl000a -> pb
         FAC1 = FAC1 * 1000d0
         write(*,*) ' ... switching to units of: pb  (for fnl000a)'
      endif

c ---------- start calculation -----------------
c - loop over scale, rap and pT bins
      do iscale=1,nscalevar
         write(*,*) ' *** now setting scales mur/pT, muf/pT to: ',
     +        mufscale(iscale)
         nbin = 0

         do irap=1,nrap         ! loop over y ranges
c            ymin = yarr(irap)
c            ymax = yarr(irap+1)
c - fnt1004: special case where 2nd y-range is identical with 1st
            ymin = yarr(1)
            ymax = yarr(2)

            do ipt=1,npt(irap)  ! loop over pT bins
               nbin = nbin+1    ! linear bin counter
               if (nbin.gt.Nrapptmax) then
                  write(*,*) '  too many bins -> increase Nrapptmax'
                  stop
               endif
               ptmin = ptarr(irap,ipt) 
               ptmax = ptarr(irap,ipt+1)

               do i=1,nptpt     ! integrate over pT inside the bin
                  pt = ptmin + dble(i-0.5d0)*(ptmax-ptmin)/dble(nptpt)
c               pt = ptmin       ! test: run at fixed pT / no integration
c --- fnt1004 - fill weigtht
                  if (irap.eq.2) ptwgt(2)=pt*pt*pt/2d0/pi

                  xt=2.*pt/rs
                  ytest=log((1.+sqrt(1.-xt**2))/xt)
                  yup=ymax
                  if(yup.gt.ytest) yup=ytest
c - normalization factor if accessible y range is smaller than measured range
                  norm = (yup-ymin)/(ymax-ymin)
c                  if (norm.lt.1.d0) write(*,*) ' reduced accessible y range:'
c     +                 ,irap,ipt,i,norm

                  if (norm.lt.0d0) then
                     write(*,*) ' no accessble y/pT range in bin:',irap,ipt
                     write(*,*) '                  in pT sub-bin:',i
                     goto 900
                  endif

                  do j=1,nypt   ! integrate over y
                     y = ymin + dble(j-0.5d0)*(yup-ymin)/dble(nypt)
c                  y = ymin      ! test: run at fixed y / no integration
                     T=-RS*PT*EXP(-Y)
                     U=-RS*PT*EXP(Y)
                     XAMIN=-U/(S+T)
                     
                     do k=1,nxpt ! integrate over xa
                        kd = dble(k)
C - w/o importance sampling
c                     xa = xamin + dble(k-0.5d0)/dble(nxpt) *(1d0-xamin)
c                     xabin = (1d0-xamin)/dble(nxpt) *norm
C - use importance sampling (quadratic)
c                     xa = xamin + (dble(k-0.5d0)/dble(nxpt))**2d0 *(1d0-xamin)
c                     xabin = (1d0-xamin)*(2d0*kd-1d0)/dble(nxpt)/dble(nxpt)*norm
C - use importance sampling (cubic)
                        xa = xamin+(dble(k-0.5d0)/dble(nxpt))**3d0*(1d0-xamin)
                        xabin = (1d0-xamin)*(3d0*kd*kd-3d0*kd+1d0)
     +                       /dble(nxpt)/dble(nxpt)/dble(nxpt) *norm
c --
                        s4max = xa*(s+t)+u
                        XB0=-XA*T/(XA*S+U)
                        SH0=XA*XB0*S
                        TH0=XA*T
                        UH0=XB0*U
                        
                        muf = mufval(irap,ipt) * mufscale(iscale)
                        mur = murval(irap,ipt) * murscale(iscale)
c - (the ren scale choice can not be changed later - aposteriori 
c    ren. scale variation for threshold corrections is not possible!!)
                        mur2=mur*mur
                        muf2=muf*muf

                        pt2 = pt*pt
                        as = alpha_s(2,mur2,ALAMBDA,NEFF)
                        FAC2 = FAC1*PT/(XA*S+U)
                        FAC = FAC2* AS*AS/4d0/PI/PI

                        CALL COEF(XA,XB0,muf,HAB0) ! gets the PDFs

c - compute Born terms
                        CALL SUB_SIG(SH0,TH0,UH0,SIG0)
                        CALL SUB_SIG(SH0,UH0,TH0,SIG0C)
                        do ipdf=1,9

c --- fnt1004 - pt-weighting
                           SIG0(ipdf) = SIG0(ipdf) *ptwgt(irap)
                           SIG0C(ipdf)= SIG0C(ipdf)*ptwgt(irap)

                           temp = HAB0(ipdf)*cf(ipdf)
     +                          *(SIG0(ipdf)+SIG0C(ipdf))*FAC
                           result(irap,ipt,iscale,1)=result(irap,ipt,iscale,1)
     +                          +temp*xabin/dble(npoint)           !*1000d0

c -                     contribution for fastNLO table
                           xsa(ipdf) = cf(ipdf)*xabin/dble(npoint)
     +                          * (SIG0(ipdf)+SIG0C(ipdf))
     +                          * FAC2 * Nevents
                        enddo


c - compute Delta(S4) 1-loop contribution
                        CALL C1CALC(SH0,TH0,UH0,mur,muf,PT,NEFF,C10)
                        CALL C1CALC(SH0,UH0,TH0,mur,muf,PT,NEFF,C10C)
                        ALS4M=LOG(S4MAX/PT2)
                        do ipdf=1,9
                           temp = HAB0(IPDF)*cf(ipdf)*(SIG0(IPDF)*(C10(1,IPDF)
     +                          +C10(2,IPDF)*ALS4M+C10(3,IPDF)*ALS4M**2/2.)
     +                          +SIG0C(IPDF)*(C10C(1,IPDF)+C10C(2,IPDF)*ALS4M
     +                          +C10C(3,IPDF)*ALS4M**2/2.))*FAC*AS/PI
                           result(irap,ipt,iscale,2)=result(irap,ipt,iscale,2)
     +                          +temp*xabin/dble(npoint)          !*1000d0

c -                     contribution for fastNLO table
                           xsb(ipdf) = cf(ipdf)*xabin/dble(npoint)
     +                          * (SIG0(ipdf)*
     +                          (C10(1,IPDF)
     +                          +C10(2,IPDF)*ALS4M
     +                          +C10(3,IPDF)*ALS4M**2/2.)
     +                          + SIG0C(ipdf)*
     +                          (C10C(1,IPDF)
     +                          +C10C(2,IPDF)*ALS4M
     +                          +C10C(3,IPDF)*ALS4M**2/2.))
     +                          * FAC2 * 2d0 * Nevents
                        enddo

c - compute Delta(S4) 2-loop contribution
                        CALL C2CALC(SH0,TH0,UH0,mur,muf,PT,NEFF,C20)
                        CALL C2CALC(SH0,UH0,TH0,mur,muf,PT,NEFF,C20C)
                        ALS4M=LOG(S4MAX/PT2)
                        do ipdf=1,9
                           temp = HAB0(ipdf)*cf(ipdf)*(SIG0(ipdf)
     +                          *(C20(1,ipdf)*ALS4M
     +                          +C20(2,ipdf)*ALS4M**2/2.
     +                          +C20(3,ipdf)*ALS4M**3/3.
     +                          +C20(4,ipdf)*ALS4M**4/4.)
     +                          +SIG0C(ipdf)*(C20C(1,ipdf)*ALS4M
     +                          +C20C(2,ipdf)*ALS4M**2/2.
     +                          +C20C(3,ipdf)*ALS4M**3/3.
     +                          +C20C(4,ipdf)*ALS4M**4/4.))*FAC*(AS/PI)**2
                           result(irap,ipt,iscale,3)=result(irap,ipt,iscale,3)
     +                          +temp*xabin/dble(npoint)          !*1000d0

c -                     contribution for fastNLO table
                           xsc(ipdf) = cf(ipdf)*xabin/dble(npoint)
     +                          * (SIG0(ipdf)*
     +                          (C20(1,ipdf)*ALS4M
     +                          +C20(2,ipdf)*ALS4M**2/2
     +                          +C20(3,ipdf)*ALS4M**3/3.
     +                          +C20(4,ipdf)*ALS4M**4/4.)
     +                          + SIG0C(ipdf)*
     +                          (C20C(1,ipdf)*ALS4M
     +                          +C20C(2,ipdf)*ALS4M**2/2.
     +                          +C20C(3,ipdf)*ALS4M**3/3.
     +                          +C20C(4,ipdf)*ALS4M**4/4.))
     +                          * FAC2 * 4d0 * Nevents 

                        enddo
c - fill all contribution into the fastNLO table (1-loop and 2-loop)
                        call FILLTABLE(0,xa,xb0,xsa,nbin,irap,ipt,iscale)
                        call FILLTABLE(1,xa,xb0,xsb,nbin,irap,ipt,iscale)
                        call FILLTABLE(2,xa,xb0,xsc,nbin,irap,ipt,iscale)



c - actual S4 integration
                        do l=1,ns4pt
                           ld = dble(l)
C - w/o importance sampling
c                        s4 = dble(l-0.5d0) * s4max /dble(ns4pt)
c                        s4bin = s4max/dble(ns4pt)
C - use importance sampling (quadratic)
c                        s4 = (dble(l-0.5d0)/dble(ns4pt))**2 * s4max 
c                        s4bin = s4max * (2d0*ld-1d0)/dble(ns4pt)/dble(ns4pt)
C - use importance sampling (cubic)
c                        s4 = (dble(l-0.5d0)/dble(ns4pt))**3 * s4max 
c                        s4bin = s4max * (3d0*ld*ld-3d0*ld+1d0)
c     +                       /dble(ns4pt)/dble(ns4pt)/dble(ns4pt)
C - use importance sampling (**4)
c                        s4 = (dble(l-0.5d0)/dble(ns4pt))**4 * s4max 
c                        s4bin = s4max* (4d0*ld*ld*ld -6d0*ld*ld +4d0*ld -1d0) 
c     +                       /dble(ns4pt)/dble(ns4pt)/dble(ns4pt)/dble(ns4pt)
C - use importance sampling (**5)
                           s4 = (dble(l-0.5d0)/dble(ns4pt))**5 * s4max 
                           s4bin = s4max * (5d0*ld*ld*ld*ld 
     +                          -10d0*ld*ld*ld+10d0*ld*ld-5d0*ld+1d0) 
     +                          /dble(ns4pt)/dble(ns4pt)/dble(ns4pt)
     +                          /dble(ns4pt)/dble(ns4pt)
c --
                           XB=(S4-XA*T)/(XA*S+U)
                           SH=XA*XB*S
                           TH=XA*T
                           UH=XB*U
                           ALS4=LOG(S4/PT2)
                           CALL SUB_SIG(SH,TH,UH,SIG1)
                           CALL SUB_SIG(SH,UH,TH,SIG1C)
                           IF (ioptimize.eq.0)  CALL COEF(XA,XB,muf,HAB1)
                           CALL C1CALC(SH,TH,UH,mur,muf,PT,NEFF,C1)
                           CALL C1CALC(SH,UH,TH,mur,muf,PT,NEFF,C1C)
                           CALL C2CALC(SH,TH,UH,mur,muf,PT,NEFF,C2)
                           CALL C2CALC(SH,UH,TH,mur,muf,PT,NEFF,C2C)
                           A1 = 0d0
                           A2 = 0d0
                           do m=1,9 ! PDF loop

c --- fnt1004 - pT-weighting
                              SIG1(m) = SIG1(m) *ptwgt(irap)
                              SIG1C(m)= SIG1C(m)*ptwgt(irap)

                              A1=A1+(HAB1(m)*cf(m)*SIG1(m)*(C1(2,m)
     +                             +ALS4*C1(3,m))
     +                             -HAB0(m)*cf(m)*SIG0(m)*(C10(2,m)
     +                             +ALS4*C10(3,m)))/S4
                              A1=A1+(HAB1(m)*cf(m)*SIG1C(m)*(C1C(2,m)
     +                             +ALS4*C1C(3,m))
     +                             -HAB0(m)*cf(m)*SIG0C(m)*(C10C(2,m)
     +                             +ALS4*C10C(3,m)))/S4
                              A2=A2+(HAB1(m)*cf(m)*SIG1(m)*(C2(1,m)
     +                             +C2(2,m)*ALS4+C2(3,m)*ALS4**2
     +                             +C2(4,m)*ALS4**3)
     +                             -HAB0(m)*cf(m)*SIG0(m)*(C20(1,m)
     +                             +C20(2,m)*ALS4+C20(3,m)*ALS4**2
     +                             +C20(4,m)*ALS4**3))/S4
                              A2=A2+(HAB1(m)*cf(m)*SIG1C(m)*(C2C(1,m)
     +                             +C2C(2,m)*ALS4+C2C(3,m)*ALS4**2
     +                             +C2C(4,m)*ALS4**3)
     +                             -HAB0(m)*cf(m)*SIG0C(m)*(C20C(1,m)
     +                             +C20C(2,m)*ALS4+C20C(3,m)*ALS4**2
     +                             +C20C(4,m)*ALS4**3))/S4

c -                 1-loop contribution for fastNLO table
c - proportional to HAB0
                              xsb1(m) =-cf(m)*s4bin*xabin/dble(npoint) / S4
     +                             *(SIG0(m)*(C10(2,m)+ALS4*C10(3,m))
     +                             + SIG0C(m)*(C10C(2,m)+ALS4*C10C(3,m)))
     +                             * FAC2 * 2d0 * Nevents 
c - proportional to HAB1
                              xsb2(m) = cf(m)*s4bin*xabin/dble(npoint) / S4
     +                             *(SIG1(m)*(C1(2,m)+ALS4*C1(3,m))
     +                             + SIG1C(m)*(C1C(2,m)+ALS4*C1C(3,m)))
     +                             * FAC2 * 2d0 * Nevents 

c -                 2-loop contribution for fastNLO table
c - proportional to HAB0
                              xsc1(m) =-cf(m)*s4bin*xabin/dble(npoint) / S4
     +                             *(SIG0(m)
     +                             *(C20(1,m) + C20(2,m)*ALS4 
     +                             + C20(3,m)*ALS4**2 + C20(4,m)*ALS4**3)
     +                             + SIG0C(m)
     +                             *(C20C(1,m) + C20C(2,m)*ALS4
     +                             + C20C(3,m)*ALS4**2 + C20C(4,m)*ALS4**3))
     +                             * FAC2 * 4d0 * Nevents 
c     +                             * AS*AS*AS*AS/4d0/4d0/PI/PI/PI/PI
c     +                             * HAB0(m)
c - proportional to HAB1
                              xsc2(m) = cf(m)*s4bin*xabin/dble(npoint) / S4
     +                             *(SIG1(m)
     +                             *(C2(1,m) + C2(2,m)*ALS4
     +                             + C2(3,m)*ALS4**2 + C2(4,m)*ALS4**3)
     +                             +SIG1C(m)
     +                             *(C2C(1,m) + C2C(2,m)*ALS4
     +                             + C2C(3,m)*ALS4**2 + C2C(4,m)*ALS4**3))
     +                             * FAC2 * 4d0 * Nevents
c     +                             * AS*AS*AS*AS/4d0/4d0/PI/PI/PI/PI
c     +                             * HAB1(m)

                           enddo ! end PDF loop

                           result(irap,ipt,iscale,2)=result(irap,ipt,iscale,2)
     +                          + A1*fac*as/pi*s4bin*xabin/dble(npoint)        !*1000d0
                           result(irap,ipt,iscale,3)=result(irap,ipt,iscale,3)
     +                          + A2*fac*(as/pi)**2
     +                          * s4bin*xabin/dble(npoint)              !*1000d0

c - fill these contributions into the fastNLO table
                        call FILLTABLE(1,xa,xb0,xsb1,nbin,irap,ipt,iscale)
                        call FILLTABLE(1,xa,xb,xsb2,nbin,irap,ipt,iscale)
                        call FILLTABLE(2,xa,xb0,xsc1,nbin,irap,ipt,iscale)
                        call FILLTABLE(2,xa,xb,xsc2,nbin,irap,ipt,iscale)
                           
                        enddo   ! end s4 integral
                     enddo      ! end xa integral
                  enddo         ! end y integral
 900              continue
               enddo            ! end pT integral

               IF (ioptimize.eq.0) then
                  write(*,*) irap,ipt,' ',result(irap,ipt,iscale,1),
     +                 result(irap,ipt,iscale,2),result(irap,ipt,iscale,3)
               else
                  write(*,*) irap,ipt,' ',result(irap,ipt,iscale,1)
               endif

c - fill histograms
               call hfill(irap,real(ptarr(irap,ipt)*1.001d0),0.0,
     +              real(result(irap,ipt,iscale,1)))
               call hfill(irap+100000+1000*iscale,real(ptarr(irap,ipt)*1.001d0),0.0,
     +              real(result(irap,ipt,iscale,2)))
               call hfill(irap+200000+1000*iscale,real(ptarr(irap,ipt)*1.001d0),0.0,
     +              real(result(irap,ipt,iscale,3)))
               
            enddo
         enddo
      enddo ! scale loop

c - close HBOOK file
      CALL HROUT (0,ICYCLE,' ')
      CALL HREND ('thrcor')

c - write fastNLO raw-table 
      CALL WRITETABLE(0,OUTTAB0)
      CALL WRITETABLE(1,OUTTAB1)
      CALL WRITETABLE(2,OUTTAB2)

      RETURN
      END


C ----------------------------------------------------------
      SUBROUTINE FILLTABLE(nloop,xa,xb,coef,nbin,nrapbin,nptbin,nscale)
c -------------------------------------------
c MW 
c fill array for fastNLO
c
c   nloop:  No of loops for which the table is filled
c                  0 Born, 1: LL, 2: NLL 
c -------------------------------------------
      implicit none
      integer nloop,nbin,nrapbin,nptbin,nscale,nf
      double precision xa,xb,coef(9)
      double precision xmin,xmax,cf(7)
      double precision reweight, hxmin, hxmax, hxone, hxlimit, hxi, hxj, 
     +     delta, deltamax, deltamin
      integer ifirst,i,j, isub,nmin,nmax
c - variables for bicubic interpolation
      Double Precision cmax(4),cmin(4), cefmax(4),cefmin(4), bicef(4,4),
     +     buffer
      Integer Imax, Imin, Di

      include 'fnthresh-v13.inc'
      data ifirst/0/,hxone/0.0/
      save ifirst

c - initialization
      if (ifirst.ne.1) then
         do i=1,nrap
            do j=1,npt(i)
               hxlim(i,j) = -sqrt(-log10(xlimit(i,j)))
            enddo
         enddo
         ifirst = 1
      endif

c - PDF reweighting
      reweight = 1d0/sqrt(xa)/sqrt(xb)*(1d0-0.99d0*xa)*(1d0-0.99d0*xb)
      reweight = reweight * reweight * reweight

c - determine coefficients for seven subprocesses
      cf(1) = reweight *(coef(9) + coef(7))
      cf(2) = reweight * coef(8)
      cf(3) = reweight * coef(8)
      cf(4) = reweight * coef(1)
      cf(5) = reweight * coef(3)
      cf(6) = reweight *(coef(4) + coef(5) + coef(6))
      cf(7) = reweight * coef(2)

c - determine xmax,xmin bin
      if (xa.lt.xb) then
         xmin = xa
         xmax = xb
      else
         xmin = xb
         xmax = xa
      endif
      hxmin = -sqrt(-log10(xmin))
      hxmax = -sqrt(-log10(xmax))
      hxlimit = hxlim(nrapbin,nptbin)
      nmin = int(nxtot * (hxmin-hxlimit)/(hxone-hxlimit))
      nmax = int(nxtot * (hxmax-hxlimit)/(hxone-hxlimit))

c - determine relative distances to bin centers
      delta = (hxone-hxlimit)/nxtot
      hxi = hxlimit+dble(nmax)/dble(nxtot)*(hxone-hxlimit)
      hxj = hxlimit+dble(nmin)/dble(nxtot)*(hxone-hxlimit)
      deltamax = (hxmax-hxi)/delta
      deltamin = (hxmin-hxj)/delta
      if (deltamax.gt.1d0 .or. deltamin.gt.1d0 .or. deltamax.lt.0d0
     +     .or. deltamin.lt.0d0) write(*,*) "error in logic"

c - convert from C++ to Fortran
      nmin = nmin + 1
      nmax = nmax + 1
      if (nmin.eq.0) write(*,*) "nmin=0!! ",nmin,nmax
      if (nmax.gt.25) write(*,*) "       nmax !! ",nmin,nmax

c === compute bicubic interpolation functions ======
      cmax(1) = deltamax + 1d0
      cmax(2) = deltamax
      cmax(3) = 1d0 - deltamax
      cmax(4) = 2d0 - deltamax

      cmin(1) = deltamin + 1d0
      cmin(2) = deltamin
      cmin(3) = 1d0 - deltamin
      cmin(4) = 2d0 - deltamin

      if (nmax.eq.1 .or. nmax.eq.nxtot) then ! linear in 1st and last bin
         cefmax(1) = 0d0
         cefmax(2) = 1d0 - cmax(2)
         cefmax(3) = 1d0 - cmax(3)
         cefmax(4) = 0d0
      else
         cefmax(2) = 1d0-2.5d0*cmax(2)*cmax(2) + 1.5d0*cmax(2)*cmax(2)*cmax(2)
         cefmax(3) = 1d0-2.5d0*cmax(3)*cmax(3) + 1.5d0*cmax(3)*cmax(3)*cmax(3)
         cefmax(1) = 2d0 - 4d0*cmax(1) + 2.5d0*cmax(1)*cmax(1)
     +        - 0.5d0*cmax(1)*cmax(1)*cmax(1)
         cefmax(4) = 2d0 - 4d0*cmax(4) + 2.5d0*cmax(4)*cmax(4)
     +        - 0.5d0*cmax(4)*cmax(4)*cmax(4)
      endif
      if (nmin.eq.1 .or. nmin.eq.nxtot) then ! linear in 1st and last bin
         cefmin(1) = 0d0
         cefmin(2) = 1d0 - cmin(2)
         cefmin(3) = 1d0 - cmin(3)
         cefmin(4) = 0d0
      else
         cefmin(2) = 1d0-2.5d0*cmin(2)*cmin(2) + 1.5d0*cmin(2)*cmin(2)*cmin(2)
         cefmin(3) = 1d0-2.5d0*cmin(3)*cmin(3) + 1.5d0*cmin(3)*cmin(3)*cmin(3)
         cefmin(1) = 2d0 - 4d0*cmin(1) + 2.5d0*cmin(1)*cmin(1)
     +        - 0.5d0*cmin(1)*cmin(1)*cmin(1)
         cefmin(4) = 2d0 - 4d0*cmin(4) + 2.5d0*cmin(4)*cmin(4)
     +        - 0.5d0*cmin(4)*cmin(4)*cmin(4)
      endif

      do i=1,4
         do j=1,4
            bicef(i,j) = cefmax(i) * cefmin(j)
         enddo
      enddo


c === here: for xb>xa - swap sbprocesses 2,3
c          --> but not needed for threshold corrections


c === fill table for seven fastNLO subprocesses
      do isub=1,7
         do i=1,4
            do j=1,4
               imax = nmax + i -2  ! the target index
               imin = nmin + j -2
c -       check if above diagonal -> project back
               if (imin.gt.imax) then
                  di = imin-imax
                  imax = imax + di
                  imin = imin - di
c -       swap subprocesses 2,3 (but not needed for thresh. corr.)                  
               endif
               if (imax.le.nxtot .and. imin.gt.0) 
     +              table(nloop,nbin,imax,imin,nscale,isub) = 
     +              table(nloop,nbin,imax,imin,nscale,isub) 
     +              + bicef(i,j) * cf(isub) 
            enddo
         enddo

      enddo

      RETURN
      END
C --------------------------------------------------------------------
      SUBROUTINE READTABLE(filename)
c -------------------------------------------
c M. Wobisch Dec 02, 2005
c read fastNLO LO/NLO table (in text format) 
c to obtain the phase space definition
c
c -------------------------------------------
c      implicit none
      CHARACTER*(*) filename
      CHARACTER*(255) label
      integer i,j, nsep, iref
      include 'fnthresh-v13.inc'

      write (*,*) 'read fastNLO table ',filename
      open(unit=2,file=filename,status='unknown')
      READ(2,*) Ireaction
      READ(2,*) Ecms
      do i=1,5
         READ(2,*) namelabel(i)
         write(*,*) ">>>   info:   ",namelabel(i)
      enddo
      READ(2,*) Iproc
      READ(2,*) Ialgo
      READ(2,*) Jetresol1
      READ(2,*) Jetresol2
      READ(2,*) Nord
      do i=1,Nord
         READ(2,*) Npow(i)
      enddo
      do i=1,Nord
         READ(2,*) label
         write(*,*) ">>>   contribution ",i,label
      enddo
      READ(2,*) Nsep
      If (nsep.ne.Iseparator) goto 999
      do i=1,Nord
         READ(2,*) Nevt(i)
      enddo
      READ(2,*) Nxtot
      READ(2,*) Ixscheme
      READ(2,*) Ipdfwgt
      READ(2,*) Iref
      IF (iref.eq.1) then
         write(*,*) ' >>> you are reading a fastNLO REFERENCE table'
         write(*,*) '   > there is no point in computing the threshold'
         write(*,*) '   > corrections for this table!'
         write(*,*) '   >  - please choose a fastNLO result table'
         stop
      endif

      READ(2,*) Nsep
      If (nsep.ne.Iseparator) goto 999
      READ(2,*) Nrap
      do i=1,Nrap+1
         READ(2,*) yarr(i)
      enddo
      do i=1,Nrap
         READ(2,*) Npt(i)
      enddo
      do i=1,Nrap
         do j=1,Npt(i)+1
            READ(2,*) ptarr(i,j)
         enddo
      enddo
      READ(2,*) Nsep
      If (nsep.ne.Iseparator) goto 999
      do i=1,Nrap
         do j=1,Npt(i)
            READ(2,*) Xlimit(i,j)
         enddo
      enddo
      READ(2,*) Nsep

      READ(2,*) scalelabel
      write(*,*) ' the scales are proportional to: ',scalelabel
      READ(2,*) NSCALEBIN
      IF (nscalebin.gt.1) then
         write(*,*) ' >>> NSCALEBIN = ',nscalebin,' not yet implemented'
         write(*,*) '   > only NSCALEBIN=1 allowed'
         stop
      endif


      If (nsep.ne.Iseparator) goto 999
      do i=1,Nrap
         do j=1,Npt(i)
            READ(2,*) Murval(i,j)
         enddo
      enddo
      READ(2,*) Nsep
      If (nsep.ne.Iseparator) goto 999
      do i=1,Nrap
         do j=1,Npt(i)
            READ(2,*) Mufval(i,j)
         enddo
      enddo
      READ(2,*) Nsep
      If (nsep.ne.Iseparator) goto 999
      READ(2,*) Nscalevar
      do i=1,Nscalevar
         READ(2,*) Murscale(i)
      enddo
      do i=1,Nscalevar
         READ(2,*) Mufscale(i)
      enddo
      READ(2,*) Nsep
      If (nsep.ne.Iseparator) goto 999
c - no need to read the sigma tilde ......
      close(unit=2)

      write(*,*)  Ireaction, Ecms,Iproc,Ialgo, Jetresol1,Jetresol2,NORD,nscalevar
      IF (ireaction.ne.2 .and. ireaction.ne.3) then
         write(*,*) '  the input table is not for hadron-hadron collisions'
         stop
      endif
      if (iproc.ne.1) then
         write(*,*) '  the input table is not for inclusive jets'
         stop
      endif



      RETURN
 999  continue
      WRITE(*,*) ' wrong table format'
      STOP
      END
C --------------------------------------------------------------------
      SUBROUTINE WRITETABLE(nloop,filename)
c -------------------------------------------
c M. Wobisch Dec 02, 2005
c write fastNLO raw-table - in text format
c
c   nloop:  No of loops for which the table is filled
c                  0 Born, 1: LL, 2: NLL 
c -------------------------------------------
      implicit none
      CHARACTER*(*) filename
      CHARACTER*(33) label(0:2)
      integer nloop, i,j,k,l,m, isub, nbin, mord,mpow,mscalevar,iref
      include 'fnthresh-v13.inc'

      nbin=0
      label(0) = 'Born'
      label(1) = 'NLO-NLL-(threshold-corrections)'
      label(2) = 'NNLO-NLL-(threshold-corrections)'
      if (nloop.eq.0) then
         mscalevar = 1
      else
         mscalevar = nscalevar
      endif

      write (*,*) 'write fastNLO raw-table ',filename

      open(unit=2,file=filename,status='unknown')
      WRITE(2,*) Ireaction
      WRITE(2,*) Ecms
      do i=1,5
         WRITE(2,*) namelabel(i)
      enddo
      WRITE(2,*) Iproc
      WRITE(2,*) Ialgo
      WRITE(2,*) Jetresol1
      WRITE(2,*) Jetresol2

c - thresh-corr are the third contribution (on top of LO, NLO)
      mord = 1 + nloop
      mpow = 2 + nloop
      WRITE(2,*) mord    ! Nord
      WRITE(2,*) mpow    ! Npow
      WRITE(2,*) label(nloop)
      WRITE(2,*) Iseparator     ! ------------------------

      WRITE(2,*) Nevents
      WRITE(2,*) Nxtot
      WRITE(2,*) Ixscheme
      WRITE(2,*) Ipdfwgt
      iref = 0
      WRITE(2,*) Iref
      WRITE(2,*) Iseparator     ! ------------------------
      WRITE(2,*) Nrap
      do i=1,Nrap+1
         WRITE(2,*) yarr(i)
      enddo
      do i=1,Nrap
         WRITE(2,*) Npt(i)
      enddo
      do i=1,Nrap
         do j=1,Npt(i)+1
            WRITE(2,*) ptarr(i,j)
         enddo
      enddo
      WRITE(2,*) Iseparator     ! ------------------------
      do i=1,Nrap
         do j=1,Npt(i)
            WRITE(2,'(f21.18)') Xlimit(i,j)
         enddo
      enddo
      WRITE(2,*) Iseparator     ! ------------------------
      WRITE(2,*) scalelabel
      WRITE(2,*) NSCALEBIN
      do i=1,Nrap
         do j=1,Npt(i)
            WRITE(2,*) Murval(i,j)
         enddo
      enddo
      WRITE(2,*) Iseparator     ! ------------------------
      do i=1,Nrap
         do j=1,Npt(i)
            WRITE(2,*) Mufval(i,j)
         enddo
      enddo
      WRITE(2,*) Iseparator     ! ------------------------
      if (nloop.gt.0) then
         mscalevar = Nscalevar
         WRITE(2,*) Nscalevar
         do i=1,Nscalevar
            WRITE(2,*) Murscale(i)
         enddo
         do i=1,Nscalevar
            WRITE(2,*) Mufscale(i)
         enddo
         WRITE(2,*) Iseparator  ! ------------------------
      else
         mscalevar=1
      endif

c - here write the sigma tilde
      nbin=0
      do i=1,Nrap
         do j=1,Npt(i)
            nbin=nbin+1
            do k=1,Nxtot
               do l=1,k
                  do isub=1,7
                     do m=1,mscalevar ! Nscalevar or one
c                        WRITE(2,*) 1d0
                        WRITE(2,*) table(nloop,nbin,k,l,m,isub)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      WRITE(2,*) Iseparator
      close(unit=2)

      write(*,*) '    ...table is written'

      RETURN
      END
C --------------------------------------------------------------------


********************************************************************
* original code by J.F.Owens
* modified by M. Wobisch for independent variations of mur, muf
* Nov 8, 2005
********************************************************************

      SUBROUTINE SUB_SIG(S,T,U,SIG)
      implicit double precision (a-h,o-z)
      DIMENSION SIG(9)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CALCULATE REDUCED BORN CROSS SECTION EXPRESSIONS. USE NICK'S NOTATION C
C  EXCEPT OVERALL ALPHA_S^2 HAS BEEN REMOVED AND 1/S INSERTED.          C
C  NUMBERING SCHEME IS ADAPTED FROM MY OLD NOTATION                      C
C  CREATED 6/23/00 BY J.F. OWENS                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      S2=S**2
      T2=T**2
      U2=U**2
C
C  Q Q' ==> Q Q'
C
      SIG(1)=4./9.*(S2+U2)/T2/S
C
C  Q QBAR' ==> Q QBAR'
C
      SIG(2)=4./9.*(S2+U2)/T2/S
C
C  Q Q ==> Q Q
C
      SIG(3)=4./9.*((S2+U2)/T2+(S2+T2)/U2-2.*S2/(3.*T*U))/S
C
C  Q QBAR ==> Q' QBAR'
C
      SIG(4)=4./9.*(T2+U2)/S2/S
C
C  Q QBAR ==> Q QBAR
C
      SIG(5)=4./9.*((T2+U2)/S2+(S2+U2)/T2-2.*U2/(3.*S*T))/S
C
C  Q QBAR ==> G G
C
      SIG(6)=8./3.*(4./9.*(T2+U2)/T/U-1.+2.*T*U/S2)/S
C
C  G G ==> Q QBAR
C
      SIG(7)=3./8.*(4./9.*(T2+U2)/T/U-1.+2.*T*U/S2)/S
C
C  Q G ==> Q G
C
      SIG(8)=(2.-1./9.-4./9.*T2/S/U-2.*S*U/T2)/S
C
C  G G ==> G G
C
      SIG(9)=(27./2.-9./2.*(S*U/T2+T*U/S2+S*T/U2))/S
C
      RETURN
      END


C **********************************************************************
C **********************************************************************
      SUBROUTINE COEF(XA,XB,Q,HAB)
*--------------------------------------------------
* Q: factorization scale 
* returns PDF linear combinations in array HAB
*
* original code by J.F.Owens
* modified by M. Wobisch  11/22/05 
*    remove "NEFF" factors - return 'pure' PDFs 
*    put "NEFF" factors in main logic
*--------------------------------------------------
      implicit double precision (a-h,o-z)
      DIMENSION HAB(9),QA(-5:5),QB(-5:5)
      COMMON/QCDTABLE/ALAMBDA,NFL,IORDER
      COMMON/TARGET/NTT
      AS=ALPHA_S(2,Q*Q,ALAMBDA,NEFF)
      DO 100 J=-5,5
c      if(j.eq.0)then
c         qa(j)=0.
c         qb(j)=0.
c         goto 100
c      else 
c      QA(J)=CTQ5PDF(J,XA,Q)
      QA(J)=CTQ6PDF(J,XA,Q)
      JJ=J
      IF(NTT.EQ.2)JJ=-J
c      QB(J)=CTQ5PDF(JJ,XB,Q)
      QB(J)=CTQ6PDF(JJ,XB,Q)
c      endif
  100 CONTINUE
      HAB(1)=0.
      HAB(2)=0.
      DO 101 J=1,5
      DO 101 K=1,5
      IF(K.EQ.J)GOTO 101     
      HAB(1)=HAB(1)+QA(J)*QB(K)+QA(-J)*QB(-K)
      HAB(2)=HAB(2)+QA(J)*QB(-K)+QA(-J)*QB(K)
  101 CONTINUE
      HAB(3)=0.
      HAB(4)=0.
      HAB(5)=0.
      DO 102 J=-5,5
      IF(J.EQ.0)GOTO 102
      HAB(3)=HAB(3)+QA(J)*QB(J)
      HAB(4)=HAB(4)+QA(J)*QB(-J)          !*(NEFF-1)
      HAB(5)=HAB(5)+QA(J)*QB(-J)
  102 CONTINUE 
      HAB(6)=HAB(5)
      HAB(7)=QA(0)*QB(0)                  !*NEFF
      A1=0.
      B1=0.
      DO 103 J=-5,5
      IF(J.EQ.0)GOTO 103
      A1=A1+QA(J)
      B1=B1+QB(J)
  103 CONTINUE
      HAB(8)=A1*QB(0)+B1*QA(0)
      HAB(9)=QA(0)*QB(0)
C
C  NOTE: DIVIDE BY 2 FOR IDENTICAL PARTICLE FINAL STATES
C  SINCE THE T <==> U SYMMETRIZATION IS DONE FOR ALL PROCESSES 
C  IN S4.
C
cMW      HAB(3)=HAB(3)/2.
cMW      HAB(6)=HAB(6)/2.
cMW      HAB(9)=HAB(9)/2.

c      do 200 j=1,9
c      if(j.ne.8)then
c         hab(j)=0.
c      endif
c  200 continue
      RETURN
      END

C **********************************************************************
C **********************************************************************
      SUBROUTINE C1CALC(S,T,U,xmur,xmuf,PT,NEFF,C1)
*--------------------------------------------------
* original code by J.F.Owens
* modified by M. Wobisch  11/22/05 
*  to use independent renormalization and 
*  factorization scales
*--------------------------------------------------
      implicit double precision (a-h,o-z)
      DIMENSION C1(3,9)
      CF=4./3.
      CA=3.
c      XMUF=Q
c      XMUR=Q
      ALF=2.*LOG(XMUF/PT)
      ALR=2.*LOG(XMUR/PT)
      ALT=LOG(-T/S)
      ALU=LOG(-U/S)
      ALTU=LOG(T/U)
      ALUT=LOG(U/T)
      ALP=LOG(PT*PT/S)
      BETA0=(33.-2.*NEFF)/3.
      S2=S*S
      T2=T*T
      U2=U*U
c  q q' --> q q' (B.2)
      C1(3,1)=2.*CF
      C1(2,1)=-2.*CF*ALF-3.*CF/2.-10./3.*ALT+2.*ALU
      C1(1,1)=-CF*(ALP+1.5)*ALF+.5*BETA0*ALR
c  q qbar' --> q qbar' (A.3)
      C1(3,2)=2.*CF
      C1(2,2)=-2.*CF*ALF-3.*CF/2.-10./3.*ALT-4./3.*ALU
      C1(1,2)=-CF*(ALP+1.5)*ALF+.5*BETA0*ALR
c  q q --> q q (B.1)
      C1(3,3)=2.*CF
      C1(2,3)=-2.*CF*ALF-3.*CF/2.-10./3.*ALT+2.*ALU
     &+(4.*CF**2/CA*ALTU*(S2+T2)/U2+8.*CF**2/CA**2*ALU*S2/T/U)
     &/(4./9.*((S2+U2)/T2+(S2+T2)/U2-2./3.*S2/T/U))
      C1(1,3)=-CF*(ALP+1.5)*ALF+.5*BETA0*ALR
c  q qbar --> q' qbar' (A.2)
      C1(3,4)=2.*CF
      C1(2,4)=-2.*CF*ALF-3.*CF/2.-2.*ALT-4./3.*ALU
      C1(1,4)=-CF*(ALP+1.5)*ALF+.5*BETA0*ALR
c  q qbar --> q qbar (A.1)
      C1(3,5)=2.*CF
      C1(2,5)=-2.*CF*ALF-3.*CF/2.-10./3.*ALT-4./3.*ALU
     &+(4.*CF**2/CA*ALT*(T2+U2)/S2-8.*CF**2/CA**2*ALU*U2/S/T)
     &/(4./9.*((S2+U2)/T2+(T2+U2)/S2-2./3.*U2/S/T))
      C1(1,5)=-CF*(ALP+1.5)*ALF+.5*BETA0*ALR
c  q qbar --> g g (C.1)
      C1(3,6)=4.*CF-2.*CA
      C1(2,6)=-2.*CF*ALF-CA*ALP-BETA0/2.
     &+(-4./9.*(T2+U2)/T/U*ALP-4.*((U2-T2)/T/U+2.*(U-T)/S)*ALUT)
     &/(8./3.*(4./9.*(T2+U2)/T/U/-1.+2.*T*U/S2))
      C1(1,6)=-CF*(ALP+1.5)*ALF+.5*BETA0*ALR
c  g g --> q qbar (C.2)
      C1(3,7)=4.*CA-2.*CF
      C1(2,7)=-3./2.*CF-(2.*CF-CA)*ALP-2.*CA*ALF
     &+(-4./9.*(T2+U2)/T/U*ALP-4.*((U2-T2)/T/U+2.*(U-T)/S)*ALUT)
     &/(3./8.*(4./9.*(T2+U2)/T/U-1.+2.*T*U/S2))
      C1(1,7)=-CA*ALP*ALF+.5*BETA0*(ALR-ALF)
c  q g --> q g (D)
      C1(3,8)=CF+CA
      C1(2,8)=-(CF+CA)*ALF+CF*(-2.*ALU-.75)-2.*CA*ALT-BETA0/4.
     &+(-(CF+CA)/9.*(T2/S/U-2.)*ALT+CA*(-1.-2.*S/T+U/2./S-S/2./U)*ALU
     &+(CF*ALT+CA/2.*ALU)*((2./9.-1.)*(T2/S/U-2.)+2.*(1.-2.*S*U/T2)))
     &/(2.-1./9.-4./9.*T2/S/U-2.*S*U/T2)
      C1(1,8)=-(CF*ALT+CA*ALU+.75*CF+BETA0/4.)*ALF+BETA0/2.*ALR
c  g g --> g g (E)
      C1(3,9)=2.*CA
      C1(2,9)=-2.*CA*ALF-2.*CA*ALP-BETA0/2.
     &+(27./8.*(2.*ALT+5.*ALU)*(1.-T*U/S2-S*T/U2+T2/S/U)-27/4.*ALU
     &*(S*T/U2-T*U/S2+U2/S/T-S2/T/U)+(3.*ALT+3./2.*ALU)*(27./4.
     &-9.*(S*U/T2+T*U/4./S2+S*T/4/U2)+9./2.*(U2/S/T+S2/T/U-T2/2./S/U)))
     &/(27./2.-9./2.*(S*U/T2+T*U/S2+S*T/U2))
C        1         2         3         4         5         6         7 *
      C1(1,9)=-CA*ALP*ALF+BETA0/2.*(ALR-ALF)
      RETURN
      END

C **********************************************************************
C **********************************************************************
      SUBROUTINE C2CALC(S,T,U,xmur,xmuf,PT,NEFF,C2)
*--------------------------------------------------
* original code by J.F.Owens
* modified by M. Wobisch  11/22/05 
*  to use independent renormalization and 
*  factorization scales
*--------------------------------------------------
      implicit double precision (a-h,o-z)
      DIMENSION C2(4,9)
      CF=4./3.
      CA=3.
c      XMUF=Q
c      XMUR=Q
      ALF=2.*LOG(XMUF/PT)
      ALR=2.*LOG(XMUR/PT)
      ALT=LOG(-T/S)
      ALU=LOG(-U/S)
      ALTU=LOG(T/U)
      ALUT=LOG(U/T)
      ALP=LOG(PT*PT/S)
      BETA0=(33.-2.*NEFF)/3.
      S2=S*S
      T2=T*T
      U2=U*U
      C2(4,1)=2.*CF**2
      C2(3,1)=3.*CF*(-2.*CF*ALF-3./2.*CF-10./3.*ALT+2.*ALU-BETA0/12.)
      C2(2,1)=4.*CF*(-CF/2.*ALP+.75*CF+10./3.*ALT-2.*ALU+CF*ALF)*ALF
     &+3./2.*CF*BETA0*ALR
      C2(1,1)=CF*(2.*CF*ALP+3.*CF+BETA0/4.)*ALF**2
     &-3./2.*CF*BETA0*ALF*ALR
      C2(4,2)=2.*CF**2
      C2(3,2)=3.*CF*(-2.*CF*ALF-3./2.*CF-10./3.*ALT-4./3.*ALU-BETA0/12.)
C        1         2         3         4         5         6         7 *
      C2(2,2)=4.*CF*(-CF/2.*ALP+.75*CF+10./3.*ALT+4./3.*ALU+CF*ALF)*ALF
     &+3./2.*CF*BETA0*ALR
      C2(1,2)=CF*(2.*CF*ALP+3.*CF+BETA0/4.)*ALF**2
     &-3./2.*CF*BETA0*ALF*ALR
      SIGB=4./9.*((S2+U2)/T2+(S2+T2)/U2-2./CA*S2/T/U)
      C2(4,3)=2.*CF**2
      C2(3,3)=3.*CF*(-2.*CF*ALF-3./2.*CF-10./3.*ALT+2.*ALU-BETA0/12.)
     &+4.*CF**3*(ALTU*(S2+T2)/U2+2./3.*ALU*S2/T/U)/SIGB
      C2(2,3)=4.*CF*(-CF/2.*ALP+3./4.*CF+10./3.*ALT-2.*ALU+CF*ALF)*ALF
     &+3./2.*CF*BETA0*ALR
     &+16./3.*CF**3*(-ALTU*(S2+T2)/U2-2./3.*ALU*S2/T/U)*ALF/SIGB
      C2(1,3)=CF*(2.*CF*ALP+3.*CF+BETA0/4.)*ALF**2
     &-3./2.*CF*BETA0*ALF*ALR
      C2(4,4)=2.*CF**2
      C2(3,4)=3.*CF*(-2.*CF*ALF-3./2.*CF+2.*ALT-4./3.*ALU-BETA0/12.)
      C2(2,4)=4.*CF*(-CF/2.*ALP+3./4.*CF-2.*ALT+4./3.*ALU+CF*ALF)*ALF
     &+3./2.*CF*BETA0*ALR
      C2(1,4)=CF*(2.*CF*ALP+3.*CF+BETA0/4.)*ALF**2
     &-3./2.*CF*BETA0*ALF*ALR
      SIGB=4./9.*((T2+U2)/S2+(S2+U2)/T2-2./3.*U2/S/T)
      C2(4,5)=2.*CF**2
      C2(3,5)=3.*CF*(-2.*CF*ALF-3./2.*CF-10./3.*ALT-4./3.*ALU-BETA0/12.)
     &+4.*CF**3*(ALT*(T2+U2)/S2-2./3.*ALU*U2/S/T)/SIGB
      C2(2,5)=4.*CF*(-CF/2.*ALP+3./4.*CF+10./3.*ALT+4./3.*ALU+CF*ALF)
     &*ALF+3./2*CF*BETA0*ALR
     &+16./3.*CF**3*(-ALT*(T2+U2)/S2+2./3.*ALU*U2/S/T)*ALF/SIGB
      C2(1,5)=CF*(2.*CF*ALP+3.*CF+BETA0/4.)*ALF**2
     &-3./2.*CF*BETA0*ALF*ALR
C        1         2         3         4         5         6         7 *
      SIGB=8./3.*(4./9.*(T2+U2)/T/U-1.+2.*T*U/S2)
      C2(4,6)=1./2.*(4.*CF-2.*CA)**2
      C2(3,6)=3.*(2.*CF-CA)*(-2.*CF*ALF-CA*ALP)+BETA0*(-4.*CF+9./4.*CA)
     &+3.*(2.*CF-CA)*(-4./9.*(T2+U2)/T/U*ALP-4.*((U2-T2)/T/U
     &+2.*(U-T)/S)*ALUT)/SIGB
      C2(2,6)=4.*CF*(-(CF-CA/2.)*(ALP+3./2.)+CA*ALP+BETA0/2.+CF*ALF)*ALF
     &+3.*BETA0*(CF-CA/2.)*ALR
     &+4.*CF*ALF*(4./9.*(T2+U2)/T/U*ALP
     &-4.*((U2-T2)/T/U+2.*(U-T)/S)*ALTU)
     &/SIGB
C        1         2         3         4         5         6         7 *
      C2(1,6)=CF*(2.*CF*ALP+3.*CF+BETA0/4.)*ALF**2
     &-3./2.*CF*BETA0*ALF*ALR
      SIGB=3./8.*(4./9.*(T2+U2)/T/U-1.+2.*T*U/S2)
      C2(4,7)=1./2.*(4.*CA-2.*CF)**2
      C2(3,7)=3.*(2.*CA-CF)*(-3./2.*CF-(2.*CF-CA)*ALP-2.*CA*ALF)
     &+BETA0*(-CA+3./4.*CF)
     &+3.*(2.*CA-CF)*(-4./9.*(T2+U2)/T/U*ALP-4.*((U2-T2)/T/U
     &+2.*(U-T)/S)*ALUT)/SIGB
      C2(2,7)=4.*CA*((5./2.*CF-2.*CA)*ALP+3./2.*CF+CA*ALF)*ALF
     &-BETA0*(2.*CA-CF)*ALF+3./2.*BETA0*(2.*CA-CF)*ALR
     &+12.*ALF*(4./9.*(T2+U2)/T/U*ALP+4.*((U2-T2)/T/U+2.*(U-T)/S)*ALUT)
     &/SIGB
      C2(1,7)=CA*(2.*CA*ALP+5.*BETA0/4.)*ALF**2
     &-3./2.*CA*BETA0*ALF*ALR
c  q g --> q g
      SIGB=2.-1./9.-4./9.*T2/S/U-2.*S*U/T2
      C2(4,8)=1./2.*(CF+CA)**2
      C2(3,8)=3./2.*(CF+CA)*(-(CF+CA)*ALF+CF*(-2.*ALU-3./4.)-2.*CA*ALT
     &-BETA0/3.)+3./2.*(CF+CA)*(-1./9.*(CF+CA)*(T2/S/U-2.)*ALT
     &+3.*(-1.-2.*S/T+U/2./S-S/2./U)*ALU+(CF*ALT+CA/2.*ALU)*((2./9.-1.)*
     &(T2/S/U-2.)+2.*(1.-2.*S*U/T2)))/SIGB
      C2(2,8)=(CF+CA)*(((4.*CF-CA)*ALU+(4.*CA-CF)*ALT+3./4.*CF+BETA0/4.
     &+(CF+CA)*ALF)*ALF+3./4.*BETA0*ALR)-2.*(CF+CA)*(-1./9.*(CF+CA)*
     &(T2/S/U-2.)*ALT+3.*(-1.-2.*S/T+U/2./S-S/2./U)*ALU+(CF*ALT
     &+CA/2.*ALU)*((2./9.-1.)*(T2/S/U-2.)+2.*(1.-2.*S*U/T2)))*ALF/SIGB
      C2(1,8)=(CF+CA)*((CF*ALT+CA*ALU+3./4.*CF+3./8.*BETA0)*ALF**2
     &-3./4.*BETA0*ALF*ALR)
c      c2(1,8)=0.
c      c2(2,8)=0.
c      c2(3,8)=0.
      SIGB=27/2.-9./2.*(S*U/T2+T*U/S2+S*T/U2)
      C2(4,9)=2.*CA**2
      C2(3,9)=3.*CA*(-2.*CA*ALF-2.*CA*ALP-7./12.*BETA0)
     &+3.*CA*(27./8.*(2.*ALT+5.*ALU)*(1.-T*U/S2-S*T/U2+T2/S/U)
     &-27./4.*ALU*(S*T/U2-T*U/S2+U2/S/T-S2/T/U)+(3.*ALT+3./2.*ALU)*
     &(27./4.-9.*(S*U/T2+T*U/4./S2+S*T/4./U2)+9./2.*(U2/S/T+S2/T/U
     &-T2/2./S/U)))/SIGB
      C2(2,9)=CA*(6.*CA*ALP+BETA0+4.*CA*ALF)*ALF+3./2.*CA*BETA0*ALR
     &-4.*CA*(27./8.*(2.*ALT+5.*ALU)*(1.-T*U/S2-S*T/U2+T2/S/U)
     &-27./4.*ALU*(S*T/U2-T*U/S2+U2/S/T-S2/T/U)+(3.*ALT+3./2.*ALU)
     &*(27./4.-9.*(S*U/T2+T*U/4./S2+S*T/4./U2)+9./2.*(U2/S/T+S2/T/U
     &-T2/2./S/U)))*ALF/SIGB
      C2(1,9)=CA*(2.*CA*ALP+5./4.*BETA0)*ALF**2-3.*BETA0/2.*CA*ALF*ALR
      RETURN
      END

