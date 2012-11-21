      subroutine diffpdf(xpom,zpom,muf,pdfs)
      implicit none
      double precision beta,muf,pt2max,pt2min,mypdf,xpom,pdfs(-6:6)
      double precision zpom
      double precision srh

      external mypdf
      
c     set consistent names
      beta = zpom

c      print * ,"zpom ", zpom
c      print * ,"muf ", muf
c      print * ,"xpom ", xpom
c      print * , "pdf-6 " , pdfs(-6)
c      print * , "pdf 0 ", pdfs(0)
c      print * , "pdf 1 ", pdfs(1)
c      print * , "pdf 6 ", pdfs(6)



c      pdfs(-6)=-3006
c      pdfs(-2)=-3002
c      pdfs(0)=-3000
c      pdfs(2)=3002
c      pdfs(6)=3006
      
c     The scale , for example, is set to the minimum accessible 
c     in VFPS analysis. 
      
c      Q2min=5d0
c      pt2max=5.5d0
c      pt2min=4.0d0
      
c     Q=(0.5d0*(pt2max+pt2min)+Q2min)**0.5
      
c      do beta=0.01d0,0.99d0,0.01d0
         
c     the bq(-nf:nf) array contains beta*distribution 



      pdfs(-6)= 0d0
      pdfs(-5)= 0d0	
      pdfs(-4)= beta*(mypdf(-4, beta, muf, xpom))
      pdfs(-3)= beta*(mypdf(-3, beta, muf, xpom))	
      pdfs(-2)= beta*(mypdf(-2, beta, muf, xpom))
      pdfs(-1)= beta*(mypdf(-1, beta, muf, xpom))	
      pdfs(0) = beta*(mypdf( 0, beta, muf, xpom))
      pdfs(1) = beta*(mypdf( 1, beta, muf, xpom))
      pdfs(2) = beta*(mypdf( 2, beta, muf, xpom))
      pdfs(3) = beta*(mypdf( 3, beta, muf, xpom))
      pdfs(4) = beta*(mypdf( 4, beta, muf, xpom))
      pdfs(5) = 0d0
      pdfs(6) = 0d0

c      print * , "dpdf(-6) ", pdfs(-6)
c      print * , "dpdf(-2) ", pdfs(-2)
c      print * , "dpdf(0) ", pdfs(0)
c      print * , "dpdf(2) ", pdfs(2)
c      print * , "dpdf(6) ", pdfs(6)
      
c      enddo
   
      return
      end

      INCLUDE 'qcd_2006.f'
      INCLUDE 'h12006flux.f'
      INCLUDE 'i_2006_fita.f'       
      INCLUDE 'i_2006_fitb.f' 	
      INCLUDE 'pion_stand_alone.f'
      
c===================================================================
      
        Double precision function mypdf (Iparton, X, Q, xpom)
        implicit none
	integer ifit,int,ipom,Iparton,calling
	double precision xpom,t,flux,qq2,q,x,dxpom,my,xpbin
	double precision XPQ(-6:6),F2(1:2),FL(1:2),C2(1:2),CL(1:2)
        double precision xpmin,xpmax,Pflux,Rflux
	double precision UPV,DNV,SEA,STR,CHM,GL
	external STROWP1

c	ifit=1		! Fit type A
	ifit=2		! Fit type B
	t=-1d0		! maximum t
	int=1		! t-integrated flux [t..tmin]

c	Get Pomeron flux in the proton. It depends on xpom & t but not on beta and Q2

	ipom=1
	call h12006flux(xpom,t,int,ifit,ipom,flux)
	Pflux=flux	    

c	Get Reggeon flux in the proton. It depends on xpom & t but not on beta and Q2 

	ipom=2
	call h12006flux(xpom,t,int,ifit,ipom,flux)
	Rflux=flux	    

	qq2=q*q

c 	call stand-alone pion structure function by Owen
c       pion pdf's are used for the reggeon term, best we can do 

        CALL STROWP1(X,Q,UPV,DNV,SEA,STR,CHM,GL)

c	get pomeron PDF's

	CALL QCD_2006(x,qq2,IFIT,XPQ,F2,FL,C2,CL)

c	various contributions are added together 
c       final DPDF(beta) = pomeron_pdf * Pomflux + reggeon_pdf * Regflux

c	mypdf return parton densities....so XPQ are devided by x (which is beta, in fact)

	if((Iparton.eq.0))then
	mypdf=XPQ(Iparton)/x*Pflux+(GL)/x*RFLUX	
	endif

	if((Iparton.eq.1))then
	mypdf=XPQ(Iparton)/x*Pflux+(UPV+SEA)/x*RFLUX	
	endif

	if((Iparton.eq.-1))then
	mypdf=XPQ(Iparton)/x*Pflux+(SEA)/x*RFLUX	
	endif

	if((Iparton.eq.2))then
	mypdf=XPQ(Iparton)/x*Pflux+(DNV+SEA)/x*RFLUX	
	endif

	if((Iparton.eq.-2))then
	mypdf=XPQ(Iparton)/x*Pflux+(SEA)/x*RFLUX	
	endif

	if((Iparton.eq.-3).or.(Iparton.eq.3))then
	mypdf=XPQ(Iparton)/x*Pflux+(STR)/x*RFLUX	
	endif

c       There is a scheme inconsistency between the massless nf=4 
c	scheme employed by NLOjet++ and the FFNS nf=3 of FitB DPDF's

c	Charm contribution are obtained from charm structure functions F2c, 
c       here is called C2(1), and returned by H1 fit 2006 
c       ( which is performed in FFNS scheme, nf=3)

	if((Iparton.eq.-4).or.(Iparton.eq.4))then
        mypdf=9.d0/8.d0*C2(1)/x*Pflux+(CHM)/x*RFLUX
	endif

c	bottom and top set to zero
	
	if((Iparton.eq.5).or.(Iparton.eq.6))then
	mypdf=0d0
	endif
	
	if((Iparton.eq.-5).or.(Iparton.eq.-6))then
        mypdf=0d0
	endif   

	return
	end

