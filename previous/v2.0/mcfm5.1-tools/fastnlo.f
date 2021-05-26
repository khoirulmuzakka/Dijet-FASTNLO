
c ------------------------------------------------------------------------
c   here are the variables used in the fastNLO interface    MW 06/16/2007
c
c   fnsum(10)    for test purposes: array to sum up the x-sect
c   fnsigma      overall factor for the coefficients
c   fncoeff      2d array of coefficients
c   fncfv2       2d array of coefficients - for virtual only
c   fncfv3       2d array of coefficients - for virtual only
c   fnpdf        2d PDF array
c   fnpdf2       2d PDF array - for virtual only
c   fnpdf3       2d PDF array - for virtual only
c   fnx(2)       the x-values where PDFs are probed
c   fnxz(2)      the additional x-values where PDFs are probed in virtual
c ------------------------------------------------------------------------
      Double Precision fnsum(10), fnsigma,
     +     fncoeff(-6:6,-6:6),fncfv2(-6:6,-6:6),fncfv3(-6:6,-6:6),
     +     fnpdf(-6:6,-6:6),fnpdf2(-6:6,-6:6),fnpdf3(-6:6,-6:6), 
     +     fnx(2),fnxz(2),
     +     fncftmp(0:9,-6:6,-6:6),fntmpsum
      Integer ifnctrb           ! 1 lo 2 virt 3 real (or directly from MCFM

      Common /cfastnlo/ fnsum, fnsigma,
     +     fncoeff,fncfv2,fncfv3,
     +     fnpdf,fnpdf2,fnpdf3,
     +     fnx,fnxz,
     +     ifnctrb


c - alphas
c   src/Inc/qcdcouple.f:      double precision gsq,as,ason2pi,ason4pi
c   src/Inc/qcdcouple.f:      common/qcdcouple/gsq,as,ason2pi,ason4pi
c .      fac=gw**4*xn
c +g     fac=gwsq**2*gsq*V


