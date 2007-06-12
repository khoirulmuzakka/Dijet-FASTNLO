c
c     parameters for fastNLO table 
c
      INTEGER NordMax, Nrapmax, Nptmax,Nscalemax,Nrapptmax,Nxmax
      PARAMETER (NordMax=10, Nrapmax=12, Nptmax=100, 
     +      Nscalemax=4, Nrapptmax=600, Nxmax=325)

      INTEGER Iseparator, Ireaction,Iproc, Ialgo, Nord, Npow(NordMax),
     +      Nxtot, Ixscheme, Ipdfwgt, 
     +      Nrap, Npt(Nrapmax),nscalevar

      DOUBLE PRECISION ECMS, JetResol1, JetResol2,
     +      Nevt(NordMax), yarr(Nrapmax+1),
     +      ptarr(Nrapmax+1,Nptmax+1), 
     +	    Xlimit(Nrapmax,Nptmax),Hxlim(Nrapmax,Nptmax),
     +      Murval(Nrapmax,Nptmax),Mufval(Nrapmax,Nptmax),
     +      Murscale(Nscalemax),Mufscale(Nscalemax), Nevents

c - question -> use 1d x binning??   (smaller but harder to fill)
c      DOUBLE PRECISION table(Nrapptmax,Nscalemax,Nxmax,7)
      DOUBLE PRECISION table(Nrapptmax,25,25,Nscalemax,7)

      COMMON /cthrcr1/
     +      ECMS, JetResol1, JetResol2, yarr, ptarr,
     +      Nevt, Nevents, Xlimit, Murval, Mufval, Murscale, Mufscale,
     +      Ireaction,Iproc, Ialgo, Nord, Npow,
     +      Nxtot, Ixscheme, Ipdfwgt, 
     +      Nrap, Npt, nscalevar

      COMMON /cthrcr2/ table



      DATA Iseparator/1234567890/