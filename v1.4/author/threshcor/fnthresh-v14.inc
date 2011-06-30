c
c     parameters for fastNLO table 
c
      INTEGER NordMax, Nrapmax, Nptmax,Nscalemax,Nrapptmax,Nxmax,
     +      NSCALEBINMAX
      PARAMETER (NordMax=10, Nrapmax=20, Nptmax=160, 
     +      Nscalemax=6, Nrapptmax=700, Nxmax=325, NSCALEBINMAX=2)
      Double Precision mu0scale
      Parameter (mu0scale=0.25) ! reference scale for scale interpolation

      INTEGER Iseparator, ITABVERSION, Ireaction,IXSECTUNITS, Iproc, Ialgo, 
     +      Nord, Npow(NordMax),
     +      Nxtot, Ixscheme, Ipdfwgt, IREF, NBINTOT,NDIMENSION,
     +      Nrap, Npt(Nrapmax),nscalevar,NSCALEBIN

      DOUBLE PRECISION ECMS, JetResol1, JetResol2,
     +      Nevt(NordMax), yarr(Nrapmax+1),
     +      ptarr(Nrapmax+1,Nptmax+1), 
     +	    Xlimit(Nrapmax,Nptmax),Hxlim(Nrapmax,Nptmax),
     +      Murval(Nrapmax,Nptmax,NSCALEBINMAX),
     +      Mufval(Nrapmax,Nptmax,NSCALEBINMAX),
     +      Murscale(Nscalemax),Mufscale(Nscalemax), 
     +      Murvaltrans(Nrapmax,Nptmax,Nscalemax,NSCALEBINMAX),
     +      Mufvaltrans(Nrapmax,Nptmax,Nscalemax,NSCALEBINMAX),
     +      Nevents

      CHARACTER*64 NAMELABEL(5)
      CHARACTER*64 POWLABEL(6)
      CHARACTER*64 CODELABEL(6)
      CHARACTER*64 DIMLABEL(2)
      CHARACTER*64 SCALELABEL


c - question -> use 1d x binning??   (smaller but harder to fill)
c      DOUBLE PRECISION table(Nrapptmax,Nscalemax,Nxmax,7)
      DOUBLE PRECISION table(0:2,Nrapptmax,25,25,Nscalemax,7,NSCALEBINMAX)

      COMMON /cthrcr1/
     +      ECMS, JetResol1, JetResol2, yarr, ptarr,
     +      Nevt, Nevents, Xlimit, Murval, Mufval, Murscale, Mufscale,
     +      Murvaltrans,Mufvaltrans,
     +      Ireaction,Iproc, Ialgo, Nord, Npow, 
     +      NAMELABEL,POWLABEL,CODELABEL,DIMLABEL,SCALELABEL, 
     +      Nxtot, Ixscheme, Ipdfwgt, IREF,NBINTOT,NDIMENSION,
     +      Nrap, Npt, nscalevar,NSCALEBIN

      COMMON /cthrcr2/ table


      DATA Iseparator/1234567890/