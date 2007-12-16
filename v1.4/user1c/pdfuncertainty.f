      PROGRAM PDFUNC
* -------------------------------------------------------------------
* M. Wobisch 02/09/2006
*
* fastNLO - example program to compute PDF uncertainties
*           using PDFs from LHAPDF
*
* -------------------------------------------------------------------
      implicit none
      CHARACTER*255 FILENAME
      CHARACTER*255 HISTOFILE
      CHARACTER*255 LHAPDF
      integer i,j,k,l1,l2,l3,l4, npdf
      integer lenocc
      INCLUDE 'fnx9999.inc'
      double precision mur,muf, diff,
     +     res0(NBINTOTMAX,NMAXSUBPROC+1,3),
     +     res1hi(NBINTOTMAX,NMAXSUBPROC+1,3),
     +     res1lo(NBINTOTMAX,NMAXSUBPROC+1,3)

c - Attention!!! - this mus be declared consistent with the
c                  definition in the commonblock!!!!!
      double precision xsect1(900,3),xsect0(900,3)


c --- parse command line
      IF ( IARGC().LT.1)  FILENAME = 'table.txt.gz'
      IF ( IARGC().LT.2)  HISTOFILE= 'fastnlo.hbk'
      IF ( IARGC().GT.0)  CALL GETARG(1,FILENAME)
      IF ( IARGC().GT.1)  CALL GETARG(2,HISTOFILE)
      IF ( IARGC().GT.2)  THEN
         write(*,*) 'fastNLO: too many arguments given. Stopping'
         RETURN
      ENDIF


c - Initialize LHAPDF
      CALL GETENV('LHAPDF',LHAPDF)
      IF (LENOCC(LHAPDF).EQ.0) THEN
         LHAPDF = '/disk2/work/wobisch/lhapdf-4.1/PDFsets/cteq61.LHgrid'
      ENDIF
      write(*,*)"Looking for LHAPDF in directory "//
     &     LHAPDF(1:LENOCC(LHAPDF))//"!"
ckr      call InitPDFset(LHAPDF(1:LENOCC(LHAPDF))//'/PDFsets/cteq61.LHgrid')
      call InitPDFset(LHAPDF(1:LENOCC(LHAPDF))//'/../PDFsets/cteq65.LHgrid')

c      call InitPDFset('/work/shootingstar-clued0/wobisch/lhapdf500/share/lhapdf/PDFsets/cteq61.LHgrid')
c      call InitPDFset('/disk2/work/wobisch/lhapdf-4.2/PDFsets/cteq61.LHgrid')
c      call InitPDFset('/disk2/work/wobisch/lhapdf-4.2/PDFsets/MRST2001E.LHgrid')
c      call InitPDFset('/disk2/work/wobisch/lhapdf-4.2/PDFsets/MRST2004nnlo.LHgrid')
c      call InitPDFset('/disk2/work/wobisch/lhapdf-4.2/PDFsets/a02m_nlo.LHgrid')


c      call InitPDFset('/h1/h1gen/lhapdf/LHAPDFv4/PDFsets/cteq61.LHgrid')
      call numberPDF(NPDF)
      WRITE(*,*) "fastNLO:   the PDF set has ",NPDF," members"

c - one initial call - to fill commonblock -> for histo-booking
      call FX9999CC(FILENAME, 1d0 , 1d0 , 0 , XSECT1)
      call PDFHIST(1,histofile)

c- new call: a single call for each scale
c         1st argument:  name of table
c         2nd argument:  xmur  prefactor for nominal ren-scale
c                              any choice is possible, but please note 
c                              that NNLO-NLL works only for xmur=xmuf
c         3rd argument:  xmuf  prefactor for nominal fact-scale
c                              only a few choices are possible
c                              (see output or table documentation)
c         4th argument:  0: no ascii output       1: print results
c         5th argument:  array to return results

c -> compute PDF uncertainties for all available scales
      do i=1,Nscalevar
         write(*,*) " now scale No.",i
         call InitPDF(0)
         mur = murscale(i)
         muf = mufscale(i)
         call FX9999CC(FILENAME, mur , muf, 0 , XSECT1)

c -    save the result array from the first call (= central result)
c      and reset result arrays   
         write(*,*) " the observable has ",NBINTOT," bins  - ",NSUBPROC," subprocesses"
         do l1=1,NBINTOT
            do l2=1,(NSUBPROC+1)
               do l3=1,NORD
                  res0(l1,l2,l3) = 0d0 
                  do l4=1,l3
                     res0(l1,l2,l3) = res0(l1,l2,l3)+result(l1,l2,l4)
                  enddo
                  res1lo(l1,l2,l3) = 0d0
                  res1hi(l1,l2,l3) = 0d0
               enddo
            enddo
         enddo

         do j=1,NPDF
            call InitPDF(j)
            call FX9999CC(FILENAME, mur , muf, 0 , XSECT0)

c - for all bins/subproc/orders: add negative/positive variations
            do l1=1,NBINTOT
               do l2=1,(NSUBPROC+1)
                  do l3=1,NORD
                     diff = - res0(l1,l2,l3)
                     do l4=1,l3
                        diff = diff + result(l1,l2,l4) 
                     enddo
                     if (diff .gt. 0d0) then
                        res1hi(l1,l2,l3) = res1hi(l1,l2,l3)+diff*diff
                     else
                        res1lo(l1,l2,l3) = res1lo(l1,l2,l3)+diff*diff
                     endif
                  enddo
               enddo
            enddo
         enddo                  ! loop over bins

c - take square-root of sum of squares
         do l1=1,NBINTOT
            do l2=1,(NSUBPROC+1)
               do l3=1,NORD
                  res1hi(l1,l2,l3) = sqrt(res1hi(l1,l2,l3))
                  res1lo(l1,l2,l3) = -sqrt(res1lo(l1,l2,l3))
               enddo
            enddo
            WRITE(*,*) l1,res0(l1,NSUBPROC+1,NORD),
     +           res1lo(l1,NSUBPROC+1,NORD)/res0(l1,NSUBPROC+1,NORD),
     +           res1hi(l1,NSUBPROC+1,NORD)/res0(l1,NSUBPROC+1,NORD)
         enddo

c - fill histograms
         call PDFFILL(i,res0,res1hi,res1lo)

      enddo                     ! loop over scales
      write(*,*) 'Bin    x-sect       lower PDF    upper PDF unc.'
      write(*,*) '(the printed uncertainties are for the highest order)'
      write(*,*) '(histograms contain results for all orders and subprocesses)'

c - close hbook file
      call PDFHIST(2,histofile)
      RETURN
      END

c
c
c ======================= do the histogramming =========================
      SUBROUTINE PDFHIST(n,histofile)
      implicit none
      Character*(*) histofile
      INTEGER IFIRST, IFILE, J,N,istat2,icycle
      Integer ix,iord,isub, iscale, irap, ihist
      INCLUDE 'fnx9999.inc'
      real pt(nptmax)

c - HBOOK common 
      INTEGER NWPAWC
      PARAMETER (NWPAWC=2500000)
      REAL HMEMOR(NWPAWC)
      COMMON /PAWC/ HMEMOR


c - open & book
      if (n.eq.1) then
         write(*,*)' ----------- book histograms -------'
         CALL HLIMIT(NWPAWC)
         CALL HROPEN(11,'fastNLO',histofile,'N',1024,ISTAT2)
         if (ISTAT2.NE.0) then
            WRITE(*,*) ' FNHBOOK: could not open histofile ',istat2
         endif
         
         do iord=0,Nord         ! order: tot, LO, NLO-corr, NNLO-corr -> Nord
            do iscale=1,NSCALEVAR ! scale variations
               do isub=0,Nsubproc      ! subprocess: 0 tot + 7 subproc
               
                  do irap=1, nrapidity
                     ihist = iord*1000000+iscale*100000+isub*10000+irap*100
                     do j=1,(npt(irap)+1)
                        pt(j) = real(PTBIN(irap,j))
                     enddo
                     call hbookb(ihist,'p?T! (GeV)' , NPT(irap) ,PT ,0)
                     call hbookb(ihist+1,'p?T! (GeV)' , NPT(irap) ,PT ,0)
                     call hbookb(ihist+2,'p?T! (GeV)' , NPT(irap) ,PT ,0)
                  enddo
               enddo
            enddo
         enddo        


c - close HBOOK file
      elseif (n.eq.2) then
         CALL HROUT (0,ICYCLE,' ')
         CALL HREND ('fastNLO')
      endif

      RETURN
      END
c ======================= do the histogramming =========================
      SUBROUTINE PDFFILL(nscale,res0,res1hi,res1lo)
      IMPLICIT NONE
      INTEGER nscale,   I,J,K,L, nbin,nx
      Integer iord,isub,isub2, iscale,ihist
      INCLUDE 'fnx9999.inc'
      double precision 
     +     res0(NBINTOTMAX,NMAXSUBPROC+1,3),
     +     res1hi(NBINTOTMAX,NMAXSUBPROC+1,3),
     +     res1lo(NBINTOTMAX,NMAXSUBPROC+1,3)
      real val0,vallo,valhi,v

      if (nscale.lt.1 .or. nscale.gt.Nscalevar) then
         write(*,*) " FNHFILL:   nscale=",nscale," is out of range"
         write(*,*) "            max: ",nscalevar
         stop
      endif
      iscale = nscale

c - fill all histograms for the given scale
      do iord=0,Nord            ! order: tot, LO, NLO-corr, 3 NNLOcorr
         do isub2=1,(Nsubproc+1) ! subprocess: Nsubproc + 1 tot
            isub=isub2
            if (isub.eq.8) isub=0
            nbin=0
            do i=1,nrapidity                   
               do j=1,npt(i)
                  nbin = nbin + 1

                  if (iord.gt.0) then
                     val0  = real(res0(nbin,isub2,iord))
                     vallo = real(res1lo(nbin,isub2,iord))
                     valhi = real(res1hi(nbin,isub2,iord))
                  else
                     val0  = real(res0(nbin,isub2,Nord))
                     vallo = real(res1lo(nbin,isub2,Nord))
                     valhi = real(res1hi(nbin,isub2,Nord))
                  endif

                  ihist = iord*1000000+iscale*100000+isub*10000+i*100
                  call hfill(ihist, real(PTBIN(i,j)+0.01) ,0.0, val0)
                  call hfill(ihist+1, real(PTBIN(i,j)+0.01) ,0.0, vallo)
                  call hfill(ihist+2, real(PTBIN(i,j)+0.01) ,0.0, valhi)
               enddo            ! pT-loop
            enddo               ! rap-loop
         enddo                  ! isub loop
      enddo                     ! iord-loop

      RETURN
      END
