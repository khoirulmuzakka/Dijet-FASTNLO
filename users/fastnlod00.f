
* -------------------------------------------------------------------
* fastNLO user package                     M. Wobisch 06/29/2005
*
* fastnlod00.f
*
* contains the following routines
*     FNCALC          computes cross sections / fills arrays
*     FNREAD(ioord)   reads coefficient tables
*     FNOUT           prints x-section results on screen 
*     FNGETPDF        gets the PDF values for all x-bins
*
* uses commonblock definitions in: fastnlo00.inc
*
* needs the following routines from fn-pdf.f 
*     FNPAPDF         computes the 7 PDF linear combinations
*                         for "PA": proton-antiproton
*     FNPPPDF         computes the 7 PDF linear combinations
*                         for "PP": proton-proton
*     FNEPPDF         computes the 3 PDF linear combinations
*                         for "EP": electron-proton
*
* needs the following routines from "fnpdfalpha.f":
*     FNALPHAS        alpha_s interface (double precision function)
*     FNPDF           PDF interface
*
* question: change 2-d x-Array to 1-dim (save factor 2 in memory)?
*          - also store PDFs in 'parallel 1d array
* -------------------------------------------------------------------



*******************************************************************
      SUBROUTINE FNCALC(FILENAME, IMUFFLAG, XSECT)
*-----------------------------------------------------------------
* TK 12/07/2005 table format contains now LO + NLO with 5 scale variations, 
* support for zipped files
* MW 04/15/2005
*
* compute cross section
*
* input:
*   FILENAME     name of input table
*   IMUFFLAG     flag to enable/disable factorization scale variations
*                 0: off   
*                 1: compute x-sect for varied factoriation scales
*                    takes 3x time for additional PDF access
* output:
*   XSECT(nbin,2,nscalevar)    array of cross sections 
*                              first dimension: Bin number
*                              second dim: 1 LO   2 NLO
*                              third dim: 1 default   2....nscalevar varaiations
*-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fastnlod00.inc'
      INTEGER IFIRST, IFILE, iord,iscale, I,J,K,L,M, IMUFFLAG, nbin,nx
      CHARACTER*(*) FILENAME
      CHARACTER*50 OLDFILENAME
      DOUBLE PRECISION binwidth
      DATA IFIRST/0/, OLDFILENAME/'xxxx'/
      SAVE IFIRST, OLDFILENAME

c - output in first call
      if (IFIRST.EQ.0) then
         write(*,*) " "
         write(*,*) "  ##################################################",
     +        "#############"
         write(*,*) "  #" 
         write(*,*) "  #   fastNLO v0.0                       June 6, 2005"
         write(*,*) "  #" 
         write(*,*) "  #   Thomas Kluge, Klaus Rabbertz, Markus Wobisch" 
         write(*,*) "  #" 
         write(*,*) "  #   user code to compute the inclusive jet ",
     +        "cross section"
         write(*,*) "  #   for the H1 measurements at HERA "
         write(*,*) "  #      hep-ex/" 
         write(*,*) "  #      hep-ex/" 
         write(*,*) "  #  " 
         write(*,*) "  #  if you use this code, please cite the ",
     +        "following reference:" 
         write(*,*) "  #  ...." 
         write(*,*) "  #  the perturbative coefficients have been computed" 
         write(*,*) "  #  using NLOJET++ by Zoltan Nagy - please cite:" 
         write(*,*) "  #  ...." 
         write(*,*) "  #  " 
      endif

c --- read the fastNLO coefficient tables
      IF (FILENAME.ne.OLDFILENAME) THEN
c         write (*,*) ' old/new table ',OLDFILENAME,'   ',FILENAME
         call FNREAD(FILENAME)  ! read combined LO & NLO table
         OLDFILENAME=FILENAME
      ENDIF

      if (IFIRST.EQ.0) then
         IFIRST = 1
         write(*,*) "  #       -- coefficient tables are now in memory --"
         write(*,*) "  #################################################",
     +        "#############"
         write(*,*) " "
      endif


c - reset result arrays
      nbin = 0                  ! continuos numbering for the final array
      do i=1,nrapidity          ! (Pseudo-)Rapidity Bins
         do j=1,NPT(i)          ! ET/pT Bins
            nbin=nbin+1
            do k=1,NSCALEVAR  ! central + scale variations
               do m=1,4         ! 7 Subproc + total
                  result(nbin,m,1,k)=0d0 ! reset LO
                  result(nbin,m,2,k)=0d0 ! reset NLO
               enddo
               xsect(nbin,1,k)=0d0 ! reset LO
               xsect(nbin,2,k)=0d0 ! reset NLO
            enddo
         enddo
      enddo

c - compute xsect for default factorization scales
c   but for different renormalization scales
c      write(*,*) '>>>>> get PDFs'
      call FNGETPDF(1d0)
c      write(*,*) '>>>>>>>start 1st multiplication'
      call FNMULT(1)

c ----------------- final touches ------------------------------------
c - sum up subprocesses / divide by pT Bins width / fill result array
c -                   multiply by 1000 to get pb / fill array xsect
      nbin=0
      do i=1,nrapidity          ! Rapidity Bins
         do j=1,NPT(i)          ! pT Bins
            nbin=nbin+1
            binwidth=(PTBIN(i,j+1)-PTBIN(i,j))*(RAPBIN(i+1)-RAPBIN(i))

            do iscale=1,NSCALEVAR
               do iord=1,2
                  result(nbin,4,iord,iscale) = 0d0
                  xsect(nbin,iord,iscale) = 0d0
                  do m=1,3      ! Subprocess
                     result(nbin,m,iord,iscale) = result(nbin,m,iord,iscale)
     +                    *1000d0
                     result(nbin,4,iord,iscale) = result(nbin,4,iord,iscale)
     +                    + result(nbin,m,iord,iscale)
                  enddo
                  xsect(nbin,iord,iscale)=result(nbin,4,iord,iscale)
               enddo
            enddo
         enddo
      enddo

      RETURN
      END

*******************************************************************
      SUBROUTINE FNMULT(iscale)
*-----------------------------------------------------------------
* MW 08/26/2005
*
* multiply the PDFs and the coefficients for a single scale setting
* 
* input:    ISCALE = 1: compute for default ren and fact scales
*                  2: compute for increased ren scale (x2)
*                  3: compute for decreased ren scale (/2)
*                  4: compute for increased fact scale (x2)
*                  5: compute for decreased fact scale (/2)
*-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fastnlod00.inc'
      INTEGER ISCALE
      INTEGER i,j,k,l,m,iord,index,    nbin,nx
      DOUBLE PRECISION  mur, as, FNALPHAS

      nbin = 0                  ! continuos numbering for the final array
      do i=1,nrapidity          ! (Pseudo-)Rapidity Bins
         do j=1,NPT(i)          ! ET/pT Bins
            nbin=nbin+1            
            
            mur = murval(i,j) * murscale(iscale)
            as = FNALPHAS(mur)  ! get alpha_s

            nx=0
            do k=1,NTOT         ! Nxmax
               nx=nx+1
               do m=1,3         ! Subprocess
                  do iord=1,2   ! order: 1 LO  2 NLO
                     if (iord.eq.1) index=1
                     if (iord.eq.2) index=iscale+1
                     if (iord.eq.1) THEN
                     ENDIF
                     result(nbin,m,iord,iscale)=result(nbin,m,iord,iscale)
     +                    + array(nbin,nx,m,index)
     +                    * as**(iord) ! multiply alpha-s
     +                    * pdf(nbin,nx,m) ! multiply PDFs
                  enddo
               enddo
            enddo
         enddo
      enddo

      RETURN
      END

*******************************************************************
      SUBROUTINE FNGETPDF(muffactor)
*-----------------------------------------------------------------
* MW 08/26/2005
*
* read the PDFs in all (rap,pT) bins - at all xmax, xmin bins
* the default factorization scale (in GeV) is multiplied by muffactor
*-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fastnlod00.inc'
      DOUBLE PRECISION MUFFACTOR
      INTEGER NBIN,NX, i,j,k,l,m
      DOUBLE PRECISION x1, x2, xlim, hx, hxlim, muf, H(3)

      nbin=0
      do i=1,nrapidity          ! Rapidity Bins
         do j=1,NPT(i)          ! pT Bins
            nbin=nbin+1

            xlim = xlimit(i,j)
            hxlim = log10(xlim)
            muf = mufval(i,j) * muffactor
            
            nx=0 
            do k=1,NTOT         ! Nx
               hx = hxlim *(1d0 - dble(k-1)/dble(ntot)) ! compute x1-value  
               x1 = 10**hx
               nx=nx+1
               call FNEPPDF(x1,muf,H) ! get PDFs at x1 at scale muf
               do m=1,3        ! Subprocess
                  pdf(nbin,nx,m) = H(m)
c     +                 / (1d0-0.99d0*x1)**3/(1d0-0.99d0*x2)**3
c     +                 * sqrt(x1*x2) ! 2 lines PDF reweighting
               enddo
            enddo
         enddo
      enddo

      RETURN 
      END

*******************************************************************
      SUBROUTINE FNOUT
*-----------------------------------------------------------------
* MW 04/18/2005
*
* output cross section
*-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IFIRST, IFILE, I,J,K,L,M, NBIN, iscale
      INCLUDE 'fastnlod00.inc'

      write(*,*) "      FNOUT:   fastNLO - results:"
      write(*,*) "            cross sections:    LO       NLOcorr     tot"

c      do iscale=1,NSCALEVAR
      do iscale=1,1
         write(*,*) " "
         write(*,*) "  ----------- scale setting No.",iscale
         
         nbin = 0               ! linear bin for the final observable
         do i=1,nrapidity       ! Rapidity Bins
            WRITE(*,*) "   Q2 from ",RAPBIN(i)," - ",RAPBIN(i+1),
     +           "         No pT bins: ",NPT(i)
            do j=1,NPT(i)       ! pT Bins
               nbin=nbin+1
               WRITE(*,*) "  pT: ",PTBIN(i,j),"-",PTBIN(i,j+1),
     +              " sigma: ",real(result(nbin,4,1,iscale)),
     +              real(result(nbin,4,2,iscale)),
     +              real(result(nbin,4,1,iscale)+result(nbin,4,2,iscale))
            enddo
         enddo
      enddo
      RETURN
      END

*******************************************************************
      SUBROUTINE FNREAD(FILENAME)
*-----------------------------------------------------------------
* MW 04/15/2005
*
* read ASCII table
* iord: order    1:LO  2:NLO table
*-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) FILENAME
      CHARACTER*255 BUFFER
      INTEGER IFIRST, IFILE, I,J,K,L,M,N,   NBIN,NX
      INTEGER FLNLEN,LENOCC
      LOGICAL ZIPPED
      INCLUDE 'fastnlod00.inc'

      DATA IFIRST/0/
      SAVE IFIRST

      FLNLEN = LENOCC (FILENAME)
      IF (FILENAME(FLNLEN-2:FLNLEN) .EQ. '.gz')  THEN
         ZIPPED = .TRUE.
         write(*,*) '  #      reading zipped table' 
      ELSE
         ZIPPED = .FALSE.
         write(*,*) '  #      reading table' 
      ENDIF

      IF(ZIPPED) THEN
         CALL GZOPEN(2,'R' , FILENAME, IFILE)
      ELSE
      OPEN(2,STATUS='OLD',
     +        FILE=FILENAME,IOSTAT=IFILE)
      ENDIF

      IF (IFILE .ne. 0) THEN
         WRITE(*,*) ' fastNLO:  table file not found ',
     +        '  -  IOSTAT = ',IFILE
         STOP
      ENDIF

      WRITE(*,*) "  #       # max bins (eta,ET,x): ",NRAPMAX,NPTMAX,NXMAX
c     -----------------------------------
      IF(ZIPPED) THEN
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) IREACTION
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) ECMS
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) IPROC
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) IALGO
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) JETRES1
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) JETRES2
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) NPOW
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) OALPHAS
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) NEVT(1)
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) NEVT(2)
         write(*,*) '  #      order ',1,': ',NEVT(1),' events'
         write(*,*) '  #      order ',2,': ',NEVT(2),' events'
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) NTOT
         WRITE(*,*) "  #       # x bins in table: ",NTOT
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) IXSCHEME
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) IPDFWGT
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) NRAPIDITY
         do i=1,nrapidity+1
            CALL GZREAD(2,BUFFER,255, IFILE)
            READ(BUFFER,*) RAPBIN(i)
         enddo
         do i=1,nrapidity
            CALL GZREAD(2,BUFFER,255, IFILE)
            READ(BUFFER,*) NPT(i)
         enddo
         do i=1,nrapidity
            do j=1,NPT(i)+1
               CALL GZREAD(2,BUFFER,255, IFILE)
               READ(BUFFER,*) PTBIN(i,j)
            enddo
         enddo
         CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,nrapidity
            do j=1,NPT(i)
               CALL GZREAD(2,BUFFER,255, IFILE)
               READ(BUFFER,"(F20.17)") XLIMIT(i,j)
            enddo
         enddo
         CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,nrapidity
            do j=1,NPT(i)
               CALL GZREAD(2,BUFFER,255, IFILE)
               READ(BUFFER,*) MURVAL(i,j)
            enddo
         enddo
               CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,nrapidity
            do j=1,NPT(i)
               CALL GZREAD(2,BUFFER,255, IFILE)
               READ(BUFFER,*) MUFVAL(i,j)
            enddo
         enddo
         CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         CALL GZREAD(2,BUFFER,255,IFILE)
         READ(BUFFER,*) NSCALEVAR
         do i=1,NSCALEVAR
            CALL GZREAD(2,BUFFER,255, IFILE)
            READ(BUFFER,*) MURSCALE(i)
         enddo
         do i=1,NSCALEVAR
            CALL GZREAD(2,BUFFER,255, IFILE)
            READ(BUFFER,*) MUFSCALE(i)
         enddo

         CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         nbin=0
         do i=1,nrapidity
            do j=1,NPT(i)
               nbin=nbin+1      ! linear numbering of rap,pT bins
               nx=0
               do k=1,NTOT      ! Nxmax
                  nx=nx+1       ! linear numbering of xmax,xmin bins
                  do m=1,3      ! Subproc
                     do n=1,1+NSCALEVAR ! LO &  NLO with  scale variations
                        CALL GZREAD(2,BUFFER,255, IFILE)
                        READ(BUFFER,*) array(nbin,nx,m,n)
c     READ(BUFFER,*) array(i,j,k,l,m,n)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         CALL GZREAD(2,BUFFER,255, IFILE)
         READ(BUFFER,*) i
         if (i.ne.iseparator) goto 999
         CALL GZCLOSE(2,IFILE)
      ELSE
         READ(2,*) IREACTION
         READ(2,*) ECMS
         READ(2,*) IPROC
         READ(2,*) IALGO
         READ(2,*) JETRES1
         READ(2,*) JETRES2
         READ(2,*) NPOW
         READ(2,*) OALPHAS
         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         READ(2,*) NEVT(1)
         READ(2,*) NEVT(2)
         write(*,*) '  #      order ',1,': ',NEVT(1),' events'
         write(*,*) '  #      order ',2,': ',NEVT(2),' events'
         READ(2,*) NTOT
         WRITE(*,*) "  #       # x bins in table: ",NTOT
         READ(2,*) IXSCHEME
         READ(2,*) IPDFWGT
         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         READ(2,*) NRAPIDITY
         do i=1,nrapidity+1
            READ(2,*) RAPBIN(i)
         enddo
         do i=1,nrapidity
            READ(2,*) NPT(i)
         enddo
         do i=1,nrapidity
            do j=1,NPT(i)+1
               READ(2,*) PTBIN(i,j)
            enddo
         enddo
         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,nrapidity
            do j=1,NPT(i)
               READ(2,*) XLIMIT(i,j)
            enddo
         enddo
         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,nrapidity
            do j=1,NPT(i)
               READ(2,*) MURVAL(i,j)
            enddo
         enddo
         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         do i=1,nrapidity
            do j=1,NPT(i)
               READ(2,*) MUFVAL(i,j)
            enddo
         enddo
         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         READ(2,*) NSCALEVAR
         do i=1,NSCALEVAR
            READ(2,*) MURSCALE(i)
         enddo
         do i=1,NSCALEVAR
            READ(2,*) MUFSCALE(i)
         enddo

         READ(2,*) i
         if (i.ne.iseparator) goto 999
c     -----------------------------------
         nbin=0
         do i=1,nrapidity
            do j=1,NPT(i)
               nbin=nbin+1      ! linear numbering of rap,pT bins
               nx=0
               do k=1,NTOT      ! Nxmax
                  nx=nx+1       ! linear numbering of xmax,xmin bins 
                  do m=1,3      ! Subproc
                     do n=1,1+NSCALEVAR ! LO &  NLO with  scale variations
                        READ(2,*) array(nbin,nx,m,n)
c     READ(2,*) array(i,j,k,l,m,n)
                     enddo
                  enddo
               enddo
            enddo
         enddo

      READ(2,*) i
      if (i.ne.iseparator) goto 999
      CLOSE(2)
      ENDIF

      goto 998
      WRITE(*,*) "ntot  ",NTOT
      do i=1,nrapidity
         WRITE(*,*) i,"  ",NPT(i)
      enddo
      goto 998
      WRITE(*,*) IREACTION
      WRITE(*,*) ECMS
      WRITE(*,*) IPROC
      WRITE(*,*) IALGO
      WRITE(*,*) JETRES1, JETRES2
      WRITE(*,*) NPOW
      WRITE(*,*) OALPHAS
      WRITE(*,*) NEVT(1)
      WRITE(*,*) NEVT(2)
      WRITE(*,*) NTOT
      WRITE(*,*) IXSCHEME
      WRITE(*,*) IPDFWGT
      WRITE(*,*) NRAPIDITY
      do i=1,nrapidity+1
         WRITE(*,*) RAPBIN(i)
      enddo
      do i=1,nrapidity
         WRITE(*,*) i,"  ",NPT(i)
      enddo
      do i=1,nrapidity
         do j=1,NPT(i)+1
            WRITE(*,*) PTBIN(i,j)
         enddo
      enddo
      do i=1,nrapidity
         do j=1,NPT(i)
            WRITE(*,*) XLIMIT(i,j)
         enddo
      enddo



      WRITE(*,*) 
      WRITE(*,*) 

 998  continue
c      write(*,*) '  *   table is in memory'
      RETURN
 999  continue
      IF(ZIPPED) THEN
         CALL GZCLOSE(2,IFILE)
      ELSE
         close (2) 
      ENDIF
      write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      write(*,*) " >>>>>   fastNLO error in table format "
c      write(*,*) " >>>>>       separator at wrong place!"
      write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      RETURN
      END

*******************************************************************
 




