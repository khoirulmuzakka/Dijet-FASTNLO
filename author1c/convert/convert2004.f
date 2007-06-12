      PROGRAM CONVERT
c ---------------------------------------------------------------------
c --- M. Wobisch 2006/05/25
c --- 
c --- CONVERT 2003+ (D=0.7,-0.5,1.0) TO 2004/2005/2006
c --- 
c ---------------------------------------------------------------------
      implicit none
      Integer  i,j,k,l,m, irap,ipt, icor, ipdf
      Double Precision kd,ld
      Integer  Iset, npoint, nypt, nptpt,nxpt,ns4pt, iscale, nbin
      Double Precision  s,t,u, xamin,xa,xabin, s4max,s4,s4bin,
     +     xb,sh,th,uh, xb0,sh0,th0,uh0, 
     +     mur,mur2, muf,muf2,
     +     rs,y,ymin,ymax,ytest,yup,xt,pt,pt2,ptmin,ptmax
      Double Precision fac,fac1,fac2,  as,alpha_s, pi
      Double Precision temp,A1,A2, als4m, als4
      Double Precision CF(9) ! coefficients removed from routine COEF
      Double Precision norm, 
     +     xsa(9),xsa1(9),xsa2(9),
     +     xsb(9),xsb1(9),xsb2(9),
     +     xsc(9),xsc1(9),xsc2(9)
      CHARACTER*100 INTABLE, OUTTAB


      write(*,*) ' '
      write(*,*) ' ***********************************************'
      write(*,*) ' *     convert table 2003+ to 2004/2005/2006'
      write(*,*) ' * '
      write(*,*) ' * Markus Wobisch May 25, 2006'
      write(*,*) ' * '
      write(*,*) ' ***********************************************'



      IF ( IARGC().LT.1)  INTABLE = 'table.tab'
      IF ( IARGC().LT.2)  OUTTAB= 'result.txt'
      IF ( IARGC().GT.0)  CALL GETARG(1,INTABLE)
      IF ( IARGC().GT.1)  CALL GETARG(2,OUTTAB)


      write(*,*) ' read  ',intable
      write(*,*) ' write ',OUTTAB


c - read fastNLO table 
      CALL READTABLE(intable)
      CALL WRITETABLE(OUTTAB)

      RETURN
      END


C ----------------------------------------------------------
      SUBROUTINE READTABLE(filename)
c -------------------------------------------
c M. Wobisch Dec 02, 2005
c read fastNLO table 
c
c -------------------------------------------
      implicit none
      CHARACTER*(*) filename
c ???
c      CHARACTER*255 BUFFER
      INTEGER IFIRST, IFILE, I,J,K,L,M,N,   NBIN,NX
      include 'convert2003.inc'



      write (*,*) 'read fastNLO table ',filename
      open(unit=2,file=filename,status='unknown')


      READ(2,*) i
      if (i.ne.iseparator) goto 999
      READ(2,*) ITABVERSION
      WRITE(*,*) "#       tableformat is version",real(itabversion)/10000d0
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c  -----------------------------------
      READ(2,*) IREACTION
      READ(2,*) ECMS
      READ(2,*) IXSECTUNITS
      do i=1,5
         READ(2,*) NAMELABEL(i)
      enddo
      READ(2,*) IPROC
      READ(2,*) IALGO
      READ(2,*) JETRES1
      READ(2,*) JETRES2
      READ(2,*) NORD
      do i=1,nord
         READ(2,*) NPOW(i)
      enddo
      do i=1,nord
         READ(2,*) POWLABEL(i)
         READ(2,*) CODELABEL(i)
      enddo
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      do i=1,nord
         READ(2,*) NEVT(i)
      enddo
      do i=1,nord
         write(*,5000) ' #      ',NEVT(i),
     +        ' events in ',POWLABEL(i)
      enddo
      READ(2,*) NXTOT
      WRITE(*,*) "#           No. of x bins: ",NXTOT
      READ(2,*) IXSCHEME
      if(IXSCHEME.NE.2) THEN
         write(*,*) 'Up to now only ixscheme=2 implemented. Stop'
         stop
      ENDIF
      READ(2,*) IPDFWGT
      READ(2,*) IREF
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      READ(2,*) NBINTOT
      READ(2,*) NDIMENSION
      do i=1,ndimension
         READ(2,*) DIMLABEL(i)
      enddo
      
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
c   -----------------------------------
      do i=1,nrapidity
         do j=1,NPT(i)
            READ(2,*) XLIMIT(i,j)
         enddo
      enddo
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      READ(2,*) SCALELABEL
      READ(2,*) NSCALEBIN
      do i=1,nrapidity
         do j=1,NPT(i)
            do k=1,NSCALEBIN
               READ(2,*) MURVAL(i,j,k)
            enddo
         enddo
      enddo
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      do i=1,nrapidity
         do j=1,NPT(i)
            do k=1,NSCALEBIN
               READ(2,*) MUFVAL(i,j,k)
            enddo
         enddo
      enddo
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      READ(2,*) NSCALEVAR
      do i=1,NSCALEVAR
         READ(2,*) MURSCALE(i)
      enddo
      do i=1,NSCALEVAR
         READ(2,*) MUFSCALE(i)
      enddo
      
      READ(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------

      if((IREACTION.eq.2).OR.(IREACTION.eq.3)) then ! pp or ppbar
         NXSUM = (NXTOT*NXTOT+NXTOT)/2
         Nsubproc = 7
      elseif(IREACTION.eq.1) then ! DIS
         NXSUM = NXTOT
         Nsubproc = 3
      endif

      nbin=0
      do i=1,nrapidity
         do j=1,NPT(i)
            nbin=nbin+1         ! linear numbering of rap,pT bins
            do k=1,NXSUM        ! tot. No. x-Bins
               do m=1,Nsubproc  ! No. of Subproc
                  do n=1,1+NSCALEVAR*(NORD-1) ! LO & NLO & w/ scale var
                     do l=1,NSCALEBIN ! No. of Subproc
                        READ(2,*) array(nbin,k,m,n,l)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo


      READ(2,*) i
      if (i.ne.iseparator) goto 999
      CLOSE(2)

      write(*,*)  " "
      write(*,*)  "      >>>>>>>>>>> table is in memory <<<<<<<<<<<<<<"
      write(*,*)  " "

      RETURN
 999  continue

      close (2) 
      write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      write(*,*) " >>>>>   fastNLO error in table format "
      write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      stop
      RETURN

 5000 FORMAT (A,I12,A,A64) 

      RETURN
      END
C --------------------------------------------------------------------
      SUBROUTINE WRITETABLE(filename)
c -------------------------------------------
c M. Wobisch Feb 22, 2006
c write fastNLO raw-table - in text format
c
c  code specific for FNT2003->FNT2004
c
c -------------------------------------------
      implicit none
      CHARACTER*(*) filename
      CHARACTER*(33) label(0:2)
      integer i,j,k,l,m,n, isub, nbin, mord,mpow,mscalevar
      Integer Nrapmin, NPTMIN(20), NBINOFFSET, NSKIP(20), nbinnew
      include 'convert2003.inc'

      WRITE(*,*) ' write: ', filename

c - define variables -> to selct the bins that we want to keep
      do i=1,20
         NPTMIN(i) = 1
      enddo
      Nrapmin = 1
c      Nbinoffset = 0           ! No of bins that are skipped at beginning

c      Nbinoffset = 7           ! No of bins that are skipped at beginning
      Nrapmin = 1
      nrapidity = 5
      NPTMIN(1) = 8
      NPT(1) = 17
      NPTMIN(2) = 8
      NPT(2) = 17
      NPTMIN(3) = 8
      NPT(3) = 16
      NPTMIN(4) = 8
      NPT(4) = 14
      NPTMIN(5) = 8
      NPT(5) = 12

      NSKIP(1) = 7 ! No skipped bins before 1st rapidity range (was Nbinoffset)
      NSKIP(2) = 8
      NSKIP(3) = 8
      NSKIP(4) = 8
      NSKIP(5) = 8

c  - update NBINTOT info
      NBINTOT=0
      do i=Nrapmin,Nrapmin-1+nrapidity
         do j=NPTMIN(i),NPTMIN(i)-1+NPT(i)
            NBINTOT=NBINTOT+1
         enddo
      enddo
      WRITE(*,*) " "
      WRITE(*,*) "    new tot No bins ",NBINTOT
      WRITE(*,*) " "




      open(unit=2,file=filename,status='unknown')
      WRITE(2,5001) Iseparator     ! ------------------------
      WRITE(2,*) ITABVERSION
      WRITE(2,5001) Iseparator     ! ------------------------
      WRITE(2,*) IREACTION
      WRITE(2,*) ECMS
      WRITE(2,*) IXSECTUNITS
      do i=1,5
         WRITE(2,5009) NAMELABEL(i)
      enddo
      WRITE(2,*) IPROC
      WRITE(2,*) IALGO
      WRITE(2,*) JETRES1
      WRITE(2,*) JETRES2
      WRITE(2,*) NORD
      do i=1,nord
         WRITE(2,*) NPOW(i)
      enddo
      do i=1,nord
         WRITE(2,5009) POWLABEL(i)
         WRITE(2,5009) CODELABEL(i)
      enddo
      WRITE(2,5001) Iseparator     ! ------------------------
      do i=1,nord
         WRITE(2,*) NEVT(i)
      enddo
      WRITE(2,*) NXTOT
      WRITE(2,*) IXSCHEME
      WRITE(2,*) IPDFWGT
      WRITE(2,*) IREF
      WRITE(2,5001) Iseparator     ! ------------------------
      WRITE(2,*) NBINTOT
      WRITE(2,*) NDIMENSION
      do i=1,ndimension
         WRITE(2,5009) DIMLABEL(i)
      enddo     
      WRITE(2,*) NRAPIDITY
      do i=Nrapmin,Nrapmin-1+nrapidity+1
         WRITE(2,*) RAPBIN(i)
      enddo
      do i=Nrapmin,Nrapmin-1+nrapidity
         WRITE(2,*) NPT(i)
      enddo
      do i=Nrapmin,Nrapmin-1+nrapidity
         do j=NPTMIN(i),NPTMIN(i)-1+NPT(i)+1
            WRITE(2,*) PTBIN(i,j)
         enddo
      enddo
      WRITE(2,5001) Iseparator     ! ------------------------
      do i=Nrapmin,Nrapmin-1+nrapidity
         do j=NPTMIN(i),NPTMIN(i)-1+NPT(i)
            WRITE(2,5002) XLIMIT(i,j)
         enddo
      enddo
      WRITE(2,5001) Iseparator     ! ------------------------
      WRITE(2,5009) SCALELABEL
      WRITE(2,*) NSCALEBIN
      do i=Nrapmin,Nrapmin-1+nrapidity
         do j=NPTMIN(i),NPTMIN(i)-1+NPT(i)
            do k=1,NSCALEBIN
               WRITE(2,*) MURVAL(i,j,k)
            enddo
         enddo
      enddo
      WRITE(2,5001) Iseparator     ! ------------------------
      do i=Nrapmin,Nrapmin-1+nrapidity
         do j=NPTMIN(i),NPTMIN(i)-1+NPT(i)
            do k=1,NSCALEBIN
               WRITE(2,*) MUFVAL(i,j,k)
            enddo
         enddo
      enddo
      WRITE(2,5001) Iseparator     ! ------------------------
      WRITE(2,*) NSCALEVAR
      do i=1,NSCALEVAR
         WRITE(2,*) MURSCALE(i)
      enddo
      do i=1,NSCALEVAR
         WRITE(2,*) MUFSCALE(i)
      enddo
      WRITE(2,5001) Iseparator     ! ------------------------

c      nbin = Nbinoffset
      nbin = 0
      nbinnew = 0
      do i=Nrapmin,Nrapmin-1+nrapidity
         nbin = nbin + nskip(i)
         do j=NPTMIN(i),NPTMIN(i)-1+NPT(i)
            nbin=nbin+1         ! linear numbering of rap,pT bins
            nbinnew = nbinnew+1
            write(*,*)' >> ',nbin,nbinnew
            do k=1,NXSUM        ! tot. No. x-Bins
               do m=1,Nsubproc  ! No. of Subproc
                  do n=1,1+NSCALEVAR*(NORD-1) ! LO & NLO & w/ scale var
                     do l=1,NSCALEBIN ! No. of Subproc
                        if (array(nbin,k,m,n,l).eq.0d0) then
                           WRITE(2,5009) '0'
                        else
                           WRITE(2,5003) array(nbin,k,m,n,l)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      WRITE(2,5001) Iseparator     ! ------------------------

      close(unit=2)

      write(*,*) '    ...table is written'

 5001 FORMAT (I10)      ! for Iseparator 
 5002 FORMAT (F20.18)   ! for long positive real 
 5003 FORMAT (D24.17)    ! for coefficients
 5009 FORMAT (A)           ! for zeros (printed as character)

      RETURN
      END

C --------------------------------------------------------------------
