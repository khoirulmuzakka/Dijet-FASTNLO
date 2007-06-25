      PROGRAM CONVERT
*********************************************************************
*   M. Wobisch      Nov 20, 2006
*   
*   convert fastNLO tables from v1.4 to v2.0
*  
*  
* contains:
*   CV20WRT  writes v2.0 tables
*   CV14RD   read v1.4 tables
*
* status:
*   block A correct  (fixed Ndim meaningful for 1d-scenarios?)
*   block B correct  (see below)
*
* check:
*   - check 1/(PDF reweighting) factors
*   - check pp reduction of subprocesses (2,3 for LO-type)
*           --> much too small table (length o.k. but too many zeros)
* to do:
*   - remove DIS subprocesses which are empty
*   - set codes for PDF Linear combinations correctly
*   - reorder pp subprocesses (move 2,3 to end)
*
* MW 06/03/2007 first version with correct format (missing reweighting)
* MW 11/22/2006 remove bug in pointer to scale variations
* MW 11/20/2006 first version
*********************************************************************
      implicit none
      CHARACTER*255 FILENAME, FILENAME2

c --- parse command line
      IF ( IARGC().LT.1)  FILENAME = 'intable.txt'
      IF ( IARGC().LT.2)  FILENAME2 = 'outtable.txt'
      IF ( IARGC().GT.0)  CALL GETARG(1,FILENAME)
      IF ( IARGC().GT.1)  CALL GETARG(2,FILENAME2)

      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) ' **************************************************'
      WRITE(*,*) ' ******* fastNLO: convert from v1.4 to v1.5 *******'
      WRITE(*,*) ' **************************************************'
      WRITE(*,*) ' *  ==> scenario name has to be edited by hand'
      WRITE(*,*) ' *'
      WRITE(*,*) ' *  read:', FILENAME
      WRITE(*,*) ' * write:', FILENAME2
      WRITE(*,*) ' *'
      WRITE(*,*) ' '

c --- read the fastNLO coefficient table from v1.4
      WRITE(*,*) ' reading v1.4 table ***************'
      call CV14RD(FILENAME)

c
c - question: write sum-table - or write single contributions????
c

      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) ' writing v2.0 table(s) ***************'
c --- write either header, single contributions, or whole table 
c      call CV20WRT('head.txt',0) ! header only (one can add hadr-cor or data)
c      call CV20WRT('contrib1.txt',1) ! LO contrib only
c      call CV20WRT('contrib2.txt',2) ! NLO contrib only
c      call CV20WRT('contrib3.txt',3) ! thresh-cor only
      call CV20WRT(FILENAME2,99) ! ------- all: v2.0 sum table


      RETURN
      END



*******************************************************************
      SUBROUTINE CV20WRT(FILENAME, IFLG)
*-----------------------------------------------------------------
* M. Wobisch   write v2.0 ASCII table
*
* input: FILENAME  name of table
*        IFLG       0: only header is written (can be used as template
*                                        for adding data or hadr.cor.) 
*                  99: whole table is written -> "v2.0 sum table"
*                else: No. of single contribution to be written (can be
*                      used in "joiner" to create sum table, together with
*                      other single contributions such as New Physics, e/w.
*
* MW 11/20/2006 
*-----------------------------------------------------------------
      implicit none
      integer Iflg
      CHARACTER*(*) filename
      CHARACTER*(33) label(0:2)
      integer i,j,k,l,m,n,n1,n2,o,p, nbin, i1,i2
      INCLUDE 'conv20.inc'

      Integer nctrmin,nctrmax, mmax
c - assume: max 500 ObsBins & max 50 xNodes
      Double Precision xarray(500,50),rewgt(500,1300), w1,w2

c - new variables that did not exist in v14
      Integer  Ncontrib, Ndim, IDiffBin(2), INormFlag, IDataFlag, IAddMultFlag,
     +     IContrFlag1,IContrFlag2,IContrFlag3, NcontrDescr, NCodeDescr,
     +     IScaleDep, NPDF, NPDFPDG(2), NPDFDim, NFragFunc, NSubpr,
     +     IPDFCoeff, IPDFdef1, IPDFdef2, IPDFdef3,
     +     NScales, NScaleDim, nxall , nscaddr, nscvar
      Double Precision XNode, hxlim,hx,a

c - pointer to subprocesses - new ordering
      Integer ipppoint(7),idispoint(3)
      Data ipppoint/1,4,5,6,7,2,3/, idispoint/3,1,2/

c - Ncontrib
      Ncontrib = Nord
c - check old 'Ndimension'
      If (Ndimension.ne.2) then
         write (*,*) '    Ndimension = ',Ndimension,'   .ne.2'
      endif
      Ndim = Ndimension
c      If  (Nrapidity.eq.1) then
c         write(*,*) ' only single rapidity range -> set Ndim=1'
c         Ndim = 1
c         DIMLABEL(1) = DIMLABEL(2)
c      endif

c - differential or binned?  -> so far always binned
      IDiffBin(1) = 2
      IDiffBin(2) = 2
c - normalization in last dimension?
      INormFlag = 0
      IDataFlag = 0
      IAddMultFlag = 0
      IContrFlag1 = 1
      IContrFlag2 = 0
      IContrFlag3 = 0
      NContrDescr = 1
      NCodeDescr = 1
      
      NPDFPDG(1) = 2212         ! PDG code proton
      NPDFPDG(2) = -2212        !          pbar
      if (ireaction.eq.1) then
         NPDF = 1
         NPDFDim = 0
         NSubpr = 3             ! see later how deal with this order-by-order
      else
         NPDF = 2
         NPDFDim = 1
         NSubpr = 7             ! deal with this later order-by-order
      endif
      if (ireaction.eq.2) NPDFPDG(2) = NPDFPDG(1)
      
      NFragFunc = 0
      NScales = 2
      NScaleDim = 1

      WRITE(*,*) ' '
      WRITE(*,*) 'CV20WRT:   write: ', filename
      WRITE(*,*) ' '

      open(unit=2,file=filename,status='unknown')
c ------------------------------------------------------------------
      WRITE(2,5001) Iseparator  ! ------------------- block A1 ------
      WRITE(2,5000) 20000       ! Itabversion
      WRITE(2,5009) 'FNX9999'   ! ScenName
      WRITE(2,4999) Ncontrib        ! Ncontrib
      WRITE(2,4999) 0        ! NMult
      WRITE(2,4999) 0        ! NData
c ------------------------------------------------------------------
      WRITE(2,5001) Iseparator  ! ------------------- block A2 ------
      WRITE(2,4999) IXSECTUNITS    ! -> Ipublunits
      WRITE(2,4999) 3           ! NScDescript 
      do i=1,3
         WRITE(2,5009) NAMELABEL(i) ! ScDescript(i)
      enddo
      WRITE(2,*) ECMS           ! Ecms
      WRITE(2,4999) NPOW(1)     ! ILOord
      WRITE(2,5000) NBINTOT     ! NObsBin
      WRITE(2,4999) Ndim        ! NDim    (so far =2 - change where needed)
      Do i=1,Ndim
         WRITE(2,5009) DIMLABEL(i)
      Enddo
      Do i=1,Ndim
         WRITE(2,4999) IDiffBin(i) ! IDiffBin
      Enddo

      i = 0
      do i1=1,Nrapidity
         do i2=1,NPT(i1)
            i = i+1             ! --- bin counter - linear loop over all bins
            do j=1,Ndim
               if (j.eq.1) then
                  WRITE(2,*) RAPBIN(i1)
                  WRITE(2,*) RAPBIN(i1+1)
               elseif (j.eq.2) then
                  WRITE(2,*) PTBIN(i1,i2)
                  WRITE(2,*) PTBIN(i1,i2+1)
               endif
            enddo
         enddo
      enddo
      WRITE(2,4999) INormFlag   ! INormFlag

c ------------------------------------------------------------------
      WRITE(2,5001) Iseparator  ! ------------------- block B ------

      If (Iflg.eq.0) then
         goto 999
      Elseif (Iflg.eq.99) then
         nctrmin = 1
         nctrmax = Ncontrib
      Else
         nctrmin = Iflg
         nctrmax = Iflg
         If (Ncontrib.lt. Iflg) then
            write(*,*) ' contribution ',iflg,' not available'
            stop
         Endif
      Endif
c ------------------------------------------------------------------
      do n=nctrmin,nctrmax           ! loop over all contributions  B01...
         WRITE(2,4999) IXSECTUNITS      ! Ixsectunits (for each block)
         WRITE(2,4999) IDataFlag     ! IDataFlag
         WRITE(2,4999) IAddMultFlag  ! IAddMultFlag
         If (n.le.2) then
            IContrFlag1 = 1
            IContrFlag2 = 0
            IContrFlag3 = 0
         elseif (n.eq.3) then
            IContrFlag1 = 2
            IContrFlag2 = 1
            IContrFlag3 = 0
         endif
         WRITE(2,4999) IcontrFlag1 ! IcontrFlag1
         WRITE(2,4999) IcontrFlag2 ! IcontrFlag2
         WRITE(2,4999) IcontrFlag3 ! IcontrFlag3
         WRITE(2,4999) NContrDescr ! NcontrDescr
         do i=1,NContrDescr
            WRITE(2,5009) POWLABEL(n) ! CtrbDescript(i)
         enddo
         WRITE(2,4999) NCodeDescr ! NCodeDescr
         do i=1,NCodeDescr
            WRITE(2,5009) CODELABEL(n) ! CodeDescript(i)
         enddo
c --- (here: neither data - nor multiplicative factors)
         WRITE(2,4999) IRef          ! IRef
         IScaleDep = 0          ! in v1.4 tables LO always 1st
         If (n.eq.2) IScaleDep = 1   ! in v1.4 tables NLO always 2nd
         If (n.eq.3) IScaleDep = 2   ! in v1.4 tables thresh always 3rd
         WRITE(2,4999) IScaleDep     ! IScaleDep
         WRITE(2,4998) Nevt(n)  ! Nevt
         WRITE(2,4999) Npow(n)  ! Npow
         WRITE(2,4999) NPDF     ! NPDF
         do i=1,NPDF
            WRITE(2,5000) NPDFPDG(i)  ! NPDFPDG
         enddo
         WRITE(2,4999) NPDFDim     ! NPDFdim
         WRITE(2,4999) NFragFunc   ! NFragFunc
         WRITE(2,4999) 0           ! NFragDim

c -   NSubProc
c         WRITE(2,*) NSubpr      ! Nsubproc   <<< modify for DIS & pp
         If (NPDFDim.eq.0) then ! ---- DIS
            If (n.eq.1) then ! --- LO
               WRITE(2,4999) 2
            Elseif (n.eq.2) Then
               WRITE(2,4999) 3
            Endif
         ElseIf(NPDFDIm.eq.1) then ! --------- pp
            If (n.eq.1 .or. n.eq.3) then    ! ------ LO or threshcor
               WRITE(2,4999) 6
            Elseif (n.eq.2) then ! ---- NLO
               Write(2,4999) 7
            Endif
         Endif

         If (ireaction.eq.1) then
            IPDFdef1 = 1
            IPDFdef2 = 1
            If (n.eq.1) then
               IPDFdef3 = 2     ! --- only Delta, G
            Elseif (n.eq.2) then
               IPDFdef3 = 3     ! --- Delta, G, Sigma
            Else
               write(*,*) 'DIS wrong n'
               STOP
            Endif
         elseif (ireaction.ge.2) then
            If (n.eq.1 .or. n.eq.3) then    ! ------ LO or threshcor
               IPDFdef1 = 2
               IPDFdef2 = 1
               IPDFdef3 = 1
               IPDFCoeff = 2000101
            Elseif (n.eq.2) then ! ---- NLO
               IPDFdef1 = 2
               IPDFdef2 = 1
               IPDFdef3 = 2
            Else
               write(*,*) " strange: n outside 1,2,3 in pp table"
            Endif
         endif
         WRITE(2,5000) IPDFdef1   ! IPDFdef1
         WRITE(2,5000) IPDFdef2   ! IPDFdef2
         WRITE(2,5000) IPDFdef3   ! IPDFdef3

c         WRITE(2,5000) Nxtot       ! Nxtot(1) - obsolete see below
         i = 0
         do i1=1,Nrapidity
            do i2=1,NPT(i1)
               i = i+1          ! --- bin counter -> linear bin numbering
               WRITE(2,5000) Nxtot ! Nxtot(1) - new for each ObsBin
               do j=1,Nxtot
                  XNode = xlimit(i1,i2)
                  if (ixscheme.eq.2) then
                     hxlim = - sqrt(-log10(xnode))
                     hx = hxlim *(1d0 - dble(j-1)/dble(Nxtot)) 
                     XNode = 10**-(hx*hx)
                  elseif (ixscheme.eq.1) then
                     hxlim = log10(xnode)
                     hx = hxlim *(1d0 - dble(j-1)/dble(Nxtot))
                     XNode = 10**(hx) 
                  endif
                  xarray(i,j) = XNode ! --- use later for unweighting
                  WRITE(2,5002) XNode
               enddo
            enddo
         enddo

         WRITE(2,*) NScales     ! NScales
         WRITE(2,*) NScaleDim   ! NScaleDim
         do i=1,Nscales
            WRITE(2,*) 1        ! Iscale(i) <- all scales = 1st Dimension
         enddo
         do i=1,Nscaledim
            WRITE(2,*) 1        ! NScaleDescript(i) <- one descriptive string
            WRITE(2,5009) SCALELABEL           
         enddo

         If (IScaleDep.eq.0) then
            nscvar = 1
         else
            nscvar = NscaleVar
         endif

         do i=1,NscaleDim
            WRITE(2,*) nscvar ! NscaleVar
            WRITE(2,*) NscaleBin ! NscaleBin
         enddo
         do i=1,Nscaledim
            do j=1,nscvar
               WRITE(2,*) MURSCALE(j) ! Scalefac(i)(j)               
            enddo
         enddo

         i = 0
         do i1=1,Nrapidity
            do i2=1,NPT(i1)
               i = i+1          ! --- bin counter
               do j=1,Nscaledim
                  do k=1,nscvar
                     do l=1,NscaleBin
                        If (n.eq.1) Then
                           WRITE(2,*)MURVAL(i1,i2,l) !ScaleNode(ijkl)
                        Else
                           WRITE(2,*)MURVAL(i1,i2,l)*MURSCALE(k)!ScaleNode(ijkl)
                        Endif
                     enddo
                  enddo
               enddo
            enddo
         enddo

c - undo PDF reweighting -> compute factors 
         i = 0
         do i1=1,Nrapidity
            do i2=1,NPT(i1)
               i = i+1          ! --- bin counter ObsBin
               if (NPDFDim.eq.0) nxall=Nxtot
               if (NPDFDim.eq.1) nxall=(Nxtot*Nxtot+Nxtot)/2
               m=0
               Do j=1,Nxtot
                  If (NPDFDim.eq.0) Then ! ---- DIS
                     mmax=1
                  Elseif (NPDFDim.eq.1) Then ! ---- pp / pp-bar
                     mmax =j
                  Else
                     WRITE(*,*) " can't deal with NPDFDim=2 here "
                     STOP
                  Endif
                  Do k=1,mmax
                     m = m+1 ! Bin No. in linear x-array
c ---------- works only for standard fastNLO PDF weighting formula
                     w1= sqrt(xarray(i,j))/(1d0-0.99d0*xarray(i,j))**3
                     w2= sqrt(xarray(i,k))/(1d0-0.99d0*xarray(i,k))**3
                     If (NPDFdim.eq.0) Then ! --- DIS
                        rewgt(i,m) = w1
                     Elseif (NPDFdim.eq.1) Then ! --- pp / pp-bar
                        rewgt(i,m) = w1*w2 ! --- undo weighting
                     Endif
                  Enddo
               Enddo
            Enddo
         Enddo

c --- write SigmaTilde array
         i = 0
         do i1=1,Nrapidity
            do i2=1,NPT(i1)
               i = i+1          ! --- bin counter ObsBin
               do j=1,Nscaledim !                 scaledim
                  do k=1,nscvar !                 scalevar
                     do l=1,NscaleBin !           scalebin
                        if (NPDFDim.eq.0) nxall=Nxtot
                        if (NPDFDim.eq.1) nxall=(Nxtot*Nxtot+Nxtot)/2
                        do m=1,nxall !            xbins
                           do n2=1,Nsubproc !     subprocesses
c - reorder subprocesses - use pointer
                              if (ireaction.eq.1) then ! DIS
                                 n1 = idispoint(n2)
                              else
                                 n1 = ipppoint(n2)
c                                 n1=n2
                              endif
                              if (n.eq.1) then
                                 nscaddr = 1
                              else
                                 nscaddr = 1+k+(n-2)*nscalevar
                              endif
                              a = array(i,m,n1,nscaddr,l)
                              a = a*rewgt(i,m) ! --- PDF unweighting

c --- special cases: in pp in LO merge 2nd & 3rd subproc for pp/ppbar ---
c ---                in DIS remove Sigma at LO
c --- to do: what is variable "n"? -> outermost loop: contributions

                              if (ireaction.eq.1) then ! --- DIS
                                 if (n.eq.1 .and. n1.eq.2) then
                                    continue
                                 else
                                    if (a.eq.0d0) then
                                       WRITE(2,5009) '0'
                                    else
                                       WRITE(2,5003) a
                                    endif
                                 endif

                              else ! ------- pp/ppbar
                                 if ((n.eq.1 .or. n.eq.3).and. n1.eq.3) a=
     +                                (array(i,m,n1-1,nscaddr,l)
     +                                +array(i,m,n1,  nscaddr,l))/2d0
                                 if ((n.eq.1 .or. n.eq.3) .and. n1.eq.2) then
c                                    if (a.eq.0d0) then  ! --- remove later
c                                       WRITE(2,5009) '0'
c                                    else
c                                       WRITE(2,5003) a
c                                    endif
                                    continue
                                 else
                                    if (a.eq.0d0) then
                                       WRITE(2,5009) '0'
                                    else
                                       WRITE(2,5003) a
                                    endif
                                 endif

c ------ old test for problem that has now diappeared ----------------
c - Problem: for scalevar 1&2 threshcor subproc 2,3 not identical
                                 if (n.eq.3 .and. n1.eq.3 .and. a.ne.0d0 .and.
     +                                (array(i,m,n1-1,nscaddr,l)
     +                                /array(i,m,n1,  nscaddr,l)).ne.1d0)
     +                                write(*,*)i,array(i,m,n1-1,nscaddr,l)
     +                                /array(i,m,n1,  nscaddr,l),k,l,m,
     +                                nscaddr,n

                              endif ! ---- DIS or pp/ppbar?
c ---------------------------------------------------------
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         WRITE(2,5001) Iseparator ! ------------------- block B ------
      enddo
      WRITE(2,5001) Iseparator  ! --------- add. check for EOF --

 999  continue
      close(unit=2)

 4998 FORMAT (I14)              ! for integers like Nevt
 4999 FORMAT (I2)               ! for integers < 100
 5000 FORMAT (I5)               ! for integers < 100k
 5001 FORMAT (I10)              ! for Iseparator 
 5002 FORMAT (F20.18)           ! for long positive real 
 5003 FORMAT (D24.17)           ! for coefficients
 5009 FORMAT (A)                ! for charactrs & zeros (printed as character)

      RETURN
      END

*******************************************************************

      SUBROUTINE CV14RD(FILENAME)
*-----------------------------------------------------------------
* M. Wobisch   read ASCII table of perturbative coefficients from v1.4
*
* input: FILENAME  name of table
*
* MW 04/15/2005
* MW 01/30/2006 -> now reading tableformat v1.4
* MW 06/02/2006 -> check table format
*-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) FILENAME
      CHARACTER*255 BUFFER
      INTEGER IFIRST, IFILE, I,J,K,L,M,N,   NBIN,NX
      INCLUDE 'conv20.inc'

      DATA IFIRST/0/
      SAVE IFIRST

      OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=IFILE)
      IF (IFILE .ne. 0) THEN
         WRITE(*,*) '          fastNLO:  table file not found ',
     +        '  -  IOSTAT = ',IFILE
         STOP
      ENDIF

      READ(2,*) i
      if (i.ne.iseparator) goto 999
      READ(2,*) ITABVERSION
      WRITE(*,*) "#       tableformat is version",real(itabversion)/10000d0
      if (ITABVERSION.ne.14000) then
         write(*,*) '#     => this conversion works only for version 1.4'
         stop
      endif
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
            nbin=nbin+1         ! linear numbering of (rap,pT) bins
            do k=1,NXSUM        ! tot. No. x-Bins
               do m=1,Nsubproc  ! No. of Subproc
                  do n=1,1+NSCALEVAR*(NORD-1) ! LO & NLO & w/ scale var
                     do l=1,NSCALEBIN ! No. of Bins in Scale
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

      RETURN
 999  continue

      close (2) 
      write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      write(*,*) " >>>>>   fastNLO error in table format "
      write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      stop
      RETURN

 5000 FORMAT (A,I12,A,A64) 
      END

*******************************************************************
