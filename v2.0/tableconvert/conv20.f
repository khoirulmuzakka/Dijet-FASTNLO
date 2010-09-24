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
*
* MW 09/23/2010 add BinWidth array to block A2
* MW 09/17/2010 remove new block A2 and put Imachine at end of block A1
* MW 09/10/2010 Add new "machine readable" block A2 (old block A2->A3)
* MW 08/17/2010 Add user variables in block A1
* KR 09/05/2009 Fixed IPDFdef1 to be 2 resp. 3 for DIS resp. pp/ppbar 
*               instead of 1,2! According to v2.0 table format 1 is
*               reserved for e+e-.
* MW 06/03/2007 first version with correct format (missing reweighting)
* MW 11/22/2006 remove bug in pointer to scale variations
* MW 11/20/2006 first version
*********************************************************************
      implicit none
      CHARACTER*255 FILENAME, FILENAME2
      LOGICAL LHEADONLY
      CHARACTER*4 SHEAD/'head'/
      CHARACTER*255 BUFFER
      

c --- parse command line
      If ( IARGC().LT.1)  FILENAME = 'intable.txt'
      If ( IARGC().LT.2)  FILENAME2 = 'outtable.txt'
      LHEADONLY = .false.
      If ( IARGC().GT.0)  Call GETARG(1,FILENAME)
      If ( IARGC().GT.1)  Call GETARG(2,FILENAME2)
      If ( IARGC().GT.2)  Call GETARG(3,BUFFER)
      
      If(INDEX(BUFFER,SHEAD).NE.0) LHEADONLY = .true.

      Write(*,*) ' '
      Write(*,*) ' '
      Write(*,*) ' '
      Write(*,*) ' **************************************************'
      Write(*,*) ' ******* fastNLO: convert from v1.4 to v2.0 *******'
      Write(*,*) ' **************************************************'
      Write(*,*) ' *  ==> scenario name has to be edited by hand'
      Write(*,*) ' *'
      Write(*,*) ' *  read:', FILENAME
      Write(*,*) ' * Write:', FILENAME2
      Write(*,*) ' *'
      Write(*,*) ' '

c --- read the fastNLO coefficient table from v1.4
      Write(*,*) ' reading v1.4 table ***************'
      Call CV14RD(FILENAME)

c
c - question: Write sum-table - or Write single contributions????
c

      Write(*,*) ' '
      Write(*,*) ' '
      Write(*,*) ' writing v2.0 table(s) ***************'
c --- Write either header, single contributions, or whole table 
      If (LHEADONLY) THEN
         Call CV20WRT(FILENAME2,0) ! header only (one can add hadr-cor or data)
      Else
         Call CV20WRT(FILENAME2,99) ! ------- all: v2.0 sum table
c      Call CV20WRT('contrib1.txt',1) ! LO contrib only
c      Call CV20WRT('contrib2.txt',2) ! NLO contrib only
c      Call CV20WRT('contrib3.txt',3) ! thresh-cor only
      Endif

      Return
      End



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
      integer i,j,k,l,m,n,n1,n2,o,p, nbin, i1,i2   !MW, iscdef
      Include 'conv20.inc'

      Integer nctrmin,nctrmax, mmax
c - assume: max 500 ObsBins & max 50 xNodes
      Double Precision xarray(500,50),rewgt(500,1300), w1,w2

c - new variables that did not exist in v14
      Integer  Ncontrib, Ndim, IDiffBin(2), INormFlag, IDataFlag, IAddMultFlag,
     +     IContrFlag1,IContrFlag2,IContrFlag3, NcontrDescr, NCodeDescr,
     +     IScaleDep, NPDF, NPDFPDG(2), NPDFDim, NFragFunc, NSubpr,
     +     IPDFCoeff, IPDFdef1, IPDFdef2, IPDFdef3,
     +     NScales, NScaleDim, nxall , nscaddr, nscvar
      Double Precision binsize,XNode, hxlim,hx,a
c - variables for diagonalization of Bernstein Scale Interpolation
      Double Precision Bnst(4)
      Double Precision t1,t2
      Double Precision Bnfactor

c - pointer to subprocesses - new ordering
      Integer ipppoint(7),idispoint(3)
      Data ipppoint/1,4,5,6,7,2,3/, idispoint/3,1,2/

c - Ncontrib
      Ncontrib = Nord
c - check old 'Ndimension'
      If (Ndimension.ne.2) then
         Write (*,*) '    Ndimension = ',Ndimension,'   .ne.2'
      endif
      Ndim = Ndimension
c      If  (Nrapidity.eq.1) then
c         Write(*,*) ' only single rapidity range -> set Ndim=1'
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
      If (ireaction.eq.1) then
         NPDF = 1
         NPDFDim = 0
         NSubpr = 3             ! see later how deal with this order-by-order
      else
         NPDF = 2
         NPDFDim = 1
         NSubpr = 7             ! deal with this later order-by-order
      endif
      If (ireaction.eq.2) NPDFPDG(2) = NPDFPDG(1)
      
      NFragFunc = 0
      NScales = 2
      NScaleDim = 1

      Write(*,*) ' '
      Write(*,*) 'CV20WRT:   Write: ', filename
      Write(*,*) ' '

      open(unit=2,file=filename,status='unknown')
c ------------------------------------------------------------------
      Write(2,5001) Iseparator  ! ------------------- block A1 ------
      Write(2,5000) 20000       ! Itabversion
      Write(2,5009) 'FNX9999'   ! ScenName
      Write(2,4999) Ncontrib        ! Ncontrib
      Write(2,4999) 0           ! NMult
      Write(2,4999) 0           ! NData
      Write(2,4999) 0           ! Nuserstrings
c      Write(2,4999) 1           ! Nuserstrings
c      Write(2,5009) 'converted_from_fastNLO_v14_table'   ! info      
      Write(2,4999) 0           ! Nuserint
      Write(2,4999) 0           ! Nuserfloat
      Write(2,4999) 0           ! Imachine
c ------------------------------------------------------------------
      Write(2,5001) Iseparator  ! ------------------- block A2 ------
      Write(2,4999) IXSECTUNITS    ! -> Ipublunits
      Write(2,4999) 3           ! NScDescript 
      Do i=1,3
         Write(2,5009) NAMELABEL(i) ! ScDescript(i)
      Enddo
      Write(2,*) ECMS           ! Ecms
      Write(2,4999) NPOW(1)     ! ILOord
      Write(2,5000) NBINTOT     ! NObsBin
      Write(2,4999) Ndim        ! NDim    (so far =2 - change where needed)
      Do i=Ndim,1,-1
c      Do i=1,Ndim
         Write(2,5009) DIMLABEL(i)
      Enddo
      Do i=Ndim,1,-1
c      Do i=1,Ndim
         Write(2,4999) IDiffBin(i) ! IDiffBin
      Enddo

      i = 0
      Do i1=1,Nrapidity
         Do i2=1,NPT(i1)
            i = i+1             ! --- bin counter - linear loop over all bins
            Do j=Ndim,1,-1
c            Do j=1,Ndim
               If (j.eq.1) then
                  Write(2,*) RAPBIN(i1)
                  Write(2,*) RAPBIN(i1+1)
               elseif (j.eq.2) then
                  Write(2,*) PTBIN(i1,i2)
                  Write(2,*) PTBIN(i1,i2+1)
               endif
            Enddo
         Enddo
      Enddo

      Write(*,*) ' Warning: bin width is computed, assuming incl. jet'
      Write(*,*) '          double differential cross section, including'
      Write(*,*) '          a factor of two (|y| -> y)'
      Do i1=1,Nrapidity
         Do i2=1,NPT(i1)
            binsize = (PTBIN(i1,i2+1)-PTBIN(i1,i2))
     +           * 2d0*(RAPBIN(i1+1)-RAPBIN(i1))
            Write(2,*) binsize
         Enddo
      Enddo
      Write(2,4999) INormFlag   ! INormFlag

c ------------------------------------------------------------------
      Write(2,5001) Iseparator  ! ------------------- block B ------

      If (Iflg.eq.0) then
         goto 999
      Elseif (Iflg.eq.99) then
         nctrmin = 1
         nctrmax = Ncontrib
      Else
         nctrmin = Iflg
         nctrmax = Iflg
         If (Ncontrib.lt. Iflg) then
            Write(*,*) ' contribution ',iflg,' not available'
            stop
         Endif
      Endif
c ------------------------------------------------------------------
      Do n=nctrmin,nctrmax           ! loop over all contributions  B01...
         Write(2,4999) IXSECTUNITS      ! Ixsectunits (for each block)
         Write(2,4999) IDataFlag     ! IDataFlag
         Write(2,4999) IAddMultFlag  ! IAddMultFlag
         If (n.lt.1) then
            Write(*,*) 'n=',n
            Stop
         Elseif (n.eq.1) then
            IContrFlag1 = 1
            IContrFlag2 = 1
            IContrFlag3 = 0
         Elseif (n.eq.2) then
            IContrFlag1 = 1
            IContrFlag2 = 2
            IContrFlag3 = 0
         Elseif (n.eq.3) then
            IContrFlag1 = 2
            IContrFlag2 = 1
            IContrFlag3 = 2
         Endif
         Write(2,4999) IcontrFlag1 ! IcontrFlag1
         Write(2,4999) IcontrFlag2 ! IcontrFlag2
         Write(2,4999) IcontrFlag3 ! IcontrFlag3
         Write(2,4999) NContrDescr ! NcontrDescr
         Do i=1,NContrDescr
            Write(2,5009) POWLABEL(n) ! CtrbDescript(i)
         Enddo
         Write(2,4999) NCodeDescr ! NCodeDescr
         Do i=1,NCodeDescr
            Write(2,5009) CODELABEL(n) ! CodeDescript(i)
         Enddo
c --- (here: neither data - nor multiplicative factors)
         Write(2,4999) IRef          ! IRef
         IScaleDep = 0          ! in v1.4 tables LO always 1st
         If (n.eq.2) IScaleDep = 1   ! in v1.4 tables NLO always 2nd
         If (n.eq.3) IScaleDep = 2   ! in v1.4 tables thresh always 3rd
         Write(2,4999) IScaleDep     ! IScaleDep
         Write(2,4998) Nevt(n)  ! Nevt
         Write(2,4999) Npow(n)  ! Npow
         Write(2,4999) NPDF     ! NPDF
         Do i=1,NPDF
            Write(2,5000) NPDFPDG(i)  ! NPDFPDG
         Enddo
         Write(2,4999) NPDFDim     ! NPDFdim
         Write(2,4999) NFragFunc   ! NFragFunc
         Write(2,4999) 0           ! NFragDim

c -   NSubProc
c         Write(2,*) NSubpr      ! Nsubproc   <<< modify for DIS & pp
         If (NPDFDim.eq.0) then ! ---- DIS
            If (n.eq.1) then ! --- LO
               Write(2,4999) 2
            Elseif (n.eq.2) Then
               Write(2,4999) 3
            Endif
         ElseIf(NPDFDim.eq.1) then ! --------- pp
            If (n.eq.1 .or. n.eq.3) then    ! ------ LO or threshcor
               Write(2,4999) 6
            Elseif (n.eq.2) then ! ---- NLO
               Write(2,4999) 7
            Endif
         Endif

         If (ireaction.eq.1) Then ! ep
            IPDFdef1 = 2        ! DIS
            IPDFdef2 = 1        ! NC DIS
            If (n.eq.1) Then
               IPDFdef3 = 2     ! --- only Delta, G
            Elseif (n.eq.2) Then
               IPDFdef3 = 3     ! --- Delta, G, Sigma
            Else
               Write(*,*) 'DIS wrong n'
               Stop
            Endif
         Elseif (ireaction.ge.2) Then ! pp:2 - ppbar:3
            If (n.eq.1 .or. n.eq.3) Then    ! ------ LO or threshcor
               IPDFdef1 = 3
               IPDFdef2 = 1
               IPDFdef3 = 1
c               IPDFCoeff = 2000101
            Elseif (n.eq.2) Then ! ---- NLO
               IPDFdef1 = 3
               IPDFdef2 = 1
               IPDFdef3 = 2
            Else
               Write(*,*) " strange: n outside 1,2,3 in pp table"
            Endif
         Endif
         Write(2,5000) IPDFdef1   ! IPDFdef1
         Write(2,5000) IPDFdef2   ! IPDFdef2
         Write(2,5000) IPDFdef3   ! IPDFdef3

         i = 0
         Do i1=1,Nrapidity
            Do i2=1,NPT(i1)
               i = i+1          ! --- bin counter -> linear bin numbering
               Write(2,5000) Nxtot ! Nxtot(1) - new for each ObsBin
               Do j=1,Nxtot
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
                  xarray(i,j) = XNode ! --- use later for PDF unweighting
                  Write(2,5002) XNode
               Enddo
            Enddo
         Enddo

         Write(2,*) NScales     ! NScales
         Write(2,*) NScaleDim   ! NScaleDim
         Do i=1,Nscales
c            Write(2,*) 1        ! Iscale(i) <- all scales = 1st Dimension
            Write(2,*) 0        ! Iscale(i) <- all scales = 1st Dimension
         Enddo
         Do i=1,Nscaledim
            Write(2,*) 1        ! NScaleDescript(i) <- one descriptive string
            Write(2,5009) SCALELABEL           
         Enddo
         
         If (IScaleDep.eq.0) then
            nscvar = 1
cMW - iscdef not needed - just set MURSCALE = 1 for the LO contribution
c            If (ireaction.eq.1) Then
c               iscdef = 2       ! in DIS v1.4 tables default scale = 2
c            Elseif (ireaction.ge.2) Then
c               iscdef = 3       ! in pp v1.4 tables default scale = 3
c            Else
c               Write(*,*) 'unknown ireaction'
c               Stop
c            Endif
         Else
            nscvar = NscaleVar
         endif
         
         Do i=1,NscaleDim
            Write(2,*) nscvar   ! NscaleVar
            Write(2,*) NscaleBin ! NscaleBin
         Enddo
         Do i=1,Nscaledim
            if (nscvar.eq.1) then
cMW               Write(2,*) MURSCALE(iscdef) ! Scalefac(i)(j)
               Write(2,*) 1d0 ! Scalefac(i)(j)
cdebug
cdebug               Write(*,*)"CONV20: IScaleDim,IScaleVar,ScaleFac: ",
cdebug     >              i,iscdef,MURSCALE(iscdef)
cdebug
            else
               Do j=1,nscvar
                  Write(2,*) MURSCALE(j) ! Scalefac(i)(j)               
cdebug
cdebug                  Write(*,*)"CONV20: IScaleDim,IScaleVar,ScaleFac: ",
cdebug     >                 i,j,MURSCALE(j)
cdebug
               Enddo
            endif
         Enddo

         i = 0
         Do i1=1,Nrapidity
            Do i2=1,NPT(i1)
               i = i+1          ! --- bin counter
               Do j=1,Nscaledim
                  Do k=1,nscvar
                     Do l=1,NscaleBin
                        If (n.eq.1) Then
                           Write(2,*)MURVAL(i1,i2,l) !ScaleNode(ijkl)
                        Else
                           Write(2,*)MURVAL(i1,i2,l)*MURSCALE(k)!ScaleNode(ijkl)
                        Endif
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo

c - undo PDF reweighting -> compute factors 
         i = 0
         Do i1=1,Nrapidity
            Do i2=1,NPT(i1)
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
                     Write(*,*) " can't deal with NPDFDim=2 here "
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


c --- diagonalize Bernstein scale interpolation in v1.4 (not in v1.6)
         If (ITabversion.eq.14000 .and. NScaleBin.gt.2) Then
            i = 0
            Do i1=1,Nrapidity
               Do i2=1,NPT(i1)
                  i = i+1       ! --- bin counter ObsBin
                  Do j=1,Nscaledim !                 scaledim
                     Do k=1,nscvar !                 scalevar
c                        Do l=1,NscaleBin !           scalebin
                        if (NPDFDim.eq.0) nxall=Nxtot
                        if (NPDFDim.eq.1) nxall=(Nxtot*Nxtot+Nxtot)/2
                        Do m=1,nxall !            xbins
                           Do n1=1,Nsubproc !     subprocesses
                              If (n.eq.1) then
                                 nscaddr = 1
                              Else
                                 nscaddr = 1+k+(n-2)*nscalevar
                              Endif

                              If (NScaleBin.eq.3) Then ! 3 Bernstein Scalebins
                                 t1 = 1./2.
                                 Bnfactor =  1./(2.*(1.-t1)*t1)
                                 Bnst(1) = array(i,m,n1,nscaddr,1)
                                 Bnst(2) = array(i,m,n1,nscaddr,2)
                                 Bnst(3) = array(i,m,n1,nscaddr,3)
                                 array(i,m,n1,nscaddr,1) = 
     +                                1d0*Bnst(1) - Bnfactor*(1.-t1)**2 *Bnst(2)             
                                 array(i,m,n1,nscaddr,2) = 
     +                                            + Bnfactor            *Bnst(2)             
                                 array(i,m,n1,nscaddr,3) = 
     +                                            - Bnfactor*(t1)**2    *Bnst(2)+ 1d0*Bnst(3)

                              Elseif (NScaleBin.eq.4) Then ! 4 Bernstein Scalebins
                                 t1 = 1./3.
                                 t2 = 2./3.
                                 Bnfactor =  1./(3.*(1.-t1)**2*t1 * 3.*(1.-t2)*(t2)**2 - 3.*(1.-t2)**2*t2 * 3.*(1.-t1)*(t1)**2)
                                 Bnst(1) = array(i,m,n1,nscaddr,1)
                                 Bnst(2) = array(i,m,n1,nscaddr,2)
                                 Bnst(3) = array(i,m,n1,nscaddr,3)
                                 Bnst(4) = array(i,m,n1,nscaddr,4)
                                        
                                 array(i,m,n1,nscaddr,1) = 
     &                                1d0*Bnst(1)+
     & Bnfactor* (-3.*(1.-t2)*(t2)**2*(1.-t1)**3+ 3.*(1.-t1)*(t1)**2*(1.-t2)**3) *Bnst(2)+ 
     & Bnfactor* (-3.*(1.-t1)**2*(t1)*(1.-t2)**3+ 3.*(1.-t2)**2*(t2)* *(1.-t1)**3) *Bnst(3)

                                 array(i,m,n1,nscaddr,2) = 
     & Bnfactor*(3.*(1.-t2)*(t2)**2) *Bnst(2)+
     & Bnfactor*(- 3.*(1.-t2)**2*(t2)) *Bnst(3)

                                 array(i,m,n1,nscaddr,3) = 
     & Bnfactor*(-3.*(1.-t1)*(t1)**2) *Bnst(2) +
     & Bnfactor*(3.*(1.-t1)**2*(t1)) *Bnst(3)

                                 array(i,m,n1,nscaddr,4) = 
     & Bnfactor* (-3.*(1.-t2)*(t2)**2*(t1)**3+ 3.*(1.-t1)*(t1)**2*(t2)**3)  *Bnst(2) +
     & Bnfactor* (-3.*(1.-t1)**2*(t1)*(t2)**3+ 3.*(1.-t2)**2*(t2)*(t1)**3)  *Bnst(3) +
     & 1d0 *Bnst(4)
 
                              Else
                                 Write(*,*)
     >                                'not implemented: NScaleBin='
     >                                ,NScaleBin
                                 Stop
                              Endif
                           Enddo
                        Enddo
                     Enddo                     
                  Enddo
               Enddo
            Enddo
         Endif                  ! Bernstein diagonalization

c --- write SigmaTilde array
         i = 0
         Do i1=1,Nrapidity
            Do i2=1,NPT(i1)
               i = i+1          ! --- bin counter ObsBin
               Do j=1,Nscaledim !                 scaledim
                  Do k=1,nscvar !                 scalevar
                     Do l=1,NscaleBin !           scalebin
                        if (NPDFDim.eq.0) nxall=Nxtot
                        if (NPDFDim.eq.1) nxall=(Nxtot*Nxtot+Nxtot)/2
                        Do m=1,nxall !            xbins
                           Do n2=1,Nsubproc !     subprocesses
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
                                       Write(2,5009) '0'
                                    else
                                       Write(2,5003) a
                                    endif
                                 endif

                              else ! ------- pp/ppbar
                                 if ((n.eq.1 .or. n.eq.3).and. n1.eq.3) a=
     +                                (array(i,m,n1-1,nscaddr,l)
     +                                +array(i,m,n1,  nscaddr,l))/2d0
     +                                *rewgt(i,m)
                                 if ((n.eq.1 .or. n.eq.3) .and. n1.eq.2) then
c                                    if (a.eq.0d0) then  ! --- remove later
c                                       Write(2,5009) '0'
c                                    else
c                                       Write(2,5003) a
c                                    endif
                                    continue
                                 else
                                    if (a.eq.0d0) then
                                       Write(2,5009) '0'
                                    else
                                       Write(2,5003) a
                                    endif
                                 endif

c ------ old test for problem that has now diappeared ----------------
c - Problem: for scalevar 1&2 threshcor subproc 2,3 not identical
                                 if (n.eq.3 .and. n1.eq.3 .and. a.ne.0d0 .and.
     +                                (array(i,m,n1-1,nscaddr,l)
     +                                /array(i,m,n1,  nscaddr,l)).ne.1d0)
     +                                Write(*,*)i,array(i,m,n1-1,nscaddr,l)
     +                                /array(i,m,n1,  nscaddr,l),k,l,m,
     +                                nscaddr,n

                              endif ! ---- DIS or pp/ppbar?
c ---------------------------------------------------------
                           Enddo
                        Enddo
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo
         Write(2,5001) Iseparator ! ------------------- block B ------
      Enddo
      Write(2,5001) Iseparator  ! --------- add. check for EOF --

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
         Write(*,*) '          fastNLO:  table file not found ',
     +        '  -  IOSTAT = ',IFILE
         STOP
      ENDIF

      Read(2,*) i
      if (i.ne.iseparator) goto 999
      Read(2,*) ITABVERSION
      Write(*,*) "#       tableformat is version",real(itabversion)/10000d0
      if (ITABVERSION.ne.14000 .and. ITABVERSION.ne.16000) then
         Write(*,*) '#     => this conversion works only for v1.4 and v1.6'
         stop
      endif
      Read(2,*) i
      if (i.ne.iseparator) goto 999
c  -----------------------------------
      Read(2,*) IREACTION
      Read(2,*) ECMS
      Read(2,*) IXSECTUNITS
      Do i=1,5
         Read(2,*) NAMELABEL(i)
      Enddo
      Read(2,*) IPROC
      Read(2,*) IALGO
      Read(2,*) JETRES1
      Read(2,*) JETRES2
      Read(2,*) NORD
      Do i=1,nord
         Read(2,*) NPOW(i)
      Enddo
      Do i=1,nord
         Read(2,*) POWLABEL(i)
         Read(2,*) CODELABEL(i)
      Enddo
      Read(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      Do i=1,nord
         Read(2,*) NEVT(i)
      Enddo
      Do i=1,nord
         Write(*,5000) ' #      ',NEVT(i),
     +        ' events in ',POWLABEL(i)
      Enddo
      Read(2,*) NXTOT
      Write(*,*) "#           No. of x bins: ",NXTOT
      Read(2,*) IXSCHEME
      Read(2,*) IPDFWGT
      Read(2,*) IREF
      Read(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      Read(2,*) NBINTOT
      Read(2,*) NDIMENSION
      Do i=1,ndimension
         Read(2,*) DIMLABEL(i)
      Enddo
      
      Read(2,*) NRAPIDITY
      Do i=1,nrapidity+1
         Read(2,*) RAPBIN(i)
      Enddo
      Do i=1,nrapidity
         Read(2,*) NPT(i)
      Enddo
      Do i=1,nrapidity
         do j=1,NPT(i)+1
            Read(2,*) PTBIN(i,j)
         Enddo
      Enddo
      Read(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      Do i=1,nrapidity
         Do j=1,NPT(i)
            Read(2,*) XLIMIT(i,j)
         Enddo
      Enddo
      Read(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      Read(2,*) SCALELABEL
      Read(2,*) NSCALEBIN
      Do i=1,nrapidity
         Do j=1,NPT(i)
            Do k=1,NSCALEBIN
               Read(2,*) MURVAL(i,j,k)
            Enddo
         Enddo
      Enddo
      Read(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      Do i=1,nrapidity
         Do j=1,NPT(i)
            Do k=1,NSCALEBIN
               Read(2,*) MUFVAL(i,j,k)
            Enddo
         Enddo
      Enddo
      Read(2,*) i
      if (i.ne.iseparator) goto 999
c   -----------------------------------
      Read(2,*) NSCALEVAR
      Do i=1,NSCALEVAR
         Read(2,*) MURSCALE(i)
      Enddo
      Do i=1,NSCALEVAR
         Read(2,*) MUFSCALE(i)
      Enddo
      
      Read(2,*) i
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
      Do i=1,nrapidity
         Do j=1,NPT(i)
            nbin=nbin+1         ! linear numbering of (rap,pT) bins
            Do k=1,NXSUM        ! tot. No. x-Bins
               Do m=1,Nsubproc  ! No. of Subproc
                  Do n=1,1+NSCALEVAR*(NORD-1) ! LO & NLO & w/ scale var
                     Do l=1,NSCALEBIN ! No. of Bins in Scale
                        Read(2,*) array(nbin,k,m,n,l)
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo
      Enddo
      Read(2,*) i
      If (i.ne.iseparator) Goto 999
      Close(2)

      Return
 999  Continue

      Close (2) 
      Write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      Write(*,*) " >>>>>   fastNLO error in table format "
      Write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      Stop
      Return

 5000 Format (A,I12,A,A64) 
      End

*******************************************************************
