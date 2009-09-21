*******************************************************************
*******************************************************************
* fastNLO user code            T. Kluge, M. Wobisch v1.4 02/01/2006      
*                                                   v1.5 31/05/2007  
*                                                   v2.0 16/07/2007 
*   >> this code should not be edited
*   >> if you have questions, please contact 
*      the authors at:  fastnlo@cedar.ac.uk
*
* fnxNNNN.f     
*     -> where "Xnnnn" is the respective fastNLO scenario name
*            "X" is an abbreviation for the collider
*                   "L" LHC
*                   "T" Tevatron
*                   "H" HERA 
*                   "R" RHIC
*            "nnnn" is a scenario number
*
* ===================== ideas ============================
* - propose to call main routine fnx9999 (just like scenario)
*           -> any reasons not to?
* - include table-checks into read/write routine - or separate check routine?
*
*
*--------- v2.0 routines - status --------------------------------
*  fx9999cc     main routine - calls other code - need cosmetics        (works)
*                                                 -> need to linearize k1,k2 
*  fx9999in     initialization (to do: check consistency)               (works)
*  fx9999pt     find pointers to contributions/scales                   (works)
*                 - interpret contribution selection/consistency (advanced)
*                 - maybe reorder contributions/optimize PDF access (not yet)
*                 - check availability of the contributions (not yet)
*  fx9999gp     get PDFs      - missing details                         (works)
*  fx9999pl     compute PDF linear combinations (done - cosmetics)      (works)
*  fx9999mt     multiply coefficients and PDFs                       (advanced)
*  fx9999pr     print results                                           (works)
*  fx9999nf     print scenario information (physics & technical)        (works)
*  fx9999nm     normalize distribution by its own integral              (to do)
*
* external:
*  fx9999rw     read (or write) table -> external                    (complete)
*  fnio.f (contains 5 routines for i/o)
*
* uses commonblock definitions in: fnXNNNN.inc
*
* needs the following routines from "fn-interface.f":
*     FNALPHAS        alpha_s interface (double precision function)
*     FNPDF           PDF interface
*
*******************************************************************
*******************************************************************



*******************************************************************
      Subroutine FX9999CC(FILENAME,XMUR,XMUF,IPRINTFLAG,XSECT)
*-----------------------------------------------------------------
* fastNLO user code v2.0 - main routine 
*
* input:
*   FILENAME    name of input table
*   XMUR        prefactor for nominal renormalization scale     
*                    any choice is possible, but please note 
*                    that 2-loop threshold corrections work
*                    only for xmur=xmuf
*   XMUF        prefactor for nominal fact-scale
*                     only a few choices are possible
*                     (see output or table documentation)
*   IPRINTFALG  =1 print results / =0 don't print
*
* output:
*   XSECT(nbin,3)    array of cross sections 
*                      first dimension: Bin number
*                      second dimension: 1 LO   2 NLO   
*                                3 NNLO / or other correction (if available)
*
*
* --- propose add flag: which order 1 LO, 2 NLO,
*
*
* MW 04/15/2005 initial version
* MW 09/02/2005 implement flexible scale variations
* TK 12/07/2005 table format contains now LO + NLO with 5 scale variations
* MW 2006/01/17 implement tableformat version 1c
* MW 2006/02/01 implement tableformat version 1.4 
* TK 2006/08/06 implement scalebins for N>2
* MW 2006/08/09 add normalization feature: 1/sigma dsigma/d[s.th.] 
* MW 2007/06/11 (@London Gatwick) implement v2.0
* KR 2009/09/22 Add filling sum of subprocesses (into 0 bin), sum of all
*               orders for isub=0 (into 0 bin) and adding up contributions
*               successively from lower to higher order according to pointers.
*               ===> result(MxObsBin,0:MxSubproc,0:MxCtrb)
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Integer IFILE, IPoint, IScPoint, I,J,K,L,M, 
     +     IPrintFlag,
     +     maxscale, nbin,nx,
     +     lenocc
      Character*(*) FILENAME
      Double Precision Xmur, Xmuf

c === Initialization: Read table
      Call FX9999IN(Filename)

c === Reset output array
      Do i=0,NContrib
         Do j=1,NObsBin
            xsect(j,i) = 0.d0
         Enddo
      Enddo

c === Determine access order to contributions and set pointers to 
c === contributions/scales
      Call FX9999PT(xmur,xmuf)

c === Loop over contributions, use pointers ordered according to FX9999PT
      Do i=1,NContrib
         IPoint   = IContrPointer(i)
         IScPoint = IScalePointer(i)
cdebug
cdebug         Write(*,*)"FX9999CC: IContr, IContrPointer, "//
cdebug     +        " IScalePointer, NSubproc(IContrPointer): ",
cdebug     +        i, IPoint, IScPoint, NSubProc(IPoint)
cdebug

c === Get PDFs
         Call FX9999GP(i)       ! Use i, pointers are derived inside. Change?
         
c === Multiply with perturbative coefficients and alphas
         Call FX9999MT(i,xmur,xmuf) ! Use i, pointers are derived inside. ?


c === Add up in result array
         Do j=1,NObsBin
            result(j,0,IPoint) = 0.D0
            Do k=1,NSubProc(IPoint)
               result(j,0,IPoint) = result(j,0,IPoint) +
     >              result(j,k,IPoint)
            Enddo
         Enddo
      Enddo

c === Fill sum of subprocesses, --> result(j,0,i), and 
c === sum of sums, --> result(j,0,0)      
      Do j=1,NObsBin
         If (NContrib.gt.1) Then
            Do i=2,NContrib
               result(j,0,IContrPointer(i)) =
     >              result(j,0,IContrPointer(i)) +
     >              result(j,0,IContrPointer(i-1))
            Enddo
         Endif
         result(j,0,0) = result(j,0,NContrib)
         Do i=0,NContrib
            xsect(j,i) = result(j,0,i)
cdebug            WRITE(*,*)"FX9999CC: IOBS,ICONTR,XSECT: ",
cdebug     >           j,i,xsect(j,i)
         Enddo
      Enddo

c === Normalization: Todo
      If (INormFlag.eq.0) Then
         Continue               ! no normalization - nothing to do
      ElseIf (INormFlag.eq.1) Then
c         Call FX9999NM          ! normalize by own integral
         continue
c      ElseIf (INormFlag.eq.2 .or INormFlag.eq.3) Then ! get denomin., divide
c         Call DX9999CC(DenomTable,Xmur,Xmuf,0,Xsect2)
c         Do i=1,NObsBin
c            Sum = 0d0 
c            Do j=IDivLoPointer(i),IDivUpPointer
c            Sum=Sum+Xsect2(j)
c         Enddo
c         Xsect(i) = Xsect(i)/Sum
c         EndDo
      EndIf

c === Print results - if requested
      If (IPrintFlag.eq.1) Call FX9999PR(xsect)
      Return

c 5000 Format (A,A64)
 5001 Format (A,A,A)
 5002 Format (A,F9.4,4X,A,F9.4)
 998  Continue
      End

*******************************************************************
*******************************************************************
      SUBROUTINE FX9999IN(Filename)
*-----------------------------------------------------------------
* MW 06/10/2007
*
* initialize fastNLO code
* 
* input:    Filename of fastNLO table
*
* 06/11/2007 MW
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Include 'strings.inc'
      Integer i,j,k,
     +     lenocc
      Character*(*) FILENAME
      Character*255 OLDFILENAME
      Data OLDFILENAME/'xxxx'/
      Save OLDFILENAME

c === reset result arrays
ckr Add back Ctrb loop for xsect as in description
      Do i=1,MxObsBin
         Xsect(i,0) = 0d0
         Xsect2(i,0)= 0d0
         Do k=1,MxCtrb
            Xsect(i,k) = 0d0
            Xsect2(i,k)= 0d0
         Enddo
         Do k=0,MxCtrb
            Do j=0,MxSubproc
               result(i,j,k) = 0d0
            Enddo
         Enddo
      Enddo

c === output in first fastNLO call
      If (IFNfirst.eq.0) Then
         Do i=1,12
            Write(*,5000) " # ",CHEADER(i)
         Enddo
         IFNfirst = 1
      Endif

c === in 1st scenario call: read fastNLO coefficient table
      If (FILENAME(1:LENOCC(FILENAME)).ne.
     >     OLDFILENAME(1:LENOCC(OLDFILENAME))) Then
c         Call FX9999RD(FILENAME)
         Call FX9999RW('read',FILENAME(1:LENOCC(FILENAME)))
                                ! new flexible version
         OLDFILENAME = FILENAME

c   - check consistency of array dimensions / commonblock parameters
c ----------> to be done
      Endif


 5000 Format (A,A64)
      Return 
      End

*******************************************************************
*******************************************************************
      SUBROUTINE FX9999PT(xmur,xmuf)
*-----------------------------------------------------------------
* MW 06/14/2007
*
* determine pointer to contributions and to scales
* 
* input:    
*           Xmur        renormalization scale factor        
*           Xmuf        factorization scale factor
*
* 06/11/2007 MW
* 19.09.2009 KR Improve contribution determination
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Integer i,j,k, i1
      Double Precision Xmur, Xmuf, XmurOld, XmufOld
      Data XMurOld/0d0/,XmufOld/0d0/
      
c --- Find particular contributions
      IContr = 0
      IContrPointer(IContr) = 0

c --- Find LO contribution
      If (PORDPTHY.ge.1) then
         j = IContr + 1
         IContrPointer(j) = -1
         Do i=1,NContrib
            If (IContrFlag1(i).eq.1.and.IContrFlag2(i).eq.1) Then
               IContr = IContr+1
               IContrPointer(IContr) = i
            Endif
         Enddo
         If (IContrPointer(j).eq.-1) Then
            Write(*,*)"FX9999PT: ERROR! Requested contribution "//
     >           "not available, stopped!"
            Write(*,*)"          PORDPTHY = ",PORDPTHY
            Stop
         Elseif (j.ne.IContr) Then
            Write(*,*)"FX9999PT: ERROR! More than one contribution "//
     >           "available, stopped!"
            Write(*,*)"          PORDPTHY = ",PORDPTHY
            Write(*,*)"          j = ",j,"IContr = ",IContr
            Stop
         Endif
      Endif
      
c --- Find NLO contribution
      If (PORDPTHY.ge.2) then
         j = IContr + 1
         IContrPointer(j) = -1
         Do i=1,NContrib
            If (IContrFlag1(i).eq.1.and.IContrFlag2(i).eq.2) Then
               IContr = IContr+1
               IContrPointer(IContr) = i
            Endif
         Enddo
         If (IContrPointer(j).eq.-1) Then
            Write(*,*)"FX9999PT: ERROR! Requested contribution "//
     >           "not available, stopped!"
            Write(*,*)"          PORDPTHY = ",PORDPTHY
            Stop
         Elseif (j.ne.IContr) Then
            Write(*,*)"FX9999PT: ERROR! More than one contribution "//
     >           "available, stopped!"
            Write(*,*)"          PORDPTHY = ",PORDPTHY
            Write(*,*)"          j = ",j,"IContr = ",IContr
            Stop
         Endif
      Endif
      
c --- Find 1-loop TC
      If (PTHRESHCOR.eq.1) then
         If (PORDPTHY.lt.1) then
            Write(*,*)"FX9999PT: ERROR! Inconsistent choice of "//
     >           "1-loop TC without LO, stopped!"
            Write(*,*)"          PORDPTHY = ",PORDPTHY,", PTHRESHCOR = "
     >           ,PTHRESHCOR
            Stop
         Endif
         If (PORDPTHY.ge.2) then
            Write(*,*)"FX9999PT: ERROR! Inconsistent choice of "//
     >           "1-loop TC with NLO, stopped!"
            Write(*,*)"          PORDPTHY = ",PORDPTHY,", PTHRESHCOR = "
     >           ,PTHRESHCOR
            Stop
         Endif
         j = IContr + 1
         IContrPointer(j) = -1
         Do i=1,NContrib
            If (IContrFlag1(i).eq.2.and.IContrFlag2(i).eq.1.and.
     >           IContrFlag3(i).eq.1) Then
               IContr = IContr+1
               IContrPointer(IContr) = i
            Endif
         Enddo
         If (IContrPointer(j).eq.-1) Then
            Write(*,*)"FX9999PT: ERROR! Requested contribution "//
     >           "not available, stopped!"
            Write(*,*)"          PORDPTHY = ",PORDPTHY,", PTHRESHCOR = "
     >           ,PTHRESHCOR
            Stop
         Elseif (j.ne.IContr) Then
            Write(*,*)"FX9999PT: ERROR! More than one contribution "//
     >           "available, stopped!"
            Write(*,*)"          PORDPTHY = ",PORDPTHY,", PTHRESHCOR = "
     >           ,PTHRESHCOR
            Write(*,*)"          j = ",j,"IContr = ",IContr
            Stop
         Endif
      Endif

c --- Find 2-loop TC
      If (PTHRESHCOR.eq.2) then
         If (PORDPTHY.lt.2) then
            Write(*,*)"FX9999PT: ERROR! Inconsistent choice of "//
     >           "2-loop TC without NLO, stopped!"
            Write(*,*)"          PORDPTHY = ",PORDPTHY,", PTHRESHCOR = "
     >           ,PTHRESHCOR
            Stop
         Endif
         If (PORDPTHY.ge.3) then
            Write(*,*)"FX9999PT: ERROR! Inconsistent choice of "//
     >           "2-loop TC with NNLO, stopped!"
            Write(*,*)"          PORDPTHY = ",PORDPTHY,", PTHRESHCOR = "
     >           ,PTHRESHCOR
            Stop
         Endif
         j = IContr + 1
         IContrPointer(j) = -1
         Do i=1,NContrib
            If (IContrFlag1(i).eq.2.and.IContrFlag2(i).eq.1.and.
     >           IContrFlag3(i).eq.2) Then
               IContr = IContr+1
               IContrPointer(IContr) = i
            Endif
         Enddo
         If (IContrPointer(j).eq.-1) Then
            Write(*,*)"FX9999PT: ERROR! Requested contribution "//
     >           "not available, stopped!"
            Write(*,*)"          PORDPTHY = ",PORDPTHY,", PTHRESHCOR = "
     >           ,PTHRESHCOR
            Stop
         Elseif (j.ne.IContr) Then
            Write(*,*)"FX9999PT: ERROR! More than one contribution "//
     >           "available, stopped!"
            Write(*,*)"          PORDPTHY = ",PORDPTHY,", PTHRESHCOR = "
     >           ,PTHRESHCOR
            Write(*,*)"          j = ",j,"IContr = ",IContr
            Stop
         Endif
      Endif

c --- Check availability of fact scale choice and assign pointer
c     (based on fact scale since renorm scale is flexible)
c           --- except for threshold corrections ---
      Do i=1,IContr
         i1 = IContrPointer(i)
         IScalePointer(i) = 0
         If (IScaleDep(i1).eq.0) Then ! Born-type w/o scale dep - use any scale
            If (NScaleVar(i1,1).ge.1) Then
               IScalePointer(i) = 1
            Else
               Write(*,*)"FX9999PT: ERROR! Not a single scale "//
     >              "available in contribution, stopped!"
               Write(*,*)"          IContr = ",IContr,
     >              ", NScaleVar(.,1) = ",NScaleVar(i1,1) 
               Stop
            Endif
            If (NScaleVar(i1,1).gt.1) Then
               Write(*,*)"FX9999PT: WARNING! Why more than one "//
     >              " scale variation for Born-type contribution?"
            Endif
         Else                   ! no Born type contribution
            Do j=1,NScaleVar(i1,1)
               If (xmuf.eq.scalefac(i1,1,j)) IScalePointer(i)=j
            Enddo
         Endif
         If (IScalePointer(i).eq.0) Then
            Write(*,*)"FX9999PT: ERROR! The requested "//
     >           " factorization scale xmuf = ",xmuf
            Write(*,*)"          is not available, stopped!"
            Do j=1,NContrDescr(i1)
               Write(*,*) "         Descriptions: ",CtrbDescript(i1,j)
            Enddo
            Stop
         Endif 
      Enddo

c --- Check if renormalization scale can be provided a posteriori if needed
      If (xmur .ne. xmuf) Then
         Do i=1,IContr
            i1 = IContrPointer(i)
            IF (IScaleDep(i1).eq.2) Then
               Write(*,*)"FX9999PT: ERROR! The requested "//
     >              " renormalization scale xmur = ",xmur
               Write(*,*)"          is not available, stopped!"
               Write(*,*)"          Only xmur=xmuf is possible."
               Do j=1,NContrDescr(i1)
                  Write(*,*) ' ',CtrbDescript(i1,j)
               Enddo
               Stop
            Endif
         Enddo
      Endif

c --- Debug print-out
c      Do i=1,IContr
c         Write(*,*) "FX9999PT: Pointer number ",i," to "//
c     >        "IContr, IScale:",IcontrPoint(i),IScalePointer(i)
c      Enddo


      Return
      End

*******************************************************************
*******************************************************************
      SUBROUTINE FX9999MT(in,xmur,xmuf)
*-----------------------------------------------------------------
* MW 06/10/2007
*
* multiply the PDFs and the perturbative coefficients 
*
*  question: give 'relative' contribution
*            IScalePointer only available for relative contrib.
* 
* input:    IN    No. of contribution in present calculation
*           XMUR  prefactor for nominal renormalization scale
*           XMUF  prefactor of nominal factorization scale setting
*
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Integer Ixmur,Ixmuf, in,ic,is, i,j,k,l,m,n, nxmax
      Double Precision xmur,xmuf,mu, FNALPHAS, scf
      Double Precision logmu, scfac,scfac2a, scfac2b, as(8),aspow(8)
      Double Precision pi, beta0, beta1, NF,CA,CF
      Parameter (PI=3.14159265358979323846, NF=5d0, CA=3d0, CF=4d0/3d0)
      Parameter (beta0=(11d0*CA-2d0*NF)/3d0) 
      Parameter (beta1=34*CA*CA/3d0-2d0*NF*(CF+5d0*CA/3d0))

c - set pointers to contribution in table and to scale variation
      ic = IContrPointer(in)
      is = IScalePointer(in)

c - the absolute order in alpha_s of the LO contribution is in ILOord
c - vary renormalization scale around the value used in orig. calculation
c      logmu = log(xmur/xmuf)    ! aposteriori change w.r.t. orig. calculation
      logmu = log(xmur/ScaleFac(ic,1,is)) ! aposteriori change wrt. orig. calc.
      scfac  = dble(ILOord)  *beta0 *logmu          ! NLO contrib.
      scfac2a= dble(ILOord+1)*beta0 *logmu          ! NNLO contrib.
      scfac2b= dble(ILOord*(ILOord+1))/2d0*beta0*beta0*logmu*logmu  
     +     + dble(ILOord)*beta1/2d0*logmu           ! NNLO contrib. continued

c - MW:  maybe simpler if we make the mur-variation later
c        for the whole contribution - instead of doing it for
c        each array element??



c - MW: strategy: - find scalevar for which muf corresponds to selected muf
c                 - check mur-factor in this scalevar
c                 - compute needed mur-variation = xmur/mur-factor
c



c - loop over coefficient array - compar with order in table storage!!!!!!
c loop: observable, scalebins,(get alphas), xbins,subproc
      Do j=1,NObsBin

         IF (NPDFdim(ic).eq.0) Then
            nxmax =  NxTot(ic,1,j) ! 1d case (DIS)
         Elseif (NPDFdim(ic).eq.1) Then
            nxmax = (NxTot(ic,1,j)*NxTot(ic,1,j)+NxTot(ic,1,j))/2 ! 2d half mat
         Elseif (NPDFdim(ic).eq.2) Then
            Write(*,*) ' fx9999mt:   2d case not yet'
            Stop
         Endif
         
         Do k=1,NScaleNode(ic,1)
            If (IScaleDep(ic).eq.0) Then ! no scale dep of coefficients
               mu = ScaleNode(ic,j,1,1,k)
            Else
               mu = ScaleNode(ic,j,1,is,k) ! ???????
            Endif
c     - get alphas
c            write(*,*) 'xmur ',xmur,mu
            as(1) =  FNALPHAS(xmur * mu) ! ??????????????
            aspow(1) = as(1)**Npow(ic)
            Do l=1,nxmax
               Do m=1,NSubProc(ic)
                  result(j,m,ic) = result(j,m,ic) + 
     +                 SigmaTilde(ic,j,1,is,k,l,m)
     +                 * aspow(1)
     +                 * pdf(j,k,l,m)
               Enddo
            Enddo
         Enddo
      Enddo

      Return 
      End

*******************************************************************
*******************************************************************
      SUBROUTINE FX9999GP(in)
*-----------------------------------------------------------------
* MW 08/26/2005
*
* read the PDFs in all (rap,pT) bins - at all xmax, xmin bins
* the default factorization scale (in GeV) is multiplied by muffactor
*
* input: IC    number of contribution
*
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Double Precision Xmuf, x, muf, 
     +     tmppdf(-6:6), xpdf1(MxNxTot,-6:6),xpdf2(MxNxTot,-6:6), H(10)
      Integer in,ic,is, i,j,k,l,m, nx, nx2limit

c - set pointers to contribution in table and to scale variation
      ic = IContrPointer(in)
      is = IScalePointer(in)

      Do i=1,NObsBin
         Do j=1,NScaleNode(ic,1)
            muf = ScaleNode(ic,i,1,is,j)
c            write(*,*) 'muf ',muf
            nx =0
            Do k=1,NxTot(ic,1,i)  ! --- fill first PDF
               x = Xnode1(ic,i,k) 
               Call FNPDF(x, muf, tmppdf)
               Do m=-6,6
                  xpdf1(k,m) = tmppdf(m)
                  If (NPDF(ic).eq.1) Then ! DIS
                     Continue   
                  Elseif (NPDF(ic).eq.2) Then ! two hadrons
                     If (NPDFPDG(ic,1).eq.NPDFPDG(ic,2)) Then ! identical
                        xpdf2(k,m) = tmppdf(m)
ckr                        write(*,*)"k,m,xpdf2",k,m,tmppdf(m)
                     Elseif (NPDFPDG(ic,1).eq.-NPDFPDG(ic,2)) Then ! h anti-h
                        xpdf2(k,-m) = tmppdf(m)
                     Else
                        Write(*,*) ' So far only the scattering of identical'
                        Write(*,*) ' hadrons or hadron&anti-hadron is'
                        Write(*,*) ' implemented -> gamma-p to be done...'
                        Stop
                     Endif
                  Else
                     write(*,*) ' neither one nor two hadrons...?',NPDF(ic)
                     Stop
                  Endif
               Enddo            ! m
            Enddo               ! k 
            If (NPDFdim(ic).eq.2) then ! --- fill second PDF
               Do k=1,NxTot(ic,2,i)
                  x = Xnode2(ic,i,l)
                  Call FNPDF(x, muf, tmppdf) ! < to be changed ->different PDF!
                  Do m=-6,6
                     xpdf2(j,m) = tmppdf(m)
                  Enddo
               Enddo
            Endif
            Do k=1,NxTot(ic,1,i) ! --- build PDF1,2 linear combinations
               If (NPDF(ic).eq.1) nx2limit = 1 !               1d case (DIS)
               If (NPDFdim(ic).eq.1) nx2limit = k !            2d half matrix
               If (NPDFdim(ic).eq.2) nx2limit = NxTot(ic,2,i) !  2d full matrix
               Do l=1,nx2limit
                  nx = nx+1
                  Call fx9999pl(IPDFdef(ic,1),IPDFdef(ic,2),
     +                 IPDFdef(ic,3),k,l,xpdf1,xpdf2,H)
                  Do m=1,NSubproc(ic)
                     pdf(i,j,nx,m) = H(m)
ckr                     write(*,*)"i,j,nx,m,pdf",i,j,nx,m,pdf(i,j,nx,m)
                  Enddo         ! m subproc
               Enddo            ! x2-loop
            Enddo               ! x1-loop
         Enddo                  ! ScaleNode-loop
      Enddo                     ! ObsBin-loop

      Return
      End

*******************************************************************
*******************************************************************
      SUBROUTINE FX9999PL(icf1,icf2,icf3,i,j,XPDF1,XPDF2,H)
* ---------------------------------------------------------------
* MW 06/13/07
* compute PDF linear combinations - for different sets of subprocesses
*
* depending on icf1,2,3, the product of the i-th and j-th entries
* of the PDF array XPDF are multiplied into the relevant linear 
* combinations in the array H
*
* input:
*    ireact             flag for reaction (1:DIS, 2:pp-jets, 3:ppbar-jets)
*    i                  x-index of first hadron        
*    j                  x-index for second hadron (if two-hadron process)
*    XPDF1(nxmax,-6:6)  PDF array for all x-bins
*    XPDF2(nxmax,-6:6)  PDF array for all x-bins
*
* output:
*    H(10)              PDF linear combinations
* ---------------------------------------------------------------
      Implicit None
      INCLUDE 'fnx9999.inc'
      Integer icf1,icf2,icf3, i,j,k
      Double Precision XPDF1(MxNxTot,-6:6),XPDF2(MxNxTot,-6:6), H(10),
     +     G1, G2,              ! gluon densities from both hadrons
     +     SumQ1, SumQ2,        ! sum of quark densities
     +     SumQB1, SumQB2,      ! sum of anti-quark densities
     +     Q1(6),Q2(6), QB1(6),QB2(6), ! arrays of 6 (anti-)quark densities
     +     S,A                  ! products S,A

c --- DIS: inclusive and jets, gammaP direct
      if (icf1.eq.2 .and. (icf2.eq.0 .or. icf2.eq.1)) then 
         H(3) = 0d0             ! Delta  at O(as^0)
         do k=1,5,2
            H(3) = H(3) + (XPDF1(i,k)+XPDF1(i,-k)+
     +           4d0*(XPDF1(i,k+1)+XPDF1(i,-k-1)))/9d0
         enddo
         H(1) = XPDF1(i,0)      ! Gluon  at O(as^1)
         H(2) = 0d0             ! Sigma  at O(as^2)
         do k=1,6
            H(2) = H(2)+XPDF1(i,k)+XPDF1(i,-k)
         enddo

c --- hadron-hadron: jets
      elseif (icf1.eq.3.and.icf2.eq.1.and.(icf3.ge.1.and.icf3.le.2))then 
         SumQ1  = 0d0
         SumQB1 = 0d0
         SumQ2  = 0d0
         SumQB2 = 0d0
         do k=1,6
            Q1(k)  = XPDF1(i,k)  ! read 1st PDF at x1
            QB1(k) = XPDF1(i,-k)
            SumQ1  = SumQ1  + Q1(k)
            SumQB1 = SumQB1 + QB1(k)
            Q2(k)  = XPDF2(j,k)  ! read 2nd PDF at x2
            QB2(k) = XPDF2(j,-k)
            SumQ2  = SumQ2  + Q2(k)
            SumQB2 = SumQB2 + QB2(k)
         enddo
         G1     = XPDF1(i,0)
         G2     = XPDF1(j,0)
c   - compute S,A
         S = 0d0
         A = 0d0
         do k=1,6
            S = S + (Q1(k)*Q2(k)) + (QB1(k)*QB2(k)) 
            A = A + (Q1(k)*QB2(k)) + (QB1(k)*Q2(k)) 
         enddo
c   - compute seven combinations
         H(1) = G1*G2
         H(2) = SumQ1*SumQ2 + SumQB1*SumQB2 - S
         H(3) = S
         H(4) = A
         H(5) = SumQ1*SumQB2 + SumQB1*SumQ2 - A
         H(6) = (SumQ1+SumQB1)*G2
         H(7) = G1*(SumQ2+SumQB2)
         if (icf3.eq.1) H(6) = H(6)+H(7) ! case: 6 subproc
c --- gammaP: direct, jets
      elseif (icf1.eq.2 .and. icf2.eq.2) then 
         write(*,*) '    gammaP to be implemented'
         stop

      else
         write(*,*) '    icf1,2,3 =',icf1,icf2,icf3
         write(*,*) '    this combination is not yet defined'
         stop
      endif

      RETURN 
      END
*******************************************************************


*******************************************************************
*******************************************************************
      SUBROUTINE FX9999PR(xsect)
*-----------------------------------------------------------------
* M. Wobisch  - print fastNLO cross section results
*
* input:
*   xsect       result array
*
* MW 04/18/2005
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Integer i,j

      write(*,*) "fastNLO results:"

ckr Add back Ctrb loop for xsect as in description
      Do j=0,NContrib
         Do i=1,NObsBin
            write(*,*) 'icontr, bin No, result ',j,i,'  ',xsect(i,j)
         Enddo
      Enddo

      Return
      End

*******************************************************************
*******************************************************************
      Subroutine FX9999NF
*-----------------------------------------------------------------
* fastNLO user code v2.0 - print scenario iNFormation
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Integer i,j

      Write(*,*)
      Write(*,*)' ######################################################'
      Write(*,*)' #    information on fastNLO scenario: ',ScenName
      Write(*,*)' #   ------------------------------------------'
      Write(*,*)' # description:'
      Do i=1,NScDescript
         Write(*,*)' #   ',SCDescript(i)  
      Enddo
      Write(*,*)' # measured at at Ecms =',Ecms,' GeV'
      Write(*,*)' #'
      Write(*,*)' # tot. No. of Observable bins: ',NObsBin,' in ',Ndim,' dimensions'
      Do i=1,NDim
         Write(*,*)' # dimension ',i,':'
         Write(*,*)' #        ',DimLabel(i)
      Enddo
      Write(*,*)' #'
      Write(*,*)' # No. of contributions: ',Ncontrib
      Do i=1,Ncontrib         
         Write(*,*)' # - contribution',i,':'
         Do j=1,NcontrDescr(i)
            Write(*,*)' #     ',CtrbDescript(i,j)
         Enddo
         Write(*,*)' #     computed by:'
         Do j=1,NcodeDescr(i)
            Write(*,*)' #     ',CodeDescript(i,j)
         Enddo
         
      Enddo
      Write(*,*)' # No. of x bins in 1st contrib.:',Nxtot(1,1,1)
      Write(*,*)' # No. of available scale variations/scale nodes:'
      Do i=1,Ncontrib         
         Write(*,*)' #   NscaleVar,NScaleNode:',NscaleVar(i,1),NScaleNode(i,1)
      eNDDO
      Write(*,*)' #'
      Write(*,*)' ######################################################'

      Return
      End

*******************************************************************
