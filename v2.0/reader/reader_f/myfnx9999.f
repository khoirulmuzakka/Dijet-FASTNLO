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
*                 - e.g. threshcore (no aposteriori mur): multiple same
*                   muf variations for different mur variations
*                   not sufficient to find first muf var. that works 
*  fx9999pm     get PDFs, multiply (does aposteriori mur variation)
*  fx9999gp     get PDFs      - missing details                         (works)
*                >>> currently under investigation
*
*  fx9999pl     compute PDF linear combinations (done - cosmetics)      (works)
*  fx9999mt     multiply coefficients and PDFs                       (advanced)
*  fx9999pr     print results                                           (works)
*  fx9999nf     print scenario information (physics & technical)        (works)
*  fx9999nm     normalize distribution by its own integral              (to do)
*
*  fx9999tb     returns information on variables from the table       (at work)
*
* external:
*  fx9999rw     read (or write) table -> external                    (complete)
*               -> need to include table consistency checks (<max?)
*
* global (scenario independent) routines:
*  fnio.f   (contains 5 routines for i/o)
*  fnset.f  (set flags which define contributions to be included)
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
*   XSECT(nbin)    array of cross sections 
*                      nbin (int): Bin number
*
*
* --- propose add flag: which order 1 LO, 2 NLO,
*
* current restrictions: (see also subroutines)
*  - can deal with single scale dim only (no DIS with muf=Q, mur=pT
*  - assumes xmur=xmuf for each contribution
*  - assumes same cross section units (nb,pb,...) for all contributions
*        -> maybe convert everything to"published units"?
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
     +     maxscale, nbin,nx
      Character*(*) FILENAME
      Double Precision Xmur, Xmuf

c === Initialization: Read table
      Call FX9999IN(Filename)

c === Reset output array
      Do i=1,NObsBin
         xsect(i) = 0d0
      Enddo
c === Determine pointers to acces order of contributions and to scales
c       >>> this should also fill new pointers to LO and NLO, needed
c       >>> for aposteriori scale variations    
      Call FX9999PT(xmur,xmuf)

c === check if aposteriori mur variation is required
c     Imurapost = 1
c     set Log term
c     fill beta0, n

c === see below: should we combine calls to GP and MT in single subroutine?
c     -> sometimes the sequence may be called twice

c === Loop over contributions, use pointers ordered according to FX9999PT
      Do i=1,Icontr
         IPoint   = IContrPointer(i)
ckr         Write(*,*)"FX9999CC: IContr, IContrPointer, "//
ckr     +        " IScalePointer, NSubproc(IContrPointer): "
ckr         Write(*,*)'   ',i, IPoint, IScalePointer(i), NSubProc(IPoint)

c - Get PDFs - Multiply with perturbative coefficients and alphas 
         Call FX9999PM(i,xmur,xmuf)

      Enddo


c - add results in output array - does not work for NP/npert cor
      Do i=1,Icontr
c      Do i=1,1   ! test - only single contribution  LO
c      Do i=2,2   ! test - only single contribution  NLO
c      Do i=3,3   ! test - only single contribution  threshcor
         IPoint   = IContrPointer(i)
         Do j=1,NObsBin
            Do k=1,NSubProc(Ipoint)
c            Do k=1,1 ! test - only gg 
c            Do k=2,2 ! test - only g (DIS)
c            Do k=2,5 ! test - only qq
               xsect(j) = xsect(j)+result(j,k,i)
            Enddo
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
c            Xsect(i) = Xsect(i)/Xsect2(IDivPointer(i))
cc            Sum = 0d0 
cc            Do j=IDivLoPointer(i),IDivUpPointer
cc              Sum=Sum+Xsect2(j)
cc            Enddo
cc         Xsect(i) = Xsect(i)/Sum
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
* current restrictions:
*     check consistency of array sizes and commonblock definitions
*        -> or do this in RW routine?? (preferred!)
*
* 06/11/2007 MW
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Include 'strings.inc'
      Integer i,j,k, ictrb(30)
      Character*(*) FILENAME
      Character*255 OLDFILENAME
      Data OLDFILENAME/'xxxx'/
      Save OLDFILENAME

c === reset result arrays
      Do i=1,MxObsBin
         Xsect(i) = 0d0
         Xsect2(i)= 0d0
         Do j=0,MxSubproc
            Do k=0,MxCtrb
               result(i,j,k) = 0d0
            Enddo
         Enddo
      Enddo

c === output in first fastNLO call
      If (IFNfirst.eq.0) Then
         Do i=1,13
            Write(*,5000)" #"//CHEADER(i)
         Enddo
         IFNfirst = 1
      Endif

c === in 1st scenario call: read fastNLO coefficient table
      If (FILENAME.ne.OLDFILENAME) Then
         Call FX9999RW('read',FILENAME,ictrb)
         OLDFILENAME = FILENAME

c   - check consistency of array dimensions / commonblock parameters
c ----------> to be done - or do this in i/o routine "RW"?
c --> better in RW so it can be used in other codes / in read and write
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
* 19.09.2009 KR Improve contribution determination, fix double comparison
*
* current restrictions:
* - muf and mur recognition:
*   So far it is assumed that all scales are stored in first dimension
*   multiple scale dimensions are not yet implemented here
* - mur recognition
*   So far it is assumed that (in one dimension) the scale factors
*   for mur and muf are the same (i.e. can not have fixed mur, but different
*   muf values (which, in principle would work with aposteriori variations)
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Integer i,j,k, i1
      Double Precision Xmur, Xmuf
      
c --- Check input
      If (xmur.lt.1d-3.or.xmuf.lt.1d-3) Then
         WRITE(*,*)"FX9999PT: ERROR! Scale factors smaller than "//
     >        "0.001 are not allowed, stopped! xmur, xmuf = ",xmur,xmuf
         STOP
      Endif

c --- Find particular contributions
      IContr = 0

c --- Find LO contribution
      If (PORDPTHY.ge.1) then
         j = IContr + 1
         IContrPointer(j) = -1
         Do i=1,NContrib
ckr            Write(*,*) 'FX9999PT Icintrflag1,2 ',IContrFlag1(i),IContrFlag2(i),Iref(i)
            If (IContrFlag1(i).eq.1.and.IContrFlag2(i).eq.1
     +           .and. Iref(i).eq.Preftab) Then 
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
c - would be strange, but why stop?
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
            If (IContrFlag1(i).eq.1.and.IContrFlag2(i).eq.2
     +           .and. Iref(i).eq.Preftab) Then
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
            Write(*,*) IContrFlag1(i),IContrFlag2(i),
     +           IContrFlag3(i),Iref(i),Preftab
            If (IContrFlag1(i).eq.2.and.IContrFlag2(i).eq.1.and.
     >           IContrFlag3(i).eq.1
     +           .and. Iref(i).eq.Preftab) Then
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
     >           IContrFlag3(i).eq.2
     +           .and. Iref(i).eq.Preftab) Then 
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

c --- Check availability of factorization scale choice and assign pointer
c     (based on factorization scale since renorm scale is flexible -
c      except for threshold corrections)
      Do i=1,IContr
         i1 = IContrPointer(i)
         IScalePointer(i) = 0
ckr
ckr         write(*,*)"AAA: ic,icp,iscaledep",i,i1,iscaledep(i1)
         If (IScaleDep(i1).eq.0) Then ! Born-type w/o scale dep - use any scale
            If (NScaleVar(i1,1).ge.1) Then
               IScalePointer(i) = 1 ! use 1st scale variation
            Else
               Write(*,*)"FX9999PT: ERROR! Not a single scale "//
     >              "available in contribution, stopped!"
               Write(*,*)"          IContr = ",IContr,
     >              ", NScaleVar(.,1) = ",NScaleVar(i1,1) 
               Stop
            Endif
            If (NScaleVar(i1,1).gt.1) Then
               Write(*,*)"FX9999PT: WARNING! Why more than one "//
     >              " scale variation for contribution"
               Write(*,*) "        with scale-independent coefficients?"
            Endif
         Else                   ! >>> to do: distinguish between IScaleDep=2,3
            Do j=1,NScaleVar(i1,1)
               If (dabs(scalefac(i1,1,j)/xmuf-1d0).lt.1d-4) Then
c - works only if IScaleDep=2 - for =3 need to find correct pair of (xmur,xmuf)
                  IScalePointer(i)=j
ckr                  Write(*,*)"iscvar,xmuf,scalfac,scalpt",j,xmuf,
ckr     >                 scalefac(i1,1,j),IScalePointer(i)
                  Exit          ! <<< what does "Exit" do?
               Endif
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

c --- Check if renormalization scale is directly available or (if not)
c     can be provided a posteriori
c
c - the current treatment works only for a single scale dimension
c   in two scale dimensions one needs to check
c    - if the chosen xmur corresponds to the ... variable for the chosen muf
c          -> can be computed 
c    - if the chosen xmur corresponds to one of the ... variables
c      for the chosen muf -> can be computed only if aposteriori is possible
c
c - for 2 scale dim: need 2nd scale pointer I2ScalePointer
c
c >>>>>>>
c > should not be necessary if we distinguish between IScaleDep=2,3 in mur code
c
      If (dabs(xmur/xmuf-1.d0).gt.1d-4) Then
         Do i=1,IContr
            i1 = IContrPointer(i)
            If (IScaleDep(i1).eq.2) Then
               Write(*,*)"FX9999PT: ERROR! The requested "//
     >              " renormalization scale xmur = ",xmur
               Write(*,*)"          is not available, stopped!"
               Write(*,*)"          Only xmur=xmuf is possible."
               Do j=1,NContrDescr(i1)
                  Write(*,*) '  ',CtrbDescript(i1,j)
               Enddo
               Stop
            Endif
         Enddo
      Endif

c --- Debug print-out
ckr      Write(*,*) "FX9999PT: No. contributions selected ",Icontr
      Do i=1,IContr
ckr         Write(*,*) "FX9999PT: Pointer number ",i," to "//
ckr     >        "IContr, IScale:",IcontrPointer(i),IScalePointer(i)
      Enddo


      Return
      End

*******************************************************************
*******************************************************************
      Subroutine FX9999PM(ictrb,xmur,xmuf)
*-----------------------------------------------------------------
* fastNLO user code v2.0 
*
*   - get PDFs & multiply with alphas and coefficients
*
* input:
*   Ictrb       Number of contribution (logical - not in table!)
*   XMUR        prefactor for nominal renormalization scale     
*                    any choice is possible, but please note 
*                    that 2-loop threshold corrections work
*                    only for xmur=xmuf
*   XMUF        prefactor for nominal fact-scale
*                     only a few choices are possible
*                     (see output or table documentation)
*
* current restrictions:
*
* MW 08/14/2010 initial version
* MW 08/24/2010 aposteriori mur variation implemented (NLO only)
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Integer Ictrb, Iaddpow, ic,is
      Double Precision Xmur, Xmuf
      Double Precision factor, beta0, beta1, NF,CA,CF, logmur
      Parameter (NF=5d0, CA=3d0, CF=4d0/3d0)
      Parameter (beta0=(11d0*CA-2d0*NF)/3d0) 
      Parameter (beta1=34*CA*CA/3d0-2d0*NF*(CF+5d0*CA/3d0))

ckr      Write(*,*) 
ckr      Write(*,*) 'FX9999PM at work ',ictrb,xmur,xmuf

c - how to identify that a posteriori mur variation is required
c - if: only one scale dim and same mur,muf scale factors used in table

c - aposteriori variation not required in the following cases
c    - if coefficient has no scale dependence
c    - if xmur=xmuf (later: if mur scalefactor in contribution is 
c                            different from xmur)   ScaleFac(j,i)


c - first: compute 'standard' contribution at chosen scale 
      Call FX9999GP(Ictrb,Xmuf)
      Call FX9999MT(Ictrb,Xmur,Xmuf,0,1d0) ! standard contribution

c - aposteriori mur variation
      ic = IContrPointer(Ictrb)
      is = IScalePointer(Ictrb)
ckr      Write(*,*) 'FX9999PM: ic,is: ',ic,is,xmur,ScaleFac(ic,1,is)

      If (IScaleDep(ic).ne.0 .and. dabs(xmur/xmuf-1d0).gt.1d-4) Then
         Write(*,*) 'FX9999PM: trying aposteriori mur contribution'
         If (IScaleDep(ic).ne.1) Then
            Write(*,*) ' a posteriori scale variation for contribution',
     +           ic,' not possible'
            Stop
         Endif

c        - abs. order in alphas:     "Npow(ic)"
c        - abs order of LO contrib:  "ILOord"
c NLO if (Npow-ILOord)=1      NNLO if (Npow(ic)-ILOord)=2

c --- if NLO 
         If ( (Npow(ic)-ILOord).eq.1) Then
            logmur = log(xmur/ScaleFac(ic,1,is))
            factor = dble(ILOord)*beta0*logmur ! n beta0 logmu

            Write(*,*) 'factor ',factor,ILOOrd,real(beta0),real(logmur)
            Call FX9999GP(Ictrb-1,Xmuf)
            Call FX9999MT(Ictrb-1,Xmur,Xmuf,1,factor) ! 1:mod NLO
c --- if NNLO
         Elseif ( (Npow(ic)-ILOord).eq.2) Then
            logmur = log(xmur/ScaleFac(ic,1,is))

            Write(*,*) ' a posteriori mur variation beyond NLO not yet ',
     +           'implemented'
            Stop
c            factor = ...dble(ILOord)*beta0*logmur ! n beta0 logmu
            Call FX9999GP(Ictrb-2,Xmuf)
            Call FX9999MT(Ictrb-2,Xmur,Xmuf,2,factor) ! 2: mod LO

c            factor = ...dble(ILOord)*beta0*logmur ! n beta0 logmu
            Call FX9999GP(Ictrb-1,Xmuf)
            Call FX9999MT(Ictrb-1,Xmur,Xmuf,1,factor) ! 1: mod NLO
         Else
            Write(*,*) ' a posteriori mur variation beyond NNLO not yet ',
     +           'implemented'
            Stop
         Endif
      Endif                     ! a posteriori mur variation

      Return
      End

*******************************************************************
*******************************************************************
      SUBROUTINE FX9999MT(in,xmur,xmuf,Iaddpow,factor)
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
*           Iaddpow add. power in alphas - for aposteriori mur variation
*                   also: increment in result-contribution index where 
*                         result is stored
*           factor   =!1 for a posteriori mur variations
*
* current restrictions:
*
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Integer Iaddpow
      Integer Ixmur,Ixmuf, in,ic,is, i,j,k,l,m,n, nxmax,
     +     istore
      Double Precision xmur,xmuf,factor, mur, FNALPHAS, scf
      Double Precision scfac,scfac2a, scfac2b, as,aspow
      Double Precision pi
      Parameter (PI=3.14159265358979323846)

c      Write(*,*)
c      Write(*,*)"FX9999MT: in,xmur,xmuf,Iaddpow,factor",in,real(xmur),
c     +     real(xmuf),Iaddpow,real(factor)

c - set pointers to contribution in table and to scale variation
      ic = IContrPointer(in)
      is = IScalePointer(in)

c - if regular contribution (i.e. no aposterior scale variation)
c   then store in slot number of logical contribution
c   otherwise increment by No. of add powers in alphas
      istore = in + Iaddpow

c - loop over coefficient array - compare with order in table storage!!!!!!
c loop: observable, scalebins,(get alphas), xbins,subproc
      Do j=1,NObsBin

         IF (NPDFdim(ic).eq.0) Then ! 1d case (DIS)
            nxmax =  NxTot(ic,1,j) 
         Elseif (NPDFdim(ic).eq.1) Then ! 2dim half-matrix
            nxmax = (NxTot(ic,1,j)*NxTot(ic,1,j)+NxTot(ic,1,j))/2 
         Elseif (NPDFdim(ic).eq.2) Then ! 2dim full matrix
            Write(*,*) ' fx9999mt:   2d case not yet'
            Stop
         Endif
         
         Do k=1,NScaleNode(ic,1)
            mur = xmur / ScaleFac(ic,1,is) * ScaleNode(ic,j,1,is,k) 

c --- get alphas
            as =  FNALPHAS(mur)
            aspow = as**(Npow(ic)+IAddPow)
cdebug
c            Write(*,*)'FX9999MT: ScaleNode',real(ScaleNode(ic,j,1,is,k))
c            Write(*,*)"FX9999MT: ScaleDep, xmur, mur, as, asp=",
c     >           IScaleDep(ic),real(xmur),real(mur),
c     >           real(as),real(aspow)
cdebug
            Do l=1,nxmax
               Do m=1,NSubProc(ic)
c                  coeff = SigmaTilde(ic,j,1,is,k,l,m)
                  If (Preftab.eq.0) Then
cMWold                     result(j,m,ic) = result(j,m,ic) + 
cMW new: result array gets first entries filled
                     result(j,m,istore) = result(j,m,istore) + 
c     +                    coeff
     +                    SigmaTilde(ic,j,1,is,k,l,m)
     +                    * aspow
     +                    * pdf(j,k,l,m)
     +                    * factor
                  Else          ! Reference Table
cMWold                     result(j,m,ic) = result(j,m,ic) + 
                     result(j,m,istore) = result(j,m,istore) + 
     +                    SigmaTilde(ic,j,1,is,k,l,m)
                  Endif
               Enddo
            Enddo
         Enddo
      Enddo

cdebug
c      Do j=1,NObsBin
c         Do m=1,NSubProc(ic)
c            Write(*,*)"FX9999MT: IOBS,IPROC,IORD,"//
c     >           "RESULT: ",
c     >           J,M,IC,RESULT(J,M,IC)
c         Enddo
c      Enddo
cdebug

      Return 
      End

*******************************************************************
*******************************************************************
      SUBROUTINE FX9999GP(ictrb,xmuf)
*-----------------------------------------------------------------
* MW 08/26/2005
*
* read the PDFs in all (rap,pT) bins - at all xmax, xmin bins
* the default factorization scale (in GeV) is multiplied by muffactor
*
* input: 
*     ictrb (integer)    number of logical contribution 
*     xmuf  (dble)       factorization scale variation factor
*
* current restrictions:
*  - not checked
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Double Precision Xmuf, x, muf, 
     +     tmppdf(-6:6), xpdf1(MxNxTot,-6:6),xpdf2(MxNxTot,-6:6), H(10)
      Integer ictrb, ic,is, i,j,k,l,m, nx, nx2limit

c - set pointers to contribution in table and to scale variation
      ic = IContrPointer(ictrb)
      is = IScalePointer(ictrb)
c      Write(*,*)'FX9999GP: ictrb,ic,is: ',ictrb,ic,is

      Do i=1,NObsBin
         Do j=1,NScaleNode(ic,1)
c            If (ic.eq.1) Then   ! <<< this is clearly insufficient!
c                                !     applies to all Born-type contrib.
c            Write(*,*) '    ic,IScaleDep(ic) ',ic,IScaleDep(ic)
            If (IScaleDep(ic).eq.0) Then
               muf = xmuf*ScaleNode(ic,i,1,is,j) ! 1: ScaleDim  is:ScaleVar
            Else
               muf = ScaleNode(ic,i,1,is,j)
            Endif

c            Write(*,*)"FX9999GP: xmuf, scalenode, muf: ",
c     >           real(xmuf),real(ScaleNode(ic,i,1,is,j)),real(muf)
            nx =0
            Do k=1,NxTot(ic,1,i)  ! --- fill first PDF
               x = Xnode1(ic,i,k) 
               Call FNPDF(x, muf, tmppdf)
               Do m=-6,6
                  xpdf1(k,m) = tmppdf(m)
c                  xpdf1(k,m) = tmppdf(m) / x  ! ???
                  If (NPDF(ic).eq.1) Then ! DIS
                     Continue   
                  Elseif (NPDF(ic).eq.2) Then ! two hadrons
                     If (NPDFPDG(ic,1).eq.NPDFPDG(ic,2)) Then ! identical hh
                        xpdf2(k,m) = tmppdf(m)
ckr                        write(*,*)"k,m,xpdf2",k,m,tmppdf(m)
                     Elseif (NPDFPDG(ic,1).eq.-NPDFPDG(ic,2)) Then ! hh-bar
                        xpdf2(k,-m) = tmppdf(m)
                     Else
                        Write(*,*)
     >                       ' So far only the scattering of identical'
                        Write(*,*) ' hadrons or hadron&anti-hadron is'
                        Write(*,*)
     >                       ' implemented -> gamma-p to be done...'
                        Stop
                     Endif
                  Else
                     Write(*,*) ' neither one nor two hadrons...?'
     >                    ,NPDF(ic)
                     Stop
                  Endif
               Enddo            ! m
            Enddo               ! k 
            If (NPDFdim(ic).eq.2) then ! --- fill second PDF
               Do k=1,NxTot(ic,2,i)
                  x = Xnode2(ic,i,l)
                  Call FNPDF(x, muf, tmppdf) ! < to be changed ->different PDF!
cMW      what does this mean 'different PDF'? - maybe for different hadrons?
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
*
* current restrictions:
*  - order of DIS subprocesses wrong
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
      If (icf1.eq.2 .and. (icf2.eq.0 .or. icf2.eq.1)) Then 
c         H(3) = 0d0             ! Delta  at O(as^0)
c         Do k=1,5,2
c            H(3) = H(3) + (XPDF1(i,k)+XPDF1(i,-k)+
c     +           4d0*(XPDF1(i,k+1)+XPDF1(i,-k-1)))/9d0
c         Enddo
c         H(1) = XPDF1(i,0)      ! Gluon  at O(as^1)
c         H(2) = 0d0             ! Sigma  at O(as^2)
c         Do k=1,6
c            H(2) = H(2)+XPDF1(i,k)+XPDF1(i,-k)
c         Enddo

c --- this is not final - just for an early DIS table
c         H(2) = 0d0             ! Delta  at O(as^0)
c         Do k=1,5,2
c            H(2) = H(2) + (XPDF1(i,k)+XPDF1(i,-k)+
c     +           4d0*(XPDF1(i,k+1)+XPDF1(i,-k-1)))/9d0
c         Enddo
c         H(1) = XPDF1(i,0)      ! Gluon  at O(as^1)
c         H(3) = 0d0             ! Sigma  at O(as^2)
c         Do k=1,6
c            H(3) = H(3)+XPDF1(i,k)+XPDF1(i,-k)
c         Enddo

c --- final: 1Delta   2Gluon  3Sigma.
         H(1) = 0d0             ! Delta  at O(as^0)
         Do k=1,5,2
            H(1) = H(1) + (XPDF1(i,k)+XPDF1(i,-k)+
     +           4d0*(XPDF1(i,k+1)+XPDF1(i,-k-1)))/9d0
         Enddo
         H(2) = XPDF1(i,0)      ! Gluon  at O(as^1)
         H(3) = 0d0             ! Sigma  at O(as^2)
         Do k=1,6
            H(3) = H(3)+XPDF1(i,k)+XPDF1(i,-k)
         Enddo


c --- hadron-hadron: jets
      Elseif (icf1.eq.3.and.icf2.eq.1.and.(icf3.ge.1.and.icf3.le.2))Then 
         SumQ1  = 0d0
         SumQB1 = 0d0
         SumQ2  = 0d0
         SumQB2 = 0d0
         Do k=1,6
            Q1(k)  = XPDF1(i,k)  ! read 1st PDF at x1
            QB1(k) = XPDF1(i,-k)
            SumQ1  = SumQ1  + Q1(k)
            SumQB1 = SumQB1 + QB1(k)
            Q2(k)  = XPDF2(j,k)  ! read 2nd PDF at x2
            QB2(k) = XPDF2(j,-k)
            SumQ2  = SumQ2  + Q2(k)
            SumQB2 = SumQB2 + QB2(k)
         Enddo
         G1     = XPDF1(i,0)
         G2     = XPDF1(j,0)
c   - compute S,A
         S = 0d0
         A = 0d0
         Do k=1,6
            S = S + (Q1(k)*Q2(k)) + (QB1(k)*QB2(k)) 
            A = A + (Q1(k)*QB2(k)) + (QB1(k)*Q2(k)) 
         Enddo
c   - compute seven combinations
         H(1) = G1*G2
         H(2) = SumQ1*SumQ2 + SumQB1*SumQB2 - S
         H(3) = S
         H(4) = A
         H(5) = SumQ1*SumQB2 + SumQB1*SumQ2 - A
         H(6) = (SumQ1+SumQB1)*G2
         H(7) = G1*(SumQ2+SumQB2)
         If (icf3.eq.1) H(6) = H(6)+H(7) ! case: 6 subproc

c --- gammaP: direct, jets
      Elseif (icf1.eq.2 .and. icf2.eq.2) Then 
         Write(*,*) '    gammaP to be implemented'
         Stop

      Else
         Write(*,*) '    icf1,2,3 =',icf1,icf2,icf3
         Write(*,*) '    this combination is not yet defined'
         Stop
      Endif

      Return
      End
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
*
* current restrictions:
* - not final - should be beautified - more info in output, like bin 
*   boundaries - non redundant in lower dimensions
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Include 'strings.inc'
      Integer i,j
      Double Precision min2old,max2old

      Write(*,*) "fastNLO results:     ndim=",Ndim
      Write(*,9003) IXsectUnits(1),Cunits(IXsectUnits(1))

      min2old=-99999d0
      max2old=-99999d0

      If (Ndim.eq.1) 
     +     Write(*,*) ' bins in ',dimlabel(Ndim)

      Do i=1,NObsBin
         If (Ndim.gt.1) Then
c            If ((min2old.ne.lobin(i,Ndim-1)).and.
c     +           (max2old.ne.upbin(i,Ndim-1))) Then
c               min2old = lobin(i,Ndim-1)
c               max2old = upbin(i,Ndim-1)
            If ((min2old.ne.lobin(i,2)).or.
     +           (max2old.ne.upbin(i,2))) Then
               min2old = lobin(i,2)
               max2old = upbin(i,2)
               Write(*,*) '------------------------------------'
c               Write(*,9005) min2old,max2old,dimlabel(Ndim-1)
               Write(*,9005) min2old,max2old,dimlabel(Ndim)
c               Write(*,*) ' bins in ',dimlabel(Ndim)
               Write(*,*) ' bins in ',dimlabel(1)
            Endif
         Endif
c         Write(*,9010) i,lobin(i,Ndim),upbin(i,Ndim),xsect(i)
         Write(*,9010) i,lobin(i,1),upbin(i,1),xsect(i)
      Enddo
 9003 Format (' results in units of 10^-',I2,' barn (',A2,')')
 9005 Format ('  range ',F10.3,' -',F10.3,' in ',A64)
 9010 Format (I4,F9.3,' -',F9.3,' :',E12.4)

      Return
      End

*******************************************************************
*******************************************************************
      Subroutine FX9999NF
*-----------------------------------------------------------------
* fastNLO user code v2.0 - print scenario information
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Integer i,j

      Write(*,*)
      WRITE(*,*)"########################################"//
     >     "################################"
      Write(*,*)"# Information on fastNLO scenario: ",ScenName
      WRITE(*,*)"# --------------------------------------"//
     >     "--------------------------------"
      Write(*,*)"# Description:"
      Do i=1,NScDescript
         Write(*,*)"#   ",SCDescript(i)  
      Enddo
      Write(*,*)"#"
      Write(*,*)"# Centre-of-mass energy Ecms: ",Ecms," GeV"
      Write(*,*)"#"
      Write(*,*)"# Tot. no. of observable bins: ",NObsBin," in ",Ndim,
     >     " dimensions:"
      Do j=1,NDim
         Write(*,*)"# Binning in dimension ",j,":",DimLabel(j)
         Do i=1,NObsBin
            If (IDiffBin(j).EQ.1) Then
               Write(*,*)"#   the bin center ",i," is at ",
     >              LoBin(i,j)
            Else
               Write(*,*)"#   the bin no. ",i," goes from ",
     >              LoBin(i,j),"up to ",UpBin(i,j)
            Endif
         Enddo
      Enddo
      Write(*,*)"#"
      Write(*,*)"# No. of contributions: ",Ncontrib
      Do i=1,Ncontrib         
         Write(*,*)"# Contribution",i,":"
         Do j=1,NcontrDescr(i)
            Write(*,*)"#   ",CtrbDescript(i,j)
         Enddo
         Write(*,*)"#   computed by:"
         Do j=1,NcodeDescr(i)
            Write(*,*)"#     ",CodeDescript(i,j)
         Enddo
      Enddo
      Write(*,*)"# No. of x bins in 1st contrib.:",Nxtot(1,1,1)
      Write(*,*)"# No. of available scale variations/scale nodes:"
      Do i=1,Ncontrib         
         Write(*,*)"#   NscaleVar,NScaleNode:",
     >        NscaleVar(i,1),NScaleNode(i,1)
      Enddo
      Write(*,*)"#"
      WRITE(*,*)"########################################"//
     >     "################################"

      Return
      End

*******************************************************************
*******************************************************************
      Subroutine FX9999TB(fnstring,ivar,dvar)
*-----------------------------------------------------------------
* fastNLO user code v2.0 - returns infortmation from table variables
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Character*(*) fnstring
      Integer ivar
      Double Precision dvar

      ivar = 0
      dvar = 0d0
      If (fnstring .eq. 'NContrib') Then
         ivar = NContrib
      Elseif (fnstring .eq. 'NObsBin') Then
         ivar = NObsBin
      Elseif (fnstring .eq. 'Nevt1') Then
         ivar = Nevt(1)
      Elseif (fnstring .eq. 'Nevt2') Then
         ivar = Nevt(2)
      Else
         Write(*,*) 'FX9999TB   unknown input - variable:',fnstring
         Stop
      Endif

      Return
      End

*******************************************************************
