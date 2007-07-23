*******************************************************************
*******************************************************************
* fastNLO user code            T. Kluge, M. Wobisch v1.4 02/01/2006      
*                                                   v1.5 31/05/2007  
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
* - include table-checks into read routine - or separate check routine?
*
*
*--------- v2.0 routines - status --------------------------------
*  fx9999cc     main routine - calls other code - need cosmetics        (works)
*  fx9999in     initialization (to do: check consistency)               (works)
*  fx9999pt     find pointers to contributions/scales                   (works)
*                 - interpret contribution selection/consistency (advanced)
*                 - maybe reorder contributions/optimize PDF access (not yet)
*                 - check availability of the contributions (not yet)
*  fx9999gp     get PDFs      - missing details                         (works)
*  fx9999pl     compute PDF linear combinations (done - cosmetics)      (works)
*  fx9999mt     multiply coefficients and PDFs                       (advanced)
*  fx9999pr     print results                                           (works)
*  fx9999rw     read (or write) table                                (complete)
*  fx9999nf     print scenario information (physics & technical)        (works)
*  fx9999nm     normalize distribution by its own integral              (to do)
*
*
* uses commonblock definitions in: fnXNNNN.inc
*
* needs the following routines from "fn-interface.f":
*     FNALPHAS        alpha_s interface (double precision function)
*     FNPDF           PDF interface
*
* also: package fnio.f (5 routines for i/o)
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
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Integer IFILE, ipoint, I,J,K,L,M, 
     +     IPrintFlag,
     +     maxscale, nbin,nx
      Character*(*) FILENAME
      Character*50 OLDFILENAME
      Double Precision Xmur, Xmuf
      Data OLDFILENAME/'xxxx'/
      Save OLDFILENAME


c === initialization: read table, set pointers to contributions
      call FX9999IN(Filename)
c === reset output array
      Do j=1,NObsBin
         xsect(j) = 0d0
      Enddo
c === determine pointers to contributions/scales
      Call FX9999PT(xmur,xmuf)

c === loop over pointers to contributions
      Do i=1,IContrib
         Ipoint = IContrPoint(i)
         write(*,*) "Ctrb",i,IContrPoint(i),IScalePoint(i),
     +        NSubproc(IContrPoint(i))
c - get PDFs
         call FX9999GP(i)
c - multiply with perturbative coefficients and alphas
         call FX9999MT(i,xmur,xmuf) ! <<< need to think about argument
c - add up in small array
         Do j=1,NObsBin
            Do k=1,NSubProc(Ipoint)
               xsect(j) = xsect(j)+result(j,k,Ipoint)
c               write(*,*) 'result ',j,k,result(j,k,Ipoint)
            Enddo
         Enddo
      Enddo

c === normalization
      If (INormFlag.eq.0) Then
         Continue               ! no normalization - nothing to do
      ElseIf (INormFlag.eq.1) Then
c         Call FX9999NM          ! normalize by own integral
         continue
c      ElseIf (INormFlag.eq.2 .or INormFlag.eq.3) Then ! get denominator, divide
c         Call DX9999CC(DenomTable,Xmur,Xmuf,0,Xsect2)
         Do i=1,NObsBin
c            Xsect(i) = Xsect(i)/Xsect2(IDivPointer(i))
         EndDo
      EndIf

c === print results - if requested
      If (IPrintFlag.eq.1) call FX9999PR(xsect)
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
      Integer i,j,k
      Character*(*) FILENAME
      Character*50 OLDFILENAME
      Data OLDFILENAME/'xxxx'/
      Save OLDFILENAME

c === reset result arrays
      Do i=1,MxObsBin
         Xsect(i) = 0d0
         Xsect2(i)= 0d0
         Do j=1,MxSubproc
            Do k=1,MxCtrb
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
      If (FILENAME.ne.OLDFILENAME) Then
c         Call FX9999RD(FILENAME)
         Call FX9999RW('read',FILENAME) ! new flexible version
         OLDFILENAME=FILENAME

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
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Integer i,j,k, i1
      Double Precision Xmur, Xmuf, XmurOld, XmufOld
      Data XMurOld/0d0/,XmufOld/0d0/

c === preliminary: assume fixed table structure: LO, NLO, threshcor
c   -->  still need to implement test if and where available
      IContrib = 0
      If (PORDPTHY.ge.1) then   ! LO contribution selected
         IContrib = IContrib+1
c         IContrib = 2           ! test: run NLO correction only
c            IContrPoint(IContrib) = function(IContrFlags: 1, 1, 0) <<<<<
         IContrPoint(IContrib) = 1 ! <<< assumes that LO comes first
      Endif
      If (PORDPTHY.ge.2) then   ! NLO contribution selected
         IContrib = IContrib+1
         IContrPoint(IContrib) = 2 ! <<< assumes that NLO comes second
      Endif
      If (PTHRESHCOR.eq.1) then ! 1-loop threshold corrections selected
         If (PORDPTHY.lt.1) then
            write(*,*)' inconsistent choice: 1-loop threshold corrections'
            write(*,*)'                      need to be matched to LO'
            Stop
         Endif
         If (PORDPTHY.ge.2) then
            write(*,*)' inconsistent choice: NLO and 1-loop threshold'
            write(*,*)'                     corrections are redundant'
            Stop
         Endif
         IContrib = IContrib+1
         IContrPoint(IContrib) = IContrib ! assume threshcor come after FO
      Endif
      If (PTHRESHCOR.eq.2) then ! 2-loop threshold corrections selected
         If (PORDPTHY.lt.2) then
            write(*,*)' inconsistent choice: 2-loop threshold corrections'
            write(*,*)'                      need to be matched to NLO'
            Stop
         Endif
         If (PORDPTHY.ge.3) then
            write(*,*)' inconsistent choice: NNLO and 2-loop threshold'
            write(*,*)'                     corrections are redundant'
            Stop
         Endif
         IContrib = IContrib+1
         IContrPoint(IContrib) = IContrib ! assume threshcor come after FO
      Endif

c ----- reorder, so that identical subprocesses are calculated in a row!!
c      Do i=1,Icontrib
c         write(*,*) 'Pointer:',IcontrPoint(i)
c      Enddo


c - check availability of fact scale choice and assign pointer
c    (based on fact scale since renorm scale is flexible
c         --- ecxept for threshold corrections ---
      Do i=1,Icontrib
         i1 = IContrPoint(i)
         IScalePoint(i) = 0
         If (IScaleDep(i1).eq.0) Then ! Born-type w/o scale dep - use any scale
            If (NScaleVar(i1,1).ge.1) Then
               IScalePoint(i) = 1
            Else
               Write(*,*) ' not a single scale variation in contrib.',i1
               Stop
            Endif
            If (NScaleVar(i1,1).gt.1) write(*,*) " why more than one",
     +           ' scale variation for Born-type contrib.',i1,'?'
         Else                   ! no Born type contribution
            Do j=1,NScaleVar(i1,1)
               If (xmuf.eq.scalefac(i1,1,j)) IScalePoint(i)=j
            Enddo
         Endif
         If (IScalePoint(i).eq.0) Then
            Write(*,*)' In fastNLO scenario ',ScenName
            Write(*,*)' the requested factorization scale xmuf of ',xmuf
            Write(*,*)' is not available for the contribution:'
            Do j=1,NContrDescr(i1)
               Write(*,*) ' ',CtrbDescript(i1,j)
            Enddo
            Stop
         Endif 
      Enddo

c - check if renormalization scale can be provided aposteriori if needed
      If (xmur .ne. xmuf) Then
         Do i=1,Icontrib
            i1 = IContrPoint(i)
            IF (IScaleDep(i1).eq.2) Then
               Write(*,*)' In fastNLO scenario ',ScenName
               Write(*,*)' the requested renormalization scale xmur of ',xmur
               Write(*,*)' is not available for the contribution:'
               Do j=1,NContrDescr(i1)
                  Write(*,*) ' ',CtrbDescript(i1,j)
               Enddo
               Write(*,*) '      (only xmur=xmuf is possible)'
               Stop
            Endif
         Enddo
      Endif

c ----- test output
      Do i=1,Icontrib
         write(*,*) 'Pointer:',i,IcontrPoint(i),IScalePoint(i)
      Enddo

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
*            IScalePoint only available for relative contrib.
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
      ic = IContrPoint(in)
      is = IScalePoint(in)

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
c               Do m=1,1 ! gg only
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
      ic = IContrPoint(in)
      is = IScalePoint(in)

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

c --- DIS: inclusive and jets
      if (icf1.eq.1 .and. icf2.eq.1 .and.(icf3.ge.1.and.icf3.le.3)) then 
         H(1) = 0d0             ! Delta  at O(as^0)
         do k=1,5,2
            H(1) = H(1) + (XPDF1(i,k)+XPDF1(i,-k)+
     +           4d0*(XPDF1(i,k+1)+XPDF1(i,-k-1)))/9d0
         enddo
         H(2) = XPDF1(i,0)      ! Gluon  at O(as^1)
         H(3) = 0d0             ! Sigma  at O(as^2)
         do k=1,6
            H(3) = H(3)+XPDF1(i,k)+XPDF1(i,-k)
         enddo

c --- hadron-hadron: jets
      elseif (icf1.eq.2.and.icf2.eq.1.and.(icf3.ge.1.and.icf3.le.2))then 
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

      Do i=1,NObsBin
         write(*,*) 'bin No/result ',i,'  ',xsect(i)
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
      Subroutine fx9999rw(crw,filename)
* ------------------------------------------------------------
*  fastNLO usercode v2.0 reads or writes a v2 table
*
*  MW 07/23/2007
*
*  character   crw:        'read' or 'write'
*  character   filename:    name of table
* ------------------------------------------------------------
      Implicit None
      Character*(*) crw,filename
      Integer Nunit, Ifile,ic,i,j,k,l,n,m, nxmax
      Include 'fnx9999.inc'

      Nunit=2

      If (crw.eq.'write') Then
         Open(unit=Nunit,file=filename,status='unknown')
      Else
         Open(Nunit,Status='OLD',File=Filename,IOSTAT=Ifile)
         If (Ifile .ne. 0) Then
            Write(*,*) '          fastNLO:  table file not found ',
     +           '  -  IOSTAT = ',Ifile
            Stop
         Endif
      Endif

c --- fastNLO table
c - block A1
      Call fnioisep(crw,nunit)
      Call fnioint(crw,nunit, Itabversion)
      Call fniochar(crw,nunit, ScenName)
      Call fnioint(crw,nunit, Ncontrib)
      Call fnioint(crw,nunit, Nmult)
      Call fnioint(crw,nunit, Ndata)
c - block A2
      Call fnioisep(crw,nunit)
      Call fnioint(crw,nunit, IpublUnits) 
      Call fnioint(crw,nunit, NscDescript)
      If (NscDescript .gt. MxScDescript) Then
         Write(*,*) ' NscDescript too large ',NscDescript,'<=',MXScDescript
      Endif
      Do i=1,NscDescript
        Call fniochar(crw,nunit, ScDescript(i))
      Enddo
      Call fniodbl(crw,nunit, Ecms)
      Call fnioint(crw,nunit, ILOord)
      Call fnioint(crw,nunit, NobsBin)
      Call fnioint(crw,nunit, NDim)
      Do i=1,NDim
         Call fniochar(crw,nunit, DimLabel(i))
      Enddo
      Do i=1,NDim
         Call fnioint(crw,nunit, IDiffBin(i))
      Enddo
      Do i=1,NObsBin
         Do j=1,NDim
            Call fniodbl(crw,nunit, LoBin(i,j))
            If (IDiffBin(j).eq.2) Call fniodbl(crw,nunit, UpBin(i,j))
         Enddo
      Enddo
      Call fnioint(crw,nunit, INormFlag)
      If (INormFlag.gt.1) Then
         Call fniochar(crw,nunit, DenomTable)
         Do i=1,NObsBin
            Call fnioint(crw,nunit, IDivPointer(i))
         Enddo
      Endif

c - block B
      Do ic=1,NContrib
         Call fnioisep(crw,nunit)
         Call fnioint(crw,nunit, IXsectUnits(ic))
         Call fnioint(crw,nunit, IDataFlag(ic))
         Call fnioint(crw,nunit, IAddMultFlag(ic))
         Call fnioint(crw,nunit, IContrFlag1(ic))
         Call fnioint(crw,nunit, IContrFlag2(ic))
         Call fnioint(crw,nunit, IContrFlag3(ic))
         Call fnioint(crw,nunit, NContrDescr(ic))
         Do i=1,NContrDescr(ic)
            Call fniochar(crw,nunit, CtrbDescript(ic,i))
         Enddo
         Call fnioint(crw,nunit, NCodeDescr(ic))
         Do i=1,NCodeDescr(ic)
            Call fniochar(crw,nunit, CodeDescript(ic,i))
         Enddo

c --------------------------- Idata ?????
         If (IDataFlag(ic).eq.1) Then
            Write(*,*) "   Data Blocks can not yet be read"
            STOP
            Goto 100
         Endif

c --------------------------- IAddMult ?????
         If (IAddMultFlag(ic).eq.1) Then
            Write(*,*) "   Multiplicative Blocks can not yet be read"
            STOP
            Goto 100
         Endif

c --- coefficient block
      Call fnioint(crw,nunit, IRef(ic))
      Call fnioint(crw,nunit, IScaleDep(ic))
      Call fniolint(crw,nunit, Nevt(ic))
      Call fnioint(crw,nunit, Npow(ic))
      Call fnioint(crw,nunit, NPDF(ic))
         Do i=1,NPDF(ic)
            Call fnioint(crw,nunit, NPDFPDG(ic,i))
         Enddo
         Call fnioint(crw,nunit, NPDFDim(ic))
         Call fnioint(crw,nunit, NFragFunc(ic))
         Do i=1,NFragFunc(ic)
            Call fnioint(crw,nunit, NFFPDG(ic,i))
         Enddo
         Call fnioint(crw,nunit, NFFDim(ic))
         Call fnioint(crw,nunit, NSubproc(ic))
         Call fnioint(crw,nunit, IPDFdef(ic,1))   
         Call fnioint(crw,nunit, IPDFdef(ic,2))  
         Call fnioint(crw,nunit, IPDFdef(ic,3))  

         IF (IPDFdef(ic,1).eq.0) then ! - no predefined set of PDF coefficients
            write(*,*) " case IPDFdef(1)=0 not yet implemented"
            STOP
         Endif
         If (NPDF(ic).gt.0) Then
            Do i=1,NObsBin
               Call fnioint(crw,nunit, Nxtot(ic,1,i))
               Do j=1,Nxtot(ic,1,i)
                  Call fniodbl(crw,nunit, XNode1(ic,i,j))
               Enddo
            Enddo 
            If (NPDFDim(ic).eq.2) Then
               Do i=1,NObsBin
                  Call fnioint(crw,nunit, Nxtot(ic,2,i))
                  Do j=1,Nxtot(ic,2,i)
                     Call fniodbl(crw,nunit, XNode2(ic,i,j))
                  Enddo
               Enddo 
            Endif
         Endif
         IF (NFragFunc(ic).gt.0) then ! - no FFs so far
            write(*,*) " fastNLO: no FragFuncs so far"
            STOP
         Endif

         Call fnioint(crw,nunit, NScales(ic))
         Call fnioint(crw,nunit, NScaleDim(ic))
         Do i=1,NScales(ic)
            Call fnioint(crw,nunit, IScale(ic,i))
         Enddo
         Do i=1,NScaleDim(ic)
            Call fnioint(crw,nunit, NScaleDescript(ic,i))
            Do j=1,NScaleDescript(ic,i)
               Call fniochar(crw,nunit, ScaleDescript(ic,i,j))
            Enddo
         Enddo
         Do i=1,NScaleDim(ic)
            Call fnioint(crw,nunit, NScaleVar(ic,i))
            Call fnioint(crw,nunit, NScaleNode(ic,i))
         Enddo
         Do i=1,NScaleDim(ic)
            Do j=1,NScaleVar(ic,i)
               Call fniodbl(crw,nunit, ScaleFac(ic,i,j))
            Enddo
         Enddo

         Do i=1,NObsBin
            Do j=1,NScaleDim(ic)
               Do k=1,NScaleVar(ic,j)
                  Do l=1,NScaleNode(ic,j)
                     Call fniodbl(crw,nunit, ScaleNode(ic,i,j,k,l))
                  Enddo
               Enddo
            Enddo
         Enddo

         Do i=1,NObsBin
            Do k=1,NScaleVar(ic,1)
               Do l=1,NScaleNode(ic,1)
c --- here we assume NFragFunc=0
                  If (NFragFunc(ic).gt.0) then
                     write(*,*) " NFragFunc>0 not yet implemented"
                     STOP
                  Endif
                  If (NPDFdim(ic).eq.0) Then
                     nxmax = Nxtot(ic,1,i)
                  Elseif (NPDFdim(ic).eq.1) Then
                     nxmax = (Nxtot(ic,1,i)**2+Nxtot(ic,1,i))/2
                  Elseif (NPDFdim(ic).eq.2) Then 
                     nxmax = 
     +                    Nxtot(ic,1,i)*Nxtot(ic,2,i)
                     write(*,*) '  NPDFdim = 2 not yet enabled'
                     Stop
                  Else
                     write(*,*) '  NPDFdim > 2 not enabled'
                     Stop
                  Endif
                  Do m=1,nxmax
                     Do n=1,NSubProc(ic)
                        Call fniodbl(crw,nunit, SigmaTilde(ic,i,1,k,l,m,n))
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo



 100     Continue
      Enddo

c      Call fnioint(crw,nunit, )
c      Call fniodbl(crw,nunit, )
c      Call fniochar(crw,nunit, )



c - end of table
      Call fnioisep(crw,nunit)
      Call fnioisep(crw,nunit)
      Close(2)

      Return
      End
