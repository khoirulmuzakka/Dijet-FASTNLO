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
*
*   xXnnnnIN        initialization 
*   xXnnnnNF        print scenario information (physics & technical)
*
* -------------------- from v1.4 ------- (to be updated)
* contains the following routines
*     fXnnnnCC        computes cross sections / fills arrays
*     fXnnnnRD        reads coefficient tables
*
*     fXnnnnGP        gets the PDF values for all x-bins
*     fXnnnnPL        computes PDF linear combinations
*     fXnnnnMT        multiply PDF array and coefficient array
*     fXnnnnPR        prints x-section results 
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
      Subroutine FX9999CC(FILENAME,IContFlag,XMUR,XMUF,IPRINTFLAG,XSECT)
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
      Integer IFILE, iord, isub, I,J,K,L,M, 
     +     IContFlag,IPrintFlag,
     +     ICtOld, 
     +     maxscale, nbin,nx, ixmur,ixmuf
      Character*(*) FILENAME
      Character*50 OLDFILENAME
      Double Precision Xmur, Xmuf, XmurOld, XmufOld
      Data OLDFILENAME/'xxxx'/
      Save OLDFILENAME
      Data ICtOld/-1/,XMurOld/0d0/,XmufOld/0d0/


c === initialization: read table, set pointers to contributions
      call FX9999IN(Filename,IContFlag,xmur,xmuf)

c === loop over pointers to contributions
      Do i=1,IContrib
         write(*,*) "Ctrb",i,IContrPoint(i),IScalePoint(i),
     +        NSubproc(IContrPoint(i))

c - get PDFs
         call FX9999GP(i)

c - multiply with perturbative coefficients and alphas
         call FX9999MT(xmur,IScalePoint(i)) ! <<< need to think about argument

c ----- maybe need no arguments at all? all info is known through commonblock

c   - add up in small array
         Do j=1,NObsBin
            Do k=1,NSubProc(IContrPoint(i))
               xsect(j) = result(j,k,i)
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
      SUBROUTINE FX9999IN(Filename,IContFlag,xmur,xmuf)
*-----------------------------------------------------------------
* MW 06/10/2007
*
* initialize fastNLO code
* 
* input:    IContFlag   defines which contributions shall be added
*           Xmur        renormalization scale factor        
*           Xmuf        factorization scale factor
*
* 06/11/2007 MW
*-----------------------------------------------------------------
      Implicit None
      Include 'fnx9999.inc'
      Integer IContFlag, ICtOld, i,j,k
c      Double Precision Xmur, Xmuf
      Character*(*) FILENAME
      Character*50 OLDFILENAME
      Double Precision Xmur, Xmuf, XmurOld, XmufOld
      Data OLDFILENAME/'xxxx'/
      Save OLDFILENAME
      Data ICtOld/-1/,XMurOld/0d0/,XmufOld/0d0/

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
         Call FX9999RD(FILENAME)
         OLDFILENAME=FILENAME
         ICtOld  = -1
         XMurOld = 0d0
         XmufOld = 0d0

c   - check consistency of array dimensions / commonblock parameters
c ----------> to be done
      Endif

c === if new combination flag: check availability - set pointers to contrib.
c === if new scale choice: check availability - set pointers to scale var.
      If (IcontFlag.ne.ICtOld .or. Xmur.ne.XmurOld .or. Xmuf.ne.XmufOld) then
         ICtOld = IContFlag
         XmurOld = Xmur
         XmufOld = Xmuf

c === preliminary: assume fixed structure LO, NLO, threshcor
c   -->  still need to implement test if and where available
         IContrib = 0
         If (PORDPTHY.ge.1) then ! LO contribution selected
            IContrib = IContrib+1
            IContrPoint(IContrib) = 1 ! <<< assumes that LO comes first
         Endif
         If (PORDPTHY.ge.2) then ! NLO contribution selected
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


c - check availability of scale choice and assign pointer
c     (based on fact scale since renorm scale flexible
c       --- ecxept for threshold corrections ---

c  --->  preliminary choice:use 3rd scale (usually =pT for pp)
         IScalePoint(1) = 1
         If (IContrib.ge.2) IScalePoint(2) = 3
         If (IContrib.ge.3) IScalePoint(3) = 3

      Endif

 5000 Format (A,A64)
      Return 
      End

*******************************************************************
*******************************************************************
      SUBROUTINE FX9999MT(xmur,ixmuf)
*-----------------------------------------------------------------
* MW 06/10/2007
*
* multiply the PDFs and the perturbative coefficients 
* 
* input:    XMUR  prefactor for nominal renormalization scale
*           IXMUF No.of factorization scale setting (as stored in table) 
*
*-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'fnx9999.inc'
      INTEGER IXMUF, i,j,k,l,m,iord, jord,    nbin,nx
      Double Precision xmur

      Return 
      End

*******************************************************************
*******************************************************************
      SUBROUTINE FX9999GP(IC)
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
     +     tmppdf(-6:6), xpdf(MxNxTot,-6:6), H(7)
      Integer ic, i,j,k,l,m, nx, nx2limit

c - temp variables for cross check
      Double Precision reweight


      Do i=1,NObsBin
         Do j=1,NScaleNode(ic,1)
            muf = ScaleNode(ic,i,1,IScalePoint(ic),j)
            nx =0
            Do k=1,NxTot(ic,1)
               x = Xnode1(ic,i,k) 
               Call FNPDF(x, muf, tmppdf)
c - temporary - just to cross check old results
               reweight = 1d0
               reweight = sqrt(x)/(1d0-0.99d0*x)**3
               do l=-6,6
                  XPDF(j,l) = tmppdf(l) * reweight
               enddo

               nx2limit = k
c - find exact condition later
c               if ( DIS: nx2limit=1
c -------------- very different treatment for full matrix w/ diff dimensions
               Do l=1,nx2limit
                  nx = nx+1
c                  call fx9999pl(ireaction,k,l,XPDF,H)
c --- Nsubproc is old -> problem: don't want to use real Nsubproc
c           since for compatibility mode we may need more than actual Nsubproc 
c                  Do m=1,Nsubproc
c                     pdf(nbin,nx,m,p) = H(m)
c                  Enddo
               Enddo            ! x2-loop
            Enddo               ! x1-loop
         Enddo                  ! ScaleNode-loop
      Enddo                     ! ObsBin-loop

      Return
      End

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
      IMPLICIT NONE
      INTEGER i,j
      INCLUDE 'fnx9999.inc'

      write(*,*) "fastNLO results:"

      Return
      End

*******************************************************************
*******************************************************************
      SUBROUTINE FX9999RD(FILENAME,IFLG)
*-----------------------------------------------------------------
* M. Wobisch   read ASCII table of perturbative coefficients
*
* input: FILENAME  name of table
*        IFLG =0 read header (A) 
*             =1 read header plus LO contribution
*             =2 read header plus LO plus NLO contribution
*             =99 read whole table
*
* To Do:
* - problem: Nevt too long?
* - read data blocks
* - read multiplicative factor-blocks
*
* MW 31/05/2007 completely rewritten for v1.5
*-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) FILENAME
      CHARACTER*255 BUFFER
      INTEGER IFLG, IFIRST, IFILE, IC,I,J,K,L,M,N,   nxmax
      INCLUDE 'fnx9999.inc'

      DATA IFIRST/0/
      SAVE IFIRST

      OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=IFILE)
      IF (IFILE .ne. 0) THEN
         WRITE(*,*) '          fastNLO:  table file not found ',
     +        '  -  IOSTAT = ',IFILE
         STOP
      ENDIF

c ---------------------------- block A1
      Read(2,*) i
      If (i.ne.Iseparator) Goto 999
      Read(2,*) Itabversion
      Read(2,*) ScenName
      Read(2,*) Ncontrib
      Read(2,*) Nmult
      Read(2,*) Ndata
c ---------------------------- block A2
      Read(2,*) i
      If (i.ne.Iseparator) Goto 999
      Read(2,*) IpublUnits
      Read(2,*) NscDescript
      If (NscDescript .gt. MxScDescript) then
         write(*,*) ' NscDescript too large ',NscDescript,'<=',MXScDescript
      Endif
      Do i=1,NscDescript
         Read(2,*) ScDescript(i)
      Enddo
      Read(2,*) Ecms
      Read(2,*) ILOord
      Read(2,*) NobsBin
      Read(2,*) NDim
      Do i=1,NDim
         Read(2,*) DimLabel(i)
      Enddo
      Do i=1,NDim
         Read(2,*) IDiffBin(i) 
      Enddo
      Do i=1,NObsBin
         Do j=1,NDim
            Read(2,*) LoBin(i,j)
            If (IDiffBin(j).eq.2) Read(2,*) Upbin(i,j)
         Enddo
      Enddo
      Read(2,*) INormFlag
      If (INormFlag.gt.1) then
         Read(2,*) DenomTable
         Do i=1,NObsBin
            Read(2,*) IDivPointer(i)
         Enddo
      Endif

c ---------------------------- block B  <<<< need to read this multiple times
      Do ic=1,NContrib
         Read(2,*) i
         If (i.ne.Iseparator) Goto 999
         Read(2,*) IXsectUnits(ic)
         Read(2,*) IDataFlag(ic)
         Read(2,*) IAddMultFlag(ic)
         Read(2,*) IContrFlag1(ic)
         Read(2,*) IContrFlag2(ic)
         Read(2,*) IContrFlag3(ic)
         Read(2,*) NContrDescr(ic)
         Do i=1,NContrDescr(ic)
            Read(2,*) CtrbDescript(ic,i)
         Enddo
         Read(2,*) NCodeDescr(ic)
         Do i=1,NCodeDescr(ic)
            Read(2,*) CodeDescript(ic,i)
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
         Read(2,*) IRef(ic)
         Read(2,*) IScaleDep(ic)
         Read(2,*) Nevt(ic)
         Read(2,*) Npow(ic)
         Read(2,*) NPDF(ic)
         Do i=1,NPDF(ic)
            Read(2,*) NPDFPDG(ic,i)  
         Enddo
         Read(2,*) NPDFDim(ic)
         Read(2,*) NFragFunc(ic)
         Do i=1,NFragFunc(ic)
            Read(2,*) NFFPDG(ic,i)
         Enddo
         Read(2,*) NFFDim(ic) 
         Read(2,*) NSubproc(ic)
         Read(2,*) IPDFCoeff(ic)   

         IF (IPDFCoeff(ic).eq.0) then ! - no predefined set of PDF coefficients
            write(*,*) " case IPDFCoeff=0 not yet implemented"
            STOP
         Endif
         If (NPDF(ic).gt.0) Then
            Read(2,*) Nxtot(ic,1)
            Do i=1,NObsBin
               Do j=1,Nxtot(ic,1)
                  Read(2,*) XNode1(ic,i,j)
               Enddo
            Enddo 
            If (NPDFDim(ic).eq.2) Then
               Read(2,*) Nxtot(ic,2)
               Do i=1,NObsBin
                  Do j=1,Nxtot(ic,2)
                     Read(2,*) XNode2(ic,i,j)
                  Enddo
               Enddo 
            Endif
         Endif
         IF (NFragFunc(ic).gt.0) then ! - no FFs so far
            write(*,*) " no FragFuncs so far"
            STOP
         Endif
         Read(2,*) NScales(ic)
         Read(2,*) NScaleDim(ic)
         Do i=1,NScales(ic)
            Read(2,*) IScale(ic,i)
         Enddo
         Do i=1,NScaleDim(ic)
            Read(2,*) NScaleDescript(ic,i)
            Do j=1,NScaleDescript(ic,i)
               Read(2,*) ScaleDescript(ic,i,j)
            Enddo
         Enddo
         Do i=1,NScaleDim(ic)
            Read(2,*) NScaleVar(ic,i)
            Read(2,*) NScaleNode(ic,i)
         Enddo
         Do i=1,NScaleDim(ic)
            Do j=1,NScaleVar(ic,i)
               Read(2,*) ScaleFac(ic,i,j)
            Enddo
         Enddo

         Do i=1,NObsBin
            Do j=1,NScaleDim(ic)
               Do k=1,NScaleVar(ic,j)
                  Do l=1,NScaleNode(ic,j)
                     Read(2,*)ScaleNode(ic,i,j,k,l)
                  Enddo
               Enddo
            Enddo
         Enddo

         Do i=1,NObsBin
            Do j=1,NScaleDim(ic)
               Do k=1,NScaleVar(ic,j)
                  Do l=1,NScaleNode(ic,j)
c --- here we assume NFragFunc=0
                     If (NFragFunc(ic).gt.0) then
                        write(*,*) " NFragFunc>0 not yet implemented"
                        STOP
                     Endif
                     IF (NPDFdim(ic).eq.0) nxmax = Nxtot(ic,1)
                     IF (NPDFdim(ic).eq.1) nxmax = 
     +                    (Nxtot(ic,1)**2+Nxtot(ic,1))/2
                     IF (NPDFdim(ic).eq.2) nxmax = 
     +                    Nxtot(ic,1)*Nxtot(ic,2)
                     Do m=1,nxmax
                        Do n=1,NSubProc(ic)
                           Read(2,*) SigmaTilde(ic,i,j,k,l,m,n)
                        Enddo
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo

 100     Continue
      Enddo

      Read(2,*) i
      If (i.ne.Iseparator) Goto 999
      Read(2,*) i
      If (i.ne.Iseparator) Goto 999
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
