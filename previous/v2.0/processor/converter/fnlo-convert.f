      SUBROUTINE FNLOCONVERT
************************************************************************
*
*     fastNLO_converter:
*     Convert fastNLO tables from v1.4 to v2.0
*
*     M. Wobisch, K. Rabbertz
*
*     Contains:
*     CV14RD   reads  v1.4 tables
*     CV20WRT  writes v2.0 tables
*
************************************************************************
      Implicit None
      Include 'strings.inc'
      Character*1 CH1TMP
      Character*2 CH2TMP
      Character*8 CH8TMP
      Character*255 INFILE,OUTFILE,SCENNAME
      Double Precision DBWID
      Integer ICONT, ISECT

*---  Parse command line
      WRITE(*,'(A)')
      WRITE(*,*)CSEPS
      WRITE(*,*)"# fnlo-convert: Program Steering"
      WRITE(*,*)LSEPS
*---  Input table name
      IF (IARGC().LT.1) THEN
         INFILE  = "intable.tab"
         OUTFILE = "outtable.tab"
         WRITE(*,*)"# fnlo-convert: WARNING! No table names given,"
         WRITE(*,*)"#   taking the defaults intable.tab and"
         WRITE(*,*)"#   outtable.tab instead!"
         WRITE(*,*)"#   For an explanation of command line "//
     >        "arguments type:"
         WRITE(*,*)"#   ./fnlo-convert -h"
      ELSE
         CALL GETARG(1,INFILE)
         IF (INFILE(1:LEN_TRIM(INFILE)).EQ."-h") THEN
            WRITE(*,*)'#'
            WRITE(*,*)'# Usage: ./fnlo-convert [arguments]'
            WRITE(*,*)'# Table input file,  def. = intable.tab'
            WRITE(*,*)'# Table output file, def. = outtable.tab'
            WRITE(*,*)'# Scenarioname, def. = fnx9999'
            WRITE(*,*)'# Table part to write, def. = 9 (all)'
            WRITE(*,*)'#    0: Table header only'
            WRITE(*,*)'#    1: LO contribution only'
            WRITE(*,*)'#    2: NLO contribution only'
            WRITE(*,*)'#    3: THC contribution only'
            WRITE(*,*)'#    9: Complete table'
            WRITE(*,*)'# Bin widths were not stored in v1.4 tables!'
            WRITE(*,*)'# ==> calculate bin widths from 2-dim. binning!'
            WRITE(*,*)'# Set additional bin width factor, def. = 1'
            WRITE(*,*)'#   0.: Ignore all binning, bin widths = 1.'
            WRITE(*,*)'#  -1.: Ignore bin widths of 1st dim. (y)'
            WRITE(*,*)'#  -2.: Ignore bin widths of 2nd dim. (pT)'
            WRITE(*,*)'#   2.: Multiply by 2,e.g. for |eta| or |y| bins'
            WRITE(*,*)'# Cross section units, def. = read from table'
            WRITE(*,*)'#'
            WRITE(*,*)'# Use "_" to skip changing a default argument.'
            WRITE(*,*)'#'
            STOP
         ELSEIF (INFILE(1:1).EQ."_") THEN
            INFILE = "intable.tab"
            WRITE(*,*)
            WRITE(*,*)"# fnlo-convert: WARNING! No intable name given,"
            WRITE(*,*)"#   taking the default outtable.tab instead!"
         ELSE
            WRITE(*,*)"# fnlo-convert: Converting table: ",
     >           INFILE(1:LEN_TRIM(INFILE))
         ENDIF
      ENDIF
*---  Output table name
      OUTFILE = "outtable.tab"
      IF (IARGC().GE.2) THEN
         CALL GETARG(2,OUTFILE)
      ENDIF
      IF (IARGC().LT.2.OR.OUTFILE(1:1).EQ."_") THEN
         OUTFILE = "outtable.tab"
         WRITE(*,*)"# fnlo-convert: WARNING! No outtable name given,"
         WRITE(*,*)"#   taking the default outtable.tab instead!"
      ENDIF
      WRITE(*,*)"# fnlo-convert: Writing to table: ",
     >     OUTFILE(1:LEN_TRIM(OUTFILE))
*---  Scenario name
      SCENNAME = "fnx9999"
      IF (IARGC().GE.3) THEN
         CALL GETARG(3,SCENNAME)
      ENDIF
      IF (IARGC().LT.3.OR.SCENNAME(1:1).EQ."_") THEN
         SCENNAME = "fnx9999"
         WRITE(*,*)"# fnlo-convert: WARNING! No scenario name given,"
         WRITE(*,*)"#   taking the default fnx9999 instead!"
      ENDIF
      WRITE(*,*)"# fnlo-convert: Setting scenario name to: ",
     >     SCENNAME(1:LEN_TRIM(SCENNAME))
*---  Table part
      CH1TMP = "X"
      IF (IARGC().GE.4) THEN
         CALL GETARG(4,CH1TMP)
      ENDIF
      IF (IARGC().LT.4.OR.CH1TMP(1:1).EQ."_") THEN
         CH1TMP = "9"
         ICONT  =  9
         WRITE(*,*)"# fnlo-convert: No table part selected,"//
     >        " writing complete table."
      ELSE
         READ(CH1TMP,'(I1)'),ICONT
         WRITE(*,*)"# fnlo-convert: Selected table part: ",ICONT
      ENDIF
*---  Bin width factor
      CH8TMP = "X"
      IF (IARGC().GE.5) THEN
         CALL GETARG(5,CH8TMP)
      ENDIF
      IF (IARGC().LT.5.OR.CH8TMP(1:1).EQ."_") THEN
         CH8TMP = "1"
         DBWID  =  1.D0
         WRITE(*,*)"# fnlo-convert: No bin width factor given,"
         WRITE(*,*)"#   using default of 1."
      ELSE
         READ(CH8TMP,'(F9.6)'),DBWID
         WRITE(*,*)"# fnlo-convert: Using bin width factor of: ",DBWID
      ENDIF
*---  Cross section units
      CH2TMP = "X"
      IF (IARGC().GE.6) THEN
         CALL GETARG(6,CH2TMP)
      ENDIF
      IF (IARGC().LT.6.OR.CH2TMP(1:1).EQ."_") THEN
         CH2TMP = "-1"
         ISECT  =  -1
         WRITE(*,*)"# fnlo-convert: No cross section unit given,"//
     >        " taking default from table."
      ELSE
         READ(CH2TMP,'(I2)'),ISECT
         WRITE(*,*)"# fnlo-convert: Setting cross section unit to: ",
     >        ISECT
      ENDIF
*---  Too many arguments
      IF (IARGC().GT.6) THEN
         WRITE(*,*)
     >        "fnlo-convert: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF
      WRITE(*,*)CSEPS

*---  Read the fastNLO v1.4 coefficient table
      WRITE(*,'(A)')""
      WRITE(*,*)CSEPS
      WRITE(*,*)"# fnlo-convert: Read v1.4 table"
      WRITE(*,*)CSEPS
      Call CV14RD(INFILE)
      WRITE(*,*)CSEPS

*---  Write the selected fastNLO v2.0 parts
      WRITE(*,'(A)')""
      WRITE(*,*)CSEPS
      WRITE(*,*)"# fnlo-convert: Write selected v2.0 parts"
      WRITE(*,*)CSEPS
      Call CV20WRT(OUTFILE,SCENNAME,ICONT,DBWID,ISECT)
      WRITE(*,*)CSEPS

*--- Finished
      WRITE(*,'(A)')""
      WRITE(*,*)CSEPS
      Write(*,*) "# fnlo-convert: Finished!"
      WRITE(*,*)CSEPS
      WRITE(*,'(A)')""

      Return
      End



      SUBROUTINE CV14RD(FILENAME)
************************************************************************
*
*     M. Wobisch, K. Rabbertz
*     Read ASCII table of perturbative coefficients from v1.4
*
*     Input: FILENAME  name of table
*
************************************************************************
      IMPLICIT NONE
      CHARACTER*(*) FILENAME
      CHARACTER*64 CHTMP
      INTEGER IFIRST, IFILE, I,J,K,L,M,N, NBIN,NX
      INCLUDE 'fnlo-convert.inc'

      DATA IFIRST/0/
      SAVE IFIRST

      OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=IFILE)
      IF (IFILE .ne. 0) THEN
         Write(*,*)'# CV14RD: ERROR! Table file not found, '//
     >        'IOSTAT = ',IFILE
         STOP
      ENDIF

*---  Table start
      Read(2,*) i
      if (i.ne.iseparator) goto 999
      Read(2,*) ITABVERSION
      Write(*,'(X,A,F4.1)') "# CV14RD: Table format is version",
     >     dble(itabversion)/10000d0
      if (ITABVERSION.ne.14000) then
         Write(*,*)
     >        '# CV14RD: ERROR! This conversion works only for v1.4!'
         stop
      endif
      Read(2,*) i
      if (i.ne.iseparator) goto 999
*---------------------------------------
      Read(2,*) IREACTION
      Read(2,*) ECMS
      Read(2,*) IXSECTUNITS
      Do i=1,5
         Read(2,'(A)') NAMELABEL(i)
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
         Read(2,'(A)') POWLABEL(i)
         Read(2,'(A)') CODELABEL(i)
      Enddo
      Read(2,*) i
      if (i.ne.iseparator) goto 999
*---------------------------------------
      Do i=1,nord
         Read(2,*) NEVT(i)
      Enddo
      Do i=1,nord
         CHTMP = POWLABEL(i)
         Write(*,'(A,I12,A,A)') ' # CV14RD: ',NEVT(i),
     >        ' events in ',CHTMP(1:LEN_TRIM(CHTMP))
      Enddo
      Read(2,*) NXTOT
      Write(*,*) "# CV14RD: No. of x bins: ",NXTOT
      Read(2,*) IXSCHEME
      Read(2,*) IPDFWGT
      Read(2,*) IREF
      Read(2,*) i
      if (i.ne.iseparator) goto 999
*---------------------------------------
      Read(2,*) NBINTOT
      Read(2,*) NDIMENSION
      Do i=1,ndimension
         Read(2,'(A)') DIMLABEL(i)
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
*---------------------------------------
      Do i=1,nrapidity
         Do j=1,NPT(i)
            Read(2,*) XLIMIT(i,j)
         Enddo
      Enddo
      Read(2,*) i
      if (i.ne.iseparator) goto 999
*---------------------------------------
      Read(2,'(A)') SCALELABEL
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
*---------------------------------------
      Do i=1,nrapidity
         Do j=1,NPT(i)
            Do k=1,NSCALEBIN
               Read(2,*) MUFVAL(i,j,k)
            Enddo
         Enddo
      Enddo
      Read(2,*) i
      if (i.ne.iseparator) goto 999
*---------------------------------------
      Read(2,*) NSCALEVAR
      Do i=1,NSCALEVAR
         Read(2,*) MURSCALE(i)
      Enddo
      Do i=1,NSCALEVAR
         Read(2,*) MUFSCALE(i)
      Enddo

      Read(2,*) i
      if (i.ne.iseparator) goto 999
*---------------------------------------
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
            nbin=nbin+1         ! Linear numbering of (rap,pT) bins
            Do k=1,NXSUM        ! Tot. no. of x-Bins
               Do m=1,Nsubproc  ! No. of sub processes
                  Do n=1,1+NSCALEVAR*(NORD-1) ! LO & NLO & w/ scale var
                     Do l=1,NSCALEBIN ! No. of bins in Scale
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
      Write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      Write(*,*) "> CV14RD: ERROR in table format!"
      Write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      Stop
      Return

      End



      SUBROUTINE CV20WRT(FILENAME,SCENNAME,ICONT,DBWID,ISECT)
************************************************************************
*
*     M. Wobisch, K. Rabbertz
*     Write ASCII table in v2.0 format
*
*     Input: FILENAME  name of table
*     Input: SCENNAME  name of scenario
*     Input: ICONT     part of table to write
*                   0: Table header only'
*                   1: LO contribution only'
*                   2: NLO contribution only'
*                   3: THC contribution only'
*                   9: Complete table'
*     Input: DBWID     factor for bin width calculation
*     Input: ISECT     Correct cross section units to ISECT
*
************************************************************************
      Implicit None
      Integer ICONT
      CHARACTER*(*) filename,scenname
      Double Precision DBWID
      Integer ISECT, i,j,k,l,m,n,n1,n2,o,p, nbin, i1,i2
      Include 'fnlo-convert.inc'
      Include 'strings.inc'

      Integer nctrmin,nctrmax, mmax
*---  Assume max 500 ObsBins & max 50 xNodes
      Double Precision xarray(500,50),rewgt(500,1300), w1,w2

*---  New variables that did not exist in v14
      Character*4 CH4TMP
      Character*64 CodeDescr(5),CHTMP
      Integer Ncontrib, Ndim, IDiffBin(2), INormFlag, IDataFlag,
     >     IAddMultFlag, IContrFlag1, IContrFlag2, NScaleDep,
     >     NcontrDescr, NCodeDescr, IScaleDep, NPDF, NPDFPDG(2),
     >     NPDFDim, NFragFunc, NSubpr, IPDFCoeff, IPDFdef1, IPDFdef2,
     >     IPDFdef3, NScales, NScaleDim, nxall, nscaddr, nscvar
      Double Precision binsize, XNode, hxlim, hx, a
*---  Variables for diagonalization of Bernstein Scale Interpolation
      Double Precision Bnst(4)
      Double Precision t1,t2
      Double Precision Bnfactor

*---  Pointer to subprocesses - new ordering
      Integer ipppoint(7),idispoint(3)
      Data ipppoint/1,4,5,6,7,2,3/, idispoint/3,1,2/

*---  Ncontrib
      Ncontrib = Nord
*---  Check old 'Ndimension'
      If (Ndimension.ne.2) then
         Write (*,*) '# CV20WRT: WARNING! Ndimension different from 2',
     >        Ndimension
      endif
      Ndim = Ndimension

*---  Differential or binned?  -> so far always binned
      IDiffBin(1) = 2
      IDiffBin(2) = 2
*---  Normalization in last dimension? -> so far never normalized
      INormFlag = 0
      IDataFlag = 0
      IAddMultFlag = 0
      IContrFlag1 = 1
      IContrFlag2 = 0
      NScaleDep = 0
      NContrDescr = 1
      NCodeDescr = 3

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
      If (ireaction.eq.2) NPDFPDG(2) = NPDFPDG(1) ! pp

      NFragFunc = 0
      NScales = 2
      NScaleDim = 1

      Write(*,*) '# CV20WRT: Write table: ',
     >     filename(1:len_trim(filename))

      open(unit=2,file=filename,status='unknown')
*---  Block A1
      Write(2,'(I10)') Iseparator
      Write(2,'(I5)')  20000    ! Itabversion
      Write(2,'(A)')   scenname(1:len_trim(scenname)) ! ScenName
      Write(2,'(I2)')  Ncontrib ! Ncontrib
      Write(2,'(I2)')  0        ! NMult
      Write(2,'(I2)')  0        ! NData
      Write(2,'(I2)')  0        ! Nuserstrings
      Write(2,'(I2)')  0        ! Nuserint
      Write(2,'(I2)')  0        ! Nuserfloat
      Write(2,'(I2)')  0        ! Imachine
*---  Block A2
      Write(2,'(I10)') Iseparator
      If (ISECT.ne.-1) Then
         IXSECTUNITS = ISECT
         Write(*,*) '# CV20WRT: Cross section units are corrected to: ',
     >        ISECT
      Endif
      Write(2,'(I2)')  IXSECTUNITS ! -> Ipublunits
*---  Compose five line scenario description out of available
*---  information and order as in v2.0
      Write(2,'(I2)')  5        ! NScDescript
*---  Observable
      CHTMP = NAMELABEL(1)
      Write(2,'(A)') CHTMP(1:LEN_TRIM(CHTMP)) ! ScDescript(1)
*---  Collaboration
      CHTMP = NAMELABEL(3)
      Write(2,'(A)') CHTMP(1:LEN_TRIM(CHTMP)) ! ScDescript(2)
*---  Process: Compose one additional line from process string
      CHTMP = CIPROC(IPROC)
      Write(2,'(A)') CHTMP(1:LEN_TRIM(CHTMP)) ! ScDescript(3)
*---  Jet algo: Compose another additional line from jet algo definition
      CHTMP = CIALGO(IALGO)
      CHTMP = CHTMP(1:LEN_TRIM(CHTMP))//"_"//CJETRES1(IALGO)
      CHTMP = CHTMP(1:LEN_TRIM(CHTMP))//"="
      Write(CH4TMP,'(F4.2)')JETRES1
      CHTMP = CHTMP(1:LEN_TRIM(CHTMP))//CH4TMP(1:LEN_TRIM(CH4TMP))
      Write(2,'(A)') CHTMP(1:LEN_TRIM(CHTMP)) ! ScDescript(4)
*---  Publication
      CHTMP = NAMELABEL(2)
      Write(2,'(A)') CHTMP(1:LEN_TRIM(CHTMP)) ! ScDescript(5)
*---
      Write(2,'(G24.17)') ECMS  ! Ecms
      Write(2,'(I2)') NPOW(1)   ! ILOord
      Write(2,'(I5)') NBINTOT   ! NObsBin
      Write(2,'(I2)') Ndim      ! NDim    (so far = 2 - change where needed)
      Do i=Ndim,1,-1
         CHTMP = DIMLABEL(i)
         Write(2,'(A)') CHTMP(1:LEN_TRIM(CHTMP))
      Enddo
      Do i=Ndim,1,-1
         Write(2,'(I2)') IDiffBin(i) ! IDiffBin
      Enddo

      i = 0
      Do i1=1,Nrapidity
         Do i2=1,NPT(i1)
            i = i+1             ! --- bin counter - linear loop over all bins
            Do j=Ndim,1,-1
               If (j.eq.1) then
                  Write(2,'(G24.17)') RAPBIN(i1)
                  Write(2,'(G24.17)') RAPBIN(i1+1)
               elseif (j.eq.2) then
                  Write(2,'(G24.17)') PTBIN(i1,i2)
                  Write(2,'(G24.17)') PTBIN(i1,i2+1)
               endif
            Enddo
         Enddo
      Enddo

      Write(*,*)'# CV20WRT: Bin width is computed using '//
     >     'a factor of:', DBWID
      Do i1=1,Nrapidity
         Do i2=1,NPT(i1)
            If (DBWID.EQ.0d0) then
               binsize = 1d0
            elseif (DBWID.EQ.-1d0) then
               binsize = (PTBIN(i1,i2+1)-PTBIN(i1,i2))
            elseif (DBWID.EQ.-2d0) then
               binsize = (RAPBIN(i1+1)-RAPBIN(i1))
            else
               binsize = (PTBIN(i1,i2+1)-PTBIN(i1,i2)) *
     >              (RAPBIN(i1+1)-RAPBIN(i1)) * DBWID
            endif
            Write(2,'(G24.17)')binsize
         Enddo
      Enddo
      Write(2,'(I2)') INormFlag ! INormFlag

*---  Block B
      Write(2,'(I10)') Iseparator

      If (ICONT.eq.0) then
         goto 999
      Elseif (ICONT.eq.9) then
         nctrmin = 1
         nctrmax = Ncontrib
      Else
         nctrmin = ICONT
         nctrmax = ICONT
         If (Ncontrib.lt.ICONT) then
            Write(*,*) '# CV20WRT: ERROR! Contribution does not exist:',
     >           ICONT
            stop
         Endif
      Endif

      Do n=nctrmin,nctrmax      ! loop over all contributions  B01...
         Write(2,'(I2)') IXSECTUNITS ! Ixsectunits (for each block)
         Write(2,'(I2)') IDataFlag ! IDataFlag
         Write(2,'(I2)') IAddMultFlag ! IAddMultFlag
         If (n.lt.1) then
            Write(*,*) '# CV20WRT: ERROR! No contribution found:',N
            Stop
         Elseif (n.eq.1) then   ! LO
            IContrFlag1 = 1
            IContrFlag2 = 1
            NScaleDep   = 0
         Elseif (n.eq.2) then   ! NLO
            IContrFlag1 = 1
            IContrFlag2 = 2
            NScaleDep   = 0
         Elseif (n.eq.3) then   ! THC 2-loop (Attention: old v2 def 2,1,2 --> new v21 2,2,0)
            IContrFlag1 = 2
            IContrFlag2 = 2
            NScaleDep   = 0
         Endif
         Write(2,'(I2)') IcontrFlag1 ! IcontrFlag1
         Write(2,'(I2)') IcontrFlag2 ! IcontrFlag2
         Write(2,'(I2)') NScaleDep ! NScaleDep
         Write(2,'(I2)') NContrDescr ! NcontrDescr = 1
         Do i=1,NContrDescr
            CHTMP = POWLABEL(n)
            Write(2,'(A)') CHTMP(1:LEN_TRIM(CHTMP)) ! CtrbDescript(i)
         Enddo
         Write(2,'(I2)') NCodeDescr ! NCodeDescr
*---  Replace with actual v2.0 type of code description
*---  N = 1, 2: NLOJet++; N = 3: Threshold Corrections
         IF (N.LE.2) THEN
            Do i=1,3
               CodeDescr(i) = CNLOJET(i)
            Enddo
         ELSEIF (N.EQ.3) THEN
            Do i=1,3
               CodeDescr(i) = CTHRCOR(i)
            Enddo
         ELSE
            CodeDescr(1) = CODELABEL(n)
         ENDIF
         Do i=1,NCodeDescr
            CHTMP = CodeDescr(i)
            Write(2,'(A)')CHTMP(1:LEN_TRIM(CHTMP))
         Enddo
*---  (here: neither data - nor multiplicative factors)
         Write(2,'(I2)') IRef   ! IRef
         IScaleDep = 0          ! in v1.4 tables LO always 1st
         If (n.eq.2) IScaleDep = 1 ! in v1.4 tables NLO always 2nd
         If (n.eq.3) IScaleDep = 2 ! in v1.4 tables thresh always 3rd
         Write(2,'(I2)') IScaleDep ! IScaleDep
         Write(2,'(I14)') Nevt(n) ! Nevt
         Write(2,'(I2)') Npow(n) ! Npow
         Write(2,'(I2)') NPDF   ! NPDF
         Do i=1,NPDF
            Write(2,'(I5)') NPDFPDG(i) ! NPDFPDG
         Enddo
         Write(2,'(I2)') NPDFDim ! NPDFdim
         Write(2,'(I2)') NFragFunc ! NFragFunc
         Write(2,'(I2)') 0      ! NFragDim

*---  NSubProc
         If (NPDFDim.eq.0) then ! ---- DIS
            If (n.eq.1) then    ! --- LO
               Write(2,'(I2)') 2
            Elseif (n.eq.2) Then
               Write(2,'(I2)') 3
            Endif
         ElseIf(NPDFDim.eq.1) then ! --------- pp
            If (n.eq.1 .or. n.eq.3) then ! ------ LO or threshcor
               Write(2,'(I2)') 6
            Elseif (n.eq.2) then ! ---- NLO
               Write(2,'(I2)') 7
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
            If (n.eq.1 .or. n.eq.3) Then ! ------ LO or threshcor
               IPDFdef1 = 3
               IPDFdef2 = 1
               IPDFdef3 = 1
            Elseif (n.eq.2) Then ! ---- NLO
               IPDFdef1 = 3
               IPDFdef2 = 1
               IPDFdef3 = 2
            Else
               Write(*,*) "# CV20WRT: WARNING! Contribution no. "
               Write(*,*) "#   outside 1,2,3 in pp table:", N
            Endif
         Endif
         Write(2,'(I5)') IPDFdef1 ! IPDFdef1
         Write(2,'(I5)') IPDFdef2 ! IPDFdef2
         Write(2,'(I5)') IPDFdef3 ! IPDFdef3

         i = 0
         Do i1=1,Nrapidity
            Do i2=1,NPT(i1)
               i = i+1          ! --- bin counter -> linear bin numbering
               Write(2,'(I5)') Nxtot ! Nxtot(1) - new for each ObsBin
               Do j=1,Nxtot
                  XNode = xlimit(i1,i2)
                  if (ixscheme.eq.2) then
                     hxlim = - sqrt(-log10(xnode))
                     hx = hxlim *(1d0 - dble(j-1)/dble(Nxtot))
                     XNode = 10**(-(hx*hx))
                  elseif (ixscheme.eq.1) then
                     hxlim = log10(xnode)
                     hx = hxlim *(1d0 - dble(j-1)/dble(Nxtot))
                     XNode = 10**(hx)
                  endif
                  xarray(i,j) = XNode ! --- use later for PDF unweighting
                  Write(2,'(G24.17)') XNode
               Enddo
            Enddo
         Enddo

         Write(2,'(I2)') NScales     ! NScales
         Write(2,'(I2)') NScaleDim   ! NScaleDim
         Do i=1,Nscales
            Write(2,'(I2)') 0        ! Iscale(i) <- all scales = 1st Dimension
         Enddo
         Do i=1,Nscaledim
            Write(2,'(I2)') 1        ! NScaleDescript(i) <- one descriptive string
            Write(2,'(A)') SCALELABEL(1:LEN_TRIM(SCALELABEL))
         Enddo

         If (IScaleDep.eq.0) then
            nscvar = 1
         Else
            nscvar = NscaleVar
         endif

         Do i=1,NscaleDim
            Write(2,'(I2)') nscvar   ! NscaleVar
            Write(2,'(I2)') NscaleBin ! NscaleBin
         Enddo
         Do i=1,Nscaledim
            if (nscvar.eq.1) then
               Write(2,'(G24.17)') 1d0   ! Scalefac(i)(j)
            else
               Do j=1,nscvar
                  Write(2,'(G24.17)') MURSCALE(j) ! Scalefac(i)(j)
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
                           Write(2,'(G24.17)')MURVAL(i1,i2,l) !ScaleNode(ijkl)
                        Else
                           Write(2,'(G24.17)')MURVAL(i1,i2,l) !ScaleNode(ijkl)
     >                          *MURSCALE(k)
                        Endif
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo

*---  undo PDF reweighting -> compute factors
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
                     Write(*,*)
     >                    "# CV20WRT: ERROR! Cannot deal with "//
     >                    "NPDFDim=2 here:", NPDFDim
                     STOP
                  Endif
                  Do k=1,mmax
                     m = m+1    ! Bin No. in linear x-array
*---  Works only for standard fastNLO PDF weighting formula
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


*---  Diagonalize Bernstein scale interpolation in v1.4
*---  Also convert the simpler linear interpolation with 2 scale bins
         If (ITabversion.eq.14000 .and. NScaleBin.ge.2) Then
            i = 0
            Do i1=1,Nrapidity
               Do i2=1,NPT(i1)
                  i = i+1       ! --- bin counter ObsBin
                  Do j=1,Nscaledim !                 scaledim
                     Do k=1,nscvar !                 scalevar
                        if (NPDFDim.eq.0) nxall=Nxtot
                        if (NPDFDim.eq.1) nxall=(Nxtot*Nxtot+Nxtot)/2
                        Do m=1,nxall !            xbins
                           Do n1=1,Nsubproc !     subprocesses
                              If (n.eq.1) then
                                 nscaddr = 1
                              Else
                                 nscaddr = 1+k+(n-2)*nscalevar
                              Endif

                              If (NScaleBin.eq.2) Then ! 2 Scalebins as in Dijet Mass
ckr Something to do here? Converted dijet mass tables of this type exhibit discrepancies!!!
                              ElseIf (NScaleBin.eq.3) Then ! 3 Bernstein Scalebins
                                 t1 = 1./2.
                                 Bnfactor =  1./(2.*(1.-t1)*t1)
                                 Bnst(1) = array(i,m,n1,nscaddr,1)
                                 Bnst(2) = array(i,m,n1,nscaddr,2)
                                 Bnst(3) = array(i,m,n1,nscaddr,3)
                                 array(i,m,n1,nscaddr,1) =
     +                                1d0*Bnst(1) - Bnfactor*(1.-t1)**2
     >                                *Bnst(2)
                                 array(i,m,n1,nscaddr,2) =
     +                                + Bnfactor            *Bnst(2)
     >
                                 array(i,m,n1,nscaddr,3) =
     +                                - Bnfactor*(t1)**2    *Bnst(2)+
     >                                1d0*Bnst(3)

                              Elseif (NScaleBin.eq.4) Then ! 4 Bernstein Scalebins
                                 t1 = 1./3.
                                 t2 = 2./3.
                                 Bnfactor =  1./(3.*(1.-t1)**2*t1 * 3.
     >                                *(1.-t2)*(t2)**2 - 3.*(1.-t2)**2
     >                                *t2 * 3.*(1.-t1)*(t1)**2)
                                 Bnst(1) = array(i,m,n1,nscaddr,1)
                                 Bnst(2) = array(i,m,n1,nscaddr,2)
                                 Bnst(3) = array(i,m,n1,nscaddr,3)
                                 Bnst(4) = array(i,m,n1,nscaddr,4)

                                 array(i,m,n1,nscaddr,1) =
     &                                1d0*Bnst(1)+
     &                                Bnfactor* (-3.*(1.-t2)*(t2)**2*(1.
     >                                -t1)**3+ 3.*(1.-t1)*(t1)**2*(1.
     >                                -t2)**3) *Bnst(2)+Bnfactor* (-3.
     >                                *(1.-t1)**2*(t1)*(1.-t2)**3+ 3.
     >                                *(1.-t2)**2*(t2)* *(1.-t1)**3)
     >                                *Bnst(3)

                                 array(i,m,n1,nscaddr,2) =
     &                                Bnfactor*(3.*(1.-t2)*(t2)**2)
     >                                *Bnst(2)+Bnfactor*(- 3.*(1.-t2)**2
     >                                *(t2)) *Bnst(3)

                                 array(i,m,n1,nscaddr,3) =
     &                                Bnfactor*(-3.*(1.-t1)*(t1)**2)
     >                                *Bnst(2) +Bnfactor*(3.*(1.-t1)**2
     >                                *(t1)) *Bnst(3)

                                 array(i,m,n1,nscaddr,4) =
     &                                Bnfactor* (-3.*(1.-t2)*(t2)**2
     >                                *(t1)**3+ 3.*(1.-t1)*(t1)**2*(t2)
     >                                **3)  *Bnst(2) +Bnfactor* (-3.*(1.
     >                                -t1)**2*(t1)*(t2)**3+ 3.*(1.-t2)
     >                                **2*(t2)*(t1)**3)  *Bnst(3) +1d0
     >                                *Bnst(4)

                              Else
                                 Write(*,*) "# CV20WRT: ERROR! "//
     >                                "Not implemented: NScaleBin = ",
     >                                NScaleBin
                                 Stop
                              Endif
                           Enddo
                        Enddo
                     Enddo
                  Enddo
               Enddo
            Enddo
         Endif                  ! Bernstein diagonalization

*---  Write SigmaTilde array
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
*---  Reorder subprocesses - use pointer
                              if (ireaction.eq.1) then ! DIS
                                 n1 = idispoint(n2)
                              else
                                 n1 = ipppoint(n2)
*---  n1=n2
                              endif
                              if (n.eq.1) then
                                 nscaddr = 1
                              else
                                 nscaddr = 1+k+(n-2)*nscalevar
                              endif
                              a = array(i,m,n1,nscaddr,l)
                              a = a*rewgt(i,m) ! --- PDF unweighting

*---  Special cases: In pp in LO merge 2nd & 3rd subproc for pp/ppbar
*---                 In DIS remove Sigma at LO

                              if (ireaction.eq.1) then ! --- DIS
                                 if (n.eq.1 .and. n1.eq.2) then
                                    continue
                                 else
                                    if (a.eq.0d0) then
                                       Write(2,'(A)') '0'
                                    else
                                       Write(2,'(G24.17)') a
                                    endif
                                 endif

                              else ! ------- pp/ppbar
                                 if ((n.eq.1 .or. n.eq.3).and. n1.eq.3)
     >                                a=(array(i,m,n1-1,nscaddr,l)
     >                                +array(i,m,n1,  nscaddr,l))/2d0
     >                                *rewgt(i,m)
                                 if ((n.eq.1 .or. n.eq.3) .and. n1.eq.2)
     >                                then
                                    continue
                                 else
                                    if (a.eq.0d0) then
                                       Write(2,'(A)') '0'
                                    else
                                       Write(2,'(G24.17)') a
                                    endif
                                 endif

*---  Old test for problem that has now diappeared
*---  Problem: for scalevar 1&2 threshcor subproc 2,3 not identical
                                 if (n.eq.3 .and. n1.eq.3 .and. a.ne.0d0
     >                                .and.(array(i,m,n1-1,nscaddr,l)
     >                                /array(i,m,n1,  nscaddr
     >                                ,l)).ne.1d0)Write(*,*)i,array(i,m
     >                                ,n1-1,nscaddr,l)/array(i,m,n1,
     >                                nscaddr,l),k,l,m,nscaddr,n

                              endif ! ---- DIS or pp/ppbar?
                           Enddo
                        Enddo
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo
         Write(2,'(I10)') Iseparator ! ------------------- Block B ------
      Enddo
      Write(2,'(I10)') Iseparator ! --------- Add. check for EOF --

 999  continue
      close(unit=2)

      RETURN
      END
