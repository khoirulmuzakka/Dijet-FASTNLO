      Subroutine fx9999rw(crw,filename,ictrb)
* ------------------------------------------------------------
*  fastNLO usercode v2.0                       MW 07/23/2007
* 
*  fx9999rw reads or writes fastNLO v2 tables
*  input:
*     character   crw:        'read' or 'write'
*     character   filename:    name of table
*     integer     ictrb(30):   pointers for writing selected contributions
*                              only contributions n with ictrb(n)<>0 are 
*                              written - reading is not affected
C
C        1         2         3         4         5         6         7 
C 3456789012345678901234567890123456789012345678901234567890123456789012
* ------------------------------------------------------------
      Implicit None
      Character*(*) crw,filename
      Integer ictrb(30)
      Integer Nunit,Ifile,ic,i,j,k,l,n,m,nxmax,nselctrb
      Integer IMult
      Character*1 CH1TMP,CH1TMP2,CH1TMP3,CH1TMP4
      Character*2 CH2TMP
      Character*3 CH3TMP
      Character*4 CH4TMP
      Character*40 CHTMP
      Character*41  CSEP41,DSEP41,LSEP41,SSEP41
      Character*82  CSEPS,DSEPS,LSEPS,SSEPS
      Character*164 CSEPL,DSEPL,LSEPL,SSEPL
      Integer iprint
      Logical lprint
      Include 'fnx9999.inc'
      Data CSEP41,DSEP41,LSEP41,SSEP41/
     >     '#########################################',
     >     '=========================================',
     >     "-----------------------------------------",
     >     "*****************************************"/

*---  Initialization
      CSEPS = CSEP41//CSEP41
      DSEPS = DSEP41//DSEP41
      LSEPS = LSEP41//LSEP41
      SSEPS = SSEP41//SSEP41
      CSEPL = CSEP41//CSEP41//CSEP41//CSEP41
      DSEPL = DSEP41//DSEP41//DSEP41//DSEP41
      LSEPL = LSEP41//LSEP41//LSEP41//LSEP41
      SSEPL = SSEP41//SSEP41//SSEP41//SSEP41
      DSEPS(1:2) = "# "
      LSEPS(1:2) = "# "
      SSEPS(1:2) = "# "
      DSEPL(1:2) = "# "
      LSEPL(1:2) = "# "
      SSEPL(1:2) = "# "
* iprint = 0: No additional printout
*          1: Print Block A1 & A2 (A1, A2)
*          2: Also print basic values of Block B (B0)
*          3: Also print x nodes of Block B for each contribution (BX)
*          4: Also print scale nodes of Block B for each contribution (BS)
*          5: Also print sigma tilde of Block B (not implemented yet)
      iprint = 2
      lprint = .false.
      Nunit=2

      If (crw.eq.'write') Then
         Open(unit=Nunit,file=filename,status='unknown')
      Else
         OPEN(Nunit,STATUS='OLD',FILE=FILENAME,IOSTAT=IFILE)
         IF (Ifile .ne. 0) THEN
            WRITE(*,*)"FX9999RW: ERROR! Table file "//
     >           FILENAME(1:LEN_TRIM(FILENAME))
            WRITE(*,*)"          not found, "//
     +           "stopped! IOSTAT = ",Ifile
            STOP
         Endif
      Endif

c --- when writing: count No. of selected contributions to be written
c                   (those with ictrb(ic)=1)
      nselctrb = 0
      If (crw .eq. 'write') Then         
         Do ic=1,NContrib
            If (ictrb(ic).ne.0) Then
               nselctrb = nselctrb + 1
               Write(*,*) '         ... writing contribution ',ic
            Endif
         Enddo
         Write(*,*) ' writing a total of ',nselctrb,' contributions',
     +        ' into: ',filename
      Endif


c --- fastNLO Table
c --- Block A1
      If (iprint.gt.0) then
         lprint = .true.
         write(*,'(A)')""
         write(*,*)SSEP41
         write(*,*)"* fastNLO Table: Block A1"
         write(*,*)SSEP41
      Endif
      Call fnioisep(crw,nunit, lprint,"  A1  ISep")
      Call fnioint(crw,nunit, Itabversion,lprint,"  A1  Itabversion")
      Call fniochar(crw,nunit, ScenName,lprint,"  A1  ScenName")
      If (crw .eq. 'write') Then
         Call fnioint(crw,nunit, nselctrb,lprint,"  A1  nselctrb")
      Else
         Call fnioint(crw,nunit, Ncontrib,lprint,"  A1  Ncontrib")
      Endif
      Call fnioint(crw,nunit, Nmult,lprint,"  A1  Nmult")
      Call fnioint(crw,nunit, Ndata,lprint,"  A1  Ndata")
      Call fnioint(crw,nunit, NuserString,lprint,"  A1  NuserString")
      If (NuserString .gt. MxUser) Then
         Write(*,*) ' NuserString too large ',NuserString,'<=',MxUser
      Endif
      Do i=1,NuserString
         Write(CH1TMP,'(I1)'),i
         CHTMP = "  A1    UserString("//CH1TMP//")"
         Call fniochar(crw,nunit, UserString(i),lprint,CHTMP)
      Enddo
      Call fnioint(crw,nunit, NuserInt,lprint,"  A1  NuserInt")
      If (NuserInt .gt. MxUser) Then
         Write(*,*) ' NuserInt too large ',NuserInt,'<=',MxUser
      Endif
      Do i=1,NuserInt
         Write(CH1TMP,'(I1)'),i
         CHTMP = "  A1    UserInt("//CH1TMP//")"
         Call fnioint(crw,nunit, UserInt(i),lprint,CHTMP)
      Enddo
      Call fnioint(crw,nunit, NuserFloat,lprint,"  A1  NuserFloat")
      If (NuserFloat .gt. MxUser) Then
         Write(*,*) ' NuserFloat too large ',NuserFloat,'<=',MxUser
      Endif
      Do i=1,NuserFloat
         Write(CH1TMP,'(I1)'),i
         CHTMP = "  A1    UserFloat("//CH1TMP//")"
         Call fniodbl(crw,nunit, UserFloat(i),lprint,CHTMP)
      Enddo
      Call fnioint(crw,nunit, Imachine,lprint,"  A1  Imachine") 
      If (Imachine.ne.0) Then
         Write(*,*) ' Imachine different from zero - ',
     +        'not yet implemented: ',Imachine
         Stop
      Endif
      If (iprint.gt.0) then
         write(*,*)CSEP41
      Endif
      lprint = .false.



c --- Block A2
      If (iprint.gt.0) then
         lprint = .true.
         write(*,'(A)')""
         write(*,*)SSEP41
         write(*,*)"* fastNLO Table: Block A2"
         write(*,*)SSEP41
      Endif
      Call fnioisep(crw,nunit,lprint,"  A2  ISep")
      Call fnioint(crw,nunit, IpublUnits,lprint,"  A2  IpublUnits") 
      Call fnioint(crw,nunit, NscDescript,lprint,"  A2  NscDescript")
      If (NscDescript .gt. MxScDescript) then
         write(*,*) ' NscDescript too large ',NscDescript,'<=',MXScDescript
      Endif
      Do i=1,NscDescript
         Write(CH1TMP,'(I1)'),i
         CHTMP = "  A2    ScDescript("//CH1TMP//")"
         Call fniochar(crw,nunit, ScDescript(i),lprint,CHTMP)
      Enddo
      Call fniodbl(crw,nunit, Ecms,lprint,"  A2  Ecms")
      Call fnioint(crw,nunit, ILOord,lprint,"  A2  ILOord")
      Call fnioint(crw,nunit, NobsBin,lprint,"  A2  NobsBin")
      Call fnioint(crw,nunit, NDim,lprint,"  A2  NDim")
      Do i=1,NDim
         Write(CH1TMP,'(I1)'),i
         CHTMP = "  A2    DimLabel("//CH1TMP//")"
         Call fniochar(crw,nunit, DimLabel(i),lprint,CHTMP)
      Enddo
      Do i=1,NDim
         Write(CH1TMP,'(I1)'),i
         CHTMP = "  A2    IDiffBin("//CH1TMP//")"
         Call fnioint(crw,nunit, IDiffBin(i),lprint,CHTMP)
      Enddo
      Do i=1,NObsBin
         Do j=1,NDim
            Write(CH3TMP,'(I3)'),i
            Write(CH1TMP,'(I1)'),j
            CHTMP = "  A2      LoBin("//
     >           CH3TMP//","//CH1TMP//")"
            Call fniodbl(crw,nunit, LoBin(i,j),lprint,CHTMP)
            If (IDiffBin(j).eq.2) then
               CHTMP = "  A2      UpBin("//
     >              CH3TMP//","//CH1TMP//")"
               Call fniodbl(crw,nunit, UpBin(i,j),lprint,CHTMP)
            Endif
         Enddo
      Enddo
      Do i=1,NObsBin
         Write(CH3TMP,'(I3)'),i
         CHTMP = "  A2    BinSize("//CH3TMP//")"
         Call fniodbl(crw,nunit, BinSize(i),lprint,CHTMP)
      Enddo
      Call fnioint(crw,nunit, INormFlag,lprint,"  A2  INormFlag")
      If (INormFlag.gt.1) then
         Call fniochar(crw,nunit, DenomTable,lprint,"  A2  DenomTable")
      Endif
      If (INormFlag.gt.0) Then
         Do i=1,NObsBin
            Write(CH3TMP,'(I3)'),i 
            CHTMP = "  A2    IDivLoPointer("//CH3TMP//")"
            Call fnioint(crw,nunit, IDivLoPointer(i),lprint,CHTMP)
            CHTMP = "  A2    IDivUpPointer("//CH3TMP//")"
            Call fnioint(crw,nunit, IDivUpPointer(i),lprint,CHTMP)
         Enddo
      Endif
      If (iprint.gt.0) then
         write(*,*)CSEP41
      Endif
      lprint = .false.



c --- Block B
      If (iprint.gt.1) then
         lprint = .true.
         write(*,'(A)')""
         write(*,*)SSEP41
         write(*,*)"* fastNLO Table: Block B"
         write(*,*)SSEP41
      Endif
c -   count number of multiplicative contributions
      IMult = 0
      Do ic=1,NContrib
c -      write only selected contributions - those with ictrb(ic)=1         
         If (crw .eq. 'write' .and. ictrb(ic).eq.0) Goto 100

         Write(CH1TMP,'(I1)'),ic
         Call fnioisep(crw,nunit,lprint,"  B0  ISep")
         CHTMP = "  B0    IXsectUnits("//CH1TMP//")"
         Call fnioint(crw,nunit, IXsectUnits(ic),lprint,CHTMP)
         CHTMP = "  B0    IDataFlag("//CH1TMP//")"
         Call fnioint(crw,nunit, IDataFlag(ic),lprint,CHTMP)
         CHTMP = "  B0    IAddMultFlag("//CH1TMP//")"
         Call fnioint(crw,nunit, IAddMultFlag(ic),lprint,CHTMP)
         CHTMP = "  B0    IContrFlag1("//CH1TMP//")"
         Call fnioint(crw,nunit, IContrFlag1(ic),lprint,CHTMP)
         CHTMP = "  B0    IContrFlag2("//CH1TMP//")"
         Call fnioint(crw,nunit, IContrFlag2(ic),lprint,CHTMP)
         CHTMP = "  B0    IContrFlag3("//CH1TMP//")"
         Call fnioint(crw,nunit, IContrFlag3(ic),lprint,CHTMP)
         CHTMP = "  B0    NContrDescr("//CH1TMP//")"
         Call fnioint(crw,nunit, NContrDescr(ic),lprint,CHTMP)
         Do i=1,NContrDescr(ic)
            Write(CH1TMP2,'(I1)'),i
            CHTMP = "  B0      CtrbDescript("
     >           //CH1TMP//","//CH1TMP2//")"
            Call fniochar(crw,nunit, CtrbDescript(ic,i),lprint,CHTMP)
         Enddo
         CHTMP = "  B0    NCodeDescr("//CH1TMP//")"
         Call fnioint(crw,nunit, NCodeDescr(ic),lprint,CHTMP)
         Do i=1,NCodeDescr(ic)
            Write(CH1TMP2,'(I1)'),i
            CHTMP = "  B0      CodeDescript("
     >           //CH1TMP//","//CH1TMP2//")"
            Call fniochar(crw,nunit, CodeDescript(ic,i),lprint,CHTMP)
         Enddo

c --------------------------- IAddMult, done (KR)
         If (IAddMultFlag(ic).eq.1) Then
            IMult = IMult + 1
            Write(CH1TMP2,'(I1)'),IMult
            CHTMP = "  B0    NMUncorrel("//CH1TMP2//")"
            Call fnioint(crw,nunit, NMUncorrel(IMult),lprint,CHTMP)
            Do i=1,NMUncorrel(IMult)
               Write(CH1TMP3,'(I1)'),i
               CHTMP = "  B0      MUncDescript("
     >              //CH1TMP2//","//CH1TMP3//")"
               Call fniochar(crw,nunit,
     >              MUncDescript(IMult,i),lprint,CHTMP)
            Enddo
            CHTMP = "  B0    NMCorrel("//CH1TMP2//")"
            Call fnioint(crw,nunit, NMCorrel(IMult),lprint,CHTMP)
            Do i=1,NMCorrel(IMult)
               Write(CH1TMP3,'(I1)'),i
               CHTMP = "  B0      MCorDescript("
     >              //CH1TMP2//","//CH1TMP3//")"
               Call fniochar(crw,nunit,
     >              MCorDescript(IMult,i),lprint,CHTMP)
            Enddo
            Do i=1,NObsbin
               Write(CH3TMP,'(I3)'),i
               CHTMP = "  B0      MFact("
     >              //CH1TMP2//","//CH3TMP//")"
               Call fniodbl(crw,nunit,
     >              MFact(IMult,i),lprint,CHTMP)
               Do j=1,NMUncorrel(IMult)
                  Write(CH1TMP4,'(I1)'),j
                  CHTMP = "  B0        MUnCorLo("
     >                 //CH1TMP2//","//CH3TMP//","//CH1TMP4//")"
                  Call fniodbl(crw,nunit,
     >                 MUnCorLo(IMult,i,j),lprint,CHTMP)
                  CHTMP = "  B0        MUnCorUp("
     >                 //CH1TMP2//","//CH3TMP//","//CH1TMP4//")"
                  Call fniodbl(crw,nunit,
     >                 MUnCorUp(IMult,i,j),lprint,CHTMP)
               Enddo
               Do j=1,NMCorrel(IMult)
                  Write(CH1TMP4,'(I1)'),j
                  CHTMP = "  B0        MCorLo("
     >                 //CH1TMP2//","//CH3TMP//","//CH1TMP4//")"
                  Call fniodbl(crw,nunit,
     >                 MCorLo(IMult,i,j),lprint,CHTMP)
                  CHTMP = "  B0        MCorUp("
     >                 //CH1TMP2//","//CH3TMP//","//CH1TMP4//")"
                  Call fniodbl(crw,nunit,
     >                 MCorUp(IMult,i,j),lprint,CHTMP)
               Enddo
            Enddo
            Goto 100
         Endif

c --------------------------- Idata, done (KR)
         If (IDataFlag(ic).eq.1) Then
            CHTMP = "  B0    NDUncorrel"
            Call fnioint(crw,nunit, NDUncorrel,lprint,CHTMP)
            Do i=1,NDUncorrel
               Write(CH1TMP3,'(I1)'),i
               CHTMP = "  B0      DUncDescript("
     >              //CH1TMP3//")"
               Call fniochar(crw,nunit,
     >              DUncDescript(i),lprint,CHTMP)
            Enddo
            CHTMP = "  B0    NDCorrel"
            Call fnioint(crw,nunit, NDCorrel,lprint,CHTMP)
            Do i=1,NDCorrel
               Write(CH1TMP3,'(I1)'),i
               CHTMP = "  B0      DCorDescript("
     >              //CH1TMP3//")"
               Call fniochar(crw,nunit,
     >              DCorDescript(i),lprint,CHTMP)
            Enddo
            Do i=1,NObsbin
               Write(CH3TMP,'(I3)'),i
               CHTMP = "  B0      DxVal("
     >              //CH3TMP//")"
               Call fniodbl(crw,nunit,
     >              DxVal(i),lprint,CHTMP)
               CHTMP = "  B0      DyVal("
     >              //CH3TMP//")"
               Call fniodbl(crw,nunit,
     >              DyVal(i),lprint,CHTMP)
ckr Use write out for changing order Up Lo to Lo Up and
ckr rel. to abs. if required for table
Comment:                write(*,'(A,G10.4)')"ABSX: ",DxVal(i)
Comment:                write(*,'(A,G10.4)')"ABSY: ",DyVal(i)
               Do j=1,NDUncorrel
                  Write(CH1TMP4,'(I1)'),j
                  CHTMP = "  B0        DUnCorLo("
     >                 //CH3TMP//","//CH1TMP4//")"
                  Call fniodbl(crw,nunit,
     >                 DUnCorLo(i,j),lprint,CHTMP)
                  CHTMP = "  B0        DUnCorUp("
     >                 //CH3TMP//","//CH1TMP4//")"
                  Call fniodbl(crw,nunit,
     >                 DUnCorUp(i,j),lprint,CHTMP)
Comment:                   write(*,'(A,G10.4)')"ABSUPU: ",DUnCorUp(i,j)*DyVal(i)/100.D0
Comment:                   write(*,'(A,G10.4)')"ABSLOU: ",DUnCorLo(i,j)*DyVal(i)/100.D0
               Enddo
               Do j=1,NDCorrel
                  Write(CH1TMP4,'(I1)'),j
                  CHTMP = "  B0        DCorLo("
     >                 //CH3TMP//","//CH1TMP4//")"
                  Call fniodbl(crw,nunit,
     >                 DCorLo(i,j),lprint,CHTMP)
                  CHTMP = "  B0        DCorUp("
     >                 //CH3TMP//","//CH1TMP4//")"
                  Call fniodbl(crw,nunit,
     >                 DCorUp(i,j),lprint,CHTMP)
Comment:                   write(*,'(A,G10.4)')"ABSUPC: ",DCorUp(i,j)*DyVal(i)/100.D0
Comment:                   write(*,'(A,G10.4)')"ABSLOC: ",DCorLo(i,j)*DyVal(i)/100.D0
               Enddo
            Enddo
            CHTMP = "  B0    NDErrMatrix"
            Call fnioint(crw,nunit, NDErrMatrix,lprint,CHTMP)
            If (NDErrMatrix.NE.0) Then
               Write(*,*)"fx9999rw: Data with error matrix not"//
     >              " yet implemented, aborted! NDErrMatrix =",
     >              NDErrMatrix
               Stop
            Endif
            Goto 100
         Endif

c --- coefficient block
         CHTMP = "  B0    IRef("//CH1TMP//")"
         Call fnioint(crw,nunit, IRef(ic),lprint,CHTMP)
         CHTMP = "  B0    IScaleDep("//CH1TMP//")"
         Call fnioint(crw,nunit, IScaleDep(ic),lprint,CHTMP)
         CHTMP = "  B0    Nevt("//CH1TMP//")"
         Call fniolint(crw,nunit, Nevt(ic),lprint,CHTMP)
         CHTMP = "  B0    Npow("//CH1TMP//")"
         Call fnioint(crw,nunit, Npow(ic),lprint,CHTMP)
         CHTMP = "  B0    NPDF("//CH1TMP//")"
         Call fnioint(crw,nunit, NPDF(ic),lprint,CHTMP)
         Do i=1,NPDF(ic)
            Write(CH1TMP2,'(I1)'),i
            CHTMP = "  B0      NPDFPDG("//CH1TMP//","//CH1TMP2//")"
            Call fnioint(crw,nunit, NPDFPDG(ic,i),lprint,CHTMP)
         Enddo
         CHTMP = "  B0    NPDFDim("//CH1TMP//")"
         Call fnioint(crw,nunit, NPDFDim(ic),lprint,CHTMP)
         CHTMP = "  B0    NFragFunc("//CH1TMP//")"
         Call fnioint(crw,nunit, NFragFunc(ic),lprint,CHTMP)
         Do i=1,NFragFunc(ic)
            Write(CH1TMP2,'(I1)'),i
            CHTMP = "  B0     NFFPDG("//CH1TMP//","//CH1TMP2//")"
            Call fnioint(crw,nunit, NFFPDG(ic,i),lprint,CHTMP)
         Enddo
         CHTMP = "  B0    NFFDim("//CH1TMP//")"
         Call fnioint(crw,nunit, NFFDim(ic),lprint,CHTMP)
         CHTMP = "  B0    NSubproc("//CH1TMP//")"
         Call fnioint(crw,nunit, NSubproc(ic),lprint,CHTMP)
         CHTMP = "  B0    IPDFdef1("//CH1TMP//")"
         Call fnioint(crw,nunit, IPDFdef(ic,1),lprint,CHTMP)
         CHTMP = "  B0    IPDFdef2("//CH1TMP//")"
         Call fnioint(crw,nunit, IPDFdef(ic,2),lprint,CHTMP)
         CHTMP = "  B0    IPDFdef3("//CH1TMP//")"
         Call fnioint(crw,nunit, IPDFdef(ic,3),lprint,CHTMP)

c --- check consistency between Nsubproc and IPDFdef1,2,
         If (IPDFdef(ic,1).eq.0) Then ! - no predefined set of PDF coefficients
            Continue
         Elseif (IPDFdef(ic,1).eq.1) Then ! --- e+e-
            Continue
         Elseif (IPDFdef(ic,1).eq.2) Then ! --- DIS
            Continue
         Elseif (IPDFdef(ic,1).eq.3) Then ! --- hh/hhbar
            If (IPDFdef(ic,2).eq.1) Then 
               If (IPDFdef(ic,3).eq.1) Then
                  If (Nsubproc(ic).ne.6) Then
                     Write(*,*) " IPDFdef(3)=1 requires ",
     +                    "6 subprocesses for hh"
                     STOP
                  Endif
               Elseif (IPDFdef(ic,3).eq.2) Then 
                  If (Nsubproc(ic).ne.7) Then
                     Write(*,*) " IPDFdef(3)=2 requires ",
     +                    "7 subprocesses for hh"
                     STOP
                  Endif
               Endif
            Else
               Write(*,*) " case IPDFdef(2)<>1 n/a for hh/hhbar"
               STOP
            Endif

         Endif

         If (IPDFdef(ic,1).eq.0) Then ! - no predefined set of PDF coefficients
            Write(*,*) " case IPDFdef(1)=0 not yet implemented"
            STOP
         Endif
         
         lprint = .false.
         If (iprint.gt.2) then
            lprint = .true.
            write(*,'(A)')""
            write(*,*)SSEP41
            write(*,*)"* fastNLO Table: Block B ic = ",ic
            write(*,*)"*    X Node Details"
            write(*,*)SSEP41
         Endif

         If (NPDF(ic).gt.0) Then
            
            Do i=1,NObsBin
               Write(CH3TMP,'(I3)'),i
               CHTMP = "  BX        Nxtot1("//
     >              CH1TMP//",1,"//CH3TMP//")"
               Call fnioint(crw,nunit, Nxtot(ic,1,i),lprint,CHTMP)
               Do j=1,Nxtot(ic,1,i)
                  Write(CH2TMP,'(I2)'),j 
                  CHTMP = "  BX        XNode1("//
     >                 CH1TMP//","//CH3TMP//","//CH2TMP//")"
                  Call fniodbl(crw,nunit, XNode1(ic,i,j),lprint,CHTMP)
               Enddo
               If (crw.eq.'read') Then
                  Hxlim1(ic,i) = -sqrt(-log10(XNode1(ic,i,1))) ! for HH
               Endif
               
            Enddo 
            If (NPDFDim(ic).eq.2) Then
               
               Do i=1,NObsBin
                  Write(CH3TMP,'(I3)'),i
                  CHTMP = "  BX        Nxtot2("//
     >                 CH1TMP//",2,"//CH3TMP//")"
                  Call fnioint(crw,nunit, Nxtot(ic,2,i),lprint,CHTMP)
                  Do j=1,Nxtot(ic,2,i)
                     Write(CH2TMP,'(I2)'),j 
                     CHTMP = "  BX        XNode2("//
     >                    CH1TMP//","//CH3TMP//","//CH2TMP//")"
                     Call fniodbl(crw,nunit, XNode2(ic,i,j),lprint
     >                    ,CHTMP)
                  Enddo
               Enddo 

            Endif
         Endif
         IF (NFragFunc(ic).gt.0) then ! - no FFs so far
            write(*,*) " fastNLO: no FragFuncs so far"
            STOP
         Endif

         If (iprint.gt.2) then
            write(*,*)CSEP41
         Endif

         lprint = .false.
         if (iprint.gt.1) lprint=.true.

         CHTMP = "  B0    NScales("//CH1TMP//")"
         Call fnioint(crw,nunit, NScales(ic),lprint,CHTMP)
         CHTMP = "  B0    NScaleDim("//CH1TMP//")"
         Call fnioint(crw,nunit, NScaleDim(ic),lprint,CHTMP)

         Do i=1,NScales(ic)
            Write(CH1TMP2,'(I1)'),i
            CHTMP = "  B0      NScales("//
     >           CH1TMP//","//CH1TMP2//")"
            Call fnioint(crw,nunit, IScale(ic,i),lprint,CHTMP)
         Enddo
         Do i=1,NScaleDim(ic)
            Write(CH1TMP2,'(I1)'),i
            CHTMP = "  B0      NScaleDescript("//
     >           CH1TMP//","//CH1TMP2//")"
            Call fnioint(crw,nunit, NScaleDescript(ic,i),lprint,CHTMP)
            Do j=1,NScaleDescript(ic,i)
               Write(CH1TMP3,'(I1)'),j
               CHTMP = "  B0        ScaleDescript("//
     >              CH1TMP//","//CH1TMP2//","//CH1TMP3//")"
               Call fniochar(crw,nunit,
     >              ScaleDescript(ic,i,j),lprint,CHTMP)
            Enddo
         Enddo
         Do i=1,NScaleDim(ic)
            Write(CH1TMP2,'(I1)'),i
            CHTMP = "  B0      NScaleVar("//
     >           CH1TMP//","//CH1TMP2//")"
            Call fnioint(crw,nunit, NScaleVar(ic,i),lprint,CHTMP)
            CHTMP = "  B0      NScaleNode("//
     >           CH1TMP//","//CH1TMP2//")"
            Call fnioint(crw,nunit, NScaleNode(ic,i),lprint,CHTMP)
         Enddo

         lprint = .false.
         If (iprint.gt.3) then
            lprint = .true.
            write(*,'(A)')""
            write(*,*)SSEP41
            write(*,*)"* fastNLO Table: Block B ic = ",ic
            write(*,*)"*    Scale Node Details"
            write(*,*)SSEP41
         Endif

         Do i=1,NScaleDim(ic)
            Do j=1,NScaleVar(ic,i)
               Write(CH1TMP2,'(I1)'),i
               Write(CH1TMP3,'(I1)'),j
               CHTMP = "  BS        ScaleFac("//
     >              CH1TMP//","//CH1TMP2//","//CH1TMP3//")"
               Call fniodbl(crw,nunit, ScaleFac(ic,i,j),lprint,CHTMP)
            Enddo
         Enddo

         Do i=1,NObsBin
            Do j=1,NScaleDim(ic)
               Do k=1,NScaleVar(ic,j)
                  Do l=1,NScaleNode(ic,j)
                     Write(CH3TMP,'(I3)'),i
                     Write(CH1TMP2,'(I1)'),j
                     Write(CH1TMP3,'(I1)'),k
                     Write(CH1TMP4,'(I1)'),l
                     CHTMP = "  BS          ScaleNode("//
Comment:      >                    CH1TMP//","//
     >                    CH3TMP//","//CH1TMP2//","//
     >                    CH1TMP3//","//CH1TMP4//")"
                     Call fniodbl(crw,nunit, ScaleNode(ic,i,j,k,l),
     >                    lprint,CHTMP)
                     If (crw.eq.'read') Then
                        HScaleNode(ic,i,j,k,l) = 
     +                       log(log(ScaleNode(ic,i,j,k,l)/0.25d0))
                     Endif
                  Enddo
               Enddo
            Enddo
         Enddo
         
         If (iprint.gt.3) then
            write(*,*)CSEP41
         Endif

         lprint = .false.
         If (iprint.gt.4) then
            lprint = .true.
            write(*,'(A)')""
            write(*,*)SSEP41
            write(*,*)"* fastNLO Table: Block B ic = ",ic
            write(*,*)"*    Sigma Tilde"
            write(*,*)"* Not implemented yet!"
            write(*,*)SSEP41
            lprint = .false.
         Endif

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
                        Call fniodbl(crw,nunit,
     +                       SigmaTilde(ic,i,1,k,l,m,n),lprint,CHTMP)
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo

         If (iprint.gt.4) then
            write(*,*)CSEP41
         Endif

         lprint = .false.
         if (iprint.gt.1) lprint=.true.
         
 100     Continue               ! - end of block
      Enddo

      If (iprint.gt.1) then
         write(*,*)CSEP41
      Endif
      lprint = .false.

c - end of table
      Call fnioisep(crw,nunit,lprint," T   ISep")
      Call fnioisep(crw,nunit,lprint," F   ISep")
      Close(2)

      Return
      End
C -------------------------------------------------------------
