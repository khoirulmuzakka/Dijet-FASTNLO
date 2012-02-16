      SUBROUTINE FX9999RW(CRW,FILENAME,ICTRB)
***********************************************************************
*     
*     Read and write complete (or partial) fastNLO v2 tables
*     
*     Input:
*     ------
*     CRW        'read' or 'write'
*     FILENAME    name of table
*     ICTRB(30)   flags for writing selected contributions
*     >           only contributions n with ictrb(n) <> 0 are 
*     >           written - reading is not affected
*     
*     
*     Current restrictions:
*     ---------------------
*     - The Info and debug screen print out is steered by the hardcoded
*     > IPRINT variable, default is IPRINT = 0, i.e. no additional
*     > print out
*     
***********************************************************************
      IMPLICIT NONE
      CHARACTER*(*) CRW,FILENAME
      INTEGER ICTRB(30)
      INTEGER NUNIT,IFILE,IC,I,J,K,L,N,M,NXMAX,NSELCTRB
      INTEGER IMULT
      CHARACTER*1 CH1TMP,CH1TMP2,CH1TMP3,CH1TMP4
      CHARACTER*2 CH2TMP
      CHARACTER*3 CH3TMP
      CHARACTER*4 CH4TMP
      CHARACTER*40 CHTMP
      INTEGER IPRINT
      LOGICAL LPRINT
      INCLUDE 'fnx9999.inc'
      INCLUDE 'strings.inc'

*---  Initialization
*     iprint = 0: No additional printout
*     >        1: Print Block A1 & A2 (A1, A2)
*     >        2: Also print basic values of Block B (B0)
*     >        3: Also print x nodes of Block B for each contr. (BX)
*     >        4: Also print scale nodes of Block B for each contr. (BS)
*     >        5: Also print sigma tilde of Block B (not implemented)
      IPRINT = 0
      LPRINT = .FALSE.
      NUNIT=2

      IF (CRW.EQ.'write') THEN
         OPEN(UNIT=NUNIT,FILE=FILENAME,STATUS='unknown')
      ELSE
         OPEN(NUNIT,STATUS='OLD',FILE=FILENAME,IOSTAT=IFILE)
         IF (IFILE .NE. 0) THEN
            WRITE(*,*)"FX9999RW: ERROR! Table file ",
     >           FILENAME(1:LEN_TRIM(FILENAME))
            WRITE(*,*)"          not found, "//
     >           "stopped! IOSTAT = ",Ifile
            STOP
         ENDIF
      ENDIF

*---  When writing count no. of selected contributions to be written
*     (those with ICTRB(IC) = 1)
      NSELCTRB = 0
      IF (CRW .EQ. 'write') THEN
         DO IC=1,NCONTRIB
            IF (ICTRB(IC).NE.0) THEN
               NSELCTRB = NSELCTRB + 1
               WRITE(*,*)'FX9999RW: ... writing contribution ',IC
            ENDIF
         ENDDO
         WRITE(*,*)'FX9999RW: Writing a total of ',NSELCTRB,
     >        ' contributions',' into: ',FILENAME
      ENDIF


*---  fastNLO Table
*---  Block A1
      IF (IPRINT.GT.0) THEN
         LPRINT = .TRUE.
         WRITE(*,'(A)')""
         WRITE(*,*)SSEP0
         WRITE(*,*)"* fastNLO Table: Block A1"
         WRITE(*,*)SSEP0
      ENDIF
      CALL FNIOISEP(CRW,NUNIT, LPRINT,"  A1  ISep")
      CALL FNIOINT(CRW,NUNIT, ITABVERSION,LPRINT,"  A1  Itabversion")
      CALL FNIOCHAR(CRW,NUNIT, SCENNAME,LPRINT,"  A1  ScenName")
*---  Ncontrib: All contributions including additive, multiplicative or
*---  data (Then: Nadd = Ncontrib - Nmult - Ndata)
      IF (CRW .EQ. 'write') THEN
         CALL FNIOINT(CRW,NUNIT, NSELCTRB,LPRINT,"  A1  nselctrb")
      ELSE
         CALL FNIOINT(CRW,NUNIT, NCONTRIB,LPRINT,"  A1  Ncontrib")
      ENDIF
      CALL FNIOINT(CRW,NUNIT, NMULT,LPRINT,"  A1  Nmult")
      CALL FNIOINT(CRW,NUNIT, NDATA,LPRINT,"  A1  Ndata")
      CALL FNIOINT(CRW,NUNIT, NUSERSTRING,LPRINT,"  A1  NuserString")
      IF (NUSERSTRING .GT. MXUSER) THEN
         WRITE(*,*)'FX9999RW: Warning! NuserString too large ',
     >        NUSERSTRING,'<=',MXUSER
      ENDIF
      DO I=1,NUSERSTRING
         WRITE(CH1TMP,'(I1)'),I
         CHTMP = "  A1    UserString("//CH1TMP//")"
         CALL FNIOCHAR(CRW,NUNIT, USERSTRING(I),LPRINT,CHTMP)
      ENDDO
      CALL FNIOINT(CRW,NUNIT, NUSERINT,LPRINT,"  A1  NuserInt")
      IF (NUSERINT .GT. MXUSER) THEN
         WRITE(*,*)'FX9999RW: Warning! NuserInt too large ',
     >        NUSERINT,'<=',MXUSER
      ENDIF
      DO I=1,NUSERINT
         WRITE(CH1TMP,'(I1)'),I
         CHTMP = "  A1    UserInt("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, USERINT(I),LPRINT,CHTMP)
      ENDDO
      CALL FNIOINT(CRW,NUNIT, NUSERFLOAT,LPRINT,"  A1  NuserFloat")
      IF (NUSERFLOAT .GT. MXUSER) THEN
         WRITE(*,*)'FX9999RW: Warning! NuserFloat too large ',
     >        NUSERFLOAT,'<=',MXUSER
      ENDIF
      DO I=1,NUSERFLOAT
         WRITE(CH1TMP,'(I1)'),I
         CHTMP = "  A1    UserFloat("//CH1TMP//")"
         CALL FNIODBL(CRW,NUNIT, USERFLOAT(I),LPRINT,CHTMP)
      ENDDO
      CALL FNIOINT(CRW,NUNIT, IMACHINE,LPRINT,"  A1  Imachine") 
      IF (IMACHINE.NE.0) THEN
         WRITE(*,*)'FX9999RW: Warning! Imachine different from zero - ',
     >        'not yet implemented: ',IMACHINE
         STOP
      ENDIF
      IF (IPRINT.GT.0) THEN
         WRITE(*,*)CSEP0
      ENDIF
      LPRINT = .FALSE.

*---  Block A2
      IF (IPRINT.GT.0) THEN
         LPRINT = .TRUE.
         WRITE(*,'(A)')""
         WRITE(*,*)SSEP0
         WRITE(*,*)"* fastNLO Table: Block A2"
         WRITE(*,*)SSEP0
      ENDIF
      CALL FNIOISEP(CRW,NUNIT,LPRINT,"  A2  ISep")
      CALL FNIOINT(CRW,NUNIT, IPUBLUNITS,LPRINT,"  A2  IpublUnits") 
      CALL FNIOINT(CRW,NUNIT, NSCDESCRIPT,LPRINT,"  A2  NscDescript")
      IF (NSCDESCRIPT .GT. MXSCDESCRIPT) THEN
         WRITE(*,*)'FX9999RW: Warning! NscDescript too large ',
     >        NSCDESCRIPT,'<=',MXSCDESCRIPT
      ENDIF
      DO I=1,NSCDESCRIPT
         WRITE(CH1TMP,'(I1)'),I
         CHTMP = "  A2    ScDescript("//CH1TMP//")"
         CALL FNIOCHAR(CRW,NUNIT, SCDESCRIPT(I),LPRINT,CHTMP)
      ENDDO
      CALL FNIODBL(CRW,NUNIT, ECMS,LPRINT,"  A2  Ecms")
      CALL FNIOINT(CRW,NUNIT, ILOORD,LPRINT,"  A2  ILOord")
      CALL FNIOINT(CRW,NUNIT, NOBSBIN,LPRINT,"  A2  NobsBin")
      CALL FNIOINT(CRW,NUNIT, NDIM,LPRINT,"  A2  NDim")
      DO I=1,NDIM
         WRITE(CH1TMP,'(I1)'),I
         CHTMP = "  A2    DimLabel("//CH1TMP//")"
         CALL FNIOCHAR(CRW,NUNIT, DIMLABEL(I),LPRINT,CHTMP)
      ENDDO
      DO I=1,NDIM
         WRITE(CH1TMP,'(I1)'),I
         CHTMP = "  A2    IDiffBin("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, IDIFFBIN(I),LPRINT,CHTMP)
      ENDDO
      DO I=1,NOBSBIN
         DO J=1,NDIM
            WRITE(CH3TMP,'(I3)'),I
            WRITE(CH1TMP,'(I1)'),J
            CHTMP = "  A2      LoBin("//
     >           CH3TMP//","//CH1TMP//")"
            CALL FNIODBL(CRW,NUNIT, LOBIN(I,J),LPRINT,CHTMP)
            IF (IDIFFBIN(J).EQ.2) THEN
               CHTMP = "  A2      UpBin("//
     >              CH3TMP//","//CH1TMP//")"
               CALL FNIODBL(CRW,NUNIT, UPBIN(I,J),LPRINT,CHTMP)
            ENDIF
         ENDDO
      ENDDO
      DO I=1,NOBSBIN
         WRITE(CH3TMP,'(I3)'),I
         CHTMP = "  A2    BinSize("//CH3TMP//")"
         CALL FNIODBL(CRW,NUNIT, BINSIZE(I),LPRINT,CHTMP)
      ENDDO
      CALL FNIOINT(CRW,NUNIT, INORMFLAG,LPRINT,"  A2  INormFlag")
      IF (INORMFLAG.GT.1) THEN
         CALL FNIOCHAR(CRW,NUNIT, DENOMTABLE,LPRINT,"  A2  DenomTable")
      ENDIF
      IF (INORMFLAG.GT.0) THEN
         DO I=1,NOBSBIN
            WRITE(CH3TMP,'(I3)'),I 
            CHTMP = "  A2    IDivLoPointer("//CH3TMP//")"
            CALL FNIOINT(CRW,NUNIT, IDIVLOPOINTER(I),LPRINT,CHTMP)
            CHTMP = "  A2    IDivUpPointer("//CH3TMP//")"
            CALL FNIOINT(CRW,NUNIT, IDIVUPPOINTER(I),LPRINT,CHTMP)
         ENDDO
      ENDIF
      IF (IPRINT.GT.0) THEN
         WRITE(*,*)CSEP0
      ENDIF
      LPRINT = .FALSE.

*---  Block B
      IF (IPRINT.GT.1) THEN
         LPRINT = .TRUE.
         WRITE(*,'(A)')""
         WRITE(*,*)SSEP0
         WRITE(*,*)"* fastNLO Table: Block B"
         WRITE(*,*)SSEP0
      ENDIF
*---  Count number of multiplicative contributions
      IMULT = 0
      DO IC=1,NCONTRIB
*---  Write only selected contributions - those with ICTRB(IC)=1
         IF (CRW .EQ. 'write' .AND. ICTRB(IC).EQ.0) GOTO 100
         WRITE(CH1TMP,'(I1)'),IC
         CALL FNIOISEP(CRW,NUNIT,LPRINT,"  B0  ISep")
         CHTMP = "  B0    IXsectUnits("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, IXSECTUNITS(IC),LPRINT,CHTMP)
         CHTMP = "  B0    IDataFlag("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, IDATAFLAG(IC),LPRINT,CHTMP)
         CHTMP = "  B0    IAddMultFlag("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, IADDMULTFLAG(IC),LPRINT,CHTMP)
         CHTMP = "  B0    IContrFlag1("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, ICONTRFLAG1(IC),LPRINT,CHTMP)
         CHTMP = "  B0    IContrFlag2("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, ICONTRFLAG2(IC),LPRINT,CHTMP)
         CHTMP = "  B0    NScaleDep("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, NSCALEDEP(IC),LPRINT,CHTMP)
         CHTMP = "  B0    NContrDescr("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, NCONTRDESCR(IC),LPRINT,CHTMP)
         DO I=1,NCONTRDESCR(IC)
            WRITE(CH1TMP2,'(I1)'),I
            CHTMP = "  B0      CtrbDescript("
     >           //CH1TMP//","//CH1TMP2//")"
            CALL FNIOCHAR(CRW,NUNIT, CTRBDESCRIPT(IC,I),LPRINT,CHTMP)
         ENDDO
         CHTMP = "  B0    NCodeDescr("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, NCODEDESCR(IC),LPRINT,CHTMP)
         DO I=1,NCODEDESCR(IC)
            WRITE(CH1TMP2,'(I1)'),I
            CHTMP = "  B0      CodeDescript("
     >           //CH1TMP//","//CH1TMP2//")"
            CALL FNIOCHAR(CRW,NUNIT, CODEDESCRIPT(IC,I),LPRINT,CHTMP)
         ENDDO
*---  IAddMult
         IF (IADDMULTFLAG(IC).EQ.1) THEN
            IMULT = IMULT + 1
            WRITE(CH1TMP2,'(I1)'),IMULT
            CHTMP = "  B0    NMUncorrel("//CH1TMP2//")"
            CALL FNIOINT(CRW,NUNIT, NMUNCORREL(IMULT),LPRINT,CHTMP)
            DO I=1,NMUNCORREL(IMULT)
               WRITE(CH1TMP3,'(I1)'),I
               CHTMP = "  B0      MUncDescript("
     >              //CH1TMP2//","//CH1TMP3//")"
               CALL FNIOCHAR(CRW,NUNIT,
     >              MUNCDESCRIPT(IMULT,I),LPRINT,CHTMP)
            ENDDO
            CHTMP = "  B0    NMCorrel("//CH1TMP2//")"
            CALL FNIOINT(CRW,NUNIT, NMCORREL(IMULT),LPRINT,CHTMP)
            DO I=1,NMCORREL(IMULT)
               WRITE(CH1TMP3,'(I1)'),I
               CHTMP = "  B0      MCorDescript("
     >              //CH1TMP2//","//CH1TMP3//")"
               CALL FNIOCHAR(CRW,NUNIT,
     >              MCORDESCRIPT(IMULT,I),LPRINT,CHTMP)
            ENDDO
            DO I=1,NOBSBIN
               WRITE(CH3TMP,'(I3)'),I
               CHTMP = "  B0      MFact("
     >              //CH1TMP2//","//CH3TMP//")"
               CALL FNIODBL(CRW,NUNIT,
     >              MFACT(IMULT,I),LPRINT,CHTMP)
               DO J=1,NMUNCORREL(IMULT)
                  WRITE(CH1TMP4,'(I1)'),j
                  CHTMP = "  B0        MUnCorLo("
     >                 //CH1TMP2//","//CH3TMP//","//CH1TMP4//")"
                  CALL FNIODBL(CRW,NUNIT,
     >                 MUNCORLO(IMULT,I,J),LPRINT,CHTMP)
                  CHTMP = "  B0        MUnCorUp("
     >                 //CH1TMP2//","//CH3TMP//","//CH1TMP4//")"
                  CALL FNIODBL(CRW,NUNIT,
     >                 MUNCORUP(IMULT,I,J),LPRINT,CHTMP)
               ENDDO
               DO J=1,NMCORREL(IMULT)
                  WRITE(CH1TMP4,'(I1)'),J
                  CHTMP = "  B0        MCorLo("
     >                 //CH1TMP2//","//CH3TMP//","//CH1TMP4//")"
                  CALL FNIODBL(CRW,NUNIT,
     >                 MCORLO(IMULT,I,J),LPRINT,CHTMP)
                  CHTMP = "  B0        MCorUp("
     >                 //CH1TMP2//","//CH3TMP//","//CH1TMP4//")"
                  CALL FNIODBL(CRW,NUNIT,
     >                 MCORUP(IMULT,I,J),LPRINT,CHTMP)
               ENDDO
            ENDDO
            GOTO 100
         ENDIF
*---  Idata
         IF (IDATAFLAG(IC).EQ.1) THEN
            CHTMP = "  B0    NDUncorrel"
            CALL FNIOINT(CRW,NUNIT, NDUNCORREL,LPRINT,CHTMP)
            DO I=1,NDUNCORREL
               WRITE(CH1TMP3,'(I1)'),i
               CHTMP = "  B0      DUncDescript("
     >              //CH1TMP3//")"
               CALL FNIOCHAR(CRW,NUNIT,
     >              DUNCDESCRIPT(I),LPRINT,CHTMP)
            ENDDO
            CHTMP = "  B0    NDCorrel"
            CALL FNIOINT(CRW,NUNIT, NDCORREL,LPRINT,CHTMP)
            DO I=1,NDCORREL
               WRITE(CH1TMP3,'(I1)'),I
               CHTMP = "  B0      DCorDescript("
     >              //CH1TMP3//")"
               CALL FNIOCHAR(CRW,NUNIT,
     >              DCORDESCRIPT(I),LPRINT,CHTMP)
            ENDDO
            DO I=1,NOBSBIN
               WRITE(CH3TMP,'(I3)'),I
               CHTMP = "  B0      DxVal("
     >              //CH3TMP//")"
               CALL FNIODBL(CRW,NUNIT,
     >              DXVAL(I),LPRINT,CHTMP)
               CHTMP = "  B0      DyVal("
     >              //CH3TMP//")"
               CALL FNIODBL(CRW,NUNIT,
     >              DYVAL(I),LPRINT,CHTMP)
C---  WRITE(*,'(A,G10.4)')"ABSX: ",DxVal(i)
C---  WRITE(*,'(A,G10.4)')"ABSY: ",DyVal(i)
               DO J=1,NDUNCORREL
                  WRITE(CH1TMP4,'(I1)'),J
                  CHTMP = "  B0        DUnCorLo("
     >                 //CH3TMP//","//CH1TMP4//")"
                  CALL FNIODBL(CRW,NUNIT,
     >                 DUNCORLO(I,J),LPRINT,CHTMP)
                  CHTMP = "  B0        DUnCorUp("
     >                 //CH3TMP//","//CH1TMP4//")"
                  CALL FNIODBL(CRW,NUNIT,
     >                 DUNCORUP(I,J),LPRINT,CHTMP)
C---  WRITE(*,'(A,G10.4)')"ABSUPU: ",DUnCorUp(i,j)*DyVal(i)/100.D0
C---  WRITE(*,'(A,G10.4)')"ABSLOU: ",DUnCorLo(i,j)*DyVal(i)/100.D0
               ENDDO
               DO J=1,NDCORREL
                  WRITE(CH1TMP4,'(I1)'),J
                  CHTMP = "  B0        DCorLo("
     >                 //CH3TMP//","//CH1TMP4//")"
                  CALL FNIODBL(CRW,NUNIT,
     >                 DCORLO(I,J),LPRINT,CHTMP)
                  CHTMP = "  B0        DCorUp("
     >                 //CH3TMP//","//CH1TMP4//")"
                  CALL FNIODBL(CRW,NUNIT,
     >                 DCORUP(I,J),LPRINT,CHTMP)
C---  WRITE(*,'(A,G10.4)')"ABSUPC: ",DCorUp(i,j)*DyVal(i)/100.D0
C---  WRITE(*,'(A,G10.4)')"ABSLOC: ",DCorLo(i,j)*DyVal(i)/100.D0
               ENDDO
            ENDDO
            CHTMP = "  B0    NDErrMatrix"
            CALL FNIOINT(CRW,NUNIT, NDERRMATRIX,LPRINT,CHTMP)
            IF (NDERRMATRIX.NE.0) THEN
               WRITE(*,*)"FX9999RW: ERROR! Data with error matrix not"//
     >              " yet implemented, aborted! NDErrMatrix =",
     >              NDERRMATRIX
               STOP
            ENDIF
            GOTO 100
         ENDIF

*---  Coefficient block
         CHTMP = "  B0    IRef("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, IREF(IC),LPRINT,CHTMP)
         CHTMP = "  B0    IScaleDep("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, ISCALEDEP(IC),LPRINT,CHTMP)
         CHTMP = "  B0    Nevt("//CH1TMP//")"
         CALL FNIOLINT(CRW,NUNIT, NEVT(IC),LPRINT,CHTMP)
         CHTMP = "  B0    Npow("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, NPOW(IC),LPRINT,CHTMP)
         CHTMP = "  B0    NPDF("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, NPDF(IC),LPRINT,CHTMP)
         DO I=1,NPDF(IC)
            WRITE(CH1TMP2,'(I1)'),I
            CHTMP = "  B0      NPDFPDG("//CH1TMP//","//CH1TMP2//")"
            CALL FNIOINT(CRW,NUNIT, NPDFPDG(IC,I),LPRINT,CHTMP)
         ENDDO
         CHTMP = "  B0    NPDFDim("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, NPDFDIM(IC),LPRINT,CHTMP)
         CHTMP = "  B0    NFragFunc("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, NFRAGFUNC(IC),LPRINT,CHTMP)
         DO I=1,NFRAGFUNC(IC)
            WRITE(CH1TMP2,'(I1)'),I
            CHTMP = "  B0     NFFPDG("//CH1TMP//","//CH1TMP2//")"
            CALL FNIOINT(CRW,NUNIT, NFFPDG(IC,I),LPRINT,CHTMP)
         ENDDO
         CHTMP = "  B0    NFFDim("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, NFFDIM(IC),LPRINT,CHTMP)
         CHTMP = "  B0    NSubproc("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, NSUBPROC(IC),LPRINT,CHTMP)
         CHTMP = "  B0    IPDFdef1("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, IPDFDEF(IC,1),LPRINT,CHTMP)
         CHTMP = "  B0    IPDFdef2("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, IPDFDEF(IC,2),LPRINT,CHTMP)
         CHTMP = "  B0    IPDFdef3("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, IPDFDEF(IC,3),LPRINT,CHTMP)

*---  Check consistency between Nsubproc and IPDFdef1,2,
         IF (IPDFDEF(IC,1).EQ.0) THEN ! - no predefined set of PDF coefficients
            CONTINUE
         ELSEIF (IPDFDEF(IC,1).EQ.1) THEN ! --- e+e-
            CONTINUE
         ELSEIF (IPDFDEF(IC,1).EQ.2) THEN ! --- DIS
            CONTINUE
         ELSEIF (IPDFDEF(IC,1).EQ.3) THEN ! --- hh/hhbar
            IF (IPDFDEF(IC,2).EQ.1) THEN 
               IF (IPDFDEF(IC,3).EQ.1) THEN
                  IF (NSUBPROC(IC).NE.6) THEN
                     WRITE(*,*)"FX9999RW: ERROR! IPDFdef(3)=1 "//
     >                    "requires 6 subprocesses for hh. Stopped."
                     STOP
                  ENDIF
               ELSEIF (IPDFDEF(IC,3).EQ.2) THEN 
                  IF (NSUBPROC(IC).NE.7) THEN
                     WRITE(*,*)"FX9999RW: ERROR! IPDFdef(3)=2 "//
     >                    "requires 7 subprocesses for hh. Stopped."
                     STOP
                  ENDIF
               ENDIF
            ELSE
               WRITE(*,*)"FX9999RW: ERROR! Case IPDFdef(2)<>1 "//
     >              "n/a for hh/hhbar. Stopped."
               STOP
            ENDIF
         ENDIF
         
         IF (IPDFDEF(IC,1).EQ.0) THEN ! - no predefined set of PDF coefficients
            WRITE(*,*)"FX9999RW: ERROR! Case IPDFdef(1)=0 "//
     >           "not yet implemented. Stopped."
            STOP
         ENDIF
         
         LPRINT = .FALSE.
         IF (IPRINT.GT.2) THEN
            LPRINT = .TRUE.
            WRITE(*,'(A)')""
            WRITE(*,*)SSEP0
            WRITE(*,*)"* fastNLO Table: Block B ic = ",IC
            WRITE(*,*)"*    X Node Details"
            WRITE(*,*)SSEP0
         ENDIF

         IF (NPDF(IC).GT.0) THEN
            
            DO I=1,NOBSBIN
               WRITE(CH3TMP,'(I3)'),I
               CHTMP = "  BX        Nxtot1("//
     >              CH1TMP//",1,"//CH3TMP//")"
               CALL FNIOINT(CRW,NUNIT, NXTOT(IC,1,I),LPRINT,CHTMP)
               DO J=1,NXTOT(IC,1,I)
                  WRITE(CH2TMP,'(I2)'),J 
                  CHTMP = "  BX        XNode1("//
     >                 CH1TMP//","//CH3TMP//","//CH2TMP//")"
                  CALL FNIODBL(CRW,NUNIT, XNODE1(IC,I,J),LPRINT,CHTMP)
               ENDDO
               IF (CRW.EQ.'read') THEN
                  HXLIM1(IC,I) = -SQRT(-LOG10(XNODE1(IC,I,1))) ! for HH
               ENDIF
            ENDDO 

            IF (NPDFDIM(IC).EQ.2) THEN
               
               DO I=1,NOBSBIN
                  WRITE(CH3TMP,'(I3)'),I
                  CHTMP = "  BX        Nxtot2("//
     >                 CH1TMP//",2,"//CH3TMP//")"
                  CALL FNIOINT(CRW,NUNIT, NXTOT(IC,2,I),LPRINT,CHTMP)
                  DO J=1,NXTOT(IC,2,I)
                     WRITE(CH2TMP,'(I2)'),J 
                     CHTMP = "  BX        XNode2("//
     >                    CH1TMP//","//CH3TMP//","//CH2TMP//")"
                     CALL FNIODBL(CRW,NUNIT, XNODE2(IC,I,J),LPRINT
     >                    ,CHTMP)
                  ENDDO
               ENDDO 

            ENDIF
         ENDIF

         IF (NFRAGFUNC(IC).GT.0) THEN ! - no FFs so far
            WRITE(*,*)"FX9999RW: ERROR! FragFuncs "//
     >           "not yet implemented. Stopped."
            STOP
         ENDIF

         IF (IPRINT.GT.2) THEN
            WRITE(*,*)CSEP0
         ENDIF

         LPRINT = .FALSE.
         IF (IPRINT.GT.1) LPRINT=.TRUE.

         CHTMP = "  B0    NScales("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, NSCALES(IC),LPRINT,CHTMP)
         CHTMP = "  B0    NScaleDim("//CH1TMP//")"
         CALL FNIOINT(CRW,NUNIT, NSCALEDIM(IC),LPRINT,CHTMP)

         DO I=1,NSCALES(IC)
            WRITE(CH1TMP2,'(I1)'),I
            CHTMP = "  B0      NScales("//
     >           CH1TMP//","//CH1TMP2//")"
            CALL FNIOINT(CRW,NUNIT, ISCALE(IC,I),LPRINT,CHTMP)
         ENDDO
         DO I=1,NSCALEDIM(IC)
            WRITE(CH1TMP2,'(I1)'),I
            CHTMP = "  B0      NScaleDescript("//
     >           CH1TMP//","//CH1TMP2//")"
            CALL FNIOINT(CRW,NUNIT, NSCALEDESCRIPT(IC,I),LPRINT,CHTMP)
            DO J=1,NSCALEDESCRIPT(IC,I)
               WRITE(CH1TMP3,'(I1)'),J
               CHTMP = "  B0        ScaleDescript("//
     >              CH1TMP//","//CH1TMP2//","//CH1TMP3//")"
               CALL FNIOCHAR(CRW,NUNIT,
     >              SCALEDESCRIPT(IC,I,J),LPRINT,CHTMP)
            ENDDO
         ENDDO
         DO I=1,NSCALEDIM(IC)
            WRITE(CH1TMP2,'(I1)'),I
            CHTMP = "  B0      NScaleVar("//
     >           CH1TMP//","//CH1TMP2//")"
            CALL FNIOINT(CRW,NUNIT, NSCALEVAR(IC,I),LPRINT,CHTMP)
            CHTMP = "  B0      NScaleNode("//
     >           CH1TMP//","//CH1TMP2//")"
            CALL FNIOINT(CRW,NUNIT, NSCALENODE(IC,I),LPRINT,CHTMP)
         ENDDO

         LPRINT = .FALSE.
         IF (IPRINT.GT.3) THEN
            LPRINT = .TRUE.
            WRITE(*,'(A)')""
            WRITE(*,*)SSEP0
            WRITE(*,*)"* fastNLO Table: Block B ic = ",ic
            WRITE(*,*)"*    Scale Node Details"
            WRITE(*,*)SSEP0
         ENDIF

         DO I=1,NSCALEDIM(IC)
            DO J=1,NSCALEVAR(IC,I)
               WRITE(CH1TMP2,'(I1)'),I
               WRITE(CH1TMP3,'(I1)'),J
               CHTMP = "  BS        ScaleFac("//
     >              CH1TMP//","//CH1TMP2//","//CH1TMP3//")"
               CALL FNIODBL(CRW,NUNIT, SCALEFAC(IC,I,J),LPRINT,CHTMP)
            ENDDO
         ENDDO

         DO I=1,NOBSBIN
            DO J=1,NSCALEDIM(IC)
               DO K=1,NSCALEVAR(IC,J)
                  DO L=1,NSCALENODE(IC,J)
                     WRITE(CH3TMP,'(I3)'),I
                     WRITE(CH1TMP2,'(I1)'),J
                     WRITE(CH1TMP3,'(I1)'),K
                     WRITE(CH1TMP4,'(I1)'),L
                     CHTMP = "  BS          ScaleNode("//
     >                    CH3TMP//","//CH1TMP2//","//
     >                    CH1TMP3//","//CH1TMP4//")"
                     CALL FNIODBL(CRW,NUNIT, SCALENODE(IC,I,J,K,L),
     >                    LPRINT,CHTMP)
                     IF (CRW.EQ.'read') THEN
                        HSCALENODE(IC,I,J,K,L) = 
     >                       LOG(LOG(SCALENODE(IC,I,J,K,L)/0.25D0))
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         
         IF (IPRINT.GT.3) THEN
            WRITE(*,*)CSEP0
         ENDIF

         LPRINT = .FALSE.
         IF (IPRINT.GT.4) THEN
            LPRINT = .TRUE.
            WRITE(*,'(A)')""
            WRITE(*,*)SSEP0
            WRITE(*,*)"* fastNLO Table: Block B ic = ",ic
            WRITE(*,*)"*    Sigma Tilde"
            WRITE(*,*)"* Not implemented yet!"
            WRITE(*,*)SSEP0
            LPRINT = .FALSE.
         ENDIF

         DO I=1,NOBSBIN
            DO K=1,NSCALEVAR(IC,1)
               DO L=1,NSCALENODE(IC,1)
*---  Here we assume NFragFunc=0
                  IF (NFRAGFUNC(IC).GT.0) THEN
                     WRITE(*,*)"FX9999RW: ERROR! NFragFunc>0 "//
     >                    "not yet implemented. Stopped."
                     STOP
                  ENDIF
                  IF (NPDFDIM(IC).EQ.0) THEN
                     NXMAX = NXTOT(IC,1,I)
                  ELSEIF (NPDFDIM(IC).EQ.1) THEN
                     NXMAX = (NXTOT(IC,1,I)**2+NXTOT(IC,1,I))/2
                  ELSEIF (NPDFDIM(IC).EQ.2) THEN 
                     NXMAX = 
     >                    NXTOT(IC,1,I)*NXTOT(IC,2,I)
                     WRITE(*,*)"FX9999RW: ERROR! NPDFdim = 2 "//
     >                    "not yet enabled. Stopped."
                     STOP
                  ELSE
                     WRITE(*,*)"FX9999RW: ERROR! NPDFdim > 2 "//
     >                    "not enabled. Stopped."
                     STOP
                  ENDIF
                  DO M=1,NXMAX
                     DO N=1,NSUBPROC(IC)
                        CALL FNIODBL(CRW,NUNIT,
     >                       SIGMATILDE(IC,I,1,K,L,M,N),LPRINT,CHTMP)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         IF (IPRINT.GT.4) THEN
            WRITE(*,*)CSEP0
         ENDIF

         LPRINT = .FALSE.
         IF (IPRINT.GT.1) LPRINT=.TRUE.
         
 100     CONTINUE               ! - end of block
      ENDDO

      IF (IPRINT.GT.1) THEN
         WRITE(*,*)CSEP0
      ENDIF
      LPRINT = .FALSE.

*---  End of table
      CALL FNIOISEP(CRW,NUNIT,LPRINT," T   ISep")
      CALL FNIOISEP(CRW,NUNIT,LPRINT," F   ISep")
      CLOSE(2)

      RETURN
      END
