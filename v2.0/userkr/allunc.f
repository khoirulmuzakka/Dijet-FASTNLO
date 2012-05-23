      PROGRAM ALLUNC
* ---------------------------------------------------------------------
* K. Rabbertz 07.09.2008 First try to integrate all uncertainties
*                        into one job
*
* ALLUNC - Program to derive all (PDF, statistical, algorithmic,
*          scale ...) uncertainties using fastNLO tables
*
* 10.06.2011 kr: Replace CERNLIB function LENOCC by f95 Standard LEN_TRIM
* 31.03.2011 kr: Some modifications to improve v2 compatibility
* ---------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      INCLUDE "uncert.inc"
      CHARACTER*255 SCENARIO,TABPATH,TABNAME,REFNAME
      CHARACTER*255 FILENAME,FILENAMES,HISTFILE
      CHARACTER*255 TABNAMN,FILENAMN
      CHARACTER*255 PDFNAM,PDFSET,PDFSET2,PDFPATH,LHAPDF,CHTMP
      CHARACTER*255 FILEBASE,LOFILE,NLOFILE
      CHARACTER*255 BORNNAME,NLONAME
      CHARACTER*8 CH8TMP
      CHARACTER*4 CH4TMP
      CHARACTER*4 NO
      INTEGER BORNN,NLON
      INTEGER I,J,MYPDF,IPDF,IPDFUD,MYPDFAS,IOPDF,IOAS
      INTEGER ITAB,NTAB,NFOUND,NFAIL
      INTEGER ISTAT,ISCL,IORD,IORD2,IBIN,NBIN,ISUB,IRAP,IPT,NTMP
      INTEGER IHIST,IPHASE,ISTEP,IETYPE
      LOGICAL LONE,LPDF,LSTAT,LSER,LSCL,LRAT,LALG,LNRM,LTAB
cnew
      INTEGER IPRINT
      LOGICAL LLO,LNLO,LTHC1L,LTHC2L,LNPC1,LDATA
      DOUBLE PRECISION ALPS,FNALPHAS,ALPHASPDF,ASMZPDF,ASUP,ASDN
      DOUBLE PRECISION XMUR,XMUF,QLAM4,QLAM5,BWGT
      DOUBLE PRECISION DSTMP(4)
      DOUBLE PRECISION WTXTMP(MXOBSBIN,MXSUBPROC+1,NMAXORD+1)
      DOUBLE PRECISION WTXLTMP(MXOBSBIN,MXSUBPROC+1,NMAXORD+1)
      DOUBLE PRECISION WTXUTMP(MXOBSBIN,MXSUBPROC+1,NMAXORD+1)
c - To unify quoted uncertainties (CL68,CL90,special)
c - Convert from CL68 to CL90 values
c - TOCL90 = 1.64485D0 ! SQRT(2.D0)/InvERF(0.9D0)
c - Convert from GJR to CTEQ (CL90) values ?
c - TOCL90GJR = 2.12766D0! 1.D0/0.47D0
      DOUBLE PRECISION TOCL90,TOCL90GJR
      PARAMETER (TOCL90 = 1.64485D0, TOCL90GJR = 2.12766D0)
      
*---  Define series of scale factor settings to test. Last and 8th entry
*---  is (0,0) and is not to be used!
      INTEGER MXSCALECOMB
      PARAMETER (MXSCALECOMB=2*MXSCALEVAR)
      INTEGER ISCLPT(MXSCALECOMB)
      DOUBLE PRECISION XMURS(MXSCALECOMB),XMUFS(MXSCALECOMB)
      DATA XMURS/1.0D0,0.5D0,2.0D0,0.5D0,1.0D0,1.D0,2.D0,0.0D0/
      DATA XMUFS/1.0D0,0.5D0,2.0D0,1.0D0,0.5D0,2.D0,1.D0,0.0D0/

c - Attention!!! This must be declared consistent with the
c                definition in the commonblock!!!!!
      DOUBLE PRECISION XSECT0(MXOBSBIN,3)
      REAL PT(NPTMAX)
      INTEGER IMODE,IWEIGHT,NRAP

ckr Old Z mass
ckr      PARAMETER (ZMASS = 91.187D0)
ckr Z mass from PDG 2006
      DOUBLE PRECISION ZMASS
      PARAMETER (ZMASS = 91.1876D0)

      CHARACTER*255 ASMODE,ASMODETMP
      DOUBLE PRECISION ASMZVAL,DASMZVAL,ASMZTMP
      INTEGER IASLOOP
      COMMON/STEER/ASMZVAL,IASLOOP,ASMODE

c --- Parse command line
      WRITE(*,*)"\n ########################################"//
     >     "################################"
      WRITE(*,*)"# ALLUNC"
      WRITE(*,*)"########################################"//
     >     "################################"
      WRITE(*,*)"# Program to derive all (PDF, statistical, "//
     >     "algorithmic, scale ...)"
      WRITE(*,*)"# uncertainties using fastNLO tables"
      WRITE(*,*)"########################################"//
     >     "################################"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      WRITE(*,*)"ALLUNC: Program Steering"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"

*---Scenario
      LRAT = .FALSE.
      LNRM = .FALSE.
      LTAB = .FALSE.
      LSCL = .FALSE.
      IF (IARGC().LT.1) THEN
         SCENARIO = "fnt2003"
         WRITE(*,*)
     >        "ALLUNC: WARNING! No scenario name given, "//
     >        "taking the default fnt2003 instead!"
         WRITE(*,*)"      For an explanation of command line "//
     >        "arguments type:"
         WRITE(*,*)"      ./allunc -h"
      ELSE
         CALL GETARG(1,SCENARIO)
         IF (SCENARIO(1:LEN_TRIM(SCENARIO)).EQ."-h") THEN
            WRITE(*,*)' '
            WRITE(*,*)'Usage: ./allunc [arguments]'
            WRITE(*,*)'  Scenario name, def. = fnt2003'
            WRITE(*,*)'  Table path, def. = .'
            WRITE(*,*)'     Table names have to be of style:'
c --- Use '...' with \", otherwise gfortran complains 
            WRITE(*,*)'     \"scenario\".tab'
            WRITE(*,*)'     \"scenario\"ref.tab'
            WRITE(*,*)'     The statistics tables named e.g.:'
            WRITE(*,*)'     \"scenario\"-hhc-nlo-2jet_nnnn.tab'
            WRITE(*,*)'     have to be in the subdir. \"stat\".'
            WRITE(*,*)'  HBOOK output file, def. = scenario.hbk'
            WRITE(*,*)'  Derive algorithmic uncertainty, def. = no'
            WRITE(*,*)'  Last LO stat. table number, def. = -1'
            WRITE(*,*)'  Last NLO stat. table number, def. = -1'
            WRITE(*,*)'  Number of pre-defined mu_r, mu_f scale '//
     >           'settings to investigate, def. = 1'
            WRITE(*,*)'  PDF set, def. = cteq66.LHgrid'
            WRITE(*,*)'  PDF path, def. = $(LHAPDF)/'//
     >           '../share/lhapdf/PDFsets'
            WRITE(*,*)'  alpha_s calc., def. = PDF (from PDF set)'
            WRITE(*,*)'    alt. = PY: 0-, 1- and 2-loop '//
     >           '(from Pythia 6.4 using Lambda_4 from PDF)'
            WRITE(*,*)'    alt. = KR: 1-, 2- and 3-loop '//
     >           '(from hep-ph/9506442)'
            WRITE(*,*)'    alt. = MW: 2- and 4-loop '//
     >           '(from hep-ph/9806404)'
            WRITE(*,*)'  alpha_s(M_Z), def. from PDF set'
            WRITE(*,*)'     (in mode PY this has to be Lambda_4/GeV!)'
            WRITE(*,*)'  alpha_s(M_Z) variation, def. = 0.D0, i.e. none'
            WRITE(*,*)'  alpha_s loop order, def. from PDF set'
            WRITE(*,*)'  PDF uncertainties in '//
     >           'eigen vector (1; CTEQ,MSTW,ABKM), '//
     >           'toy MC (2; NNPDF) or mixed method (3; HERAPDF), '//
     >           'def. = 1'
            WRITE(*,*)' '
            STOP
ckr To be cross-checked for each new scenario
         ELSEIF (SCENARIO(1:8).EQ."fnl2332c") THEN
            LNRM = .TRUE.
            LTAB = .FALSE.
            WRITE(*,*)
     >           "ALLUNC: Deriving x section ratios"
         ELSEIF (SCENARIO(1:7).EQ."fnl2442") THEN
ckr Original version for fnl2442: Works fine, trivial division in rap 3
            LRAT = .FALSE.
ckr New norm. version for fnl2442: Works fine, trivial division in rap 4
            LNRM = .TRUE.
            LTAB = .FALSE.
            WRITE(*,*)
     >           "ALLUNC: Deriving x section ratios"
         ELSEIF (SCENARIO(1:11).EQ."fnl2522diff") THEN
            LNRM = .TRUE.
            LTAB = .TRUE.
            WRITE(*,*)
     >           "ALLUNC: Deriving normalized distributions"
         ELSEIF (SCENARIO(1:7).EQ."fnl2622e".OR.
     >           SCENARIO(1:7).EQ."fnl2652") THEN
            LNRM = .TRUE.
            LTAB = .FALSE.
            WRITE(*,*)
     >           "ALLUNC: Deriving normalized distributions"
         ELSEIF (SCENARIO(1:10).EQ."fnl2722num".OR.
     >           SCENARIO(1:10).EQ."fnl2732num".OR.
     >           SCENARIO(1:10).EQ."fnl2742num") THEN
            LNRM = .TRUE.
            LTAB = .TRUE.
            WRITE(*,*)
     >           "ALLUNC: Deriving normalized distributions"
         ENDIF
         WRITE(*,*)"ALLUNC: Evaluating scenario: ",
     >        SCENARIO(1:LEN_TRIM(SCENARIO))
      ENDIF
ckr      LNRM = .FALSE.
      TABNAME = SCENARIO(1:LEN_TRIM(SCENARIO))//".tab"
      IF (LNRM.AND.LTAB) THEN
         TABNAMN = SCENARIO(1:7)//"norm"//".tab"
      ENDIF
      REFNAME = SCENARIO(1:LEN_TRIM(SCENARIO))//"ref.tab"

*---Path to tables
      TABPATH = "X"
      IF (IARGC().GE.2) THEN
         CALL GETARG(2,TABPATH)
      ENDIF
      IF (IARGC().LT.2.OR.TABPATH(1:1).EQ."_") THEN
         TABPATH = "."
         WRITE(*,*)
     >        "ALLUNC: WARNING! No table path given, "//
     >        "taking . instead!"
      ENDIF
      WRITE(*,*)"ALLUNC: Using table path: ",
     >     TABPATH(1:LEN_TRIM(TABPATH))
      FILENAME = TABPATH(1:LEN_TRIM(TABPATH))//"/"//
     >     TABNAME(1:LEN_TRIM(TABNAME))
      WRITE(*,*)"ALLUNC: Taking primary table ",
     >     FILENAME(1:LEN_TRIM(FILENAME))
      IF (LNRM.AND.LTAB) THEN
         FILENAMN = TABPATH(1:LEN_TRIM(TABPATH))//"/"//
     >        TABNAMN(1:LEN_TRIM(TABNAMN))
         WRITE(*,*)"ALLUNC: Taking normalization table ",
     >        FILENAMN(1:LEN_TRIM(FILENAMN))
      ENDIF

*---HBOOK filename
      HISTFILE = "X"
      IF (IARGC().GE.3) THEN
         CALL GETARG(3,HISTFILE)
      ENDIF
      IF (IARGC().LT.3.OR.HISTFILE(1:1).EQ."_") THEN
         HISTFILE = SCENARIO(1:LEN_TRIM(SCENARIO))//".hbk"
         WRITE(*,*)
     >        "ALLUNC: WARNING! No output filename given, "//
     >        "taking scenario.hbk instead!"
      ELSE
         WRITE(*,*)"ALLUNC: Creating output file: ",
     >        HISTFILE(1:LEN_TRIM(HISTFILE))
      ENDIF

*---Derive algorithmic uncertainty?
      CH4TMP = "X"
      IF (IARGC().GE.4) THEN
         CALL GETARG(4,CH4TMP)
      ENDIF
      IF (IARGC().LT.4.OR.CH4TMP(1:1).EQ."_".OR.
     >     CH4TMP(1:2).EQ."no") THEN
         LALG   =  .FALSE.
         WRITE(*,*)
     >        "ALLUNC: Do NOT derive algorithmic uncertainty!"
      ELSE
         IF (CH4TMP(1:3).EQ."yes") THEN
            LALG = .TRUE.
            WRITE(*,*)"ALLUNC: Derive algorithmic uncertainty."
         ELSE
            LALG = .FALSE.
            WRITE(*,*)
     >           "ALLUNC: WARNING! "//
     >           "Do NOT derive algorithmic uncertainty!"
         ENDIF
      ENDIF

*---Last LO table to use 
      CH4TMP = "X"
      IF (IARGC().GE.5) THEN
         CALL GETARG(5,CH4TMP)
      ENDIF
      IF (IARGC().LT.5.OR.CH4TMP(1:1).EQ."_") THEN
         CH4TMP = "-1"
         BORNN  =  -1
         WRITE(*,*)
     >        "ALLUNC: WARNING! Last number of LO tables not given, "//
     >        "using -1 instead ==> no stat. uncertainty!"
      ELSE
         READ(CH4TMP,'(I4)'),BORNN
         WRITE(*,*)"ALLUNC: Last LO table number: ",BORNN
      ENDIF

*---Last NLO table to use 
      CH4TMP = "X"
      IF (IARGC().GE.6) THEN
         CALL GETARG(6,CH4TMP)
      ENDIF
      IF (IARGC().LT.6.OR.CH4TMP(1:1).EQ."_") THEN
         CH4TMP = "-1"
         NLON   =  -1
         WRITE(*,*)
     >        "ALLUNC: WARNING! Last number of NLO tables not given, "//
     >        "using -1 instead ==> no stat. uncertainty!"
      ELSE
         READ(CH4TMP,'(I4)'),NLON
         WRITE(*,*)"ALLUNC: Last NLO table number: ",NLON
      ENDIF

*---No. of pre-defined scale settings to investigate
      CH4TMP = "X"
      IF (IARGC().GE.7) THEN
         CALL GETARG(7,CH4TMP)
      ENDIF
      IF (IARGC().LT.7.OR.CH4TMP(1:1).EQ."_") THEN
         CH4TMP = "1"
         NSCLS   = 1
         WRITE(*,*)
     >        "ALLUNC: WARNING! No. of scale settings not given, "//
     >        "using 1 instead ==> no scale uncertainties!"
      ELSE
         READ(CH4TMP,'(I1)'),NSCLS
         IF (NSCLS.LT.1) THEN
            WRITE(*,*)
     >           "ALLUNC: ERROR! No scale setting "//
     >           "or even less??? Aborting! NSCLS = ",
     >           NSCLS
            STOP
         ELSEIF (NSCLS.GT.MXSCALECOMB-1) THEN
            WRITE(*,*)
     >           "ALLUNC: ERROR! Too many scale settings "//
     >           "requested, aborting! NSCLS = ",
     >           NSCLS
            STOP
         ELSE
            IF (NSCLS.GE.3) LSCL = .TRUE.
            WRITE(*,*)"ALLUNC: No. of scale settings: ",NSCLS
         ENDIF
      ENDIF
      
*---PDF set
      PDFSET = "X"
      IF (IARGC().GE.8) THEN
         CALL GETARG(8,PDFSET)
      ENDIF
      IF (IARGC().LT.8.OR.PDFSET(1:1).EQ."_") THEN
         PDFSET = "cteq66.LHgrid"
         WRITE(*,*)
     >        "ALLUNC: WARNING! No PDF set given, "//
     >        "taking cteq66.LHgrid instead!"
      ELSE
         WRITE(*,*)"ALLUNC: Using PDF set: ",
     >        PDFSET(1:LEN_TRIM(PDFSET))
      ENDIF
      PDFNAM = PDFSET(1:LEN_TRIM(PDFSET))

*---Path to PDF sets
      CHTMP = "X"
      IF (IARGC().GE.9) THEN
         CALL GETARG(9,CHTMP)
      ENDIF
      IF (IARGC().LT.9.OR.CHTMP(1:1).EQ."_") THEN
         PDFPATH = "/../share/lhapdf/PDFsets"
         WRITE(*,*)
     >        "ALLUNC: No PDF path given, "//
     >        "assuming: $(LHAPDF)"//PDFPATH(1:LEN_TRIM(PDFPATH))
*---  Initialize path to LHAPDF libs
         CALL GETENV("LHAPDF",LHAPDF)
         IF (LEN_TRIM(LHAPDF).EQ.0) THEN
            WRITE(*,*)"\nALLUNC: ERROR! $LHAPDF not set, aborting!"
            STOP
         ENDIF
         PDFPATH = LHAPDF(1:LEN_TRIM(LHAPDF))//
     >        PDFPATH(1:LEN_TRIM(PDFPATH))
      ELSE
         PDFPATH = CHTMP(1:LEN_TRIM(CHTMP)) 
      ENDIF
      WRITE(*,*)"ALLUNC: Looking for LHAPDF PDF sets in path: "//
     >     PDFPATH(1:LEN_TRIM(PDFPATH))
      PDFSET = PDFPATH(1:LEN_TRIM(PDFPATH))//"/"//PDFSET
      WRITE(*,*)"ALLUNC: Taking PDF set "
     >     //PDFSET(1:LEN_TRIM(PDFSET))

*---alpha_s mode
      CHTMP = "X"
      IF (IARGC().GE.10) THEN
         CALL GETARG(10,CHTMP)
      ENDIF
      IF (IARGC().LT.10.OR.CHTMP(1:1).EQ."_") THEN
         ASMODE = "PDF"
         WRITE(*,*)
     >        "ALLUNC: No alpha_s mode given, "//
     >        "using alpha_s according to PDF set"
      ELSE
         ASMODE = CHTMP(1:LEN_TRIM(CHTMP))
         WRITE(*,*)"ALLUNC: Using alpha_s mode: ",
     >        ASMODE(1:LEN_TRIM(ASMODE))
      ENDIF

*---alpha_s(M_Z)
      CH8TMP = "X"
      IF (IARGC().GE.11) THEN
         CALL GETARG(11,CH8TMP)
      ENDIF
      IF (IARGC().LT.11.OR.CH8TMP(1:1).EQ."_") THEN
         ASMZVAL = -1D0
         WRITE(*,*)
     >        "ALLUNC: No alpha_s(M_Z) value given, "//
     >        "using alpha_s according to PDF set"
      ELSE
         READ(CH8TMP,'(F9.6)'),ASMZVAL
         WRITE(*,*)"ALLUNC: Using alpha_s(M_Z):",ASMZVAL
      ENDIF

*---alpha_s(M_Z) +- variation
      CH8TMP = "X"
      IF (IARGC().GE.12) THEN
         CALL GETARG(12,CH8TMP)
      ENDIF
      IF (IARGC().LT.12.OR.CH8TMP(1:1).EQ."_") THEN
         DASMZVAL = 0.D0
         WRITE(*,*)
     >        "ALLUNC: No alpha_s(M_Z) variation given, "//
     >        "using central value only"
      ELSE
         READ(CH8TMP,'(F9.6)'),DASMZVAL
         WRITE(*,*)"ALLUNC: Varying PDF alpha_s set member or "//
     >        "alpha_s(M_Z) value by +-:",DASMZVAL
      ENDIF

*---alpha_s loop order in evolution
      CH4TMP = "X"
      IF (IARGC().GE.13) THEN
         CALL GETARG(13,CH4TMP)
      ENDIF
      IF (IARGC().LT.13.OR.CH4TMP(1:1).EQ."_") THEN
         IASLOOP = -1
         WRITE(*,*)
     >        "ALLUNC: No alpha_s loop order given, "//
     >        "using alpha_s loop order according to PDF set"
      ELSE
         READ(CH4TMP,'(I4)'),IASLOOP
         WRITE(*,*)"ALLUNC: Using alpha_s loop order:",IASLOOP
      ENDIF

*---Use eigen vector (CTEQ/MSTW) or toy MC (NNPDF) method for
*---PDF uncertainties
      CH4TMP = "X"
      IF (IARGC().GE.14) THEN
         CALL GETARG(14,CH4TMP)
      ENDIF
      IF (IARGC().LT.14.OR.CH4TMP(1:1).EQ."_".OR.
     >     CH4TMP(1:1).EQ."1") THEN
         IETYPE = 1
         WRITE(*,*)
     >        "ALLUNC: Use eigen vector (CTEQ/MSTW/ABKM) method."
      ELSE
         READ(CH4TMP,'(I4)'),IETYPE
         IF (IETYPE.EQ.1) THEN
            WRITE(*,*)
     >           "ALLUNC: Use eigen vector (CTEQ/MSTW/ABKM) method."
         ELSEIF (IETYPE.EQ.2) THEN
            WRITE(*,*)"ALLUNC: Use toy MC (NNPDF) method."
         ELSEIF (IETYPE.EQ.3) THEN
            WRITE(*,*)"ALLUNC: Use mixed (HERAPDF) method."
         ELSE
            WRITE(*,*)
     >           "ALLUNC: ERROR! Undefined PDF error method, aborted!"//
     >           " IETYPE = ",IETYPE
         ENDIF
      ENDIF

*---Too many arguments
      IF (IARGC().GT.15) THEN
         WRITE(*,*)"\nALLUNC: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF

*---  Initialize LHAPDF, no PDF set printout after first call
      CALL INITPDFSET(PDFSET(1:LEN_TRIM(PDFSET)))
      CALL SETLHAPARM('SILENT')
      
*---  Initialize one member, 0=best fit member
      CALL INITPDF(0)
      
*---  Write out some info on best fit member      
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      WRITE(*,*)"ALLUNC: PDF Set Info"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      CALL NUMBERPDF(MYPDF)
      CALL GETORDERPDF(IOPDF)
      CALL GETORDERAS(IOAS)
      CALL GETLAM4(0,QLAM4)
      CALL GETLAM5(0,QLAM5)
      WRITE(*,*) "ALLUNC: The PDF set has",MYPDF+1," members"
      WRITE(*,*) "ALLUNC: The PDF is of order",IOPDF+1
      WRITE(*,*) "ALLUNC: alpha_s was used in",IOAS+1,
     >     "-loop order in the PDF"
      WRITE(*,*) "ALLUNC: The lambda_4 value for member 0 is",QLAM4
      WRITE(*,*) "ALLUNC: The lambda_5 value for member 0 is",QLAM5
      
*---  Check primary table existence
      FILENAME = TABPATH(1:LEN_TRIM(TABPATH))//"/"//TABNAME
      WRITE(*,*)"ALLUNC: Checking primary table: "//
     >     FILENAME(1:LEN_TRIM(FILENAME))
      OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=ISTAT)
      IF (ISTAT.NE.0) THEN
         WRITE(*,*)"ALLUNC: ERROR! Primary table not found, "//
     >        "aborting! IOSTAT = ",ISTAT
         STOP
      ELSE
         CLOSE(2)
      ENDIF

*---  Check uncertainties to derive 
      LONE  = MYPDF.LE.1
      LSTAT = BORNN.GE.2.OR.NLON.GE.2
      IF (LALG) THEN
         FILENAME = TABPATH(1:LEN_TRIM(TABPATH))//"/"//REFNAME
         WRITE(*,*)"ALLUNC: Checking reference table: "//
     >        FILENAME(1:LEN_TRIM(FILENAME))
         OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=ISTAT)
         IF (ISTAT.NE.0) THEN
            WRITE(*,*)"ALLUNC: WARNING! Reference table not found, "//
     >           "skipped! IOSTAT = ",ISTAT
            LALG = .FALSE.
         ELSE
            CLOSE(2)
         ENDIF
      ENDIF
ckr To be implemented
ckr      LSER  = .NOT.LONE.AND.MYPDF.LT.10.AND..NOT.LSTAT.AND..NOT.LALG
      LSER  = .FALSE.
      IF (DASMZVAL.GT.0.D0) LSER = .TRUE.
      LPDF  = .NOT.LONE.AND..NOT.LSER

      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      WRITE(*,*)"ALLUNC: Uncertainties"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      IF (LONE) THEN
         WRITE(*,*)"ALLUNC: Only central PDF available."
      ENDIF
      IF (LPDF) THEN
         WRITE(*,*)"ALLUNC: Deriving PDF uncertainties"
      ENDIF
      IF (LSTAT) THEN
         WRITE(*,*)"ALLUNC: Deriving statistical uncertainties"
      ENDIF
      IF (LSER) THEN
         WRITE(*,*)"ALLUNC: Deriving uncertainty due to a "//
     >        "variation of alpha_s(M_Z)"
         IF (LSTAT.OR.LALG) THEN
            WRITE(*,*)"ALLUNC: WARNING! Uncertainties for"//
     >           " alpha_s(M_Z) cannot be done at the"//
     >           " same time as statistical or algorithmic ones,"//
     >           " switched off !"
            LSER = .FALSE.
         ENDIF
      ENDIF
      IF (LSCL) THEN
         WRITE(*,*)"ALLUNC: Deriving scale uncertainties"
      ENDIF
      IF (LALG) THEN
         WRITE(*,*)"ALLUNC: Deriving algorithmic uncertainties"
         IF (LRAT.OR.LNRM) THEN
            WRITE(*,*)"ALLUNC: WARNING! Uncertainties for ratios or"//
     >           " normalized distributions cannot be done at the"//
     >           " same time as algorithmic ones, switched off"//
     >           " both!"
            LRAT = .FALSE.
            LNRM = .FALSE.
         ENDIF
      ENDIF
      IF (LRAT) THEN
         WRITE(*,*)"ALLUNC: Deriving uncertainties for ratios"
         IF (LNRM) THEN
            WRITE(*,*)"ALLUNC: WARNING! Uncertainties for ratios and"
     >           //" normalized distributions cannot be done at the "
     >           //"same time, switched off normalization!"
            LNRM = .FALSE.
         ENDIF
      ENDIF
      IF (LNRM) THEN
         WRITE(*,*)"ALLUNC: Deriving uncertainties for normalized "//
     >        "distributions"
      ENDIF
      
*---  One initial call - to fill commonblock -> for histo-booking
*---  Use primary table for this (recall: ref. table has 2 x rap. bins)
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      WRITE(*,*)"ALLUNC: Initialize Table, Book Histograms"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      FILENAME = TABPATH(1:LEN_TRIM(TABPATH))//"/"//TABNAME
      
*---  Initialize table
      Call FX9999IN(FILENAME)
      
*---  Initial call to alpha_s interface
      ALPS = FNALPHAS(91.1876D0)
      
*---  Print out contribution list
      Call FX9999CL

*---  Print out scenario information
      Call FX9999NF

*---  Check on existence of LO, NLO, THC, NPC and DATA contributions
      LLO    = .FALSE.
      LNLO   = .FALSE.
      LTHC1L = .FALSE.
      LTHC2L = .FALSE.
      LNPC1  = .FALSE.
      LDATA  = .FALSE.
      DO I=1,NCONTRIB
         IF (ICONTRFLAG1(I).EQ.1.AND.ICONTRFLAG2(I).EQ.1)LLO    = .TRUE.
         IF (ICONTRFLAG1(I).EQ.1.AND.ICONTRFLAG2(I).EQ.2)LNLO   = .TRUE.
         IF (ICONTRFLAG1(I).EQ.2.AND.ICONTRFLAG2(I).EQ.1)LTHC1L = .TRUE.
         IF (ICONTRFLAG1(I).EQ.2.AND.ICONTRFLAG2(I).EQ.2)LTHC2L = .TRUE.
         IF (ICONTRFLAG1(I).EQ.4.AND.ICONTRFLAG2(I).EQ.1.AND.
     >        IADDMULTFLAG(I).EQ.1)
     >        LNPC1 = .TRUE.
         IF (ICONTRFLAG1(I).EQ.0.AND.ICONTRFLAG2(I).EQ.0.AND.
     >        IDATAFLAG(I).EQ.1)
     >        LDATA = .TRUE.
      ENDDO

*---  Determine dimensional subdivisions (NDIM=1,2 only!)
      IF (NDIM.GT.2) THEN
         WRITE(*,*)"ALLUNC: ERROR! Cannot deal yet with "//
     >        "more than 2 dimensions. NDIM = ",NDIM
         STOP
      ENDIF
      DO I=1,NOBSBIN
*---  Very first bin; initialize
         IF (I.EQ.1) THEN
            NTMP = 1
            NBINCOUNTER(1) = 0
            NBINCOUNTER(2) = 1
*---  Bin of 2nd dimension identical 
         ELSEIF ((LOBIN(I-1,2).EQ.LOBIN(I,2)).AND.
     >           (UPBIN(I-1,2).EQ.UPBIN(I,2))) THEN
            NTMP = NTMP + 1
*---  Bin of 2nd dimension changed 
         ELSE
            NDIVCOUNTER(NBINCOUNTER(2)) = NTMP
            NBINCOUNTER(1) = NBINCOUNTER(1) + NTMP
            NBINCOUNTER(2) = NBINCOUNTER(2) + 1
            NTMP = 1
         ENDIF
      ENDDO
      NDIVCOUNTER(NBINCOUNTER(2)) = NTMP
      NBINCOUNTER(1) = NBINCOUNTER(1) + NTMP
      IBIN = 0
      DO J=1,NBINCOUNTER(2)
         DO I=1,NDIVCOUNTER(J)
            IBIN = IBIN + 1
            IOBSPOINTER(I,J) = IBIN
         ENDDO
      ENDDO            

*---  Define v14 variables
      NRAPIDITY = NBINCOUNTER(2)
      DO I=1,NRAPIDITY
         NPT(I) = NDIVCOUNTER(I)
         DO J=1,NPT(I)
            PTBIN(I,J) = LOBIN(IOBSPOINTER(J,I),1)
         ENDDO
         PTBIN(I,NPT(I)+1) = UPBIN(IOBSPOINTER(NPT(I),I),1)
      ENDDO
C---  Set no. of subprocesses (old: 7) to the one for NLO contribution 2
C---  (new: 6 for 2-3 partons, 7 for 3-4 partons)
      NSBPRC = NSUBPROC(2)
C---  Set NORD to 2 for LO & NLO
      NORD = 2
C---  Contribution 2 (NLO), only 1 scale dimension
ckr      NSCLS = NSCALEVAR(2,1)

*---  Define output formats
ckr 900     FORMAT(1P,I5,3(3X,E21.14))
 900  FORMAT(1P,I5,3(6X,E18.11))
 901  FORMAT(1P,I5,2(6X,E18.11))
 902  FORMAT(3I6,3E16.5,5(F10.3,3X))
      
*---  Initial settings
      CALL FNSET("P_RESET",0)   ! Reset all selections to zero
      IF (LLO) THEN
         CALL FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
      ELSE 
         WRITE(*,*)"ALLUNC: ERROR! No LO found, stopped."
         STOP
      ENDIF
      IF (LNLO) THEN
         CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
      ELSE 
         WRITE(*,*)"ALLUNC: ERROR! No NLO found, stopped."
         STOP
      ENDIF

*---  Look for pointers matching requested scale settings
      IPRINT = 0
      DO ISCL=1,NSCLS
         XMUR = XMURS(ISCL)
         XMUF = XMUFS(ISCL)
         CALL FX9999PT(XMUR,XMUF,IPRINT)
         IF (ISCALEPOINTER(INLO).GT.0) THEN
            ISCLPT(ISCL) = ISCALEPOINTER(INLO)
         ELSE 
            WRITE(*,*)"ALLUNC: ERROR! Required scale no. ",ISCL,
     >           "not found, stopped."
            STOP
         ENDIF
      ENDDO

*---  Book histograms
      CALL PDFHIST(1,HISTFILE,LONE,LPDF,LSTAT,LALG,LSER,MYPDF,
     >     LRAT.OR.LNRM,LSCL,ISCLPT(1))
      WRITE(*,*)"ALLUNC: The observable has",NOBSBIN," bins -",
     >     NSBPRC," subprocesses"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      


*---  PDF part
*     Use primary table
*     New call: a single call for each scale
*     1st argument:  name of table
*     2nd argument:  xmur  prefactor for nominal ren-scale
*     -              any choice is possible, but please note 
*     -              that NNLO-NLL works only for xmur=xmuf
*     3rd argument:  xmuf  prefactor for nominal fact-scale
*     -              only a few choices are possible
*                    (see output or table documentation)
*     4th argument:  0: no ascii output       1: print results
*     5th argument:  array to return results

*---  Compute central result (LO+NLO) and PDF uncertainties
      IF (LONE.OR.LPDF) THEN
         WRITE(*,*)"****************************************"//
     >        "********************************"
         WRITE(*,*)"ALLUNC: Evaluating central result"
         IF (LPDF) WRITE(*,*)"ALLUNC: Evaluating PDF uncertainties"
         WRITE(*,*)"****************************************"//
     >        "********************************"
         WRITE(*,*)"========================================"//
     >        "================================"
         IF (LONE) THEN
            WRITE(*,*)"No PDF Uncertainties possible"
         ELSEIF (LPDF) THEN
            WRITE(*,*)"Relative PDF Uncertainties"
            WRITE(*,*)"- the printed values are for the total "//
     >           "cross section summed over all subprocesses"
            WRITE(*,*)"- histograms contain more detailed results"
         ENDIF
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         WRITE(*,*)" bin       cross section           "//
     >        "lower PDF uncertainty   upper PDF uncertainty"

*---  Only primary scale
         DO I=1,1
            ISCL = ISCLPT(I)
            XMUR = XMURS(I)
            XMUF = XMUFS(I)
*---  Part 1: Basic PDF uncertainty using a single PDF set (and only one
*---         initialization call!). For PDF uncertainty calculation from
*---         multiple PDF sets like full HERAPDF make sure to switch back
*---         to original one here
            IF (IETYPE.EQ.3) THEN
               CALL INITPDFSET(PDFSET(1:LEN_TRIM(PDFSET)))
            ENDIF
*---  NLO contribution 2 (INLO), scale dimension 1
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            WRITE(*,*)"ALLUNC: Now scale no.",i,"; mur, muf = ",
     >           xmur,xmuf
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            CALL INITPDF(0)
            IPHASE  = 1
            IMODE   = 1
            IWEIGHT = 0
            NRAP  = NRAPIDITY
            IF (IETYPE.EQ.2) THEN
               IMODE = 2
            ENDIF
            IF (LRAT.OR.LNRM) THEN
               NRAP = 2*NRAPIDITY
            ENDIF
            CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
            ISTEP = 0
C---  WRITE(*,*)"AAAAA: ALLUNC STEP = ",ISTEP
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
            IF (LNRM) THEN
*--- Load normalization table with potentially different binning!
               IF (LTAB) THEN
                  CALL FX9999IN(FILENAMN)
                  CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
               ENDIF
               ISTEP = 1
C---  WRITE(*,*)"BBBBB: ALLUNC STEP = ",ISTEP
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))
               IF (LTAB) THEN
                  CALL FX9999IN(FILENAME)
                  CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
               ENDIF
            ENDIF
            ISTEP = 2
C---  WRITE(*,*)"CCCCC: ALLUNC STEP = ",ISTEP
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
            
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
            IPHASE = 2

*--- Do loop runs once even if MYPDF=0! => Avoid with IF statement
            IF (LPDF) THEN
               DO J=1,MYPDF
ckr               DO J=1,3
                  CALL INITPDF(J)
                  CALL FX9999IN(FILENAME)
                  CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                  IF (LNRM) THEN
                     ISTEP = 3
C---  WRITE(*,*)"DDDDD: ALLUNC STEP = ",ISTEP
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
*--- Load normalization table with potentially different binning!
                     IF (LTAB) THEN
                        CALL FX9999IN(FILENAMN)
                        CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                     ENDIF
                     ISTEP = 4
ckr                     WRITE(*,*)"EEEEE: ALLUNC STEP = ",ISTEP
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                     IF (LTAB) THEN
                        CALL FX9999IN(FILENAME)
                        CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                     ENDIF
                     ISTEP = 5
ckr                     WRITE(*,*)"FFFFF: ALLUNC STEP = ",ISTEP
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                  ENDIF
                  CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
               ENDDO
            ENDIF
            
            IPHASE = 3
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

*--- Part 2: Additional PDF uncertainty parts using a separate PDF set
            IF (IETYPE.EQ.3) THEN
               IF (PDFNAM(1:LEN_TRIM(PDFNAM)).EQ.
     >              "HERAPDF10_EIG.LHgrid") THEN
                  PDFSET2 = PDFPATH(1:LEN_TRIM(PDFPATH))//
     >                 "/HERAPDF10_VAR.LHgrid"
               ELSEIF (PDFNAM(1:LEN_TRIM(PDFNAM)).EQ.
     >                 "HERAPDF15_EIG.LHgrid") THEN
                  PDFSET2 = PDFPATH(1:LEN_TRIM(PDFPATH))//
     >                 "/HERAPDF15_VAR.LHgrid"
ckr                  WRITE(*,*)"ALLUNC: Taking second HERAPDF set: "//
ckr     >                 PDFSET2(1:LEN_TRIM(PDFSET2))
               ELSE
                  WRITE(*,*)"ALLUNC: Illegal HERAPDF set, aborted! "//
     >                 "PDFNAM: ",PDFNAM(1:LEN_TRIM(PDFNAM))
                  STOP
               ENDIF
*--- Back up central result and relevant uncertainty
               IBIN = 0
               DO IRAP=1,NRAP
                  DO IPT=1,NPT(IRAP)
                     IBIN = IBIN+1
                     DO IORD=1,NORD+1
                        DO ISUB=1,NSBPRC+1
*--- Central result
                           WTXTMP(IBIN,ISUB,IORD) =
     >                          MYRESN(IBIN,ISUB,IORD)
*--- Rel. lower uncertainty
                           WTXLTMP(IBIN,ISUB,IORD) =
     >                          WTDXL2(IBIN,ISUB,IORD)
*--- Rel. upper uncertainty
                           WTXUTMP(IBIN,ISUB,IORD) =
     >                          WTDXU2(IBIN,ISUB,IORD)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO

*--- HERAPDF1.0: 2nd PDF set for additional PDF uncertainty
*---             with quadratic addition of
*---             all lower/upper deviations 
               CALL INITPDFSET(PDFSET2(1:LEN_TRIM(PDFSET2)))
               CALL INITPDF(0)
               IPHASE  = 1
               IMODE   = 1
               IWEIGHT = 0
               NRAP  = NRAPIDITY
               IF (LRAT.OR.LNRM) THEN
                  NRAP = 2*NRAPIDITY
               ENDIF
               CALL FX9999IN(FILENAME)
               CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
               ISTEP = 0
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))
               IF (LNRM) THEN
*--- Load normalization table with potentially different binning!
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAMN)
                     CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                  ENDIF
                  ISTEP = 1
                  CALL CENRES(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAME)
                     CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                  ENDIF
               ENDIF
               ISTEP = 2
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))
               
               CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
               IPHASE = 2
               
*--- HERAPDF1.0: Do loop runs from 1 - 8 for this part
               DO J=1,8
                  CALL INITPDF(J)
                  CALL FX9999IN(FILENAME)
                  CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                  IF (LNRM) THEN
                     ISTEP = 3
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                     IF (LTAB) THEN
                        CALL FX9999IN(FILENAMN)
                        CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                     ENDIF
                     ISTEP = 4
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                     IF (LTAB) THEN
                        CALL FX9999IN(FILENAME)
                        CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                     ENDIF
                     ISTEP = 5
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                  ENDIF
                  CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
               ENDDO
               
               IPHASE = 3
               CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

*--- Add quadratically to previously backed-up result
               IBIN = 0
               DO IRAP=1,NRAP
                  DO IPT=1,NPT(IRAP)
                     IBIN = IBIN+1
                     DO IORD=1,NORD+1
                        DO ISUB=1,NSBPRC+1
*--- Central result
*--- Rel. lower uncertainty
                           WTXLTMP(IBIN,ISUB,IORD) =
     >                          -SQRT(WTXLTMP(IBIN,ISUB,IORD)**2.+
     >                          WTDXL2(IBIN,ISUB,IORD)**2.)
*--- Rel. upper uncertainty
                           WTXUTMP(IBIN,ISUB,IORD) =
     >                          +SQRT(WTXUTMP(IBIN,ISUB,IORD)**2.+
     >                          WTDXU2(IBIN,ISUB,IORD)**2.)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO

*--- HERAPDF1.0: 2nd PDF set for additional PDF uncertainty
*---             with quadratic addition of
*---             minimal/maximal lower/upper deviations 
               CALL INITPDF(0)
               IPHASE  = 1
               IMODE   = 3
               IWEIGHT = 0
               NRAP  = NRAPIDITY
               IF (LRAT.OR.LNRM) THEN
                  NRAP = 2*NRAPIDITY
               ENDIF
               CALL FX9999IN(FILENAME)
               CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
               ISTEP = 0
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))
               IF (LNRM) THEN
*--- Load normalization table with potentially different binning!
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAMN)
                     CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                  ENDIF
                  ISTEP = 1
                  CALL CENRES(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAME)
                     CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                  ENDIF
               ENDIF
               ISTEP = 2
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))
               
               CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
               IPHASE = 2

*--- HERAPDF1.0: Do loop runs from 9 - 13 for this part
               DO J=9,13
                  CALL INITPDF(J)
                  CALL FX9999IN(FILENAME)
                  CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                  IF (LNRM) THEN
                     ISTEP = 3
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                     IF (LTAB) THEN
                        CALL FX9999IN(FILENAMN)
                        CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                     ENDIF
                     ISTEP = 4
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                     IF (LTAB) THEN
                        CALL FX9999IN(FILENAME)
                        CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                     ENDIF
                     ISTEP = 5
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                  ENDIF
                  CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
               ENDDO
               
               IPHASE = 3
               CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
               
*--- Add quadratically to previously backed-up result
               IBIN = 0
               DO IRAP=1,NRAP
                  DO IPT=1,NPT(IRAP)
                     IBIN = IBIN+1
                     DO IORD=1,NORD+1
                        DO ISUB=1,NSBPRC+1
*--- Central result
*--- Rel. lower uncertainty
                           WTXLTMP(IBIN,ISUB,IORD) =
     >                          -SQRT(WTXLTMP(IBIN,ISUB,IORD)**2.+
     >                          WTDXLM(IBIN,ISUB,IORD)**2.)
*--- Rel. upper uncertainty
                           WTXUTMP(IBIN,ISUB,IORD) =
     >                          +SQRT(WTXUTMP(IBIN,ISUB,IORD)**2.+
     >                          WTDXUM(IBIN,ISUB,IORD)**2.)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO

*--- Copy back to original arrays as for single PDF set use
               IBIN = 0
               DO IRAP=1,NRAP
                  DO IPT=1,NPT(IRAP)
                     IBIN = IBIN+1
                     DO IORD=1,NORD+1
                        DO ISUB=1,NSBPRC+1
*--- Central result
                           MYRESN(IBIN,ISUB,IORD) =
     >                          WTXTMP(IBIN,ISUB,IORD)
*--- Rel. lower uncertainty
                           WTDXL2(IBIN,ISUB,IORD) =
     >                          WTXLTMP(IBIN,ISUB,IORD)
*--- Rel. upper uncertainty
                           WTDXU2(IBIN,ISUB,IORD) =
     >                          WTXUTMP(IBIN,ISUB,IORD)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF



c - Give some standard output, fill histograms
            IBIN = 0
            DO IRAP=1,NRAP
               DO IPT=1,NPT(IRAP)
                  IBIN = IBIN+1
                  WRITE(*,900) IBIN,MYRESN(IBIN,NSBPRC+1,NORD+1),
     >                 WTDXL2(IBIN,NSBPRC+1,NORD+1),
     >                 WTDXU2(IBIN,NSBPRC+1,NORD+1)
               ENDDO
            ENDDO
            
c - Fill histograms
            CALL PDFFILL(NRAP,0,-1,I,MYRESN)
            CALL PDFFILL(NRAP,1,-1,I,WTDXL2)
            CALL PDFFILL(NRAP,2,-1,I,WTDXU2)
         ENDDO                     ! Loop over scales
      ENDIF
*---  End of PDF uncertainties



ckr The following should not have been changed in PDF part!
*---  Central scale
ckr      ISCL = ISCLPR
*--- NLO contribution 2, scale dimension 1
ckr      XMUR = SCALEFAC(2,1,ISCL)
ckr      XMUF = SCALEFAC(2,1,ISCL)
ckr      XMUR = XMURS(1)
ckr      XMUF = XMUFS(1)
c - Make sure to use again the central PDF!
      CALL INITPDFSET(PDFSET(1:LEN_TRIM(PDFSET)))
      CALL INITPDF(0)
      CALL FX9999IN(FILENAME)
      CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)



c - Statistics part
c - Use statistics tables
c - Call statistical error-code for scenario
      IF (LSTAT) THEN
         WRITE(*,*)"****************************************"//
     >        "********************************"
         WRITE(*,*)"ALLUNC: Evaluating statistical uncertainties"
         WRITE(*,*)"****************************************"//
     >        "********************************"
         BORNNAME = TABPATH(1:LEN_TRIM(TABPATH))//"/stat/"//
     >        SCENARIO(1:LEN_TRIM(SCENARIO))//"-hhc-born-"
         NLONAME  = TABPATH(1:LEN_TRIM(TABPATH))//"/stat/"//
     >        SCENARIO(1:LEN_TRIM(SCENARIO))//"-hhc-nlo-"
*---  Central scale
         ISCL = ISCLPT(1)
ckr NLO contribution 2, scale dimension 1
ckr         XMUR = SCALEFAC(2,1,ISCL)
ckr         XMUF = SCALEFAC(2,1,ISCL)
         CALL INITPDF(0)
         DO IORD=0,MIN(2,NORD)
            IORD2 = IORD
            IF (IORD.EQ.0) THEN
               IORD2 = NORD+1
            ENDIF
            NFOUND = 0
            NFAIL  = 0
            IPHASE  = 1
            IMODE   = 1
ckr            IWEIGHT = 1
            IWEIGHT = 0
c - This is the normal table to fill the default values
            CALL FX9999IN(FILENAME)
            CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
            IPHASE = 2

c - Loop over files
            IF (IORD.EQ.0) THEN
c - Total x section
               NTAB = NLON 
               BWGT = 1D0/1D8
               FILEBASE = NLONAME(1:LEN_TRIM(NLONAME))
            ELSEIF (IORD.EQ.1) THEN
c - LO
               NTAB = BORNN
               BWGT = 1D0/1D9
               FILEBASE = BORNNAME(1:LEN_TRIM(BORNNAME))
            ELSEIF (IORD.EQ.2) THEN
c - NLO
               NTAB = NLON
               BWGT = 1D0/1D8
               FILEBASE = NLONAME(1:LEN_TRIM(NLONAME))
            ELSE
               WRITE(*,*)"ALLUNC: ERROR! "//
     >              "Illegal order for stat. calc:",IORD
               STOP
            ENDIF
            DO ITAB=0,NTAB
               WRITE(NO,'(I4.4)'),ITAB
               FILENAMES = FILEBASE(1:LEN_TRIM(FILEBASE))//"2jet_"//NO
     >              //".tab"
               OPEN(2,STATUS='OLD',FILE=FILENAMES,IOSTAT=ISTAT)
               IF (ISTAT.NE.0) THEN
                  FILENAMES = FILEBASE(1:LEN_TRIM(FILEBASE))//"3jet_"
     >                 //NO//".tab"
                  OPEN(2,STATUS='OLD',FILE=FILENAMES,IOSTAT=ISTAT)
                  IF (ISTAT.NE.0) THEN
C---  WRITE(*,*)"Filename for order",IORD,":",FILENAMES
                     NFAIL = NFAIL + 1
                     IF (NFAIL.LT.2) THEN
                        WRITE(*,*)"ALLUNC: WARNING! Table file "//
     >                       "not found, skipped ! IOSTAT = ",ISTAT
                     ENDIF
ckr While using f90 DO-ENDDO ... one could also use the EXIT statement
ckr instead of GOTO. However, EXIT leaves the loop!
ckr To continue the loop use CYCLE instead! 
                     GOTO 20
ckr            CYCLE
                  ENDIF
               ENDIF
               CLOSE(2)
               NFOUND = NFOUND + 1
               IF (NFOUND.LT.2) THEN
                  WRITE(*,*)"Filename for order",IORD,":",FILENAMES
               ENDIF
               
               NRAP  = NRAPIDITY
               IMODE = 2
               IF (LRAT) THEN
                  NRAP = 2*NRAPIDITY
               ENDIF
c - These are the stat. tables to get the deviations
               CALL FX9999IN(FILENAMES)
               CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
               CALL UNCERT(IPHASE,IMODE,IWEIGHT,ITAB,LRAT,LNRM)
 20            CONTINUE
            ENDDO
            
            WRITE(*,*)"########################################"//
     >           "################################"
            WRITE(*,*)"No. of filenames skipped:",NFAIL
            WRITE(*,*)"No. of filenames analyzed:",NFOUND
            WRITE(*,*)"########################################"//
     >           "################################"
            
            IPHASE = 3
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

c - Special statistics output

ckr            WRITE(*,*
ckr     >           )"\n *************************************************"
ckr            WRITE(*,*)"ALLUNC: No looping over scales for order:",IORD
ckr            WRITE(*,*
ckr     >           )"*************************************************"

            WRITE(*,*)"========================================"//
     >           "===================================="//
     >           "=========================="//
     >           "==========================="
            WRITE(*,*)"#IBIN #IMIN #IMAX         <s>           "//
     >           "s_min           s_max       ds/<s>/%"//
     >           " ds_min/<s>/% ds_max/<s>/%"//
     >           "     ds_min/ds    ds_max/ds"
            WRITE(*,*)"----------------------------------------"//
     >           "------------------------------------"//
     >           "--------------------------"//
     >           "---------------------------"

            IBIN = 0
            DO IRAP=1,NRAP
               DO IPT=1,NPT(IRAP)
                  IBIN = IBIN+1
                  DSTMP(3) = -1D0
                  DSTMP(4) = -1D0
                  IF (MYRES(IBIN,NSBPRC+1,IORD2).GE.0.D0) THEN
                     DSTMP(1) = (WTXMIN(IBIN,NSBPRC+1,IORD2)-1D0)
     >                    *100D0
                     DSTMP(2) = (WTXMAX(IBIN,NSBPRC+1,IORD2)-1D0)
     >                    *100D0
                     IF (WTDXU2(IBIN,NSBPRC+1,IORD2).GT.1D-99) THEN
                        DSTMP(3) = 
     >                       (WTXMIN(IBIN,NSBPRC+1,IORD2)-1D0) /
     >                       WTDXU2(IBIN,NSBPRC+1,IORD2)
                        DSTMP(4) = 
     >                       (WTXMAX(IBIN,NSBPRC+1,IORD2)-1D0) /
     >                       WTDXU2(IBIN,NSBPRC+1,IORD2)
                     ENDIF
                  ELSE
                     DSTMP(1) = (WTXMIN(IBIN,NSBPRC+1,IORD2)+1D0)
     >                    *100D0
                     DSTMP(2) = (WTXMAX(IBIN,NSBPRC+1,IORD2)+1D0)
     >                    *100D0
                     IF (WTDXU2(IBIN,NSBPRC+1,IORD2).GT.1D-99) THEN
                        DSTMP(3) = 
     >                       (WTXMIN(IBIN,NSBPRC+1,IORD2)+1D0) /
     >                       WTDXU2(IBIN,NSBPRC+1,IORD2)
                        DSTMP(4) = 
     >                       (WTXMAX(IBIN,NSBPRC+1,IORD2)+1D0) /
     >                       WTDXU2(IBIN,NSBPRC+1,IORD2)
                     ENDIF
                  ENDIF
                  WRITE(*,902) IBIN,
     >                 IJMIN(IBIN),
     >                 IJMAX(IBIN),
     >                 MYRES(IBIN,NSBPRC+1,IORD2),
     >                 WTXMIN(IBIN,NSBPRC+1,IORD2) *
     >                 DABS(MYRES(IBIN,NSBPRC+1,IORD2)),
     >                 WTXMAX(IBIN,NSBPRC+1,IORD2) *
     >                 DABS(MYRES(IBIN,NSBPRC+1,IORD2)),
     >                 WTDXMN(IBIN,NSBPRC+1,IORD2)*100D0,
     >                 DSTMP(1),DSTMP(2),
     >                 DSTMP(3),DSTMP(4)
               ENDDO
            ENDDO

c - End special statistics output



c - Give some standard output, fill histograms
            WRITE(*,*)"========================================"//
     >           "================================"
            WRITE(*,*)"Statistical Uncertainties"
            WRITE(*,*)"- the printed values are for the total "//
     >           "cross section, primary scale (3 in v1.4, 1 in v2), "//
     >           "summed over all subprocesses"
            WRITE(*,*)"- histograms contain more detailed results"
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            WRITE(*,*)" bin       cross section           "//
     >           "average error of the mean"
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            WRITE(*,*)"ALLUNC: Now order no.",iord
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            IBIN = 0
            DO IRAP=1,NRAP
               DO IPT=1,NPT(IRAP)
                  IBIN = IBIN+1
                  WRITE(*,901) IBIN,MYRES(IBIN,NSBPRC+1,IORD2),
     >                 WTDXMN(IBIN,NSBPRC+1,IORD2)
               ENDDO
            ENDDO
c - Fill histograms
            CALL PDFFILL(NRAP,3,IORD,ISCL,WTDXMN)
            CALL PDFFILL(NRAP,4,IORD,ISCL,WTDXUL)
            
         ENDDO
      ENDIF
*---  End of statistical uncertainty
      
      
      
Comment: c - alpha_s uncertainty part
Comment: c - Use primary table
Comment: c - Check that FILENAME is still the primary table here ...!!!
Comment: c - alpha_s(M_Z) variations, either implicit via PDF set or explicitly
Comment: c - Uncertainty filled at central scale no. 3
Comment: c - (ISCL=3 in FORTRAN, refscale=2 in C++ parlance of author code)
Comment:       IF (LSER) THEN
Comment:          WRITE(*,*)"****************************************"//
Comment:      >        "********************************"
Comment:          WRITE(*,*)"ALLUNC: Evaluating alpha_s uncertainties"
Comment:          WRITE(*,*)"****************************************"//
Comment:      >        "********************************"
Comment: ckr Attention: Assumption here is that delta(alpha_s) = +- 0.00n
Comment: ckr            equals a change in PDF set member number by n
Comment: ckr            This is true for:
Comment: ckr                              cteq66alphas:           2 +- up to  2
Comment: ckr                              CT10as:                 5 +- up to  5
Comment: ckr                              MSTW2008nlo_asmzrange: 11 +- up to 10
Comment: ckr Central member number  : Integer part of DASMZVAL
Comment: ckr Up/down member numbers : Integer part of (DASMZVAL-NINT(DASMZVAL)) times 1000 
Comment: ckr alpha_s variation alone: (DASMZVAL-NINT(DASMZVAL))
Comment:          CALL INITPDF(0)
Comment:          CALL NUMBERPDF(MYPDFAS)
Comment:          IPDF = NINT(ABS(DASMZVAL))
Comment:          WRITE(*,*)"ALLUNC: PDF member for central alpha_s value:"
Comment:      >        ,IPDF
Comment:          WRITE(*,*)"ALLUNC: Number of alpha_s variations in set:"
Comment:      >        ,MYPDFAS
Comment:          IF (ASMODE.EQ."PDF".AND.IPDF.GT.MYPDFAS) THEN
Comment:             WRITE(*,*)"ALLUNC: ERROR! Central PDF member for,"//
Comment:      >           " alpha_s variation does not exist, aborting!"//
Comment:      >           "IPDF = ",IPDF
Comment:             STOP
Comment:          ENDIF
Comment:          IPDFUD = NINT(1000.D0*(ABS(DASMZVAL)-DBLE(IPDF)))
Comment:          WRITE(*,*)"ALLUNC: alpha_s variation * 1000 here:"
Comment:      >        ,IPDFUD
Comment:          IF (ASMODE.EQ."PDF".AND.
Comment:      >        (IPDF-IPDFUD.LT.0.OR.IPDF+IPDFUD.GT.MYPDFAS)) THEN
Comment:             WRITE(*,*)"ALLUNC: ERROR! Varied PDF members for,"//
Comment:      >           " alpha_s variation do not exist, aborting!"//
Comment:      >           "IPDFUD = ",IPDFUD
Comment:             STOP
Comment:          ENDIF
Comment: 
Comment:          WRITE(*,*)"========================================"//
Comment:      >        "================================"
Comment:          WRITE(*,*)"Relative alpha_s Uncertainties"
Comment:          WRITE(*,*)"- the printed values are for the total "//
Comment:      >        "cross section summed over all subprocesses"
Comment:          WRITE(*,*)"- histograms contain more detailed results"
Comment:          WRITE(*,*)"----------------------------------------"//
Comment:      >        "--------------------------------"
Comment:          WRITE(*,*)" bin       cross section           "//
Comment:      >        "lower a_s uncertainty   upper a_s uncertainty"
Comment:          WRITE(*,*)"----------------------------------------"//
Comment:      >        "--------------------------------"
Comment:          
Comment:          ISCL = 3
Comment:          XMUR = MURSCALE(ISCL)
Comment:          XMUF = MUFSCALE(ISCL)
Comment: 
Comment:          IF (ASMODE.EQ."PDF") THEN
Comment: ckr alpha_s variations in PDF set
Comment:             IPHASE  = 1
Comment:             IMODE   = 3
Comment:             IWEIGHT = 0
Comment:             NRAP  = NRAPIDITY
Comment:             IF (LRAT.OR.LNRM) THEN
Comment:                NRAP = 2*NRAPIDITY
Comment:             ENDIF
Comment:             CALL INITPDF(IPDF)
Comment: cnew
Comment:             ASMZPDF = ALPHASPDF(ZMASS)
Comment: cnew
Comment:             CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:             ISTEP = 0
Comment:             CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             IF (LNRM) THEN
Comment: ckr Load normalization table with potentially different binning!
Comment:                IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
Comment:                ISTEP = 1
Comment:                CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:                IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:             ENDIF
Comment:             ISTEP = 2
Comment:             CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment: 
Comment:             CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
Comment:             IPHASE = 2
Comment:             
Comment:             CALL INITPDF(IPDF-IPDFUD)
Comment: cnew
Comment:             ASMZTMP = ALPHASPDF(ZMASS)
Comment:             ASDN = ASMZTMP - ASMZPDF
Comment: cnew
Comment:             CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:             IF (LNRM) THEN
Comment:                ISTEP = 3
Comment:                CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:                IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
Comment:                ISTEP = 4
Comment:                CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:                IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:                ISTEP = 5
Comment:                CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             ENDIF
Comment:             CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
Comment:             
Comment:             CALL INITPDF(IPDF+IPDFUD)
Comment: cnew
Comment:             ASMZTMP = ALPHASPDF(ZMASS)
Comment:             ASUP = ASMZTMP - ASMZPDF
Comment: cnew
Comment:             CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:             IF (LNRM) THEN
Comment:                ISTEP = 3
Comment:                CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:                IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
Comment:                ISTEP = 4
Comment:                CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:                IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:                ISTEP = 5
Comment:                CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             ENDIF
Comment:             CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
Comment:             
Comment:             IPHASE = 3
Comment:             CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
Comment: 
Comment: c - Give some standard output, fill histograms
Comment:             WRITE(*,*)"ALLUNC: Uncertainties from"//
Comment:      >           " alpha_s variations in PDF members"
Comment:             WRITE(*,*)"----------------------------------------"//
Comment:      >           "--------------------------------"
Comment:             ISCL = 3
Comment:             IBIN   = 0
Comment:             DO IRAP=1,NRAP
Comment:                DO IPT=1,NPT(IRAP)
Comment:                   IBIN = IBIN+1
Comment:                   WRITE(*,900) IBIN,MYRESN(IBIN,NSBPRC+1,NORD+1),
Comment:      >                 WTDXLM(IBIN,NSBPRC+1,NORD+1),
Comment:      >                 WTDXUM(IBIN,NSBPRC+1,NORD+1)
Comment:                ENDDO
Comment:             ENDDO
Comment:             CALL PDFFILL(NRAP,0,-1,ISCL,MYRESN)
Comment:             CALL PDFFILL(NRAP,1,-1,ISCL,WTDXLM)
Comment:             CALL PDFFILL(NRAP,2,-1,ISCL,WTDXUM)
Comment:          ENDIF
Comment: 
Comment: ckr Standalone alpha_s variations
Comment:          IPHASE  = 1
Comment:          IMODE   = 3
Comment:          IWEIGHT = 0
Comment:          NRAP  = NRAPIDITY
Comment:          IF (LRAT.OR.LNRM) THEN
Comment:             NRAP = 2*NRAPIDITY
Comment:          ENDIF
Comment:          CALL INITPDF(IPDF)
Comment:          ASMODETMP = ASMODE
Comment:          ASMODE    = "KR"
Comment:          ASMZTMP   = ASMZVAL
Comment:          ASMZVAL   = ALPHASPDF(ZMASS)
Comment: 
Comment:          CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:          ISTEP = 0
Comment:          CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:          IF (LNRM) THEN
Comment: ckr Load normalization table with potentially different binning!
Comment:             IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
Comment:             ISTEP = 1
Comment:             CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:          ENDIF
Comment:          ISTEP = 2
Comment:          CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:          
Comment:          CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
Comment:          IPHASE = 2
Comment:          
Comment: cnew         ASMZVAL = ASMZVAL - DBLE(IPDFUD)/1000.
Comment:          ASMZVAL = ASMZVAL + ASDN
Comment:          CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:          IF (LNRM) THEN
Comment:             ISTEP = 3
Comment:             CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >           SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
Comment:             ISTEP = 4
Comment:             CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >           SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:             ISTEP = 5
Comment:             CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >           SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:          ENDIF
Comment:          CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
Comment:          
Comment: ckr Add twice to correct previous subtraction!
Comment: cnew         ASMZVAL = ASMZVAL + DBLE(2*IPDFUD)/1000.
Comment:          ASMZVAL = ASMZVAL - ASDN + ASUP
Comment:          CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:          IF (LNRM) THEN
Comment:             ISTEP = 3
Comment:             CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >           SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
Comment:             ISTEP = 4
Comment:             CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >           SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:             ISTEP = 5
Comment:             CALL CENRES(ISTEP,LRAT,LNRM,
Comment:      >           SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:          ENDIF
Comment:          CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
Comment:          
Comment:          IPHASE = 3
Comment:          CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
Comment: 
Comment:          ASMODE  = ASMODETMP
Comment:          ASMZVAL = ASMZTMP
Comment: 
Comment: c - Give some standard output, fill histograms
Comment:          WRITE(*,*)"ALLUNC: Uncertainties from"//
Comment:      >        " standalone alpha_s variations"
Comment:          WRITE(*,*)"----------------------------------------"//
Comment:      >        "--------------------------------"
Comment:          ISCL = 3
Comment:          IBIN   = 0
Comment:          DO IRAP=1,NRAP
Comment:             DO IPT=1,NPT(IRAP)
Comment:                IBIN = IBIN+1
Comment:                WRITE(*,900) IBIN,MYRESN(IBIN,NSBPRC+1,NORD+1),
Comment:      >              WTDXLM(IBIN,NSBPRC+1,NORD+1),
Comment:      >              WTDXUM(IBIN,NSBPRC+1,NORD+1)
Comment:             ENDDO
Comment:          ENDDO
Comment:          IF (.NOT.ASMODE.EQ."PDF") THEN
Comment:             CALL PDFFILL(NRAP,0,-1,ISCL,MYRESN)
Comment:          ENDIF
Comment:          CALL PDFFILL(NRAP,3,-1,ISCL,WTDXLM)
Comment:          CALL PDFFILL(NRAP,4,-1,ISCL,WTDXUM)
Comment:       ENDIF
      


c - Scale uncertainty part
c - Use primary table
c - Check that FILENAME is still the primary table here ...!!!
c - Two schemes are implemented for hadron-hadron:
c - 1. NSCALEVAR = 4 with standard mur, muf factors of
c      (1/4,1/4), (1/2,1/2), (  1,  1), (  2,  2)
c - 2. NSCALEVAR = 8 with additional mur, muf combinations
c      (1/4,1/4), (1/2,1/2), (  1,  1), (  2,  2) as before plus
c      (  1,1/2), (  1,  2), (1/2,  1), (  2,  1)
c - Uncertainty filled at central scale no. 3 (v14) or 1 (v2)
c - (ISCL=3 in FORTRAN, refscale=2 in C++ parlance of author code)
      IF (LSCL) THEN
         WRITE(*,*)"****************************************"//
     >        "********************************"
         WRITE(*,*)"ALLUNC: Evaluating scale uncertainties"
         WRITE(*,*)"****************************************"//
     >        "********************************"
c - 2-point scheme
         IF (NSCLS.GE.3) THEN
            CALL INITPDF(0)
*---  Central scale
            ISCL = ISCLPT(1)
            XMUR = XMURS(1)
            XMUF = XMUFS(1)
            IPHASE  = 1
            IMODE   = 3
            IWEIGHT = 0
            NRAP  = NRAPIDITY
            IF (LRAT.OR.LNRM) THEN
               NRAP = 2*NRAPIDITY
            ENDIF
            CALL FX9999IN(FILENAME)
            CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
            ISTEP = 0
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
            IF (LNRM) THEN
ckr Load normalization table with potentially different binning!
               IF (LTAB) THEN
                  CALL FX9999IN(FILENAMN)
                  CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
               ENDIF
               ISTEP = 1
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))
               IF (LTAB) THEN
                  CALL FX9999IN(FILENAME)
                  CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
               ENDIF
            ENDIF
            ISTEP = 2
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
            
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
            IPHASE = 2

            DO I=2,3
               ISCL = ISCLPT(I)
               XMUR = XMURS(I)
               XMUF = XMUFS(I)
               CALL FX9999IN(FILENAME)
               CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
               IF (LNRM) THEN
                  ISTEP = 3
                  CALL CENRES(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
ckr Load normalization table with potentially different binning!
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAMN)
                     CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                  ENDIF
                  ISTEP = 4
                  CALL CENRES(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAME)
                     CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                  ENDIF
                  ISTEP = 5
                  CALL CENRES(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
               ENDIF
               CALL UNCERT(IPHASE,IMODE,IWEIGHT,ISCL,LRAT,LNRM)
            ENDDO
            IPHASE = 3
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

*---  Give some standard output, fill histograms
            WRITE(*,*)"========================================"//
     >           "================================"
            WRITE(*,*)"Relative Scale Uncertainties (2-point)"
            WRITE(*,'(A,I2,A)')
     >           " - the printed values are for the total "/
     >           /"cross section, scale no. ",ISCLPT(1),","
            WRITE(*,*)"  summed over all subprocesses"
            WRITE(*,*)"- histograms contain more detailed results"
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            WRITE(*,*)" bin       cross section           "//
     >           "lower scale uncertainty upper scale uncertainty"
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            WRITE(*,*)"ALLUNC: Uncertainties from"//
     >           " symmetric (2-point) scale variations"
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            IORD   = 0
            ISCL   = ISCLPT(1)
            ISUB   = 0
            IBIN   = 0
            DO IRAP=1,NRAP
               DO IPT=1,NPT(IRAP)
                  IBIN = IBIN+1
                  WRITE(*,900) IBIN,MYRESN(IBIN,NSBPRC+1,NORD+1),
     >                 WTDXLM(IBIN,NSBPRC+1,NORD+1),
     >                 WTDXUM(IBIN,NSBPRC+1,NORD+1)
               ENDDO
            ENDDO

*---  Central scale
            ISCL = ISCLPT(1)
            CALL PDFFILL(NRAP,6,-1,ISCL,WTDXLM)
            CALL PDFFILL(NRAP,7,-1,ISCL,WTDXUM)
         ENDIF

c - 6-point scheme
         IF (NSCLS.GE.7) THEN
            CALL INITPDF(0)
*---  Central scale
            ISCL = ISCLPT(1)
            XMUR = XMURS(1)
            XMUF = XMUFS(1)
            IPHASE  = 1
            IMODE   = 3
            IWEIGHT = 0
            NRAP  = NRAPIDITY
            IF (LRAT.OR.LNRM) THEN
               NRAP = 2*NRAPIDITY
            ENDIF
            CALL FX9999IN(FILENAME)
            CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
            ISTEP = 0
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
            IF (LNRM) THEN
ckr Load normalization table with potentially different binning!
               IF (LTAB) THEN
                  CALL FX9999IN(FILENAMN)
                  CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
               ENDIF
               ISTEP = 1
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))
               IF (LTAB) THEN
                  CALL FX9999IN(FILENAME)
                  CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
               ENDIF
            ENDIF
            ISTEP = 2
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))

            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
            IPHASE = 2

            DO I=2,7
               ISCL = ISCLPT(I)
               XMUR = XMURS(I)
               XMUF = XMUFS(I)
               CALL FX9999IN(FILENAME)
               CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
               IF (LNRM) THEN
                  ISTEP = 3
                  CALL CENRES(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
ckr Load normalization table with potentially different binning!
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAMN)
                     CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                  ENDIF
                  ISTEP = 4
                  CALL CENRES(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAME)
                     CALL FX9999CC(XMUR,XMUF,XSECT0,XSUNCOR,XSCOR)
                  ENDIF
                  ISTEP = 5
                  CALL CENRES(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
               ENDIF
               CALL UNCERT(IPHASE,IMODE,IWEIGHT,ISCL,LRAT,LNRM)
            ENDDO
            IPHASE = 3
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

c - Give some standard output, fill histograms
            WRITE(*,*)"========================================"//
     >           "================================"
            WRITE(*,*)"Relative Scale Uncertainties (6-point)"
            WRITE(*,'(A,I2,A)')
     >           " - the printed values are for the total "/
     >           /"cross section, scale no. ",ISCLPT(1),","
            WRITE(*,*)"  summed over all subprocesses"
            WRITE(*,*)"- histograms contain more detailed results"
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            WRITE(*,*)" bin       cross section           "//
     >           "lower scale uncertainty upper scale uncertainty"
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            WRITE(*,*)"ALLUNC: Uncertainties from"//
     >           " asymmetric (6-point) scale variations"
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            IORD   = 0
            ISCL   = ISCLPT(1)
            ISUB   = 0
            IBIN   = 0
            DO IRAP=1,NRAP
               DO IPT=1,NPT(IRAP)
                  IBIN = IBIN+1
                  WRITE(*,900) IBIN,MYRESN(IBIN,NSBPRC+1,NORD+1),
     >                 WTDXLM(IBIN,NSBPRC+1,NORD+1),
     >                 WTDXUM(IBIN,NSBPRC+1,NORD+1)
               ENDDO
            ENDDO
            CALL PDFFILL(NRAP,8,-1,ISCL,WTDXLM)
            CALL PDFFILL(NRAP,9,-1,ISCL,WTDXUM)
         ENDIF
      ENDIF
*---  End of scale uncertainties



Comment: c - Algorithmic part
Comment: c - Use reference table
Comment: c - Reference result is always stored in scale variation no. 1
Comment: c - Reference table has doubled number of rapidity bins and
Comment: c -   the reference result is stored in upper half of these bins
Comment: c - Reference result is evaluated with CTEQ61 PDFs
Comment: c - Default scale in hadron-hadron usually is the third variation:
Comment: c -   C++ no. 2 --> Fortran no. 3 ==> normal result (but compare to author code!)
Comment: c - (Use other scale in FORTRAN ONLY if refscale <> 2 in author code)
Comment: c - Attention: From now on ref. table loaded ==> rap. bins doubled,
Comment: c                  no ratio calcs below ...
Comment:       LRAT = .FALSE.
Comment:       IF (LALG) THEN
Comment:          WRITE(*,*)"****************************************"//
Comment:      >        "********************************"
Comment:          WRITE(*,*)"ALLUNC: Evaluating algorithmic uncertainties"
Comment:          WRITE(*,*)"****************************************"//
Comment:      >        "********************************"
Comment:          FILENAME = TABPATH(1:LEN_TRIM(TABPATH))//"/"//REFNAME
Comment:          WRITE(*,*)"----------------------------------------"//
Comment:      >        "--------------------------------"
Comment:          WRITE(*,*)"ALLUNC: Taking reference table: "//
Comment:      >        FILENAME(1:LEN_TRIM(FILENAME))
Comment: c - Initialize CTEQ61 reference PDFs
Comment:          PDFSET = PDFPATH(1:LEN_TRIM(PDFPATH))//"/cteq61.LHgrid"
Comment:          WRITE(*,*)"ALLUNC: Taking reference PDF: "//
Comment:      >        PDFSET(1:LEN_TRIM(PDFSET))
Comment:          WRITE(*,*)"----------------------------------------"//
Comment:      >        "--------------------------------"
Comment:          CALL INITPDFSET(PDFSET(1:LEN_TRIM(PDFSET)))
Comment:          CALL INITPDF(0)
Comment: 
Comment: c - Reference result
Comment: c - Need scale variation no. 1 for this, upper half of rap. bins
Comment: c - NOTE1: The values of mur and muf are internally used to select the
Comment: c          the corresponding table number, so they have to be set
Comment: c          according to scale variation 1 although the factors are wrong
Comment: c          for the comparison ==>
Comment: c            compare to scale variation 3, lower half of rap. bins
Comment: c - NOTE2: The reference result should not depend on the mur or muf values
Comment: c          apart from the table choice, BUT if mur is not equal to muf
Comment: c          additional factors are taken into account which must be
Comment: c          avoided for the reference calculation!
Comment:          ISCL = 1
Comment:          XMUR = MURSCALE(ISCL)
Comment:          XMUF = MUFSCALE(ISCL)
Comment: ckr         CALL FX9999CC(FILENAME,XMUR,XMUF,1,XSECT0)
Comment:          CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:          IPHASE  = 1
Comment:          IMODE   = 4
Comment:          IWEIGHT = 0
Comment:          CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:          CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
Comment: 
Comment: c - Normal result
Comment:          ISCL = 3
Comment:          XMUR = MURSCALE(ISCL)
Comment:          XMUF = MUFSCALE(ISCL)
Comment:          CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
Comment:          IPHASE = 2
Comment:          CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
Comment:          IPHASE = 3
Comment:          CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
Comment: 
Comment: c - Give some standard output, fill histograms
Comment:          WRITE(*,*)"========================================"//
Comment:      >        "================================"
Comment:          WRITE(*,*)"Relative Algorithmic Error"
Comment:          WRITE(*,*)"- the printed values are for the total "//
Comment:      >        "cross section summed over all subprocesses"
Comment:          WRITE(*,*)"- histograms contain more detailed results"
Comment:          WRITE(*,*)"----------------------------------------"//
Comment:      >        "--------------------------------"
Comment:          WRITE(*,*)" bin       ref. cross section      "//
Comment:      >        "algorithmic deviation"
Comment:          WRITE(*,*)"----------------------------------------"//
Comment:      >        "--------------------------------"
Comment:          ISCL = 3
Comment:          IBIN   = 0
Comment:          DO IRAP=1,INT(NRAPIDITY/2)
Comment:             DO IPT=1,NPT(IRAP)
Comment:                IBIN = IBIN+1
Comment:                WRITE(*,901) IBIN,WTX(IBIN,NSBPRC+1,NORD+1),
Comment:      >              WTDXUM(IBIN,NSBPRC+1,NORD+1)
Comment:             ENDDO
Comment:          ENDDO
Comment:          CALL PDFFILL(NRAPIDITY,5,-1,ISCL,WTDXUM)
Comment:       ENDIF



c - Close hbook file
      CALL PDFHIST(2,HISTFILE,LONE,LPDF,LSTAT,LALG,LSER,MYPDF,
     >     LRAT.OR.LNRM,LSCL,ISCLPT(1))
      END



c
c ======================= Book the histograms ========================
c
      SUBROUTINE PDFHIST(N,HISTFILE,
     &     LONE,LPDF,LSTAT,LALG,LSER,MYPDF,LRAT,LSCL,ISCLPR)
      IMPLICIT NONE
      CHARACTER*(*) HISTFILE
      CHARACTER*255 CSTRNG,CBASE1,CBASE2,CTMP
      INTEGER N,IPDF,MYPDF,IPTMAX,NRAP,MYNSCLS,ISCLPR
      LOGICAL LONE,LPDF,LSTAT,LALG,LSER,LRAT,LSCL
      
      INTEGER I,J,ISTAT2,ICYCLE
      INTEGER IORD,ISUB,ISCL,IRAP,IPT,IHIST,NHIST
      INCLUDE "fnx9999.inc"
      INCLUDE "strings.inc"
      INCLUDE "uncert.inc"
      REAL PT(NPTMAX)

c - HBOOK common 
      INTEGER NWPAWC
      PARAMETER (NWPAWC=2500000)
      REAL HMEMOR(NWPAWC)
      COMMON /PAWC/ HMEMOR

c - Open & book
      IF (N.EQ.1) THEN
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         WRITE(*,*)"PDFHIST: Book Histograms"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         CALL HLIMIT(NWPAWC)
         CALL HROPEN(11,"fastNLO",HISTFILE,"N",1024,ISTAT2)
         IF (ISTAT2.NE.0) THEN
            WRITE(*,*)"\nPDFHIST: ERROR! Could not open histofile: ",
     >           ISTAT2," Aborted!"
            STOP
         ENDIF
         
         NRAP = NRAPIDITY
         IF (LRAT) THEN
            NRAP = 2*NRAPIDITY
            DO I=NRAPIDITY+1,NRAP
               NPT(I) = NPT(I-NRAPIDITY)
               DO J=1,(NPT(I)+1)
                  PTBIN(I,J) = PTBIN(I-NRAPIDITY,J)
               ENDDO
            ENDDO
         ENDIF
         NHIST = 0
         CBASE1 = "CIPROC"
ckr         CBASE1 = CIPROC(IPROC)
         
         CBASE1 = CBASE1(1:LEN_TRIM(CBASE1))//"_"
     >        //"NAMELABEL1"
ckr     >        //NAMELABEL(1)
         CBASE2 = CBASE1(1:LEN_TRIM(CBASE1))//"_"
     >        //"CIALGO"
ckr     >        //CIALGO(IALGO)
         CBASE2 = CBASE2(1:LEN_TRIM(CBASE2))//"_"
     >        //"CJETRES1"
ckr     >        //CJETRES1(IALGO)
ckr         WRITE(CTMP,'(F3.1)'),JETRES1
         CTMP = "JETRES1"
         CBASE2 = CBASE2(1:LEN_TRIM(CBASE2))//"="
     >        //CTMP
ckr         DO IORD=0,NORD         ! Order: tot, LO, NLO-corr, NNLO-corr
         DO IORD=0,2            ! Order: tot, LO, NLO-corr, NNLO-corr
            MYNSCLS = NSCLS
            IF (NSCALEMAX.GE.8.AND.NSCLS.EQ.4) THEN
               MYNSCLS = 8
            ENDIF
            DO ISCL=1,MYNSCLS   ! Scale variations
ckr NLO contribution 2
               DO ISUB=0,NSBPRC ! Subprocesses: 0 tot + 7 subproc
                  DO IRAP=1, NRAP
                     IHIST = IORD*1000000 + ISCL*100000 +
     >                    ISUB*10000 + IRAP*100
ckr                     write(*,*)"ALL: iord,isc,isub,irap,ihist",
ckr     >                    iord,iscale,isub,irap,ihist
                     DO J=1,(NPT(IRAP)+1)
                        PT(J) = REAL(PTBIN(IRAP,J))
                     ENDDO
                     CSTRNG = CBASE2
                     WRITE(CTMP,'(I1)'),IORD
                     CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//"_iord="
     >                    //CTMP
                     WRITE(CTMP,'(I1)'),ISCL
                     CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//"_imu="
     >                    //CTMP
                     WRITE(CTMP,'(I1)'),IRAP
                     CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//"_iy="
     >                    //CTMP
                     CALL HBOOKB(IHIST,
     >                    CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                    NPT(IRAP),PT,0)
ckr                     CALL HBARX(IHIST)
                     NHIST = NHIST+1
ckr                     write(*,*)"1. Booked histo #",nhist
Comment:                      IF (LSER) THEN
Comment:                         DO IPDF=1,MYPDF
Comment:                            CALL HBOOKB(IHIST+IPDF,
Comment:      >                          CSTRNG(1:LEN_TRIM(CSTRNG)),
Comment:      >                          NPT(IRAP),PT,0)
Comment:                            NHIST = NHIST+1
Comment:                         ENDDO
Comment:                      ENDIF
                     IF (LPDF) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_PDF_low/xsect"
                        CALL HBOOKB(IHIST+1,
     >                       CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_PDF_up/xsect"
                        CALL HBOOKB(IHIST+2,
     >                       CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+2
ckr                        write(*,*)"2. Booked histo #",nhist
                     ELSEIF (LSER.AND.ISCL.EQ.ISCLPR) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_asPDF_low/xsect"
                        CALL HBOOKB(IHIST+1,
     >                       CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_asPDF_up/xsect"
                        CALL HBOOKB(IHIST+2,
     >                       CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+2
ckr                        write(*,*)"2. Booked histo #",nhist
                     ENDIF
                     IF (LSTAT) THEN
                        IF (IORD.LE.2.AND.
     >                       ISCL.EQ.ISCLPR.AND.ISUB.EQ.0) THEN
                           CSTRNG = CBASE1
                           CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                          "_dstat/xsect"
                           CALL HBOOKB(IHIST+3,
     >                          CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                          NPT(IRAP),PT,0)
                           CSTRNG = CBASE1
                           CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                          "_dmax/2/xsect"
                           CALL HBOOKB(IHIST+4,
     >                          CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                          NPT(IRAP),PT,0)
                           NHIST = NHIST+2
ckr                        write(*,*)"3. Booked histo #",nhist
                        ENDIF
                     ELSEIF (LSER.AND.ISCL.EQ.ISCLPR) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_as_low/xsect"
                        CALL HBOOKB(IHIST+3,
     >                       CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_as_up/xsect"
                        CALL HBOOKB(IHIST+4,
     >                       CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+2
ckr                        write(*,*)"3. Booked histo #",nhist
                     ENDIF
                     IF (LALG.AND.IORD.LE.2.AND.
     >                    ISCL.EQ.ISCLPR) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_xsect/xsect_ref"
                        CALL HBOOKB(IHIST+5,
     >                       CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
ckr                        write(*,*)"4. Booked histo #",nhist
                     ENDIF
                     IF (LSCL.AND.IORD.LE.2.AND.
     >                    ISCL.EQ.ISCLPR) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_SCL2_low/xsect"
                        CALL HBOOKB(IHIST+6,
     >                       CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_SCL2_up/xsect"
                        CALL HBOOKB(IHIST+7,
     >                       CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_SCL6_low/xsect"
                        CALL HBOOKB(IHIST+8,
     >                       CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_SCL6_up/xsect"
                        CALL HBOOKB(IHIST+9,
     >                       CSTRNG(1:LEN_TRIM(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
ckr                        write(*,*)"4b. Booked histo #",nhist
                     ENDIF
                     IF (LSTAT.AND.ISUB.EQ.0) THEN
                        IHIST = IORD*1000000 + ISCL*100000 +
     >                       ISUB*10000 + IRAP*100
                        CSTRNG = CBASE1(1:LEN_TRIM(CBASE1))
                        WRITE(CTMP,'(I1)'),IORD
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//"_iord="
     >                       //CTMP
                        WRITE(CTMP,'(I1)'),ISCL
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//"_imu="
     >                       //CTMP
                        WRITE(CTMP,'(I1)'),IRAP
                        CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//"_iy="
     >                       //CTMP
                        CTMP = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_norm_mu_sig_all_pt"
                        CALL HBOOK1(IHIST + 10,
     >                       CTMP(1:LEN_TRIM(CTMP)),
     >                       63,-10.5,10.5,0)
                        CALL HIDOPT(IHIST + 10,'STAT')
                        CTMP = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                       "_norm_mu_dmax_all_pt"
                        CALL HBOOK1(IHIST + 11,
     >                       CTMP(1:LEN_TRIM(CTMP)),
     >                       30,-1.5,1.5,0)
                        CALL HIDOPT(IHIST + 11,'STAT')
                        NHIST = NHIST+2
ckr                        write(*,*)"5. Booked histo #, IHIST",nhist,
ckr     >                       IHIST + 11
ckr IHIST limit before next rapidity bin: IHIST+2*IPT+11 < IHIST+100
ckr => Maximal IPT = IPTMAX < 45
                        IPTMAX = 44
                        DO IPT=1,MIN(NPT(IRAP),IPTMAX)
                           CSTRNG = CBASE1(1:LEN_TRIM(CBASE1))
                           WRITE(CTMP,'(I1)'),IORD
                           CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//"_iord="
     >                          //CTMP
                           WRITE(CTMP,'(I1)'),ISCL
                           CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))/
     >                          /"_imu="//CTMP
                           WRITE(CTMP,'(I1)'),IRAP
                           CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//"_iy="
     >                          //CTMP
                           WRITE(CTMP,'(I2)'),IPT
                           CSTRNG = CSTRNG(1:LEN_TRIM(CSTRNG))//"_ipt="
     >                          //CTMP
                           CTMP = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                          "_norm_mu_sig"
                           CALL HBOOK1(IHIST + 2*IPT + 10,
     >                          CTMP(1:LEN_TRIM(CTMP)),
     >                          63,-10.5,10.5,0)
                           CALL HIDOPT(IHIST + 2*IPT + 10,'STAT')
                           CTMP = CSTRNG(1:LEN_TRIM(CSTRNG))//
     >                          "_norm_mu_dmax"
                           CALL HBOOK1(IHIST + 2*IPT + 11,
     >                          CTMP(1:LEN_TRIM(CTMP)),
     >                          30,-1.5,1.5,0)
                           CALL HIDOPT(IHIST + 2*IPT + 11,'STAT')
                           NHIST = NHIST+2
ckr                           write(*,*)"6. Booked histo #, IHIST",nhist,
ckr     >                          IHIST + 2*IPT + 11
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO               ! End od ISCL loop, IORD still on
         ENDDO
         WRITE(*,*)"PDFHIST: Number of histograms booked:",NHIST
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"



c - Close HBOOK file
      ELSEIF (N.EQ.2) THEN
         CALL HROUT(0,ICYCLE," ")
         CALL HREND("fastNLO")
      ENDIF

      RETURN
      END

c
c ======================= Fill the histograms =========================
c
      SUBROUTINE PDFFILL(NRAP,IOFF,IORDFIL,ISCL,DVAL)
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      INCLUDE "uncert.inc"
      INTEGER NRAP,IOFF,IORDFIL,ISCL
      DOUBLE PRECISION DVAL(MXOBSBIN,MXSUBPROC+1,4)
      INTEGER I,J,IBIN,IORD,IORD2,ISUB,ISUB2,IHIST
      REAL RVAL
      
      IF (ISCL.LT.1 .OR. ISCL.GT.NSCLS) THEN
         WRITE(*,*) "\nPDFFILL: ERROR! ISCL ",ISCL,
     >        " is out of range, aborted!"
         WRITE(*,*) "PDFFILL: Max. ISCL: ",NSCLS
         STOP
      ENDIF
      
c - Fill all histograms for the given scale
c - Fill sums over subprocesses and/or orders into zero factors for IHIST
      DO IORD2=1,NORD+1         ! Order: LO, NLO-corr, NNLO-corr, ... , tot
         IORD = IORD2
         IF (IORD2.EQ.NORD+1) IORD = 0 
         IF (IORDFIL.EQ.-1.OR.IORDFIL.EQ.IORD) THEN
            DO ISUB2=1,NSBPRC+1 ! Subprocesses hh: 6/7 subprocesses + total
               ISUB = ISUB2
               IF (ISUB2.EQ.NSBPRC+1) ISUB=0
               IBIN=0
               DO I=1,NRAP
                  DO J=1,NPT(I)
                     IBIN = IBIN + 1
ckr Recall: HBOOK understands only single precision
                     RVAL  = REAL(DVAL(IBIN,ISUB2,IORD2))
                     IHIST = IORD*1000000+ISCL*100000+ISUB*10000+I*100
                     IHIST = IHIST+IOFF
                     CALL HFILL(IHIST,REAL(PTBIN(I,J)+0.01),0.0,RVAL)
                  ENDDO         ! pT-loop
               ENDDO            ! rap-loop
            ENDDO               ! isub-loop
         ENDIF
      ENDDO                     ! iord-loop
      
      RETURN
      END
