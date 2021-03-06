      SUBROUTINE ALLUNC
* ---------------------------------------------------------------------
* K. Rabbertz 07.09.2008 First try to integrate all uncertainties
*                        into one job
*
* ALLUNC - Program to derive all (PDF, statistical, algorithmic,
*          scale ...) uncertainties using fastNLO tables
*
* 09.04.2013 kr: Write out killfile for largest fluctuations
* 14.02.2013 kr: Improve PDF uncertainty setup
* 13.08.2012 kr: Implement symmetrized eigen vector method,
*                small bug fix in asymm. eigen vector method
* 10.06.2011 kr: Replace CERNLIB function LENOCC by f95 Standard LEN_TRIM
* 31.03.2011 kr: Some modifications to improve v2 compatibility
* ---------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      INCLUDE "uncert.inc"
      INCLUDE "v14unc.inc"
      CHARACTER*255 SCENARIO,TABPATH,TABNAME,REFNAME
      CHARACTER*255 FILENAME,FILENAMES,HISTFILE,KILLFILE
      CHARACTER*255 TABNAMN,FILENAMN
      CHARACTER*255 PDFNAM,PDFSET,PDFSET2,PDFPATH,LHAPDF,CHTMP
      CHARACTER*255 FILEBASE,LOFILE,NLOFILE
      CHARACTER*255 BORNNAME,NLONAME
      CHARACTER*8 CH8TMP
      CHARACTER*4 CH4TMP
      CHARACTER*4 NO
      INTEGER BORNN,NLON
      INTEGER I,J,MYPDF,IPDF,IPDFUD,MYPDFAS,IOPDF,IOAS
      INTEGER ITAB,NTAB,NFOUND,NFAIL,IOPEN
      INTEGER ISTAT,ISCL,IORD,IORD2,IBIN,NBIN,ISUB,IRAP,IPT,NTMP
      INTEGER IHIST,IPHASE,ISTEP,IETYPE,NPDFMOD,NPDFPAR
      LOGICAL LONE,LPDF,LSTAT,LSER,LSCL,LRAT,LALG,LNRM,LTAB

      INTEGER IPRINT
      LOGICAL LLO,LNLO,LTHC1L,LTHC2L,LNPC1,LDATA
      DOUBLE PRECISION ALPS,FNALPHAS,ALPHASPDF,ASMZPDF,ASUP,ASDN
      DOUBLE PRECISION XMUR,XMUF,QLAM4,QLAM5,BWGT
      DOUBLE PRECISION DSTMP(MXOBSBIN,4)
      DOUBLE PRECISION WTXTMP(MXOBSBIN,MXSUBPROC+1,NMAXORD+1)
      DOUBLE PRECISION WTXLTMP(MXOBSBIN,MXSUBPROC+1,NMAXORD+1)
      DOUBLE PRECISION WTXUTMP(MXOBSBIN,MXSUBPROC+1,NMAXORD+1)

*---  Arrays for crosscheck on PDF uncertainties
      INTEGER K,L,ii
      INTEGER IDim0Bin,IDim1Bin,IDim2Bin
      INTEGER NDim0Bins,NDim1Bins,NDim2Bins
      DOUBLE PRECISION XS0(2*MXOBSBIN),DXU(2*MXOBSBIN),DXL(2*MXOBSBIN)
      DOUBLE PRECISION DXU2(2*MXOBSBIN),DXL2(2*MXOBSBIN)
      DOUBLE PRECISION DXS0U,DXS0L,DELTA1,DELTA2,DELTA3

c - To unify quoted uncertainties (CL68,CL90,special)
c - Convert from CL68 to CL90 values
c - TOCL90 = 1.64485D0 ! SQRT(2.D0)*InvERF(0.9D0)
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
*---  Use argument counter to simplify logic
      INTEGER NARG
      DATA NARG/0/

*---  ATTENTION: This is the most likely source of Fortran problems!
*---  For each scenario, the result array must be declared at least
*---  as large as in the definition in the common block of the
*---  corresponding scenario.
*---  See the scenario info printout "Tot. no. of observable bins"
*---  and compare to the value of the parameter MxObsBin
*---  in the file fnx9999.inc or [scenario].inc.
*---  Adapt the following to your scenario:
*---  Integer MxObsBin
*---  Parameter (MxObsBin = nnn)
*---  We recommend to name the array according to the scenario.
      DOUBLE PRECISION XSLO(MXOBSBIN),XSCLLO(MXOBSBIN)
      DOUBLE PRECISION XSNLO(MXOBSBIN),XSCLNLO(MXOBSBIN)
      DOUBLE PRECISION XSTHC(MXOBSBIN),XSCLTHC(MXOBSBIN)
      DOUBLE PRECISION DXSUCTMP(MXOBSBIN,2),DXSCORTMP(MXOBSBIN,2)
      DOUBLE PRECISION XSNPC(MXOBSBIN),XSCLNPC(MXOBSBIN)
      DOUBLE PRECISION DXSUCNPC(MXOBSBIN,2),DXSCORNPC(MXOBSBIN,2)
      DOUBLE PRECISION XSDAT(MXOBSBIN),XSCLDAT(MXOBSBIN)
      DOUBLE PRECISION DXSUCDATA(MXOBSBIN,2),DXSCORDATA(MXOBSBIN,2)
      DOUBLE PRECISION KFAC(MXOBSBIN),KTHC(MXOBSBIN),KNPC(MXOBSBIN)

      REAL PT(NPTMAX)
      INTEGER IMODE,IWEIGHT,NRAP

ckr Old Z mass
ckr      PARAMETER (ZMASS = 91.187D0)
ckr Z mass from PDG 2012
      DOUBLE PRECISION ZMASS
      PARAMETER (ZMASS = 91.1876D0)

      CHARACTER*255 ASMODE,ASMODETMP
      DOUBLE PRECISION ASMZVAL,DASMZVAL,ASMZTMP
      INTEGER IASLOOP
      COMMON/STEER/ASMZVAL,IASLOOP,ASMODE

c --- Set debug printout level
      IDEBUG = 0

c --- Parse command line
      WRITE(*,*)"########################################"//
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
      NARG = NARG + 1
      IF (IARGC().LT.NARG) THEN
         SCENARIO = "fnt2003"
         WRITE(*,*)
     >        "ALLUNC: WARNING! No scenario name given, "//
     >        "taking the default fnt2003 instead!"
         WRITE(*,*)"      For an explanation of command line "//
     >        "arguments type:"
         WRITE(*,*)"      ./allunc -h"
      ELSE
         CALL GETARG(NARG,SCENARIO)
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
     >           'share/lhapdf/PDFsets'
            WRITE(*,*)'  alpha_s calc., def. = PDF (from PDF set)'
            WRITE(*,*)'    alt. = PY: 0-, 1- and 2-loop '//
     >           '(from Pythia 6.4 using Lambda_4 from PDF)'
            WRITE(*,*)'    alt. = KR: 1-, 2- and 3-loop '//
     >           '(from hep-ph/9506442)'
            WRITE(*,*)'    alt. = MW: 2-, 3- and 4-loop '//
     >           '(from hep-ph/9806404)'
            WRITE(*,*)'  alpha_s(M_Z), def. from PDF set'
            WRITE(*,*)'     (in mode PY this has to be Lambda_4/GeV!)'
            WRITE(*,*)'  alpha_s(M_Z) variation, def. = 0.D0, i.e. none'
            WRITE(*,*)'  alpha_s loop order, def. from PDF set'
            WRITE(*,*)'  PDF uncertainties in '
            WRITE(*,*)'    1. asymmetric pairwise eigen vector method'
            WRITE(*,*)'       (CTEQ/MSTW),'
            WRITE(*,*)'    2. asymmetric pairwise eigen vector method'
            WRITE(*,*)'       adding also same-side deviations,'
            WRITE(*,*)'    3. symmetric pairwise eigen vector method'
            WRITE(*,*)'       (HERAPDF),'
            WRITE(*,*)'    4. toy MC method (NNPDF),'
            WRITE(*,*)'    5. symmetrized eigen vector method'
            WRITE(*,*)'       (AB(K)M),'
            WRITE(*,*)'    6. asymmetric eigen vector method,'
            WRITE(*,*)'       identical to method 2,'
            WRITE(*,*)'    def. = 1.'
            WRITE(*,*)'    Add 10 to include HERAPDF parameterization '
            WRITE(*,*)'    uncertainties.'
            WRITE(*,*)' '
            STOP
         ELSE
            WRITE(*,*)"ALLUNC: Evaluating scenario: ",
     >           SCENARIO(1:LEN_TRIM(SCENARIO))
ckr To be cross-checked for each new scenario
            IF (SCENARIO(1:8).EQ."fnl2332c") THEN
               LNRM = .TRUE.
               LTAB = .FALSE.
               WRITE(*,*)
     >              "ALLUNC: Deriving x section ratios"
            ELSEIF (SCENARIO(1:7).EQ."fnl2380") THEN
               LNRM = .TRUE.
               LTAB = .TRUE.
               WRITE(*,*)
     >              "ALLUNC: Deriving x section ratios"
            ELSEIF (SCENARIO(1:7).EQ."fnl2442") THEN
ckr Original version for fnl2442: Works fine, trivial division in rap 3
               LRAT = .FALSE.
ckr New norm. version for fnl2442: Works fine, trivial division in rap 4
               LNRM = .TRUE.
               LTAB = .FALSE.
               WRITE(*,*)
     >              "ALLUNC: Deriving x section ratios"
            ELSEIF (SCENARIO(1:11).EQ."fnl2522diff") THEN
               LNRM = .TRUE.
               LTAB = .TRUE.
               WRITE(*,*)
     >              "ALLUNC: Deriving normalized distributions"
            ELSEIF (SCENARIO(1:7).EQ."fnl2622".OR.
     >              SCENARIO(1:7).EQ."fnl3622".OR.
     >              SCENARIO(1:7).EQ."fnl2652") THEN
               LNRM = .TRUE.
               LTAB = .FALSE.
               WRITE(*,*)
     >              "ALLUNC: Deriving normalized distributions"
            ELSEIF (SCENARIO(1:10).EQ."fnl2722num".OR.
     >              SCENARIO(1:10).EQ."fnl2732num".OR.
     >              SCENARIO(1:10).EQ."fnl2742num") THEN
               LNRM = .TRUE.
               LTAB = .TRUE.
               WRITE(*,*)
     >              "ALLUNC: Deriving normalized distributions"
            ENDIF
         ENDIF
      ENDIF
      TABNAME = SCENARIO(1:LEN_TRIM(SCENARIO))//".tab"
      IF (LNRM.AND.LTAB) THEN
         TABNAMN = SCENARIO(1:7)//"norm"//".tab"
      ENDIF
      REFNAME = SCENARIO(1:LEN_TRIM(SCENARIO))//"ref.tab"

*---Path to tables
      NARG = NARG + 1
      TABPATH = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,TABPATH)
      ENDIF
      IF (IARGC().LT.NARG.OR.TABPATH(1:1).EQ."_") THEN
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
      NARG = NARG + 1
      HISTFILE = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,HISTFILE)
      ENDIF
      IF (IARGC().LT.NARG.OR.HISTFILE(1:1).EQ."_") THEN
         HISTFILE = SCENARIO(1:LEN_TRIM(SCENARIO))//".hbk"
         WRITE(*,*)
     >        "ALLUNC: WARNING! No output filename given, "//
     >        "taking scenario.hbk instead!"
      ELSE
         WRITE(*,*)"ALLUNC: Creating output file: ",
     >        HISTFILE(1:LEN_TRIM(HISTFILE))
      ENDIF

*---Derive algorithmic uncertainty?
      NARG = NARG + 1
      CH4TMP = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,CH4TMP)
      ENDIF
      IF (IARGC().LT.NARG.OR.CH4TMP(1:1).EQ."_".OR.
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
      NARG = NARG + 1
      CH4TMP = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,CH4TMP)
      ENDIF
      IF (IARGC().LT.NARG.OR.CH4TMP(1:1).EQ."_") THEN
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
      NARG = NARG + 1
      CH4TMP = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,CH4TMP)
      ENDIF
      IF (IARGC().LT.NARG.OR.CH4TMP(1:1).EQ."_") THEN
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
      NARG = NARG + 1
      CH4TMP = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,CH4TMP)
      ENDIF
      IF (IARGC().LT.NARG.OR.CH4TMP(1:1).EQ."_") THEN
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
      NARG = NARG + 1
      PDFSET = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,PDFSET)
      ENDIF
      IF (IARGC().LT.NARG.OR.PDFSET(1:1).EQ."_") THEN
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
      NARG = NARG + 1
      CHTMP = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,CHTMP)
      ENDIF
      IF (IARGC().LT.NARG.OR.CHTMP(1:1).EQ."_") THEN
         PDFPATH = "/share/lhapdf/PDFsets"
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
      NARG = NARG + 1
      CHTMP = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,CHTMP)
      ENDIF
      IF (IARGC().LT.NARG.OR.CHTMP(1:1).EQ."_") THEN
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
      NARG = NARG + 1
      CH8TMP = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,CH8TMP)
      ENDIF
      IF (IARGC().LT.NARG.OR.CH8TMP(1:1).EQ."_") THEN
         ASMZVAL = -1D0
         WRITE(*,*)
     >        "ALLUNC: No alpha_s(M_Z) value given, "//
     >        "using alpha_s according to PDF set"
      ELSE
         READ(CH8TMP,'(F9.6)'),ASMZVAL
         WRITE(*,*)"ALLUNC: Using alpha_s(M_Z):",ASMZVAL
      ENDIF

*---alpha_s(M_Z) +- variation
      NARG = NARG + 1
      CH8TMP = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,CH8TMP)
      ENDIF
      IF (IARGC().LT.NARG.OR.CH8TMP(1:1).EQ."_") THEN
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
      NARG = NARG + 1
      CH4TMP = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,CH4TMP)
      ENDIF
      IF (IARGC().LT.NARG.OR.CH4TMP(1:1).EQ."_") THEN
         IASLOOP = -1
         WRITE(*,*)
     >        "ALLUNC: No alpha_s loop order given, "//
     >        "using alpha_s loop order according to PDF set"
      ELSE
         READ(CH4TMP,'(I4)'),IASLOOP
         WRITE(*,*)"ALLUNC: Using alpha_s loop order:",IASLOOP
      ENDIF

*---Determine method to use for PDF uncertainties
      NARG = NARG + 1
      CH4TMP = "X"
      IF (IARGC().GE.NARG) THEN
         CALL GETARG(NARG,CH4TMP)
      ENDIF
      IF (IARGC().LT.NARG.OR.CH4TMP(1:1).EQ."_") THEN
         IETYPE = 1
         WRITE(*,*)
     >        "ALLUNC: Use asymmetric eigen vector "//
     >        "(CTEQ/MSTW) method."
      ELSE
         READ(CH4TMP,'(I4)'),IETYPE
         IF (MOD(IETYPE,10).EQ.1) THEN
            WRITE(*,*)
     >           "ALLUNC: Use asymmetric pairwise eigen vector "//
     >           "method (CTEQ/MSTW)."
         ELSEIF (MOD(IETYPE,10).EQ.2) THEN
            WRITE(*,*)
     >           "ALLUNC: Use asymmetric pairwise eigen vector "//
     >           "method adding also same-side deviations."
         ELSEIF (MOD(IETYPE,10).EQ.3) THEN
            WRITE(*,*)
     >           "ALLUNC: Use symmetric pairwise eigen vector "//
     >           "method (HERAPDF)."
         ELSEIF (MOD(IETYPE,10).EQ.4) THEN
            WRITE(*,*)
     >           "ALLUNC: Use toy MC method (NNPDF)."
         ELSEIF (MOD(IETYPE,10).EQ.5) THEN
            WRITE(*,*)
     >           "ALLUNC: Use symmetrized eigen vector "//
     >           "method (AB(K)M)."
         ELSEIF (MOD(IETYPE,10).EQ.6) THEN
            WRITE(*,*)
     >           "ALLUNC: Use asymmetric eigen vector "//
     >           "method, identical to method 2."
         ELSE
            WRITE(*,*)
     >           "ALLUNC: ERROR! Undefined PDF error method, aborted!"//
     >           " IETYPE = ",IETYPE
         ENDIF
         IF (INT(IETYPE/10).EQ.0) THEN
         ELSEIF (INT(IETYPE/10).EQ.1) THEN
            WRITE(*,*)
     >           "ALLUNC: Add quadratically parameterization "//
     >           "uncertainties of HERAPDF."
         ELSE
            WRITE(*,*)
     >           "ALLUNC: ERROR! Undefined PDF error method, aborted!"//
     >           " IETYPE = ",IETYPE
         ENDIF
      ENDIF

*---Too many arguments
      NARG = NARG + 1
      IF (IARGC().GT.NARG) THEN
         WRITE(*,*)"\nALLUNC: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF

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

*---  Check uncertainties to derive
      LONE  = MYPDF.LE.1
      LSTAT = BORNN.GE.2.OR.NLON.GE.2
      IF ((LALG.OR.LNRM.OR.LRAT).AND.LSTAT) THEN
         LALG = .FALSE.
         LNRM = .FALSE.
         LRAT = .FALSE.
         WRITE(*,*)
     >        "ALLUNC: WARNING! Statistics mode!"
         WRITE(*,*)"        Algorithmic, normalization "//
     >        "or ratio uncertainties are not possible!"
      ENDIF
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

*---  Determine pointers to contributions and print out contribution list
      Call FX9999PT(1D0,1D0,1)
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
ctmp         IF (ICONTRFLAG1(I).EQ.2.AND.ICONTRFLAG2(I).EQ.1)LTHC1L = .TRUE.
ctmp         IF (ICONTRFLAG1(I).EQ.2.AND.ICONTRFLAG2(I).EQ.2)LTHC2L = .TRUE.
         IF (ICONTRFLAG1(I).EQ.4.AND.ICONTRFLAG2(I).EQ.1.AND.
     >        IADDMULTFLAG(I).EQ.1)
     >        LNPC1 = .TRUE.
         IF (ICONTRFLAG1(I).EQ.0.AND.ICONTRFLAG2(I).EQ.0.AND.
     >        IDATAFLAG(I).EQ.1)
     >        LDATA = .TRUE.
      ENDDO

*---  Determine dimensional subdivisions (NDIM=1,2 only!)
ckr      IF (NDIM.GT.2) THEN
Comment:       j = NDim0Bins()
Comment:       if (ndim.gt.1) then
Comment:          do i=1,j
Comment:             k = NDim1Bins(i)
Comment:             if (ndim.gt.2) then
Comment:                do ii=1,k
Comment:                   l = NDim2Bins(i,ii)
Comment:                   write(*,*)"NDim0Bins,NDim1Bins,NDim2Bins",j,k,l
Comment:                enddo
Comment:             else
Comment:                write(*,*)"NDim0Bins,NDim1Bins",j,k
Comment:             endif
Comment:          enddo
Comment:       else
Comment:          write(*,*)"NDim0Bins",j
Comment:       endif

Comment:          DO I=1,NOBSBIN
Comment:             j = IDim0Bin(i)
Comment:             k = IDim1Bin(i)
Comment:             l = IDim2Bin(i)
Comment:             write(*,*)"iobsbin,IDim0Bin,IDim1Bin,IDim2Bin",i,j,k,l
Comment:          enddo
Comment:          WRITE(*,*)"ALLUNC: ERROR! Cannot deal yet with "//
Comment:      >        "more than 2 dimensions. NDIM = ",NDIM
ckr         STOP
ckr      ENDIF
      DO I=1,NOBSBIN
*---  Very first bin; initialize
         IF (I.EQ.1) THEN
            NTMP = 1
            NBINCOUNTER(1) = 0
            NBINCOUNTER(2) = 1
*---  Bin of 2nd dimension identical
         ELSEIF ((LOBIN(I-1,1).EQ.LOBIN(I,1)).AND.
     >           (UPBIN(I-1,1).EQ.UPBIN(I,1))) THEN
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
c      write(*,*)"NBINCOUNTER(1)",NBINCOUNTER(1)
c      write(*,*)"NBINCOUNTER(2),NDim0Bins()",NBINCOUNTER(2),NDim0Bins()
      DO J=1,NBINCOUNTER(2)
c         write(*,*)"J,NDIVCOUNTER(J),NDim1Bins(J)",
c     >        J,NDIVCOUNTER(J),NDim1Bins(J)
         DO I=1,NDIVCOUNTER(J)
            IBIN = IBIN + 1
            IOBSPOINTER(I,J) = IBIN
c            write(*,*)"I,J,IBIN,IOBSPOINTER",I,J,
c     >           IBIN,IOBSPOINTER(I,J)
         ENDDO
      ENDDO

*---  Define v14 variables
ckr      NRAPIDITY = NBINCOUNTER(2)
      NRAPIDITY = NDim0Bins()
      DO I=1,NRAPIDITY
         NPT(I) = NDIVCOUNTER(I)
         DO J=1,NPT(I)
            PTBIN(I,J) = LOBIN(IOBSPOINTER(J,I),2)
         ENDDO
         PTBIN(I,NPT(I)+1) = UPBIN(IOBSPOINTER(NPT(I),I),1)
         RAPBIN(I) = LOBIN(IOBSPOINTER(1,I),1)
      ENDDO
      RAPBIN(NRAPIDITY+1) = UPBIN(IOBSPOINTER(1,NRAPIDITY),2)

      IF (IDEBUG.GT.1) THEN
         WRITE(*,*)"DEBUG2: QQQ1 NBINCOUNTER 1, 2  = ",
     >        NBINCOUNTER(1),NBINCOUNTER(2)
         DO I=1,NOBSBIN
            WRITE(*,*)"DEBUG2: QQQ2: IOBS, LOBIN I,1-2 = ",
     >           I,LOBIN(I,1),LOBIN(I,2)
            WRITE(*,*)"DEBUG2: QQQ3: IOBS, UPBIN I,1-2 = ",
     >           I,UPBIN(I,1),UPBIN(I,2)
         ENDDO
         DO I=1,NRAPIDITY
            NPT(I) = NDIVCOUNTER(I)
            WRITE(*,*)"DEBUG2: QQQ4: IRAP, RAPBIN, NPT = ",
     >           I,RAPBIN(I),NPT(I)
            DO J=1,NPT(I)+1
               WRITE(*,*)"DEBUG2: QQQ5: IRAP, IPT, PTBIN  = ",
     >              I,J,PTBIN(I,J)
            ENDDO
         ENDDO
         WRITE(*,*)"DEBUG2: QQQ4: IRAP, RAPBIN      = ",
     >        NRAPIDITY+1,RAPBIN(NRAPIDITY+1)
      ENDIF

*---  Store no. of subprocesses in NSBPRC
*---  v14 NSUBPROC       : always 7
*---  v20 NSUBPROC(iCtrb): 2-parton processes: 6   (e.g. LO dijet)
*---  .    3-parton processes: 7   (e.g. NLO dijet or LO 3-jet)
*---  .    4-parton processes: 7   (e.g. NLO 3-jet)
*---  ATTENTION: Uncertainty derivation always assumes summed-up
*---  .          subprocesses to be found in IORD bin 7 + 1 !
*---  .          In case of LO dijet with 6 subprocesses only this has
*---  .          to be properly initialized to 0 !
      NSBPRC = NSUBPROC(2)

*---  Define output formats
ckr 900     FORMAT(1P,I5,3(3X,E21.14))
 900  FORMAT(1P,I5,3(6X,E18.11))
 901  FORMAT(1P,I5,2(6X,E18.11))
 902  FORMAT(3I6,3E16.5,5(F10.3,3X))
 910  FORMAT(1P,I5,3(6X,E18.11),0P,2X,1(3X,F9.5))
 920  FORMAT(1P,I5,3(6X,E18.11),0P,2X,2(3X,F9.5))

*---  Initial settings
      CALL FNSET("P_RESET",0)   ! Reset all selections to zero

*---  Set v14 variable NORD to maximal numbering value of
*---  ILO, INLO, ITHC1L, ITHC2L
*---  Make sure that:
*---  1. NORD is equal to maximal ICONT counter
*---  2. Unwanted contributions in result() array are filled with zeros
*---  ==> avoid double counting of 1- and 2-loop threshold corrections!
      IF (.NOT.LLO) THEN
         WRITE(*,*)"ALLUNC: ERROR! No LO found, stopped!"
         STOP
      ELSE
         NORD = ILO
         CALL FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
      ENDIF
      IF (LNLO) THEN
         NORD = INLO
         CALL FNSET("P_ORDPTHY",2) ! select order pert. theory: 1=LO, 2=NLO
         WRITE(*,*)"ALLUNC: INFO! Evaluating NLO."
      ENDIF
      IF (LNLO.AND.LTHC2L) THEN
         NORD = ITHC2L
         CALL FNSET("P_THRESHCOR",2) ! select 2-loop threshold corrections
         WRITE(*,*)"ALLUNC: INFO! "//
     >        "Evaluating 2-loop threshold corrections."
      ENDIF
      IF (.NOT.LNLO.AND.LTHC1L) THEN
         NORD = ITHC1L
         CALL FNSET("P_THRESHCOR",1) ! select 1-loop threshold corrections
         WRITE(*,*)"ALLUNC: INFO! "//
     >        "Evaluating 1-loop threshold corrections."
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
Comment:       CALL PDFHIST(1,HISTFILE,LONE,LPDF,LSTAT,LALG,LSER,MYPDF,
Comment:      >     LRAT.OR.LNRM,LSCL,ISCLPT(1))
      WRITE(*,'(A,I4,A)')" ALLUNC: The observable has",NOBSBIN," bins."
      WRITE(*,'(A,I2,A)')" ALLUNC: The  LO is subdivided into ",
     >     NSUBPROC(1)," subprocesses and"
      WRITE(*,'(A,I2,A)')"         the NLO is subdivided into ",
     >     NSUBPROC(2)," subprocesses."
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

*---  Compute central result and PDF uncertainties
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
         CHTMP = " bin       cross section           "//
     >        "lower PDF uncertainty   upper PDF uncertainty"
         IF (LNLO) THEN
            CHTMP = CHTMP(1:LEN_TRIM(CHTMP))//"   KNLO"
         ENDIF
         IF (LNLO.AND.LTHC2L) THEN
            CHTMP = CHTMP(1:LEN_TRIM(CHTMP))//"        KTHC"
         ELSEIF (LLO.AND.LTHC1L) THEN
            CHTMP = CHTMP(1:LEN_TRIM(CHTMP))//"   KTHC"
         ENDIF
         WRITE(*,*)CHTMP(1:LEN_TRIM(CHTMP))

*---  Only primary scale
         DO I=1,1
            ISCL = ISCLPT(I)
            XMUR = XMURS(I)
            XMUF = XMUFS(I)
*---  Part 1: Basic PDF uncertainty using a single PDF set (and only one
*---          initialization call!). For PDF uncertainty calculation from
*---          multiple PDF sets like full HERAPDF make sure to switch back
*---          to original one here.
*---          Only necessary when running for multiple scales.
            IF (INT(IETYPE/10).NE.0) THEN
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
*--- IMODE=1: Eigen vectors, IMODE=2: Toy MC replicas
            IF (MOD(IETYPE,10).EQ.4) THEN
               IMODE = 2
            ENDIF
            IF (LRAT.OR.LNRM) THEN
               NRAP = 2*NRAPIDITY
            ENDIF
            CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)

*--- Initialize arrays for crosscheck on PDF uncertainties, part 1 (unnormalized)
            DO K=1,NOBSBIN
               XS0(K)  = XSNLO(K)
               DXU(K)  = 0D0
               DXL(K)  = 0D0
               DXU2(K) = 0D0
               DXL2(K) = 0D0
            ENDDO

            ISTEP = 0
            IF (IDEBUG.GT.0)
     >           WRITE(*,*)"DEBUG1: AAA ALLUNC STEP = ",
     >           ISTEP
            CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))

            IF (LNRM) THEN
*--- Load normalization table with potentially different binning!
               IF (LTAB) THEN
                  CALL FX9999IN(FILENAMN)
                  CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
               ENDIF
               ISTEP = 1
               IF (IDEBUG.GT.0)
     >              WRITE(*,*)"DEBUG1: BBB ALLUNC STEP = ",
     >              ISTEP
               CALL CENRES2(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))
               IF (LTAB) THEN
                  CALL FX9999IN(FILENAME)
                  CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
               ENDIF
            ENDIF

            ISTEP = 2
            IF (IDEBUG.GT.0)
     >           WRITE(*,*)"DEBUG1: CCC ALLUNC STEP = ",
     >           ISTEP
            CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
            CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

*--- Initialize arrays for crosscheck on PDF uncertainties, part 2 (normalized)
            IF (LNRM) THEN
               DO K=NOBSBIN+1,2*NOBSBIN
*--- After CENRES STEP2: MYRESN = MYRES
                  XS0(K)  = MYRES(K,NSBPRC+1,NORD+1)
                  DXU(K)  = 0D0
                  DXL(K)  = 0D0
                  DXU2(K) = 0D0
                  DXL2(K) = 0D0
               ENDDO
            ENDIF

*--- Do loop runs once even if MYPDF=0! => Avoid with IF statement
            IF (LPDF) THEN

*--- Pairwise +/- eigen vector evaluation (CTEQ/MSTW, HERAPDF)
*--- (Assumes EV1+-,EV2+-,EV3+- order of PDF members!)
               IF (MOD(IETYPE,10).EQ.1.OR.
     >              MOD(IETYPE,10).EQ.2.OR.
     >              MOD(IETYPE,10).EQ.3) THEN
                  DO J=1,MYPDF-1,2
                     CALL INITPDF(J)
                     CALL FX9999IN(FILENAME)
                     CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                    XSUNCOR,XSCOR)

*--- Fill EV+ into arrays for crosscheck on PDF uncertainties, part 1 (unnormalized)
                     DO K=1,NOBSBIN
                        DXU(K) = MAX(XSNLO(K)-XS0(K),0D0)
                        DXL(K) = MIN(XSNLO(K)-XS0(K),0D0)
                        IF (MOD(IETYPE,10).EQ.2) THEN
                           DXU2(K) = DXU2(K) + DXU(K)*DXU(K)
                           DXL2(K) = DXL2(K) + DXL(K)*DXL(K)
                        ENDIF
                     ENDDO

                     IF (LNRM) THEN
                        ISTEP = 3
                        IF (IDEBUG.GT.0)
     >                       WRITE(*,*)"DEBUG1: DD1A ALLUNC STEP = ",
     >                       ISTEP
                        CALL CENRES2(ISTEP,LRAT,LNRM,
     >                       SCENARIO(1:LEN_TRIM(SCENARIO)))
*--- Load normalization table with potentially different binning!
                        IF (LTAB) THEN
                           CALL FX9999IN(FILENAMN)
                           CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
                        ENDIF
                        ISTEP = 4
                        IF (IDEBUG.GT.0)
     >                       WRITE(*,*)"DEBUG1: EE1A ALLUNC STEP = ",
     >                       ISTEP
                        CALL CENRES2(ISTEP,LRAT,LNRM,
     >                       SCENARIO(1:LEN_TRIM(SCENARIO)))
                        IF (LTAB) THEN
                           CALL FX9999IN(FILENAME)
                           CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
                        ENDIF
                        ISTEP = 5
                        IF (IDEBUG.GT.0)
     >                       WRITE(*,*)"DEBUG1: FF1A ALLUNC STEP = ",
     >                       ISTEP
                        CALL CENRES2(ISTEP,LRAT,LNRM,
     >                       SCENARIO(1:LEN_TRIM(SCENARIO)))
                     ENDIF

*--- Fill EV+ into arrays for crosscheck on PDF uncertainties, part 2 (normalized)
*--- Needs to be MYRES here after CENRES STEP 5 to give result for the
*--- actual PDF
                     IF (LNRM) THEN
                        DO K=NOBSBIN+1,2*NOBSBIN
                           DXU(K) = MAX(MYRES(K,NSBPRC+1,NORD+1)-XS0(K),
     >                          0D0)
                           DXL(K) = MIN(MYRES(K,NSBPRC+1,NORD+1)-XS0(K),
     >                          0D0)
                           IF (MOD(IETYPE,10).EQ.2) THEN
                              DXU2(K) = DXU2(K) + DXU(K)*DXU(K)
                              DXL2(K) = DXL2(K) + DXL(K)*DXL(K)
                           ENDIF
                        ENDDO
                     ENDIF

                     IPHASE = 2
                     CALL UNCERT2(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)

                     CALL INITPDF(J+1)
                     CALL FX9999IN(FILENAME)
                     CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                    XSUNCOR,XSCOR)

*--- Fill EV+/EV- difference into arrays for crosscheck on PDF uncertainties, part 1 (unnormalized)
                     DO K=1,NOBSBIN
                        IF (MOD(IETYPE,10).EQ.1) THEN
                           DXU(K) = MAX(XSNLO(K)-XS0(K),DXU(K),0D0)
                           DXL(K) = MIN(XSNLO(K)-XS0(K),DXL(K),0D0)
                        ELSEIF (MOD(IETYPE,10).EQ.2) THEN
                           DXU(K) = MAX(XSNLO(K)-XS0(K),0D0)
                           DXL(K) = MIN(XSNLO(K)-XS0(K),0D0)
                        ELSEIF (MOD(IETYPE,10).EQ.3) THEN
                           IF (DXU(K).GT.ABS(DXL(K))) THEN
                              DXU(K) = ABS(DXU(K)-(XSNLO(K)-XS0(K)))
     >                             /2D0
                              DXL(K) = -DXU(K)
                           ELSE
                              DXL(K) = -ABS(DXL(K)-(XSNLO(K)-XS0(K)))
     >                             /2D0
                              DXU(K) = -DXL(K)
                           ENDIF
                        ENDIF
                        DXU2(K) = DXU2(K) + DXU(K)*DXU(K)
                        DXL2(K) = DXL2(K) + DXL(K)*DXL(K)
                     ENDDO

                     IF (LNRM) THEN
                        ISTEP = 3
                        IF (IDEBUG.GT.0)
     >                       WRITE(*,*)"DEBUG1: DD1B ALLUNC STEP = ",
     >                       ISTEP
                        CALL CENRES2(ISTEP,LRAT,LNRM,
     >                       SCENARIO(1:LEN_TRIM(SCENARIO)))
*--- Load normalization table with potentially different binning!
                        IF (LTAB) THEN
                           CALL FX9999IN(FILENAMN)
                           CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                          XSUNCOR,XSCOR)
                        ENDIF
                        ISTEP = 4
                        IF (IDEBUG.GT.0)
     >                       WRITE(*,*)"DEBUG1: EE1B ALLUNC STEP = ",
     >                       ISTEP
                        CALL CENRES2(ISTEP,LRAT,LNRM,
     >                       SCENARIO(1:LEN_TRIM(SCENARIO)))
                        IF (LTAB) THEN
                           CALL FX9999IN(FILENAME)
                           CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                          XSUNCOR,XSCOR)
                        ENDIF
                        ISTEP = 5
                        IF (IDEBUG.GT.0)
     >                       WRITE(*,*)"DEBUG1: FF1B ALLUNC STEP = ",
     >                       ISTEP
                        CALL CENRES2(ISTEP,LRAT,LNRM,
     >                       SCENARIO(1:LEN_TRIM(SCENARIO)))
                     ENDIF

*--- Fill EV+/EV- difference into arrays for crosscheck on PDF uncertainties, part 2 (normalized)
*--- Needs to be MYRES here after CENRES STEP 5 to give result for the
*--- actual PDF
                     IF (LNRM) THEN
                        DO K=NOBSBIN+1,2*NOBSBIN
                           IF (MOD(IETYPE,10).EQ.1) THEN
                              DXU(K) = MAX(MYRES(K,NSBPRC+1,NORD+1)-
     >                             XS0(K),DXU(K),0D0)
                              DXL(K) = MIN(MYRES(K,NSBPRC+1,NORD+1)-
     >                             XS0(K),DXL(K),0D0)
                           ELSEIF (MOD(IETYPE,10).EQ.2) THEN
                              DXU(K) = MAX(MYRES(K,NSBPRC+1,NORD+1)-
     >                             XS0(K),0D0)
                              DXL(K) = MIN(MYRES(K,NSBPRC+1,NORD+1)-
     >                             XS0(K),0D0)
                           ELSEIF (MOD(IETYPE,10).EQ.3) THEN
                              IF (DXU(K).GT.ABS(DXL(K))) THEN
                                 DXU(K) = ABS(DXU(K)-
     >                                (MYRES(K,NSBPRC+1,NORD+1) -
     >                                XS0(K)))/2D0
                                 DXL(K) = -DXU(K)
                              ELSE
                                 DXL(K) = -ABS(DXL(K)-
     >                                (MYRES(K,NSBPRC+1,NORD+1) -
     >                                XS0(K)))/2D0
                                 DXU(K) = -DXL(K)
                              ENDIF
                           ENDIF
                           DXU2(K) = DXU2(K) + DXU(K)*DXU(K)
                           DXL2(K) = DXL2(K) + DXL(K)*DXL(K)
                        ENDDO
                     ENDIF

                     IPHASE = 3
                     IF (MOD(IETYPE,10).EQ.2) THEN
                        IPHASE = 4
                     ELSEIF (MOD(IETYPE,10).EQ.3) THEN
                        IPHASE = 5
                     ENDIF
                     CALL UNCERT2(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
                  ENDDO

*--- Separate treatment of each PDF member for deviations
*--- (NNPDF, AB(K)M)
               ELSEIF (MOD(IETYPE,10).EQ.4.OR.
     >                 MOD(IETYPE,10).EQ.5.OR.
     >                 MOD(IETYPE,10).EQ.6) THEN
                  IPHASE = 6
                  DO J=1,MYPDF
                     CALL INITPDF(J)
                     CALL FX9999IN(FILENAME)
                     CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                    XSUNCOR,XSCOR)

*--- Fill member into arrays for crosscheck on PDF uncertainties, part 1 (unnormalized)
                     IF (MOD(IETYPE,10).NE.4) THEN
                        DO K=1,NOBSBIN
                           DXU(K) = MAX(XSNLO(K)-XS0(K),0D0)
                           DXL(K) = MIN(XSNLO(K)-XS0(K),0D0)
                           DXU2(K) = DXU2(K) + DXU(K)*DXU(K)
                           DXL2(K) = DXL2(K) + DXL(K)*DXL(K)
                        ENDDO
                     ELSE
                        DO K=1,NOBSBIN
                           DXL(K)  = DXL(K)  + 1D0
                           DXU2(K) = DXU2(K) + XSNLO(K)
                           DXL2(K) = DXL2(K) + XSNLO(K)*XSNLO(K)
                        ENDDO
                     ENDIF

                     IF (J.EQ.MYPDF) THEN
                        IF (MOD(IETYPE,10).EQ.4) THEN
                           DO K=1,NOBSBIN
                              XS0(K)  = DXU2(K)/DXL(K)
                              DXU2(K) = DXL2(K)/DXL(K) - XS0(K)*XS0(K)
                              DXL2(K) = DXU2(K)
                           ENDDO
                        ELSEIF (MOD(IETYPE,10).EQ.5) THEN
                           DO K=1,NOBSBIN
                              DXU2(K) = (DXU2(K)+DXL2(K))
                              DXL2(K) = DXU2(K)
                           ENDDO
                        ENDIF
                     ENDIF

                     IF (LNRM) THEN
                        ISTEP = 3
                        IF (IDEBUG.GT.0)
     >                       WRITE(*,*)"DEBUG1: DD2 ALLUNC STEP = ",
     >                       ISTEP
                        CALL CENRES2(ISTEP,LRAT,LNRM,
     >                       SCENARIO(1:LEN_TRIM(SCENARIO)))
*--- Load normalization table with potentially different binning!
                        IF (LTAB) THEN
                           CALL FX9999IN(FILENAMN)
                           CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                          XSUNCOR,XSCOR)
                        ENDIF
                        ISTEP = 4
                        IF (IDEBUG.GT.0)
     >                       WRITE(*,*)"DEBUG1: EE2 ALLUNC STEP = ",
     >                       ISTEP
                        CALL CENRES2(ISTEP,LRAT,LNRM,
     >                       SCENARIO(1:LEN_TRIM(SCENARIO)))
                        IF (LTAB) THEN
                           CALL FX9999IN(FILENAME)
                           CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                          XSUNCOR,XSCOR)
                        ENDIF
                        ISTEP = 5
                        IF (IDEBUG.GT.0)
     >                       WRITE(*,*)"DEBUG1: FF2 ALLUNC STEP = ",
     >                       ISTEP
                        CALL CENRES2(ISTEP,LRAT,LNRM,
     >                       SCENARIO(1:LEN_TRIM(SCENARIO)))
                     ENDIF

*--- Fill member into arrays for crosscheck on PDF uncertainties, part 2 (normalized)
*--- Needs to be MYRES here after CENRES STEP 5 to give result for the
*--- actual PDF
                     IF (LNRM) THEN
                        IF (MOD(IETYPE,10).NE.4) THEN
                           DO K=NOBSBIN+1,2*NOBSBIN
                              DXU(K) = MAX(MYRES(K,NSBPRC+1,NORD+1)-
     >                             XS0(K),0D0)
                              DXL(K) = MIN(MYRES(K,NSBPRC+1,NORD+1)-
     >                             XS0(K),0D0)
                              DXU2(K) = DXU2(K) + DXU(K)*DXU(K)
                              DXL2(K) = DXL2(K) + DXL(K)*DXL(K)
                           ENDDO
                        ELSE
                           DO K=NOBSBIN+1,2*NOBSBIN
                              DXL(K)  = DXL(K) +1D0
                              DXU2(K) = DXU2(K)+MYRES(K,NSBPRC+1,NORD+1)
                              DXL2(K) = DXL2(K)+MYRES(K,NSBPRC+1,NORD+1)
     >                             **2D0
                           ENDDO
                        ENDIF
                        IF (J.EQ.MYPDF) THEN
                           IF (MOD(IETYPE,10).EQ.4) THEN
                              DO K=NOBSBIN+1,2*NOBSBIN
                                 XS0(K)  = DXU2(K)/DXL(K)
                                 DXU2(K) = DXL2(K)/DXL(K)-XS0(K)*XS0(K)
                                 DXL2(K) = DXU2(K)
                              ENDDO
                           ELSEIF (MOD(IETYPE,10).EQ.5) THEN
                              DO K=NOBSBIN+1,2*NOBSBIN
                                 DXU2(K) = (DXU2(K)+DXL2(K))
                                 DXL2(K) = DXU2(K)
                              ENDDO
                           ENDIF
                        ENDIF
                     ENDIF

                     CALL UNCERT2(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
                  ENDDO
               ELSE
                  WRITE(*,*)"ALLUNC: ERROR! Undefined PDF error "//
     >                 "method, aborted! IETYPE = ",IETYPE
               ENDIF

            ENDIF

            IPHASE = 9
            CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

*--- Symmetrize uncertainties as required for AB(K)M
            IF (MOD(IETYPE,10).EQ.5) THEN
               IBIN = 0
               DO IRAP=1,NRAP
                  DO IPT=1,NPT(IRAP)
                     IBIN = IBIN+1
                     DO IORD=1,NORD+1
                        DO ISUB=1,NSBPRC+1
Comment:                            write(*,*)"AAAA: ib,is,io,up,lo,",
Comment:      >                          ibin,isub,iord,
Comment:      >                          WTDXU2(IBIN,ISUB,IORD),
Comment:      >                          WTDXL2(IBIN,ISUB,IORD)
*--- Symmetrized upper uncertainty
                           WTDXU2(IBIN,ISUB,IORD) =
     >                          + SQRT(
     >                          WTDXL2(IBIN,ISUB,IORD)**2D0 +
     >                          WTDXU2(IBIN,ISUB,IORD)**2D0 )
*--- Symmetrized lower uncertainty
                           WTDXL2(IBIN,ISUB,IORD) =
     >                          - WTDXU2(IBIN,ISUB,IORD)
Comment:                            write(*,*)"BBBB: ib,is,io,up,lo,",
Comment:      >                          ibin,isub,iord,
Comment:      >                          WTDXU2(IBIN,ISUB,IORD),
Comment:      >                          WTDXL2(IBIN,ISUB,IORD)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

*--- Part 2: Additional PDF uncertainty parts using a separate PDF set
            IF (INT(IETYPE/10).EQ.1) THEN
               IF (PDFNAM(1:LEN_TRIM(PDFNAM)).EQ.
     >              "HERAPDF10_EIG.LHgrid") THEN
                  PDFSET2 = PDFPATH(1:LEN_TRIM(PDFPATH))//
     >                 "/HERAPDF10_VAR.LHgrid"
                  NPDFMOD = 8
                  NPDFPAR = 13
               ELSEIF (PDFNAM(1:LEN_TRIM(PDFNAM)).EQ.
     >                 "HERAPDF15_EIG.LHgrid") THEN
                  PDFSET2 = PDFPATH(1:LEN_TRIM(PDFPATH))//
     >                 "/HERAPDF15_VAR.LHgrid"
                  NPDFMOD = 8
                  NPDFPAR = 12
               ELSEIF (PDFNAM(1:LEN_TRIM(PDFNAM)).EQ.
     >                 "HERAPDF15NLO_EIG.LHgrid") THEN
                  PDFSET2 = PDFPATH(1:LEN_TRIM(PDFPATH))//
     >                 "/HERAPDF15NLO_VAR.LHgrid"
                  NPDFMOD = 8
                  NPDFPAR = 12
               ELSEIF (PDFNAM(1:LEN_TRIM(PDFNAM)).EQ.
     >                 "HERAPDF15NNLO_EIG.LHgrid") THEN
                  PDFSET2 = PDFPATH(1:LEN_TRIM(PDFPATH))//
     >                 "/HERAPDF15NNLO_VAR.LHgrid"
                  NPDFMOD = 8
                  NPDFPAR = 10
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
               CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
               ISTEP = 0
               CALL CENRES2(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))
               IF (LNRM) THEN
*--- Load normalization table with potentially different binning!
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAMN)
                     CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                    XSUNCOR,XSCOR)
                  ENDIF
                  ISTEP = 1
                  CALL CENRES2(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAME)
                     CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                    XSUNCOR,XSCOR)
                  ENDIF
               ENDIF
               ISTEP = 2
               CALL CENRES2(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))

               CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
               IPHASE = 6

*--- HERAPDF1.0 & 1.5: Do loop runs from 1 - 8 (NPDFMOD) for this part
               DO J=1,NPDFMOD
                  CALL INITPDF(J)
                  CALL FX9999IN(FILENAME)
                  CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                 XSUNCOR,XSCOR)
                  IF (LNRM) THEN
                     ISTEP = 3
                     CALL CENRES2(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                     IF (LTAB) THEN
                        CALL FX9999IN(FILENAMN)
                        CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                       XSUNCOR,XSCOR)
                     ENDIF
                     ISTEP = 4
                     CALL CENRES2(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                     IF (LTAB) THEN
                        CALL FX9999IN(FILENAME)
                        CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                       XSUNCOR,XSCOR)
                     ENDIF
                     ISTEP = 5
                     CALL CENRES2(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                  ENDIF
                  CALL UNCERT2(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
               ENDDO

               IPHASE = 9
               CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

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
               CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >              XSUNCOR,XSCOR)
               ISTEP = 0
               CALL CENRES2(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))
               IF (LNRM) THEN
*--- Load normalization table with potentially different binning!
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAMN)
                     CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                    XSUNCOR,XSCOR)
                  ENDIF
                  ISTEP = 1
                  CALL CENRES2(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAME)
                     CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                    XSUNCOR,XSCOR)
                  ENDIF
               ENDIF
               ISTEP = 2
               CALL CENRES2(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))

               CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
               IPHASE = 6

*--- HERAPDF1.0:               Do loop runs from 9 - 13 for this part (NPDFMOD+1,NPDFPAR)
*--- HERAPDF1.5,HERAPDF1.5NLO: Do loop runs from 9 - 12 for this part
*--- HERAPDF1.5NNLO:           Do loop runs from 9 - 10 for this part
               DO J=NPDFMOD+1,NPDFPAR
                  CALL INITPDF(J)
                  CALL FX9999IN(FILENAME)
                  CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                 XSUNCOR,XSCOR)
                  IF (LNRM) THEN
                     ISTEP = 3
                     CALL CENRES2(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                     IF (LTAB) THEN
                        CALL FX9999IN(FILENAMN)
                        CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                       XSUNCOR,XSCOR)
                     ENDIF
                     ISTEP = 4
                     CALL CENRES2(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                     IF (LTAB) THEN
                        CALL FX9999IN(FILENAME)
                        CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                       XSUNCOR,XSCOR)
                     ENDIF
                     ISTEP = 5
                     CALL CENRES2(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LEN_TRIM(SCENARIO)))
                  ENDIF
                  CALL UNCERT2(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
               ENDDO

               IPHASE = 9
               CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

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
                  IF (NORD.EQ.1) THEN
                     WRITE(*,900) IBIN,MYRESN(IBIN,NSBPRC+1,NORD+1),
     >                    WTDXL2(IBIN,NSBPRC+1,NORD+1),
     >                    WTDXU2(IBIN,NSBPRC+1,NORD+1)
                  ELSEIF (NORD.EQ.2) THEN
                     WRITE(*,910) IBIN,MYRESN(IBIN,NSBPRC+1,NORD+1),
     >                    WTDXL2(IBIN,NSBPRC+1,NORD+1),
     >                    WTDXU2(IBIN,NSBPRC+1,NORD+1),
     >                    MYRESN(IBIN,NSBPRC+1,NORD+1) /
     >                    (MYRESN(IBIN,NSBPRC+1,1) +
     >                    MYRESN(IBIN,NSBPRC+1,2))
                  ELSEIF (NORD.GE.3) THEN
                     WRITE(*,920) IBIN,MYRESN(IBIN,NSBPRC+1,NORD+1),
     >                    WTDXL2(IBIN,NSBPRC+1,NORD+1),
     >                    WTDXU2(IBIN,NSBPRC+1,NORD+1),
     >                    (MYRESN(IBIN,NSBPRC+1,1) +
     >                    MYRESN(IBIN,NSBPRC+1,2)) /
     >                    MYRESN(IBIN,NSBPRC+1,1),
     >                    MYRESN(IBIN,NSBPRC+1,NORD+1) /
     >                    (MYRESN(IBIN,NSBPRC+1,1) +
     >                    MYRESN(IBIN,NSBPRC+1,2))
                  ENDIF
c - Debug: Output should be identical in single PDF set case!
                  IF (INT(IETYPE/10).EQ.0) THEN
                     DXS0L = -SQRT(DXL2(IBIN))/XS0(IBIN)
                     DXS0U = +SQRT(DXU2(IBIN))/XS0(IBIN)
                     DELTA1 = ABS(MYRESN(IBIN,NSBPRC+1,NORD+1) -
     >                    XS0(IBIN)) / XS0(IBIN)
                     DELTA2 = ABS(WTDXL2(IBIN,NSBPRC+1,NORD+1) -
     >                    DXS0L) / (-DXS0L)
                     DELTA3 = ABS(WTDXU2(IBIN,NSBPRC+1,NORD+1) -
     >                    DXS0U) / DXS0U
                     IF (DELTA1.GT.1D-10.OR.
     >                    DELTA2.GT.1D-10.OR.
     >                    DELTA3.GT.1D-10) THEN
                        WRITE(*,*)"ALLUNC: WARNING! Potential mismatch"/
     >                       /" in PDF uncertainty calculation!",
     >                       DELTA1,DELTA2,DELTA3
                        WRITE(*,900) IBIN,XS0(IBIN),DXS0L,DXS0U
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO

c - Fill histograms
Comment:             CALL PDFFILL(NRAP,0,-1,I,MYRESN)
Comment:             CALL PDFFILL(NRAP,1,-1,I,WTDXL2)
Comment:             CALL PDFFILL(NRAP,2,-1,I,WTDXU2)
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
      CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)



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
ckr Account for different event numbers in the merged tables
            IWEIGHT = 1
ckr            IWEIGHT = 0
c - This is the normal table to fill the default values
            CALL FX9999IN(FILENAME)
            CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
            CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
            CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
            IPHASE = 6

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
                     FILENAMES = FILEBASE(1:LEN_TRIM(FILEBASE))//
     >                    "2jet_v22_"
     >                    //NO//".tab"
                     OPEN(2,STATUS='OLD',FILE=FILENAMES,IOSTAT=ISTAT)
                     IF (ISTAT.NE.0) THEN
C---  WRITE(*,*)"Filename for order",IORD,":",FILENAMES
                     NFAIL = NFAIL + 1
                     IF (NFAIL.LT.2) THEN
                        WRITE(*,*)"ALLUNC: WARNING! Table file "//
     >                       "not found, skipped ! Filename = ",
     >                       FILENAMES," IOSTAT = ",ISTAT
                     ENDIF
ckr While using f90 DO-ENDDO ... one could also use the EXIT statement
ckr instead of GOTO. However, EXIT leaves the loop!
ckr To continue the loop use CYCLE instead!
                     GOTO 20
ckr            CYCLE
                  ENDIF
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
               CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
               CALL UNCERT2(IPHASE,IMODE,IWEIGHT,ITAB,LRAT,LNRM)
 20            CONTINUE
            ENDDO

            WRITE(*,*)"########################################"//
     >           "################################"
            WRITE(*,*)"No. of filenames skipped:",NFAIL
            WRITE(*,*)"No. of filenames analyzed:",NFOUND
            WRITE(*,*)"########################################"//
     >           "################################"

            IPHASE = 9
            CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

c - Special statistics output
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
                  DSTMP(IBIN,3) = -1D0
                  DSTMP(IBIN,4) = -1D0
                  IF (MYRES(IBIN,NSBPRC+1,IORD2).GE.0.D0) THEN
                     DSTMP(IBIN,1) = (WTXMIN(IBIN,NSBPRC+1,IORD2)-1D0)
     >                    *100D0
                     DSTMP(IBIN,2) = (WTXMAX(IBIN,NSBPRC+1,IORD2)-1D0)
     >                    *100D0
                     IF (WTDXU2(IBIN,NSBPRC+1,IORD2).GT.1D-99) THEN
                        DSTMP(IBIN,3) =
     >                       (WTXMIN(IBIN,NSBPRC+1,IORD2)-1D0) /
     >                       WTDXU2(IBIN,NSBPRC+1,IORD2)
                        DSTMP(IBIN,4) =
     >                       (WTXMAX(IBIN,NSBPRC+1,IORD2)-1D0) /
     >                       WTDXU2(IBIN,NSBPRC+1,IORD2)
                     ENDIF
                  ELSE
                     DSTMP(IBIN,1) = (WTXMIN(IBIN,NSBPRC+1,IORD2)+1D0)
     >                    *100D0
                     DSTMP(IBIN,2) = (WTXMAX(IBIN,NSBPRC+1,IORD2)+1D0)
     >                    *100D0
                     IF (WTDXU2(IBIN,NSBPRC+1,IORD2).GT.1D-99) THEN
                        DSTMP(IBIN,3) =
     >                       (WTXMIN(IBIN,NSBPRC+1,IORD2)+1D0) /
     >                       WTDXU2(IBIN,NSBPRC+1,IORD2)
                        DSTMP(IBIN,4) =
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
     >                 DSTMP(IBIN,1),DSTMP(IBIN,2),
     >                 DSTMP(IBIN,3),DSTMP(IBIN,4)
               ENDDO
            ENDDO

c - Write out list of files with largest fluctuations
* ---
* Could replace the DO 10 CONTINUE and GOTO with named DO-loop and CYCLE
* here, but then the fixed-format source code doesn't work with emacs
* f90-mode indentation :-(
* ---
            IOPEN = 0
            KILLFILE = SCENARIO(1:LEN_TRIM(SCENARIO))//"_kill.txt"
            IF (IORD.EQ.0) THEN
               DO 10 ITAB=0,NTAB
                  WRITE(NO,'(I4.4)'),ITAB
                  FILENAMES = FILEBASE(1:LEN_TRIM(FILEBASE))//"2jet_"/
     >                 /NO//".tab"
                  OPEN(2,STATUS='OLD',FILE=FILENAMES,IOSTAT=ISTAT)
                  IF (ISTAT.NE.0) THEN
                     FILENAMES = FILEBASE(1:LEN_TRIM(FILEBASE))//"3jet_"
     >                    //NO//".tab"
                  ENDIF
                  CLOSE(2)
                  IBIN = 0
                  DO IRAP=1,NRAP
                     DO IPT=1,NPT(IRAP)
                        IBIN = IBIN+1
                        IF ((ITAB.EQ.IJMIN(IBIN).AND.
     >                       DSTMP(IBIN,3).LE.-10D0).OR.
     >                       (ITAB.EQ.IJMAX(IBIN).AND.
     >                       DSTMP(IBIN,4).GE.10D0)) THEN
                           IF (IOPEN.EQ.0) THEN
                              OPEN(99,
     >                             FILE=KILLFILE(1:LEN_TRIM(KILLFILE)),
     >                             IOSTAT=ISTAT)
                              IF (ISTAT.EQ.0) THEN
                                 IOPEN = 1
                              ELSE
                                 WRITE(*,*)"ALLUNC: ERROR! Could not "//
     >                                "open killfile, aborted!"
                                 STOP
                              ENDIF
                           ENDIF
                           WRITE(99,'(A)')
     >                          FILENAMES(1:LEN_TRIM(FILENAMES))
                           GOTO 10
                           WRITE(99,'(A)')
     >                          FILENAMES(1:LEN_TRIM(FILENAMES))
                           GOTO 10
                        ENDIF
                     ENDDO
                  ENDDO
 10            CONTINUE
            ENDIF
            IF (IOPEN.EQ.1) CLOSE(99)
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
Comment:             CALL PDFFILL(NRAP,3,IORD,ISCL,WTDXMN)
Comment:             CALL PDFFILL(NRAP,4,IORD,ISCL,WTDXUL)

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
Comment:             CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:             ISTEP = 0
Comment:             CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             IF (LNRM) THEN
Comment: ckr Load normalization table with potentially different binning!
Comment:                IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:                ISTEP = 1
Comment:                CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:                IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:             ENDIF
Comment:             ISTEP = 2
Comment:             CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:
Comment:             CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
Comment:             IPHASE = 6
Comment:
Comment:             CALL INITPDF(IPDF-IPDFUD)
Comment: cnew
Comment:             ASMZTMP = ALPHASPDF(ZMASS)
Comment:             ASDN = ASMZTMP - ASMZPDF
Comment: cnew
Comment:             CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:             IF (LNRM) THEN
Comment:                ISTEP = 3
Comment:                CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:                IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:                ISTEP = 4
Comment:                CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:                IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:                ISTEP = 5
Comment:                CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             ENDIF
Comment:             CALL UNCERT2(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
Comment:
Comment:             CALL INITPDF(IPDF+IPDFUD)
Comment: cnew
Comment:             ASMZTMP = ALPHASPDF(ZMASS)
Comment:             ASUP = ASMZTMP - ASMZPDF
Comment: cnew
Comment:             CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:             IF (LNRM) THEN
Comment:                ISTEP = 3
Comment:                CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:                IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:                ISTEP = 4
Comment:                CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:                IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:                ISTEP = 5
Comment:                CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >              SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             ENDIF
Comment:             CALL UNCERT2(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
Comment:
Comment:             IPHASE = 9
Comment:             CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
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
Comment:          CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:          ISTEP = 0
Comment:          CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:          IF (LNRM) THEN
Comment: ckr Load normalization table with potentially different binning!
Comment:             IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:             ISTEP = 1
Comment:             CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:          ENDIF
Comment:          ISTEP = 2
Comment:          CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:
Comment:          CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
Comment:          IPHASE = 6
Comment:
Comment: cnew         ASMZVAL = ASMZVAL - DBLE(IPDFUD)/1000.
Comment:          ASMZVAL = ASMZVAL + ASDN
Comment:          CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:          IF (LNRM) THEN
Comment:             ISTEP = 3
Comment:             CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >           SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:             ISTEP = 4
Comment:             CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >           SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:             ISTEP = 5
Comment:             CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >           SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:          ENDIF
Comment:          CALL UNCERT2(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
Comment:
Comment: ckr Add twice to correct previous subtraction!
Comment: cnew         ASMZVAL = ASMZVAL + DBLE(2*IPDFUD)/1000.
Comment:          ASMZVAL = ASMZVAL - ASDN + ASUP
Comment:          CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:          IF (LNRM) THEN
Comment:             ISTEP = 3
Comment:             CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >           SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:             ISTEP = 4
Comment:             CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >           SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:             IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:             ISTEP = 5
Comment:             CALL CENRES2(ISTEP,LRAT,LNRM,
Comment:      >           SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:          ENDIF
Comment:          CALL UNCERT2(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
Comment:
Comment:          IPHASE = 9
Comment:          CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
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
            CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
            ISTEP = 0
            CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
            IF (LNRM) THEN
ckr Load normalization table with potentially different binning!
               IF (LTAB) THEN
                  CALL FX9999IN(FILENAMN)
                  CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
               ENDIF
               ISTEP = 1
               CALL CENRES2(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))
               IF (LTAB) THEN
                  CALL FX9999IN(FILENAME)
                  CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
               ENDIF
            ENDIF
            ISTEP = 2
            CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))

            CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
            IPHASE = 6

            DO I=2,3
               ISCL = ISCLPT(I)
               XMUR = XMURS(I)
               XMUF = XMUFS(I)
               CALL FX9999IN(FILENAME)
               CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
               IF (LNRM) THEN
                  ISTEP = 3
                  CALL CENRES2(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
ckr Load normalization table with potentially different binning!
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAMN)
                     CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                    XSUNCOR,XSCOR)
                  ENDIF
                  ISTEP = 4
                  CALL CENRES2(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAME)
                     CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                    XSUNCOR,XSCOR)
                  ENDIF
                  ISTEP = 5
                  CALL CENRES2(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
               ENDIF
               CALL UNCERT2(IPHASE,IMODE,IWEIGHT,ISCL,LRAT,LNRM)
            ENDDO
            IPHASE = 9
            CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

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
cnew            DO IRAP=1,NRAP
cnew               DO IPT=1,NPT(IRAP)
cnew                  IBIN = IBIN+1
            DO IBIN=1,NOBSBIN
                  WRITE(*,900) IBIN,MYRESN(IBIN,NSBPRC+1,NORD+1),
     >                 WTDXLM(IBIN,NSBPRC+1,NORD+1),
     >                 WTDXUM(IBIN,NSBPRC+1,NORD+1)
cnew               ENDDO
            ENDDO

*---  Central scale
            ISCL = ISCLPT(1)
Comment:             CALL PDFFILL(NRAP,6,-1,ISCL,WTDXLM)
Comment:             CALL PDFFILL(NRAP,7,-1,ISCL,WTDXUM)
         ENDIF

c - 6-point scheme, only if possible, e.g. not in case of threshold corrections
         IF (NSCLS.GE.7.AND.NORD.LT.3) THEN
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
            CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
            ISTEP = 0
            CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
            IF (LNRM) THEN
ckr Load normalization table with potentially different binning!
               IF (LTAB) THEN
                  CALL FX9999IN(FILENAMN)
                  CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
               ENDIF
               ISTEP = 1
               CALL CENRES2(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LEN_TRIM(SCENARIO)))
               IF (LTAB) THEN
                  CALL FX9999IN(FILENAME)
                  CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
               ENDIF
            ENDIF
            ISTEP = 2
            CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))

            CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
            IPHASE = 6

            DO I=2,7
               ISCL = ISCLPT(I)
               XMUR = XMURS(I)
               XMUF = XMUFS(I)
               CALL FX9999IN(FILENAME)
               CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,XSUNCOR,XSCOR)
               IF (LNRM) THEN
                  ISTEP = 3
                  CALL CENRES2(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
ckr Load normalization table with potentially different binning!
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAMN)
                     CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                    XSUNCOR,XSCOR)
                  ENDIF
                  ISTEP = 4
                  CALL CENRES2(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
                  IF (LTAB) THEN
                     CALL FX9999IN(FILENAME)
                     CALL FX9999CC(XMUR,XMUF,XSNLO,XSCLNLO,
     >                    XSUNCOR,XSCOR)
                  ENDIF
                  ISTEP = 5
                  CALL CENRES2(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LEN_TRIM(SCENARIO)))
               ENDIF
               CALL UNCERT2(IPHASE,IMODE,IWEIGHT,ISCL,LRAT,LNRM)
            ENDDO
            IPHASE = 9
            CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

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
Comment:             CALL PDFFILL(NRAP,8,-1,ISCL,WTDXLM)
Comment:             CALL PDFFILL(NRAP,9,-1,ISCL,WTDXUM)
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
Comment: ckr         CALL FX9999CC(FILENAME,XMUR,XMUF,1,XSNLO,XSCLNLO)
Comment:          CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:          IPHASE  = 1
Comment:          IMODE   = 4
Comment:          IWEIGHT = 0
Comment:          CALL CENRES2(ISTEP,LRAT,LNRM,SCENARIO(1:LEN_TRIM(SCENARIO)))
Comment:          CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
Comment:
Comment: c - Normal result
Comment:          ISCL = 3
Comment:          XMUR = MURSCALE(ISCL)
Comment:          XMUF = MUFSCALE(ISCL)
Comment:          CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSNLO,XSCLNLO)
Comment:          IPHASE = 6
Comment:          CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
Comment:          IPHASE = 9
Comment:          CALL UNCERT2(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
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
Comment:       CALL PDFHIST(2,HISTFILE,LONE,LPDF,LSTAT,LALG,LSER,MYPDF,
Comment:      >     LRAT.OR.LNRM,LSCL,ISCLPT(1))

      RETURN
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
      INCLUDE "v14unc.inc"
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
      INCLUDE "v14unc.inc"
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
