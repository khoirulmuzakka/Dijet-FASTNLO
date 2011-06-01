      PROGRAM ALLUNC
* ---------------------------------------------------------------------
* K. Rabbertz 07.09.2008 First try to integrate all uncertainties
*                        into one job
*
* ALLUNC - Program to derive all (PDF, statistical, algorithmic,
*          scale ...) uncertainties using fastNLO tables
*
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
      INTEGER BORNN,NLON,LENOCC
      INTEGER I,J,MYPDF,IPDF,IPDFUD,MYPDFAS,IOPDF,IOAS,NSCLS
      INTEGER ITAB,NTAB,NFOUND,NFAIL
      INTEGER ISTAT,ISCL,IORD,IORD2,IBIN,NBIN,ISUB,IRAP,IPT
      INTEGER IHIST,IPHASE,ISTEP,IETYPE
      LOGICAL LONE,LPDF,LSTAT,LSER,LSCL,LRAT,LALG,LNRM,LTAB
cnew
      DOUBLE PRECISION ALPHASPDF,ASMZPDF,ASUP,ASDN
      DOUBLE PRECISION XMUR,XMUF,QLAM4,QLAM5,BWGT
      DOUBLE PRECISION DSTMP(4)
      DOUBLE PRECISION WTXTMP(NBINTOTMAX,NMAXSUBPROC+1,NMAXORD+1)
      DOUBLE PRECISION WTXLTMP(NBINTOTMAX,NMAXSUBPROC+1,NMAXORD+1)
      DOUBLE PRECISION WTXUTMP(NBINTOTMAX,NMAXSUBPROC+1,NMAXORD+1)
c - To unify quoted uncertainties (CL68,CL90,special)
c - Convert from CL68 to CL90 values
c - TOCL90 = 1.64485D0 ! SQRT(2.D0)/InvERF(0.9D0)
c - Convert from GJR to CTEQ (CL90) values ?
c - TOCL90GJR = 2.12766D0! 1.D0/0.47D0
      DOUBLE PRECISION TOCL90,TOCL90GJR
      PARAMETER (TOCL90 = 1.64485D0, TOCL90GJR = 2.12766D0)
      
c - Attention!!! This must be declared consistent with the
c                definition in the commonblock!!!!!
      DOUBLE PRECISION XSECT0(NBINTOTMAX,3)
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
      LSCL = .TRUE.
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
         IF (SCENARIO(1:LENOCC(SCENARIO)).EQ."-h") THEN
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
         ELSEIF (SCENARIO(1:7).EQ."fnl2622".OR.
     >           SCENARIO(1:7).EQ."fnl2652") THEN
            LNRM = .TRUE.
            LTAB = .FALSE.
            WRITE(*,*)
     >           "ALLUNC: Deriving normalized distributions"
         ELSEIF (SCENARIO(1:10).EQ."fnl2722num") THEN
            LNRM = .TRUE.
            LTAB = .TRUE.
            WRITE(*,*)
     >           "ALLUNC: Deriving normalized distributions"
         ENDIF
         WRITE(*,*)"ALLUNC: Evaluating scenario: ",
     >        SCENARIO(1:LENOCC(SCENARIO))
      ENDIF
ckr      lnrm = .false.
      TABNAME = SCENARIO(1:LENOCC(SCENARIO))//".tab"
      IF (LNRM.AND.LTAB) THEN
         TABNAMN = SCENARIO(1:7)//"norm"//".tab"
      ENDIF
      REFNAME = SCENARIO(1:LENOCC(SCENARIO))//"ref.tab"

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
     >     TABPATH(1:LENOCC(TABPATH))
      FILENAME = TABPATH(1:LENOCC(TABPATH))//"/"//
     >     TABNAME(1:LENOCC(TABNAME))
      WRITE(*,*)"ALLUNC: Taking primary table ",
     >     FILENAME(1:LENOCC(FILENAME))
      IF (LNRM.AND.LTAB) THEN
         FILENAMN = TABPATH(1:LENOCC(TABPATH))//"/"//
     >        TABNAMN(1:LENOCC(TABNAMN))
         WRITE(*,*)"ALLUNC: Taking normalization table ",
     >        FILENAMN(1:LENOCC(FILENAMN))
      ENDIF

*---HBOOK filename
      HISTFILE = "X"
      IF (IARGC().GE.3) THEN
         CALL GETARG(3,HISTFILE)
      ENDIF
      IF (IARGC().LT.3.OR.HISTFILE(1:1).EQ."_") THEN
         HISTFILE = SCENARIO(1:LENOCC(SCENARIO))//".hbk"
         WRITE(*,*)
     >        "ALLUNC: WARNING! No output filename given, "//
     >        "taking scenario.hbk instead!"
      ELSE
         WRITE(*,*)"ALLUNC: Creating output file: ",
     >        HISTFILE(1:LENOCC(HISTFILE))
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

*---PDF set
      PDFSET = "X"
      IF (IARGC().GE.7) THEN
         CALL GETARG(7,PDFSET)
      ENDIF
      IF (IARGC().LT.7.OR.PDFSET(1:1).EQ."_") THEN
         PDFSET = "cteq66.LHgrid"
         WRITE(*,*)
     >        "ALLUNC: WARNING! No PDF set given, "//
     >        "taking cteq66.LHgrid instead!"
      ELSE
         WRITE(*,*)"ALLUNC: Using PDF set: ",
     >        PDFSET(1:LENOCC(PDFSET))
      ENDIF
      PDFNAM = PDFSET(1:LENOCC(PDFSET))

*---Path to PDF sets
      CHTMP = "X"
      IF (IARGC().GE.8) THEN
         CALL GETARG(8,CHTMP)
      ENDIF
      IF (IARGC().LT.8.OR.CHTMP(1:1).EQ."_") THEN
         PDFPATH = "/../share/lhapdf/PDFsets"
         WRITE(*,*)
     >        "ALLUNC: No PDF path given, "//
     >        "assuming: $(LHAPDF)"//PDFPATH(1:LENOCC(PDFPATH))
c - Initialize path to LHAPDF libs
         CALL GETENV("LHAPDF",LHAPDF)
         IF (LENOCC(LHAPDF).EQ.0) THEN
            WRITE(*,*)"\nALLUNC: ERROR! $LHAPDF not set, aborting!"
            STOP
         ENDIF
         PDFPATH = LHAPDF(1:LENOCC(LHAPDF))//
     >        PDFPATH(1:LENOCC(PDFPATH))
      ELSE
         PDFPATH = CHTMP(1:LENOCC(CHTMP)) 
      ENDIF
      WRITE(*,*)"ALLUNC: Looking for LHAPDF PDF sets in path: "//
     >     PDFPATH(1:LENOCC(PDFPATH))
      PDFSET = PDFPATH(1:LENOCC(PDFPATH))//"/"//PDFSET
      WRITE(*,*)"ALLUNC: Taking PDF set "
     >     //PDFSET(1:LENOCC(PDFSET))

*---alpha_s mode
      CHTMP = "X"
      IF (IARGC().GE.9) THEN
         CALL GETARG(9,CHTMP)
      ENDIF
      IF (IARGC().LT.9.OR.CHTMP(1:1).EQ."_") THEN
         ASMODE = "PDF"
         WRITE(*,*)
     >        "ALLUNC: No alpha_s mode given, "//
     >        "using alpha_s according to PDF set"
      ELSE
         ASMODE = CHTMP(1:LENOCC(CHTMP))
         WRITE(*,*)"ALLUNC: Using alpha_s mode: ",
     >        ASMODE(1:LENOCC(ASMODE))
      ENDIF

*---alpha_s(M_Z)
      CH8TMP = "X"
      IF (IARGC().GE.10) THEN
         CALL GETARG(10,CH8TMP)
      ENDIF
      IF (IARGC().LT.10.OR.CH8TMP(1:1).EQ."_") THEN
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
      IF (IARGC().GE.11) THEN
         CALL GETARG(11,CH8TMP)
      ENDIF
      IF (IARGC().LT.11.OR.CH8TMP(1:1).EQ."_") THEN
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
      IF (IARGC().GE.12) THEN
         CALL GETARG(12,CH4TMP)
      ENDIF
      IF (IARGC().LT.12.OR.CH4TMP(1:1).EQ."_") THEN
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
      IF (IARGC().GE.13) THEN
         CALL GETARG(13,CH4TMP)
      ENDIF
      IF (IARGC().LT.13.OR.CH4TMP(1:1).EQ."_".OR.
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
      IF (IARGC().GT.13) THEN
         WRITE(*,*)"\nALLUNC: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF



c - Initialize LHAPDF, no PDF set printout after first call
      CALL INITPDFSET(PDFSET(1:LENOCC(PDFSET)))
      CALL SETLHAPARM('SILENT')

c - Initialize one member, 0=best fit member
      CALL INITPDF(0)

c - Write out some info on best fit member      
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
      
c - Check primary table existence
      FILENAME = TABPATH(1:LENOCC(TABPATH))//"/"//TABNAME
      WRITE(*,*)"ALLUNC: Checking primary table: "//
     >     FILENAME(1:LENOCC(FILENAME))
      OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=ISTAT)
      IF (ISTAT.NE.0) THEN
         WRITE(*,*)"ALLUNC: ERROR! Primary table not found, "//
     >        "aborting! IOSTAT = ",ISTAT
         STOP
      ELSE
         CLOSE(2)
      ENDIF

c - Check uncertainties to derive 
      LONE  = MYPDF.LE.1
      LSTAT = BORNN.GE.2.OR.NLON.GE.2
      IF (LALG) THEN
         FILENAME = TABPATH(1:LENOCC(TABPATH))//"/"//REFNAME
         WRITE(*,*)"ALLUNC: Checking reference table: "//
     >        FILENAME(1:LENOCC(FILENAME))
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
      

      
c - One initial call - to fill commonblock -> for histo-booking
c - Use primary table for this (recall: ref. table has 2 x rap. bins)
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      WRITE(*,*)"ALLUNC: Initialize Table, Book Histograms"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      FILENAME = TABPATH(1:LENOCC(TABPATH))//"/"//TABNAME
ckr      CALL FX9999CC(FILENAME,1D0,1D0,1,XSECT0)
      CALL FX9999CC(FILENAME,1D0,1D0,0,XSECT0)
      CALL PDFHIST(1,HISTFILE,LONE,LPDF,LSTAT,LALG,LSER,MYPDF,
     >     LRAT.OR.LNRM,LSCL)
      WRITE(*,*)"ALLUNC: The observable has",NBINTOT," bins -",
     >     NSUBPROC," subprocesses"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
      


c - Define output formats
ckr 900     FORMAT(1P,I5,3(3X,E21.14))
 900        FORMAT(1P,I5,3(6X,E18.11))
 901        FORMAT(1P,I5,2(6X,E18.11))
 902        FORMAT(3I6,3E16.5,5(F10.3,3X))
      

      
c - PDF part
c - Use primary table
c - New call: a single call for each scale
c         1st argument:  name of table
c         2nd argument:  xmur  prefactor for nominal ren-scale
c                              any choice is possible, but please note 
c                              that NNLO-NLL works only for xmur=xmuf
c         3rd argument:  xmuf  prefactor for nominal fact-scale
c                              only a few choices are possible
c                              (see output or table documentation)
c         4th argument:  0: no ascii output       1: print results
c         5th argument:  array to return results

c - Compute central result and, if possible, PDF uncertainties for
c - all precalculated scale variations
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
         DO I=1,NSCALEVAR
ckr         DO I=1,1
ckr Part 1: Basic PDF uncertainty using a single PDF set (and only one
ckr         initialization call!). For PDF uncertainty calculation from
ckr         multiple PDF sets like full HERAPDF make sure to switch back
ckr         to original one here
            IF (IETYPE.EQ.3) THEN
               CALL INITPDFSET(PDFSET(1:LENOCC(PDFSET)))
            ENDIF
            XMUR = MURSCALE(I)
            XMUF = MUFSCALE(I)
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
            CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
            ISTEP = 0
ckr            WRITE(*,*)"AAAAA: ALLUNC STEP = ",ISTEP
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
            IF (LNRM) THEN
ckr Load normalization table with potentially different binning!
               IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
               ISTEP = 1
ckr               WRITE(*,*)"BBBBB: ALLUNC STEP = ",ISTEP
               CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
               IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
            ENDIF
            ISTEP = 2
ckr            WRITE(*,*)"CCCCC: ALLUNC STEP = ",ISTEP
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
            
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
            IPHASE = 2

ckr Do loop runs once even if MYPDF=0! => Avoid with IF statement
            IF (LPDF) THEN
               DO J=1,MYPDF
ckr               DO J=1,3
                  CALL INITPDF(J)
                  CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
                  IF (LNRM) THEN
                     ISTEP = 3
ckr                     WRITE(*,*)"DDDDD: ALLUNC STEP = ",ISTEP
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LENOCC(SCENARIO)))
ckr Load normalization table with potentially different binning!
                     IF (LTAB)
     >                    CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
                     ISTEP = 4
ckr                     WRITE(*,*)"EEEEE: ALLUNC STEP = ",ISTEP
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LENOCC(SCENARIO)))
                     IF (LTAB)
     >                    CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
                     ISTEP = 5
ckr                     WRITE(*,*)"FFFFF: ALLUNC STEP = ",ISTEP
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LENOCC(SCENARIO)))
                  ENDIF
                  CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
               ENDDO
            ENDIF
            
            IPHASE = 3
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

ckr Part 2: Additional PDF uncertainty parts using a separate PDF set
            IF (IETYPE.EQ.3) THEN
               IF (PDFNAM(1:LENOCC(PDFNAM)).EQ.
     >              "HERAPDF10_EIG.LHgrid") THEN
                  PDFSET2 = PDFPATH(1:LENOCC(PDFPATH))//
     >                 "/HERAPDF10_VAR.LHgrid"
ckr                  WRITE(*,*)"ALLUNC: Taking second HERAPDF set: "//
ckr     >                 PDFSET2(1:LENOCC(PDFSET2))
               ELSE
                  WRITE(*,*)"ALLUNC: Illegal HERAPDF set, aborted! "//
     >                 "PDFNAM: ",PDFNAM(1:LENOCC(PDFNAM))
                  STOP
               ENDIF
ckr Back up central result and relevant uncertainty
               IBIN = 0
               DO IRAP=1,NRAP
                  DO IPT=1,NPT(IRAP)
                     IBIN = IBIN+1
                     DO IORD=1,NORD+1
                        DO ISUB=1,NSUBPROC+1
ckr Central result
                           WTXTMP(IBIN,ISUB,IORD) =
     >                          MYRESN(IBIN,ISUB,IORD)
ckr Rel. lower uncertainty
                           WTXLTMP(IBIN,ISUB,IORD) =
     >                          WTDXL2(IBIN,ISUB,IORD)
ckr Rel. upper uncertainty
                           WTXUTMP(IBIN,ISUB,IORD) =
     >                          WTDXU2(IBIN,ISUB,IORD)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO

ckr HERAPDF1.0: 2nd PDF set for additional PDF uncertainty
ckr             with quadratic addition of
ckr             all lower/upper deviations 
               CALL INITPDFSET(PDFSET2(1:LENOCC(PDFSET2)))
               CALL INITPDF(0)
               IPHASE  = 1
               IMODE   = 1
               IWEIGHT = 0
               NRAP  = NRAPIDITY
               IF (LRAT.OR.LNRM) THEN
                  NRAP = 2*NRAPIDITY
               ENDIF
               CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
               ISTEP = 0
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
               IF (LNRM) THEN
ckr Load normalization table with potentially different binning!
                  IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,
     >                 XSECT0)
                  ISTEP = 1
                  CALL CENRES(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LENOCC(SCENARIO)))
                  IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,
     >                 XSECT0)
               ENDIF
               ISTEP = 2
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
               
               CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
               IPHASE = 2
               
ckr HERAPDF1.0: Do loop runs from 1 - 8 for this part
               DO J=1,8
                  CALL INITPDF(J)
                  CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
                  IF (LNRM) THEN
                     ISTEP = 3
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LENOCC(SCENARIO)))
                     IF (LTAB)
     >                    CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
                     ISTEP = 4
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LENOCC(SCENARIO)))
                     IF (LTAB)
     >                    CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
                     ISTEP = 5
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LENOCC(SCENARIO)))
                  ENDIF
                  CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
               ENDDO
               
               IPHASE = 3
               CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

ckr Add quadratically to previously backed-up result
               IBIN = 0
               DO IRAP=1,NRAP
                  DO IPT=1,NPT(IRAP)
                     IBIN = IBIN+1
                     DO IORD=1,NORD+1
                        DO ISUB=1,NSUBPROC+1
ckr Central result
ckr Rel. lower uncertainty
                           WTXLTMP(IBIN,ISUB,IORD) =
     >                          -SQRT(WTXLTMP(IBIN,ISUB,IORD)**2.+
     >                          WTDXL2(IBIN,ISUB,IORD)**2.)
ckr Rel. upper uncertainty
                           WTXUTMP(IBIN,ISUB,IORD) =
     >                          +SQRT(WTXUTMP(IBIN,ISUB,IORD)**2.+
     >                          WTDXU2(IBIN,ISUB,IORD)**2.)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO

ckr HERAPDF1.0: 2nd PDF set for additional PDF uncertainty
ckr             with quadratic addition of
ckr             minimal/maximal lower/upper deviations 
               CALL INITPDF(0)
               IPHASE  = 1
               IMODE   = 3
               IWEIGHT = 0
               NRAP  = NRAPIDITY
               IF (LRAT.OR.LNRM) THEN
                  NRAP = 2*NRAPIDITY
               ENDIF
               CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
               ISTEP = 0
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
               IF (LNRM) THEN
ckr Load normalization table with potentially different binning!
                  IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,
     >                 XSECT0)
                  ISTEP = 1
                  CALL CENRES(ISTEP,LRAT,LNRM,
     >                 SCENARIO(1:LENOCC(SCENARIO)))
                  IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,
     >                 XSECT0)
               ENDIF
               ISTEP = 2
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
               
               CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
               IPHASE = 2

ckr HERAPDF1.0: Do loop runs from 9 - 13 for this part
               DO J=9,13
                  CALL INITPDF(J)
                  CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
                  IF (LNRM) THEN
                     ISTEP = 3
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LENOCC(SCENARIO)))
                     IF (LTAB)
     >                    CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
                     ISTEP = 4
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LENOCC(SCENARIO)))
                     IF (LTAB)
     >                    CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
                     ISTEP = 5
                     CALL CENRES(ISTEP,LRAT,LNRM,
     >                    SCENARIO(1:LENOCC(SCENARIO)))
                  ENDIF
                  CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
               ENDDO
               
               IPHASE = 3
               CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
               
ckr Add quadratically to previously backed-up result
               IBIN = 0
               DO IRAP=1,NRAP
                  DO IPT=1,NPT(IRAP)
                     IBIN = IBIN+1
                     DO IORD=1,NORD+1
                        DO ISUB=1,NSUBPROC+1
ckr Central result
ckr Rel. lower uncertainty
                           WTXLTMP(IBIN,ISUB,IORD) =
     >                          -SQRT(WTXLTMP(IBIN,ISUB,IORD)**2.+
     >                          WTDXLM(IBIN,ISUB,IORD)**2.)
ckr Rel. upper uncertainty
                           WTXUTMP(IBIN,ISUB,IORD) =
     >                          +SQRT(WTXUTMP(IBIN,ISUB,IORD)**2.+
     >                          WTDXUM(IBIN,ISUB,IORD)**2.)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO

ckr Copy back to original arrays as for single PDF set use
               IBIN = 0
               DO IRAP=1,NRAP
                  DO IPT=1,NPT(IRAP)
                     IBIN = IBIN+1
                     DO IORD=1,NORD+1
                        DO ISUB=1,NSUBPROC+1
ckr Central result
                           MYRESN(IBIN,ISUB,IORD) =
     >                          WTXTMP(IBIN,ISUB,IORD)
ckr Rel. lower uncertainty
                           WTDXL2(IBIN,ISUB,IORD) =
     >                          WTXLTMP(IBIN,ISUB,IORD)
ckr Rel. upper uncertainty
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
                  WRITE(*,900) IBIN,MYRESN(IBIN,NSUBPROC+1,NORD+1),
     >                 WTDXL2(IBIN,NSUBPROC+1,NORD+1),
     >                 WTDXU2(IBIN,NSUBPROC+1,NORD+1)
               ENDDO
            ENDDO
            
c - Fill histograms
            CALL PDFFILL(NRAP,0,-1,I,MYRESN)
            CALL PDFFILL(NRAP,1,-1,I,WTDXL2)
            CALL PDFFILL(NRAP,2,-1,I,WTDXU2)
         ENDDO                     ! Loop over scales
      ENDIF
c - Make sure to use again the central PDF!
      CALL INITPDFSET(PDFSET(1:LENOCC(PDFSET)))
      ISCL = 3
      XMUR = MURSCALE(ISCL)
      XMUF = MUFSCALE(ISCL)
      CALL INITPDF(0)
      CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)



c - Statistics part
c - Use statistics tables
c - Call statistical error-code for scenario
      IF (LSTAT) THEN
         WRITE(*,*)"****************************************"//
     >        "********************************"
         WRITE(*,*)"ALLUNC: Evaluating statistical uncertainties"
         WRITE(*,*)"****************************************"//
     >        "********************************"
ckr Replace old and crosschecked code
ckr         CALL STATCODE(TABPATH,SCENARIO,BORNN,NLON)
         BORNNAME = TABPATH(1:LENOCC(TABPATH))//"/stat/"//
     >        SCENARIO(1:LENOCC(SCENARIO))//"-hhc-born-"
         NLONAME  = TABPATH(1:LENOCC(TABPATH))//"/stat/"//
     >        SCENARIO(1:LENOCC(SCENARIO))//"-hhc-nlo-"
         ISCL = 3
         XMUR = MURSCALE(ISCL)
         XMUF = MUFSCALE(ISCL)
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
            CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
            IPHASE = 2

c - Loop over files
            IF (IORD.EQ.0) THEN
c - Total x section
               NTAB = NLON 
               BWGT = 1D0/1D8
               FILEBASE = NLONAME(1:LENOCC(NLONAME))
            ELSEIF (IORD.EQ.1) THEN
c - LO
               NTAB = BORNN
               BWGT = 1D0/1D9
               FILEBASE = BORNNAME(1:LENOCC(BORNNAME))
            ELSEIF (IORD.EQ.2) THEN
c - NLO
               NTAB = NLON
               BWGT = 1D0/1D8
               FILEBASE = NLONAME(1:LENOCC(NLONAME))
            ELSE
               WRITE(*,*)"ALLUNC: ERROR! "//
     >              "Illegal order for stat. calc:",IORD
               STOP
            ENDIF
            DO ITAB=0,NTAB
               WRITE(NO,'(I4.4)'),ITAB
               FILENAMES = FILEBASE(1:LENOCC(FILEBASE))//"2jet_"//NO
     >              //".tab"
               OPEN(2,STATUS='OLD',FILE=FILENAMES,IOSTAT=ISTAT)
               IF (ISTAT.NE.0) THEN
                  FILENAMES = FILEBASE(1:LENOCC(FILEBASE))//"3jet_"
     >                 //NO//".tab"
                  OPEN(2,STATUS='OLD',FILE=FILENAMES,IOSTAT=ISTAT)
                  IF (ISTAT.NE.0) THEN
ckr                     WRITE(*,*)"Filename for order",IORD,":",FILENAMES
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
               CALL FX9999CC(FILENAMES,XMUR,XMUF,0,XSECT0)
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
                  IF (MYRES(IBIN,NSUBPROC+1,IORD2).GE.0.D0) THEN
                     DSTMP(1) = (WTXMIN(IBIN,NSUBPROC+1,IORD2)-1D0)
     >                    *100D0
                     DSTMP(2) = (WTXMAX(IBIN,NSUBPROC+1,IORD2)-1D0)
     >                    *100D0
                     IF (WTDXU2(IBIN,NSUBPROC+1,IORD2).GT.1D-99) THEN
                        DSTMP(3) = 
     >                       (WTXMIN(IBIN,NSUBPROC+1,IORD2)-1D0) /
     >                       WTDXU2(IBIN,NSUBPROC+1,IORD2)
                        DSTMP(4) = 
     >                       (WTXMAX(IBIN,NSUBPROC+1,IORD2)-1D0) /
     >                       WTDXU2(IBIN,NSUBPROC+1,IORD2)
                     ENDIF
                  ELSE
                     DSTMP(1) = (WTXMIN(IBIN,NSUBPROC+1,IORD2)+1D0)
     >                    *100D0
                     DSTMP(2) = (WTXMAX(IBIN,NSUBPROC+1,IORD2)+1D0)
     >                    *100D0
                     IF (WTDXU2(IBIN,NSUBPROC+1,IORD2).GT.1D-99) THEN
                        DSTMP(3) = 
     >                       (WTXMIN(IBIN,NSUBPROC+1,IORD2)+1D0) /
     >                       WTDXU2(IBIN,NSUBPROC+1,IORD2)
                        DSTMP(4) = 
     >                       (WTXMAX(IBIN,NSUBPROC+1,IORD2)+1D0) /
     >                       WTDXU2(IBIN,NSUBPROC+1,IORD2)
                     ENDIF
                  ENDIF
                  WRITE(*,902) IBIN,
     >                 IJMIN(IBIN),
     >                 IJMAX(IBIN),
     >                 MYRES(IBIN,NSUBPROC+1,IORD2),
     >                 WTXMIN(IBIN,NSUBPROC+1,IORD2) *
     >                 DABS(MYRES(IBIN,NSUBPROC+1,IORD2)),
     >                 WTXMAX(IBIN,NSUBPROC+1,IORD2) *
     >                 DABS(MYRES(IBIN,NSUBPROC+1,IORD2)),
     >                 WTDXMN(IBIN,NSUBPROC+1,IORD2)*100D0,
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
     >           "cross section, scale no. 3, "//
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
                  WRITE(*,901) IBIN,MYRES(IBIN,NSUBPROC+1,IORD2),
     >                 WTDXMN(IBIN,NSUBPROC+1,IORD2)
               ENDDO
            ENDDO
c - Fill histograms
ckr Replaces STATCODE
            CALL PDFFILL(NRAP,3,IORD,ISCL,WTDXMN)
ckr Put replacement for STATCODE here
            CALL PDFFILL(NRAP,4,IORD,ISCL,WTDXUL)
            
         ENDDO
      ENDIF
      
      
      
c - alpha_s uncertainty part
c - Use primary table
c - Check that FILENAME is still the primary table here ...!!!
c - alpha_s(M_Z) variations, either implicit via PDF set or explicitly
c - Uncertainty filled at central scale no. 3
c - (ISCL=3 in FORTRAN, refscale=2 in C++ parlance of author code)
      IF (LSER) THEN
         WRITE(*,*)"****************************************"//
     >        "********************************"
         WRITE(*,*)"ALLUNC: Evaluating alpha_s uncertainties"
         WRITE(*,*)"****************************************"//
     >        "********************************"
ckr Attention: Assumption here is that delta(alpha_s) = +- 0.00n
ckr            equals a change in PDF set member number by n
ckr            This is true for:
ckr                              cteq66alphas:           2 +- up to  2
ckr                              CT10as:                 5 +- up to  5
ckr                              MSTW2008nlo_asmzrange: 11 +- up to 10
ckr Central member number  : Integer part of DASMZVAL
ckr Up/down member numbers : Integer part of (DASMZVAL-NINT(DASMZVAL)) times 1000 
ckr alpha_s variation alone: (DASMZVAL-NINT(DASMZVAL))
         CALL INITPDF(0)
         CALL NUMBERPDF(MYPDFAS)
         IPDF = NINT(ABS(DASMZVAL))
         WRITE(*,*)"ALLUNC: PDF member for central alpha_s value:"
     >        ,IPDF
         WRITE(*,*)"ALLUNC: Number of alpha_s variations in set:"
     >        ,MYPDFAS
         IF (ASMODE.EQ."PDF".AND.IPDF.GT.MYPDFAS) THEN
            WRITE(*,*)"ALLUNC: ERROR! Central PDF member for,"//
     >           " alpha_s variation does not exist, aborting!"//
     >           "IPDF = ",IPDF
            STOP
         ENDIF
         IPDFUD = NINT(1000.D0*(ABS(DASMZVAL)-DBLE(IPDF)))
         WRITE(*,*)"ALLUNC: alpha_s variation * 1000 here:"
     >        ,IPDFUD
         IF (ASMODE.EQ."PDF".AND.
     >        (IPDF-IPDFUD.LT.0.OR.IPDF+IPDFUD.GT.MYPDFAS)) THEN
            WRITE(*,*)"ALLUNC: ERROR! Varied PDF members for,"//
     >           " alpha_s variation do not exist, aborting!"//
     >           "IPDFUD = ",IPDFUD
            STOP
         ENDIF

         WRITE(*,*)"========================================"//
     >        "================================"
         WRITE(*,*)"Relative alpha_s Uncertainties"
         WRITE(*,*)"- the printed values are for the total "//
     >        "cross section summed over all subprocesses"
         WRITE(*,*)"- histograms contain more detailed results"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         WRITE(*,*)" bin       cross section           "//
     >        "lower a_s uncertainty   upper a_s uncertainty"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         
         ISCL = 3
         XMUR = MURSCALE(ISCL)
         XMUF = MUFSCALE(ISCL)

         IF (ASMODE.EQ."PDF") THEN
ckr alpha_s variations in PDF set
            IPHASE  = 1
            IMODE   = 3
            IWEIGHT = 0
            NRAP  = NRAPIDITY
            IF (LRAT.OR.LNRM) THEN
               NRAP = 2*NRAPIDITY
            ENDIF
            CALL INITPDF(IPDF)
cnew
            ASMZPDF = ALPHASPDF(ZMASS)
cnew
            CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
            ISTEP = 0
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
            IF (LNRM) THEN
ckr Load normalization table with potentially different binning!
               IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
               ISTEP = 1
               CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
               IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
            ENDIF
            ISTEP = 2
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))

            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
            IPHASE = 2
            
            CALL INITPDF(IPDF-IPDFUD)
cnew
            ASMZTMP = ALPHASPDF(ZMASS)
            ASDN = ASMZTMP - ASMZPDF
cnew
            CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
            IF (LNRM) THEN
               ISTEP = 3
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
               IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
               ISTEP = 4
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
               IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
               ISTEP = 5
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
            ENDIF
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
            
            CALL INITPDF(IPDF+IPDFUD)
cnew
            ASMZTMP = ALPHASPDF(ZMASS)
            ASUP = ASMZTMP - ASMZPDF
cnew
            CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
            IF (LNRM) THEN
               ISTEP = 3
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
               IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
               ISTEP = 4
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
               IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
               ISTEP = 5
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
            ENDIF
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
            
            IPHASE = 3
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

c - Give some standard output, fill histograms
            WRITE(*,*)"ALLUNC: Uncertainties from"//
     >           " alpha_s variations in PDF members"
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            ISCL = 3
            IBIN   = 0
            DO IRAP=1,NRAP
               DO IPT=1,NPT(IRAP)
                  IBIN = IBIN+1
                  WRITE(*,900) IBIN,MYRESN(IBIN,NSUBPROC+1,NORD+1),
     >                 WTDXLM(IBIN,NSUBPROC+1,NORD+1),
     >                 WTDXUM(IBIN,NSUBPROC+1,NORD+1)
               ENDDO
            ENDDO
            CALL PDFFILL(NRAP,0,-1,ISCL,MYRESN)
            CALL PDFFILL(NRAP,1,-1,ISCL,WTDXLM)
            CALL PDFFILL(NRAP,2,-1,ISCL,WTDXUM)
         ENDIF

ckr Standalone alpha_s variations
         IPHASE  = 1
         IMODE   = 3
         IWEIGHT = 0
         NRAP  = NRAPIDITY
         IF (LRAT.OR.LNRM) THEN
            NRAP = 2*NRAPIDITY
         ENDIF
         CALL INITPDF(IPDF)
         ASMODETMP = ASMODE
         ASMODE    = "KR"
         ASMZTMP   = ASMZVAL
         ASMZVAL   = ALPHASPDF(ZMASS)

         CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
         ISTEP = 0
         CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
         IF (LNRM) THEN
ckr Load normalization table with potentially different binning!
            IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
            ISTEP = 1
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
            IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
         ENDIF
         ISTEP = 2
         CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
         
         CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
         IPHASE = 2
         
cnew         ASMZVAL = ASMZVAL - DBLE(IPDFUD)/1000.
         ASMZVAL = ASMZVAL + ASDN
         CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
         IF (LNRM) THEN
            ISTEP = 3
            CALL CENRES(ISTEP,LRAT,LNRM,
     >           SCENARIO(1:LENOCC(SCENARIO)))
            IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
            ISTEP = 4
            CALL CENRES(ISTEP,LRAT,LNRM,
     >           SCENARIO(1:LENOCC(SCENARIO)))
            IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
            ISTEP = 5
            CALL CENRES(ISTEP,LRAT,LNRM,
     >           SCENARIO(1:LENOCC(SCENARIO)))
         ENDIF
         CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
         
ckr Add twice to correct previous subtraction!
cnew         ASMZVAL = ASMZVAL + DBLE(2*IPDFUD)/1000.
         ASMZVAL = ASMZVAL - ASDN + ASUP
         CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
         IF (LNRM) THEN
            ISTEP = 3
            CALL CENRES(ISTEP,LRAT,LNRM,
     >           SCENARIO(1:LENOCC(SCENARIO)))
            IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
            ISTEP = 4
            CALL CENRES(ISTEP,LRAT,LNRM,
     >           SCENARIO(1:LENOCC(SCENARIO)))
            IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
            ISTEP = 5
            CALL CENRES(ISTEP,LRAT,LNRM,
     >           SCENARIO(1:LENOCC(SCENARIO)))
         ENDIF
         CALL UNCERT(IPHASE,IMODE,IWEIGHT,J,LRAT,LNRM)
         
         IPHASE = 3
         CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

         ASMODE  = ASMODETMP
         ASMZVAL = ASMZTMP

c - Give some standard output, fill histograms
         WRITE(*,*)"ALLUNC: Uncertainties from"//
     >        " standalone alpha_s variations"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         ISCL = 3
         IBIN   = 0
         DO IRAP=1,NRAP
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               WRITE(*,900) IBIN,MYRESN(IBIN,NSUBPROC+1,NORD+1),
     >              WTDXLM(IBIN,NSUBPROC+1,NORD+1),
     >              WTDXUM(IBIN,NSUBPROC+1,NORD+1)
            ENDDO
         ENDDO
         IF (.NOT.ASMODE.EQ."PDF") THEN
            CALL PDFFILL(NRAP,0,-1,ISCL,MYRESN)
         ENDIF
         CALL PDFFILL(NRAP,3,-1,ISCL,WTDXLM)
         CALL PDFFILL(NRAP,4,-1,ISCL,WTDXUM)
      ENDIF
      


c - Scale uncertainty part
c - Use primary table
c - Check that FILENAME is still the primary table here ...!!!
c - Two schemes are implemented for hadron-hadron:
c - 1. NSCALEVAR = 4 with standard mur, muf factors of
c      (1/4,1/4), (1/2,1/2), (  1,  1), (  2,  2)
c - 2. NSCALEVAR = 8 with additional mur, muf combinations
c      (1/4,1/4), (1/2,1/2), (  1,  1), (  2,  2) as before plus
c      (  1,1/2), (  1,  2), (1/2,  1), (  2,  1)
c - Uncertainty filled at central scale no. 3
c - (ISCL=3 in FORTRAN, refscale=2 in C++ parlance of author code)
      IF (LSCL) THEN
         WRITE(*,*)"****************************************"//
     >        "********************************"
         WRITE(*,*)"ALLUNC: Evaluating scale uncertainties"
         WRITE(*,*)"****************************************"//
     >        "********************************"
c - 2-point scheme
         NSCLS = NSCALEVAR
         IF (NSCALEVAR.GT.4) THEN
            NSCLS = 4
         ENDIF
         CALL INITPDF(0)
         ISCL = 3
         XMUR = MURSCALE(ISCL)
         XMUF = MUFSCALE(ISCL)
         IPHASE  = 1
         IMODE   = 3
         IWEIGHT = 0
         NRAP  = NRAPIDITY
         IF (LRAT.OR.LNRM) THEN
            NRAP = 2*NRAPIDITY
         ENDIF
         CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
         ISTEP = 0
         CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
         IF (LNRM) THEN
ckr Load normalization table with potentially different binning!
            IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
            ISTEP = 1
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
            IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
         ENDIF
         ISTEP = 2
         CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))

         CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
         IPHASE = 2

         DO ISCL=1,NSCLS
ckr Do neither use scale 1 with factor of 1/4 nor default scale 3
ckr Ugly goto construction avoidable with f90 CYCLE command
            IF (ISCL.EQ.1.OR.ISCL.EQ.3) GOTO 10
            XMUR = MURSCALE(ISCL)
            XMUF = MUFSCALE(ISCL)
            CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
            IF (LNRM) THEN
               ISTEP = 3
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
ckr Load normalization table with potentially different binning!
               IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
               ISTEP = 4
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
               IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
               ISTEP = 5
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
            ENDIF
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,ISCL,LRAT,LNRM)
 10         CONTINUE
         ENDDO
         IPHASE = 3
         CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
ckr No log file output, but store in histograms
         ISCL = 3
         CALL PDFFILL(NRAP,6,-1,ISCL,WTDXLM)
         CALL PDFFILL(NRAP,7,-1,ISCL,WTDXUM)

c - 6-point scheme
         NSCLS = NSCALEVAR
         IF (NSCALEMAX.GE.8.AND.NSCALEVAR.EQ.4) THEN
            NSCLS = 8
            MURSCALE(5) = 1.0D0
            MUFSCALE(5) = 0.5D0
            MURSCALE(6) = 1.0D0
            MUFSCALE(6) = 2.0D0
            MURSCALE(7) = 0.5D0
            MUFSCALE(7) = 1.0D0
            MURSCALE(8) = 2.0D0
            MUFSCALE(8) = 1.0D0
         ENDIF
         CALL INITPDF(0)
         ISCL = 3
         XMUR = MURSCALE(ISCL)
         XMUF = MUFSCALE(ISCL)
         IPHASE  = 1
         IMODE   = 3
         IWEIGHT = 0
         NRAP  = NRAPIDITY
         IF (LRAT.OR.LNRM) THEN
            NRAP = 2*NRAPIDITY
         ENDIF
         CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
         ISTEP = 0
         CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
         IF (LNRM) THEN
ckr Load normalization table with potentially different binning!
            IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
            ISTEP = 1
            CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
            IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
         ENDIF
         ISTEP = 2
         CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))

         CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
         IPHASE = 2

         DO ISCL=1,NSCLS
ckr Do neither use scale 1 with factor of 1/4 nor default scale 3
ckr Ugly goto construction avoidable with f90 CYCLE command
            IF (ISCL.EQ.1.OR.ISCL.EQ.3) GOTO 11
            XMUR = MURSCALE(ISCL)
            XMUF = MUFSCALE(ISCL)
            CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
            IF (LNRM) THEN
               ISTEP = 3
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
ckr Load normalization table with potentially different binning!
               IF (LTAB) CALL FX9999CC(FILENAMN,XMUR,XMUF,0,XSECT0)
               ISTEP = 4
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
               IF (LTAB) CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
               ISTEP = 5
               CALL CENRES(ISTEP,LRAT,LNRM,
     >              SCENARIO(1:LENOCC(SCENARIO)))
            ENDIF
            CALL UNCERT(IPHASE,IMODE,IWEIGHT,ISCL,LRAT,LNRM)
 11         CONTINUE
         ENDDO
         IPHASE = 3
         CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

c - Give some standard output, fill histograms
         WRITE(*,*)"========================================"//
     >        "================================"
         WRITE(*,*)"Relative Scale Uncertainties"
            WRITE(*,*)"- the printed values are for the total "//
     >           "cross section, scale no. 3, "//
     >           "summed over all subprocesses"
         WRITE(*,*)"- histograms contain more detailed results"
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
         WRITE(*,*)" bin       cross section           "//
     >        "lower scale uncertainty upper scale uncertainty"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         WRITE(*,*)"ALLUNC: Uncertainties from",NSCLS-2,
     >        " scale variations"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         IORD   = 0
         ISCL = 3
         ISUB   = 0
         IBIN   = 0
         DO IRAP=1,NRAP
Comment:             IHIST = IORD*1000000 + ISCL*100000 +
Comment:      >           ISUB*10000 + IRAP*100
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
Comment:                PT(IBIN) = REAL(PTBIN(IRAP,IPT))
               WRITE(*,900) IBIN,MYRESN(IBIN,NSUBPROC+1,NORD+1),
     >              WTDXLM(IBIN,NSUBPROC+1,NORD+1),
     >              WTDXUM(IBIN,NSUBPROC+1,NORD+1)
            ENDDO
         ENDDO
         CALL PDFFILL(NRAP,8,-1,ISCL,WTDXLM)
         CALL PDFFILL(NRAP,9,-1,ISCL,WTDXUM)
      ENDIF



c - Algorithmic part
c - Use reference table
c - Reference result is always stored in scale variation no. 1
c - Reference table has doubled number of rapidity bins and
c -   the reference result is stored in upper half of these bins
c - Reference result is evaluated with CTEQ61 PDFs
c - Default scale in hadron-hadron usually is the third variation:
c -   C++ no. 2 --> Fortran no. 3 ==> normal result (but compare to author code!)
c - (Use other scale in FORTRAN ONLY if refscale <> 2 in author code)
c - Attention: From now on ref. table loaded ==> rap. bins doubled,
c                  no ratio calcs below ...
      LRAT = .FALSE.
      IF (LALG) THEN
         WRITE(*,*)"****************************************"//
     >        "********************************"
         WRITE(*,*)"ALLUNC: Evaluating algorithmic uncertainties"
         WRITE(*,*)"****************************************"//
     >        "********************************"
         FILENAME = TABPATH(1:LENOCC(TABPATH))//"/"//REFNAME
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         WRITE(*,*)"ALLUNC: Taking reference table: "//
     >        FILENAME(1:LENOCC(FILENAME))
c - Initialize CTEQ61 reference PDFs
         PDFSET = PDFPATH(1:LENOCC(PDFPATH))//"/cteq61.LHgrid"
         WRITE(*,*)"ALLUNC: Taking reference PDF: "//
     >        PDFSET(1:LENOCC(PDFSET))
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         CALL INITPDFSET(PDFSET(1:LENOCC(PDFSET)))
         CALL INITPDF(0)

c - Reference result
c - Need scale variation no. 1 for this, upper half of rap. bins
c - NOTE1: The values of mur and muf are internally used to select the
c          the corresponding table number, so they have to be set
c          according to scale variation 1 although the factors are wrong
c          for the comparison ==>
c            compare to scale variation 3, lower half of rap. bins
c - NOTE2: The reference result should not depend on the mur or muf values
c          apart from the table choice, BUT if mur is not equal to muf
c          additional factors are taken into account which must be
c          avoided for the reference calculation!
         ISCL = 1
         XMUR = MURSCALE(ISCL)
         XMUF = MUFSCALE(ISCL)
ckr         CALL FX9999CC(FILENAME,XMUR,XMUF,1,XSECT0)
         CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
         IPHASE  = 1
         IMODE   = 4
         IWEIGHT = 0
         CALL CENRES(ISTEP,LRAT,LNRM,SCENARIO(1:LENOCC(SCENARIO)))
         CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

c - Normal result
         ISCL = 3
         XMUR = MURSCALE(ISCL)
         XMUF = MUFSCALE(ISCL)
         CALL FX9999CC(FILENAME,XMUR,XMUF,0,XSECT0)
         IPHASE = 2
         CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)
         IPHASE = 3
         CALL UNCERT(IPHASE,IMODE,IWEIGHT,0,LRAT,LNRM)

c - Give some standard output, fill histograms
         WRITE(*,*)"========================================"//
     >        "================================"
         WRITE(*,*)"Relative Algorithmic Error"
         WRITE(*,*)"- the printed values are for the total "//
     >        "cross section summed over all subprocesses"
         WRITE(*,*)"- histograms contain more detailed results"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         WRITE(*,*)" bin       ref. cross section      "//
     >        "algorithmic deviation"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         ISCL = 3
         IBIN   = 0
         DO IRAP=1,INT(NRAPIDITY/2)
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               WRITE(*,901) IBIN,WTX(IBIN,NSUBPROC+1,NORD+1),
     >              WTDXUM(IBIN,NSUBPROC+1,NORD+1)
            ENDDO
         ENDDO
         CALL PDFFILL(NRAPIDITY,5,-1,ISCL,WTDXUM)
      ENDIF



c - Close hbook file
      CALL PDFHIST(2,HISTFILE,LONE,LPDF,LSTAT,LALG,LSER,MYPDF,
     >     LRAT.OR.LNRM,LSCL)
      END



c
c ======================= Book the histograms ========================
c
      SUBROUTINE PDFHIST(N,HISTFILE,
     &     LONE,LPDF,LSTAT,LALG,LSER,MYPDF,LRAT,LSCL)
      IMPLICIT NONE
      CHARACTER*(*) HISTFILE
      CHARACTER*255 CSTRNG,CBASE1,CBASE2,CTMP
      INTEGER N,LENOCC,IPDF,MYPDF,IPTMAX,NRAP
      LOGICAL LONE,LPDF,LSTAT,LALG,LSER,LRAT,LSCL

      INTEGER I,J,ISTAT2,ICYCLE,NSCLS
      INTEGER IORD,ISUB,ISCL,IRAP,IPT,IHIST,NHIST
      INCLUDE "fnx9999.inc"
      INCLUDE "strings.inc"
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
         CBASE1 = CIPROC(IPROC)
         CBASE1 = CBASE1(1:LENOCC(CBASE1))//"_"
     >        //NAMELABEL(1)
         CBASE2 = CBASE1(1:LENOCC(CBASE1))//"_"
     >        //CIALGO(IALGO)
         CBASE2 = CBASE2(1:LENOCC(CBASE2))//"_"
     >        //CJETRES1(IALGO)
         WRITE(CTMP,'(F3.1)'),JETRES1
         CBASE2 = CBASE2(1:LENOCC(CBASE2))//"="
     >        //CTMP
         DO IORD=0,NORD         ! Order: tot, LO, NLO-corr, NNLO-corr
            NSCLS = NSCALEVAR
            IF (NSCALEMAX.GE.8.AND.NSCALEVAR.EQ.4) THEN
               NSCLS = 8
            ENDIF
            DO ISCL=1,NSCLS ! Scale variations
               DO ISUB=0,NSUBPROC ! Subprocesses: 0 tot + 7 subproc
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
                     CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//"_iord="
     >                    //CTMP
                     WRITE(CTMP,'(I1)'),ISCL
                     CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//"_imu="
     >                    //CTMP
                     WRITE(CTMP,'(I1)'),IRAP
                     CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//"_iy="
     >                    //CTMP
                     CALL HBOOKB(IHIST,
     >                    CSTRNG(1:LENOCC(CSTRNG)),
     >                    NPT(IRAP),PT,0)
ckr                     CALL HBARX(IHIST)
                     NHIST = NHIST+1
ckr                     write(*,*)"1. Booked histo #",nhist
Comment:                      IF (LSER) THEN
Comment:                         DO IPDF=1,MYPDF
Comment:                            CALL HBOOKB(IHIST+IPDF,
Comment:      >                          CSTRNG(1:LENOCC(CSTRNG)),
Comment:      >                          NPT(IRAP),PT,0)
Comment:                            NHIST = NHIST+1
Comment:                         ENDDO
Comment:                      ENDIF
                     IF (LPDF) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_PDF_low/xsect"
                        CALL HBOOKB(IHIST+1,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_PDF_up/xsect"
                        CALL HBOOKB(IHIST+2,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+2
ckr                        write(*,*)"2. Booked histo #",nhist
                     ELSEIF (LSER.AND.ISCL.EQ.3) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_asPDF_low/xsect"
                        CALL HBOOKB(IHIST+1,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_asPDF_up/xsect"
                        CALL HBOOKB(IHIST+2,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+2
ckr                        write(*,*)"2. Booked histo #",nhist
                     ENDIF
                     IF (LSTAT) THEN
                        IF (IORD.LE.2.AND.
     >                       ISCL.EQ.3.AND.ISUB.EQ.0) THEN
                           CSTRNG = CBASE1
                           CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                          "_dstat/xsect"
                           CALL HBOOKB(IHIST+3,
     >                          CSTRNG(1:LENOCC(CSTRNG)),
     >                          NPT(IRAP),PT,0)
                           CSTRNG = CBASE1
                           CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                          "_dmax/2/xsect"
                           CALL HBOOKB(IHIST+4,
     >                          CSTRNG(1:LENOCC(CSTRNG)),
     >                          NPT(IRAP),PT,0)
                           NHIST = NHIST+2
ckr                        write(*,*)"3. Booked histo #",nhist
                        ENDIF
                     ELSEIF (LSER.AND.ISCL.EQ.3) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_as_low/xsect"
                        CALL HBOOKB(IHIST+3,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_as_up/xsect"
                        CALL HBOOKB(IHIST+4,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+2
ckr                        write(*,*)"3. Booked histo #",nhist
                     ENDIF
                     IF (LALG.AND.IORD.LE.2.AND.
     >                    ISCL.EQ.3) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_xsect/xsect_ref"
                        CALL HBOOKB(IHIST+5,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
ckr                        write(*,*)"4. Booked histo #",nhist
                     ENDIF
                     IF (LSCL.AND.IORD.LE.2.AND.
     >                    ISCL.EQ.3) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_SCL2_low/xsect"
                        CALL HBOOKB(IHIST+6,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_SCL2_up/xsect"
                        CALL HBOOKB(IHIST+7,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_SCL6_low/xsect"
                        CALL HBOOKB(IHIST+8,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_SCL6_up/xsect"
                        CALL HBOOKB(IHIST+9,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
ckr                        write(*,*)"4b. Booked histo #",nhist
                     ENDIF
                     IF (LSTAT.AND.ISUB.EQ.0) THEN
                        IHIST = IORD*1000000 + ISCL*100000 +
     >                       ISUB*10000 + IRAP*100
                        CSTRNG = CBASE1(1:LENOCC(CBASE1))
                        WRITE(CTMP,'(I1)'),IORD
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//"_iord="
     >                       //CTMP
                        WRITE(CTMP,'(I1)'),ISCL
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//"_imu="
     >                       //CTMP
                        WRITE(CTMP,'(I1)'),IRAP
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//"_iy="
     >                       //CTMP
                        CTMP = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_norm_mu_sig_all_pt"
                        CALL HBOOK1(IHIST + 10,
     >                       CTMP(1:LENOCC(CTMP)),
     >                       63,-10.5,10.5,0)
                        CALL HIDOPT(IHIST + 10,'STAT')
                        CTMP = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_norm_mu_dmax_all_pt"
                        CALL HBOOK1(IHIST + 11,
     >                       CTMP(1:LENOCC(CTMP)),
     >                       30,-1.5,1.5,0)
                        CALL HIDOPT(IHIST + 11,'STAT')
                        NHIST = NHIST+2
ckr                        write(*,*)"5. Booked histo #, IHIST",nhist,
ckr     >                       IHIST + 11
ckr IHIST limit before next rapidity bin: IHIST+2*IPT+11 < IHIST+100
ckr => Maximal IPT = IPTMAX < 45
                        IPTMAX = 44
                        DO IPT=1,MIN(NPT(IRAP),IPTMAX)
                           CSTRNG = CBASE1(1:LENOCC(CBASE1))
                           WRITE(CTMP,'(I1)'),IORD
                           CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//"_iord="
     >                          //CTMP
                           WRITE(CTMP,'(I1)'),ISCL
                           CSTRNG = CSTRNG(1:LENOCC(CSTRNG))/
     >                          /"_imu="//CTMP
                           WRITE(CTMP,'(I1)'),IRAP
                           CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//"_iy="
     >                          //CTMP
                           WRITE(CTMP,'(I2)'),IPT
                           CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//"_ipt="
     >                          //CTMP
                           CTMP = CSTRNG(1:LENOCC(CSTRNG))//
     >                          "_norm_mu_sig"
                           CALL HBOOK1(IHIST + 2*IPT + 10,
     >                          CTMP(1:LENOCC(CTMP)),
     >                          63,-10.5,10.5,0)
                           CALL HIDOPT(IHIST + 2*IPT + 10,'STAT')
                           CTMP = CSTRNG(1:LENOCC(CSTRNG))//
     >                          "_norm_mu_dmax"
                           CALL HBOOK1(IHIST + 2*IPT + 11,
     >                          CTMP(1:LENOCC(CTMP)),
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
      INTEGER NRAP,IOFF,IORDFIL,ISCL
      DOUBLE PRECISION DVAL(NBINTOTMAX,NMAXSUBPROC+1,4)
      INTEGER I,J,IBIN,IORD,IORD2,ISUB,ISUB2,IHIST
      REAL RVAL
      
      IF (ISCL.LT.1 .OR. ISCL.GT.NSCALEVAR) THEN
         WRITE(*,*) "\nPDFFILL: ERROR! ISCL ",ISCL,
     >        " is out of range, aborted!"
         WRITE(*,*) "PDFFILL: Max. ISCL: ",NSCALEVAR
         STOP
      ENDIF
      
c - Fill all histograms for the given scale
c - Fill sums over subprocesses and/or orders into zero factors for IHIST
      DO IORD2=1,NORD+1         ! Order: LO, NLO-corr, NNLO-corr, ... , tot
         IORD = IORD2
         IF (IORD2.EQ.NORD+1) IORD = 0 
         IF (IORDFIL.EQ.-1.OR.IORDFIL.EQ.IORD) THEN
            DO ISUB2=1,NSUBPROC+1 ! Subprocesses hh: 7 subprocesses + total
               ISUB = ISUB2
               IF (ISUB2.EQ.NSUBPROC+1) ISUB=0
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
