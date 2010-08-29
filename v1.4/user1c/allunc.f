      PROGRAM ALLUNC
* ---------------------------------------------------------------------
* K. Rabbertz 07.09.2008 First try to integrate all uncertainties
*                        into one job
*
* ALLUNC - Program to derive the algorithmic, statistical and PDF
*          uncertainties using fastNLO tables
*
* ---------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      INCLUDE "uncert.inc"
      CHARACTER*255 HISTFILE
      CHARACTER*255 SCENARIO,FILENAME,TABPATH,TABNAME,REFNAME
      CHARACTER*255 PDFSET,PDFPATH,LHAPDF,CHTMP
      CHARACTER*8 CH8TMP
      CHARACTER*4 CH4TMP
      INTEGER BORNN,NLON,LENOCC
      INTEGER I,J,L1,L2,L3,L4,NPDF,IOPDF,IOAS,NSCALES
      INTEGER ISTAT,ISCALE,IORD,IBIN,NBIN,ISUB,IRAP,IPT,IHIST
      LOGICAL LONE,LPDF,LSTAT,LSER,LSCL,LRAT,LALG
      LOGICAL LTOY
      DOUBLE PRECISION MUR,MUF,DIFF,SUMM,QLAM4,QLAM5
      DOUBLE PRECISION
     >     RES0(NBINTOTMAX,NMAXSUBPROC+1,0:3),
     >     RES1HI(NBINTOTMAX,NMAXSUBPROC+1,0:3),
     >     RES1LO(NBINTOTMAX,NMAXSUBPROC+1,0:3),
     >     RESLO,RESHI,DREF
c - To unify quoted uncertainties (CL68,CL90,special)
c - Convert from CL68 to CL90 values
c - TOCL90 = 1.64485D0 ! SQRT(2.D0)/InvERF(0.9D0)
c - Convert from GJR to CTEQ (CL90) values ?
c - TOCL90GJR = 2.12766D0! 1.D0/0.47D0
      DOUBLE PRECISION TOCL90,TOCL90GJR
      PARAMETER (TOCL90 = 1.64485D0, TOCL90GJR = 2.12766D0)
      
c - Attention!!! This must be declared consistent with the
c                definition in the commonblock!!!!!
      DOUBLE PRECISION XSECT0(NBINTOTMAX,3),XSECT1(NBINTOTMAX,3)
      REAL PT(NPTMAX)
      INTEGER IMODE,NRAP

      CHARACTER*255 ASMODE
      DOUBLE PRECISION ASMZVAL
      INTEGER IASLOOP
      COMMON/STEER/ASMZVAL,IASLOOP,ASMODE

c --- Parse command line
      WRITE(*,*)"\n #################################################"
      WRITE(*,*)"# ALLUNC"
      WRITE(*,*)"#################################################"
      WRITE(*,*)"# Program to derive the algorithmic, statistical, "
      WRITE(*,*)"# and PDF uncertainties using fastNLO tables"
      WRITE(*,*)"#################################################"
      WRITE(*,*)"#"

*---Scenario
      LRAT = .FALSE.
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
C --- Use '...' with \", otherwise gfortran complains 
            WRITE(*,*)'     \"scenario\".tab'
            WRITE(*,*)'     \"scenario\"ref.tab'
            WRITE(*,*)'     The statistics tables named e.g.:'
            WRITE(*,*)'     \"scenario\"-hhc-nlo-2jet_nnnn.tab'
            WRITE(*,*)'     have to be in the subdir. \"stat\".'
            WRITE(*,*)'  HBOOK output file, def. = scenario.hbk'
            WRITE(*,*)'  Derive algorithmic uncertainty, def. = no'
            WRITE(*,*)'  Last LO stat. table number, def. = -1'
            WRITE(*,*)'  Last NLO stat. table number, def. = -1'
            WRITE(*,*)'  PDF set, def. = cteq65.LHgrid'
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
            WRITE(*,*)'  alpha_s loop order, def. from PDF set'
            WRITE(*,*)'  Use MC sampling method for PDF uncertainty,'//
     >           ' def. = no'
            WRITE(*,*)' '
            STOP
         ELSEIF (SCENARIO(1:LENOCC(SCENARIO)).EQ."fnl2442".OR.
     >           SCENARIO(1:LENOCC(SCENARIO)).EQ."fnl2442a") THEN
            LRAT = .TRUE.
            WRITE(*,*)
     >           "ALLUNC: Deriving x section ratios"
         ENDIF
         WRITE(*,*)"ALLUNC: Evaluating scenario: ",
     >        SCENARIO(1:LENOCC(SCENARIO))
      ENDIF
      TABNAME = SCENARIO(1:LENOCC(SCENARIO))//".tab"
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
      FILENAME = TABPATH(1:LENOCC(TABPATH))//"/"//TABNAME
      WRITE(*,*)"ALLUNC: Taking primary table ",
     >     FILENAME(1:LENOCC(FILENAME))

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
     >         "ALLUNC: WARNING! Do NOT derive algorithmic uncertainty!"
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
         PDFSET = "cteq65.LHgrid"
         WRITE(*,*)
     >        "ALLUNC: WARNING! No PDF set given, "//
     >        "taking cteq65.LHgrid instead!"
      ELSE
         WRITE(*,*)"ALLUNC: Using PDF set: ",
     >        PDFSET(1:LENOCC(PDFSET))
      ENDIF

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

*---alpha_s loop order in evolution
      CH4TMP = "X"
      IF (IARGC().GE.11) THEN
         CALL GETARG(11,CH4TMP)
      ENDIF
      IF (IARGC().LT.11.OR.CH4TMP(1:1).EQ."_") THEN
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
      IF (IARGC().GE.12) THEN
         CALL GETARG(12,CH4TMP)
      ENDIF
      IF (IARGC().LT.12.OR.CH4TMP(1:1).EQ."_".OR.
     >     CH4TMP(1:2).EQ."no") THEN
         LTOY   =  .FALSE.
         WRITE(*,*)
     >        "ALLUNC: Use eigen vector (CTEQ/MSTW) method."
      ELSE
         IF (CH4TMP(1:3).EQ."yes") THEN
            LTOY = .TRUE.
            WRITE(*,*)"ALLUNC: Use toy MC (NNPDF) method."
         ELSE
            LTOY = .FALSE.
            WRITE(*,*)
     >         "ALLUNC: WARNING! Use eigen vector (CTEQ/MSTW) method!"
         ENDIF
      ENDIF

*---Too many arguments
      IF (IARGC().GT.12) THEN
         WRITE(*,*)"\nALLUNC: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF
      WRITE(*,*)" "



c - Initialize LHAPDF
      CALL INITPDFSET(PDFSET(1:LENOCC(PDFSET)))

c - Initialize one member, 0=best fit member
      CALL INITPDF(0)

c - Write out some info on best fit member      
      CALL NUMBERPDF(NPDF)
      CALL GETORDERPDF(IOPDF)
      CALL GETORDERAS(IOAS)
      CALL GETLAM4(0,QLAM4)
      CALL GETLAM5(0,QLAM5)
      WRITE(*,*) "ALLUNC: The PDF set has",NPDF+1," members"
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
      LONE  = NPDF.LE.1
      LSER  = .NOT.LONE.AND.NPDF.LT.10.AND..NOT.LSTAT.AND..NOT.LALG
      LPDF  = .NOT.LONE.AND..NOT.LSER
      IF (LONE) THEN
         WRITE(*,*)"ALLUNC: Only central PDF available."
      ENDIF
      IF (LPDF) THEN
         WRITE(*,*)"ALLUNC: Deriving PDF uncertainties."
      ENDIF
      IF (LSTAT) THEN
         WRITE(*,*)"ALLUNC: Deriving statistical uncertainties."
      ENDIF
      IF (LSER) THEN
         WRITE(*,*)"ALLUNC: Deriving cross sections of series variation"
      ENDIF
      IF (LSCL) THEN
         WRITE(*,*)"ALLUNC: Deriving scale uncertainties"
      ENDIF
      IF (LRAT) THEN
         WRITE(*,*)"ALLUNC: Deriving uncertainties for ratios"
      ENDIF
      IF (LALG) THEN
         WRITE(*,*)"ALLUNC: Deriving algorithmic uncertainties."
      ENDIF
      WRITE(*,*)" "


      
c - One initial call - to fill commonblock -> for histo-booking
c - Use primary table for this (recall: ref. table has 2 x rap. bins)
      FILENAME = TABPATH(1:LENOCC(TABPATH))//"/"//TABNAME
      CALL FX9999CC(FILENAME,1D0,1D0,0,XSECT1)
      CALL PDFHIST(1,HISTFILE,LONE,LPDF,LSTAT,LALG,LSER,NPDF,LRAT,LSCL)
      WRITE(*,*)"ALLUNC: The observable has",NBINTOT," bins -",
     >     NSUBPROC," subprocesses"
      
      
      
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

c - Compute PDF uncertainties for all available scales
c - Check that FILENAME is still the primary table here ...!!!
      IF (LPDF.AND..NOT.LSER) THEN
         WRITE(*,*)"========================================"//
     >        "================================"
         WRITE(*,*)"Relative PDF Uncertainties"
         WRITE(*,*)"(the printed values are for the total)"
         WRITE(*,*)"(histograms contain results for all orders and "//
     >        "subprocesses)"
         WRITE(*,*)" bin       cross section           "//
     >        "lower PDF uncertainty   upper PDF uncertainty"
         DO I=1,NSCALEVAR
            CALL INITPDF(0)
            MUR = MURSCALE(I)
            MUF = MUFSCALE(I)
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            WRITE(*,*)"ALLUNC: Now scale no.",i,"; mur, muf = ",mur,muf
            WRITE(*,*)"----------------------------------------"//
     >           "--------------------------------"
            CALL FX9999CC(FILENAME,MUR,MUF,0,XSECT1)
            
            ISTAT = 1
            IMODE = 1
            NRAP  = NRAPIDITY
            IF (LTOY) THEN
               IMODE = 2
            ENDIF
            IF (LRAT) THEN
               IMODE = 5
               NRAP = 2*NRAPIDITY
               write(*,*)"nrapidity,nrap",nrapidity,nrap
            ENDIF
            CALL CENRES(IMODE)
            CALL UNCERT(ISTAT,IMODE)

ckr Do loop runs once even if NPDF=0! => Avoid with IF statement
            IF (LPDF) THEN
               DO J=1,NPDF
                  CALL INITPDF(J)
                  CALL FX9999CC(FILENAME,MUR,MUF,0,XSECT0)
                  ISTAT = 2
                  CALL UNCERT(ISTAT,IMODE)
               ENDDO
            ENDIF
            
            ISTAT = 3
            CALL UNCERT(ISTAT,IMODE)

c - Give some standard output, fill histograms
            IBIN = 0
            write(*,*)"nrapidity,nrap",nrapidity,nrap
            DO IRAP=1,NRAP
               write(*,*)"irap,npt",irap,npt(irap)
               DO IPT=1,NPT(IRAP)
                  IBIN = IBIN+1
                  WRITE(*,900) IBIN,MYRES(IBIN,NSUBPROC+1,NORD+1),
     >                 WTDXL2(IBIN,NSUBPROC+1,NORD+1)-1D0,
     >                 WTDXU2(IBIN,NSUBPROC+1,NORD+1)-1D0
                  WRITE(*,900) IBIN,WTX(IBIN,NSUBPROC+1,NORD+1),
     >                 WTDXL2(IBIN,NSUBPROC+1,NORD+1)-1D0,
     >                 WTDXU2(IBIN,NSUBPROC+1,NORD+1)-1D0
               ENDDO
            ENDDO
ckr 900     FORMAT(1P,I5,3(3X,E21.14))
 900        FORMAT(1P,I5,3(6X,E18.11))

c - Fill histograms
            CALL PDFFILL(NRAP,0,I,MYRES)
            CALL PDFFILL(NRAP,1,I,WTDXL2)
            CALL PDFFILL(NRAP,2,I,WTDXU2)
         ENDDO                  ! Loop over scales
      ENDIF



c - Statistics part
c - Use statistics tables
c - Call statistical error-code for scenario
      IF (LSTAT) THEN
         WRITE(*,*)"ALLUNC: Evaluating statistical uncertainties"
         CALL STATCODE(TABPATH,SCENARIO,BORNN,NLON)
      ENDIF
      
      
      
c - PDF series of variations (e.g. alpha_s(M_Z))
c - Use primary table
c - Check that FILENAME is still the primary table here ...!!!
c - Fill only central scale
c - (ISCALE=3 in FORTRAN, refscale=2 in C++ parlance of author code)
      IF (LSER) THEN
         ISCALE = 3
         MUR = MURSCALE(ISCALE)
         MUF = MUFSCALE(ISCALE)
         WRITE(*,*)"ALLUNC: For PDF series fill only scale no.",ISCALE,
     >        "; mur, muf = ",mur,muf

         CALL INITPDF(0)
         CALL FX9999CC(FILENAME,MUR,MUF,0,XSECT0)
         ISTAT = 1
         IMODE = 3
         CALL CENRES(IMODE)
         CALL UNCERT(ISTAT,IMODE)
         ISTAT = 2
         DO J=1,NPDF
            CALL INITPDF(J)
            CALL FX9999CC(FILENAME,MUR,MUF,0,XSECT0)
            CALL UNCERT(ISTAT,IMODE)
         ENDDO
         ISTAT = 3
         CALL UNCERT(ISTAT,IMODE)

c - Give some standard output, fill histograms
         WRITE(*,*)"========================================"//
     >        "================================"
         WRITE(*,*)"Relative Series Uncertainties"
         WRITE(*,*)"(the printed values are for the total)"
         WRITE(*,*)"(histograms contain results for all orders and "//
     >        "subprocesses)"
         WRITE(*,*)" bin       cross section           "//
     >        "lower ser. uncertainty  upper ser. uncertainty"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         WRITE(*,*)"ALLUNC: Uncertainties from"//
     >        " series variations"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         ISCALE = 3
         IBIN   = 0
         DO IRAP=1,NRAPIDITY
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               WRITE(*,900) IBIN,MYRES(IBIN,NSUBPROC+1,NORD+1),
     >              WTDXLM(IBIN,NSUBPROC+1,NORD+1)-1D0,
     >              WTDXUM(IBIN,NSUBPROC+1,NORD+1)-1D0
            ENDDO
         ENDDO
         CALL PDFFILL(NRAPIDITY,1,ISCALE,WTDXLM)
         CALL PDFFILL(NRAPIDITY,2,ISCALE,WTDXUM)
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
c - (ISCALE=3 in FORTRAN, refscale=2 in C++ parlance of author code)
      IF (LSCL) THEN
         NSCALES = NSCALEVAR
         IF (NSCALEMAX.GE.8.AND.NSCALEVAR.EQ.4) THEN
            NSCALES = 8
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
         ISCALE = 3
         MUR = MURSCALE(ISCALE)
         MUF = MUFSCALE(ISCALE)
         CALL FX9999CC(FILENAME,MUR,MUF,0,XSECT0)
         ISTAT = 1
         IMODE = 3
         CALL CENRES(IMODE)
         CALL UNCERT(ISTAT,IMODE)
         ISTAT = 2
         DO ISCALE=1,NSCALES
ckr Do neither use scale 1 with factor of 1/4 nor default scale 3
ckr Ugly goto construction avoidable with f90 CYCLE command
            IF (ISCALE.EQ.1.OR.ISCALE.EQ.3) GOTO 10
            MUR = MURSCALE(ISCALE)
            MUF = MUFSCALE(ISCALE)
            CALL FX9999CC(FILENAME,MUR,MUF,0,XSECT0)
            CALL UNCERT(ISTAT,IMODE)
 10         CONTINUE
         ENDDO
         ISTAT = 3
         CALL UNCERT(ISTAT,IMODE)

c - Give some standard output, fill histograms
         WRITE(*,*)"========================================"//
     >        "================================"
         WRITE(*,*)"Relative Scale Uncertainties"
         WRITE(*,*)"(the printed values are for the total)"
         WRITE(*,*)"(histograms contain results for all orders and "//
     >        "subprocesses)"
         WRITE(*,*)" bin       cross section           "//
     >        "lower scale uncertainty upper scale uncertainty"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         WRITE(*,*)"ALLUNC: Uncertainties from",NSCALES,
     >        " scale variations"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         IORD   = 0
         ISCALE = 3
         ISUB   = 0
         IBIN   = 0
         DO IRAP=1,NRAPIDITY
            IHIST = IORD*1000000 + ISCALE*100000 +
     >           ISUB*10000 + IRAP*100
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               PT(IBIN) = REAL(PTBIN(IRAP,IPT))
               WRITE(*,900) IBIN,MYRES(IBIN,NSUBPROC+1,NORD+1),
     >              WTDXLM(IBIN,NSUBPROC+1,NORD+1)-1D0,
     >              WTDXUM(IBIN,NSUBPROC+1,NORD+1)-1D0
            ENDDO
         ENDDO
         CALL PDFFILL(NRAPIDITY,6,ISCALE,WTDXLM)
         CALL PDFFILL(NRAPIDITY,7,ISCALE,WTDXUM)
      ENDIF



c - Algorithmic part
c - Use reference table
c - Reference scale is always no. 1
c - Reference result is in nrap/2 ++ bins
c - Reference result is evaluated with CTEQ61 PDFs
c - Default scale (C++ 2, Fortran 3) ==> normal result
c - (Use other ISCALE in FORTRAN ONLY if refscale <> 2 in author code)
c - Attention: From now on ref. table loaded ==> rap. bins doubled
      IF (LALG) THEN
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
         ISCALE = 1
         MUR = MURSCALE(ISCALE)
         MUF = MUFSCALE(ISCALE)
         CALL FX9999CC(FILENAME,MUR,MUF,1,XSECT)
         ISTAT = 1
         IMODE = 4
         CALL CENRES(IMODE)
         CALL UNCERT(ISTAT,IMODE)

c - Normal result
         ISCALE = 3
         MUR = MURSCALE(ISCALE)
         MUF = MUFSCALE(ISCALE)
         CALL FX9999CC(FILENAME,MUR,MUF,1,XSECT)
         ISTAT = 2
         CALL UNCERT(ISTAT,IMODE)
         ISTAT = 3
         CALL UNCERT(ISTAT,IMODE)

c - Give some standard output, fill histograms
         WRITE(*,*)"========================================"//
     >        "================================"
         WRITE(*,*)"Relative Algorithmic Uncertainty"
         WRITE(*,*)"(the printed values are for the total)"
         WRITE(*,*)"(histograms contain results for all orders and "//
     >        "subprocesses)"
         WRITE(*,*)" bin       cross section           "//
     >        "algorithmic uncertainty"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         WRITE(*,*)"ALLUNC: Algorithmic Uncertainty with "//
     >        "reference to CTEQ61 PDF set"
         WRITE(*,*)"----------------------------------------"//
     >        "--------------------------------"
         ISCALE = 3
         IBIN   = 0
         DO IRAP=1,INT(NRAPIDITY/2)
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               WRITE(*,901) IBIN,WTX(IBIN,NSUBPROC+1,NORD+1),
     >              WTDXUM(IBIN,NSUBPROC+1,NORD+1)-1D0
 901           FORMAT(1P,I5,2(6X,E18.11))
            ENDDO
         ENDDO
         CALL PDFFILL(NRAPIDITY,5,ISCALE,WTDXUM)
      ENDIF



c - Close hbook file
      CALL PDFHIST(2,HISTFILE,LONE,LPDF,LSTAT,LALG,LSER,NPDF,LRAT,LSCL)
      END



c
c ======================= Book the histograms ========================
c
      SUBROUTINE PDFHIST(N,HISTFILE,
     &     LONE,LPDF,LSTAT,LALG,LSER,NPDF,LRAT,LSCL)
      IMPLICIT NONE
      CHARACTER*(*) HISTFILE
      CHARACTER*255 CSTRNG,CBASE1,CBASE2,CTMP
      INTEGER N,LENOCC,IPDF,NPDF,IPTMAX,NRAP
      LOGICAL LONE,LPDF,LSTAT,LALG,LSER,LRAT,LSCL

      INTEGER I,J,ISTAT2,ICYCLE,NSCALES
      INTEGER IORD,ISUB,ISCALE,IRAP,IPT,IHIST,NHIST
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
         WRITE(*,*)"----------- Book histograms -------"
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
            NSCALES = NSCALEVAR
            IF (NSCALEMAX.GE.8.AND.NSCALEVAR.EQ.4) THEN
               NSCALES = 8
            ENDIF
            DO ISCALE=1,NSCALES ! Scale variations
               DO ISUB=0,NSUBPROC ! Subprocesses: 0 tot + 7 subproc
                  DO IRAP=1, NRAP
                     IHIST = IORD*1000000 + ISCALE*100000 +
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
                     WRITE(CTMP,'(I1)'),ISCALE
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
                     IF (LSER) THEN
                        DO IPDF=1,NPDF
                           CALL HBOOKB(IHIST+IPDF,
     >                          CSTRNG(1:LENOCC(CSTRNG)),
     >                          NPT(IRAP),PT,0)
                           NHIST = NHIST+1
                        ENDDO
                     ENDIF
                     IF (LPDF) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_dPDF_low"
                        CALL HBOOKB(IHIST+1,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_dPDF_up"
                        CALL HBOOKB(IHIST+2,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+2
ckr                        write(*,*)"2. Booked histo #",nhist
                     ENDIF
                     IF (LSTAT.AND.IORD.LE.2.AND.
     >                    ISCALE.EQ.3.AND.ISUB.EQ.0) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_dstat/xsect_%"
                        CALL HBOOKB(IHIST+3,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_dmax/2/xsect_%"
                        CALL HBOOKB(IHIST+4,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+2
ckr                        write(*,*)"3. Booked histo #",nhist
                     ENDIF
                     IF (LALG.AND.IORD.LE.2.AND.
     >                    ISCALE.EQ.3) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_dref/xsect_%"
                        CALL HBOOKB(IHIST+5,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
ckr                        write(*,*)"4. Booked histo #",nhist
                     ENDIF
                     IF (LSCL.AND.IORD.LE.2.AND.
     >                    ISCALE.EQ.3) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_dscl_low/xsect_%"
                        CALL HBOOKB(IHIST+6,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_dscl_up/xsect_%"
                        CALL HBOOKB(IHIST+7,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
ckr                        write(*,*)"4b. Booked histo #",nhist
                     ENDIF
                     IF (LSTAT.AND.ISUB.EQ.0) THEN
                        IHIST = IORD*1000000 + ISCALE*100000 +
     >                       ISUB*10000 + IRAP*100
                        CSTRNG = CBASE1(1:LENOCC(CBASE1))
                        WRITE(CTMP,'(I1)'),IORD
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//"_iord="
     >                       //CTMP
                        WRITE(CTMP,'(I1)'),ISCALE
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
                           WRITE(CTMP,'(I1)'),ISCALE
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
            ENDDO               ! End od ISCALE loop, IORD still on
         ENDDO
         WRITE(*,*)"Number of histograms booked:",NHIST
         WRITE(*,*)"-----------------------------------"



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
      SUBROUTINE PDFFILL(NRAP,IOFF,ISCALE,DVAL)
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      INTEGER NRAP,IOFF,ISCALE
      DOUBLE PRECISION DVAL(NBINTOTMAX,NMAXSUBPROC+1,4)
      INTEGER I,J,IBIN,IORD,IORD2,ISUB,ISUB2,IHIST
      REAL RVAL
      
      IF (ISCALE.LT.1 .OR. ISCALE.GT.NSCALEVAR) THEN
         WRITE(*,*) "\nPDFFILL: ERROR! ISCALE ",ISCALE,
     >        " is out of range, aborted!"
         WRITE(*,*) "PDFFILL: Max. ISCALE: ",NSCALEVAR
         STOP
      ENDIF

c - Fill all histograms for the given scale
c - Fill sums over subprocesses and/or orders into zero factors for IHIST
      DO IORD2=1,NORD+1         ! Order: LO, NLO-corr, NNLO-corr, ... , tot
         IORD = IORD2
         IF (IORD2.EQ.NORD+1) IORD = 0 
         DO ISUB2=1,NSUBPROC+1  ! Subprocesses hh: 7 subprocesses + total
            ISUB = ISUB2
            IF (ISUB2.EQ.NSUBPROC+1) ISUB=0
            IBIN=0
            DO I=1,NRAP
               DO J=1,NPT(I)
                  IBIN = IBIN + 1
ckr Recall: HBOOK understands only single precision
                  RVAL  = REAL(DVAL(IBIN,ISUB2,IORD2))
                  IHIST = IORD*1000000+ISCALE*100000+ISUB*10000+I*100
                  IHIST = IHIST+IOFF
                  CALL HFILL(IHIST,REAL(PTBIN(I,J)+0.01),0.0,RVAL)
               ENDDO            ! pT-loop
            ENDDO               ! rap-loop
         ENDDO                  ! isub-loop
      ENDDO                     ! iord-loop
      
      RETURN
      END
