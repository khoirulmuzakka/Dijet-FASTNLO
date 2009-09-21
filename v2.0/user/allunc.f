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
      CHARACTER*255 HISTFILE
      CHARACTER*255 SCENARIO,FILENAME,TABPATH,TABNAME,REFNAME
      CHARACTER*255 PDFSET,PDFPATH,LHAPDF,CHTMP
      CHARACTER*8 CH8TMP
      CHARACTER*4 CH4TMP
      INTEGER BORNN,NLON,LENOCC
      INTEGER I,J,L1,L2,L3,L4,MYNPDF,IOPDF,IOAS
      INTEGER ISTAT,MYISCALE,IORD,IBIN,NBIN,ISUB,IRAP,IPT,IHIST
      LOGICAL LONE,LALG,LSTAT,LPDF,LSER,LTOY
      DOUBLE PRECISION MUR,MUF,DIFF,SUMM,QLAM4,QLAM5
      DOUBLE PRECISION
cnew5     >     RES0(MXOBSBIN,MXSUBPROC+1,0:3),
cnew5     >     RES1HI(MXOBSBIN,MXSUBPROC+1,0:3),
cnew5     >     RES1LO(MXOBSBIN,MXSUBPROC+1,0:3),
cnew5     >     RESLO,RESHI,DREF
     >     RES0  (MXOBSBIN,0:MXSUBPROC,0:MXCTRB),
     >     RES1HI(MXOBSBIN,0:MXSUBPROC,0:MXCTRB),
     >     RES1LO(MXOBSBIN,0:MXSUBPROC,0:MXCTRB),
     >     RESLO,RESHI,DREF
c - NNPDF method
Comment:       DOUBLE PRECISION WGT(MXOBSBIN,MXSUBPROC+1,0:3)
Comment:       DOUBLE PRECISION WGT2(MXOBSBIN,MXSUBPROC+1,0:3)
Comment:       DOUBLE PRECISION WGTX(MXOBSBIN,MXSUBPROC+1,0:3)
Comment:       DOUBLE PRECISION WGTX2(MXOBSBIN,MXSUBPROC+1,0:3)
      DOUBLE PRECISION WGT(MXOBSBIN,0:MXSUBPROC,0:MXCTRB)
      DOUBLE PRECISION WGT2(MXOBSBIN,0:MXSUBPROC,0:MXCTRB)
      DOUBLE PRECISION WGTX(MXOBSBIN,0:MXSUBPROC,0:MXCTRB)
      DOUBLE PRECISION WGTX2(MXOBSBIN,0:MXSUBPROC,0:MXCTRB)
c - To unify quoted uncertainties (CL68,CL90,special)
      DOUBLE PRECISION CL90,CLGJR
c - Attention!!! This must be declared consistent with the
c                definition in the commonblock!!!!!
      DOUBLE PRECISION XS0(MXOBSBIN,0:MxCtrb),XS1(MXOBSBIN,0:MxCtrb)
      DOUBLE PRECISION XS2(MXOBSBIN,0:MxCtrb)
      REAL PT(MXOBSBIN)

      CHARACTER*255 ASMODE
      DOUBLE PRECISION ASMZVAL
      INTEGER IASLOOP
      COMMON/STEER/ASMZVAL,IASLOOP,ASMODE

c --- Mods for v2.0
      INTEGER IOBS,IPROC,ICP
      INTEGER IC,MYICONTR,NORD

      MYICONTR = 2

c --- Parse command line
      WRITE(*,*)"\n #################################################"
      WRITE(*,*)"# ALLUNC"
      WRITE(*,*)"#################################################"
      WRITE(*,*)"# Program to derive the algorithmic, statistical, "
      WRITE(*,*)"# and PDF uncertainties using fastNLO tables"
      WRITE(*,*)"#################################################"
      WRITE(*,*)"#"

*---Scenario
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
     >           'share/lhapdf/PDFsets'
            WRITE(*,*)'  alpha_s calc., def. = PDF (from PDF set)'
            WRITE(*,*)'    alt. = PY: 0-, 1- and 2-loop '//
     >           '(from Pythia 6.4 using Lambda_5 from PDF)'
            WRITE(*,*)'    alt. = KR: 1-, 2- and 3-loop '//
     >           '(from hep-ph/9506442)'
            WRITE(*,*)'    alt. = MW: 2- and 4-loop '//
     >           '(from hep-ph/9806404)'
            WRITE(*,*)'  alpha_s(M_Z), def. from PDF set'
            WRITE(*,*)'     (in mode PY this has to be Lambda_5!)'
            WRITE(*,*)'  alpha_s loop order, def. from PDF set'
            WRITE(*,*)' '
            STOP
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
         PDFPATH = "/share/lhapdf/PDFsets"
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
      CALL NUMBERPDF(MYNPDF)
      CALL GETORDERPDF(IOPDF)
      CALL GETORDERAS(IOAS)
      CALL GETLAM4(0,QLAM4)
      CALL GETLAM5(0,QLAM5)
      WRITE(*,*) "ALLUNC: The PDF set has",MYNPDF+1," members"
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
ckr      LALG  = .TRUE.
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
      LONE  = MYNPDF.LE.1
      LSER  = .NOT.LONE.AND.MYNPDF.LT.10.AND..NOT.LSTAT.AND..NOT.LALG
      LPDF  = .NOT.LONE.AND..NOT.LSER
ctmp
      lalg  = .false.
ctmp      lpdf  = .false.
      lser  = .false.
      lstat = .false.
ctmp
      IF (LONE) THEN
         WRITE(*,*)"ALLUNC: Only central PDF available."
      ENDIF
      IF (LALG) THEN
         WRITE(*,*)"ALLUNC: Deriving algorithmic uncertainties."
      ENDIF
      IF (LSTAT) THEN
         WRITE(*,*)"ALLUNC: Deriving statistical uncertainties."
      ENDIF
      IF (LPDF) THEN
         WRITE(*,*)"ALLUNC: Deriving PDF uncertainties."
      ENDIF
      IF (LSER) THEN
         WRITE(*,*)"ALLUNC: Deriving cross sections of series variation"
      ENDIF
      WRITE(*,*)" "


      
c - One initial call - to fill commonblock -> for histo-booking
c - Use primary table for this (recall: ref. table has 2 x rap. bins)
      FILENAME = TABPATH(1:LENOCC(TABPATH))//"/"//TABNAME
ckr Markus selection scheme
ckr      call FNSET("P_ORDPTHY",1) ! select order pert. theory: 1=LO, 2=NLO
ckr Select order pert. theory: 1=LO, 2=NLO
      CALL FNSET("P_ORDPTHY",2)
ckr Select no. of loops in threshold corrections
      CALL FNSET("P_THRESHCOR",0)
      CALL FX9999CC(FILENAME,1D0,1D0,0,XS1)

ckr Add some initializations after first table read
      NORD   = NPOW(MYICONTR) - ILOORD +1
      write(*,*)"icontr,npow,iloord,nord",icontr,npow(icontr),iloord
     >     ,nord
ckr Determine npt,nrap,rapindex; attention, simple workaround,
ckr works only for two-dim. pt,y or y,pt binning!
      IF (NDimCounter(1).EQ.NObsBin) Then ! 1st dim.: pT, 2nd: rap
         IPTDIM  = 1
         IRAPDIM = 2
      ELSEIF (NDimCounter(2).EQ.NObsBin) Then ! 1st dim.: rap, 2nd: pT
         IPTDIM  = 2
         IRAPDIM = 1
      ELSE
         WRITE(*,*)"ALLUNC: Cannot decrypt pt, rap binning, stopped!"
         STOP
      ENDIF
      NRAPIDITY = NDimCounter(IRAPDIM)
      DO I=1,NRAPIDITY-1
         NPT(I) = IDimPointer(I+1,IRAPDIM) - IDimPointer(I,IRAPDIM)
      ENDDO
      NPT(NRAPIDITY) = NDimCounter(IPTDIM) -
     >     IDimPointer(NRAPIDITY,IRAPDIM) + 1
      CALL PDFHIST(1,HISTFILE,LONE,LPDF,LSTAT,LALG,LSER,MYNPDF)



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
      IF (.NOT.LSER) THEN
ckr      DO I=1,NSCALEVAR(MYICONTR,NSCALEDIM(MYICONTR))
      DO I=3,3
         CALL INITPDF(0)
         MUR = SCALEFAC(MYICONTR,NSCALEDIM(MYICONTR),I)
         MUF = SCALEFAC(MYICONTR,NSCALEDIM(MYICONTR),I)
         WRITE(*,*)"ALLUNC: Now scale no.",i,"; mur, muf = ",mur,muf
         CALL FX9999CC(FILENAME,MUR,MUF,0,XS1)

c - Save the result array from the first call (= central result)
c   and reset result arrays   
         myicontr = 1 
         WRITE(*,*)"ALLUNC: The observable has",NOBSBIN," bins -",
     >        NSUBPROC(MYICONTR)," subprocesses"
         
         DO IOBS=1,NOBSBIN
            DO ICONTR=0,NContrib
               ICP = IContrPointer(ICONTR)
               DO IPROC = 0,NSUBPROC(ICP)
                  RES0(IOBS,IPROC,ICONTR)   = RESULT(IOBS,IPROC,ICP)
                  RES1LO(IOBS,IPROC,ICONTR) = 0D0
                  RES1HI(IOBS,IPROC,ICONTR) = 0D0
                  WGT(IOBS,IPROC,ICONTR)    = 0D0
                  WGT2(IOBS,IPROC,ICONTR)   = 0D0
                  WGTX(IOBS,IPROC,ICONTR)   = 0D0
                  WGTX2(IOBS,IPROC,ICONTR)  = 0D0
               ENDDO
cdebug               WRITE(*,*)"ALLUNC: IOBS,ICONTR,XS: ",
cdebug     >              IOBS,ICONTR,RES0(IOBS,0,ICONTR)
            ENDDO
         ENDDO
         
         
ckr Do loop runs once even if MYNPDF=0! => Avoid with IF statement
ckr         IF (MYNPDF.GT.1) THEN
         IF (LPDF) THEN
*---Convert from CL68 to CL90 values
            CL90  = 1.64485D0   ! SQRT(2.D0)/InvERF(0.9D0)
*---Convert from GJR to CTEQ (CL90) values
            CLGJR = 1.D0/0.47D0
            DO J=1,MYNPDF
               CALL INITPDF(J)
               CALL FX9999CC(FILENAME,MUR,MUF,0,XS0)
c - For all bins/subproc/orders: Add negative/positive variations
               DO IOBS=1,NOBSBIN
                  DO ICONTR=0,NContrib
                     ICP = IContrPointer(ICONTR)
                     DO IPROC = 0,NSUBPROC(ICP)
                        DIFF = RESULT(IOBS,IPROC,ICP) -
     >                       RES0(IOBS,IPROC,ICONTR)
                        IF (DIFF.GT.0D0) THEN
                           RES1HI(IOBS,IPROC,ICONTR) =
     >                          RES1HI(IOBS,IPROC,ICONTR)+DIFF*DIFF
                        ELSE
                           RES1LO(IOBS,IPROC,ICONTR) =
     >                          RES1LO(IOBS,IPROC,ICONTR)+DIFF*DIFF
                        ENDIF
                        SUMM = RESULT(IOBS,IPROC,ICP)
                        WGT(IOBS,IPROC,ICONTR)   =
     >                       WGT(IOBS,IPROC,ICONTR)   + 1.D0 
                        WGT2(IOBS,IPROC,ICONTR)  =
     >                       WGT2(IOBS,IPROC,ICONTR)  + 1.D0 
                        WGTX(IOBS,IPROC,ICONTR)  =
     >                       WGTX(IOBS,IPROC,ICONTR)  + SUMM 
                        WGTX2(IOBS,IPROC,ICONTR) =
     >                       WGTX2(IOBS,IPROC,ICONTR) + SUMM*SUMM 
                     ENDDO
                  ENDDO
               ENDDO
Comment:                DO L1=1,NOBSBIN
Comment:                   DO L2=1,(NSUBPROC(MYICONTR)+1)
Comment:                      DO L3=1,NORD
Comment:                         SUMM = 0.D0
Comment:                         DIFF = - RES0(L1,L2,L3)
Comment:                         DO L4=1,L3
Comment:                            SUMM = SUMM + RESULT(L1,L2,L4)
Comment:                            DIFF = DIFF + RESULT(L1,L2,L4) 
Comment:                         ENDDO
Comment:                         WGT(L1,L2,L3)   = WGT(L1,L2,L3)   + 1.D0 
Comment:                         WGT2(L1,L2,L3)  = WGT2(L1,L2,L3)  + 1.D0 
Comment:                         WGTX(L1,L2,L3)  = WGTX(L1,L2,L3)  + SUMM 
Comment:                         WGTX2(L1,L2,L3) = WGTX2(L1,L2,L3) + SUMM*SUMM 
Comment:                         IF (DIFF.GT.0D0) THEN
Comment:                            RES1HI(L1,L2,L3) =
Comment:      >                          RES1HI(L1,L2,L3)+DIFF*DIFF
Comment:                         ELSE
Comment:                            RES1LO(L1,L2,L3) =
Comment:      >                          RES1LO(L1,L2,L3)+DIFF*DIFF
Comment:                         ENDIF
Comment:                      ENDDO
Comment:                   ENDDO
Comment:                ENDDO
            ENDDO               ! Loop over bins
         ENDIF                  ! Not done for npdf <= 1

c - Take square-root of sum of squares
c - or apply Toy MC method for NNPDF
         DO IOBS=1,NOBSBIN
            IF (LPDF) THEN
               DO ICONTR=0,NContrib
                  ICP = IContrPointer(ICONTR)
                  DO IPROC = 0,NSUBPROC(ICP)
                     IF (.NOT.LTOY) THEN
c     kr                        RES1HI(IOBS,IPROC,ICONTR) = CL90
c     *SQRT(RES1HI(IOBS,IPROC,ICONTR))
c     kr                        RES1LO(IOBS,IPROC,ICONTR) = -CL90
c     *SQRT(RES1LO(IOBS,IPROC,ICONTR))
c     kr                        RES1HI(IOBS,IPROC,ICONTR) = CLGJR
c     *SQRT(RES1HI(IOBS,IPROC,ICONTR))
c     kr                        RES1LO(IOBS,IPROC,ICONTR) = -CLGJR
c     *SQRT(RES1LO(IOBS,IPROC,ICONTR))
                        RES1HI(IOBS,IPROC,ICONTR) =  SQRT(RES1HI(IOBS
     >                       ,IPROC,ICONTR))
                        RES1LO(IOBS,IPROC,ICONTR) = -SQRT(RES1LO(IOBS
     >                       ,IPROC,ICONTR))
                     ELSE
                        RES0(IOBS,IPROC,ICONTR) = WGTX(IOBS,IPROC
     >                       ,ICONTR)/WGT(IOBS,IPROC,ICONTR)
                        RES1HI(IOBS,IPROC,ICONTR) = 
     >                       (WGTX2(IOBS,IPROC,ICONTR)/WGT(IOBS,IPROC
     >                       ,ICONTR) -RES0(IOBS,IPROC,ICONTR)*RES0(IOBS
     >                       ,IPROC,ICONTR))
c     kr                        RES1HI(IOBS,IPROC,ICONTR) = CL90
c     *SQRT(RES1HI(IOBS,IPROC,ICONTR))
                        RES1HI(IOBS,IPROC,ICONTR) = SQRT(RES1HI(IOBS
     >                       ,IPROC,ICONTR))
                        RES1LO(IOBS,IPROC,ICONTR) = -RES1HI(IOBS,IPROC
     >                       ,ICONTR)
                     ENDIF
                  ENDDO
               ENDDO
               RESLO = RES1LO(IOBS,0,0)/RES0(IOBS,0,0)
               RESHI = RES1HI(IOBS,0,0)/RES0(IOBS,0,0)
ckr 30.01.2008: Change output format for better comp. with C++ version
            ELSE
               RESLO = 0D0
               RESHI = 0D0
            ENDIF
            WRITE(*,900) IOBS,RES0(IOBS,0,0),RESLO,RESHI
         ENDDO
ckr 900     FORMAT(1P,I5,3(3X,E21.14))
 900     FORMAT(1P,I5,3(6X,E18.11))

c - Fill histograms
         CALL PDFFILL(I,RES0,RES1HI,RES1LO)
         
      ENDDO                     ! Loop over scales

      WRITE(*,*)"Bin    x-sect       lower PDF    upper PDF unc."
      WRITE(*,*)"(the printed uncertainties are for "//
     >     "the highest order)"
      WRITE(*,*)"(histograms contain results for all orders and "//
     >     "subprocesses)"
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
c - (MYISCALE=3 in FORTRAN, refscale=2 in C++ parlance of author code)
      IF (LSER) THEN
         MYISCALE = 3
         MUR = SCALEFAC(MYICONTR,NSCALEDIM(MYICONTR),MYISCALE)
         MUF = SCALEFAC(MYICONTR,NSCALEDIM(MYICONTR),MYISCALE)
         WRITE(*,*)"ALLUNC: For PDF series fill only scale no.",MYISCALE,
     >        "; mur, muf = ",mur,muf
         DO J=0,MYNPDF
            CALL INITPDF(J)
            CALL FX9999CC(FILENAME,MUR,MUF,0,XS0)
            IBIN = 0
            DO IRAP=1,NRAPIDITY
               DO IPT=1,NPT(IRAP)
                  IBIN = IBIN+1
ckr                  PT(IBIN) = REAL(PTBIN(IRAP,IPT))
                  PT(IBIN) = REAL(LOBIN(IBIN,IPTDIM))
                  DO ISUB=0,NSUBPROC(MYICONTR)
                     RES0(IBIN,ISUB+1,0) = 0.D0
                     DO IORD=1,NORD
                        RES0(IBIN,ISUB+1,IORD) = XS0(IBIN,IORD)
                        RES0(IBIN,ISUB+1,0) = RES0(IBIN,ISUB+1,0) +
     >                       RES0(IBIN,ISUB+1,IORD)
                     ENDDO
c - Fill histograms
                     DO IORD=0,NORD
                        IHIST = IORD*1000000 + MYISCALE*100000 +
     >                       ISUB*10000 + IRAP*100
                        CALL HFILL(IHIST+J,PT(IBIN),0.,
     >                       REAL(RES0(IBIN,ISUB+1,IORD)))
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            NBIN = IBIN
         ENDDO
      ENDIF



c - Algorithmic part
c - Use reference table
      IF (LALG) THEN
         FILENAME = TABPATH(1:LENOCC(TABPATH))//"/"//REFNAME
         WRITE(*,*)"ALLUNC: Taking reference table: "//
     >        FILENAME(1:LENOCC(FILENAME))
c - Initialize CTEQ61 reference PDFs
         PDFSET = PDFPATH(1:LENOCC(PDFPATH))//"/cteq61.LHgrid"
         WRITE(*,*)"ALLUNC: Taking reference PDF: "//
     >        PDFSET(1:LENOCC(PDFSET))
         CALL INITPDFSET(PDFSET(1:LENOCC(PDFSET)))
         CALL INITPDF(0)

c - Default scale (C++ 2, Fortran 3) ==> normal result
         MYISCALE = 3
c - (Use other MYISCALE in FORTRAN ONLY if refscale <> 2 in author code)
c         MYISCALE = 2
         MUR = SCALEFAC(MYICONTR,NSCALEDIM(MYICONTR),MYISCALE)
         MUF = SCALEFAC(MYICONTR,NSCALEDIM(MYICONTR),MYISCALE)
         CALL FX9999CC(FILENAME,MUR,MUF,1,XS2)
c - Attention: From now on ref. table loaded ==> rap. bins doubled
c - Get the normal result, scale 3, from first half of doubled rap bins
         IBIN = 0
         DO IRAP=1,INT(NRAPIDITY/2)
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               DO ISUB=0,NSUBPROC(MYICONTR)
                  RES0(IBIN,ISUB+1,0) = 0.D0
                  DO IORD=1,NORD
                     RES0(IBIN,ISUB+1,IORD) = XS2(IBIN,IORD)
                     RES0(IBIN,ISUB+1,0) = RES0(IBIN,ISUB+1,0) +
     >                    RES0(IBIN,ISUB+1,IORD)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         NBIN = IBIN

c - Reference scale no. always 1
         MYISCALE = 1
         MUR = SCALEFAC(MYICONTR,NSCALEDIM(MYICONTR),MYISCALE)
         MUF = SCALEFAC(MYICONTR,NSCALEDIM(MYICONTR),MYISCALE)
         CALL FX9999CC(FILENAME,MUR,MUF,1,XS2)
c - Get the reference result, scale 1, nrap/2 ++ bins
         DO IRAP=INT(NRAPIDITY/2)+1,NRAPIDITY
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               DO ISUB=0,NSUBPROC(MYICONTR)
                  RES0(IBIN,isub+1,0) = 0.D0
                  DO IORD=1,NORD
                     RES0(IBIN,isub+1,IORD) = XS2(IBIN,IORD)
                     RES0(IBIN,isub+1,0) = RES0(IBIN,isub+1,0) +
     >                    RES0(IBIN,isub+1,IORD)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

c - Compare results and fill histos
         MYISCALE = 3
         ISUB   = 0
         IBIN   = 0
         DO IRAP=1,INT(NRAPIDITY/2)
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
ckr               PT(IBIN) = REAL(PTBIN(IRAP,IPT))
               PT(IBIN) = REAL(LOBIN(IBIN,IPTDIM))
               DO ISUB=0,NSUBPROC(MYICONTR)
                  DO IORD=0,NORD
                     DREF = RES0(IBIN,isub+1,IORD)/
     >                    RES0(IBIN+NBIN,isub+1,IORD) - 1.D0
                     IHIST = IORD*1000000 + MYISCALE*100000 +
     >                    ISUB*10000 + IRAP*100
                     CALL HFILL(IHIST+5,PT(IBIN),0.,REAL(100D0*DREF))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF



c - Close hbook file
      CALL PDFHIST(2,HISTFILE,LONE,LPDF,LSTAT,LALG,LSER,MYNPDF)
      END

c
c ======================= Book the histograms ========================
c
      SUBROUTINE PDFHIST(N,HISTFILE,LONE,LPDF,LSTAT,LALG,LSER,MYNPDF)
      IMPLICIT NONE
      CHARACTER*(*) HISTFILE
      CHARACTER*255 CSTRNG,CBASE1,CBASE2,CTMP
      INTEGER N,LENOCC,IPDF,MYNPDF,IPTMAX
      LOGICAL LONE,LPDF,LSTAT,LALG,LSER

      INTEGER J,ISTAT2,ICYCLE
      INTEGER IORD,ISUB,MYISCALE,IRAP,IPT,IHIST,NHIST
      INCLUDE "fnx9999.inc"
      INCLUDE "strings.inc"
      REAL PT(MXOBSBIN)
      

c - HBOOK common 
      INTEGER NWPAWC
      PARAMETER (NWPAWC=2500000)
      REAL HMEMOR(NWPAWC)
      COMMON /PAWC/ HMEMOR

c --- Mods for v2.0
      INTEGER MYICONTR,NORD

      MYICONTR = 2
      NORD   = NPOW(MYICONTR) - ILOORD + 1
ckr      write(*,*)"PDFHIST: icontr,npow,iloord,nord",icontr,npow(icontr)
ckr     >     ,iloord,nord

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
         
         NHIST = 0
         IF (NSCDESCRIPT.LT.3) THEN
            WRITE(*,*)"\nPDFHIST: ERROR! Insufficient number of "//
     >           "table descriptors for histogram booking: ",
     >           NSCDESCRIPT
            WRITE(*,*)"         Please add descriptive lines, stopped."
            STOP
         ENDIF
         
*--- Process name
ckr         CBASE1 = CIPROC(IPROC)
         IF (NSCDESCRIPT.GT.3) THEN
            CBASE1 = SCDESCRIPT(4)
         ELSEIF (NSCDESCRIPT.GT.2) THEN
            CBASE1 = SCDESCRIPT(3)
         ELSE
            CBASE1 = "Undefined_Process"
         ENDIF
*--- Differential formula
ckr         CBASE1 = CBASE1(1:LENOCC(CBASE1))//"_"
ckr     >        //NAMELABEL(1)
         CBASE1 = CBASE1(1:LENOCC(CBASE1))//"_"
     >        //SCDESCRIPT(1)
*--- Jet algorithm & parameters
ckr         CBASE2 = CBASE1(1:LENOCC(CBASE1))//"_"
ckr     >        //CIALGO(IALGO)
         IF (NSCDESCRIPT.GT.4) THEN
            CBASE2 = CBASE1(1:LENOCC(CBASE1))//"_"//SCDESCRIPT(5)
         ELSE
            CBASE2 = CBASE1(1:LENOCC(CBASE1))//"_Undefined_JetAlgo"
         ENDIF
ckr         CBASE2 = CBASE2(1:LENOCC(CBASE2))//"_"
ckr     >        //CJETRES1(IALGO)
ckr         WRITE(CTMP,'(F3.1)'),JETRES1
ckr         CBASE2 = CBASE2(1:LENOCC(CBASE2))//"="
ckr     >        //CTMP
ckr         write(*,*)"BBB nord",nord
         DO IORD=0,NORD         ! Order: tot, LO, NLO-corr, NNLO-corr
ckr            write(*,*)"BBB scales",NSCALEVAR(MYICONTR,NSCALEDIM(MYICONTR))
            DO MYISCALE=1,NSCALEVAR(MYICONTR,NSCALEDIM(MYICONTR)) ! Scale variations
ckr               write(*,*)"BBB nsubp",NSUBPROC(MYICONTR)
               DO ISUB=0,NSUBPROC(MYICONTR) ! Subprocesses: 0 tot + 7 subproc
ckr                  write(*,*)"BBB nrap",nrapidity
                  DO IRAP=1,NRAPIDITY
                     IHIST = IORD*1000000 + MYISCALE*100000 +
     >                    ISUB*10000 + IRAP*100
ckr                     write(*,*)"ALL: iord,isc,isub,irap,ihist",
ckr     >                    iord,iscale,isub,irap,ihist
                     DO J=1,NPT(IRAP)
ckr                        PT(J) = REAL(PTBIN(IRAP,J))
ckr                        PT(J) = REAL(LOBIN(RAPINDEX(IRAP)+J-1,1))
                        PT(J) = REAL(LOBIN(
     >                       IDimPointer(IRAP,IRAPDIM)+J-1,IPTDIM
     >                       ))
ckr                        write(*,*)"irapdim,iptdim,idimpointer",
ckr     >                       irapdim,iptdim,IDimPointer(IRAP,IRAPDIM)
ckr                        write(*,*)"j,irap,npt(irap),pt",
ckr     >                       j,irap,npt(irap),pt(j)
                     ENDDO
                     PT(NPT(IRAP)+1) = REAL(UPBIN(
     >                    IDimPointer(IRAP,IRAPDIM)+NPT(IRAP)-1,IPTDIM
     >                    ))
ckr                     write(*,*)"irap,npt(irap),pt",
ckr     >                    irap,npt(irap),pt(npt(irap)+1)
                     CSTRNG = CBASE2
                     WRITE(CTMP,'(I1)'),IORD
                     CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//"_iord="
     >                    //CTMP
                     WRITE(CTMP,'(I1)'),MYISCALE
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
                        DO IPDF=1,MYNPDF
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
     >                    MYISCALE.EQ.3.AND.ISUB.EQ.0) THEN
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
     >                    MYISCALE.EQ.3) THEN
                        CSTRNG = CBASE1
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//
     >                       "_dref/xsect_%"
                        CALL HBOOKB(IHIST+5,
     >                       CSTRNG(1:LENOCC(CSTRNG)),
     >                       NPT(IRAP),PT,0)
                        NHIST = NHIST+1
ckr                        write(*,*)"4. Booked histo #",nhist
                     ENDIF
                     IF (LSTAT.AND.ISUB.EQ.0) THEN
                        IHIST = IORD*1000000 + MYISCALE*100000 +
     >                       ISUB*10000 + IRAP*100
                        CSTRNG = CBASE1(1:LENOCC(CBASE1))
                        WRITE(CTMP,'(I1)'),IORD
                        CSTRNG = CSTRNG(1:LENOCC(CSTRNG))//"_iord="
     >                       //CTMP
                        WRITE(CTMP,'(I1)'),MYISCALE
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
                           WRITE(CTMP,'(I1)'),MYISCALE
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
            ENDDO               ! End od MYISCALE loop, IORD still on
         ENDDO
         WRITE(*,*)"Number of histograms booked:",NHIST



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
      SUBROUTINE PDFFILL(NSCALE,RES0,RES1HI,RES1LO)
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      INTEGER NSCALE
      DOUBLE PRECISION 
     >     RES0  (MXOBSBIN,0:MXSUBPROC,0:MXCTRB),
     >     RES1HI(MXOBSBIN,0:MXSUBPROC,0:MXCTRB),
     >     RES1LO(MXOBSBIN,0:MXSUBPROC,0:MXCTRB)
      INTEGER I,J,NBIN,IORD,ISUB,MYISCALE,IHIST
      REAL VAL0,VALLO,VALHI
      
c --- Mods for v2.0
      REAL PT(MXOBSBIN)
      INTEGER MYICONTR,NORD
      
      MYICONTR = 2
      NORD   = NPOW(MYICONTR) - ILOORD + 1
ckr      write(*,*)"PDFFILL: icontr,npow,iloord,nord",icontr,npow(icontr)
ckr     >     ,iloord,nord


      IF (NSCALE.LT.1 .OR.
     >     NSCALE.GT.NSCALEVAR(MYICONTR,NSCALEDIM(MYICONTR))) THEN
         WRITE(*,*) "\nPDFFILL: ERROR! NSCALE ",NSCALE,
     >        " is out of range, aborted!"
         WRITE(*,*) "PDFFILL: Max. NSCALE: ",
     >        NSCALEVAR(MYICONTR,NSCALEDIM(MYICONTR))
         STOP
      ENDIF
      MYISCALE = NSCALE

c - Fill all histograms for the given scale
      DO IORD=0,NORD            ! Order: tot, LO, NLO-corr, NNLO-corr
         DO ISUB=0,NSUBPROC(MYICONTR) ! Subprocesses: 0 tot + 7 subproc
            NBIN=0
            DO I=1,NRAPIDITY                   
               DO J=1,NPT(I)
                  NBIN = NBIN + 1
                  IF (IORD.GT.0) THEN
                     VAL0  = REAL(RES0(NBIN,ISUB,IORD))
                     VALLO = REAL(RES1LO(NBIN,ISUB,IORD))
                     VALHI = REAL(RES1HI(NBIN,ISUB,IORD))
                  ELSE
ckr                     VAL0  = REAL(RES0(NBIN,ISUB,NORD))
ckr                     VALLO = REAL(RES1LO(NBIN,ISUB,NORD))
ckr                     VALHI = REAL(RES1HI(NBIN,ISUB,NORD))
                     VAL0  = REAL(RES0(NBIN,ISUB,0))
                     VALLO = REAL(RES1LO(NBIN,ISUB,0))
                     VALHI = REAL(RES1HI(NBIN,ISUB,0))
                  ENDIF
ckr Recall: HBOOK understands only single precision
                  IHIST = IORD*1000000+MYISCALE*100000+ISUB*10000+I*100
ckr                  write(*,*)"iord,iscale,isub/2,i",iord,iscale,
ckr     >                 isub,isub2,i
ckr                  write(*,*)"ihist,val0,vallo,valhi",
ckr     >                 ihist,val0,vallo,valhi
ckr                  write(*,*)"pt",REAL(LOBIN(
ckr     >                 IDimPointer(I,IRAPDIM)+J-1,IPTDIM
ckr     >                 )+0.01)
                  CALL HFILL(IHIST,
     >                 REAL(LOBIN(
     >                 IDimPointer(I,IRAPDIM)+J-1,IPTDIM
     >                 )+0.01),0.0,VAL0)
ckr     >                 REAL(LOBIN(RAPINDEX(I)+J-1,1)+0.01),0.0,VAL0)
ckr     >                 REAL(PTBIN(I,J)+0.01),0.0,VAL0)
                  CALL HFILL(IHIST+1,
     >                 REAL(LOBIN(
     >                 IDimPointer(I,IRAPDIM)+J-1,IPTDIM
     >                 )+0.01),0.0,VALLO)
ckr     >                 REAL(LOBIN(RAPINDEX(I)+J-1,1)+0.01),0.0,VALLO)
ckr     >                 REAL(PTBIN(I,J)+0.01),0.0,VALLO)
                  CALL HFILL(IHIST+2,
     >                 REAL(LOBIN(
     >                 IDimPointer(I,IRAPDIM)+J-1,IPTDIM
     >                 )+0.01),0.0,VALHI)
ckr     >                 REAL(LOBIN(RAPINDEX(I)+J-1,1)+0.01),0.0,VALHI)
ckr     >                 REAL(PTBIN(I,J)+0.01),0.0,VALHI)
                  PT(J) = REAL(LOBIN(
     >                 IDimPointer(I,IRAPDIM)+J-1,IPTDIM
     >                 ))
ckr                        PT(J) = REAL(LOBIN(RAPINDEX(I)+J-1,1))
               ENDDO            ! pT-loop
            ENDDO               ! rap-loop
         ENDDO                  ! isub-loop
      ENDDO                     ! iord-loop
      
      RETURN
      END
