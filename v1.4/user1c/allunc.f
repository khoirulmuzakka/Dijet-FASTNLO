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
      CHARACTER*255 BORNNAME,NLONAME
      CHARACTER*255 HISTFILE
      CHARACTER*255 SCENARIO,FILENAME,TABPATH,TABNAME,REFNAME
      CHARACTER*255 PDFSET,PDFPATH,LHAPDF,CHTMP
      CHARACTER*8 CH8TMP
      CHARACTER*4 CH4TMP
      INTEGER BORNN,NLON,LENOCC
      INTEGER I,J,L1,L2,L3,L4,NPDF,IOPDF,IOAS
      INTEGER ISTAT,ISCALE,IORD,IBIN,NBIN,ISUB,IRAP,IPT,IHIST
      LOGICAL LALG,LSTAT,LPDF,LSER
      DOUBLE PRECISION MUR,MUF,DIFF,QLAM4,QLAM5
      DOUBLE PRECISION
     >     RES0(NBINTOTMAX,NMAXSUBPROC+1,0:3),
     >     RES1HI(NBINTOTMAX,NMAXSUBPROC+1,0:3),
     >     RES1LO(NBINTOTMAX,NMAXSUBPROC+1,0:3),
     >     RESLO,RESHI,DREF
c - Attention!!! This must be declared consistent with the
c                definition in the commonblock!!!!!
      DOUBLE PRECISION XSECT0(NBINTOTMAX,3),XSECT1(NBINTOTMAX,3)
      REAL PT(NPTMAX)

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
            WRITE(*,*)'  alpha_s calc., def. from PDF set'
            WRITE(*,*)'  alpha_s(M_Z), def. from PDF set'
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
     >        "using alpha_s according to PDF set (mode PDF) or"
         WRITE(*,*)
     >        "        PDG book (mode MW and KR)"
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

*---Too many arguments
      IF (IARGC().GT.11) THEN
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
      LSER  = NPDF.LT.10.AND..NOT.LSTAT.AND..NOT.LALG
      LPDF  = .NOT.LSER
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
      CALL FX9999CC(FILENAME,1D0,1D0,0,XSECT1)
      CALL PDFHIST(1,HISTFILE,LPDF,LSTAT,LALG,LSER,NPDF)



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
      DO I=1,NSCALEVAR
         CALL INITPDF(0)
         MUR = MURSCALE(I)
         MUF = MUFSCALE(I)
         WRITE(*,*)"ALLUNC: Now scale no.",i,"; mur, muf = ",mur,muf
         CALL FX9999CC(FILENAME,MUR,MUF,0,XSECT1)

c - Save the result array from the first call (= central result)
c   and reset result arrays   
         WRITE(*,*)"ALLUNC: The observable has",NBINTOT," bins -",
     >        NSUBPROC," subprocesses"
         DO L1=1,NBINTOT
            DO L2=1,(NSUBPROC+1)
               DO L3=1,NORD
                  RES0(L1,L2,L3) = 0D0 
                  DO L4=1,L3
                     RES0(L1,L2,L3) = RES0(L1,L2,L3)+RESULT(L1,L2,L4)
                  ENDDO
                  RES1LO(L1,L2,L3) = 0D0
                  RES1HI(L1,L2,L3) = 0D0
               ENDDO
            ENDDO
         ENDDO

ckr Do loop runs once even if NPDF=0! => Avoid with IF statement
ckr         IF (NPDF.GT.1) THEN
         IF (LPDF) THEN
            DO J=1,NPDF
               CALL INITPDF(J)
               CALL FX9999CC(FILENAME,MUR,MUF,0,XSECT0)
c - For all bins/subproc/orders: Add negative/positive variations
               DO L1=1,NBINTOT
                  DO L2=1,(NSUBPROC+1)
                     DO L3=1,NORD
                        DIFF = - RES0(L1,L2,L3)
                        DO L4=1,L3
                           DIFF = DIFF + RESULT(L1,L2,L4) 
                        ENDDO
                        IF (DIFF.GT.0D0) THEN
                           RES1HI(L1,L2,L3) =
     >                          RES1HI(L1,L2,L3)+DIFF*DIFF
                        ELSE
                           RES1LO(L1,L2,L3) =
     >                          RES1LO(L1,L2,L3)+DIFF*DIFF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO               ! Loop over bins
         ENDIF                  ! Not done for npdf <= 1

c - Take square-root of sum of squares
         DO L1=1,NBINTOT
ckr            IF (NPDF.GT.1) THEN
            IF (LPDF) THEN
               DO L2=1,(NSUBPROC+1)
                  DO L3=1,NORD
                     RES1HI(L1,L2,L3) =  SQRT(RES1HI(L1,L2,L3))
                     RES1LO(L1,L2,L3) = -SQRT(RES1LO(L1,L2,L3))
                  ENDDO
               ENDDO
               RESLO = RES1LO(L1,NSUBPROC+1,NORD)/
     >              RES0(L1,NSUBPROC+1,NORD)
               RESHI = RES1HI(L1,NSUBPROC+1,NORD)/
     >              RES0(L1,NSUBPROC+1,NORD)
ckr 30.01.2008: Change output format for better comp. with C++ version
            ELSE
               RESLO = 0D0
               RESHI = 0D0
            ENDIF
            WRITE(*,900) L1,RES0(L1,NSUBPROC+1,NORD),RESLO,RESHI
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
ckr         BORNNAME = TABPATH(1:LENOCC(TABPATH))//"/"//
ckr     >        SCENARIO(1:LENOCC(SCENARIO))//"-hhc-born-2jet_"
ckr         NLONAME  = TABPATH(1:LENOCC(TABPATH))//"/"//
ckr     >        SCENARIO(1:LENOCC(SCENARIO))//"-hhc-nlo-2jet_"
         BORNNAME = TABPATH(1:LENOCC(TABPATH))//"/stat/"//
     >        SCENARIO(1:LENOCC(SCENARIO))//"-hhc-born-2jet_"
         NLONAME  = TABPATH(1:LENOCC(TABPATH))//"/stat/"//
     >        SCENARIO(1:LENOCC(SCENARIO))//"-hhc-nlo-2jet_"
         CALL STATCODE(BORNN,BORNNAME,NLON,NLONAME)
      ENDIF
      
      
      
c - PDF series of variations (e.g. alpha_s(M_Z))
c - Use primary table
c - Check that FILENAME is still the primary table here ...!!!
c - Fill only central scale
      IF (LSER) THEN
         ISCALE = 3
         MUR = MURSCALE(ISCALE)
         MUF = MUFSCALE(ISCALE)
         WRITE(*,*)"ALLUNC: For PDF series fill only scale no.",ISCALE,
     >        "; mur, muf = ",mur,muf
         DO J=0,NPDF
            CALL INITPDF(J)
            CALL FX9999CC(FILENAME,MUR,MUF,0,XSECT0)
            IBIN = 0
            DO IRAP=1,NRAPIDITY
               DO IPT=1,NPT(IRAP)
                  IBIN = IBIN+1
                  PT(IBIN) = REAL(PTBIN(IRAP,IPT))
                  DO ISUB=0,NSUBPROC
                     RES0(IBIN,ISUB+1,0) = 0.D0
                     DO IORD=1,NORD
                        RES0(IBIN,ISUB+1,IORD) = XSECT0(IBIN,IORD)
                        RES0(IBIN,ISUB+1,0) = RES0(IBIN,ISUB+1,0) +
     >                       RES0(IBIN,ISUB+1,IORD)
                     ENDDO
c - Fill histograms
                     DO IORD=0,NORD
                        IHIST = IORD*1000000 + ISCALE*100000 +
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
         ISCALE = 3
         MUR = MURSCALE(ISCALE)
         MUF = MUFSCALE(ISCALE)
         CALL FX9999CC(FILENAME,MUR,MUF,1,XSECT)
c - Attention: From now on ref. table loaded ==> rap. bins doubled
c - Get the normal result, scale 3, from first half of doubled rap bins
         IBIN = 0
         DO IRAP=1,INT(NRAPIDITY/2)
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               DO ISUB=0,NSUBPROC
                  RES0(IBIN,ISUB+1,0) = 0.D0
                  DO IORD=1,NORD
                     RES0(IBIN,ISUB+1,IORD) = XSECT(IBIN,IORD)
                     RES0(IBIN,ISUB+1,0) = RES0(IBIN,ISUB+1,0) +
     >                    RES0(IBIN,ISUB+1,IORD)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         NBIN = IBIN

c - Reference scale no. always 1
         ISCALE = 1
         MUR = MURSCALE(ISCALE)
         MUF = MUFSCALE(ISCALE)
         CALL FX9999CC(FILENAME,MUR,MUF,1,XSECT)
c - Get the reference result, scale 1, nrap/2 ++ bins
         DO IRAP=INT(NRAPIDITY/2)+1,NRAPIDITY
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               DO ISUB=0,NSUBPROC
                  RES0(IBIN,isub+1,0) = 0.D0
                  DO IORD=1,NORD
                     RES0(IBIN,isub+1,IORD) = XSECT(IBIN,IORD)
                     RES0(IBIN,isub+1,0) = RES0(IBIN,isub+1,0) +
     >                    RES0(IBIN,isub+1,IORD)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

c - Compare results and fill histos
         ISCALE = 3
         ISUB   = 0
         IBIN   = 0
         DO IRAP=1,INT(NRAPIDITY/2)
            DO IPT=1,NPT(IRAP)
               IBIN = IBIN+1
               PT(IBIN) = REAL(PTBIN(IRAP,IPT))
               DO ISUB=0,NSUBPROC
                  DO IORD=0,NORD
                     DREF = RES0(IBIN,isub+1,IORD)/
     >                    RES0(IBIN+NBIN,isub+1,IORD) - 1.D0
                     IHIST = IORD*1000000 + ISCALE*100000 +
     >                    ISUB*10000 + IRAP*100
                     CALL HFILL(IHIST+5,PT(IBIN),0.,REAL(100D0*DREF))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF



c - Close hbook file
      CALL PDFHIST(2,HISTFILE,LPDF,LSTAT,LALG,LSER,NPDF)
      END

c
c ======================= Book the histograms ========================
c
      SUBROUTINE PDFHIST(N,HISTFILE,LPDF,LSTAT,LALG,LSER,NPDF)
      IMPLICIT NONE
      CHARACTER*(*) HISTFILE
      CHARACTER*255 CSTRNG,CBASE1,CBASE2,CTMP
      INTEGER N,LENOCC,IPDF,NPDF
      LOGICAL LPDF,LSTAT,LALG,LSER

      INTEGER J,ISTAT2,ICYCLE
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
            DO ISCALE=1,NSCALEVAR ! Scale variations
               DO ISUB=0,NSUBPROC ! Subprocesses: 0 tot + 7 subproc
                  DO IRAP=1, NRAPIDITY
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
ckr                        write(*,*)"5. Booked histo #",nhist
                        DO IPT=1,NPT(IRAP)
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
ckr                           write(*,*)"6. Booked histo #",nhist
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO               ! End od ISCALE loop, IORD still on
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
     >     RES0(NBINTOTMAX,NMAXSUBPROC+1,0:3),
     >     RES1HI(NBINTOTMAX,NMAXSUBPROC+1,0:3),
     >     RES1LO(NBINTOTMAX,NMAXSUBPROC+1,0:3)
      INTEGER I,J,NBIN,IORD,ISUB,ISUB2,ISCALE,IHIST
      REAL VAL0,VALLO,VALHI
      
      IF (NSCALE.LT.1 .OR. NSCALE.GT.NSCALEVAR) THEN
         WRITE(*,*) "\nPDFFILL: ERROR! NSCALE ",NSCALE,
     >        " is out of range, aborted!"
         WRITE(*,*) "PDFFILL: Max. NSCALE: ",NSCALEVAR
         STOP
      ENDIF
      ISCALE = NSCALE

c - Fill all histograms for the given scale
      DO IORD=0,NORD            ! Order: tot, LO, NLO-corr, NNLO-corr
         DO ISUB2=1,(NSUBPROC+1) ! Subprocesses: 0 tot + 7 subproc
            ISUB=ISUB2
            IF (ISUB.EQ.8) ISUB=0
            NBIN=0
            DO I=1,NRAPIDITY                   
               DO J=1,NPT(I)
                  NBIN = NBIN + 1
                  IF (IORD.GT.0) THEN
                     VAL0  = REAL(RES0(NBIN,ISUB2,IORD))
                     VALLO = REAL(RES1LO(NBIN,ISUB2,IORD))
                     VALHI = REAL(RES1HI(NBIN,ISUB2,IORD))
                  ELSE
                     VAL0  = REAL(RES0(NBIN,ISUB2,NORD))
                     VALLO = REAL(RES1LO(NBIN,ISUB2,NORD))
                     VALHI = REAL(RES1HI(NBIN,ISUB2,NORD))
                  ENDIF
ckr Recall: HBOOK understands only single precision
                  IHIST = IORD*1000000+ISCALE*100000+ISUB*10000+I*100
ckr                  write(*,*)"iord,iscale,isub/2,i",iord,iscale,
ckr     >                 isub,isub2,i
ckr                  write(*,*)"ihist,val0,vallo,valhi",
ckr     >                 ihist,val0,vallo,valhi
                  CALL HFILL(IHIST,  REAL(PTBIN(I,J)+0.01),0.0,VAL0)
                  CALL HFILL(IHIST+1,REAL(PTBIN(I,J)+0.01),0.0,VALLO)
                  CALL HFILL(IHIST+2,REAL(PTBIN(I,J)+0.01),0.0,VALHI)
               ENDDO            ! pT-loop
            ENDDO               ! rap-loop
         ENDDO                  ! isub-loop
      ENDDO                     ! iord-loop
      
      RETURN
      END
