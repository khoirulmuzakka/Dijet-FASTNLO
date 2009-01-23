      PROGRAM PDFUNC
* ---------------------------------------------------------------------
* M. Wobisch  02/09/2006
* K. Rabbertz 01/06/2008 Restructured and cleaned up version
*
* PDFUNC - example program to compute PDF uncertainties
*          using a fastNLO table and PDFs from LHAPDF
*
* ---------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      CHARACTER*255 FILENAME,HISTFILE,PDFSET,PDFPATH,LHAPDF,ASMODE
      INTEGER I,J,L1,L2,L3,L4,NPDF,IOPDF,IOAS
      INTEGER LENOCC
      DOUBLE PRECISION MUR,MUF,DIFF,QLAM4,QLAM5
      DOUBLE PRECISION
     >     RES0(NBINTOTMAX,NMAXSUBPROC+1,3),
     >     RES1HI(NBINTOTMAX,NMAXSUBPROC+1,3),
     >     RES1LO(NBINTOTMAX,NMAXSUBPROC+1,3),
     >     RESLO,RESHI
c - Attention!!! This mus be declared consistent with the
c                definition in the commonblock!!!!!
      DOUBLE PRECISION XSECT0(NBINTOTMAX,3),XSECT1(NBINTOTMAX,3)
      COMMON/STEER/ASMODE

c --- Parse command line
ckr 30.01.2008: Some more checks on input arguments
      WRITE(*,*)"\n ##############################################"
      WRITE(*,*)"# PDFUNC"
      WRITE(*,*)"##############################################"
      WRITE(*,*)"# Example program to compute PDF uncertainties"
      WRITE(*,*)"# using a fastNLO table and PDFs from LHAPDF"
      WRITE(*,*)"##############################################"
      WRITE(*,*)"#"
      IF (IARGC().LT.1) THEN
         FILENAME = "table.txt"
         WRITE(*,*)
     >        "PDFUNC: WARNING! No input table given, "//
     >        "taking table.txt instead!"
         WRITE(*,*)"      For an explanation of command line "//
     >        "arguments type:"
         WRITE(*,*)"      ./pdfunc -h"
      ELSE
         CALL GETARG(1,FILENAME)
         IF (FILENAME(1:LENOCC(FILENAME)).EQ."-h") THEN
            WRITE(*,*)" "
            WRITE(*,*)"Usage: ./pdfunc [arguments]"
            WRITE(*,*)"  NLO input table, def. = table.txt"
            WRITE(*,*)"  HBOOK output file, def. = fastnlo.hbk"
            WRITE(*,*)"  PDF set, def. = cteq65.LHgrid"
            WRITE(*,*)"  PDF path, def. = $(LHAPDF)/"//
     >           "../share/lhapdf/PDFsets"
            WRITE(*,*)"  alpha_s calc., def. from PDF set"
            WRITE(*,*)" "
            STOP
         ENDIF
         WRITE(*,*)"PDFUNC: Using input table: ",
     >        FILENAME(1:LENOCC(FILENAME))
      ENDIF
      IF (IARGC().LT.2) THEN
         HISTFILE = "fastnlo.hbk"
         WRITE(*,*)
     >        "PDFUNC: WARNING! No output filename given, "//
     >        "taking fastnlo.hbk instead!"
      ELSE
         CALL GETARG(2,HISTFILE)
         WRITE(*,*)"PDFUNC: Creating output file: ",
     >        HISTFILE(1:LENOCC(HISTFILE))
      ENDIF
      IF (IARGC().LT.3) THEN
         PDFSET = "cteq65.LHgrid"
         WRITE(*,*)
     >        "PDFUNC: WARNING! No PDF set given, "//
     >        "taking cteq65.LHgrid instead!"
      ELSE
         CALL GETARG(3,PDFSET)
         WRITE(*,*)"PDFUNC: Using PDF set: ",
     >        PDFSET(1:LENOCC(PDFSET))
      ENDIF
      IF (IARGC().LT.4) THEN
         PDFPATH = "/../share/lhapdf/PDFsets"
         WRITE(*,*)
     >        "PDFUNC: No PDF path given, "//
     >        "assuming: $(LHAPDF)"//PDFPATH
c - Initialize path to LHAPDF libs
         CALL GETENV("LHAPDF",LHAPDF)
         IF (LENOCC(LHAPDF).EQ.0) THEN
            WRITE(*,*)"\nPDFUNC: ERROR! $LHAPDF not set, aborting!"
            STOP
         ENDIF
         PDFPATH = LHAPDF(1:LENOCC(LHAPDF))//
     >        PDFPATH(1:LENOCC(PDFPATH))
      ELSE
         CALL GETARG(4,PDFPATH)
      ENDIF
      WRITE(*,*)"PDFUNC: Looking for LHAPDF PDF sets in path: ",
     >     PDFPATH(1:LENOCC(PDFPATH))
      PDFSET = PDFPATH(1:LENOCC(PDFPATH))//"/"//PDFSET
      WRITE(*,*)"PDFUNC: Taking PDF set "
     >     //PDFSET(1:LENOCC(PDFSET))
      IF (IARGC().LT.5) THEN
         ASMODE = "PDF"
         WRITE(*,*)
     >        "PDFUNC: No alpha_s mode given, "//
     >        "using alpha_s according to PDF set"
      ELSE
         CALL GETARG(5,ASMODE)
         WRITE(*,*)"PDFUNC: Using alpha_s mode:",ASMODE
      ENDIF
      IF (IARGC().GT.5) THEN
         WRITE(*,*)"\nPDFUNC: ERROR! Too many arguments, aborting!"
         STOP
      ENDIF

c - Initialization      
      CALL INITPDFSET(PDFSET(1:LENOCC(PDFSET)))
      CALL NUMBERPDF(NPDF)
      CALL GETORDERPDF(IOPDF)
      CALL GETORDERAS(IOAS)
      CALL GETLAM4(0,QLAM4)
      CALL GETLAM5(0,QLAM5)
      WRITE(*,*) "PDFUNC: The PDF set has",NPDF+1," members"
      WRITE(*,*) "PDFUNC: The PDF is of order",IOPDF+1
      WRITE(*,*) "PDFUNC: alpha_s was used in",IOAS+1,
     >     "-loop order in the PDF"
      WRITE(*,*) "PDFUNC: The lambda_4 value for member 0 is",QLAM4
      WRITE(*,*) "PDFUNC: The lambda_5 value for member 0 is",QLAM5

c - One initial call - to fill commonblock -> for histo-booking
      CALL INITPDF(0)
      CALL FX9999CC(FILENAME,1D0,1D0,0,XSECT1)
      CALL PDFHIST(1,HISTFILE)

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
      DO I=1,NSCALEVAR
         CALL INITPDF(0)
         MUR = MURSCALE(I)
         MUF = MUFSCALE(I)
         WRITE(*,*)"PDFUNC: Now scale no.",i,"; mur, muf = ",mur,muf
         CALL FX9999CC(FILENAME,MUR,MUF,0,XSECT1)

c - Save the result array from the first call (= central result)
c   and reset result arrays   
         WRITE(*,*)"PDFUNC: The observable has",NBINTOT," bins -",
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
         IF (NPDF.GT.1) THEN
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
            IF (NPDF.GT.1) THEN
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

c - Close hbook file
      CALL PDFHIST(2,HISTFILE)
      END

c
c ======================= Book the histograms ========================
c
      SUBROUTINE PDFHIST(N,HISTFILE)
      IMPLICIT NONE
      CHARACTER*(*) HISTFILE
      INTEGER N

      INTEGER J,ISTAT2,ICYCLE
      INTEGER IORD,ISUB,ISCALE,IRAP,IHIST,NHIST
      INCLUDE "fnx9999.inc"
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
         DO IORD=0,NORD         ! Order: tot, LO, NLO-corr, NNLO-corr
            DO ISCALE=1,NSCALEVAR ! Scale variations
               DO ISUB=0,NSUBPROC ! Subprocesses: 0 tot + 7 subproc
                  DO IRAP=1, NRAPIDITY
                     IHIST = IORD*1000000 + ISCALE*100000 +
     >                    ISUB*10000 + IRAP*100
                     DO J=1,(NPT(IRAP)+1)
                        PT(J) = REAL(PTBIN(IRAP,J))
                     ENDDO
                     CALL HBOOKB(IHIST,
     >                    "P?T! (GEV)",
     >                    NPT(IRAP),PT,0)
                     CALL HBOOKB(IHIST+1,
     >                    "P?T! (GEV)",
     >                    NPT(IRAP),PT,0)
                     CALL HBOOKB(IHIST+2,
     >                    "P?T! (GEV)",
     >                    NPT(IRAP),PT,0)
                     NHIST = NHIST+3
                  ENDDO
               ENDDO
            ENDDO
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
     >     RES0(NBINTOTMAX,NMAXSUBPROC+1,3),
     >     RES1HI(NBINTOTMAX,NMAXSUBPROC+1,3),
     >     RES1LO(NBINTOTMAX,NMAXSUBPROC+1,3)
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
                  IHIST = IORD*1000000+ISCALE*100000+ISUB*10000+I*100
                  CALL HFILL(IHIST,  REAL(PTBIN(I,J)+0.01),0.0,VAL0)
                  CALL HFILL(IHIST+1,REAL(PTBIN(I,J)+0.01),0.0,VALLO)
                  CALL HFILL(IHIST+2,REAL(PTBIN(I,J)+0.01),0.0,VALHI)
               ENDDO            ! pT-loop
            ENDDO               ! rap-loop
         ENDDO                  ! isub-loop
      ENDDO                     ! iord-loop
      
      RETURN
      END
