      PROGRAM EXAMPLE
* -------------------------------------------------------------------
* M. Wobisch 06/29/2005
*
* fastNLO - example program to compute Run I cross section
*           using PDFs from LHAPDF
*
* 10.06.2011 kr: Replace CERNLIB function LENOCC by f95 Standard LEN_TRIM
* 04.03.2006 kr: Get LHAPDF path from environment variable
* -------------------------------------------------------------------
      implicit none
      CHARACTER*255 FILENAME
      CHARACTER*255 HISTOFILE
      CHARACTER*255 LHAPDF
      integer i

ckr Added to make compatible with modified fn-interface.f
      CHARACTER*255 ASMODE
      DOUBLE PRECISION ASMZVAL
      INTEGER IASLOOP
      COMMON/STEER/ASMZVAL,IASLOOP,ASMODE

c - Attention!!! - this mus be declared consistent with its 
c                  definition in the commonblock!!!!!
      double precision xst1001(700,3) 

ckr Default: alpha_s from PDF
      ASMODE = "PDF"

c --- parse command line
      IF ( IARGC().LT.1)  FILENAME = 'table.txt.gz'
      IF ( IARGC().LT.2)  HISTOFILE= 'fastnlo.hbk'
      IF ( IARGC().GT.0)  CALL GETARG(1,FILENAME)
      IF ( IARGC().GT.1)  CALL GETARG(2,HISTOFILE)
      IF ( IARGC().GT.2)  THEN
         write(*,*) 'fastNLO: too many arguments given. Stopping'
         STOP
      ENDIF

c - Initialize LHAPDF    - for CTEQ6.1M   
      CALL GETENV('LHAPDF',LHAPDF)
      IF (LEN_TRIM(LHAPDF).EQ.0) THEN
         LHAPDF = '/disk2/work/wobisch/lhapdf-4.1/PDFsets/cteq61.LHgrid'
      ENDIF
      write(*,*)"Looking for LHAPDF in directory "//
     &     LHAPDF(1:LEN_TRIM(LHAPDF))//"!"
      call InitPDFset(LHAPDF(1:LEN_TRIM(LHAPDF))/
     >     /'/../share/lhapdf/PDFsets/cteq61.LHgrid')

c - initialize one member, 0=best fit member
      call InitPDF(0)

c - .... or initialize CTEQ code - needs table file:  ctq61.00.tbl
c      Call SetCtq6(200)         ! for CTEQ6.1M


c - compute the cross sections
      write(*,*) '     fastNLO: compute the cross section'

c- new call: a single call for each scale
c         1st argument:  name of table
c         2nd argument:  xmur  prefactor for nominal ren-scale
c                              any choice is possible, but please note 
c                              that NNLO-NLL works only for xmur=xmuf
c         3rd argument:  xmuf  prefactor for nominal fact-scale
c                              only a few choices are possible
c                              (see output or table documentation)
c         4th argument:  0: no ascii output       1: print results
c         5th argument:  array to return results

      call FX9999CC(FILENAME, 0.25d0, 0.25d0,1 , XST1001)
      call FX9999CC(FILENAME, 0.5d0 , 0.5d0, 1 , XST1001)
      call FX9999CC(FILENAME, 1.0d0 , 1.0d0, 1 , XST1001)
      call FX9999CC(FILENAME, 2.0d0 , 2.0d0, 1 , XST1001)

c - the results of the last call can be accessed in the array:  xst1001(n,iord)
c      n: continuous bin number for all D0/CDF bins (see documentation)
c   iord: order  1 LO, 
c                2 NLO correction (add 1,2 to get the NLO x-section)
c                3 NNLO-NLL correction (add 1,2,3 to get full prediction)




c- check speed
c      do i=1,100
c         call FX9999CC(FILENAME, 0.5d0 , 0.5d0, 0 , XST1001)
c      enddo



c - book, fill and store histograms - works now for v1.4! (Feb 8, 2006 MW)
c     Histogram-code books & fills histos only for all default scales, 
c     independent of previous FX9999CC calls.
c     -> But at least one previous call of  FX9999CC is required
c        because the commonblock needs to be filled before the call
c
      call FNHIST(FILENAME,HISTOFILE,1)  ! to plot x-sect histos
c      call FNHIST(FILENAME,HISTOFILE,2)  ! to plot x-sect histos + x-histos

      END
