      PROGRAM EXAMPLE_REF
* -------------------------------------------------------------------
* M. Wobisch 01/18/2006
*
* fastNLO - example program to evaluate reference jobs
*           using PDFs from CTEQ6 code (to be consistent 
*           reference table)
* -------------------------------------------------------------------
      implicit none
      CHARACTER*255 FILENAME
      CHARACTER*255 HISTOFILE

      double precision h1incl(4,4) , h1dijet(4,5) ! for FNH0001

c - Attention!!! - this mus be declared consistent with its 
c                  definition in the commonblock!!!!!
      double precision xst1001(600,3,4) 


c --- parse command line
      IF ( IARGC().LT.1)  FILENAME = 'table.txt'
      IF ( IARGC().LT.2)  HISTOFILE= 'fastnlo.hbk'
      IF ( IARGC().GT.0)  CALL GETARG(1,FILENAME)
      IF ( IARGC().GT.1)  CALL GETARG(2,HISTOFILE)
      IF ( IARGC().GT.2)  THEN
         write(*,*) 'fastNLO: too many arguments given. Stopping'
         RETURN
      ENDIF

      
c - Initialize CTEQ code - needs table file:  ctq61.00.tbl
      Call SetCtq6(200)         ! for CTEQ6.1M

c - Initialize LHAPDF    - for CTEQ6.1M   
c      call InitPDFset('/disk2/work/wobisch/lhapdf-4.1/PDFsets/cteq61.LHgrid')
c
c      call InitPDFset('/h1/h1gen/lhapdf/LHAPDFv4/PDFsets/cteq61.LHgrid')
c - initialize one member, 0=best fit member
c      call InitPDF(0)

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

c      call FX9999CC(FILENAME, 0.25d0, 0.25d0,1 , XST1001)
c      call FX9999CC(FILENAME, 0.5d0 , 0.5d0, 1 , XST1001)
      call FX9999CC(FILENAME, 1.0d0 , 1.0d0, 1 , XST1001)
c      call FX9999CC(FILENAME, 2.0d0 , 2.0d0, 1 , XST1001)


c - book, fill and store histograms - works now for v1.4! (Feb 8, 2006 MW)
c     Histogram-code books & fills histos only for all default scales, 
c     independent of previous FX9999CC calls.
c     -> But at least one previous call of  FX9999CC is required
c        because the commonblock needs to be filled before the call
c
      call FNHIST(FILENAME,HISTOFILE,1)  ! to plot x-sect histos
c      call FNHIST(FILENAME,HISTOFILE,2)  ! to plot x-sect histos + x-histos


      END
