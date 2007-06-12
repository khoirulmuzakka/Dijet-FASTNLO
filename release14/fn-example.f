      PROGRAM EXAMPLE
* -------------------------------------------------------------------
* M. Wobisch 04/05/2006
*
* fastNLO - example program to call fastNLO calculations
*           for different scenarios using PDFs from LHAPDF
*
* the user needs to set the correct path to the PDFs in LHAPDF
* 
* -------------------------------------------------------------------
      implicit none

c - Attention!!! - the following arrays must be declared consistent 
c                  with their definitions in the commonblocks 
c                  in the scenario-specific include-files: fn[xxxxx].inc
      double precision 
     +     xst1001(33,3),
     +     xst1002(90,3),
     +     xst1003(51,3),
     +     xst1004(20,3),
     +     xst1005(20,3),
     +     xst1007(18,3),
     +     xst1008(15,3),
     +     xst2001(94,3),
     +     xst2002(20,3),
     +     xst2003(17,3),
     +     xsl0004(103,3),
     +     xsr0001(12,3),
     +     xsh1001(60,3),
     +     xsh1002(30,3),
     +     xsh1003(19,3),
     +     xsh1004(16,3)



c - Initialize LHAPDF    - for CTEQ6.1M      MRST2004NLO:0.1205
c      call InitPDFset('pdfpath/MRST2004nlo.LHpdf')
      call InitPDFset('pdfpath/cteq61.LHgrid')
c      call InitPDFset('pdfpath/a02m_nlo.LHgrid')

c - initialize one member
      call InitPDF(0)      !  0: best fit member
c      call InitPDF(29)    ! 29: EV 15+ for CTEQ6.1M
c      call InitPDF(30)    ! 30: EV 15- for CTEQ6.1M


c - compute the cross sections
      write(*,*) ' ******************************************************'
      write(*,*) '     fastNLO: call user codes for different scenarios'
      write(*,*) ' ******************************************************'
      write(*,*) ' '

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
c
c - the results of the last call can be accessed  
c   in the array:  xst1001(n,iord)
c      n: continuous bin number 
c   iord: order  1 LO, 
c                2 NLO correction (add 1,2 to get the NLO x-section)
c                3 2-loop threshold correction
c         (add 1,2,3 to get full prediction)





c =================================================================
c ======== Tevatron Run I =========================================
c =================================================================

c --- fnt1001 - hep-ph/0102074 - Run I CDF incl. jets (ET)
      call FT1001CC('tablepath/fnt1001midp.tab',1.0d0 ,1.0d0, 1,XST1001)
c      call FT1001CC('tablepath/fnt1001rsep.tab',1.0d0 ,1.0d0, 1,XST1001)

c --- fnt1002 - hep-ex/0011036 - Run I D0 incl. jets (ET, eta)
      call FT1002CC('tablepath/fnt1002midp.tab',1.0d0 ,1.0d0, 1,XST1002)
c      call FT1002CC('tablepath/fnt1002rsep.tab',1.0d0 ,1.0d0, 1,XST1002)

c --- fnt1003 - hep-ex/0012013 - Run I CDF dijets (ET, eta)
      call FT1003CC('tablepath/fnt1003midp.tab',1.0d0 ,1.0d0, 1,XST1003)
c      call FT1003CC('tablepath/fnt1003rsep.tab',1.0d0 ,1.0d0, 1,XST1003)

c --- fnt1005 - hep-ex/0012046 - Run I D0 ratio  incl. jets 630/1800GeV (xT)
c -   -> special interface for ratio (2 tables)
      call FT1005CC('tablepath/fnt1005a-midp.tab',
     +     'tablepath/fnt1005b-midp.tab',1.0d0 ,1.0d0, 1,XST1005)
c      call FT1005CC('tablepath/fnt1005a-rsep.tab',
c     +     'tablepath/fnt1005b-rsep.tab',1.0d0 ,1.0d0, 1,XST1005)

c --- fnt1004 - hep-ex/0012046 - Run I 630GeV D0 incl. jets (ET)
      call FT1004CC('tablepath/fnt1004midp.tab',1.0d0 ,1.0d0, 1,XST1004)
c      call FT1004CC('tablepath/fnt1004rsep.tab',1.0d0 ,1.0d0, 1,XST1004)

c --- fnt1007 - hep-ex/9912022 - Run I CDF dijets (mass)
      call FT1007CC('tablepath/fnt1007midp.tab',1.0d0 ,1.0d0, 1,XST1007)
c      call FT1007CC('tablepath/fnt1007rsep.tab',1.0d0 ,1.0d0, 1,XST1007)

c --- fnt1008 - hep-ex/0012046 - Run I D0 dijets (mass)
      call FT1008CC('tablepath/fnt1008midp.tab',1.0d0 ,1.0d0, 1,XST1008)
c      call FT1008CC('tablepath/fnt1008rsep.tab',1.0d0 ,1.0d0, 1,XST1008)



c =================================================================
c ======== Tevatron Run II ========================================
c =================================================================

c --- fnt2001 - hep-ex/0409040 - Run II D0 dijets DeltaPhi
c -   -> special interface for ratio (use 2 tables)
      call FT2001CC('tablepath/fnt2001a.tab',
     +     'tablepath/fnt2001b.tab',1.0d0 ,1.0d0, 1,XST2001)

c --- fnt2002 - hep-ex/0512020 - CDF incl. jets cone (pT)
      call FT2002CC('tablepath/fnt2002midp.tab',1.0d0 ,1.0d0, 1,XST2002)
      call FT2002CC('tablepath/fnt2002rsep.tab',1.0d0 ,1.0d0, 1,XST2002)

c --- fnt2003 - hep-ex/0512062 - CDF incl. jets kT (pT)
      call FT2003CC('tablepath/fnt2003.tab',1.0d0 ,1.0d0, 1,XST2003)


c =================================================================
c ======== RHIC ===================================================
c =================================================================
c --- fnr0001 - STAR preliminary - incl. jets (pT)
      call FR0001CC('tablepath/fnr0001midp.tab',1.0d0 ,1.0d0, 1,XSR0001)


c =================================================================
c ======== LHC ====================================================
c =================================================================
c - fnl0004 - test scenario incl jets (pT,y) optimized for ATLAS
      call FL0004CC('tablepath/fnl0004.tab',1.0d0 ,1.0d0, 1,XSL0004)


c =================================================================
c ======== HERA ===================================================
c =================================================================
c - fnh1001 - hep-ex/0010054 - H1 incl. jets (ET, Q2)
      call FH1001CC('tablepath/fnh1001.tab',1.0d0 ,1.0d0, 1,XSH1001)

c - fnh1002 - hep-ex/0208037 - ZEUS incl. jets (ET, Q2) 
      call FH1002CC('tablepath/fnh1002.tab',1.0d0 ,1.0d0, 1,XSH1002)

c - fnh1003 - hep-ex/0206029, HERA, H1 incl. jets  (ET) @low Q2
      call FH1003CC('tablepath/fnh1003.tab',1.0d0 ,1.0d0, 1,XSH1003)

c - fnh1004 - hep-ex/0010054, HERA, H1 dijets (ET, Q2) 
      call FH1004CC('tablepath/fnh1004.tab',1.0d0 ,1.0d0, 1,XSH1004)



      END

