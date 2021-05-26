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
     +     xst2004(76,3),
     +     xst2007(76,3),
     +     xst2008(21,3),
     +     xst2009(110,3),
     +     xst2010(120,3),
     +     xst2011(71,3),
     +     xsh1001(60,3),
     +     xsh1002(30,3),
     +     xsh1003(19,3),
     +     xsh1004(16,3),
     +     xsh2001(30,3),
     +     xsh2003(24,3),
     +     xsr0001(12,3),
     +     xsl0004(103,3),
     +     xsl1013(176,3),
     +     xsl1014(125,3)



c - Initialize LHAPDF
      call InitPDFset('pdfpath/CT10.LHgrid')

c - initialize one member
      call InitPDF(0)      !  0: best fit member
c      call InitPDF(29)    ! 29: member 29


c - compute the cross sections

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

      WRITE(*,*)"\n ########################################"//
     >     "################################"
      WRITE(*,*)"# fastnlo: Call user codes for different scenarios"

c =================================================================
c ======== Tevatron Run I =========================================
c =================================================================
Comment:       WRITE(*,*)"########################################"//
Comment:      >     "################################"
Comment:       WRITE(*,*)"# Tevatron Run I"
Comment:       WRITE(*,*)"########################################"//
Comment:      >     "################################"
Comment: 
Comment: c --- fnt1001 - hep-ph/0102074 - Run I CDF incl. jets (ET)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt1001"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT1001CC('tablepath/fnt1001midp.tab',1.0d0 ,1.0d0, 1,XST1001)
Comment: c      call FT1001CC('tablepath/fnt1001rsep.tab',1.0d0 ,1.0d0, 1,XST1001)
Comment: 
Comment: c --- fnt1002 - hep-ex/0011036 - Run I D0 incl. jets (ET, eta)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt1002"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT1002CC('tablepath/fnt1002midp.tab',1.0d0 ,1.0d0, 1,XST1002)
Comment: c      call FT1002CC('tablepath/fnt1002rsep.tab',1.0d0 ,1.0d0, 1,XST1002)
Comment: 
Comment: c --- fnt1003 - hep-ex/0012013 - Run I CDF dijets (ET, eta)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt1003"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT1003CC('tablepath/fnt1003midp.tab',1.0d0 ,1.0d0, 1,XST1003)
Comment: c      call FT1003CC('tablepath/fnt1003rsep.tab',1.0d0 ,1.0d0, 1,XST1003)
Comment: 
Comment: c --- fnt1004 - hep-ex/0012046 - Run I 630GeV D0 incl. jets (ET)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt1004"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT1004CC('tablepath/fnt1004midp.tab',1.0d0 ,1.0d0, 1,XST1004)
Comment: c      call FT1004CC('tablepath/fnt1004rsep.tab',1.0d0 ,1.0d0, 1,XST1004)
Comment: 
Comment: c --- fnt1005 - hep-ex/0012046 - Run I D0 ratio  incl. jets 630/1800GeV (xT)
Comment: c -   -> special interface for ratio (2 tables)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt1005"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT1005CC('tablepath/fnt1005a-midp.tab',
Comment:      +     'tablepath/fnt1005b-midp.tab',1.0d0 ,1.0d0, 1,XST1005)
Comment: c      call FT1005CC('tablepath/fnt1005a-rsep.tab',
Comment: c     +     'tablepath/fnt1005b-rsep.tab',1.0d0 ,1.0d0, 1,XST1005)
Comment: 
Comment: c --- fnt1006 todo
Comment: 
Comment: c --- fnt1007 - hep-ex/9912022 - Run I CDF dijets (mass)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt1007"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT1007CC('tablepath/fnt1007midp.tab',1.0d0 ,1.0d0, 1,XST1007)
Comment: c      call FT1007CC('tablepath/fnt1007rsep.tab',1.0d0 ,1.0d0, 1,XST1007)
Comment: 
Comment: c --- fnt1008 - hep-ex/0012046 - Run I D0 dijets (mass)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt1008"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT1008CC('tablepath/fnt1008midp.tab',1.0d0 ,1.0d0, 1,XST1008)
Comment: c      call FT1008CC('tablepath/fnt1008rsep.tab',1.0d0 ,1.0d0, 1,XST1008)
Comment: 
Comment: 
Comment: 
Comment: c =================================================================
Comment: c ======== Tevatron Run II ========================================
Comment: c =================================================================
Comment:       WRITE(*,*)"########################################"//
Comment:      >     "################################"
Comment:       WRITE(*,*)"# Tevatron Run II"
Comment:       WRITE(*,*)"########################################"//
Comment:      >     "################################"
Comment: 
Comment: c --- fnt2001 - hep-ex/0409040 - Run II D0 dijets DeltaPhi
Comment: c -   -> special interface for ratio (use 2 tables)
Comment: ckr Files not on HepForge
Comment: ckr      call FT2001CC('tablepath/fnt2001a.tab',
Comment: ckr     +     'tablepath/fnt2001b.tab',1.0d0 ,1.0d0, 1,XST2001)
Comment: 
Comment: c --- fnt2002 - hep-ex/0512020 - CDF incl. jets cone (pT)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt2002"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT2002CC('tablepath/fnt2002midp.tab',1.0d0 ,1.0d0, 1,XST2002)
Comment:       call FT2002CC('tablepath/fnt2002rsep.tab',1.0d0 ,1.0d0, 1,XST2002)
Comment: 
Comment: c --- fnt2003 - hep-ex/0512062 - CDF incl. jets kT (pT)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt2003"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT2003CC('tablepath/fnt2003.tab',1.0d0 ,1.0d0, 1,XST2003)
Comment: 
Comment: c --- fnt2004 - hep-ex/0701051 - CDF incl. jets, kT algo D=0.7 (pT, y)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt2004"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT2004CC('tablepath/fnt2004.tab',1.0d0 ,1.0d0, 1,XST2004)
Comment: 
Comment: c --- fnt2007 - arXiv:0807.2204 [hep-ex] - CDF incl. jets, cone algo (pT, y)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt2007"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment: ckr Doesn't work: fastNLO: error occured - parameter outside range - check commonblock
Comment: ckr      call FT2007CC('tablepath/fnt2007midp.tab',1.0d0 ,1.0d0, 1,XST2007)
Comment: 
Comment: c --- fnt2008 - arXiv:0812.4036 [hep-ex] - CDF dijets, cone algo R=0.7 (Mjj)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt2008"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT2008CC('tablepath/fnt2008midp.tab',1.0d0 ,1.0d0, 1,XST2008)
Comment: 
Comment: c --- fnt2009 - arXiv:0802.2400 [hep-ex] - D0 inclusive jets, cone algo R=0.7 (pT,y)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt2009"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT2009CC('tablepath/fnt2009midp.tab',1.0d0 ,1.0d0, 1,XST2009)
Comment: 
Comment: c --- fnt2010 - arXiv:0906.4819 [hep-ex] - D0 dijets, cone algo R=0.7 (chi, Mjj)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt2010"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT2010CC('tablepath/fnt2010midp.tab',1.0d0 ,1.0d0, 1,XST2010)
Comment: 
Comment: c --- fnt2011 - arXiv:1002.4594 [hep-ex] - D0 dijets, cone algo R=0.7 (Mjj, |y|-max)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnt2011"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FT2011CC('tablepath/fnt2011midp.tab',1.0d0 ,1.0d0, 1,XST2011)
Comment: 
Comment: 
Comment: 
Comment: c =================================================================
Comment: c ======== HERA ===================================================
Comment: c =================================================================
Comment:       WRITE(*,*)"########################################"//
Comment:      >     "################################"
Comment:       WRITE(*,*)"# HERA"
Comment:       WRITE(*,*)"########################################"//
Comment:      >     "################################"
Comment: c - fnh1001 - hep-ex/0010054 - H1 incl. jets (ET, Q2)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnh1001"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FH1001CC('tablepath/fnh1001.tab',1.0d0 ,1.0d0, 1,XSH1001)
Comment: 
Comment: c - fnh1002 - hep-ex/0208037 - ZEUS incl. jets (ET, Q2) 
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnh1002"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FH1002CC('tablepath/fnh1002.tab',1.0d0 ,1.0d0, 1,XSH1002)
Comment: 
Comment: c - fnh1003 - hep-ex/0206029, HERA, H1 incl. jets  (ET) @low Q2
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnh1003"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FH1003CC('tablepath/fnh1003.tab',1.0d0 ,1.0d0, 1,XSH1003)
Comment: 
Comment: c - fnh1004 - hep-ex/0010054, HERA, H1 dijets (ET, Q2) 
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnh1004"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FH1004CC('tablepath/fnh1004.tab',1.0d0 ,1.0d0, 1,XSH1004)
Comment: 
Comment: c - fnh2001 - hep-ex/0608048, HERA, ZEUS incl. jets, kT algo (ET, Q2)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnh2001"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FH2001CC('tablepath/fnh2001.tab',1.0d0 ,1.0d0, 1,XSH2001)
Comment: 
Comment: c - fnh2003 - arXiv:0706.3722 [hep-ex], H1 incl. jets, kT algo (ET, Q2)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnh2003"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FH2003CC('tablepath/fnh2003.tab',1.0d0 ,1.0d0, 1,XSH2003)
Comment: 
Comment: 
Comment: 
Comment: c =================================================================
Comment: c ======== RHIC ===================================================
Comment: c =================================================================
Comment:       WRITE(*,*)"########################################"//
Comment:      >     "################################"
Comment:       WRITE(*,*)"# HERA"
Comment:       WRITE(*,*)"########################################"//
Comment:      >     "################################"
Comment: 
Comment: c --- fnr0001 - STAR preliminary - incl. jets (pT)
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnr0001"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FR0001CC('tablepath/fnr0001midp.tab',1.0d0 ,1.0d0, 1,XSR0001)



c =================================================================
c ======== LHC ====================================================
c =================================================================
      WRITE(*,*)"########################################"//
     >     "################################"
      WRITE(*,*)"# LHC"
      WRITE(*,*)"########################################"//
     >     "################################"

c - fnl0004 - test scenario incl jets (pT,y) optimized for ATLAS
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       WRITE(*,*)"Scenario: fnl0004"
Comment:       WRITE(*,*)"----------------------------------------"//
Comment:      >     "--------------------------------"
Comment:       call FL0004CC('tablepath/fnl0004.tab',1.0d0 ,1.0d0, 1,XSL0004)

c - fnl1013 - arXiv:1104.1693 - LHC, CMS dijet mass (M_JJ,|y|_max) 
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      WRITE(*,*)"Scenario: fnl1013"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      call FL1013CC('tablepath/fnl1013.tab',1.0d0 ,1.0d0, 1,XSL1013)

c - fnl1014 - arXiv:1106.0208 - LHC, CMS inclusive jets (pT,|y|) 
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      WRITE(*,*)"Scenario: fnl1014"
      WRITE(*,*)"----------------------------------------"//
     >     "--------------------------------"
      call FL1014CC('tablepath/fnl1014.tab',1.0d0 ,1.0d0, 1,XSL1013)



      END
