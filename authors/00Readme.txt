***********************************************************
* June 24, 2005        T. Kluge, K. Rabbertz, M. Wobisch
*
* FNd0run1.tar
* package with author's code for D0 Run I measurement of 
* the inclusive jet cross section  (hep-ex/0011036)
*
* This file is part of the package  FNd0run1.tar.gz
***********************************************************


* ---------------------------------------------------------
* The following files are included in this package:
* ---------------------------------------------------------
 00Readme.txt           (this file)
 01ChangeLog.txt        history of modifications
 Makefile               Makefile
 fn-d0run1v00.cc        NLOJET++ user routine for D0 Run I measurement 
                         of the incl jet x-sect
*** PDF code:
 cteq6.cc               NLOJET++ interface to CTEQ code
 cteq6.h                    -''-
 Cteq61Pdf.f            original CTEQ code for the 6.1 fit
 ctq61.00.tbl           original CTEQ table for the CTEQ6.1M parametrization

*** jet algorithms:
 d0run1cone.cc          implementation of the Run I algorithm 
                          (ET-scheme, Rsep=1.3) for 2-jet NLO 
 d0run1cone.h               -''-
 --- the following are not yet used - but already provided:
 d0run2cone.cc          implementation of the Run II algorithm 
                          (E-scheme, no Rsep) for 2-jet NLO 
 d0run2cone.h               -''-
 d0run2cone3.cc         interface to PXCONE (needed for 3-jet NLO in Run II)
 d0run2cone3.h              -''-
 pxcone_mod_new.f       PXCONE algorithm - modified for Run II (E-scheme)
 px_interface.f         interface to PXCONE


* ---------------------------------------------------------
* step 1: installing a slightly modified NLOJET++ version
* ---------------------------------------------------------
To use all features of the new fastNLO code you need to
install a slightly modified version of NLOJET++
(includes only one small modification - see below):
* create a new directory: NLOJET
* get the code from Zoltan Nagy's webpage  
  http://www.cpt.dur.ac.uk/~nagyz/nlo++/
  (so far fastNLO is using NLOJET++ version 2.0.1)
* copy the file nlojet++-2.0.1.tar.gz into the new directory NLOJET 
  and unpack the code using:
       cat nlojet++-2.0.1.tar.gz | gunzip | tar xvf -
* a minor modification of the NLOJET++ code is now needed:
  change line 40 in the file:  include/bits/nlo-basic_user.h
  from:        void phys_output(const std::basic_..... 
  to:          virtual void phys_output(const std::basic_..... 
* compile NLOJET++ (may take 15 minutes) using:
       mkdir bld-nlojet 
       cd bld-nlojet
       ../nlojet++-2.0.1/configure --prefix [complete path of your NLOJET directory]
       make install CFLAGS="-O3 -Wall" CXXFLAGS="-O3 -Wall"


* ---------------------------------------------------------
* step 2: installing the user code for the D0 Run I data
* ---------------------------------------------------------

* create a new subdirectory "fastNLO" in your NLOJET directory:
       mkdir fastNLO
* Get the FNd0run1.tar.gz package from the fastNLO author's webpage
     http://www-d0.fnal.gov/~wobisch/fastNLO/authorsonly/
  save it in the "fastNLO" subdirectory and do the following:
        gunzip FNd0run1.tar.gz
        tar -xvf FNd0run1.tar
* compile the fstNLO user code
        make d0run1v00
* run the code:
  ../bin/nlojet++ --save-after 100000 -cborn -n d0run1v00-00 
          -d [directory for the tables] -u fnd0run1v00.la

* comments:
  - "fnd0run1v00.la" is the library that you created using "make"
  - "d0run1v00-00" is the prefix for the name of the output-tables 
  - "born" is for the LO jobs 0- to be replaced 
    by "nlo" is for NLO jobs
  - Recommended numbers for the "--save-after" option are
    50000000 for LO and 5000000 for NLO jobs 
    (this stores the table every 40 minutes).
