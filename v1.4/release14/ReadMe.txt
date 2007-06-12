***********************************************************************
***********************************************************************
*
*                                         M. Wobisch - June 05, 2006
*
* instructions how to setup and run fastNLO  
*
* fastNLO allows fast computations of cross sections in hadron-induced 
* processes for arbitrary PDFs and alphas values. The calculations
* are based on tables of pre-computed perturbative coefficients.
* Tables and corresponding user-code are available for a large number
* of datasets which can be very helpful in PDF fits to constrain
* the parton distributions.
*
* All code can be obtained from the webpage
*             http://hepforge.cedar.ac.uk/fastnlo
*
*
***********************************************************************
***********************************************************************



***********************************************************************
*         the concept
***********************************************************************

*** fastNLO is structured in "scenarios". 
One scenario computes the results of one measurement, usually 
corresponding to the results of one publication.


*** There is general code - needed for all scenarios:
This are the interfaces to the PDFs and to alphas - contained
in the file "fn-interface.f".
 - If you want to use your own PDF/alphas code, you need 
   to edit "fn-interface.f". (the subroutines and the interface
   are well-desribed)
 - If you want to compute cross sections, based on existing PDFs
   you can use the existing code, but you have to install
   LHAPDF. This is very easy - instructions can be found at
   http://hepforge.cedar.ac.uk/lhapdf/


*** There is scenario-specific code
The subroutines for each scenario are included in the files
which are named based on the scenario-names. In addition
there is one (or in some cases two) include files.
Example:
The usercode for scenario "fnt1001" is contained in the 
files "fnt1001.f" and "fnt1001.inc"



***********************************************************************
*         getting started -> running the example
***********************************************************************

*** install LHAPDF
  The fastNLO example requires LHAPDF to be installed. 
  If this is not yet installed on your system, simply
  follow the instructions at http://hepforge.cedar.ac.uk/lhapdf/


*** get the fastNLO code from: http://hepforge.cedar.ac.uk/fastnlo
  Make sure you get the common routines in "fn-common.tar.gz".
  This file contains:
     ReadMe.txt         (this file)
     fn-alphas-demo.f   (a simple alphas routine which works only
                         above the b-quark mass-threshold)
     fn-interface.f     (the interface to the PDFs and alphas -
                         as default LHAPDF and fn-alphas-demo.f are used)
     Makefile
     fn-example.f       (an example routine which demonstrates how to
                        call the usercode)

  In addition, you need the code for at least one scenario
  "ft[xnnnn]-code.tar.gz" where [xnnnn] is scenario-specific, as
  documented on the webpage. This contains:
     fn[xnnnn].f        - code
     fn[xnnnn].inc      - commonblock definition
 
  Also get the coefficient table for the desired scenario(s)
  "fn[xnnn].tab.gz"


*** installation
  - create a new directory (e.g.: "fastnlo") and copy
    your fastNLO code into this directory

  - inside the directory "fastnlo" create two links:
    - one link "tablepath" pointing to a directory where
      you store the (large!) coefficient tables
    - one link "pdfpath" pointing to the LHAPDF directory
      where the .grid and .pdf files from LHAPDF are stored
      (/share/lhapdf/PDFsets/ -> inside your LHAPDF directory) 

  - edit the "Makefile" to set the two paths to LHAPDF and CERNLIB
    (and delete the filenames for the scenarios for which you have 
    not downloaded the code)

  - unless you downloaded all usercode and tables from the webpage:
    uncomment in "fn-example.f" the scenarios that you don't want.

  - type "make"

  - run "./fastnlo"
    -> the numbers that you get should be identical to the numbers
       that you can compute online using the webinterface at
       http://hepforge.cedar.ac.uk/fastnlo


***********************************************************************
*         getting started -> interface your own PDFs
***********************************************************************

*** if you want to use your own PDFs and/or alphas routine,
  you can edit "fn-interface.f". There is one function
      DOUBLE PRECISION FUNCTION FNALPHAS(MUR)
  which returns the value of alphas/(2pi)
  and a subroutine
      SUBROUTINE FNPDF(X,MUF,XPDF)
  which returns the parton densities, multiplied by x
  using the LHAPDF numbering convention:
     tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t
      -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 , 0, 1, 2, 3, 4, 5, 6


***********************************************************************
* questions? please contact the authors at:   fastnlo@cedar.ac.uk
***********************************************************************




