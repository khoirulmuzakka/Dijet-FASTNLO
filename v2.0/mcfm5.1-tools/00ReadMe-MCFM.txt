
*******************************************************************
*  fastNLO interface to MCFM                 M. Wobisch 06/19/2007
*
*  these fastNLO routines were made for MCFM v5.1 
*                       (version from June 1, 2006) 
*
*
*******************************************************************

********************* Installation ********************************
* get the MCFM file from http://mcfm.fnal.gov  
  (and hope that you have uudecode on your machine...)

  Text from the MCFM webpage:
  The source is available as a tar'ed, gzip'ed, uu-encoded package. 
  This should be placed in a new directory and then extracted. The 
  source can be initialized by running the Install command and then 
  compiled with make. ... Execution should be carried out in the Bin 
  directory, whilst some documentation (for example, instructions 
  regarding the format for the input file) can be found in Doc.

  comment MW: the installation worked fine for me once I replaced
  the compilerflag to "FC = f95" - but might not be necessary for 
  most users. If emacs tells you that there are suspicious line in 
  the makefile just ignore this.


* add the following routines to the MCFM /Inc and /User directories:
  cp svn/fastNLO/trunk/v2.0/mcfm5.1-tools/fastnlo.f     mcfm/src/Inc/.
  cp svn/fastNLO/trunk/v2.0/mcfm5.1-tools/fn-author.f   mcfm/src/User/.
  (fastnlo.f defines the fastNLO commonblock / fn-author.f has the code)

* replace the following MCFM routines with those from the fastNLO 
  SVN repository (we assume that the SVN files are in the local 
  directory "svn" and the MCFM installation is in the local 
  directory "mcfm"

  cp svn/fastNLO/trunk/v2.0/mcfm5.1-tools/nplotter.f        mcfm/src/User/.
  cp svn/fastNLO/trunk/v2.0/mcfm5.1-tools/lowint_incldip.f  mcfm/src/Need/.
  cp svn/fastNLO/trunk/v2.0/mcfm5.1-tools/realint.f         mcfm/src/Need/.
  cp svn/fastNLO/trunk/v2.0/mcfm5.1-tools/virtint_incldip.f mcfm/src/Need/.
  cp svn/fastNLO/trunk/v2.0/mcfm5.1-tools/mcfm_exit.f       mcfm/src/Need/.

  all additional lines in the modified MCFM files are enclosed in:
   c === fastNLO
         ...
   c ====================================================

  nplotter.f is the MCFM user interface -> added calls to FnInit & FnInterface
  mcfm_exit.f make the MCFM output -> added call to FnTerm
  lowint_incldip.f computes the LO x-sect -> copy coefficients in fN array
  realint.f computes the real corrections -> copy coefficients in fN array
  virtint_incldip.f computes the virtual corrections -> copy coefficients in fN array


* edit makefile: 
  we have used "indent" in emacs to make the logical structure of
  "virtint_incldip.f" more clearly visible - this requires to allow
  a longer line-length for the Fortran compiler: 
    #FFLAGS = -fno-automatic -fno-f2c -O0 -g -I$(INCPATH)
    FFLAGS  = -fno-automatic -ffixed-line-length-132 -fno-f2c -O0 -g -I$(INCPATH)
 
  under USERFILES (line 340) add the following entries in the makefile:
  fn-author.o


********************* Implementation: fn-author.f *************************

* Subroutine FnInit


* Subroutine FnTerm
  - output of reference values

* Subroutine FnInterface(wt)
  - sum up cross-sections to check if we got all contributions



********************** Examples & References *****************************
*** obtained with MFMC's internal CTEQ6M PDFs

* subproc 1 (W) for Tevatron Run II
lord     9860371.575 +/-  8547.403 fb
real      175183.404 +/-  2890.431 fb
virt    11730783.538 +/- 11299.466 fb
tota    11929392.742 +/- 11406.520 fb       = 11.9nb
--> Sum of Virt & Real gives total 
 -> virtual includes LO 
 -> how can we disentangle contributions with different orders in alphas?


