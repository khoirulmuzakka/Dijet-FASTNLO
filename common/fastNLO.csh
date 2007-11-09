#
# Set some environment variables necessary for fastNLO
#
# Location of CERNLIB, LHAPDF and NLOJET
#
if ( $1 == "kr" || $1 == "grid" ) then
    setenv FASTNLO $HOME/local
    setenv CERNLIB $FASTNLO/cernlib 
    setenv LHAPDF  $FASTNLO/lhapdf/lib
    setenv NLOJET  $FASTNLO/nlojet
    setenv CXXFLAGS "-O3 -I ."
endif
if ( $1 == "mw" ) then  
    setenv CERNLIB /D0/ups/cern/Linux-2-4/2002/lib
    setenv LHAPDF /disk2/work/wobisch/LHAPDFv4
    setenv NLOJET /home/wobisch/NLOJET
    setenv CXXFLAGS "-O3 -I ."
endif
if ( $1 == "tk" ) then  
    setenv CERNLIB /opt/products/cernlib/pro/lib
    setenv LHAPDF /h1/h1gen/lhapdf/lib
    setenv NLOJET /scratch/usr/kluget/nlojet++-2.0.1
    setenv CXXFLAGS "-O3 -I ."
endif
