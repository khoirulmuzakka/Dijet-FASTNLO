#
# Set some environment variables necessary for fastNLO
#
# Location of CERNLIB, LHAPDF and NLOJET
#
if ( $1 == "kr" ) then
    setenv CERNLIB /opt/cern/2003/lib 
    setenv LHAPDF /home/rabbertz/LHAPDFv4
    setenv NLOJET /home/rabbertz/nlolib/src/Nlojet++
    setenv CXXFLAGS "-O3 -I ."
endif
if ( $1 == "mw" ) then  
    setenv CERNLIB /D0/ups/cern/Linux-2-4/2002/lib
#    setenv LHAPDF /disk2/work/wobisch/LHAPDFv4
    setenv LHAPDF /disk2/work/wobisch/LHAPDFv4.1/lib
    setenv NLOJET /home/wobisch/NLOJET
    setenv CXXFLAGS "-O3 -I ."
endif
if ( $1 == "tk" ) then  
    setenv CERNLIB /opt/products/cernlib/pro/lib
    setenv LHAPDF /h1/h1gen/lhapdf/lib
    setenv NLOJET /scratch/usr/kluget/nlojet++-2.0.1
    setenv CXXFLAGS "-O3 -I ."
endif
