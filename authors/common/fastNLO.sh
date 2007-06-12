#
# Set some environment variables necessary for fastNLO
#
# Location of CERNLIB, LHAPDF and NLOJET
#
if [ "$1" = "kr" ]; then
    export CERNLIB=/opt/cern/2003/lib
    export LHAPDF=$HOME/LHAPDFv4
    export NLOJET=$HOME/nlolib
    export CXXFLAGS="-O3 -I ."
fi
if [ "$1" = "mw" ]; then
    export CERNLIB=/D0/ups/cern/Linux-2-4/2002/lib
    export LHAPDF=/disk2/work/wobisch/LHAPDFv4
    export NLOJET=/home/wobisch/NLOJET
    export CXXFLAGS="-O3 -I ."
fi
if [ "$1" = "tk" ]; then
    export CERNLIB=/opt/products/cernlib/pro/lib
    export LHAPDF=/h1/h1gen/lhapdf/lib
    export NLOJET=/scratch/usr/kluget/nlojet++-2.0.1
    export CXXFLAGS="-O3 -I ."
fi
