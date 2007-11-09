#
# Set some environment variables necessary for fastNLO
#
# Location of CERNLIB, LHAPDF and NLOJET
#
if [ "$1" = "kr" -o "$1" = "grid" ] ; then
    export FASTNLO=$HOME/local
    export CERNLIB=$FASTNLO/cernlib 
    export LHAPDF=$FASTNLO/lhapdf/lib
    export NLOJET=$FASTNLO/nlojet
    export CXXFLAGS="-O3 -I ."
fi
if [ "$1" = "mw" ] ; then  
    export SVN_EDITOR="emacs -nw"
    export CERNLIB=/usr/lib/cernlib/2005/lib
    export LHAPDF=/usr/local/lib
    export NLOJET=/home/wobisch/NLOJET
    export CXXFLAGS="-O3 -I ."
fi
if [ "$1" = "tk" ] ; then  
    export CERNLIB=/opt/products/cernlib/pro/lib
    export LHAPDF=/h1/h1gen/lhapdf/lib
    export NLOJET=/scratch/usr/kluget/nlojet++-2.0.1
    export CXXFLAGS="-O3 -I ."
fi
