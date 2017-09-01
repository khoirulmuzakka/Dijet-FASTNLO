#!/bin/csh -f
###############################################################################
#
# Installing any fastNLO package is, in principal, quite easy and
# requires only three steps:
#
# 1. Unpacking the tarball in some source directory:
#    tar xzf fastnlo_toolkit-v.v.vpre-rrrr.tar.gz
# 2. cd into the created directory and run configure:
#    cd fastnlo_toolkit-v.v.vpre-rrrr
#    ./configure --prefix=/your/install/directory
# 3. make -j ncores install
#
###############################################################################
#
# The configure step checks whether all pre-conditions are met and
# here it is, where things can get complicated. There are two categories
# of items which might be missing:
#
# 1. Essential development tools like a C++ compiler are missing from
#    from the system installation.
#    ==> Solution: Use your system tools to install the required packages.
#
# 2. Special requirements of a fastNLO package are not met, e.g. the
#    interface package to NLOJet++ needs NLOJet++ to be installed first.
#    ==> Solution: Download the necessary packages either from this web page
#        http://fastnlo.hepforge.org/code/v23/code-v23.html
#        or from the respective projects web page, frequently also hosted
#        at http://www.hepforge.org
#
#        Then read on for recommended install options. For these we assume
#        that all 'local' packages of this step are unpacked in
#        $HOME/local/src
#        and that the common install directory will be
#        $HOME/local
#
#
#
#==============================================================================
# Mandatory and sincerely recommended packages
#==============================================================================
#
# fastjet (any version >= 3 should work):
#------------------------------------------------------------------------------
tar xzf fastjet-3.1.2.tar.gz
cd fastjet-3.1.2
./configure --enable-shared --enable-allplugins --prefix=$HOME/local --bindir=$HOME/local/bin
make -j4 install
cd ..
#
# LHAPDF (last v5 or v6; if you also want the old fastNLO_reader package use v5):
#------------------------------------------------------------------------------
tar xzf lhapdf-5.9.1.tar.gz
cd lhapdf-5.9.1
./configure --prefix=$HOME/local
make -j4 install
cd ..
#
#
#==============================================================================
# Additional packages required for some of fastNLO optional features
# Skip any of these parts that is not desired.
#==============================================================================
#
# HepMC (e.g. v2.06.09; needed by Rivet):
#------------------------------------------------------------------------------
tar xzf HepMC-2.06.09.tar.gz
cd HepMC-2.06.09
./configure --prefix=$HOME/local --with-momentum=GEV --with-length=MM
make -j4 install
cd ..
#
# YODA (e.g. v1.3.1; needed by Rivet):
#------------------------------------------------------------------------------
tar xzf YODA-1.3.1.tar.gz
cd YODA-1.3.1
./configure --prefix=$HOME/local
make -j4 install
cd ..
#
# Rivet (versions >=2 are necessary; v1 is too old):
#------------------------------------------------------------------------------
tar xzf Rivet-2.2.1.tar.gz
cd Rivet-2.2.1
./configure --prefix=$HOME/local
make -j4 install
cd ..
#
# ==> This enables the option --with-yoda of the fastNLO_toolkit to produce
#     YODA formatted e.g. NLO predictions for comparison to data with Rivet
#
#
# HOPPET (e.g. v1.1.5):
#------------------------------------------------------------------------------
tar xzf hoppet-1.1.5.tgz
cd hoppet-1.1.5
./configure --prefix=$HOME/local
make -j4 install
cd ..
#
# ==> This enables the option --with-hoppet of the fastNLO_toolkit to use
#     the alpha_s evolutions within HOPPET
#     fastNLO comes already with alpha_s evolutions from GRV or CRunDec, or
#     uses the one from the PDFs in LHAPDF
#
#
# QCDNUM (patched version from fastNLO web page is required: QCDNUM-17.00.06-patched):
#------------------------------------------------------------------------------
# (Alternatively to using the patched version, all makelib commands have to be
#  changed to use the -fPIC option in the gfortran compile statement.)
tar xzf qcdnum-17.00.06-patched.tar.gz
cd qcdnum-17-00-06
./makelibs
cp -p lib/lib*.a $HOME/local/lib
cd ..
#
# ==> This enables the option --with-qcdnum of the fastNLO_toolkit to use
#     the alpha_s evolutions within QCDNUM
#
#
#==============================================================================
# The fastNLO Toolkit
#==============================================================================
#
# fastNLO Toolkit:
#------------------------------------------------------------------------------
#tar xzf fastnlo_toolkit-v.v.vpre-rrrr.tar.gz
#cd fastnlo_toolkit-v.v.vpre-rrrr
#./configure --prefix=$HOME/local
# options depending on previous choices: --with-yoda --with-hoppet --with-qcdnum
# option to use python interface to library: --enable-pyext
#make -j4 install
#
#
#==============================================================================
# fastNLO use with NLOJet++
#==============================================================================
#
# NLOJet++ (patched version from fastNLO web page is required: NLOJet-4.1.3-patched):
#------------------------------------------------------------------------------
tar xzf nlojet++-4.1.3-patched.tar.gz
cd nlojet++-4.1.3
./configure --prefix=$HOME/local
make -j4 install
cd ..
#
# fastNLO Interface NLOJet++:
#------------------------------------------------------------------------------
#tar xfz fastnlo_interface_nlojet-v.v.vpre-rrrr.tar.gz
#cd fastnlo_interface_nlojet-v.v.vpre-rrrr
#./configure --prefix=$HOME/local
#make -j4 install
#
#
#==============================================================================
# Old fastNLO Reader package including Fortran and deprecated C++ code
#==============================================================================
#
# fastNLO Reader:
#------------------------------------------------------------------------------
tar xzf fastnlo_reader-2.1.0-2066.tar.gz
cd fastnlo_reader-2.1.0-2066
./configure --prefix=$HOME/local --enable-pyext --with-hoppet --with-qcdnum
make -j4 install
cd ..
#
#
#==============================================================================
# fastNLO use with Sherpa & MCgrid
#==============================================================================
#
# MCgrid not released yet, please be patient.
#
