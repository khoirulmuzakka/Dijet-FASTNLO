#!/bin/csh -efx
###############################################################################
# QUICK INSTALL
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
# PRECONDITIONS
#
# This scripts expects all properly named (see the "set arc=" statements below)
# source code packages to be available in the current working directory.
# In case of a successful install a correspondingly named empty file is created.
# When re-running this script the presumably successful installs are skipped.
#
# The configure step checks whether all pre-conditions are met and
# here it is, where things can get complicated. There are two categories
# of items which might be missing:
#
# 1. Essential development tools are missing from the system installation or
#    are not available in a compatible or recent enough version.
#    Examples are the GCC compiler collection, an MPI installation like OpenMPI or
#    MPICH(2), the GNU scientific library GSL, or the Boost C++ libraries.
#    Most problems encountered here are NOT because of fastNLO itself, but
#    because of requirements of other packages!
#    ==> Possible solutions:
#        a) Use your system tools to install and configure the required packages.
#        b) Use environment modules available on your system to activate
#           adequate software versions via "module load module/name/version".
#           (Use "module avail" to check on availability.)
#
# 2. Special requirements of a fastNLO package are not met, e.g. the interface
#    package to NLOJet++ needs, of course, NLOJet++ to be installed first.
#    ==> Solution: Download the necessary packages either from this web page
#        http://fastnlo.hepforge.org/code/v23/code-v23.html
#        or from the respective projects web page, frequently also hosted
#        at http://www.hepforge.org
#        But please note that for NLOJet++ to work with fastNLO the patched
#        version from our HepForge web page is mandatory!
#
#        Then read on for recommended install options. For these we assume
#        that all 'local' packages of this step are unpacked in
#        $HOME/local/src
#        and that the common install directory will be
#        $HOME/local
#        Do NOT install a source package into the directory, where it is built,
#        i.e. src in this case. This violates principles of system building.
#        Otherwise, any local directory with sufficient space is perfectly acceptable.
#        If it is intended to use this installation from the worker nodes of
#        a batch system, then the install directory must either be visible from these
#        worker nodes, or the installation must somehow be transferred. The latter
#        depends, of course, on the individual grid/batch system and cannot be
#        described here.
#
#==============================================================================
# IMPORTANT REMARKS: Read before installing.
#
# 1. The install directory assumed here is $HOME/local with binaries to be
#    installed in $HOME/local/bin. Make sure that
#    - $HOME/local/bin is in the search path to find commands like "toolname"-config,
#    - the right versions are used, either the desired system ones or the ones from
#      $HOME/local/bin.
#    Usually, this is set up via defining the PATH environment variable appropriately.
#
# 2. If you want to install BlackHat 0.9.9 for use within Sherpa-->MCGrid-->fastNLO:
#    - Compilation with gcc-4.8.x gives an error -->
#      Use older version, i.e. gcc, g++, and gfortran 4.7.x
#      In Ubuntu, this can be achieved via "sudo update-alternatives --config gcc" etc.,
#      or switch to version 4.7.x using the corresponding "module load" command.
#    - Install the libssl-dev (Ubuntu) to resolve potentially missing dependency
#    - The qd-2.3.17 package required by BlackHat does not necessarily choose your
#      system default compiler. To avoid e.g. using the Intel compiler configure
#      this package with explicit options "CXX=g++ CC=gcc FC=gfortran".
#      If you DO use the Intel compiler, then the file qd_real.cpp has to be
#      made with "make CXXFLAGS=-O1" to avoid infinite loops.
#
# 3. If you want to install OpenLoops 1.1.1 for use within Sherpa-->MCGrid-->fastNLO:
#    - gcc >= 4.6 is required!
#
# 4. If you switch to a newer compiler via the "module" command,
#    "module load compiler/gnu/4.7"
#    it might be necessary to also activate a newer boost library (>= 1.48.0) used e.g.
#    by the LHAPDF, YODA (1.4.x < YODA < 1.6.x), and RIVET (2.3.x < Rivet < 2.5.x) packages:
#    "module load lib/boost/1.56.0"
#    In that case configure the packages by specifying additional include
#    paths with CPPFLAGS="-I${INCLUDE}" or more explicitly
#    by CPPFLAGS="-I/add/include/path1 -I/add/include/path2".
#
# 5. The latest versions YODA 1.6.x and Rivet 2.5.x are not yet supported by fastNLO!
#    If the optional YODA/Rivet support of the fastNLO Toolkit is requested to be
#    installed by specifying "--with-yoda" below, it must be ensured that previous
#    versions of these tools are used.
#
# 6. SHERPA can be used with MPI for parallel processing if desired. In that case,
#    the proper MPICH or OpenMPI software needs to be installed or activated, e.g.
#    via "module load mpi/openmpi/1.6.5-gnu-4.7". Do not forget to specify corresponding
#    configuration options.
#
#==============================================================================
#
#
#==============================================================================
# Check and set defaults for some specific environment variables
#==============================================================================
# Environment settings e.g. via "module load" might set a global INCLUDE variable,
# for example to include BOOST headers! By default it does not exist.
setenv MYCPPFLAGS
if ( $?INCLUDE ) then
   if ( ! $?CPPFLAGS ) then
      setenv MYCPPFLAGS "-I${INCLUDE}"
   else
      setenv MYCPPFLAGS "${CPPFLAGS} -I${INCLUDE}"
   endif
endif
# If MPI is used, the installation path might be set here
if ( $?WITH_MPI ) then
  if ( $?MPI_HOME ) then
    setenv WITH_MPI "--enable-mpi=${MPI_HOME}"
  else
    setenv WITH_MPI "--enable-mpi"
  endif
else
  setenv WITH_MPI
endif
#
#==============================================================================
# Mandatory and sincerely recommended packages
#==============================================================================
#
# fastjet (any version >= 3 should work):
# Use the "--enable-allplugins" options to enable use of e.g. older Tevatron jet algorithms.
#------------------------------------------------------------------------------
set arc="fastjet-3.1.3"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --enable-shared --enable-allplugins --prefix=$HOME/local --bindir=$HOME/local/bin
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
# LHAPDF (last v5 or v6; if you also want the old fastNLO_reader package use v5):
# In the near future, support for LHAPDF v5 will be stopped!
#------------------------------------------------------------------------------
# set arc="lhapdf-5.9.1"
set arc="LHAPDF-6.1.6"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=$HOME/local CPPFLAGS="${MYCPPFLAGS}"
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
#
#==============================================================================
# Additional packages required for some of fastNLO optional features
# Skip any of these parts that is not desired.
#==============================================================================
#
# ROOT (e.g. v5.34.25; optional use by YODA):
# If Python support is desired, use "--enable-python".
#------------------------------------------------------------------------------
set arc="root-5.34.25"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  mv root ${arc}
  cd ${arc}
  ./configure --prefix=$HOME/local --etcdir=$HOME/local/etc --enable-python --enable-minuit2
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
# ==> This enables the option --with-root of the fastNLO_toolkit to produce
#     ROOT histograms from the calculated cross sections and uncertainties
#
#
# HepMC (e.g. v2.06.09; needed by Rivet):
#------------------------------------------------------------------------------
set arc="HepMC-2.06.09"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=$HOME/local --with-momentum=GEV --with-length=MM
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
# YODA (must be < 1.6.x; needed by Rivet):
# Python is enabled by default. If Python with ROOT interfacing is desired, use "--enable-root".
#------------------------------------------------------------------------------
set arc="YODA-1.5.9"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=$HOME/local --enable-root CPPFLAGS="${MYCPPFLAGS}"
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
# Rivet (versions >=2 are necessary; v1 is too old; must be < 2.5.x):
# Python is enabled by default.
#------------------------------------------------------------------------------
set arc="Rivet-2.3.0"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=$HOME/local CPPFLAGS="${MYCPPFLAGS}"
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
# ==> This enables the option --with-yoda of the fastNLO_toolkit to produce
#     YODA formatted e.g. NLO predictions for comparison to data with Rivet
#
#
# HOPPET (e.g. v1.1.5):
#------------------------------------------------------------------------------
set arc="hoppet-1.1.5"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=$HOME/local
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
# ==> This enables the option --with-hoppet of the fastNLO_toolkit to use
#     the alpha_s evolutions within HOPPET
#     fastNLO comes already with alpha_s evolutions from GRV or CRunDec, or
#     uses the one from the PDFs in LHAPDF
#
#
# QCDNUM:
# Use the newer autotools-enabled version 17.01.12.
# Older QCDNUM versions including 17.00.06 have an error in the time-like evolution
# of the singlet fragmentation function at NLO, see arXiv:1602.08383.
# Non-autotools-enabled versions are not supported anymore.
#------------------------------------------------------------------------------
set arc="qcdnum-17-01-12"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=$HOME/local
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
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
set arc="fastnlo_toolkit-2.3.1pre-2212"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=$HOME/local --with-yoda --with-hoppet --with-qcdnum --with-root --enable-pyext
# options depending on previous choices: --with-yoda --with-hoppet --with-qcdnum --with-root
# option to use python interface to library: --enable-pyext
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
#
#==============================================================================
# fastNLO use with NLOJet++
#==============================================================================
#
# NLOJet++ (patched version from fastNLO web page is required: NLOJet-4.1.3-patched):
#------------------------------------------------------------------------------
set arc="nlojet++-4.1.3"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}-patched.tar.gz
  cd ${arc}
  ./configure --prefix=$HOME/local
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
# fastNLO Interface NLOJet++:
#------------------------------------------------------------------------------
set arc="fastnlo_interface_nlojet-2.3.1pre-2125"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=$HOME/local
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
#
#==============================================================================
# Old fastNLO Reader package including Fortran and deprecated C++ code
#==============================================================================
#
# fastNLO Reader:
# Only works with LHAPDF v5!
#------------------------------------------------------------------------------
# set arc="fastnlo_reader-2.1.0-2066"
# if ( ! -e ${arc}_installed  ) then
#   tar xzf ${arc}.tar.gz
#   cd ${arc}
#   ./configure --prefix=$HOME/local --enable-pyext --with-hoppet --with-qcdnum
#   make -j4 install
#   cd ..
#   touch ${arc}_installed
# endif
#
#
#==============================================================================
# fastNLO Converter package to convert old v1.4 tables to new v2 format
#==============================================================================
#
# fastNLO Converter:
#------------------------------------------------------------------------------
# set arc="fastnlo_converter-2.1.0-2107"
# if ( ! -e ${arc}_installed  ) then
#   tar xzf ${arc}.tar.gz
#   cd ${arc}
#   ./configure --prefix=$HOME/local
#   make -j4 install
#   cd ..
#   touch ${arc}_installed
# endif
#
#
#==============================================================================
# fastNLO use with Sherpa & MCgrid
#==============================================================================
#
# QD library required by BlackHat:
#------------------------------------------------------------------------------
#set arc="qd-2.3.17"
#if ( ! -e ${arc}_installed  ) then
#  tar xzf ${arc}.tar.gz
#  cd ${arc}
#  ./configure --prefix=$HOME/local --enable-shared CXX=g++ CC=gcc FC=gfortran
#  make -j4 install
#  cd ..
#  touch ${arc}_installed
#endif
#
# BlackHat for use within Sherpa:
# Compile error for gcc > 4.7!
#------------------------------------------------------------------------------
# Avoid redoing BlackHat when the directory exists already assuming
# the installation was successful.
# If not treat the blackhat heavy stuff by hand!
#set arc="blackhat-0.9.9"
#if ( ! -e ${arc}_installed  ) then
#  tar xzf ${arc}.tar.gz
#  cd ${arc}
#  ./configure --prefix=$HOME/local
## Make without -j4 multicore, might be too heavy on some machines
#  make install
#  cd ..
#  touch ${arc}_installed
#endif
#
# NJet for use within Sherpa:
#------------------------------------------------------------------------------
#set arc="njet-2.0.0"
#if ( ! -e ${arc}_installed  ) then
#  tar xzf ${arc}.tar.gz
#  cd ${arc}
#  ./configure --prefix=$HOME/local
#  make -j4 install
#  cd ..
#  touch ${arc}_installed
#endif
#
# OpenLoops for use within Sherpa:
# Needs gcc >= 4.6!
#------------------------------------------------------------------------------
# set arc="OpenLoops-1.1.1"
#   if ( ! -e ${arc}_installed  ) then
#   tar xzf ${arc}.tar.gz
#   cd ${arc}
#   ./scons
#   cd ..
#   touch ${arc}_installed
# endif
#
# Sherpa for use with MCgrid & fastNLO:
# Version >= 2.2.0 is required!
# For use of MPI with Sherpa corresponding system packages are required (MPICH, OpenMPI).
#------------------------------------------------------------------------------
#set arc="SHERPA-MC-2.2.0"
#if ( ! -e ${arc}_installed  ) then
#  tar xzf ${arc}.tar.gz
#  cd ${arc}
#  ./configure --prefix=$HOME/local --with-sqlite3=install --enable-gzip --enable-lhole ${WITH_MPI} --enable-fastjet=$HOME/local --enable-lhapdf=$HOME/local --enable-hepmc2=$HOME/local --enable-rivet=$HOME/local --enable-root=$HOME/local --enable-blackhat=$HOME/local --enable-openloops=$HOME/local
#  make -j4 install
#  cd ..
#  touch ${arc}_installed
#endif
#
# MCGrid for use with fastNLO:
# (Seems to be fixed: Requires CXXFLAGS+=-fpermissive flag in make step with gcc-4.7.x)
#------------------------------------------------------------------------------
#set arc="mcgrid-2.0"
#if ( ! -e ${arc}_installed  ) then
#  tar xzf ${arc}-fixed.tar.gz
#  cd ${arc}
#  ./configure --prefix=$HOME/local
#  make -j4 install
## CXXFLAGS+=-fpermissive
#  cd ..
#  touch ${arc}_installed
#endif
#
# For usage, make sure to set PATH, PYTHON, LHAPDF, RIVET, PKG_CONFIG_PATH, ROOT environment variables appropriately!!!
#
# setenv PKG_CONFIG_PATH $HOME/local/lib/pkgconfig
# source rivetanalysis.csh
# add mcgrid example dirs
#
#
#==============================================================================
# fastNLO use with Herwig7 (NOT YET IMPLEMENTED)
#==============================================================================
#
# ThePEG:
#--------
#autoreconf -i
#./configure --prefix=$HOME/local --with-fastjet=$HOME/local --with-hepmc=$HOME/local --with-lhapdf=$HOME/local --with-rivet=$HOME/local
#make -j8
#make check
#make install
#
# NJet-matchbox:
#---------------
#./configure --prefix=$HOME/local --enable-5jet
#make -j8
#make check
#make install
#
# Herwig++-matchbox:
#-------------------
#install gengetopt from distro
#autoreconf -i
#./configure --prefix=$HOME/local --with-nlojet=$HOME/local --with-njet=$HOME/local --with-fastjet=$HOME/local
#./configure --prefix=$HOME/local --with-fastjet=$HOME/local --with-madgraph=$HOME/local/MadGraph-matchbox --with-njet=$HOME/local
#make -j8
#make check
#make install
#make  check-local
