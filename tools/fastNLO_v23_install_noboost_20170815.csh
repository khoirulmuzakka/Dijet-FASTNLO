#!/bin/csh -efx
###############################################################################
# QUICK INSTALL, KR 18.08.2017
#
# Depending on the purpose of an installation, up to two fastNLO packages are
# necessary:
#
# 1) The fastNLO Toolkit
#       fastnlo_toolkit-v.v.vpre-rrrr.tar.gz
#    for the evaluation of any fastNLO table
# 2) The fastNLO Interface to NLOJet++
#       fastnlo_interface_nlojet-v.v.vpre-rrrr.tar.gz
#    for the creation of LO or NLO fastNLO tables with NLOJet++
#
# Further programs with other processes up to NNLO have been linked to fastNLO.
# The corresponding tables can be evaluated with this Toolkit. Their creation
# will be dealt with in more detailed instructions in the near future.
#
# Installing any fastNLO package is, in principal, quite easy and
# requires only three steps:
#
# 1. Unpacking the tarball in some source directory:
#    tar xzf fastnlo_*-v.v.vpre-rrrr.tar.gz
# 2. cd into the created directory and run configure:
#    cd fastnlo_*-v.v.vpre-rrrr
#    ./configure --prefix=/your/install/directory
# 3. make -j ncores install
# 4. Optional: make check
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
#    Starting with revision 2411 of the fastNLO packages, full C++11 compatibility
#    of the compiler is enforced, which requires at least version 4.8.1 of the
#    GNU compiler collection.
#    To avoid recurring issues with incompatible versions of the Boost C++ libraries,
#    it is strongly recommended to use LHAPDF version 6.2.0 or newer and, if desired,
#    YODA >= version 1.6.7 and Rivet >= version 2.5.4. Using these or newer versions
#    the use of the Boost C++ libraries can be avoided altogether!
#
#    Most problems encountered here are NOT because of fastNLO itself, but
#    because of requirements of other packages! For example, fastNLO itself
#    never needed Boost ...
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
#        But please note that for NLOJet++ to work with fastNLO, the patched
#        version from our HepForge web page is mandatory!
#
#        Then read on for recommended install options. For these we assume
#        that all 'local' packages of this step are unpacked in
#        $HOME/${local}/src
#        and that the common install directory will be
#        $HOME/${local}
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
# 1. The install directory assumed here is $HOME/${local} with binaries to be
#    installed in $HOME/${local}/bin. Make sure that
#    - $HOME/${local}/bin is in the search path to find commands like "toolname"-config,
#    - the right version of "toolname"-config is used, either the desired system one or
#      the one from $HOME/${local}/bin.
#    Usually, this is set up via defining the PATH environment variable appropriately.
#
# 2. If you want to install BlackHat 0.9.9 for use within Sherpa-->MCGrid-->fastNLO:
#    THIS IS CURRENTLY NOT SUPPORTED ANYMORE!
#    Please wait for news or contact the fastNLO authors.
##    - Compilation with gcc-4.8.x gives an error -->
##      You must use older versions of fastNLO and of gcc, g++, gfortran, e.g. 4.7.x
##      In Ubuntu, the latter can be achieved via "sudo update-alternatives --config gcc" etc.,
##      or switch to version 4.7.x using the corresponding "module load" command.
##    - Install the libssl-dev (Ubuntu) to resolve potentially missing dependency
##    - The qd-2.3.17 package required by BlackHat does not necessarily choose your
##      system default compiler. To avoid e.g. using the Intel compiler configure
##      this package with explicit options "CXX=g++ CC=gcc FC=gfortran".
##      If you DO use the Intel compiler, then the file qd_real.cpp has to be
##      made with "make CXXFLAGS=-O1" to avoid infinite loops.
#
# 3. If you want to install OpenLoops 1.1.1 for use within Sherpa-->MCGrid-->fastNLO:
#    THIS IS CURRENTLY NOT SUPPORTED ANYMORE!
#    Please wait for news or contact the fastNLO authors.
##    - gcc >= 4.8 is required!
#
# 4. If you switch to a newer compiler via the "module" command,
#    "module load compiler/gnu/4.8"
#    and you still use Boost-dependent packages it might be necessary to
#    also activate a newer boost library (>= 1.48.0) used e.g. by older LHAPDF (< 6.2.0),
#    YODA (1.4.x < YODA < 1.6.x), or RIVET (2.3.x < Rivet < 2.5.x) packages:
#    "module load lib/boost/1.56.0"
#    In that case configure the packages by specifying additional include
#    paths with CPPFLAGS="-I${INCLUDE}" or more explicitly
#    by CPPFLAGS="-I/add/include/path1 -I/add/include/path2".
#
# 5. SHERPA can be used with MPI for parallel processing if desired. In that case,
#    the proper MPICH or OpenMPI software needs to be installed or activated, e.g.
#    via "module load mpi/openmpi/1.6.5-gnu-4.7". Do not forget to specify corresponding
#    configuration options.
#
#==============================================================================



# NEW!!!!!!!!!!!!!!!!!!
#==============================================================================
# Check command line arguments
#==============================================================================
if ($#argv != 2 && $#argv != 3) then
   echo "Usage: $0 basedir subdir [path to cvmfs software]"
   echo "1st argument: Base dir for installations"
   echo "2nd argument: Subdir in base dir for this installation"
   echo "3rd optional argument: Path to newer gcc in cvmfs software distribution,"
   echo "                       e.g. /cvmfs/cms.cern.ch/slc6_amd64_gcc481"
   exit 0
endif
#==============================================================================
# Set installation location
#==============================================================================
set base=$1
set local=$2
# NNLOJET
#set revision=3738
#
#==============================================================================
# Check and set defaults for some specific environment variables
#==============================================================================
# Environment settings e.g. via "module load" might set a global INCLUDE variable.
# By default it does not exist.
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
# Use CVMFS software repository, if necessary
if ($3 != "") then
   setenv MYCVMFS $3
   echo "MYCVMFS is $MYCVMFS"
else
   unsetenv MYCVMFS
   echo "MYCVMFS is not set"
endif
# Find root-config
# Do not set ROOTBIN to crappy CVMFS ROOT installation!
#    setenv ROOTBIN ${MYCVMFS}/cms/cmssw/CMSSW_7_2_4/external/slc6_amd64_gcc481/bin
setenv ROOTBIN
# NNLOJET executable
#setenv NNLOJETBINPATH ${base}/${local}/src/NNLOJET_rev${revision}/driver
# PATH adaptation
if ( $?PATH ) then
#    setenv PATH ${MYCVMFS}/external/lhapdf/6.1.6/bin:${NNLOJETBINPATH}:${base}/${local}/bin:${ROOTBIN}:${PATH}
    setenv PATH ${base}/${local}/bin:${ROOTBIN}:${PATH}
else
#    setenv PATH ${MYCVMFS}/external/lhapdf/6.1.6/bin:${NNLOJETBINPATH}:${base}/${local}/bin:${ROOTBIN}
    setenv PATH ${base}/${local}/bin:${ROOTBIN}
endif
# Set additional (MY)CPPFLAGS for later use
setenv MYCPPFLAGS
if ( $?INCLUDE ) then
   if ( $?CPPFLAGS ) then
      setenv MYCPPFLAGS "${CPPFLAGS} -I${INCLUDE}"
   else
      setenv MYCPPFLAGS "-I${INCLUDE}"
   endif
else
   if ( $?CPPFLAGS ) then
      setenv MYCPPFLAGS "${CPPFLAGS}"
   endif
endif
echo "MYCPPFLAGS is $MYCPPFLAGS"
# Use modern gcc; gcc 4.4 from slc6 is too antique
if ( $?MYCVMFS ) then
   source ${MYCVMFS}/external/gcc/4.8.1/etc/profile.d/init.csh
endif
# LD_LIBRARY_PATH adaptation
if ( $?LD_LIBRARY_PATH ) then
#  setenv LD_LIBRARY_PATH ${MYCVMFS}/external/lhapdf/6.1.6/lib:${base}/${local}/lib:${base}/${local}/lib/root:${LD_LIBRARY_PATH}
  setenv LD_LIBRARY_PATH ${base}/${local}/lib:${base}/${local}/lib/root:${LD_LIBRARY_PATH}
else
#  setenv LD_LIBRARY_PATH ${MYCVMFS}/external/lhapdf/6.1.6/lib:${base}/${local}/lib:${base}/${local}/lib/root
  setenv LD_LIBRARY_PATH ${base}/${local}/lib:${base}/${local}/lib/root
endif



#==============================================================================
# Mandatory and sincerely recommended packages
#==============================================================================
#
# fastjet (any version >= 3 should work):
# Use the "--enable-allplugins" options to enable use of e.g. older Tevatron jet algorithms.
#------------------------------------------------------------------------------
#set local="local"
set arc="fastjet-3.3.0"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --enable-shared --enable-allplugins --prefix=${base}/${local} --bindir=${base}/${local}/bin
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
# LHAPDF (>= 6.2.0 is recommended):
# LHAPDF v5 might still work, but support for it has been stopped!
#------------------------------------------------------------------------------
# set arc="lhapdf-5.9.1"
set arc="LHAPDF-6.2.0"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=${base}/${local} CPPFLAGS="${MYCPPFLAGS}"
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
# ==> This enables the option --with-root of the fastNLO_toolkit to produce
#     ROOT histograms from the calculated cross sections and uncertainties
#
set arc="root-5.34.25"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  mv root ${arc}
  cd ${arc}
  ./configure --prefix=${base}/${local} --etcdir=${base}/${local}/etc --enable-python --enable-minuit2
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
#
# HepMC, YODA, and Rivet:
#------------------------------------------------------------------------------
# ==> These enable the option --with-yoda of the fastNLO_toolkit to produce
#     YODA formatted e.g. NLO predictions for direct comparison to data with Rivet
#
# HepMC (should be < 3.0.0, e.g. v2.06.09; needed by Rivet):
# Version 3.0.0 uses cmake and gives errors while compiling with ROOT v5.34.25!
#------------------------------------------------------------------------------
set arc="HepMC-2.06.09"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=${base}/${local} --with-momentum=GEV --with-length=MM
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
# YODA (>= 1.6.7 is recommended; needed by Rivet):
# Python is enabled by default. If Python with ROOT interfacing is desired, use "--enable-root".
#------------------------------------------------------------------------------
set arc="YODA-1.6.7"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=${base}/${local} --enable-root CPPFLAGS="${MYCPPFLAGS}"
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
# Rivet (>= 2.5.4 is recommended; v1 is too old):
# Python is enabled by default.
#------------------------------------------------------------------------------
set arc="Rivet-2.5.4"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=${base}/${local} CPPFLAGS="${MYCPPFLAGS}"
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
#
# HOPPET (e.g. v1.1.5):
#------------------------------------------------------------------------------
# ==> This enables the option --with-hoppet of the fastNLO_toolkit to use
#     the alpha_s evolutions within HOPPET
#     fastNLO comes already with alpha_s evolutions from GRV or CRunDec, or
#     uses the one from the PDFs in LHAPDF
#
set arc="hoppet-1.1.5"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=${base}/${local}
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
#
# QCDNUM (Use the newer autotools-enabled version 17.01.12):
# Older QCDNUM versions including 17.00.06 have an error in the time-like evolution
# of the singlet fragmentation function at NLO, see arXiv:1602.08383.
# Non-autotools-enabled versions are not supported anymore.
#------------------------------------------------------------------------------
# ==> This enables the option --with-qcdnum of the fastNLO_toolkit to use
#     the alpha_s evolutions within QCDNUM
#
set arc="qcdnum-17-01-12"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=${base}/${local}
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
#
#==============================================================================
# The fastNLO Toolkit
#==============================================================================
#
# fastNLO Toolkit:
#------------------------------------------------------------------------------
set arc="fastnlo_toolkit-2.3.1pre-2411"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
#  ./configure --prefix=${base}/${local} --enable-pyext
  ./configure --prefix=${base}/${local} --with-yoda --with-hoppet --with-qcdnum --with-root --enable-pyext
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
  ./configure --prefix=${base}/${local}
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
# fastNLO Interface NLOJet++:
#------------------------------------------------------------------------------
set arc="fastnlo_interface_nlojet-2.3.1pre-2411"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=${base}/${local}
  make -j4 install
  cd ..
  touch ${arc}_installed
endif
#
#
#==============================================================================
# fastNLO use with Sherpa & MCgrid (TO BE REENABLED)
#==============================================================================
#
# QD library required by BlackHat:
#------------------------------------------------------------------------------
#set arc="qd-2.3.17"
#if ( ! -e ${arc}_installed  ) then
#  tar xzf ${arc}.tar.gz
#  cd ${arc}
#  ./configure --prefix=${base}/${local} --enable-shared CXX=g++ CC=gcc FC=gfortran
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
#  ./configure --prefix=${base}/${local}
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
#  ./configure --prefix=${base}/${local}
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
#  ./configure --prefix=${base}/${local} --with-sqlite3=install --enable-gzip --enable-lhole ${WITH_MPI} --enable-fastjet=${base}/${local} --enable-lhapdf=${base}/${local} --enable-hepmc2=${base}/${local} --enable-rivet=${base}/${local} --enable-root=${base}/${local} --enable-blackhat=${base}/${local} --enable-openloops=${base}/${local}
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
#  ./configure --prefix=${base}/${local}
#  make -j4 install
## CXXFLAGS+=-fpermissive
#  cd ..
#  touch ${arc}_installed
#endif
#
# For usage, make sure to set PATH, PYTHON, LHAPDF, RIVET, PKG_CONFIG_PATH, ROOT environment variables appropriately!!!
#
# setenv PKG_CONFIG_PATH ${base}/${local}/lib/pkgconfig
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
#./configure --prefix=${base}/${local} --with-fastjet=${base}/${local} --with-hepmc=${base}/${local} --with-lhapdf=${base}/${local} --with-rivet=${base}/${local}
#make -j8
#make check
#make install
#
# NJet-matchbox:
#---------------
#./configure --prefix=${base}/${local} --enable-5jet
#make -j8
#make check
#make install
#
# Herwig++-matchbox:
#-------------------
#install gengetopt from distro
#autoreconf -i
#./configure --prefix=${base}/${local} --with-nlojet=${base}/${local} --with-njet=${base}/${local} --with-fastjet=${base}/${local}
#./configure --prefix=${base}/${local} --with-fastjet=${base}/${local} --with-madgraph=${base}/${local}/MadGraph-matchbox --with-njet=${base}/${local}
#make -j8
#make check
#make install
#make  check-local
#
#
#==============================================================================
# DEPRECATED
#==============================================================================
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
#   ./configure --prefix=${base}/${local} --enable-pyext --with-hoppet --with-qcdnum
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
#   ./configure --prefix=${base}/${local}
#   make -j4 install
#   cd ..
#   touch ${arc}_installed
# endif
#
