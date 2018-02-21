#!/bin/csh -ef
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
#    MPICH(2), the GNU scientific library GSL, or SWIG for Python interfaces to
#    the C++ libraries.
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
#           adequate software versions e.g. via "module load module/name/version".
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
#        that all 'local' packages of this step are unpacked and built in
#        ${base}/${local}/src
#        and that the common install directory will be
#        ${base}/${local}
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
# 0. Three commandline arguments are supported:
#     - The usage is: $0 ${base} ${local} [path to cvmfs software]
#     - 1st argument (mandatory): Base dir for installation
#     - 2nd argument (mandatory): Subdir in base dir for this installation
#     - 3rd argument [optional]: Path to newer gcc in cvmfs software distribution
#
#       In case of CMS e.g. "/cvmfs/cms.cern.ch/slc6_amd64_gcc481" can be used with
#          source ${MYCVMFS}/external/gcc/4.8.1/etc/profile.d/init.csh
#       later on.
#       The LCG versions like in "/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8.4"
#       have not been tested. They would require e.g.
#          source ${MYCVMFS}/x86_64-cc7/setup.csh
#       and do set the FC, CXX, and CC environment variables, which is not advisable,
#       because these variables should be left for the end user to modify as is done
#       in some command lines here. Beware, since this might lead to conflicts!
#
# 1. The install directory assumed here is ${base}/${local} with binaries to be
#    installed in ${base}/${local}/bin. Make sure that right from the start
#    - ${base}/${local}/bin is in the search path to find commands like "toolname"-config,
#    - the right version of "toolname"-config is used, either the desired system one or
#      the one from ${base}/${local}/bin.
#    Usually, this is set up via defining the PATH environment variable appropriately.
#
# 2. If you want to install BlackHat 0.9.9 for use within Sherpa-->MCGrid-->fastNLO:
#    - Compilation with gcc-4.8.x or newer gives compile errors -->
#      you must use the patched version blackhat-0.9.9-patched.tar.gz you should have
#      received together with this script.
#    - Install the libssl-dev (Ubuntu) to resolve potentially missing dependency
#    - The qd-2.3.17 package required by BlackHat does not necessarily choose your
#      system default compiler. To avoid e.g. using the Intel compiler configure
#      this package with explicit options "CXX=g++ CC=gcc FC=gfortran".
#      If you DO use the Intel compiler, then the file qd_real.cpp has to be
#      made with "make CXXFLAGS=-O1" to avoid infinite loops.
#
# 3. If you want to install OpenLoops 1.3.1 for use within Sherpa-->MCGrid-->fastNLO:
#    - The use of OpenLoops is currently tested and the package is installed.
#      By default only the ppjj process is downloaded, see below.
#      ATTENTION: You need to be online to download and install processes!
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
#    via "module load mpi/openmpi/version". In that case set below withmpi to "1".
#    Please note that by default CXX support is not available in standard
#    installations of OpenMPI version 2. OpenMPI v2 is usable only if compiled
#    with --enable-mpi-cxx!
#
#==============================================================================
# Check command line arguments
#==============================================================================
if ( $#argv < 2 ) then
   echo ""
   echo "#=============================================================================="
   echo "Usage: $0 basedir subdir [path to cvmfs software]"
   echo "1st argument: Base dir for installations, absolute path to e.g. $HOME"
   echo "2nd argument: Subdir in base dir for this installation, relative path to e.g. local"
   echo "3rd optional argument: Path to newer gcc in cvmfs software distribution,"
   echo "                       e.g. /cvmfs/cms.cern.ch/slc6_amd64_gcc481"
   echo "                       The minimally required version is gcc 4.8.1!"
   echo "                       If set, please check also that LHAPDF can be found in:"
   echo '                       $MYCVMFS/external/lhapdf'
   echo "                       Set to _ to skip this setting."
   echo "4th optional argument: Include grid creation with Sherpa+MCgrid? def.=0"
   echo "5th optional argument: Include grid creation with NNLOJET? def.=0"
   echo "6th optional argument: Include optional packages for grid evaluation? def.=0"
   echo "7th optional argument: Include optional python extensions to packages? def.=0"
   echo "8th optional argument: Include ROOT 5 (gcc<v5) or 6? def.=0, alt.=[5,6]"
   echo "9th optional argument: No. of cores to be used, def.=8"
   echo "#=============================================================================="
   echo ""
   exit 0
endif
#==============================================================================
# Set installation location
#==============================================================================
set base=$1
set local=$2
set srcdir=$cwd
#==============================================================================
# Define scope of installation
#==============================================================================
# Use CVMFS software repository for newer gcc and for LHAPDF if necessary
if ( $#argv > 2 && $3 != "_" ) then
   setenv MYCVMFS $3
# Do NOT use 6.2.1. This CVMFS installation is SHIT and does not find the PDF sets!
#   set lhapdfbasepath=${MYCVMFS}/external/lhapdf/6.2.1
   set lhapdfbasepath=${MYCVMFS}/external/lhapdf/6.1.6
#   echo "MYCVMFS is set to: $MYCVMFS"
else
   unsetenv MYCVMFS
   set lhapdfbasepath=${base}/${local}
#   echo "MYCVMFS is not set"
endif
# With interface to Sherpa via MCGrid?
set withsherpa=0
if ( $#argv > 3 ) then
    set withsherpa=$4
endif
# With interface to NNLOJET? Attention! NNLOJET is not yet publically available!
set withnnlojet=0
if ( $#argv > 4 ) then
    set withnnlojet=$5
# Previous: set revision=3738
    set revision=4585
endif
# With optional packages for grid evaluation?
set withoptional=0
if ( $#argv > 5 ) then
    set withoptional=$6
endif
# With optional Python extensions? On some systems compile errors occur!
# Note: Python is quite useful for evaluating or preparing results, but
#       not required for mass production on compute clusters.
# BUT: BlackHat needs Python!
set withpython=0
if ( $#argv > 6 ) then
    set withpython=$7
endif
# Install optional ROOT extensions? On some systems compile errors occur!
# Note: ROOT is quite useful for evaluating or preparing results, but
#       not required for mass production on compute clusters.
# BUT:  APPLgrid requires ROOT!
set withroot=0
if ( $#argv > 7 ) then
    set withroot=$8
endif

# Default is with OpenMPI if either NNLOJET or Sherpa are requested
set withmpi=0
set withmpiinstall=0
if ( $withsherpa || $withnnlojet ) then
    set withmpi=1
    # Default is with OpenMPI installation to ensure proper configuration
    set withmpiinstall=1 # def.=1
    # If not in the search path, set MPI_HOME to specify where i.a. bin/mpicxx is found.
    set mpihome=""
    # set mpihome=${base}/${local} # If not already in search path as recommended
    # If MPI available, decide whether to run NNLOJET in multithreaded mode or not
    set mpinnlo=$withmpi # def.=$withmpi
    # set mpinnlo=0 # quattro, zen, ekpcms006
endif

# Number of cores to be used
set cores=8
if ( $#argv > 8 ) then
    set cores=$9
endif

# With interface to HERWIG7? Not yet implemented!
set withherwig=0
#==============================================================================
# Check and set defaults for some specific environment variables
#==============================================================================
echo "#=============================================================================="
echo "# Preparing installation:"
echo "#=============================================================================="
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
# Use Python options/extensions?
set pyextopt="--disable-pyext"
set pythonopt="--disable-python"
if ( $withpython ) then
    set pyextopt="--enable-pyext"
    set pythonopt="--enable-python"
endif
# Use Root options/extensions?
set rootopt="--without-root"
set rootenable="--disable-root"
set rootenablepath="--disable-root"
if ( $withroot > 1 ) then
    set rootopt="--with-root"
    set rootenable="--enable-root"
    set rootenablepath="--enable-root=${base}/${local}"
endif
# Use Sherpa with MPI
# Attention! Buggy configure.ac in Sherpa:
#   Do not specify any path with --enable-mpi=..., because that would switch off MPI.
#   Instead, one MUST specify, where the executables for compilation are; CXX= etc.
if ( $withmpi ) then
  if ( $?MPI_HOME ) then
    set mpiopt="--enable-mpi CC=${MPI_HOME}/bin/mpicc CXX=${MPI_HOME}/bin/mpicxx MPICXX=${MPI_HOME}/bin/mpicxx FC=${MPI_HOME}/bin/mpifort"
  else if ( $mpihome != "" ) then
    set mpiopt="--enable-mpi CC=${mpihome}/bin/mpicc CXX=${mpihome}/bin/mpicxx MPICXX=${mpihome}/bin/mpicxx FC=${mpihome}/bin/mpifort"
  else # try
    set mpiopt="--enable-mpi CC=mpicc CXX=mpicxx MPICXX=mpicxx FC=mpifort"
  endif
else
  set mpiopt=""
endif
#
# PATH adaptation
#
echo '#\!/bin/csh -f'  >! fnlosrc_source.csh
echo '#\!/bin/bash -f' >! fnlosrc_source.sh
if ( $?PATH ) then
   setenv PATH ${base}/${local}/bin:${PATH}
   echo 'setenv PATH '"${base}/${local}/bin:"'${PATH}' >> fnlosrc_source.csh
   echo 'export PATH='"${base}/${local}/bin:"'${PATH}' >> fnlosrc_source.sh
else
   setenv PATH ${base}/${local}/bin
   echo 'setenv PATH '"${base}/${local}/bin" >> fnlosrc_source.csh
   echo 'export PATH='"${base}/${local}/bin" >> fnlosrc_source.sh
endif
# $PATH set from now on ...
#
# If $MYCVMFS is set, use the LHAPDF installation and PDF sets from CVMFS.
if ( $?MYCVMFS ) then
   setenv PATH ${lhapdfbasepath}/bin:${PATH}
   echo 'setenv PATH '"${lhapdfbasepath}/bin:"'${PATH}' >> fnlosrc_source.csh
   echo 'export PATH='"${lhapdfbasepath}/bin:"'${PATH}' >> fnlosrc_source.sh
endif
# If $MYCVMFS is set, one could also use the ROOT installation from CVMFS, BUT
# because of frequent issues with preinstalled ROOT versions we will either install
# our own version below, if desired, or nothing at all, see command line options.
# In particular for producing fastNLO grids on compute clusters ROOT is not required at all!
# APPLgrid needs ROOT though, but is not installed currently.
#  (Uncomment to try using ROOT installation from CVMFS)
#if ( ($withroot > 1) && $?MYCVMFS ) then
#   setenv ROOTBINPATH ${MYCVMFS}/cms/cmssw/CMSSW_7_2_4/external/slc6_amd64_gcc481/bin
#   setenv PATH ${ROOTBINPATH}:${PATH}
#   echo 'setenv ROOTBINPATH '"${MYCVMFS}/cms/cmssw/CMSSW_7_2_4/external/slc6_amd64_gcc481/bin"'' >> fnlosrc_source.csh
#   echo 'export ROOTBINPATH='"${MYCVMFS}/cms/cmssw/CMSSW_7_2_4/external/slc6_amd64_gcc481/bin"'' >> fnlosrc_source.sh
#   echo 'setenv PATH '"${ROOTBINPATH}:"'${PATH}' >> fnlosrc_source.csh
#   echo 'export PATH='"${ROOTBINPATH}:"'${PATH}' >> fnlosrc_source.sh
#endif
# Just in case: NNLOJET executable
if ( $withnnlojet ) then
   setenv NNLOJET_BIN_PATH ${base}/${local}/src/NNLOJET_rev${revision}/driver
   setenv PATH ${NNLOJET_BIN_PATH}:${PATH}
   echo 'setenv NNLOJET_BIN_PATH '"${NNLOJET_BIN_PATH}" >> fnlosrc_source.csh
   echo 'export NNLOJET_BIN_PATH='"${NNLOJET_BIN_PATH}" >> fnlosrc_source.sh
   echo 'setenv PATH '"${NNLOJET_BIN_PATH}:"'${PATH}'   >> fnlosrc_source.csh
   echo 'export PATH='"${NNLOJET_BIN_PATH}:"'${PATH}'   >> fnlosrc_source.sh
endif
echo ""
echo "ATTENTION: PATH environment complemented!"
echo "   PATH has been set to:"
echo "   $PATH"
echo ""
#
# Set additional (MY)CPPFLAGS for later use
#
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
if ( $MYCPPFLAGS != "" ) then
    echo "MYCPPFLAGS is $MYCPPFLAGS"
endif
# Use modern gcc; gcc 4.4 from slc6 is too antique
if ( $?MYCVMFS ) then
   source ${MYCVMFS}/external/gcc/4.8.1/etc/profile.d/init.csh
   echo 'source '"${MYCVMFS}"/external/gcc/4.8.1/etc/profile.d/init.csh'' >> fnlosrc_source.csh
   echo 'source '"${MYCVMFS}"/external/gcc/4.8.1/etc/profile.d/init.sh'' >> fnlosrc_source.sh
   echo ""
   echo "ATTENTION: Environment has been complemented by source'ing"
   echo "   ${MYCVMFS}/external/gcc/4.8.1/etc/profile.d/init.[c]sh!"
   echo ""
endif
#
# LD_LIBRARY_PATH adaptation
#
if ( $?LD_LIBRARY_PATH ) then
  setenv LD_LIBRARY_PATH ${base}/${local}/lib:${LD_LIBRARY_PATH}
  echo 'setenv LD_LIBRARY_PATH '"${base}/${local}/lib:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.csh
  echo 'export LD_LIBRARY_PATH='"${base}/${local}/lib:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.sh
else
  setenv LD_LIBRARY_PATH ${base}/${local}/lib
  echo 'setenv LD_LIBRARY_PATH '"${base}/${local}/lib" >> fnlosrc_source.csh
  echo 'export LD_LIBRARY_PATH='"${base}/${local}/lib" >> fnlosrc_source.sh
endif
# $LD_LIBRARY_PATH set from now on ...
#
# If $MYCVMFS is set, use the LHAPDF installation and PDF sets from CVMFS.
if ( $?MYCVMFS ) then
  setenv LD_LIBRARY_PATH ${lhapdfbasepath}/lib:${LD_LIBRARY_PATH}
  echo 'setenv LD_LIBRARY_PATH '"${lhapdfbasepath}/lib:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.csh
  echo 'export LD_LIBRARY_PATH='"${lhapdfbasepath}/lib:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.sh
  setenv LHAPDF_DATA_PATH ${lhapdfbasepath}/share/LHAPDF
  echo 'setenv LHAPDF_DATA_PATH '"${LHAPDF_DATA_PATH}" >> fnlosrc_source.csh
  echo 'export LHAPDF_DATA_PATH='"${LHAPDF_DATA_PATH}" >> fnlosrc_source.sh
endif
# If $MYCVMFS is set, one could also use the ROOT installation from CVMFS, BUT
# because of frequent issues with preinstalled ROOT versions we will either install
# our own version below, if desired, or nothing at all, see command line options.
# In particular for producing fastNLO grids on compute clusters ROOT is not required at all!
# APPLgrid needs ROOT though, but is not installed currently.
#  (Uncomment to try using ROOT installation from CVMFS)
if ( $withroot > 1 ) then
#  if ( $?MYCVMFS ) then
#    setenv LD_LIBRARY_PATH ${MYCVMFS}/cms/cmssw/CMSSW_7_2_4/external/slc6_amd64_gcc481/lib/root:${LD_LIBRARY_PATH}
#    echo 'setenv LD_LIBRARY_PATH '"${MYCVMFS}/cms/cmssw/CMSSW_7_2_4/external/slc6_amd64_gcc481/lib/root:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.csh
#    echo 'export LD_LIBRARY_PATH='"${MYCVMFS}/cms/cmssw/CMSSW_7_2_4/external/slc6_amd64_gcc481/lib/root:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.sh
#  else
    setenv LD_LIBRARY_PATH ${base}/${local}/lib/root:${LD_LIBRARY_PATH}
    echo 'setenv LD_LIBRARY_PATH '"${base}/${local}/lib/root:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.csh
    echo 'export LD_LIBRARY_PATH='"${base}/${local}/lib/root:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.sh
#  endif
endif
echo ""
echo "ATTENTION: LD_LIBRARY_PATH environment complemented!"
echo "   LD_LIBRARY_PATH has been set to:"
echo "   $LD_LIBRARY_PATH"
echo ""
#
#
#==============================================================================
# Mandatory and sincerely recommended packages
#==============================================================================
#
# fastjet (any version >= 3 should work):
# Use the "--enable-allplugins" options to enable use of e.g. older Tevatron jet algorithms.
#------------------------------------------------------------------------------
set arc="fastjet-3.3.0"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --enable-shared --enable-allplugins --prefix=${base}/${local} --bindir=${base}/${local}/bin
  make -j${cores} install
  cd ..
  touch ${arc}_installed
endif
#
# LHAPDF (>= 6.2.0 is recommended):
#
#------------------------------------------------------------------------------
if ( ! $?MYCVMFS ) then
  set arc="LHAPDF-6.2.1"
  if ( ! -e ${arc}_installed  ) then
    tar xzf ${arc}.tar.gz
    cd ${arc}
    ./configure --prefix=${base}/${local} CPPFLAGS="${MYCPPFLAGS}" ${pythonopt}
    make -j${cores} install
    cd ..
# Download default PDF sets, only possible if installed with Python support!
    if ( $withpython ) then
      ${base}/${local}/bin/lhapdf install NNPDF31_nlo_as_0118 NNPDF31_nnlo_as_0118 CT14nlo CT14nnlo
    else
      echo ""
      echo "ATTENTION: LHAPDF has been installed without any PDF sets!"
      echo "           Either install with Python support or make such sets accessible by other means."
      echo ""
    endif
    touch ${arc}_installed
  endif
endif
#
# OpenMPI (version 2.1.2 has been tested to work, if it was configured with --enable-mpi-cxx):
#
#------------------------------------------------------------------------------
if ( $withmpiinstall ) then
    set arc="openmpi-2.1.2"
    if ( ! -e ${arc}_installed  ) then
    tar xzf ${arc}.tar.gz
    cd ${arc}
    ./configure --prefix=${base}/${local} --enable-mpi-cxx
    make -j${cores} install
    cd ..
    touch ${arc}_installed
    endif
endif
#
#
#==============================================================================
# Additional packages required for some of fastNLO optional features
# Skip any of these parts that is not desired.
#==============================================================================
#
# Don't try for now to use cvmfs installation of ROOT
#
# ROOT (e.g. v5.34.25; optional use by YODA):
# If Python support is desired, use "--enable-python".
#------------------------------------------------------------------------------
# ==> This enables the option --with-root of the fastNLO_toolkit to produce
#     ROOT histograms from the calculated cross sections and uncertainties
#
if ( $withroot > 1 ) then
    if ( $withroot == 5 ) then
#        set arc="root-5.34.25" # OK for gcc versions < 5
        set arc="root-5.34.26" # Patched thanks to A. Wisecarver
        if ( ! -e ${arc}_installed  ) then
            tar xzf ${arc}-patched.tar.gz # Needed for older versions 5
            mv root ${arc}
        endif
    else if ( $withroot == 6 ) then
        set arc="root-6.08.06" # OK for gcc >= 5. Needs CMake >= 3.4.3!
#        set arc="root-6.10.06" # Don't use this one! Buggy -lImt library dependence ...
        if ( ! -e ${arc}_installed  ) then
            tar xzf ${arc}.tar.gz
        endif
    else
        echo "ERROR! Unknown ROOT version selected: $withroot"
        exit(1)
    endif
    if ( ! -e ${arc}_installed  ) then
        cd ${arc}
        ./configure --prefix=${base}/${local} --etcdir=${base}/${local}/etc ${pythonopt} --enable-minuit2
        make -j${cores} install
        cd ..
        touch ${arc}_installed
    endif
  endif
endif

#
#
# HepMC, YODA, and Rivet:
#------------------------------------------------------------------------------
# ==> These enable the option --with-yoda of the fastNLO_toolkit to produce
#     YODA formatted e.g. NLO predictions for direct comparison to data with Rivet
#     Also these are required for use with Sherpa+MCgrid.
#
# HepMC (should be < 3.0.0, e.g. v2.06.09; needed by Rivet):
# Version 3.0.0 uses cmake and gives errors while compiling with ROOT v5.34.25!
#------------------------------------------------------------------------------
set arc="HepMC-2.06.09"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=${base}/${local} --with-momentum=GEV --with-length=MM
  make -j${cores} install
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
  ./configure --prefix=${base}/${local} ${pyextopt} ${rootenable} CPPFLAGS="${MYCPPFLAGS}"
  make -j${cores} install
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
  ./configure --prefix=${base}/${local} ${pyextopt} CPPFLAGS="${MYCPPFLAGS}"
  make -j${cores} install
  cd ..
  touch ${arc}_installed
endif
#
#
# HOPPET (the latest, v1.1.5, needs to be patched for a perl issue on newer systems):
# (Problem with perl script in systems with gcc >= 5! Patch needed: hoppet-1.1.5-patched.)
#------------------------------------------------------------------------------
# ==> This enables the option --with-hoppet of the fastNLO_toolkit to use
#     the alpha_s evolutions within HOPPET
#     fastNLO comes already with alpha_s evolutions from GRV or CRunDec, or
#     uses the one from the PDFs in LHAPDF
#
if ( $withoptional ) then
    set arc="hoppet-1.1.5"
    if ( ! -e ${arc}_installed  ) then
    tar xzf ${arc}-patched.tar.gz # Need patched version on newer systems like Ubuntu 16.04
    cd ${arc}
    ./configure --prefix=${base}/${local}
    make -j${cores} install
    cd ..
    touch ${arc}_installed
    endif
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
#set arc="qcdnum-17-01-12"
if ( $withoptional ) then
    set arc="qcdnum-17-01-14"
    if ( ! -e ${arc}_installed  ) then
    tar xzf ${arc}.tar.gz
    cd ${arc}
    ./configure --prefix=${base}/${local}
    make -j${cores} install
    cd ..
    touch ${arc}_installed
    endif
endif
#
#
#==============================================================================
# The fastNLO Toolkit
#==============================================================================
#
# fastNLO Toolkit:
#------------------------------------------------------------------------------
set arc="fastnlo_toolkit-2.3.1pre-2441"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
# options depending on previous choices: --with-yoda --with-hoppet --with-qcdnum --with-root
# option to use python interface to library: --enable-pyext
#  ./configure --prefix=${base}/${local} --enable-pyext
  if ( $withoptional ) then
    ./configure --prefix=${base}/${local} --with-yoda --with-hoppet --with-qcdnum ${rootopt} ${pyextopt}
  else
    ./configure --prefix=${base}/${local} --with-yoda ${rootopt} ${pyextopt}
  endif
  make -j${cores} install
  cd ..
  touch ${arc}_installed
endif
#
#
#==============================================================================
# fastNLO use with NLOJet++
#==============================================================================
#
# NLOJet++ (patched version from fastNLO web page is required: NLOJet-4.1.3-patched2):
#------------------------------------------------------------------------------
set arc="nlojet++-4.1.3"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}-patched2.tar.gz
  cd ${arc}
  ./configure --prefix=${base}/${local}
  make -j${cores} install
  cd ..
  touch ${arc}_installed
endif
#
# fastNLO Interface NLOJet++:
#------------------------------------------------------------------------------
set arc="fastnlo_interface_nlojet-2.3.1pre-2424"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=${base}/${local}
  make -j${cores} install
  cd ..
  touch ${arc}_installed
endif
#
#
#==============================================================================
# fastNLO use with NNLOJET, if available
#==============================================================================
if ( $withnnlojet ) then
#
# APPLgrid for use in nnlo-bridge to NNLOJET:
#------------------------------------------------------------------------------
# APPLgrid requires ROOT!
# Do not install for now. On Centos 7 the APPLgrid strategy of using
# 'gfortran -print-file-name=libgfortran.a' does not work because of
# buggy gfortran installation.
# Anyway, a newer version from our APPLgrid colleagues will be necessary!!
#    if ( $withroot > 1 ) then
#        set arc="applgrid-1.4.93-rev1594M"
#        if ( ! -e ${arc}_installed  ) then
#            tar xzf ${arc}.tar.gz
#            cd ${arc}
#            ./configure --prefix=${base}/${local}
# Attention: No concurrent compilation with -j here!
#            make install
#            cd ..
#            touch ${arc}_installed
#        endif
#    endif
#
# nnlo-bridge to NNLOJet:
#------------------------------------------------------------------------------
    set arc="nnlo-bridge-0.0.36"
# Previous buggy: set rev="rev1683M3"
    set rev="rev1683M4"
    if ( ! -e ${arc}_installed  ) then
    tar xzf ${arc}-${rev}.tar.gz
    cd ${arc}
    ./configure --prefix=${base}/${local}
    make -j${cores} install
    cd ..
    touch ${arc}_installed
    endif
#
# NNLOJET
#------------------------------------------------------------------------------
# Set for single-thread usage of NNLOJET
    if ( $withmpi ) then
        setenv OMP_STACKSIZE 999999
        setenv OMP_NUM_THREADS 1
        echo 'setenv OMP_STACKSIZE 999999' >> fnlosrc_source.csh
        echo 'setenv OMP_NUM_THREADS 1' >> fnlosrc_source.csh
        echo 'export OMP_STACKSIZE=999999' >> fnlosrc_source.sh
        echo 'export OMP_NUM_THREADS=1' >> fnlosrc_source.sh
        echo ""
        echo "ATTENTION: OpenMP environment set to default for NNLOJET!"
        echo "   OMP_STACKSIZE and OMP_NUM_THREADS have been set to:"
        echo "   $OMP_STACKSIZE and $OMP_NUM_THREADS"
        echo ""
    endif
    set arc="NNLOJET_rev"${revision}
    if ( ! -e ${arc}_installed  ) then
    tar xzf ${arc}.tar.gz
    cd ${arc}/driver
    make depend
    make -j${cores} fillgrid=1
    cd ../..
    touch ${arc}_installed
    endif
endif
#
#
#==============================================================================
# fastNLO use with Sherpa & MCgrid (testing)
#==============================================================================
if ( $withsherpa ) then
#
# QD library required by BlackHat:
#------------------------------------------------------------------------------
set arc="qd-2.3.17"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}.tar.gz
  cd ${arc}
  ./configure --prefix=${base}/${local} --enable-shared CXX=g++ CC=gcc FC=gfortran
  make -j${cores} install
  cd ..
  touch ${arc}_installed
endif
#
# BlackHat for use within Sherpa:
# (Compile error for gcc > 4.7! Patch needed: blackhat-0.9.9-patched.)
#------------------------------------------------------------------------------
set arc="blackhat-0.9.9"
if ( ! -e ${arc}_installed  ) then
  tar xzf ${arc}-patched.tar.gz
  cd ${arc}
  ./configure --prefix=${base}/${local}
## Make without -j${cores} ... might be too heavy on some machines
##  make install
# Switch off deprecated header warnings
# Permissive setting for newer gcc5 compiler
    make -j${cores} install CXXFLAGS="-Wno-deprecated -fpermissive"
  cd ..
  touch ${arc}_installed
endif
#
# NJet for use within Sherpa:
#------------------------------------------------------------------------------
    set arc="njet-2.0.0"
    if ( ! -e ${arc}_installed  ) then
        tar xzf ${arc}.tar.gz
        cd ${arc}
        ./configure --prefix=${base}/${local}
        make -j${cores} install
        cd ..
        touch ${arc}_installed
    endif
#
# OpenLoops for use within Sherpa:
# Needs gcc >= 4.6!
#------------------------------------------------------------------------------
    set arc="OpenLoops-1.3.1"
    if ( ! -e ${arc}_installed  ) then
        tar xzf ${arc}.tar.gz
        cd ${arc}
        ./scons
# Install process lib for pp --> jj at NLO
# One can add other processes later on
        echo "Trying to download and install the ppjj process for OpenLoops."
        echo "This must fail if you are offline!"
        echo "In this case comment out the following command of this script and do the libinstall later."
        ./openloops libinstall ppjj
        cd ..
        touch ${arc}_installed
    endif
# Environment variables for OpenLoops
# Properly installed these should in principal not be necessary ...
# Probably not correct; workaround to download further processes?
#    setenv OL_PREFIX ${cwd}/OpenLoops-1.3.1
#    setenv LD_LIBRARY_PATH ${cwd}/OpenLoops-1.3.1/proclib:${cwd}/OpenLoops-1.3.1/lib:${LD_LIBRARY_PATH}
#
# Sherpa for use with MCgrid & fastNLO:
# Version >= 2.2.0 is required!
# For use of MPI with Sherpa corresponding system packages are required (MPICH, OpenMPI).
# (Repacked version required for compatibility with C++11 standard enforcement.)
#------------------------------------------------------------------------------
    set arc="SHERPA-MC-2.2.4"
    if ( ! -e ${arc}_installed  ) then
        tar xzf ${arc}-repacked.tar.gz
        cd ${arc}
        ./configure --prefix=${base}/${local} --with-sqlite3=install --enable-gzip --enable-lhole --enable-fastjet=${base}/${local} --enable-lhapdf=${lhapdfbasepath} --enable-hepmc2=${base}/${local} --enable-rivet=${base}/${local} ${rootenablepath} --enable-openloops=${base}/${local}/src/OpenLoops-1.3.1 ${mpiopt} --enable-blackhat=${base}/${local}
        make -j${cores} install
        cd ..
        touch ${arc}_installed
    endif
#
# MCGrid for use with fastNLO:
# (Seems to be fixed: Requires CXXFLAGS=-fpermissive flag in make step with gcc-4.7.x)
# (Repacked version required for compatibility with C++11 standard enforcement.)
#------------------------------------------------------------------------------
    set arc="mcgrid-2.0.2"
    if ( ! -e ${arc}_installed  ) then
        tar xzf ${arc}-patched.tar.gz
        cd ${arc}
        ./configure --prefix=${base}/${local}
        make -j${cores} install # CXXFLAGS=-fpermissive
        cd ..
        touch ${arc}_installed
    endif
#
# Environment variables for Rivet with Sherpa/MCgrid/fastNLO
#------------------------------------------------------------------------------
    setenv RIVET_ANALYSIS_PATH ${base}/${local}/share/Rivet:${base}/${local}/share/fastnlo_interface_nlojet/RivetAdditions
    setenv RIVET_INFO_PATH ${base}/${local}/share/Rivet:${base}/${local}/share/fastnlo_interface_nlojet/RivetAdditions
    setenv RIVET_REF_PATH ${base}/${local}/share/Rivet:${base}/${local}/share/fastnlo_interface_nlojet/RivetAdditions
    setenv PKG_CONFIG_PATH ${base}/${local}/lib/pkgconfig
    echo 'setenv RIVET_ANALYSIS_PATH '"${base}/${local}/share/Rivet:${base}/${local}/share/fastnlo_interface_nlojet/RivetAdditions" >> fnlosrc_source.csh
    echo 'setenv RIVET_INFO_PATH '"${base}/${local}/share/Rivet:${base}/${local}/share/fastnlo_interface_nlojet/RivetAdditions" >> fnlosrc_source.csh
    echo 'setenv RIVET_REF_PATH '"${base}/${local}/share/Rivet:${base}/${local}/share/fastnlo_interface_nlojet/RivetAdditions" >> fnlosrc_source.csh
    echo 'setenv PKG_CONFIG_PATH '"${base}/${local}/lib/pkgconfig" >> fnlosrc_source.csh
    echo 'export RIVET_ANALYSIS_PATH='"${base}/${local}/share/Rivet:${base}/${local}/share/fastnlo_interface_nlojet/RivetAdditions" >> fnlosrc_source.sh
    echo 'export RIVET_INFO_PATH='"${base}/${local}/share/Rivet:${base}/${local}/share/fastnlo_interface_nlojet/RivetAdditions" >> fnlosrc_source.sh
    echo 'export RIVET_REF_PATH='"${base}/${local}/share/Rivet:${base}/${local}/share/fastnlo_interface_nlojet/RivetAdditions" >> fnlosrc_source.sh
    echo 'export PKG_CONFIG_PATH='"${base}/${local}/lib/pkgconfig" >> fnlosrc_source.sh
    echo ""
    echo "ATTENTION: Rivet environment has been adapted!"
    echo "   RIVET_ANALYSIS_PATH has been set to:"
    echo "   $RIVET_ANALYSIS_PATH"
    echo "   RIVET_INFO_PATH has been set to:"
    echo "   $RIVET_INFO_PATH"
    echo "   RIVET_REF_PATH has been set to:"
    echo "   $RIVET_REF_PATH"
    echo "   PKG_CONFIG_PATH has been set to:"
    echo "   $PKG_CONFIG_PATH"
    echo ""
#
# Examples for MCGrid use with fastNLO:
# Will be installed in ${base}/${local}/share/MCgrid2Examples-2.2.0
#------------------------------------------------------------------------------
    set arc="MCgrid2Examples-2.2.0"
    if ( ! -e ${arc}_installed  ) then
        tar xzf ${arc}.tar.gz -C ${base}/${local}/share
        if ( ! -d ${base}/${local}/share/${arc} ) then
           mv ${base}/${local}/share/examples ${base}/${local}/share/${arc}
        else
           echo "Warning: Could not move ${base}/${local}/share/examples!"
           echo "Directory ${base}/${local}/share/${arc} exists already."
           echo "Testing already present MCgrid plugin CMS_2011_S9086218."
        endif
        cd ${base}/${local}/share/${arc}/CMS_2011_S9086218
        make plugin-fastnlo install
        cd ${srcdir}
        touch ${arc}_installed
    endif
endif
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
#
# PYTHONPATH adaptation
#
if ( $withpython ) then
    setenv PYTHONPATHADD `find ${base}/${local}/lib* -name site-packages`
    if ( $?PYTHONPATH ) then
    setenv PYTHONPATH ${PYTHONPATHADD}:${PYTHONPATH}
    echo 'setenv PYTHONPATH '"${PYTHONPATHADD}:"'${PYTHONPATH}' >> fnlosrc_source.csh
    echo 'export PYTHONPATH='"${PYTHONPATHADD}:"'${PYTHONPATH}' >> fnlosrc_source.sh
    else
    setenv PYTHONPATH ${PYTHONPATHADD}
    echo 'setenv PYTHONPATH '"${PYTHONPATHADD}" >> fnlosrc_source.csh
    echo 'export PYTHONPATH='"${PYTHONPATHADD}" >> fnlosrc_source.sh
    endif
    echo ""
    echo "ATTENTION: PYTHONPATH environment complemented!"
    echo "   PYTHONPATH has been set to:"
    echo "   $PYTHONPATH"
    echo ""
endif
#
# Clean up unused variables and finish
#------------------------------------------------------------------------------
if ( $MYCPPFLAGS == "" ) then
    unsetenv MYCPPFLAGS
endif
if ( $withpython ) then
    if ( $PYTHONPATHADD == "" ) then
        unsetenv PYTHONPATHADD
    endif
endif
#
echo "ATTENTION: The required changes to your environment have been written to:"
echo "           fnlosrc_source.[c]sh"
echo "Please check that environment variables have been set appropriately,"
echo "e.g. PATH, LD_LIBRARY_PATH, ..."
echo "and 'source' the corresponding file before using this installation."
#
echo ""
echo "#=============================================================================="
echo "# All done! Have fun!"
echo "#=============================================================================="
echo ""
exit
#
# Please ignore!
# xfitter needs: libblasdevel, liblapackdevel, libyamldevel
#
