#!/bin/csh -ef
###############################################################################
# FULL INSTALL, KR 18.08.2017
#
# Depending on the purpose of an installation, up to two fastNLO packages are
# necessary plus additional nnlo-bridge and theory code packages.
#
# 1) The fastNLO Toolkit for the evaluation of any fastNLO table
# 2) The fastNLO Interface to NLOJet++ for the creation of LO or NLO
#    fastNLO tables with NLOJet++
# 3) The nnlo-bridge common interface with APPLgrid to NNLOJET
#    for the creation of LO, NLO, or NNLO APPLgrid/fastNLO grids
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
#    never needed Boost or ROOT ...
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
#        ${base}/src
#        and that the common install directory will be
#        ${base}
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
# 0. If no argument is given, the supported command line arguments are printed out,
#    see also below.
#
# 1. The install directory assumed here is ${base} with binaries to be
#    installed in ${base}/bin. Make sure that right from the start
#    - ${base}/bin is in the search path to find commands like "toolname"-config,
#    - the right version of "toolname"-config is used, either the desired system one or
#      the one from ${base}/bin.
#    Usually, this is set up via defining the PATH environment variable appropriately.
#
# 2. If the system compiler is not recent enough, is not C++11 compatible, or
#    otherwise incompatible with the selected packages, another GCC compiler
#    collection from /cvmfs can be used:
#       In case of CMS e.g. "/cvmfs/cms.cern.ch/slc6_amd64_gcc481" can be used with
#          source ${MYCVMFS}/external/gcc/4.8.1/etc/profile.d/init.csh
#
#       The LCG versions like in "/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8.4"
#       have not been tested. They would require e.g.
#          source ${MYCVMFS}/x86_64-cc7/setup.csh
#       and do set the FC, CXX, and CC environment variables, which is not advisable,
#       because these variables should be left for the end user to modify as is done
#       in some command lines here. Beware, since this might lead to conflicts!
#
# 3. If you want to install BlackHat 0.9.9 for use within Sherpa-->MCGrid-->fastNLO:
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
# 4. If you want to install OpenLoops 1.3.1 for use within Sherpa-->MCGrid-->fastNLO:
#    - The use of OpenLoops is currently tested and the package is installed.
#      By default only the ppjj process is downloaded, see below.
#      ATTENTION: You need to be online to download and install processes!
#
# 5. If you switch to a newer compiler via the "module" command,
#    "module load compiler/gnu/4.8"
#    and you still use Boost-dependent packages it might be necessary to
#    also activate a newer boost library (>= 1.48.0) used e.g. by older LHAPDF (< 6.2.0),
#    YODA (1.4.x < YODA < 1.6.x), or RIVET (2.3.x < Rivet < 2.5.x) packages:
#    "module load lib/boost/1.56.0"
#    In that case configure the packages by specifying additional include
#    paths with CPPFLAGS="-I${INCLUDE}" or more explicitly
#    by CPPFLAGS="-I/add/include/path1 -I/add/include/path2".
#
# 6. SHERPA can be used with MPI for parallel processing if desired. In that case,
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
   echo "Usage: $0 basedir [optional_argument_2] [optional_argument_3] ..."
   echo "  1st argument: Base dir for installations, absolute path to e.g. ${HOME}/local"
   echo "  2nd optional argument: Base path to additional software to be taken from cvmfs,"
   echo "                         e.g. /cvmfs/cms.cern.ch/slc6_amd64_gcc700"
   echo "                         or   /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8"
   echo "                         Set to _ to only use your local system and skip this setting."
   echo "  3rd optional argument: Sub path to GCC setup script in cvmfs software distribution."
   echo "                         The minimally required version is gcc 4.8.1(!)"
   echo "                         e.g. external/gcc/7.0.0-omkpbe2/etc/profile.d/init"
   echo "                         or   x86_64-slc6/setup"
   echo "                         (.csh or .sh are added automatically)."
   echo "                         Set to _ to use your system compiler and skip this setting."
   echo "  4th optional argument: Include LHAPDF6 from CVMFS? def.=0"
   echo "                         By default this script installs its own version of LHAPDF."
   echo "                         Give path to desired bin/lhapdf-config in cvmfs to try using another one,"
   echo "                         e.g. /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/lhapdf6/6.1.5"
   echo "  5th optional argument: Include grid creation with NLOJet++? def.=1"
   echo "  6th optional argument: Include grid creation with NNLOJET? def.=0"
   echo "  7th optional argument: Include grid creation with Sherpa+MCgrid? def.=0"
   echo "  8th optional argument: Include optional packages for grid evaluation? def.=0"
   echo "  9th optional argument: Include optional python extensions to packages? def.=0"
   echo " 10th optional argument: Include ROOT extensions to packages (def.=0, alt.=[5,6] or"
   echo "                         give path to bin/root-config of preinstalled ROOT)"
   echo "                         0: no ROOT, 5: try src install of ROOT 5 (gcc <= v5), 6: try src install of ROOT 6"
   echo "                         (requires cmake), path: try preinstalled ROOT by giving path to bin/root-config"
   echo "                         e.g. /cvmfs/sft.cern.ch/lcg/releases/ROOT/5.34.25-8ef6d/x86_64-slc6-gcc48-opt"
   echo " 11th optional argument: No. of cores to be used, def.=8"
   echo " 12th optional argument: Activate multithread integrations for NNLOJET standalone installation, def. = 0"
   echo "#=============================================================================="
   echo ""
   exit 0
endif

#==============================================================================
# Set installation location
#==============================================================================
echo "#=============================================================================="
echo "# Installation settings:"
echo "#=============================================================================="
set base=$1
set srcdir=$cwd
set tab = "`printf '\t'`"
echo "Installation location: $tab $1"
#==============================================================================
# Define scope of installation
#==============================================================================
# Use CVMFS software repository?
if ( $#argv > 1 && $2 != "_" ) then
   setenv MYCVMFS $2
   echo "MYCVMFS is set to: $tab $MYCVMFS"
else
   if ( $?MYCVMFS ) then
      unsetenv MYCVMFS
   endif
   echo "MYCVMFS is undefined: $tab $2"
endif
# gcc from CMVFS (e.g. external/gcc/4.8.1/etc/profile.d/init.[c]sh)?
set gccsetup="_"
if ( $#argv > 2 && $3 != "_" ) then
   set gccsetup=$3
endif
echo "gcc setup: $tab$tab $gccsetup"
# Use CVMFS software repository also for LHAPDF if requested
# Do NOT use 6.2.1. This CVMFS installation is buggy and does not find the PDF sets!
# Dare to try LHAPDF from cvmfs, e.g. external/lhapdf6/6.1.5/bin?
set lhapdfbasepath="${base}"
set lhapdfbinpath="${base}/bin"
#set lhapdfdatapath="${base}"
set withcvmfslhapdf=0
if ( $#argv > 3 && $4 != "_" ) then
   set withcvmfslhapdf=$4
endif
echo "LHAPDF from CVMFS: $tab $withcvmfslhapdf"
#
# Interfaces
#
# With interface to NLOJet++?
set withnlojetpp=1
if ( $#argv > 4 && $5 != "_" ) then
   set withnlojetpp=$5
endif
echo "NLOJet++ usage: $tab$tab $withnlojetpp"
# With interface to NNLOJET? Attention! NNLOJET is not yet publically available!
set withnnlojet=0
if ( $#argv > 5 && $6 != "_" ) then
   set withnnlojet=$6
# Previous: set revision=3738
# Used for fnl2332d:       set revision=4585
# Bug in Z+jet RV channel: set revision=4708
# Bug in DIS & pp jets VV: set revision=5088
   set revision=5918
endif
echo "NNLOJET usage: $tab$tab $withnnlojet"
# With interface to Sherpa via MCGrid?
set withsherpa=0
if ( $#argv > 6 && $7 != "_" ) then
   set withsherpa=$7
endif
echo "Sherpa usage: $tab$tab $withsherpa"
#
# Optional
#
# With optional packages for grid evaluation?
set withoptional=0
if ( $#argv > 7 && $8 != "_" ) then
   set withoptional=$8
endif
echo "Optional packages: $tab $withoptional"

# With optional Python extensions? On some systems compile errors occur!
# Note: Python is quite useful for evaluating or preparing results, but
#       not required for mass production on compute clusters.
# BUT: BlackHat needs Python!
set withpython=0
if ( $#argv > 8 && $9 != "_" ) then
   set withpython=$9
endif
echo "Python support: $tab $withpython"

# With optional ROOT extensions? On some systems compile errors occur!
# Note: ROOT is quite useful for evaluating or preparing results, but
#       not required for mass production on compute clusters.
# BUT:  APPLgrid requires ROOT!
# Dare to try ROOT from cvmfs, e.g. /cvmfs/sft.cern.ch/lcg/releases/ROOT/5.34.25-8ef6d/x86_64-slc6-gcc48-opt ?
# Be prepared for trouble when using ROOT!
set rootbasepath="${base}/root"
set rootbinpath="${base}/root/bin"
set withroot=0
set withextroot=0
if ( $#argv > 9 && $10 != "_" ) then
   set withroot=$10
   if ( $withroot != 0 && $withroot != 5 && $withroot != 6 ) then
      set withextroot=$withroot
   endif
endif
echo "ROOT support: $tab$tab $withroot"
# Do we have an extra path to ROOT?
if ( $withextroot != "0" ) then
   set rootbasepath="${withextroot}"
   set rootbinpath="${rootbasepath}/bin"
endif

# BUT LHAPDF version 6.1.6 still uses BOOST, uaargh!
#        set lhapdfbasepath=${MYCVMFS}/external/lhapdf/6.1.6
#        set lhapdfdatapath=${MYCVMFS}/external/lhapdf/6.1.6
# AND version 6.2.1 does not find its PDF sets, uaargh!
if ( $withcvmfslhapdf != "0" ) then
   set lhapdfbasepath="${withcvmfslhapdf}"
   set lhapdfbinpath="${lhapdfbasepath}/bin"
#   set lhapdfdatapath=${MYCVMFS}/external/lhapdf/6.2.1
endif

# Default is with OpenMPI if Sherpa is requested
set withmpi=0
set withmpiinstall=0
if ( $withsherpa ) then
   set withmpi=1
   # Default is with OpenMPI installation to ensure proper configuration
   set withmpiinstall=1 # def.=1
   # If not in the search path, set MPI_HOME to specify where i.a. bin/mpicxx is found.
   set mpihome=""
   # set mpihome=${base} # If not already in search path as recommended
endif

# Number of cores to be used
set cores=8
if ( $#argv > 10 && $11 != "_" ) then
   set cores=$11
endif
echo "No. of cores: $tab$tab $cores"

# Do multithread integrations with NNLOJET?
# NNLOJET uses OpenMP, which must be supported by the compiler. It does NOT use OpenMPI.
set mpnnlo=0
if ( $#argv > 11 && $12 != "_" ) then
   set mpnnlo=$12
endif
echo "Multithreaded NNLOJET: $tab $mpnnlo"

if ( $mpnnlo ) then
   if ( ! $withnnlojet ) then
      echo "ERROR! Please select NNLOJET installation!"
      exit 1
   endif
   set withnlojetpp=0
   set withsherpa=0
   set withoptional=0
   set withpython=0
   set withroot=0
   echo "ATTENTION: Only LHAPDF and NNLOJET in multithreaded mode will be installed!"
   echo "           NNLOJET multithread operations do not (yet) work together with fastNLO grid production."
endif

# With interface to HERWIG7? Not yet implemented!
set withherwig=0
#==============================================================================
# Check and set defaults for some specific environment variables
#==============================================================================
echo ""
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
set rootoptpath="--without-root"
set rootenable="--disable-root"
set rootenablepath="--disable-root"
if ( $withroot != "0" ) then
   set rootopt="--with-root"
   set rootoptpath="--with-root=${rootbasepath}"
   set rootenable="--enable-root"
   set rootenablepath="--enable-root=${rootbasepath}"
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
   setenv PATH ${base}/bin:${PATH}
   echo 'setenv PATH '"${base}/bin:"'${PATH}' >> fnlosrc_source.csh
   echo 'export PATH='"${base}/bin:"'${PATH}' >> fnlosrc_source.sh
else
   setenv PATH ${base}/bin
   echo 'setenv PATH '"${base}/bin" >> fnlosrc_source.csh
   echo 'export PATH='"${base}/bin" >> fnlosrc_source.sh
endif
# $PATH set from now on ...
#
# If $withroot equals 5 or 6, do local ROOT src installation
if ( $withroot == 5 || $withroot == 6 ) then
   setenv PATH ${rootbinpath}:${PATH}
   echo 'setenv PATH '"${rootbinpath}:"'${PATH}' >> fnlosrc_source.csh
   echo 'export PATH='"${rootbinpath}:"'${PATH}' >> fnlosrc_source.sh
endif
# If $withextroot is set, use the preinstalled ROOT.
if ( $withextroot != "0" ) then
   set rootoptpath="--with-root=${rootbasepath}"
   set rootenablepath="--enable-root=${rootbasepath}"
   setenv PATH ${rootbinpath}:${PATH}
   echo 'setenv PATH '"${rootbinpath}:"'${PATH}' >> fnlosrc_source.csh
   echo 'export PATH='"${rootbinpath}:"'${PATH}' >> fnlosrc_source.sh
endif
# If $withcvmfslhapdf is set, use the LHAPDF installation and PDF sets from CVMFS.
if ( $withcvmfslhapdf != "0" ) then
   setenv PATH ${lhapdfbinpath}:${PATH}
   echo 'setenv PATH '"${lhapdfbinpath}:"'${PATH}' >> fnlosrc_source.csh
   echo 'export PATH='"${lhapdfbinpath}:"'${PATH}' >> fnlosrc_source.sh
endif
# Just in case: NNLOJET executables
if ( $withnnlojet ) then
   setenv NNLOJET_BIN_PATH ${base}/src/NNLOJET_rev${revision}/driver
   setenv PATH ${NNLOJET_BIN_PATH}:${NNLOJET_BIN_PATH}/bin:${PATH}
   echo 'setenv NNLOJET_BIN_PATH '"${NNLOJET_BIN_PATH}" >> fnlosrc_source.csh
   echo 'export NNLOJET_BIN_PATH='"${NNLOJET_BIN_PATH}" >> fnlosrc_source.sh
   echo 'setenv PATH ${NNLOJET_BIN_PATH}:${NNLOJET_BIN_PATH}/bin:${PATH}' >> fnlosrc_source.csh
   echo 'export PATH=${NNLOJET_BIN_PATH}:${NNLOJET_BIN_PATH}/bin:${PATH}' >> fnlosrc_source.sh
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
   source ${MYCVMFS}/${gccsetup}.csh
   echo 'source '"${MYCVMFS}/${gccsetup}.csh"'' >> fnlosrc_source.csh
   echo 'source '"${MYCVMFS}/${gccsetup}.sh"''  >> fnlosrc_source.sh
   echo ""
   echo "ATTENTION: Environment has been complemented by source'ing"
   echo "   ${MYCVMFS}/${gccsetup}"
   echo ""
endif
#
# LD_LIBRARY_PATH adaptation
#
if ( $?LD_LIBRARY_PATH ) then
   setenv LD_LIBRARY_PATH ${base}/lib:${LD_LIBRARY_PATH}
   echo 'setenv LD_LIBRARY_PATH '"${base}/lib:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.csh
   echo 'export LD_LIBRARY_PATH='"${base}/lib:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.sh
else
   setenv LD_LIBRARY_PATH ${base}/lib
   echo 'setenv LD_LIBRARY_PATH '"${base}/lib" >> fnlosrc_source.csh
   echo 'export LD_LIBRARY_PATH='"${base}/lib" >> fnlosrc_source.sh
endif
# $LD_LIBRARY_PATH set from now on ...
#
if ( $withroot == 5 || $withroot == 6 ) then
   setenv LD_LIBRARY_PATH ${rootbasepath}/lib:${LD_LIBRARY_PATH}
   echo 'setenv LD_LIBRARY_PATH '"${rootbasepath}/lib:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.csh
   echo 'export LD_LIBRARY_PATH='"${rootbasepath}/lib:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.sh
endif
# If $withextroot is set, use the preinstalled ROOT.
if ( $withextroot != "0" ) then
    setenv LD_LIBRARY_PATH ${rootbasepath}/lib:${LD_LIBRARY_PATH}
    echo 'setenv LD_LIBRARY_PATH '"${rootbasepath}/lib:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.csh
    echo 'export LD_LIBRARY_PATH='"${rootbasepath}/lib:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.sh
endif
# If $withcvmfslhapdf is set, use the LHAPDF installation and PDF sets from CVMFS.
if ( $withcvmfslhapdf != "0" ) then
   setenv LD_LIBRARY_PATH ${lhapdfbasepath}/lib:${LD_LIBRARY_PATH}
   echo 'setenv LD_LIBRARY_PATH '"${lhapdfbasepath}/lib:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.csh
   echo 'export LD_LIBRARY_PATH='"${lhapdfbasepath}/lib:"'${LD_LIBRARY_PATH}' >> fnlosrc_source.sh
#   setenv LHAPDF_DATA_PATH ${lhapdfdatapath}/share/LHAPDF
#   echo 'setenv LHAPDF_DATA_PATH '"${LHAPDF_DATA_PATH}" >> fnlosrc_source.csh
#   echo 'export LHAPDF_DATA_PATH='"${LHAPDF_DATA_PATH}" >> fnlosrc_source.sh
endif
echo ""
echo "ATTENTION: LD_LIBRARY_PATH environment complemented!"
echo "   LD_LIBRARY_PATH has been set to:"
echo "   $LD_LIBRARY_PATH"
echo ""
#
# PYTHON_PATH backup
#
if ( $?PYTHONPATH ) then
   setenv PYTHONPATHORIG ${PYTHONPATH}
endif
#
#
#==============================================================================
# Mandatory and/or sincerely recommended packages
#==============================================================================
#
# fastjet (any version >= 3 should work):
# Use the "--enable-allplugins" options to enable use of e.g. older Tevatron jet algorithms.
#------------------------------------------------------------------------------
if ( ($withnlojetpp || $withsherpa) && ! $mpnnlo ) then
#   set arc="fastjet-3.3.0" # Last one ok; test 3.3.4
   set arc="fastjet-3.3.4"
   if ( ! -e ${arc}_installed  ) then
      tar xzf ${arc}.tar.gz
      cd ${arc}
      ./configure --enable-shared --enable-allplugins --prefix=${base} --bindir=${base}/bin
      make -j${cores} install
      cd ..
      touch ${arc}_installed
   endif
endif
#
# LHAPDF (>= 6.2.0 is recommended):
#
#------------------------------------------------------------------------------
if ( $withcvmfslhapdf == "0" ) then
#   set arc="LHAPDF-6.2.1"# Test 6.3.0
   set arc="LHAPDF-6.3.0"
   if ( ! -e ${arc}_installed  ) then
      tar xzf ${arc}.tar.gz
      cd ${arc}
      ./configure --prefix=${base} CPPFLAGS="${MYCPPFLAGS}" ${pythonopt}
      make -j${cores} install
      cd ..
# Download default PDF sets, only possible if installed with Python support!
      if ( $withpython ) then
# Temporarily complement PYTHONPATH here
         setenv PYTHONPATHADD `find ${base}/lib* -name site-packages | tr '[:space:]' ':'`
         if ( $?PYTHONPATHORIG ) then
            setenv PYTHONPATH ${PYTHONPATHADD}:${PYTHONPATHORIG}
         else
            setenv PYTHONPATH ${PYTHONPATHADD}
         endif
         ${base}/bin/lhapdf install NNPDF23_nnlo_as_0118 NNPDF31_nlo_as_0118 NNPDF31_nnlo_as_0118 CT14nlo CT14nnlo
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
      ./configure --prefix=${base} --enable-mpi-cxx
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
# ROOT (e.g. v5.34.25; optional use by YODA):
# If Python support is desired, use "--enable-python".
#------------------------------------------------------------------------------
# ==> This enables the option --with-root of the fastNLO_toolkit to produce
#     ROOT histograms from the calculated cross sections and uncertainties
#
if ( $withroot != 0 && $withroot != 5 && $withroot != 6 ) then
   echo ""
   echo "ATTENTION: Trying ROOT preinstalled to ${rootbasepath} for optional ROOT extensions!"
   echo "           For full ROOT usage you might want to source the thisroot.(csh) environment."
   echo ""
else if ( $withroot != 0 ) then
   if ( $withroot == 5 ) then
      set arc="root-5.34.25" # OK for gcc versions < 5
#      set arc="root-5.34.26" # Patched thanks to A. Wisecarver
#      set arc="root-5.34.38" # Gives errors
      if ( ! -e ${arc}_installed  ) then
#         tar xzf ${arc}-patched.tar.gz # Needed for older versions 5
         tar xzf ${arc}.tar.gz # Needed for older versions 5
         mv root ${arc}
      endif
   else if ( $withroot == 6 ) then
#      set arc="root-6.08.06" # OK for gcc >= 5. Needs CMake >= 3.4.3!
      set arc="root-6.14.06" # Not sure that it works!
#      set arc="root-6.20.00" # Needs CMake >= 3.9
#      set arc="root-6.22.02" # Needs CMake >= 3.9
      if ( ! -e ${arc}_installed  ) then
         tar xzf ${arc}.tar.gz
      endif
   endif
   if ( ! -e ${arc}_installed ) then
      cd ${arc}
      if ( $withroot == 5 ) then
         ./configure --prefix=${rootbasepath} --etcdir=${rootbasepath}/etc ${pythonopt} --enable-minuit2 --disable-xrootd
         make -j${cores} install
         cd ..
         touch ${arc}_installed
      else if ( $withroot == 6 ) then
         mkdir -p mybuild
         cd mybuild
         cmake .. -DCMAKE_INSTALL_PREFIX=${rootbasepath} -Dgnuinstall=ON -Ddavix=OFF # davix tries to download nonexisting stuff, very bad!
         make -j${cores} install
#         cmake ..
#         cmake --build . -- -j${cores}
         # As usual ROOT is buggy: CMAKE_INSTALL_BINDIR is not respected and everything ends in CMAKE_INSTALL_PREFIX!
         #      cmake -DCMAKE_INSTALL_PREFIX=${base}/root -DCMAKE_INSTALL_BINDIR=${base}/bin -DCMAKE_INSTALL_DATAROOTDIR=${base}/share -P cmake_install.cmake
         # Bug adapted version using addition of ${rootbasepath}/bin to PATH
#         cmake -DCMAKE_INSTALL_PREFIX=${rootbasepath} -P cmake_install.cmake
         cd ..
         cd ..
         touch ${arc}_installed
      else
         cd ..
      endif
   endif
endif

#
#
# YODA, HepMC, fastjet contrib, and Rivet:
#------------------------------------------------------------------------------
# ==> These enable the option --with-yoda of the fastNLO_toolkit to produce
#     YODA formatted e.g. NLO predictions for direct comparison to data with Rivet
#     Also these are required for use with Sherpa+MCgrid.
#
# HepMC (should be < 3.0.0, e.g. v2.06.09; needed by Rivet-2.x.x):
# Version 3.0.0 uses cmake and gives errors while compiling with ROOT v5.34.25!
#------------------------------------------------------------------------------
if ( $withoptional ) then
#   set arc="HepMC-2.06.09"
#   if ( ! -e ${arc}_installed  ) then
#      tar xzf ${arc}.tar.gz
#      cd ${arc}
#      ./configure --prefix=${base} --with-momentum=GEV --with-length=MM
#      make -j${cores} install
#      cd ..
#      touch ${arc}_installed
#   endif
#
# HepMC3 (needed by Rivet-3.x.x):
# Python is enabled by default. Requires cmake version 3.
#------------------------------------------------------------------------------
   set arc="HepMC3-3.2.2"
   if ( ! -e ${arc}_installed  ) then
      tar xzf ${arc}.tar.gz
      cd ${arc}
      mkdir -p hepmc3-build
      cd hepmc3-build
# Try to adapt the following line when Python support is desired AND cmake version is >= 3.7 (UNDOCUMENTED requirement, grrr)
#      cmake -DHEPMC3_ENABLE_ROOTIO=OFF -DCMAKE_INSTALL_PREFIX=${base} -DHEPMC3_Python_SITEARCH27=${base}/lib/python2.7/site-packages -DHEPMC3_Python_SITEARCH36=${base}/lib/python3.6/site-packages -DHepMC3_DIR=${base} ..
# HepMC3: Python bindings request a cmake policy of version 3.7 --> switch off Python support by default.
      cmake -DCMAKE_INSTALL_PREFIX=${base} -DHepMC3_DIR=${base} -DHEPMC3_ENABLE_PYTHON=OFF -DHEPMC3_ENABLE_ROOTIO=OFF ..
      make -j${cores} install
      cd ../..
      touch ${arc}_installed
   endif
#
# fastjet contrib (fjcontrib; needed by Rivet-3.x.x):
#
#------------------------------------------------------------------------------
   set arc="fjcontrib-1.045"
   if ( ! -e ${arc}_installed  ) then
      tar xzf ${arc}.tar.gz
      cd ${arc}
      ./configure --prefix=${base} --fastjet-config=${base}/bin/fastjet-config
      make -j${cores} install
      make -j${cores} fragile-shared-install # From fjcontrib: Dirty hack to provide a shared library
      cd ..
      touch ${arc}_installed
   endif
#
# YODA (>= 1.6.7 is recommended; needed by Rivet):
# Python is enabled by default. If Python with ROOT interfacing is desired, use "--enable-root".
#------------------------------------------------------------------------------
#   set arc="YODA-1.6.7"
   set arc="YODA-1.8.3"
   if ( ! -e ${arc}_installed  ) then
      tar xzf ${arc}.tar.gz
      cd ${arc}
      ./configure --prefix=${base} ${pyextopt} ${rootenable} CPPFLAGS="${MYCPPFLAGS}"
      make -j${cores} install
      cd ..
      touch ${arc}_installed
   endif
#
# Rivet (>= 2.5.4 is recommended; v1 is too old):
# Python is enabled by default.
#------------------------------------------------------------------------------
#   set arc="Rivet-2.5.4"
   set arc="Rivet-3.1.4" # Needs HepMC3 and fastjet contrib
   if ( ! -e ${arc}_installed  ) then
      tar xzf ${arc}.tar.gz
      cd ${arc}
#      ./configure --prefix=${base} ${pyextopt} CPPFLAGS="${MYCPPFLAGS}"
# ATTENTION: Remove analysis giving compile errors with gcc5
      rm -f analyses/pluginATLAS/*1790439*
      rehash
      ./configure --prefix=${base} ${pyextopt} --with-hepmc3=`HepMC3-config --prefix` CPPFLAGS="${MYCPPFLAGS}"
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
#   set arc="hoppet-1.1.5"
   set arc="hoppet-1.2.0"
   if ( ! -e ${arc}_installed  ) then
#      tar xzf ${arc}-patched.tar.gz # Need patched version on newer systems like Ubuntu 16.04
      tar xzf ${arc}.tar.gz
      cd ${arc}
      ./configure --prefix=${base}
      make -j${cores} install
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
#   set arc="qcdnum-17-01-12"
   set arc="qcdnum-17-01-14"
#   set arc="qcdnum-17-01-15" # Provokes error in fastNLO Toolkit check
   if ( ! -e ${arc}_installed  ) then
      tar xzf ${arc}.tar.gz
      cd ${arc}
      ./configure --prefix=${base}
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
if ( ! $mpnnlo ) then
#   set arc="fastnlo_toolkit-2.3.1-2585"
#   set arc="fastnlo_toolkit-2.3.1-2657"
#   set arc="fastnlo_toolkit-2.3.1-2753"
#   set arc="fastnlo_toolkit-2.3.1-2771"
   set arc="fastnlo_toolkit-2.5.0-2826"
   set rev=""
   if ( ! -e ${arc}_installed  ) then
      tar xzf ${arc}${rev}.tar.gz
      cd ${arc}
# options depending on previous choices: --with-yoda --with-hoppet --with-qcdnum --with-root
# option to use python interface to library: --enable-pyext
#     ./configure --prefix=${base} --enable-pyext
      if ( $withoptional ) then
         ./configure --prefix=${base} --with-yoda --with-hoppet --with-qcdnum ${rootoptpath} ${pyextopt}
      else
        ./configure --prefix=${base} ${rootoptpath} ${pyextopt}
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
   if ( $withnlojetpp ) then
#
# NLOJet++ (patched version from fastNLO web page is required: NLOJet-4.1.3-patched2):
#------------------------------------------------------------------------------
      set arc="nlojet++-4.1.3"
      if ( ! -e ${arc}_installed  ) then
         tar xzf ${arc}-patched2.tar.gz
         cd ${arc}
         ./configure --prefix=${base}
         make -j${cores} install
         cd ..
         touch ${arc}_installed
      endif
#
# fastNLO Interface NLOJet++:
#------------------------------------------------------------------------------
#      set arc="fastnlo_interface_nlojet-2.3.1pre-2424"
#      set arc="fastnlo_interface_nlojet-2.3.1pre-2657"
#      set arc="fastnlo_interface_nlojet-2.3.1pre-2771"
      set arc="fastnlo_interface_nlojet-2.5.0-2819"
      if ( ! -e ${arc}_installed  ) then
         tar xzf ${arc}.tar.gz
         cd ${arc}
         ./configure --prefix=${base}
         make -j${cores} install
         cd ..
         touch ${arc}_installed
      endif
   endif
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
# APPLgrid also requires the static library libgfortran.a not present in numerous cases ...
#   if ( $withroot != "0" && ! $mpnnlo ) then
#      set arc="applgrid-1.5.6"
#      if ( ! -e ${arc}_installed  ) then
#         tar xzf ${arc}.tar.gz
#         cd ${arc}
#         ./configure --prefix=${base}
# Attention: No concurrent compilation with -j here!
#         make install
#         cd ..
#         touch ${arc}_installed
#      endif
#   endif
#
# nnlo-bridge to NNLOJet:
#------------------------------------------------------------------------------
   if ( ! $mpnnlo ) then
#      set arc="nnlo-bridge-0.0.40"# updated scale settings for jetpt scale
      set arc="nnlo-bridge-0.0.46"
# Previous buggy:    set rev="-rev1683M3"       ; fixed in 0.0.39
# Improve interface --> M7: set rev="-rev1683M5"; fixed in 0.0.39
#      set rev="-rev1683M7"                       ; fixed in 0.0.39
# M1: Two printout fixes
#      set rev="M1"
      set rev=""
      if ( ! -e ${arc}_installed  ) then
         tar xzf ${arc}${rev}.tar.gz
         cd ${arc}
         ./configure --prefix=${base}
         make -j${cores} install
         cd ..
         touch ${arc}_installed
      endif
   endif
#
# NNLOJET
#------------------------------------------------------------------------------
# Single-thread usage of NNLOJET
   setenv OMP_STACKSIZE 999999
   setenv OMP_NUM_THREADS 1
# or multi-thread usage of NNLOJET; values should be adapted to respective system
   if ( $mpnnlo ) then
      setenv OMP_STACKSIZE 999999
      setenv OMP_NUM_THREADS 4
   endif
   echo 'setenv OMP_STACKSIZE '"$OMP_STACKSIZE" >> fnlosrc_source.csh
   echo 'setenv OMP_NUM_THREADS '"$OMP_NUM_THREADS" >> fnlosrc_source.csh
   echo 'export OMP_STACKSIZE='"$OMP_STACKSIZE" >> fnlosrc_source.sh
   echo 'export OMP_NUM_THREADS='"$OMP_NUM_THREADS" >> fnlosrc_source.sh
   echo ""
   echo "ATTENTION: OpenMP environment set for NNLOJET!"
   echo "   OMP_STACKSIZE and OMP_NUM_THREADS have been set to:"
   echo "   $OMP_STACKSIZE and $OMP_NUM_THREADS"
   echo ""
   set arc="NNLOJET_rev"${revision}
   if ( ! -e ${arc}_installed  ) then
   tar xzf ${arc}.tar.gz
   cd ${arc}/driver
   make depend
   if ( $mpnnlo ) then
      make -j${cores}
   else
      make -j${cores} fillgrid=1
   endif
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
      ./configure --prefix=${base} --enable-shared CXX=g++ CC=gcc FC=gfortran
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
      ./configure --prefix=${base}
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
      ./configure --prefix=${base}
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
# (Repacked version required for compatibility with C++11 standard enforcement.
#  Plus additional patch in ./ATOOLS/Org/Gzip_Stream.C)
#------------------------------------------------------------------------------
   set arc="SHERPA-MC-2.2.4"
#   set arc="SHERPA-MC-2.2.5" # Does not work, not even after repacking
   if ( ! -e ${arc}_installed  ) then
      tar xzf ${arc}-patched.tar.gz
#      tar xzf ${arc}.tar.gz
      cd ${arc}
      ./configure --prefix=${base} --with-sqlite3=install --enable-gzip --enable-lhole --enable-fastjet=${base} --enable-lhapdf=${lhapdfbasepath} --enable-hepmc2=${base} --enable-rivet=${base} ${rootenablepath} --enable-openloops=${base}/src/OpenLoops-1.3.1 ${mpiopt} --enable-blackhat=${base}
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
      tar xzf ${arc}-patched2.tar.gz
      cd ${arc}
      ./configure --prefix=${base}
      make -j${cores} install # CXXFLAGS=-fpermissive
      cd ..
      touch ${arc}_installed
   endif
#
# Environment variables for Rivet with Sherpa/MCgrid/fastNLO
#------------------------------------------------------------------------------
   setenv RIVET_ANALYSIS_PATH ${base}/share/Rivet:${base}/share/fastnlo_interface_nlojet/RivetAdditions
   setenv RIVET_INFO_PATH ${base}/share/Rivet:${base}/share/fastnlo_interface_nlojet/RivetAdditions
   setenv RIVET_REF_PATH ${base}/share/Rivet:${base}/share/fastnlo_interface_nlojet/RivetAdditions
   setenv PKG_CONFIG_PATH ${base}/lib/pkgconfig
   echo 'setenv RIVET_ANALYSIS_PATH '"${base}/share/Rivet:${base}/share/fastnlo_interface_nlojet/RivetAdditions" >> fnlosrc_source.csh
   echo 'setenv RIVET_INFO_PATH '"${base}/share/Rivet:${base}/share/fastnlo_interface_nlojet/RivetAdditions" >> fnlosrc_source.csh
   echo 'setenv RIVET_REF_PATH '"${base}/share/Rivet:${base}/share/fastnlo_interface_nlojet/RivetAdditions" >> fnlosrc_source.csh
   echo 'setenv PKG_CONFIG_PATH '"${base}/lib/pkgconfig" >> fnlosrc_source.csh
   echo 'export RIVET_ANALYSIS_PATH='"${base}/share/Rivet:${base}/share/fastnlo_interface_nlojet/RivetAdditions" >> fnlosrc_source.sh
   echo 'export RIVET_INFO_PATH='"${base}/share/Rivet:${base}/share/fastnlo_interface_nlojet/RivetAdditions" >> fnlosrc_source.sh
   echo 'export RIVET_REF_PATH='"${base}/share/Rivet:${base}/share/fastnlo_interface_nlojet/RivetAdditions" >> fnlosrc_source.sh
   echo 'export PKG_CONFIG_PATH='"${base}/lib/pkgconfig" >> fnlosrc_source.sh
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
# Will be installed in ${base}/share/MCgrid2Examples-2.2.0
#------------------------------------------------------------------------------
   set arc="MCgrid2Examples-2.2.0"
   if ( ! -e ${arc}_installed  ) then
      tar xzf ${arc}.tar.gz -C ${base}/share
      if ( ! -d ${base}/share/${arc} ) then
         mv ${base}/share/examples ${base}/share/${arc}
      else
         echo "Warning: Could not move ${base}/share/examples!"
         echo "Directory ${base}/share/${arc} exists already."
         echo "Testing already present MCgrid plugin CMS_2011_S9086218."
      endif
      cd ${base}/share/${arc}/CMS_2011_S9086218
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
#./configure --prefix=${base} --with-fastjet=${base} --with-hepmc=${base} --with-lhapdf=${base} --with-rivet=${base}
#make -j8
#make check
#make install
#
# NJet-matchbox:
#---------------
#./configure --prefix=${base} --enable-5jet
#make -j8
#make check
#make install
#
# Herwig++-matchbox:
#-------------------
#install gengetopt from distro
#autoreconf -i
#./configure --prefix=${base} --with-nlojet=${base} --with-njet=${base} --with-fastjet=${base}
#./configure --prefix=${base} --with-fastjet=${base} --with-madgraph=${base}/MadGraph-matchbox --with-njet=${base}
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
#   ./configure --prefix=${base} --enable-pyext --with-hoppet --with-qcdnum
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
#   ./configure --prefix=${base}
#   make -j4 install
#   cd ..
#   touch ${arc}_installed
# endif
#
#
# PYTHONPATH adaptation
#
if ( $withpython ) then
   setenv PYTHONPATHADD `find ${base}/lib* -name site-packages | tr '[:space:]' ':'`
   if ( ($withroot == 5 || $withroot == 6) && $pythonopt == "--enable-python" ) then
      setenv PYTHONPATHADD ${base}/root/lib:${PYTHONPATHADD}
   endif
   if ( $?PYTHONPATHORIG ) then
   setenv PYTHONPATH ${PYTHONPATHADD}:${PYTHONPATHORIG}
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
