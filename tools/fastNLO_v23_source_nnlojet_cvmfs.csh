#!/bin/csh -ef
#==============================================================================
# Set installation location
#==============================================================================
set base=/cvmfs/etp.kit.edu
set local="nnlo"
set revision=3738
#==============================================================================
# Check and set defaults for some specific environment variables
#==============================================================================
# Use CVMFS software repository, if necessary
setenv MYCVMFS /cvmfs/cms.cern.ch/slc6_amd64_gcc481
echo "MYCVMFS is $MYCVMFS"
# Use newer boost headers; necessary for newer gcc
# Check correct location and version!
if ( $?MYCVMFS ) then
   setenv INCLUDE ${MYCVMFS}/external/boost/1.57.0/include
endif
# Find root-config
# Do not set ROOTBIN to crappy CVMFS ROOT installation!
#    setenv ROOTBIN ${MYCVMFS}/cms/cmssw/CMSSW_7_2_4/external/slc6_amd64_gcc481/bin
setenv ROOTBIN
# NNLOJET executable
setenv NNLOJETBINPATH ${base}/${local}/src/NNLOJET_rev${revision}/driver
# PATH adaptation
if ( $?PATH ) then
    setenv PATH ${NNLOJETBINPATH}:${base}/${local}/bin:${ROOTBIN}:${PATH}
else
    setenv PATH ${NNLOJETBINPATH}:${base}/${local}/bin:${ROOTBIN}
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
  setenv LD_LIBRARY_PATH ${base}/${local}/lib:${base}/${local}/lib/root:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH ${base}/${local}/lib:${base}/${local}/lib/root
endif
#
# NNLOJET
#------------------------------------------------------------------------------
# Set for single-thread usage of NNLOJET
setenv OMP_STACKSIZE 999999
setenv OMP_NUM_THREADS 1
#
# Finish
#------------------------------------------------------------------------------
echo "#=============================================================================="
echo "# All done! Have fun!"
echo "#=============================================================================="
