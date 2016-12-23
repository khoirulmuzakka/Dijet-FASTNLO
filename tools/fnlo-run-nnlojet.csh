#!/bin/csh -f
echo "######################################";
echo "# fnlo-run-nnlojet.csh: Sourcing environment for NNLOJET with grids.";
echo "######################################";
# Source environment at Freiburg for more recent GCC than SLC6 standard
set gccbase=/cvmfs/cms.cern.ch/slc6_amd64_gcc481
source ${gccbase}/external/gcc/4.8.1/etc/profile.d/init.csh
# Then call perl execution script
echo "Now calling perl script with arguments $1, $2, $3"
./fnlo-run-nnlojet.pl $1 $2 $3
