#!/bin/csh -f
echo "######################################";
echo "# fnlo-run-nnlojet.csh: Sourcing environment for NNLOJET with grids.";
echo "######################################";
# Source environment for local or NEMO cluster jobs with cvmfs software installation
source /cvmfs/etp.kit.edu/fnlo/fnlosrc_source.csh
# Then call perl execution script
echo "Now calling perl script with arguments $1, $2, $3, $4"
./fnlo-run-nnlojet.pl $1 $2 $3 $4
