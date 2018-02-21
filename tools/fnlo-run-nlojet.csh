#!/bin/csh -f
echo "######################################";
echo "# fnlo-run-nlojet.csh: Sourcing environment for NLOJet++ with grids.";
echo "######################################";
# Source environment for local or NEMO cluster jobs with cvmfs software installation
source /cvmfs/etp.kit.edu/fnlo/fnlosrc_source.csh
# Then call perl execution script
echo "Now calling perl script with arguments $1, $2, $3, $4, $5, $6"
./fnlo-run-nlojet.pl -b CVMFS -e $1 -o $2 -n $3 -t . -v2.3 -x $4 $5 $6
