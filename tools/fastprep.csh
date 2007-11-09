#!/bin/csh
echo "Preparing fastNLO archive for submission in directory $FASTNLO"
echo "Only hh collisions for now."
cd $FASTNLO
tar cfz fastNLO-exe.tar.gz gcclocal.csh nlojet/bin lib/lib*.so* nlojet/lib fastNLO-rev260/trunk/common fastNLO-rev260/trunk/v1.4/author1c/hadron/*.la fastNLO-rev260/trunk/v1.4/author1c/hadron/.libs
