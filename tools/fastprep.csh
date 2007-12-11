#!/bin/csh
echo "Preparing fastNLO archive for submission in directory $FASTNLO"
echo "Only hh collisions for now."
cd $FASTNLO/..
tar cfz fastNLO-bin.tgz lib/lib*.so* lib64/lib*.so* nlojet/bin nlojet/lib fastNLO/trunk/common fastNLO/trunk/tools fastNLO/trunk/v1.4/author1c/hadron/*.la fastNLO/trunk/v1.4/author1c/hadron/.libs
