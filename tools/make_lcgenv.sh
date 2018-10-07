#PLATFORM=x86_64-slc6-gcc62-opt   
#PLATFORM=x86_64-slc6-gcc49-opt
#PATH=/cvmfs/sft.cern.ch/lcg/releases/automake/1.14-c6378/${PLATFORM}/bin:$PATH
#PATH=/cvmfs/sft.cern.ch/lcg/releases/autoconf/2.69-2b964/${PLATFORM}/bin:$PATH
#PATH=/cvmfs/sft.cern.ch/lcg/releases/libtool/2.4.2-9ad34/${PLATFORM}/bin:$PATH
#PATH=/cvmfs/sft.cern.ch/lcg/releases/m4/1.4.17-3bfe7/${PLATFORM}/bin:$PATH
#LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/libtool/2.4.2-9ad34/${PLATFORM}/lib:$LD_LIBRARY_PATH

####
echo ""
echo "------------------------"
echo "  generating lcgenv.sh"
echo "------------------------"
echo ""
LCGENV_PATH=/cvmfs/sft.cern.ch/lcg/releases
#cat $LCGENV_PATH/HEAD/Readme.md
#LCG_VERSION=LCG_93
#LCG_VERSION=LCG_94rc1
LCG_VERSION=LCG_93c
PLATFORM=x86_64-slc6-gcc62-opt   
#PLATFORM=${PLATFORM}
#$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} gcc >> lcgenv.sh


echo "# setup new env" > lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} m4 >> lcgenv.sh  # included in automake
#$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} libtool >> lcgenv.sh  # included in root
$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} automake >> lcgenv.sh 
#$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} autoconf >> lcgenv.sh  # included in automake
$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} ROOT >> lcgenv.sh
$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} doxygen >> lcgenv.sh
$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} swig >> lcgenv.sh
$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} zlib >> lcgenv.sh
$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} fastjet >> lcgenv.sh
$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} yoda 1.7.0 >> lcgenv.sh
$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} Boost >> lcgenv.sh
$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} yamlcpp >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} lhapdf 6.2.0 >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} lhapdf 6.1.6 >> lcgenv.sh
#$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} lhapdf 6.2.1 >> lcgenv.sh
$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} lhapdf 6.1.6.cxxstd >> lcgenv.sh
echo "export LHAPDF_DATA_PATH=\"\$LHAPDF_DATA_PATH:\$(lhapdf-config --prefix)/share/LHAPDF\"" >> lcgenv.sh # add path manually, such that .conf file is found
#$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} lhapdfsets >> lcgenv.sh
$LCGENV_PATH/lcgenv/latest/lcgenv -p $LCG_VERSION ${PLATFORM} hoppet >> lcgenv.sh
##source lcgenv.sh

echo "done. see: lcgenv.sh"
