#!/bin/sh

# Source this script via 'source setup_fastNLO.sh' in a bash
# It installs/configures all the neccesary software packages and 
# sets environment variables

CURRENT_DIR=$PWD

echo "#####################################"
echo "Adding NLOJET++ to the PATH env"
export PATH=$HOME/nlojet++-4.1.3/bin:$PATH
echo 'PATH=$HOME/nlojet++-4.1.3/bin:$PATH' > $HOME/fastNLO_env.sh


echo "#####################################"
echo "Installing Imagemagick"
echo "Required for plotting"

# Installing imagemagick and skipping sudo ahem
echo "desi1234" | sudo -S apt-get install -y imagemagick

if [ $? -eq 0 ]; then
    echo "Imagemagick succesfully installed"
else
    echo "Failed to install Imagemagick"
    return 1
fi

echo "#####################################"
echo "Reconfiguring fastjet with all plugins"
cd $HOME/Programs/fastjet-3.0.6 || return 1

./configure --prefix=/home/school/local/ --enable-allplugins || return 1
make -j2 || return 1
make install || return 1

echo "#####################################"
echo "Uninstall old fastNLO Toolkit"

cd "$HOME/Programs/fastnlo_toolkit-2.3.1pre-1871"
make uninstall || return



echo "#####################################"
echo "Installing fastNLO Toolkit"

echo "Downloading new fastNLO version"
cd "$HOME/Programs" || return 1

FASTNLO_RELEASE="fastnlo_toolkit-2.3.1pre-1896"
wget -N "http://fastnlo.hepforge.org/code/v23/${FASTNLO_RELEASE}.tar.gz" || return 1

echo "Unpacking fastNLO tarball"
tar xzf ${FASTNLO_RELEASE}.tar.gz || return 1

echo "cd to fastNLO"
cd "${FASTNLO_RELEASE}" || return 1

echo "Run configure"
./configure --prefix=/home/school/local/ --with-lhapdf=/home/school/lhapdf591 --with-yoda

echo "Run make"
make -j2 || return 1
echo "Install fastNLO"
make install || return 1

echo "Succesfully installed fastNLO toolkit"

echo "#####################################"
echo "Installing fastNLO Toolkit Generator Interface for NLOJet++"

FASTNLO_INTERFACE_RELEASE="fastnlo_interface_nlojet-2.3.1pre-1898"
cd $HOME/Programs || return 1

echo "Downloading new fastNLO version"

wget -N "http://fastnlo.hepforge.org/code/v23/${FASTNLO_INTERFACE_RELEASE}.tar.gz" || return 1

echo "Unpacking fastNLO tarball"
tar xzf ${FASTNLO_INTERFACE_RELEASE}.tar.gz || return 1

echo "cd to fastNLO"
cd "${FASTNLO_INTERFACE_RELEASE}" || return 1

echo "Run configure"
./configure --prefix=${HOME}/local --with-nlojet=$HOME/nlojet++-4.1.3 || return 1

echo "Run make"
make -j2 || return 1
echo "Install fastNLO"
make install || return 1

echo "Succesfully installed fastNLO toolkit generator interface"

echo "#####################################"
echo "Setting up Rivet"

#Rivet needs to be fixed to use system fastjet version
alias rivet="LD_LIBRARY_PATH='' rivet"
echo 'alias rivet="LD_LIBRARY_PATH='\'\'' rivet"' >> $HOME/fastNLO_env.sh
alias rivet-mkhtml="LD_LIBRARY_PATH='' rivet-mkhtml"
echo 'alias rivet-mkhtml="LD_LIBRARY_PATH='\'\'' rivet-mkhtml"' >> $HOME/fastNLO_env.sh

echo "#####################################"
echo "Downloading PDFs which will be used"

mkdir -p $HOME/lhapdf591/share/lhapdf/PDFsets
cd $HOME/lhapdf591/share/lhapdf/PDFsets || return 1

echo "Downloading CT10 NLO PDF set"
wget -N "http://www.hepforge.org/archive/lhapdf/pdfsets/current/CT10nlo.LHgrid"
wget -N "http://www.hepforge.org/archive/lhapdf/pdfsets/current/CT10nlo_as_0124.LHgrid"

echo "Downloading HERAPDF1.5 NLO set"
wget -N "http://www.hepforge.org/archive/lhapdf/pdfsets/current/HERAPDF15NLO_EIG.LHgrid"

echo "#####################################"
echo "Create scratch space for fastNLO and download sample files"

mkdir -p "$HOME/fastNLO_scratch" || return 1
cd "$HOME/fastNLO_scratch" || return 1

wget -N "http://fastnlo.hepforge.org/scenarios/tables-v21/fnl1014_I902309.tab.gz"
gunzip fnl1014_I902309.tab.gz
mv fnl1014_I902309.tab fnl2342b_I902309_HepForge.tab
wget -N "http://fastnlo.hepforge.org/PDFSchool.tar.gz"
tar xzvf PDFSchool.tar.gz

echo "#####################################"
echo "Prepare scratch area for table creation"

cd "$HOME/fastNLO_scratch"
cp -p $HOME/local/share/fastnlo_interface_nlojet/fnl2342b_v23_fix.str .
ln -s fnl2342b_v23_fix.str InclusiveJets.str 

echo "#####################################"
echo "Write out file with all environment variables"



echo "#####################################"
echo "Test table creation with NLOJet++"

cp -p InclusiveJets_fnl2342b_v23_fix_warmup.txt fastNLO-warmup.txt
nlojet++ --calculate -cborn --max-event=10000 -n Test -u ../local/lib/fastnlo_interface_nlojet/libInclusiveJets.la > /dev/null

if [ $? -eq 0 ]; then
    echo "fastNLO table creation with NLOJet++ is working."
else
    echo "There was an error running fastNLO with NLOJet++!"
    return 1
fi

fnlo-tk-cppread fnl2342b_I902309_HepForge.tab CT10nlo.LHgrid > /dev/null

if [ $? -eq 0 ]; then
    echo "fastNLO table evaluation is working."
else
    echo "There was an error running fastNLO toolkit reader."
    return 1
fi
