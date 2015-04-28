#
# Get update script from HepForge
#
wget fastnlo.hepforge.org/setup_fastNLO.sh
# or from KIT Karlsruhe as backup
# wget www.ekp.kit.edu/~rabbertz/PDFSchool/setup_fastNLO_KIT.sh
# source script to execute and set proper environment variables
source setup_fastNLO.sh
cd
#
# What has been installed ?
#
cd Programs
#
# Let´s have a look the fastNLO Toolkit ...
#
cd fastnlo_toolkit-2.3.1pre-1896
./configure -h
head config.log
# make -j2
make install
#
# What has been installed ?
# a libary
cd ~/local/lib
# some executables
cd ../bin
# header files reside in ../include/fastnlotk
#
# Now to the interface to NLOJet++
#
cd ~/Programs/fastnlo_interface_nlojet-2.3.1pre-1897
./configure -h
head config.log
# What kind of code needed ?
cd interface
#
# A sample table production needs two ingredients
# First let´s have a look at the steering file ...
#
# emacs, vim, nano ... fnl2342b_v23_fix.str
cd data
emacs fnl2342b_v23_fix.str &
#
# And now the C++ code to be given to NLOJet++
#
# In emacs (or otherwise), go to ../hadron/InclusiveJets.cc
#
# The example can be steered to do any hh inclusive jets table
# How to include another type of analysis ?
#
# Certainly you do not want to edit Makefiles ...
# In emacs (or otherwise), go to ../hadron/Makefile.am
# Installing an autotools package only requires make and the compilers/libs
# Adding a scenario by modifying the Makefile.am and reconfiguring requires the autotools to be available
# (Don´t try in this VM)
# But if available it´s basically a simple copy/paste renaming the example.
#
# Leaving the source code directory ...
# What has been installed by the interface package
# libraries and the steering file (data)
cd ~/local/lib/fastnlo_interface_nlojet
cd ~/local/share/fastnlo_interface_nlojet
#
# In this case there is NO executable here.
#
# The NLO calculation is done by NLOJet++ (Z. Nagy), which has been installed
# already for us in the VM
#
# We also downloaded some precalculated stuff for demonstration
# Let´s go to the scratch area ...
#
cd ~/fastNLO_scratch
#
# We have defined the binning for our observable, so we can repeat
# NLOJet++ calculations all the time ... but this is too long to
# be included in fits with different PDFs, scales, strong coupling constants
#
# We have not yet defined our interpolation grid to facilitate the
# REPEATED evaluations of higher-order calculations. Some relevant settings
# we have seen in the "normally don´t change part" of the steering file
#
# The interpolation grid(s) live in x and mu space for each relevant subprocess
# of a reaction. What "relevant" means is an interesting reaction dependent question.
# For jet production we optimimized this a long time ago when we developed a proof
# of principle implementation in NLOJet++, so no worry here.
# Since this optimization is quite a lot of work, APPLGRID have an interface to MCFM for
# example, recent work is going into automizing this step (MCGrids,aMCfast@NLO,...)
#
# Let´s get back to the remaining part. We have to find out which phase space in
# fractional momenta x and scales mu is accessed for each bin of our observable
# Running a scenario for the first time a warmup run is required with sufficient no.
# of events (@NLO) to protocol the x and mu values ... then write min/max to file
# for optimization.
#
# Let´s see such a file:
#
emacs fastNLO-warmup.txt
#
# OK, now let´s call nlojet++. Here´s the command line:
#
nlojet++ --calculate -cborn --max-event=1000000 -n PDFSchool -u ../local/lib/fastnlo_interface_nlojet/libInclusiveJets.la
#
# This is a LO run
#
# Now NLO:
#
nlojet++ --calculate -cnlo --max-event=100000 -n PDFSchool -u ../local/lib/fastnlo_interface_nlojet/libInclusiveJets.la
#
# Takes a lot longer and is BY FAR not sufficient for any reasonable NLO cross section!
#
# Let´s have a look at such a table in our format (extension tab, simple txt file)
# The tables are written by nlojet++ into the directory output
#
cd output
#
# emacs ...
#
# Now what you really want to do is to run hundreds of these jobs in parallel and combine
#
cd ..
#
# I brought some production tables with me ...
# Important remark: Make sure to have sufficiently random seeds for each job!
# You can use the "-s seednumber" option for this. The inherent randomization of NLOJet++
# is not sufficient, if you have hundreds of jobs starting at the same time ...
#
# OK, let´s combine the three tables I have brought with me. This is done with
# fnlo-tk-merge
fnlo-tk-merge InclusiveJets_fnl2342b_v23_fix-hhc-born-2jet_0000.tab InclusiveJets_fnl2342b_v23_fix-hhc-nlo-2jet_0050.tab InclusiveJets_fnl2342b_v23_fix-hhc-nlo-2jet_0399.tab InclusiveJets_fnl2342b_v23_fix_merged.tab
#
# The table fnl2342b_I902309_HepForge.tab is one of our final tables with many more events.
# Let´s have a look at the predicted cross sections from the three tables:
#
cd output
# fnlo-tk-merge lo.tab nlo.tab InclusiveJets_fnl2342b_v23_fix_school.tab
mv InclusiveJets_fnl2342b_v23_fix_school.tab ..
cd ..
fnlo-tk-cppread InclusiveJets_fnl2342b_v23_fix_school.tab CT10nlo.LHgrid > school.log
fnlo-tk-cppread InclusiveJets_fnl2342b_v23_fix_merged.tab CT10nlo.LHgrid > combined.log
fnlo-tk-cppread fnl2342b_I902309_HepForge.tab CT10nlo.LHgrid             > hepforge.log
#
# At last a little comparison
#
kdiff3 school.log combined.log &
kdiff3 combined.log hepforge.log &
#
# Questions ? Problems ?
#
# And now over to Georg
#
