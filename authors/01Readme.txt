

*****************************************************************
* July 5, 2005  M. Wobisch
*
* Run I D0/CDF measurements of incl. jet cross section
*****************************************************************
run author code:

first:
make 

then create a link "/path" to where your tables should be stored

then run:
../bin/nlojet++ --save-after 50000000 -cborn -n run1incl01a00 -d /path -u fnrun1incl01a.la
../bin/nlojet++ --save-after 5000000 -cnlo -n run1incl01a00 -d /path -u fnrun1incl01a.la

(stores tables under /path)




in /path       -> add/convert tables  (use new version of nlofsast-add)

./nlofast-add run1incl01a00-hhc-born-2jet.raw   run1incl01-mw-born.tab
./nlofast-add run1incl01a01-hhc-nlo-2jet.raw    run1incl01-mw-nlo.tab

(or use tables, provided by MW -> put in path)
****************************************************************
run user code (in package: fnrun1incl01-user.tar  )

set paths for CERNLIB and LHAPDF (v4) in make file

in example01.f: set path to CTEQ6.1M grid in LHAPDF directory:
    call InitPDFset("/disk2/work/wobisch/LHAPDFv4/PDFsets/cteq61.LHgrid")

make -> then run: 
./example1


(stores hbookfile under /path)

run:
paw -b fastnlo

to create plots in file: fastnlo.ps