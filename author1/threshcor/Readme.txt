************************************************************************
****** calculation of threshold corrections for fastNLO
************************************************************************

*** edit fnthresh.f
* set the No. of integration points
* good precision (but takes very long)
      nypt  = 8
      nptpt = 16
      nxpt  = 90
      ns4pt = 60
* for faster test (but should not be taken too seriously)
      nypt  = 5
      nptpt = 12
      nxpt  = 70
      ns4pt = 40

* full print output (including the NLO-NLL contribution) - but 4x slower
      IOPTIMIZE = 0            
* only LO output + fast table of NNLO-NLL coefficients
      IOPTIMIZE = 1            



*** make:
make fnthresh



*** run
./fnthresh path/[input].tab path/[output].txt  [hist.hbk]
* - where input is an existing fastNLO sum (LO+NLO) table from which
*   phase space, the cuts. the scales and the xlimits are extracted
* - "output" is the result table in text format
* - "hist.hbk" ist the name for the HBOOK histogram file
*    (only meaningful for  IOPTIMIZE = 0 )


*** convert table using "txt2raw"
 ./txt2raw [mytexttable].txr [outputrawtable].raw


*** add with LO and NLO raw tables to sum table
 ./nlofast-add [myLOtable].raw [myNLOtable].raw [myTHRESHOLDtable].raw sum.tab
* -> produces summary table "sum.tab" which is read by the usercode


