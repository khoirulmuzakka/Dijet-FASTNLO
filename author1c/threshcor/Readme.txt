************************************************************************
****** calculation of threshold corrections for fastNLO
************************************************************************

*** edit fnthresh-v14.f
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
make fnthresh-v14



*** run
./fnthresh-v14  path/[input].tab path/[output0].txt path/[output1].txt path/[output2].txt  [hist.hbk]
*
* - where input is an existing fastNLO sum (LO+NLO) table from which
*   the phase space, the cuts. the scales and the xlimits are extracted
* - "output0", "output1" "output2" are the result tables in text format
*        0 - Born result
*        1 - 1-loop correction
*        2 - 2-loop correction
*   -> in the fastNLO tables we only add the 2-loop correction table!!!
*      ... the Born table can be used to check the fastNLO LO table
* 
* - "hist.hbk" ist the name for the HBOOK histogram file
*    (only meaningful for  IOPTIMIZE = 0 )


*** convert tables using "txt2raw" to raw format
 ./txt2raw [mytexttable].txr [outputrawtable].raw


*** add 2-loop table with LO and NLO raw tables to sum table
 ./nlofast-add [myLOtable].raw [myNLOtable].raw [myTHRESHOLDtable].raw sum.tab
* -> produces summary table "sum.tab" which is read by the usercode



**** don't add Born and 1-loop tables to the fastNLO table!!!!

