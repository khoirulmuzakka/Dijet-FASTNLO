#!/bin/sh
# This gets the environment
. /etc/bashrc


###########################################################################
#  this shell script creates sum tables 
#    - for all single LO tables - together with one 'dummy' NLO table
#    - for all single NLO tables - together with *all* LO tables
#
#  the sum tables will have the ending ".stc" (for "statistic"
#  so far the script can handle 30 LO and 100 NLO files, assuming
#  that file counting starts at 00
#
#  put this script in the folder where you store the raw tables and 
#  the nlofast-add
#  the user needs to edit the next *four* lines where the filnames
#  of the input-files are defined
#    --- filenames for NLO raw tables
startnlo=fn01mn
endnlo=-hhc-nlo-2jet.raw
#
#    --- filenames for LO raw tables
startborn=fn01ml
endborn=-hhc-born-2jet.raw
#
#
#
# --------------------------------------------------------------------
# ----- end of user-edits - don't change the following!!!
# --------------------------------------------------------------------
#
# single NLO file  (used as dummy in LO tables)
nlofile=$startnlo'01'$endnlo
#
# a few variables
outstartborn=$startborn'-born-'
outstartnlo=$startnlo'-nlo-'
outend=.stc
#
bornfiles=' '
#
# now:
# loop over all LO tables (combined with a single NLO table)
#
echo
echo ' #######################################################################'
echo ' ######################################################################'
echo ' ## create single LO/NLO tables to study statistical fluctuations'
echo ' ## -------- now LO ----------'
echo ' ##'
for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
#
do
  outfile=$outstartborn$i$outend
  infile=$startborn$i$endborn
  if [ -f $infile ]
  then
    bornfiles=$bornfiles' '$infile
    echo ' ##          ---------'
    echo ' ##  ' $infile ' -> ' $outfile  
    if [ -f $outfile ] 
    then
      rm $outfile
    fi
    ./nlofast-add $infile $nlofile $outfile
  fi
done
#
#
echo ' ## '
echo ' ## ============================================================='
echo ' ## '
echo ' ## -------- now NLO ----------'
echo ' ## ... use the following LO files with all single NLO tables:'
echo ' ## ' $bornfiles
echo ' ## '

#
# loop over all NLO tables (combined with all LO tables)
#

for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99

do
  outfile=$outstartnlo$i$outend
  infile=$startnlo$i$endnlo
# - check if NLO raw table file exists
  if [ -f $infile ]
  then
    echo ' ##          ---------'
    echo ' ##  ' $infile ' -> ' $outfile 
    if [ -f $outfile ] 
    then
      rm $outfile
    fi
    ./nlofast-add $bornfiles $infile $outfile
  fi
done


