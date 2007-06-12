#!/bin/sh
# This gets the environment
. /etc/bashrc


###########################################################################
#  this shell script creates the sum table from LO and NLO tables
#
#  the sum table will have the ending ".tab" 
#
#  The script can handle up 100 LO and NLO tables, assuming
#  that counting starts at 00
#
#  put this script in the folder where you store the raw tables and 
#  the nlofast-add
#  the user needs to edit the next *five* lines where the filnames
#  of the input-files are defined
#    --- filenames for NLO raw tables

startnlo=fn01mn
endnlo=-hhc-nlo-2jet.raw

#
#    --- filenames for LO raw tables

startborn=fn01ml
endborn=-hhc-born-2jet.raw

# result table name

result=fnt1001mid.tab


# ===================================================================
# ============ no user-editing below this line
#
echo ' '
echo ' '
echo '  ******************************************************************'
echo '  ******  starting the shell script to add the fastNLO tables ******'
echo '  ******************************************************************'
#
bornfiles=' '
nlofiles=' '
nnlofiles=' '


for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99

do
  infileborn=$startborn$i$endborn
  infilenlo=$startnlo$i$endnlo
  if [ -f $infileborn ]
  then
    bornfiles=$bornfiles' '$infileborn
  fi

  if [ -f $infilenlo ]
  then
    nlofiles=$nlofiles' '$infilenlo
  fi


done

echo '>>> the LO files are:'
echo $bornfiles
echo '--------------------' 
echo '>>> the NLO files are:'
echo $nlofiles 
echo '--------------------' 
echo '>>> the NNLO files are:'
echo $nnlofiles
echo '--------------------'
echo '>>> the result will be stored in:'
echo $result
if [ -f $result ]
then
  echo ' >>>> found existing file '$result'   - move to '$result'-old'
  mv $result $result'-old'
fi

echo ' '
echo '  ========================================================'
echo '  ===============    now run nlofast-add   ===============' 
echo '  ========================================================'
./nlofast-add $bornfiles $nlofiles $nnlofiles $result


