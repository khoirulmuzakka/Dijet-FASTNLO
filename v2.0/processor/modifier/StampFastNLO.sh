#/sh/bash
#
#  Stamp FastNLO file as official table
#  Call
#  StampFastNLO.sh table.tab
#
#  a new table in the same directory as table.tab is created
#  with the subscript '_Stamp'  
# 

Scenario=`echo "$1" | sed 's/.tab//g'`
fnlo-modify InTable=$Scenario.tab OutTable=${Scenario}_Stamp.tab steerfile=SteerFastNLOStamp.str
echo 'StampFastNLO.sh.  Replacing initial file with stamped file.'
mv -f ${Scenario}_Stamp.tab $Scenario.tab