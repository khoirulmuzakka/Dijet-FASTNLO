
 ************************************************************
 **********      author tools for fastNLO     ***************
 ************************************************************



add.sh
  MW 2005-12-27
   - flexible shell script to add raw tables into sum-table
   - needs to reside in directory together with raw tables and nlofast-add
   - user needs to edit 5 lines to define
     * first part of LO table name
     * last part of LO table name
     * first part of NLO table name
     * last part of NLO table name
     * name of result table
   --> the LO and NLO table names are expected to be of the form
           [1st part]nn[last part]
       where nn are two digits in the range 00-99 



addstaterr.sh 
  MW 2005-12-27
   - flexible shell script to create one sum-table for each individual
     raw table, as needed to determine the statistical errors from the
     fluctuations of the single results
   - needs to reside in directory together with raw tables and nlofast-add
   - user needs to edit 4 lines to define
     * first part of LO table name
     * last part of LO table name
     * first part of NLO table name
     * last part of NLO table name
   --> the LO and NLO table names are expected to be of the form
           [1st part]nn[last part]
       where nn are two digits in the range 00-99 
   --> for each LO/NLO raw table a summary table is created which 
       ends with '.stc' (as in "statistical")
       - for the LO tables a single NLO table is added (as a dummy) -
         this is the NLO table No. 01  (make sure it exists!!)
       - for each NLO table *all* LO tables are added - because we don't
         determine the stat. errors for the NLO correction, but for the
         sum of (LO+NLO) - therefore we use the most precise result
         for LO.



staterr.f       
  MW 2005-12-27
    - This is a Fortran routine that can be copied into the user-directory.
      If "example01.f" is replaced by this routine, one can compute the
      statistical errors using the tables created by "addstaterr.sh"
    - Four lines need to be edited by the user:
      * the number of LO tables
      * the num er of NLO tables
      * the filename of the LO tables
      * the filename of the NLO tables
    - Output (for each scale for each bin):
         Bin Number 
         Mean Value (= cross section)
         Stat Error (in percent)
         largest single lower deviation (in percent)
         largest single upper deviation (in percent)
    >>>>
    >>>> at present this routine does not consider that the tables
         will have different statistics!!! They are all treated on 
         the same footing.
         Please note that this will result in a larger quoted error
         (so we are safe!!)
    >>>> to be added later: 
         - code to fill the results into histograms for plotting
    >>>> the logic should be checked again (It seems that the errors
         are quite small - on the other hand I used high statistics! 
         MW-051227)



######################################
# fastrun.pl: Starting run of fastNLO: FASTRUN0_15052008_223606
######################################


fastrun.pl
Usage: fastrun.pl [switches/options] ([scenario])
  -b batch        Batch system used: LOCAL (def.), GRID or PBS
  -d dir          Installation directory (def.=.)
  -e max-events   Maximal number of events (def.=0 => 4,294,967,295)
  -f rev          fastNLO revision to use (def.=187)
  -h              Print this text
  -j jobnr        Job number to attach (def.=0001)
  -m mode         Job mode: 0 do all (def.), 1 install only, 2 make only, 3 run only
  -o order        LO (def.) or NLO calculation
  -p pdf          CTEQ parton densities (def.) or LHAPDF
  -r              Reference calculation incl. pdf access
  -s dir          Archive source directory (def.=.)
  -t dir          Output target directory: (def.= {scen}{ref}_{jobnr} with
                  ref. to working directory in fastNLO installation)
  -v              Switch verbose mode on

Examples:
1) Install only (to install with LHAPDF use option -p):
   ./fastrun.pl [-d .|installdir] [-f 187|rev] -m 1 [-p CTEQ|LHAPDF] [-s .|sdir]

2) Make only scenario (to make scenario for reference mode use option -r):
   ./fastrun.pl [-d .|installdir] [-f 187|rev] -m 2 [-p CTEQ|LHAPDF] [-r] scenarioname

3) Run only (to run scenario in reference mode use option -r):
   ./fastrun.pl [-b LOCAL|GRID|batch] [-d .|installdir] [-e max-events] [-f 187|rev] -m 3 [-p CTEQ|LHAPDF] [-r] [-t ./{scen}{ref}_{jobnr}|tdir] scenarioname



#####################################################
# fastprep.pl: Starting archive creation for fastNLO: FASTPREP_15052008_223611
#####################################################


fastprep.pl
Usage: fastprep.pl [switches/options]
  -h              Print this text
  -p pdf          Add CTEQ parton densities or LHAPDF



######################################################
# fastidcheck.pl: Starting table id check for fastNLO: FASTIDCHECK_15052008_223624
######################################################


fastidcheck.pl
Usage: fastidcheck.pl glob (selects all files matching glob)
  -h              Print this text



##################################################
# fastadd.pl: Starting table addition for fastNLO: FASTADD_15052008_223633
##################################################


fastadd.pl
Usage: fastadd.pl [switches/options] scenario
  -h              Print this text
  -l dir          Directory for LO tables, (def.=scenario)
  -n dir          Directory for NLO tables, (def.=scenario)
  -s              Produce tables for statistical evaluation,
                  i.e. combinations of each LO with 1 NLO table and
                  all LO with each NLO table
  -v              Verbose output
