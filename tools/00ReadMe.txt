
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



gridrun.pl
  KR 2006-02-02
################################################
# gridrun.pl: Starting grid run of fastNLO: 02022006_1918
################################################


gridrun.pl
Usage: gridrun.pl [switches/options] scenario
  -h              Print this text
  -o order        LO (def.) or NLO calculation
  -p pdf          CTEQ parton densities (def.) or LHAPDF
  -r              Reference calculation incl. pdf access

    gridrun.pl is a perl program to install and run in one go all
    necessary ingredients of cernlib, lhapdf, nlojet and fastNLO
    as suitable for running on the grid.
    It expects to find in the directory where it is called from
    either the archives cernlib-2003.tar.gz, lhapdf-4.2.tar.gz,
    nlojet++-2.0.1.tar.gz, nlojet++-2.0.1-fix.tar.gz and 
    fastNLO-rev125.tar.gz or the corresponding subdirectories
    cernlib-1003, lhapdf-4.2, nlojet++-2.0.1 and fastNLO-rev125,
    assuming that the archives have already been installed there.

    The scenario to calculate is given as argument, e.g. fnt1003rsep.
    The default order for option -o is LO for leading order.
    By default, option -p is set to CTEQ so that LHAPDF use and
    installation is skipped.
    Option -r can be set to run the reference computations. According
    to this setting the files scenario.cc and Makefile are automatically
    adapted to use iref = 0/1 resp. -o scenario/scenarioref as name.

    For the moment, this has been set up for author1c scenarios only.
    
    If this is used as executable on the grid, it has to be accompanied
    by a jdl (job description language) file, specifying the files to
    transfer with the input sandbox, the program to run with which arguments,
    where to write STDOUT, STDERR and what files to return to the submitter
    with the output sandbox. An example is given as gridrun.jdl



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

