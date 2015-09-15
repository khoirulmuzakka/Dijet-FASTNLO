///********************************************************************
///
///     fnlo-tk-merge
///     Program to merge several fastNLO files with different or
///     identical contributions into one table.
///
///********************************************************************
// DB, 13.11.13, (re)write fnlo-merge for the toolkit
// KR, 01.06.15, adapt command line treatment to our standard

#include <cstdlib>
#include <iostream>
#include <vector>
#include <unistd.h>
#include "fastnlotk/fastNLOTable.h"
#include "fastnlotk/speaker.h"


//__________________________________________________________________________________________________________________________________
int main(int argc, char** argv) {
   //! ----------------------------------------------------------------------//
   //! Program to merge fastNLO tables.
   //! Input can be LO, NLO or other contributions (non-pert., data, etc)
   //! for the identical scenario.
   //! ----------------------------------------------------------------------//

   //! --- namespaces
   using namespace std;
   using namespace say;       //! namespace for 'speaker.h'-verbosity levels
   using namespace fastNLO;   //! namespace for fastNLO constants

   //! --- Set verbosity level
   SetGlobalVerbosity(INFO);

   //! --- Parse commmand line
   cout << _CSEPSC << endl;
   shout["fnlo-tk-merge"] << "Program Steering" << endl;
   cout << _SSEPSC << endl;
   //! --- Usage info
   string tablename;
   if (argc <= 1) {
      error["fnlo-tk-merge"] << "No table names given, but need at least three!" << endl;
      shout["fnlo-tk-merge"] << "For an explanation of command line arguments type:" << endl;
      shout["fnlo-tk-merge"] << "./fnlo-tk-merge -h" << endl;
      exit(1);
   } else {
      tablename = (const char*) argv[1];
      if (tablename == "-h") {
         shout << "" << endl;
         shout << "Usage: ./fnlo-tk-merge [arguments]" << endl;
         shout << "List of blank-separated table input files, at least two," << endl;
         shout << "   plus output file for merged table." << endl;
         shout << "" << endl;
         cout  << _CSEPSC << endl;
         return 0;
      }
   }

   //! --- Check no. of file names
   if (argc <= 3) {
      error["fnlo-tk-merge"] << "Not enough table names given, need at least three!" << endl;
      exit(1);
   }
   //! --- Check if output file already exists
   int nFiles = argc - 1;
   if (access(argv[nFiles], R_OK) == 0) {
      error["fnlo-tk-merge"]<<"Output file " << argv[nFiles] << " exists already!" << endl;
      shout["fnlo-tk-merge"]<<"Please remove it first." << endl;
      return 1;
   }

   //! --- Initialize output table
   fastNLOTable* resultTable = NULL;

   //! --- Loop over argument list and check existence of files
   int nValidTables = 0;
   for (int idxFile=0; idxFile<nFiles-1; idxFile++) {
      string path = argv[idxFile+1];
      //! --- File there?
      if (access(path.c_str(), R_OK) != 0) {
         warn["fnlo-tk-merge"]<<"Unable to access file, skipping "<<path<<endl;
      }
      //! --- OK, file exists
      else {
         //! --- Reading table
         info["fnlo-tk-merge"]<<"Reading table "<<path<<endl;
         fastNLOTable tab(path);
         // Todo: check validity of table, here!
         {
            //! Check, whether an additive contribution has number of events == -1.
            //! This indicates, that this contribution cannot be merged anymore,
            //! because the number of event normalization information got lost.
            //TODO: KR: Note, the -1 indicator is not active in fnlo-tk-append.
            //          It again writes "1" for the moment to avoid complications.
            const int nc = tab.GetNcontrib() + tab.GetNdata();
            for ( int ic=0 ; ic<nc; ic++ ) {
               if ( tab.GetCoeffTable(ic)->GetIAddMultFlag()==0) {
                  fastNLOCoeffAddBase* cadd = (fastNLOCoeffAddBase*)tab.GetCoeffTable(ic);
                  if ( cadd->GetNevt(0,0) == -1 ) {
                     error["fnlo-tk-merge"]<<"Contribution #"<<ic<<" in table "<<path<<endl;
                     error>>"     cannot be merged, because no valid number-of-events information"<<endl;
                     error>>"     is available: Nevt = " << cadd->GetNevt(0,0)<<endl;
                     //                     error>>"     Please use program fnlo-tk-append instead."<<endl;
                     error>>"     Exiting."<<endl;
                     exit(1);
                  }
               }
            }
         }
         if ( !resultTable ) {
            resultTable = new fastNLOTable(tab);
            nValidTables++;
            //! Hint to use fnlo-tk-merge correctly.
            warn["fnlo-tk-merge"]<<"The user has to ensure, that merged input tables are statistically independent."<<endl;
            warn["fnlo-tk-merge"]<<"This is not checked by the program!"<<endl;
         }
         //! adding table to result table
         else {
            //! check if 'scenario' is compatible
            if ( !resultTable->IsCompatible(tab) )
               warn["fnlo-tk-merge"]<<"Table '"<<path<<"' is not compatible with initial table '"<<resultTable->GetFilename()<<"'. Skipping table."<<endl;
            //! adding tables
            else {
               resultTable->AddTable(tab);
               nValidTables++;
            }
         }
      }
   }
   info["fnlo-tk-merge"]<<"Found "<<nValidTables<<" table file(s)."<<endl;
   if (nValidTables < 2) {
      error["fnlo-tk-merge"]<<"Found less than two valid tables, no merging possible!"<<endl;
      exit(1);
   }

   //! Write result
   string outfile = argv[nFiles];
   resultTable->SetFilename(outfile);
   info["fnlo-tk-merge"]<<"Write merged results to file "<<resultTable->GetFilename()<<"."<<endl;
   resultTable->WriteTable();
   return 0;
}
