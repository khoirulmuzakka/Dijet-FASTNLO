///********************************************************************
///
///     fastNLO_toolkit: fnlo-tk-merge
///     Tool to merge fastNLO tables with different contributions or
///     to combine identical statistically independent contributions
///
///     For more explanations type:
///     ./fnlo-tk-merge -h
///
///     D. Britzger, K. Rabbertz
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

   //! --- namespaces
   using namespace std;
   using namespace say;          //! namespace for 'speaker.h'-verbosity levels
   using namespace fastNLO;      //! namespace for fastNLO constants

   //! --- Set verbosity level
   SetGlobalVerbosity(INFO);

   //! --- Print program purpose
   yell << _CSEPSC << endl;
   info["fnlo-tk-merge"] << "Tool to merge fastNLO tables with different contributions or" << endl;
   info["fnlo-tk-merge"] << "to combine identical statistically independent contributions" << endl;
   yell << _SSEPSC << endl;
   info["fnlo-tk-merge"] << "For more explanations type:" << endl;
   info["fnlo-tk-merge"] << "./fnlo-tk-merge -h" << endl;
   yell << _CSEPSC << endl;

   //! --- Parse commmand line
   yell << "" << endl;
   yell << _CSEPSC << endl;
   shout["fnlo-tk-merge"] << "fastNLO Table Merger"<<endl;
   yell << _SSEPSC << endl;
   string tablename;
   if (argc <= 1) {
      error["fnlo-tk-merge"] << "No table names given, but need at least three!" << endl;
      shout["fnlo-tk-merge"] << "For an explanation of command line arguments type:" << endl;
      shout["fnlo-tk-merge"] << "./fnlo-tk-merge -h" << endl;
      yell << _CSEPSC << endl;
      exit(1);
   } else {
      tablename = (const char*) argv[1];
      //! --- Usage info
      if (tablename == "-h") {
         yell << " #" << endl;
         info["fnlo-tk-merge"] << "The purpose of this tool is to merge fastNLO tables with different contributions or" << endl;
         info["fnlo-tk-merge"] << "to combine identical statistically independent additive contributions to improve" << endl;
         info["fnlo-tk-merge"] << "the statistical precision." << endl;
         info["fnlo-tk-merge"] << "The statistical information of each additive contribution is checked." << endl;
         info["fnlo-tk-merge"] << "An event number of unity indicates that this contribution" << endl;
         info["fnlo-tk-merge"] << "has been combined from multiple contributions losing the" << endl;
         info["fnlo-tk-merge"] << "the event normalisation information that is stored contribution-" << endl;
         info["fnlo-tk-merge"] << "and not bin-wise. Further merging of such tables is not possible." << endl;
         man << "" << endl;
         man << "Usage: ./fnlo-tk-merge <InTable_1.tab> <InTable_2.tab> [InTable_n.tab] <OutTable.tab>" << endl;
         man << "       Specification: <> mandatory; [] optional." << endl;
         man << "       List of blank-separated table files, at least three!" << endl;
         man << "       Mandatory are:" << endl;
         man << "<InTable_1.tab>:   First table input file to be merged" << endl;
         man << "<InTable_2.tab>:   Second table input file to be merged" << endl;
         man << "<OutTable.tab>:    Output filename, to which the merged table is written" << endl;
         yell << " #" << endl;
         yell << _CSEPSC << endl;
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
      exit(1);
   }

   //! --- Initialize output table
   fastNLOTable* resultTable = NULL;
   string outfile = argv[nFiles];

   //! --- Loop over argument list and check existence of files
   int nValidTables = 0;
   for (int idxFile=0; idxFile<nFiles-1; idxFile++) {
      string path = argv[idxFile+1];
      //! --- File there?
      if (access(path.c_str(), R_OK) != 0) {
         warn["fnlo-tk-merge"]<<"Unable to access filec'" << path << "', skipped!" << endl;
      }
      //! --- OK, file exists
      else {
         //! --- Reading table
         info["fnlo-tk-merge"]<<"Reading table '" << path << "'" << endl;
         yell << _CSEPSC << endl;
         fastNLOTable tab(path);

         //! --- Check statistical information of additive contributions
         int nc = tab.GetNcontrib() + tab.GetNdata();
         for ( int ic=0 ; ic<nc; ic++ ) {
            bool quiet = true;
            fastNLOCoeffBase* cnew = (fastNLOCoeffBase*)tab.GetCoeffTable(ic);
            // Identify type of new coeff table
            // Only additive ones have event numbers
            // Additive?
            if ( fastNLOCoeffAddBase::CheckCoeffConstants(cnew,quiet) ) {
               fastNLOCoeffAddBase* cadd = (fastNLOCoeffAddBase*)tab.GetCoeffTable(ic);
               if ( cadd->GetNevt() == 1 ) {
                  error["fnlo-tk-merge"]<<"Contribution #" << ic << " in table " << path << endl;
                  error["fnlo-tk-merge"]<<"has no valid number-of-events information and cannot be merged. Aborted!" << endl;
                  error["fnlo-tk-merge"]<<"Nevt = " << cadd->GetNevt() << endl;
                  exit(1);
               }
            }
         }

         //! --- Initialising result with first read table
         if ( !resultTable ) {
            info["fnlo-tk-merge"]<<"Initialising result table '" << outfile << "'" << endl;
            resultTable = new fastNLOTable(tab);
            nValidTables++;
         }
         //! --- Adding further tables to result table
         else {
            //! check if 'scenario' is compatible
            if ( !resultTable->IsCompatible(tab) )
               warn["fnlo-tk-merge"]<<"Table '" << path << "' is not compatible with initial table '" << resultTable->GetFilename() << "', skipped!"<<endl;
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
      error["fnlo-tk-merge"]<<"Found less than two valid tables, no merging possible. Aborted!"<<endl;
      exit(1);
   }

   //! Write result
   resultTable->SetFilename(outfile);
   info["fnlo-tk-merge"]<<"Write merged results to file "<<resultTable->GetFilename()<<"."<<endl;
   resultTable->WriteTable();
   return 0;
}
