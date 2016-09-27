///********************************************************************
///
///     fastNLO_toolkit: fnlo-tk-cat
///     Program to catenate observable bins of compatible fastNLO v2 tables
///
///     For more explanations type:
///     ./fnlo-tk-cat -h
///
///     K. Rabbertz
///
///********************************************************************
// KR, 26.09.2016, first instance of fnlo-tk-cat

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
   SetGlobalVerbosity(DEBUG);

   //! --- Print program purpose
   cout << _CSEPSC << endl;
   info["fnlo-tk-cat"] << "Program to catenate observable bins of compatible fastNLO v2 tables" << endl;
   cout << _SSEPSC << endl;
   info["fnlo-tk-cat"] << "For more explanations type:" << endl;
   info["fnlo-tk-cat"] << "./fnlo-tk-cat -h" << endl;
   cout << _CSEPSC << endl;

   //! --- Parse commmand line
   cout << endl;
   cout << _CSEPSC << endl;
   shout["fnlo-tk-cat"] << "fastNLO Toolkit"<<endl;
   cout << _SSEPSC << endl;
   string tablename;
   if (argc <= 1) {
      error["fnlo-tk-cat"] << "No table names given, but need at least three!" << endl;
      shout["fnlo-tk-cat"] << "For an explanation of command line arguments type:" << endl;
      shout["fnlo-tk-cat"] << "./fnlo-tk-cat -h" << endl;
      cout << _CSEPSC << endl;
      exit(1);
   } else {
      tablename = (const char*) argv[1];
      //! --- Usage info
      if (tablename == "-h") {
         cout << " #" << endl;
         shout << "Usage: ./fnlo-tk-cat <InTable_1.tab> <InTable_2.tab> [InTable_n.tab] <OutTable.tab>" << endl;
         shout << "       Specification: <> mandatory; [] optional." << endl;
         shout << "       List of blank-separated table files, at least three!" << endl;
         shout << "       Mandatory are:" << endl;
         shout << "<InTable_1.tab>:   First table input file, to which observable bins are catenated" << endl;
         shout << "<InTable_2.tab>:   Second table input file, from which observable bins are catenated" << endl;
         shout << "<OutTable.tab>:    Output filename, to which the table with catenated observable bins is written" << endl;
         cout << " #" << endl;
         cout  << _CSEPSC << endl;
         return 0;
      }
   }

   //! --- Check no. of file names
   if (argc <= 3) {
      error["fnlo-tk-cat"] << "Not enough table names given, need at least three!" << endl;
      exit(1);
   }

   //! --- Check if output file already exists
   int nFiles = argc - 1;
   if (access(argv[nFiles], R_OK) == 0) {
      error["fnlo-tk-cat"]<<"Output file " << argv[nFiles] << " exists already!" << endl;
      shout["fnlo-tk-cat"]<<"Please remove it first." << endl;
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
      if ( access(path.c_str(), R_OK) != 0 ) {
         warn["fnlo-tk-cat"]<<"Unable to access file '" << path << "', skipped!" << endl;
      }
      //! --- OK, file exists
      else {
         //! --- Reading table
         info["fnlo-tk-cat"]<<"Reading table '" << path << "'" << endl;
         fastNLOTable tab(path);

         //! --- Initialising result with first read table
         if ( !resultTable ) {
            info["fnlo-tk-cat"]<<"Initialising result table '" << outfile << "'" << endl;
            resultTable = new fastNLOTable(tab);
            cout << "First result table NObsbin = " << resultTable->GetNObsBin() << endl;
            nValidTables++;
            //! Hint to use fnlo-tk-cat correctly
            info["fnlo-tk-cat"]<<"Numerous technical checks for compatibility are performed."<<endl;
            info["fnlo-tk-cat"]<<"However, it is the user's responsibility to ensure that the same observable definition" << endl;
            info["fnlo-tk-cat"]<<"was used to produce the various observable bins to catenate."<<endl;
            info["fnlo-tk-cat"]<<"This cannot be checked by the program!"<<endl;
         }
         //! --- Catenating further tables to result table
         else {
            //! check if 'scenario' is compatible
            if ( !(resultTable->IsCatenable(tab)) )
               warn["fnlo-tk-cat"]<<"Table '" << path << "' is not catenable with initial table '" << resultTable->GetFilename() << "', skipped!" << endl;
            //! catenating tables
            else {
               cout << "AAA Following result table NObsbin = " << resultTable->GetNObsBin() << endl;
               resultTable->CatenateTable(tab);
               cout << "BBB Following result table NObsbin = " << resultTable->GetNObsBin() << endl;
               nValidTables++;
            }
         }
      }
   }
   info["fnlo-tk-cat"]<<"Found "<<nValidTables<<" table file(s)."<<endl;
   if (nValidTables < 2) {
      error["fnlo-tk-cat"]<<"Found less than two valid tables, no catenating possible!"<<endl;
      exit(1);
   }

   //! Write result
   cout << "New NObsbin = " << resultTable->GetNObsBin() << endl;


   resultTable->SetFilename(outfile);
   info["fnlo-tk-cat"]<<"Write catenated results to file '" << resultTable->GetFilename() << "'"<<endl;
   resultTable->WriteTable();
   return 0;
}
