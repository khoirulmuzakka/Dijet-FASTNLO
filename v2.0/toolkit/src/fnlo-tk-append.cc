///********************************************************************
///
///     fastNLO_toolkit: fnlo-tk-append
///     Tool to add several fastNLO tables with different parts
///     (real, virtual, ...) of the SAME additive contribution.
///
///     For more explanations type:
///     ./fnlo-tk-append -h
///
///     D. Britzger, K. Rabbertz
///
///********************************************************************
// DB, 19.04.14, first instance of fnlo-tk-append

#include <cstdlib>
#include <vector>
#include <iostream>
#include <unistd.h>
#include "fastnlotk/fastNLOTable.h"
#include "fastnlotk/fastNLOCoeffAddBase.h"
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
   info["fnlo-tk-append"] << "Tool to add fastNLO tables with different parts" << endl;
   info["fnlo-tk-append"] << "(real, virtual, ...) of the SAME additive contribution." << endl;
   yell << _SSEPSC << endl;
   info["fnlo-tk-append"] << "For more explanations type:" << endl;
   info["fnlo-tk-append"] << "./fnlo-tk-append -h" << endl;
   yell << _CSEPSC << endl;

   //! --- Parse commmand line
   yell << "" << endl;
   yell << _CSEPSC << endl;
   shout["fnlo-tk-append"] << "fastNLO Table Contribution Appender"<<endl;
   yell << _SSEPSC << endl;
   string tablename;
   if (argc <= 1) {
      error["fnlo-tk-append"] << "No table names given, but need at least three!" << endl;
      shout["fnlo-tk-append"] << "For an explanation of command line arguments type:" << endl;
      shout["fnlo-tk-append"] << "./fnlo-tk-append -h" << endl;
      yell << _CSEPSC << endl;
      exit(1);
   } else {
      tablename = (const char*) argv[1];
      //! --- Usage info
      if (tablename == "-h") {
         yell << " #" << endl;
         info["fnlo-tk-append"] << "The purpose of this tool is the summation of different parts" << endl;
         info["fnlo-tk-append"] << "(real, virtual, ...) of the SAME additive contribution, which may" << endl;
         info["fnlo-tk-append"] << "have been calculated separately for computational optimisation purposes." << endl;
         info["fnlo-tk-append"] << "It is the user's responsibility to ensure, that the input tables contain" << endl;
         info["fnlo-tk-append"] << "the final results of the different parts that are to be added."<<endl;
         info["fnlo-tk-append"] << "This can not be checked by the program."<<endl;
         info["fnlo-tk-append"] << "Furthermore, it is assumed that statistically independent tables for the SAME parts" << endl;
         info["fnlo-tk-append"] << "have been combined beforehand using 'fnlo-tk-merge', because" << endl;
         info["fnlo-tk-append"] << "differing event numbers/normalisation between the corresponding parts" << endl;
         info["fnlo-tk-append"] << "lead to a loss of the statistical information required for further merging." << endl;
         info["fnlo-tk-append"] << "Event numbers, which are stored contribution- and not bin-wise, are set to 1."<<endl;
         man << "" << endl;
         man << "Usage: ./fnlo-tk-append <InTable_1.tab> <InTable_2.tab> [InTable_n.tab] <OutTable.tab>" << endl;
         man << "       Specification: <> mandatory; [] optional." << endl;
         man << "       List of blank-separated table files, at least three!" << endl;
         man << "       Mandatory are:" << endl;
         man << "<InTable_1.tab>:   First table input file, to which observable bins are catenated" << endl;
         man << "<InTable_2.tab>:   Second table input file, from which observable bins are catenated" << endl;
         man << "<OutTable.tab>:    Output filename, to which the table with catenated observable bins is written" << endl;
         yell << " #" << endl;
         yell << _CSEPSC << endl;
         return 0;
      }
   }

   //! --- Check no. of file names
   if (argc <= 3) {
      error["fnlo-tk-append"] << "Not enough table names given, need at least three!" << endl;
      exit(1);
   }

   //! --- Check if output file already exists
   int nFiles = argc - 1;
   if (access(argv[nFiles], R_OK) == 0) {
      error["fnlo-tk-append"]<<"Output file " << argv[nFiles] << " exists already!" << endl;
      shout["fnlo-tk-append"]<<"Please remove it first." << endl;
      exit(1);
   }

   //! --- Initialize output table
   fastNLOTable* resultTable = NULL;
   string outfile = argv[nFiles];

   //! loop over argument list and check existence of files
   int nValidTables = 0;
   for (int idxFile=0; idxFile<nFiles-1; idxFile++) {
      string path = argv[idxFile+1];
      //! File there?
      if (access(path.c_str(), R_OK) != 0) {
         warn["fnlo-tk-append"]<<"Unable to access file '" << path << "', skipped!" << endl;
      }
      //! --- OK, file exists
      else {
         //! Reading table
         info["fnlo-tk-append"]<<"Reading table '" << path << "'" << endl;
         yell << _CSEPSC << endl;
         fastNLOTable tab(path);

         //! --- Always normalise additive contributions from 'new'-table
         for ( int ic=0; ic<tab.GetNcontrib()+tab.GetNdata(); ic++ ) {
            bool quiet = true;
            fastNLOCoeffBase* cnew = (fastNLOCoeffBase*)tab.GetCoeffTable(ic);
            // Identify type of new coeff table
            // Only additive ones have event numbers
            // Additive?
            if ( fastNLOCoeffAddBase::CheckCoeffConstants(cnew,quiet) ) {
               fastNLOCoeffAddBase* cadd = (fastNLOCoeffAddBase*)tab.GetCoeffTable(ic);
               cadd->NormalizeCoefficients();
            }
         }

         //! --- Initialising result with first read table
         if ( !resultTable ) {
            info["fnlo-tk-append"]<<"Initialising normalised result table '" << outfile << "'" << endl;
            resultTable = new fastNLOTable(tab);
            nValidTables++;
         }
         //! --- Adding table to result table
         else {
            //! check if 'scenario' is compatible
            if ( !resultTable->IsCompatible(tab) )
               warn["fnlo-tk-append"]<<"Table '"<<path<<"' is not compatible with initial table '"<<resultTable->GetFilename()<<"', skipped!"<<endl;
            //! adding tables
            else {
               resultTable->AddTable(tab);
               nValidTables++;
            }
         }
         //! Set no. of events for all additive result contributions
         for ( int jc=0; jc<resultTable->GetNcontrib()+resultTable->GetNdata(); jc++) {
            bool quiet = true;
            fastNLOCoeffBase* cres = (fastNLOCoeffBase*)resultTable->GetCoeffTable(jc);
            // Additive?
            if ( fastNLOCoeffAddBase::CheckCoeffConstants(cres,quiet) ) {
               fastNLOCoeffAddBase* cadd = (fastNLOCoeffAddBase*)resultTable->GetCoeffTable(jc);
               cadd->SetNevt(1);
            }
         }
      }
   }
   info["fnlo-tk-append"]<<"Found "<<nValidTables<<" table file(s)."<<endl;
   if (nValidTables < 2) exit(1);

   //! Write result
   resultTable->SetFilename(outfile);
   info["fnlo-tk-append"]<<"Write added results to file "<<resultTable->GetFilename()<<endl;
   resultTable->WriteTable();
   return 0;
}
