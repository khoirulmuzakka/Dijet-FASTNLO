//********************************************************************
//
//     fnlo-append
//     Program to add several fastNLO files with different parts
//     to one contribution. This program also merges tables
//     with different contribtuions.
//     Whenever fnlo-append is used, the resulting table cannot be 
//     'merged' with other tables to accumulate higher statistics
//     since the 'normalization' of the single parts is lost.
//
//********************************************************************
// DB, 19.04.14, frist instance of fnlo-append

#include <cstdlib>
#include <vector>
#include <iostream>
#include <unistd.h>
#include "fastnlotk/speaker.h"
#include "fastnlotk/fastNLOTable.h"
#include "fastnlotk/fastNLOCoeffAddBase.h"

using namespace std;
using namespace say;

int main(int argc, char** argv)
{
   // ----------------------------------------------------------------------//
   //! Program to 'add' fastNLO tables (in contrast to 'merge' tables).
   //! Input can be LO, NLO or other contributions (non-pert., data, etc)
   //! for the identical scenario.
   // ----------------------------------------------------------------------//

   // check parameters
   if (argc < 3) {
      error["fnlo-append"]<<"Usage: fnlo-append file1.txt [filex.txt]+ result.txt"<<endl;
      return 1;
   }

   // check if output file already exists
   int nFiles = argc - 1;
   if (access(argv[nFiles], R_OK) == 0) {
      error["fnlo-append"]<<"Error: Output file " << argv[nFiles] << " already exists!" << endl;
      return 1;
   }

   // init output table
   fastNLOTable* resultTable = NULL;

   // loop over arguments and check existence of files
   int nValidTables = 0;
   for (int idxFile=0; idxFile<nFiles-1; idxFile++) {
      string path = argv[idxFile+1];
      // File there?
      if (access(path.c_str(), R_OK) != 0) {
         warn["fnlo-append"]<<"Unable to access file. Skipping "<<path<<endl;
      }
      else { //ok, file exists
         // reading table
         fastNLOTable tab(path);
	 // normalize coefficients
	 const int nc = tab.GetNcontrib() + tab.GetNdata();
	 for ( int ic=0 ; ic<nc; ic++ ) {
	    // is additive?
	    if ( tab.GetCoeffTable(ic)->GetIAddMultFlag()==0) {
	       fastNLOCoeffAddBase* cadd = (fastNLOCoeffAddBase*)tab.GetCoeffTable(ic);
	       cadd->NormalizeCoefficients();
	    }
	 }
	 // Todo: check validity of table, here!
         if ( !resultTable ) {
            resultTable = new fastNLOTable(tab);
	    nValidTables++;
	    // Hint to use fnlo-append correctly.
	    info>>"\n";
	    info["fnlo-append"]<<"\n";
	    info>>"     The user has to ensure, that all input tables contain different parts"<<endl;
	    info>>"     to one contribution or are different contributions."<<endl;
	    info>>"     This is not checked by the program."<<endl;
	    info>>"     The 'number of event' information is lost in the resulting tables."<<endl;
	    info>>"     Therefore, the resulting table cannot be 'merged' with other statistically"<<endl;
	    info>>"     independent tables any longer. If needed, this 'merging' has to be "<<endl;
	    info>>"     performed beforehand for each single input table seperately using fnlo-merge."<<endl;
	    info>>"     The user has to ensure, that the single input parts all contain sufficiently high"<<endl;
	    info>>"     numerical accuracy/statistics."<<endl;
	    info>>"\n";
	 }
         else { // adding table to result table
            // check if 'scenario' is compatible
            if ( !resultTable->IsCompatible(tab) )
               warn["fnlo-append"]<<"Table '"<<path<<"' is not compatible with initial table '"<<resultTable->GetFilename()<<"'. Skipping table."<<endl;
            else { // adding tables
               resultTable->AddTable(tab);
               nValidTables++;
            }
         }
	 // normalize all contributions
	 //resultTable->SetNumberOfEvents(1); (for all contributions)
	 const int nc2 = resultTable->GetNcontrib() + resultTable->GetNdata();
	 for ( int ic=0 ; ic<nc2; ic++ ) {
	    if ( resultTable->GetCoeffTable(ic)->GetIAddMultFlag()==0) {
	       fastNLOCoeffAddBase* cadd = (fastNLOCoeffAddBase*)resultTable->GetCoeffTable(ic);
	       cadd->SetNevt(1);
	    }
	 }
      }
   }
   info["fnlo-append"]<<"Found "<<nValidTables<<" table file(s)."<<endl;
   if (nValidTables < 1) exit(1);


   // Write result
   string outfile = argv[nFiles];
   resultTable->SetFilename(outfile);
   info["fnlo-append"]<<"Write results to file "<<resultTable->GetFilename()<<endl;
   resultTable->WriteTable();
   return 0;
}
