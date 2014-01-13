//********************************************************************
//
//     fnlo-merge
//     Program to merge several fastNLO files with different or
//     identical contributions into one table.
//
//********************************************************************
// DB, 13.11.13, (re)write fnlo-merge for the toolkit 

#include <cstdlib>
#include <vector>
#include <iostream>
#include <unistd.h>
#include "fastnlotk/speaker.h"
#include "fastnlotk/fastNLOTable.h"


using namespace std;
using namespace say;

int main(int argc, char** argv) 
{
   // ----------------------------------------------------------------------//
   // Program to merge fastNLO tables.
   // Input can be LO, NLO or other contributions (non-pert., data, etc)
   // for the identical scenario.
   // ----------------------------------------------------------------------//
   
   // check parameters
   if (argc < 3) {
      error["fnlo-merge"]<<"Usage: fnlo-merge file1.txt [filex.txt]+ result.txt"<<endl;
      return 1;
   }
   
   // check if output file already exists
   int nFiles = argc - 1;
   if (access(argv[nFiles], R_OK) == 0) {
      error["fnlo-merge"]<<"Error: Output file " << argv[nFiles] << " already exists!" << endl;
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
	 warn["fnlo-merge"]<<"Unable to access file. Skipping "<<path<<endl;
      }
      else { //ok, file exists
	 // reading table
	 fastNLOTable tab(path);
	 //! Todo: check validity of table, here!
	 if ( !resultTable )  
	    resultTable = new fastNLOTable(tab);
	 else { // adding table to result table 
	    // check if 'scenario' is compatible
	    if ( !resultTable->IsCompatible(tab) ) 
	       warn["fnlo-merge"]<<"Table '"<<path<<"' is not compatible with initial table '"<<resultTable->GetFilename()<<"'. Skipping table."<<endl;
	    else { // adding tables
	       resultTable->AddTable(tab);
	       nValidTables++;
	    }
	 }
      }
   }
   info["fnlo-merge"]<<"Found "<<nValidTables<<" table file(s)."<<endl;
   if (nValidTables < 1) exit(1);

   // Write result
   string outfile = argv[nFiles];
   resultTable->SetFilename(outfile);
   info["fnlo-merge"]<<"Write merged results to file "<<resultTable->GetFilename()<<"."<<endl;
   resultTable->WriteTable();
   return 0;
}
