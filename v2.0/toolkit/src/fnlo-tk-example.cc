//********************************************************************
//
//     fastNLO
//     This is your playgroud to calculate cross sections
//     with fastNLO.
//
//********************************************************************
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "fastnlotk/fastNLOLHAPDF.h"


//______________________________________________________________________________________________________________
int main(int argc, char** argv) {
   // usage: fastNLO <fastNLO-table.tab> [<LHAPDFfile>]
   using namespace std;
   using namespace fastNLO;      // namespace for fastNLO constants

   // ---  Parse commmand line
   if (argc <= 1) {
      cout<<" Program fastNLO\n   usage:\n   fastNLO <fastNLO-table.tab> [<LHAPDF-file>]"<<endl;
      exit(1);
   }
   // --- fastNLO table
   string tablename =  (const char*) argv[1];
   //---  PDF set
   string PDFFile = "MSTW2008lo68cl.LHgrid";
   if (argc > 2)    PDFFile = (const char*) argv[2];


   //--- give some output
   cout<<" fastNLO: Evaluating table: " << tablename << endl;
   cout<<" fastNLO: Using PDF set   : " << PDFFile << endl;


   // --- this is your playgroud to use fastNLO
   //  Calculate cross setions and/or test some options
   //  For a documentation and function calls, please see
   //    './src/fnlo-cppread.cc'

   //--- example calculation
   fastNLOLHAPDF fnlo(tablename,"cteq6m.LHpdf",0);
   fnlo.CalcCrossSection();
   fnlo.PrintCrossSections();
   fnlo.PrintCrossSectionsDefault();

   return 0;
}
