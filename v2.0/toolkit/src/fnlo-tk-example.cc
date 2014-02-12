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

double Function_Mu(double s1, double s2);

//______________________________________________________________________________________________________________
int main(int argc, char** argv) {
   // usage: fastNLO <fastNLO-table.tab> [<LHAPDFfile>]
   using namespace std;
   using namespace fastNLO;      // namespace for fastNLO constants

   // ---  Parse commmand line
   if (argc <= 1) {
      cout<<" Program fastNLO"<<endl;
      cout<<"   usage:"<<endl;
      cout<<"   fastNLO <fastNLO-table.tab> [<LHAPDF-file>]"<<endl;
      return 0;
   }
   // --- fastNLO table
   string tablename =  (const char*) argv[1];
   //---  PDF set
   string PDFFile = "CT10nnlo.LHgrid";
   if (argc > 2)    PDFFile = (const char*) argv[2];

   //--- give some output
   cout<<" fastNLO: Evaluating table: " << tablename << endl;
   cout<<" fastNLO: Using PDF set   : " << PDFFile << endl;


   // --- this is your playgroud to use fastNLO
   //  Calculate cross setions and/or test some options
   //  For a documentation and function calls, please see
   //    './src/fnlo-cppread.cc'

   //--- example calculation
   fastNLOLHAPDF fnlo(tablename,PDFFile,0);
   fnlo.PrintFastNLOTableConstants();		// print even more information
   fnlo.PrintTableInfo();			// print some valuable information
   fnlo.SetContributionON(fastNLO::kFixedOrder,0,true); // switch contributions on/off. By default LO and NLO.
   fnlo.SetContributionON(fastNLO::kFixedOrder,1,true);
   fnlo.SetContributionON(fastNLO::kFixedOrder,2,false); // NNLO must be switched on explicitly
   fnlo.SetExternalFuncForMuR( &Function_Mu );	// set function to calculate renormalization scale. See manual for default settings.
   fnlo.SetExternalFuncForMuF( &Function_Mu );
   fnlo.CalcCrossSection();
   fnlo.PrintCrossSections();
   //fnlo.PrintCrossSectionsDefault();
   vector<double > xs = fnlo.GetCrossSection();
   return 0;

}


double Function_Mu(double s1, double s2) {
   // --- fastNLO user: This is an example function
   //     to demonstrate how you might perform the
   //     definition of the scales using a
   //     'flexible-scale'-table, where a function
   //     of s1 and s2 can be used.
   //     Which variables s1 and s2 stand for are
   //     coded in the fastNLO table.
   double mu = 173.;
   return mu;
}

