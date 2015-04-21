///********************************************************************
///
///     fnlo-tk-example
///     This is your playgroud to calculate cross sections
///     with fastNLO.
///
///********************************************************************

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "fastnlotk/fastNLOLHAPDF.h"

double Function_Mu(double s1, double s2);

//______________________________________________________________________________________________________________
int main(int argc, char** argv) {
   //! usage: fastNLO <fastNLO-table.tab> [<LHAPDFfile>]
   using namespace std;
   using namespace fastNLO;      //! namespace for fastNLO constants

   //! ---  Parse commmand line
   if (argc <= 1) {
      cout<<" Program fnlo-tk-example"<<endl;
      cout<<"   usage:"<<endl;
      cout<<"   fastNLO <fastNLO-table.tab> [<LHAPDF-file>]"<<endl;
      return 0;
   }
   //! --- fastNLO table
   string tablename =  (const char*) argv[1];
   //---  PDF set
   string PDFFile = "CT10nnlo.LHgrid";
   if (argc > 2)    PDFFile = (const char*) argv[2];

   //--- give some output
   cout<<" fnlo-tk-example: Evaluating table: " << tablename << endl;
   cout<<" fnlo-tk-example: Using PDF set   : " << PDFFile << endl;


   //! --- this is your playgroud to use fastNLO
   //!  Calculate cross setions and/or test some options
   //!  For a documentation and function calls, please see
   //!    './src/fnlo-cppread.cc'

   //--- example calculation
   fastNLOLHAPDF fnlo(tablename,PDFFile,0);     //! initialize a fastNLO instance with interface to LHAPDF.
   fnlo.PrintTableInfo();                       //! print some valuable information
   //fnlo.PrintFastNLOTableConstants();         //! print even more information
   //fnlo.SetUnits(kAbsoluteUnits);             //! Use units as specified in the publication or in barns.
   //fnlo.SetContributionON(fastNLO::kFixedOrder,0,false); //! switch contributions on/off. By default LO and NLO.
   //fnlo.SetContributionON(fastNLO::kFixedOrder,1,true);
   //fnlo.SetContributionON(fastNLO::kFixedOrder,2,true); //! NNLO must be switched on explicitly
   fnlo.CalcCrossSection();                     //! Calculate the cross section
   fnlo.PrintCrossSections();                   //! Print cross section to screen


   vector<double> cs = fnlo.GetCrossSection();  //! Access cross sections for later usage

   //! finish:
   return 0;


   //! example code how to loop over all PDF eigenvectors
   cout<<"\n fnlo-tk-example: Now we want to loop over the eigenvectors of "<<PDFFile<<"."<<endl<<endl;
   fnlo.SetLHAPDFFilename(PDFFile); //! we use again the 'nominal' PDF-file
   int nEig = fnlo.GetNPDFMembers(); //! How many eigenvectors are there?
   cout<<" fnlo-tk-example: There are "<<nEig<<" Eigenvalue sets in "<<PDFFile<<endl;
   for ( int i = 0 ; i<nEig ; i++ ) { //! start with 0
      cout<<" fnlo-tk-example: Setting PDF member: "<<i<<" ***"<<endl;
      fnlo.SetLHAPDFMember(i);  //! specify the PDF member
      fnlo.CalcCrossSection();  //! redo cross section calculation
      fnlo.PrintCrossSections(); //! print new cross sections to screen
      //! write tot file
      vector<double> cs = fnlo.GetCrossSection(); //! get cross setions for further usage (i.e. printing to file)
   }

   //! finish
   return 0;

   //! example to use different scale factors
   cout<<"\n fnlo-tk-example: Now we use a scale factor of 2."<<endl;
   fnlo.SetScaleFactorsMuRMuF(2.0,2.0);
   fnlo.CalcCrossSection();
   fnlo.PrintCrossSections();


   cout<<"\n fnlo-tk-example: Now we use a scale factor of 0.5."<<endl;
   fnlo.SetScaleFactorsMuRMuF(0.5,0.5);
   fnlo.CalcCrossSection();
   fnlo.PrintCrossSections();


   cout<<"\n fnlo-tk-example: Now we go back to the nominal result: 1."<<endl;
   fnlo.SetScaleFactorsMuRMuF(1.0,1.0);
   fnlo.CalcCrossSection();
   fnlo.PrintCrossSections();

   //! finish
   return 0;

}


double Function_Mu(double s1, double s2) {
   //! --- fastNLO user: This is an example function
   //!     to demonstrate how you might perform the
   //!     definition of the scales using a
   //!     'flexible-scale'-table, where a function
   //!     of s1 and s2 can be used.
   //!     Which variables s1 and s2 stand for are
   //!     coded in the fastNLO table.
   double mu = 173.;
   return mu;
}
