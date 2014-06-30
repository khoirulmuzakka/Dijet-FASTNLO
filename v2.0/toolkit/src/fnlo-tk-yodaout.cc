//********************************************************************
//
//     fnlo-tk-yodaout
//     This is your playgroud to calculate cross sections
//     with fastNLO.
//
//********************************************************************
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include "fastnlotk/fastNLOLHAPDF.h"
#include <YODA/Histo1D.h>
#include <YODA/HistoBin1D.h>
#include <YODA/WriterYODA.h>

double Function_Mu(double s1, double s2);

//______________________________________________________________________________________________________________
int main(int argc, char** argv) {
   // usage: fastNLO <fastNLO-table.tab> [<LHAPDFfile>]
   using namespace std;
   using namespace fastNLO;      // namespace for fastNLO constants

   // ---  Parse commmand line
   if (argc <= 1) {
      cout<<" Program fnlo-tk-yodaout"<<endl;
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
   cout<<" fnlo-tk-yodaout: Evaluating table: " << tablename << endl;
   cout<<" fnlo-tk-yodaout: Using PDF set   : " << PDFFile << endl;


   // --- this is your playgroud to use fastNLO
   //  Calculate cross setions and/or test some options
   //  For a documentation and function calls, please see
   //    './src/fnlo-cppread.cc'

   //--- example calculation
   fastNLOLHAPDF fnlo(tablename,PDFFile,0);     // initialize a fastNLO instance with interface to LHAPDF.
   //fnlo.PrintTableInfo();                     // print some valuable information
   //fnlo.PrintFastNLOTableConstants();         // print even more information
   //fnlo.SetUnits(kAbsoluteUnits);             // Use units as specified in the publication or in barns.
   //fnlo.SetContributionON(fastNLO::kFixedOrder,0,false); // switch contributions on/off. By default LO and NLO.
   //fnlo.SetContributionON(fastNLO::kFixedOrder,1,true);
   //fnlo.SetContributionON(fastNLO::kFixedOrder,2,true); // NNLO must be switched on explicitly
   fnlo.CalcCrossSection();                     // Calculate the cross section
   fnlo.PrintCrossSections();                 // Print cross section to screen






//--------------------------- yodaout code - start

// First we find the point of change of the rapidity value, then for each one we create a histogram which we store in a vector and
// in the end save it in a yoda file


   std::vector<YODA::HistoBin1D> bins;                             // vector that will accept the pT bins


   YODA::Writer & writer = YODA::WriterYODA::create();                  // creat the writer for the yoda fiel
   std::vector< YODA::AnalysisObject * > ao;                            // vector that will accept the pointers of the histograms for each rapidity value

  

   for (unsigned int i=0; i<fnlo.GetNBinDimI(); i++) {                // for all rapidity bins
	std::vector<YODA::HistoBin1D> bins;                             // vector that will accept the pT bins

        stringstream histno;                                            // just to make i+1 from int
        histno << i+1;                                                  // to a string for the naming

   
 	 for (unsigned int k = 0 ; k< fnlo.GetNBinDimII(i) ; k++) {                // starting from the first pT bin of each rapidity
                bins.push_back(YODA::HistoBin1D( fnlo.GetBinDimII(i)[k].first , fnlo.GetBinDimII(i)[k].second ) );            // insert pT bin into the vector
     	 }



        YODA::Histo1D * hist = new YODA::Histo1D(bins, "/" + fnlo.GetRivetId() + "/d0" + histno.str() + "-x01-y01", "fastNLO" );  
	//create histogram pointer
        // pointer in order not to be deleted after we exit the loop, so we can then save them into the yoda file



     	for (unsigned int k =0 ; k< fnlo.GetNBinDimII(i) ; k++) { 
	      hist->fill( (fnlo.GetBinDimII(i)[k].first + fnlo.GetBinDimII(i)[k].second)/2.0 ,
                           fnlo.GetCrossSection2Dim()[i][k]*(fnlo.GetBinDimII(i)[k].second - fnlo.GetBinDimII(i)[k].first) ); 
									// fill in the histogram (*length as we fill area,not height)
        }


        ao.push_back(hist);                                             // insert the histogram pointer into the vector
   }



   writer.write( "fastNLO_Histograms.yoda", ao );                                // save histograms into the yoda file

   cout << "fastNLO_Histograms.yoda was succesfully produced" << "\n \n";



   //delete ao;                                                         // need to find a way to delete hists from memory now
//-------------------------- yodaout code - end







//    vector<double> cs = fnlo.GetCrossSection(); // Access cross sections for later usage

   // finish:
//   return 0;


   // example code how to loop over all PDF eigenvectors
//   cout<<"\n fnlo-tk-yodaout: Now we want to loop over the eigenvectors of "<<PDFFile<<"."<<endl<<endl;
//   fnlo.SetLHAPDFFilename(PDFFile); // we use again the 'nominal' PDF-file
//   int nEig = fnlo.GetNPDFMembers(); // How many eigenvectors are there?
//   cout<<" fnlo-tk-yodaout: There are "<<nEig<<" Eigenvalue sets in "<<PDFFile<<endl;
//   for ( int i = 0 ; i<nEig ; i++ ) { // start with 0
//      cout<<" fnlo-tk-yodaout: Setting PDF member: "<<i<<" ***"<<endl;
//      fnlo.SetLHAPDFMember(i);  // specify the PDF member
//      fnlo.CalcCrossSection();  // redo cross section calculation
//      fnlo.PrintCrossSections(); // print new cross sections to screen
      // write tot file
//      vector<double> cs = fnlo.GetCrossSection(); // get cross setions for further usage (i.e. printing to file)
//   }

   // finish
 //  return 0;

   // example to use different scale factors
//   cout<<"\n fnlo-tk-yodaout: Now we use a scale factor of 2."<<endl;
//   fnlo.SetScaleFactorsMuRMuF(2.0,2.0);
//   fnlo.CalcCrossSection();
//   fnlo.PrintCrossSections();


//   cout<<"\n fnlo-tk-yodaout: Now we use a scale factor of 0.5."<<endl;
//   fnlo.SetScaleFactorsMuRMuF(0.5,0.5);
//   fnlo.CalcCrossSection();
//   fnlo.PrintCrossSections();


//   cout<<"\n fnlo-tk-yodaout: Now we go back to the nominal result: 1."<<endl;
//   fnlo.SetScaleFactorsMuRMuF(1.0,1.0);
//   fnlo.CalcCrossSection();
 //  fnlo.PrintCrossSections();

   // finish
//   return 0;

//}


//double Function_Mu(double s1, double s2) {
   // --- fastNLO user: This is an example function
   //     to demonstrate how you might perform the
   //     definition of the scales using a
   //     'flexible-scale'-table, where a function
   //     of s1 and s2 can be used.
   //     Which variables s1 and s2 stand for are
   //     coded in the fastNLO table.
//   double mu = 173.;
//   return mu;
}
