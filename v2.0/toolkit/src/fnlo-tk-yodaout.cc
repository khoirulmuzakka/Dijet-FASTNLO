//********************************************************************
//
//     fnlo-tk-yodaout
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
#include "fastnlotk/speaker.h"
#include <YODA/Histo1D.h>
#include <YODA/HistoBin1D.h>
#include <YODA/WriterYODA.h>

double Function_Mu(double s1, double s2);

//   ------------------------ use ------------------------
//   ./fnlo-tk-yodaout fnlo_table.tab PDF.LHgrid xmur,xmuf
//   -----------------------------------------------------

//______________________________________________________________________________________________________________
int main(int argc, char** argv) {
   // usage: fastNLO <fastNLO-table.tab> [<LHAPDFfile>]
   using namespace std;
   using namespace say;          // namespace for 'speaker.h'-verbosity levels
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
   string PDFFile = "CT10nlo.LHgrid";
   if (argc > 2)    PDFFile = (const char*) argv[2];

   //--- give some output
   cout<<" fnlo-tk-yodaout: Evaluating table: " << tablename << endl;
   cout<<" fnlo-tk-yodaout: Using PDF set   : " << PDFFile << endl;




   // First we find the point of change of the rapidity value, then for each one we create a histogram which we store in a vector and
   // in the end save it in a yoda file

   float xmur = 1.0;
   float xmuf = 1.0;
   if (argc > 3){                           //getting the scale factors from input in form xmur,xmuf
      string scalefactors = argv[3];
      xmur = (float)::atof(scalefactors.substr( 0, scalefactors.find(",") ).c_str());
      xmuf = (float)::atof(scalefactors.substr( scalefactors.find(",")+1, scalefactors.size()-scalefactors.find(",")-1 ).c_str());
   }
   cout << " fnlo-tk-yodaout: Renormalization and factorization scale factors: " << xmur << " and " << xmuf << "\n";



   fastNLOLHAPDF fnlo(tablename,PDFFile,0);                    // initialize a fastNLO instance with interface to LHAPDF.
   fnlo.PrintTableInfo();                                       // print some valuable information
   if ( xmur != 1.0 || xmuf != 1.0 ){ fnlo.SetScaleFactorsMuRMuF(xmur,xmuf); }     //set the desired scale factors (1.0 is default)
   fnlo.CalcCrossSection();                                    // Calculate the cross section
   //fnlo.PrintCrossSections();                                // Print cross section to screen


   std::vector<YODA::HistoBin1D> bins;                                      // vector that will accept the pT bins
   YODA::Writer & writer = YODA::WriterYODA::create();                        // creat the writer for the yoda file
   std::vector< YODA::AnalysisObject * > ao;                                  // vector that will accept the pointers of the histograms for each rapidity value


   //Get RivetID and running number in histogram name by spotting the capital letter in "RIVETID=" in the fnlotable
   size_t capital_pos = 0;                           //RivetId capital letter where the histogram counter runs
   string RivetId = fnlo.GetRivetId();
   if (RivetId.empty()) {
      error["fnlo-tk-yodaout"] << "No Rivet ID found in fastNLO Table." << endl;
      exit(1);
   }

   for (size_t i = fnlo.GetRivetId().find("/"); i < fnlo.GetRivetId().size(); i++){
      if (isupper(fnlo.GetRivetId()[i])){                  //find capital letter
         RivetId[i] = tolower(RivetId[i]);               //and lower it
         capital_pos = i;
         break;
      }
   }
   if (capital_pos == 0) {
      error["fnlo-tk-yodaout"] << "Rivet ID in fastNLO table is not valid." << endl;
      exit(1);
   }

   //naming the file with PDF and scaling factors
   string PDFName = PDFFile.substr(0, PDFFile.size() - 7);
   stringstream Xmur, Xmuf;
   Xmur << xmur;
   Xmuf << xmuf;
   string FileName = PDFName + "_" + Xmur.str() + "_" + Xmuf.str();



   //histogram creation and storage
   size_t offset = atoi(RivetId.substr(capital_pos +1, 2).c_str());
   for (size_t i=0; i < fnlo.GetNBinDimI(); i++) {                                        // for all rapidity bins
      std::vector<YODA::HistoBin1D> bins;                                               // vector that will accept the pT bins

      stringstream histno;                                                              // just to make i+1 from int
      histno << offset+i;                                                                    // to a string for the naming

      RivetId.replace( capital_pos +3 - histno.str().size(), histno.str().size(), histno.str());         // next histogram name

      for (int k = 0 ; k< fnlo.GetNBinDimII(i) ; k++) {                                  // starting from the first pT bin of each rapidity
         bins.push_back(YODA::HistoBin1D( fnlo.GetBinDimII(i)[k].first , fnlo.GetBinDimII(i)[k].second ) );       // insert pT bin into the vector
      }

      YODA::Histo1D * hist = new YODA::Histo1D(bins, "/" + RivetId, FileName   );               //create histogram pointer
      // pointer in order not to be deleted after we exit the loop, so we can then save them into the yoda file

      for (int k =0 ; k< fnlo.GetNBinDimII(i) ; k++) {                        // fill in the histogram (*length as we fill area)
         hist->fill( (fnlo.GetBinDimII(i)[k].first + fnlo.GetBinDimII(i)[k].second)/2.0 ,
               fnlo.GetCrossSection2Dim()[i][k]*(fnlo.GetBinDimII(i)[k].second - fnlo.GetBinDimII(i)[k].first) );
      }
      ao.push_back(hist);                                                               // insert the histogram pointer into the vector
   }

   writer.write( FileName + ".yoda", ao );                                                  // save histograms into the yoda file
   cout << "###################################################################################" << "\n \n";
   cout << FileName + ".yoda was succesfully produced" << "\n \n";

}
