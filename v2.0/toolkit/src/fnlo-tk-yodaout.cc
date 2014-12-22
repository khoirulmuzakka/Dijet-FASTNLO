///********************************************************************
///
///     fastNLO_reader: FNLOYODAOUT
///     Program to read fastNLO v2 tables and write out
///     QCD cross sections in YODA format for use with Rivet
///
///     K. Rabbertz, G. Sieber, S. Tyros
///
///********************************************************************
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

/// Function prototype for flexible-scale function
double Function_Mu(double s1, double s2);

//__________________________________________________________________________________________________________________________________
int main(int argc, char** argv) {

   // namespaces
   using namespace std;
   using namespace say;          // namespace for 'speaker.h'-verbosity levels
   using namespace fastNLO;      // namespace for fastNLO constants

   // ---  Parse commmand line
   cout << _CSEPSC << endl;
   shout["fnlo-tk-yodaout"] << "Program Steering" << endl;
   cout << _SSEPSC << endl;
   // --- fastNLO table
   string tablename;
   if (argc <= 1) {
      error["fnlo-tk-yodaout"] << "No table name given!" << endl;
      shout["fnlo-tk-yodaout"] << "For an explanation of command line arguments type:" << endl;
      shout["fnlo-tk-yodaout"] << "./fnlo-tk-yodaout -h" << endl;
      exit(1);
   } else {
      tablename = (const char*) argv[1];
      if (tablename == "-h") {
         shout << "" << endl;
         shout << "Usage: ./fnlo-tk-yodaout [arguments]" << endl;
         shout << "Table input file, mandatory, e.g. fnl2342b.tab" << endl;
         shout << "PDF set, def. = CT10nlo.LHgrid" << endl;
         shout << "   For LHAPDF5: Give full path(s), if the PDF file is not in the cwd." << endl;
         shout << "   For LHAPDF6: Drop filename extensions and give PDF directory instead." << endl;
         shout << "Values of MuR,MuF scale factors to investigate, if possible, def. = 1.0,1.0" << endl;
         shout << "" << endl;
         shout << "Use \"_\" to skip changing a default argument." << endl;
         shout << "" << endl;
         cout  << _CSEPSC << endl;
         return 0;
      } else {
         shout["fnlo-tk-yodaout"] << "Evaluating table: "  <<  tablename << endl;
      }
   }
   // --- PDF set
   string PDFFile = "CT10nlo.LHgrid";
   if (argc > 2) {
      PDFFile = (const char*) argv[2];
   }
   if (argc <= 2 || PDFFile == "_") {
      PDFFile = "CT10nlo.LHgrid";
      warn ["fnlo-tk-yodaout"] << "No PDF set given," << endl;
      shout["fnlo-tk-yodaout"] << "taking CT10nlo.LHgrid instead!" << endl;
   } else {
      shout["fnlo-tk-yodaout"] << "Using PDF set   : " << PDFFile << endl;
   }
   // --- (MuR,MuF) scale factors
   double xmur = 1.0;
   double xmuf = 1.0;
   string chtmp;
   if (argc > 3) {
      chtmp = (const char*) argv[3];
   }
   if (argc <= 3 || chtmp == "_") {
      warn ["fnlo-tk-yodaout"] << "No request given for scale factors," << endl;
      shout["fnlo-tk-yodaout"] << "investigating default factors of MuR,MuF = 1.0,1.0" << endl;
   } else {
      xmur = (double)::atof(chtmp.substr( 0, chtmp.find(",") ).c_str());
      xmuf = (double)::atof(chtmp.substr( chtmp.find(",")+1, chtmp.size()-chtmp.find(",")-1 ).c_str());
      shout["fnlo-tk-yodaout"] << "Using scale factors of MuR,MuF = " << xmur << "," << xmuf << endl;
   }
   cout << _CSEPSC << endl;

   // --- fastNLO initialisation, read table
   fastNLOLHAPDF fnlo(tablename,PDFFile,0);                                     // Initialise a fastNLO instance with interface to LHAPDF.
   fnlo.PrintTableInfo();                                                       // Print some table information
   if ( xmur != 1.0 || xmuf != 1.0 ){ fnlo.SetScaleFactorsMuRMuF(xmur,xmuf);}   // Set the desired scale factors (1.0 is default)
   fnlo.CalcCrossSection();                                                     // Calculate the cross section
   //fnlo.PrintCrossSections();                                                 // Print cross section to screen

   // --- Get RivetID and running number in histogram name by spotting the capital letter in "RIVETID=" in the fnlo table
   size_t capital_pos = 0;                                                      // RivetId capital letter where the histogram counter runs
   string RivetId = fnlo.GetRivetId();
   if (RivetId.empty()) {
      error["fnlo-tk-yodaout"] << "No Rivet ID found in fastNLO Table, aborted!" << endl;
      exit(1);
   }
   for (size_t i = fnlo.GetRivetId().find("/"); i < fnlo.GetRivetId().size(); i++) {
      if (isupper(fnlo.GetRivetId()[i])) {                                      // Find capital letter
         RivetId[i] = tolower(RivetId[i]);                                      // and lower it
         capital_pos = i;
         break;
      }
   }
   if (capital_pos == 0) {
      error["fnlo-tk-yodaout"] << "Rivet ID found in fastNLO table does not indicate the histogram counter, aborted." << endl;
      exit(1);
   }

   // --- Naming the file with PDF and scaling factors
   string PDFName = PDFFile.substr(0, PDFFile.size() - 7);
   stringstream Xmur, Xmuf;
   Xmur << xmur;
   Xmuf << xmuf;
   string FileName = PDFName + "_" + Xmur.str() + "_" + Xmuf.str();

   // --- Histogram creation and storage
   YODA::Writer & writer = YODA::WriterYODA::create();                          // Create the writer for the yoda file
   std::vector< YODA::AnalysisObject * > ao;                                    // Vector of pointers to each of multiple histograms

   // --- Determine binning/dimensioning of observable bins in table
   const int NDim = fnlo.GetNumDiffBin();
   int NDimBins[NDim];
   size_t offset = atoi(RivetId.substr(capital_pos +1, 2).c_str());

   // --- 1D
   if (NDim == 1) {
      std::vector<YODA::HistoBin1D> bins;                                          // Vector for 1D histogram binning
      stringstream histno;                                                                              // To make i+1 from int
      histno << offset+(NDim-1);                                                                        // to a string for the naming
      RivetId.replace( capital_pos +3 - histno.str().size(), histno.str().size(), histno.str());        // Next histogram name
      NDimBins[0] = fnlo.GetNDim0Bins();
      for (int k = 0; k<NDimBins[0]; k++) {
         bins.push_back(YODA::HistoBin1D( fnlo.GetDim0BinBoundaries()[k].first , fnlo.GetDim0BinBoundaries()[k].second ) ); // Insert bin into the vector
      }
      YODA::Histo1D * hist = new YODA::Histo1D(bins, "/" + RivetId, FileName);                          // Create histogram pointer
      for (int k =0 ; k<NDimBins[0] ; k++) {                                                            // Fill in the histogram (* bin size factor as we fill area)
         hist->fill( (fnlo.GetDim0BinBoundaries()[k].first + fnlo.GetDim0BinBoundaries()[k].second)/2.0 ,
                     fnlo.GetCrossSection()[k]*(fnlo.GetDim0BinBoundaries()[k].second - fnlo.GetDim0BinBoundaries()[k].first) );
         //                     fnlo.GetCrossSection()[k]*fnlo.GetBinSize(k) );
      }
      ao.push_back(hist);                                                                                               // insert the histogram pointer into the vector
   } else if (NDim == 2) {
      unsigned int iobs = 0;
      for (size_t i=0; i < fnlo.GetNDim0Bins(); i++) {                                               // For all bins in outer (first) dimension
         std::vector<YODA::HistoBin1D> bins;                                                        // Vector for 1D histogram binning
         stringstream histno;                                                                       // To make i+1 from int
         histno << offset+i;                                                                        // to a string for the naming
         RivetId.replace( capital_pos +3 - histno.str().size(), histno.str().size(), histno.str()); // Next histogram name
         NDimBins[1] = fnlo.GetNDim1Bins(i);
         cout << "nbindim1 = " << fnlo.GetNDim0Bins() << ", nbindim2 = " << fnlo.GetNDim1Bins(i) << endl;
         for (int k = 0; k<NDimBins[1]; k++) {                                                      // Starting from the first bin in outer (first) dimension
            bins.push_back(YODA::HistoBin1D( fnlo.GetDim1BinBoundaries(i)[k].first , fnlo.GetDim1BinBoundaries(i)[k].second ) ); // Insert bin into the vector
         }
         // Pointer in order not to be deleted after we exit the loop, so we can then save them into the yoda file
         YODA::Histo1D * hist = new YODA::Histo1D(bins, "/" + RivetId, FileName);                              // Create histogram pointer
         for (int k = 0; k<NDimBins[1]; k++) {                                                      // Fill in the histogram (* bin size factor as we fill area)
            iobs++;
            cout << "iobs = " << iobs << ", first = " << fnlo.GetDim1BinBoundaries(i)[k].first << ", second = " << fnlo.GetDim1BinBoundaries(i)[k].second << endl;
            hist->fill( (fnlo.GetDim1BinBoundaries(i)[k].first + fnlo.GetDim1BinBoundaries(i)[k].second)/2.0 ,
                        fnlo.GetCrossSection2Dim()[i][k]*(fnlo.GetDim1BinBoundaries(i)[k].second - fnlo.GetDim1BinBoundaries(i)[k].first) );
            //                        fnlo.GetCrossSection2Dim()[i][k]*fnlo.GetBinSize(iobs) );
         }
         ao.push_back(hist);                                                                        // Insert the histogram pointer into the vector
      }
   } else {
      error["fnlo-tk-yodaout"] << "More than 2 dimensions easily possible in table, but not yet implemented here, aborted!" << endl;
      exit(1);
   }

   // --- Output
   writer.write( FileName + ".yoda", ao );                                                          // Save histograms into the yoda file
   cout << endl;
   cout << _CSEPSC << endl;
   shout << "" << endl;
   shout << FileName + ".yoda was succesfully produced" << endl;
   shout << "" << endl;
   cout << _CSEPSC << endl << endl;
}
