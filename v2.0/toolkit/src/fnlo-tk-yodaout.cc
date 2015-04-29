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
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "fastnlotk/fastNLOLHAPDF.h"
#include "fastnlotk/speaker.h"
#include "YODA/Scatter2D.h"
#include "YODA/WriterYODA.h"


//__________________________________________________________________________________________________________________________________
int main(int argc, char** argv) {

   //! --- namespaces
   using namespace std;
   using namespace say;       //! namespace for 'speaker.h'-verbosity levels
   using namespace fastNLO;   //! namespace for fastNLO constants

   //! --- Parse commmand line
   char buffer[1024];
   cout << _CSEPSC << endl;
   shout["fnlo-tk-yodaout"] << "Program Steering" << endl;
   cout << _SSEPSC << endl;
   //! --- fastNLO table and usage info
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
         shout << "Uncertainty to show, def. = none" << endl;
         shout << "   Alternatives: NN (none, but correct MC sampling average value --> NNPDF PDFs)" << endl;
         shout << "                 2P (symmetric 2-point scale factor variation)" << endl;
         shout << "                 6P (asymmetric 6-point scale factor variation)" << endl;
         shout << "                 HS (symmetric Hessian PDF uncertainty --> ABM PDFs)" << endl;
         shout << "                 HA (asymmetric Hessian PDF uncertainty)" << endl;
         shout << "                 HP (pairwise asymmetric Hessian PDF uncertainty --> CTEQ|MSTW PDFs)" << endl;
         shout << "                 HC (pairwise asymmetric Hessian PDF uncertainty rescaled to CL68 --> CTEQ PDFs)" << endl;
         shout << "                 MC (MC sampling PDF uncertainty --> NNPDF PDFs)" << endl;
         shout << "" << endl;
         shout << "Use \"_\" to skip changing a default argument." << endl;
         shout << "" << endl;
         cout  << _CSEPSC << endl;
         return 0;
      } else {
         shout["fnlo-tk-yodaout"] << "Evaluating table: "  <<  tablename << endl;
      }
   }
   //! --- PDF choice
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
   //! --- Uncertainty choice
   EScaleUncertaintyStyle eScaleUnc = kScaleNone;
   EPDFUncertaintyStyle   ePDFUnc   = kPDFNone;
   string chunc = "none";
   if (argc > 3) {
      chunc = (const char*) argv[3];
   }
   if (argc <= 3 || chunc == "_") {
      info["fnlo-tk-yodaout"] << "No request given for uncertainty, none evaluated." << endl;
   } else {
      if ( chunc == "NN" ) {
         shout["fnlo-tk-yodaout"] << "No uncertainty, but correct MC sampling average value as needed for NNPDF." << endl;
      } else if ( chunc == "2P" ) {
         eScaleUnc = kSymmetricTwoPoint;
         shout["fnlo-tk-yodaout"] << "Showing uncertainty from symmetric 2-point scale factor variation." << endl;
      } else if ( chunc == "6P" ) {
         eScaleUnc = kAsymmetricSixPoint;
         shout["fnlo-tk-yodaout"] << "Showing uncertainty from asymmetric 6-point scale factor variation." << endl;
      } else if ( chunc == "HS" ) {
         ePDFUnc = kHessianSymmetric;
         shout["fnlo-tk-yodaout"] << "Showing symmetric Hessian PDF uncertainty (--> ABM PDFs)." << endl;
      } else if ( chunc == "HA" ) {
         ePDFUnc = kHessianAsymmetric;
         shout["fnlo-tk-yodaout"] << "Showing asymmetric Hessian PDF uncertainty." << endl;
      } else if ( chunc == "HP" ) {
         ePDFUnc = kHessianAsymmetricMax;
         shout["fnlo-tk-yodaout"] << "Showing pairwise asymmetric Hessian PDF uncertainty (--> CTEQ|MSTW PDFs)." << endl;
      } else if ( chunc == "HC" ) {
         ePDFUnc = kHessianCTEQCL68;
         shout["fnlo-tk-yodaout"] << "Showing pairwise asymmetric Hessian PDF uncertainty rescaled to CL68 (--> CTEQ PDFs)." << endl;
      } else if ( chunc == "MC" ) {
         ePDFUnc = kMCSampling;
         shout["fnlo-tk-yodaout"] << "Showing MC sampling PDF uncertainty (--> NNPDF PDFs)." << endl;
      } else {
         error["fnlo-tk-yodaout"] << "Illegal choice of uncertainty, " << chunc << ", aborted!" << endl;
         exit(1);
      }
   }
   cout << _CSEPSC << endl;

   //! --- fastNLO initialisation, read & evaluate table
   //! Select verbosity level
   SetGlobalVerbosity(WARNING);
   //! Initialise a fastNLO instance with interface to LHAPDF
   //! Note: This also initializes the cross section to the LO/NLO one!
   fastNLOLHAPDF fnlo(tablename,PDFFile,0);
   //! Print essential table information
   fnlo.PrintTableInfo();
   // //! Check on existence of non-perturbative corrections from LO MC
   // int inpc1 = fnlo.ContrId(kNonPerturbativeCorrection, kLeading);
   // if (inpc1 > -1) {
   //    info["fnlo-read"] << "Found non-perturbative correction factors. Switch on." << endl;
   //    bool SetOn = fnlo.SetContributionON(kNonPerturbativeCorrection, inpc1, true);
   //    if (!SetOn) {
   //       error["fnlo-tk-yodaout"] << "NPC1 not found, nothing to be done!" << endl;
   //       error["fnlo-tk-yodaout"] << "This should have been caught before!" << endl;
   //       exit(1);
   //    }
   // }

   //! --- Determine dimensioning of observable bins in table
   const int NDim = fnlo.GetNumDiffBin();
   if (NDim < 1 || NDim > 2) {
      error["fnlo-tk-yodaout"] << "Found " << NDim << "-dimensional observable binning in table." << endl;
      error["fnlo-tk-yodaout"] << "Only up to two dimensions currently possible with YODA/Rivet, aborted!" << endl;
      exit(1);
   }

   //! --- Get all required info from table
   //! Get binning
   vector < pair < double, double > > bins = fnlo.GetObsBinsBounds(NDim-1);
   //! Re-calculate cross sections to potentially include the above-selected non-perturbative factors
   fnlo.CalcCrossSection();
   //! Get cross sections
   vector < double > xs = fnlo.GetCrossSection();
   //! If required get uncertainties (only for additive perturbative contributions)
   vector < pair < double, pair < double, double > > > xsdxs;
   vector < pair < double, double > > dxs;
   if ( chunc == "2P" || chunc == "6P" ) {
      xsdxs = fnlo.GetScaleUncertainty(eScaleUnc);
      snprintf(buffer, sizeof(buffer), " # Relative Scale Uncertainties (%s)",chunc.c_str());
   } else if ( chunc != "none" ) {
      xsdxs = fnlo.GetPDFUncertainty(ePDFUnc);
      snprintf(buffer, sizeof(buffer), " # Relative PDF Uncertainties (%s)",chunc.c_str());
   }
   if ( xsdxs.size() ) {
      cout << _CSEPSC << endl;
      cout << " # fnlo-tk-yodaout: Evaluating uncertainties" << endl;
      cout << _CSEPSC << endl;
      cout << _DSEPSC << endl;
      cout <<  buffer << endl;
      cout << _SSEPSC << endl;
      cout << " # bin      cross section           lower uncertainty       upper uncertainty" << endl;
      cout << _SSEPSC << endl;
   }
   for ( unsigned int iobs=0;iobs<xs.size();iobs++ ) {
      if ( xsdxs.size() ) {
         xs[iobs] = xsdxs[iobs].first;
         dxs.push_back(xsdxs[iobs].second);
         printf("%5.i      %#18.11E      %#18.11E      %#18.11E\n",iobs+1,xs[iobs],dxs[iobs].second,dxs[iobs].first);
      } else {
         dxs.push_back(make_pair(0.,0.));
      }
   }
   if ( xsdxs.size() ) {
      cout << _SSEPSC << endl;
   }

   //! --- Get RivetID and running number in Rivet plot name by spotting the capital letter in "RIVET_ID=" in the fnlo table
   size_t capital_pos = 0;
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

   //! --- Naming the file (and the legend line!) according to calculation order, PDF, and uncertainty choice
   string TabName = tablename.substr(0, tablename.size() - 4);
   string PDFName = PDFFile.substr(0, PDFFile.size() - 7);
   string FileName = "NLO_" + PDFName + "_" + chunc;
   string LineName = "NLO-" + PDFName + "_" + "dscale";

   //! --- YODA analysis object creation and storage
   YODA::Writer & writer = YODA::WriterYODA::create();                          //! Create the writer for the yoda file
   vector< YODA::AnalysisObject * > aos;                                   //! Vector of pointers to each of multiple analysis objects
   size_t offset = atoi(RivetId.substr(capital_pos +1, 2).c_str());

   //! --- Initialize dimension bin and continuous observable bin counter
   unsigned int NDimBins[NDim];
   NDimBins[0] = fnlo.GetNDim0Bins();
   unsigned int iobs = 0;
   //! Vector of 2D scatter plots
   vector<YODA::Scatter2D> plots;

   //! --- 1D
   if (NDim == 1) {
      //! Vectors to fill 2D scatter plot
      vector < double > x;
      vector < double > y;
      vector < double > exminus;
      vector < double > explus;
      vector < double > eyminus;
      vector < double > eyplus;
      //! Loop over bins in outer (1st) dimension
      for (unsigned int k =0 ; k<NDimBins[0] ; k++) {
         x.push_back((bins[iobs].second + bins[iobs].first)/2.0);
         explus.push_back((bins[iobs].second - bins[iobs].first)/2.0);
         exminus.push_back((bins[iobs].second - bins[iobs].first)/2.0);
         y.push_back(xs[iobs]);
         eyplus.push_back( xs[iobs]*dxs[iobs].first);
         eyminus.push_back(xs[iobs]*abs(dxs[iobs].second));
         iobs++;
      }
      stringstream plotno;                                                                         // To make i+1 from int
      plotno << offset;                                                                            // to a string for the naming
      RivetId.replace( capital_pos +3 - plotno.str().size(), plotno.str().size(), plotno.str());   // Next plot name
      // Pointer in order not to be deleted after we exit the loop, so we can then save them into the yoda file
      YODA::Scatter2D * plot = new YODA::Scatter2D(x,y,exminus,explus,eyminus,eyplus,"/" + RivetId,LineName);
      // Insert the plot pointer into the vector of analysis object pointers
      aos.push_back(plot);
   }
   //! --- 2D
   else if (NDim == 2) {
      //! Loop over bins in outer (1st) dimension
      for (unsigned int j=0; j<NDimBins[0]; j++) {
         //! Vectors to fill 2D scatter plot
         vector < double > x;
         vector < double > y;
         vector < double > exminus;
         vector < double > explus;
         vector < double > eyminus;
         vector < double > eyplus;
         //! Loop over bins in inner (2nd) dimension
         NDimBins[1] = fnlo.GetNDim1Bins(j);
         for (unsigned int k = 0; k<NDimBins[1]; k++) {
            x.push_back((bins[iobs].second + bins[iobs].first)/2.0);
            explus.push_back((bins[iobs].second - bins[iobs].first)/2.0);
            exminus.push_back((bins[iobs].second - bins[iobs].first)/2.0);
            y.push_back(xs[iobs]);
            eyplus.push_back( xs[iobs]*dxs[iobs].first);
            eyminus.push_back(xs[iobs]*abs(dxs[iobs].second));
            iobs++;
         }
         stringstream plotno;                                                                         // To make i+1 from int
         plotno << offset+j;                                                                          // to a string for the naming
         RivetId.replace( capital_pos +3 - plotno.str().size(), plotno.str().size(), plotno.str());   // Next plot name
         // Pointer in order not to be deleted after we exit the loop, so we can then save them into the yoda file
         YODA::Scatter2D * plot = new YODA::Scatter2D(x,y,exminus,explus,eyminus,eyplus,"/" + RivetId,LineName);
         // Insert the plot pointer into the vector of analysis object pointers
         aos.push_back(plot);
      }
   }

   //! --- Output
   //! Save histograms into the yoda file
   writer.write( FileName + ".yoda", aos );
   cout << endl;
   cout << _CSEPSC << endl;
   shout << "" << endl;
   shout << FileName + ".yoda was successfully produced" << endl;
   shout << "" << endl;
   cout << _CSEPSC << endl << endl;
}
