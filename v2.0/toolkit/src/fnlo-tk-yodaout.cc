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
   cout << endl;
   cout << _CSEPSC << endl;
   shout["fnlo-tk-yodaout"] << "fastNLO YODA Writer" << endl;
   cout << _SSEPSC << endl;
   //! --- fastNLO table and usage info
   string tablename;
   if (argc <= 1) {
      error["fnlo-tk-yodaout"] << "No fastNLO table specified!" << endl;
      shout["fnlo-tk-yodaout"] << "For an explanation of command line arguments type:" << endl;
      shout["fnlo-tk-yodaout"] << "./fnlo-tk-yodaout -h" << endl;
      cout << _CSEPSC << endl;
      exit(1);
   } else {
      tablename = (const char*) argv[1];
      if (tablename == "-h") {
         shout << "" << endl;
         shout << "Usage: ./fnlo-tk-yodaout <fastNLOtable.tab> [PDF] [uncertainty]" << endl;
         shout << "       Arguments: <> mandatory; [] optional." << endl;
         shout << "<fastNLOtable.tab>: Table input file, e.g. fnl2342b.tab" << endl;
         shout << "[PDF]: PDF set, def. = CT10nlo" << endl;
         shout << "   For LHAPDF5: Specify set names WITH filename extension, e.g. \".LHgrid\"." << endl;
         shout << "   For LHAPDF6: Specify set names WITHOUT filename extension." << endl;
         shout << "   If the PDF set still is not found, then:" << endl;
         shout << "   - Check, whether the LHAPDF environment variable is set correctly." << endl;
         shout << "   - Specify the PDF set including the absolute path." << endl;
         shout << "   - Download the desired PDF set from the LHAPDF web site." << endl;
         shout << "[uncertainty]: Uncertainty to show, def. = none" << endl;
         shout << "   Alternatives: NN (none, but correct MC sampling average value --> NNPDF PDFs)" << endl;
         shout << "                 2P (symmetric 2-point scale factor variation)" << endl;
         shout << "                 6P (asymmetric 6-point scale factor variation)" << endl;
         shout << "                 HS (symmetric Hessian PDF uncertainty --> ABM PDFs)" << endl;
         shout << "                 HA (asymmetric Hessian PDF uncertainty)" << endl;
         shout << "                 HP (pairwise asymmetric Hessian PDF uncertainty --> CTEQ|MSTW PDFs)" << endl;
         shout << "                 HC (pairwise asymmetric Hessian PDF uncertainty rescaled to CL68 --> CTEQ PDFs)" << endl;
         shout << "                 MC (MC sampling PDF uncertainty --> NNPDF PDFs)" << endl;
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
         shout << "                 L6 (LHAPDF6 PDF uncertainty --> LHAPDF6 PDFs)" << endl;
#endif
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
   string PDFFile = "X";
   if (argc > 2) {
      PDFFile = (const char*) argv[2];
   }
   if (argc <= 2 || PDFFile == "_") {
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
      PDFFile = "CT10nlo";
#else
      PDFFile = "CT10nlo.LHgrid";
#endif
      shout["fnlo-tk-yodaout"] << "No PDF set given, taking " << PDFFile << " instead!" << endl;
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
      shout["fnlo-tk-yodaout"] << "No request given for uncertainty, none evaluated." << endl;
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
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
      } else if ( chunc == "L6" ) {
         ePDFUnc = kLHAPDF6;
         shout["fnlo-tk-yodaout"] << "Showing LHAPDF6 PDF uncertainty (--> LHAPDF6 PDFs)." << endl;
#endif
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
   fastNLOReader::XsUncertainty XsUnc;
   string LineName;
   if ( chunc == "2P" || chunc == "6P" ) {
      XsUnc = fnlo.GetScaleUncertainty(eScaleUnc);
      snprintf(buffer, sizeof(buffer), " # Relative Scale Uncertainties (%s)",chunc.c_str());
      LineName += "_dxscl";
// #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
//    } else if ( chunc != "L6" ) {
//       vector<LHAPDF::PDFUncertainty> L6PDFUnc = fnlo.GetPDFUncertaintyLHAPDF();
//       xs = fnlo.CalcPDFUncertaintyCentral(L6PDFUnc);
//       vector<double> errdn = fnlo.CalcPDFUncertaintyRelMinus(L6PDFUnc);
//       vector<double> errup = fnlo.CalcPDFUncertaintyRelPlus(L6PDFUnc);
//       for ( unsigned int iobs=0;iobs<xs.size();iobs++ ) {
//          xsdxs.push_back(make_pair(xs[iobs],make_pair(errdn[iobs],errup[iobs])));
//       }
//       snprintf(buffer, sizeof(buffer), " # Relative PDF Uncertainties (%s)",chunc.c_str());
// #endif
   } else if ( chunc != "none" ) {
      XsUnc = fnlo.GetPDFUncertainty(ePDFUnc);
      snprintf(buffer, sizeof(buffer), " # Relative PDF Uncertainties (%s)",chunc.c_str());
      LineName += "_dxpdf";
   }

   if ( XsUnc.xs.size() ) {
      cout << _CSEPSC << endl;
      cout << " # fnlo-tk-yodaout: Evaluating uncertainties" << endl;
      cout << _CSEPSC << endl;
      cout << _DSEPSC << endl;
      cout <<  buffer << endl;
      cout << _SSEPSC << endl;
      cout << " # bin      cross section           lower uncertainty       upper uncertainty" << endl;
      cout << _SSEPSC << endl;
   }
   vector < double > dxsu;
   vector < double > dxsl;
   for ( unsigned int iobs=0;iobs<xs.size();iobs++ ) {
      if ( XsUnc.xs.size() ) {
         printf("%5.i      %#18.11E      %#18.11E      %#18.11E\n",iobs+1,XsUnc.xs[iobs],XsUnc.dxsl[iobs],XsUnc.dxsu[iobs]);
         xs[iobs] = XsUnc.xs[iobs];
         dxsu.push_back(XsUnc.xs[iobs]*XsUnc.dxsu[iobs]);
         dxsl.push_back(XsUnc.xs[iobs]*XsUnc.dxsl[iobs]);
      } else {
         dxsu.push_back(0);
         dxsl.push_back(0);
      }
   }
   if ( XsUnc.xs.size() ) {
      cout << _SSEPSC << endl;
   }

   //! --- Get RivetID
   //!     For 2-dimensions determine running number in Rivet plot name by spotting the capital letter in "RIVET_ID=" in the fnlo table
   size_t capital_pos = 0;
   string RivetId = fnlo.GetRivetId();
   if (RivetId.empty()) {
      error["fnlo-tk-yodaout"] << "No Rivet ID found in fastNLO Table, aborted!" << endl;
      exit(1);
   }
   if ( NDim == 2 ) {
      for (size_t i = fnlo.GetRivetId().find("/"); i < fnlo.GetRivetId().size(); i++) {
         if (isupper(fnlo.GetRivetId()[i])) {                                      // Find capital letter
            RivetId[i] = tolower(RivetId[i]);                                      // and lower it
            capital_pos = i;
            break;
         }
      }
      if (capital_pos == 0) {
         error["fnlo-tk-yodaout"] << "Rivet ID found in fastNLO table does not indicate the 2-dim. histogram counter, aborted." << endl;
         exit(1);
      }
   }

   //! --- Naming the file (and the legend line!) according to calculation order, PDF, and uncertainty choice
   //   string TabName = tablename.substr(0, tablename.size() - 4);
   string PDFName  = PDFFile.substr(0, min(11,(int)PDFFile.size()) );
   string FileName = "NLO_" + PDFName + "_" + chunc;
   LineName = "NLO_" + PDFName + LineName;

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
         eyplus.push_back(dxsu[iobs]);
         eyminus.push_back(abs(dxsl[iobs]));
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
            eyplus.push_back(dxsu[iobs]);
            eyminus.push_back(abs(dxsl[iobs]));
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
