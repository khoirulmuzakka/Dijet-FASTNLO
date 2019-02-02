///********************************************************************
///
///     fastNLO_toolkit: fnlo-tk-yodaout
///     Program to read fastNLO tables and write out
///     QCD cross sections in YODA format for use with Rivet
///
///     K. Rabbertz, G. Sieber, S. Tyros
///
///********************************************************************

// Precompiler variables for conditional compilation are generated and
// stored automatically in config.h via AC_DEFINE statements in configure.ac.
// To enable conditional compilation, e.g. using HAVE_LIBZ, this config file
// MUST be the very first one to be included with
#include <config.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "fastnlotk/fastNLOAlphas.h"
#include "fastnlotk/fastNLOLHAPDF.h"
#include "fastnlotk/speaker.h"
#ifdef WITH_YODA
#include "YODA/Scatter2D.h"
#include "YODA/WriterYODA.h"
#endif

//__________________________________________________________________________________________________________________________________
int main(int argc, char** argv) {

   //! --- namespaces
   using namespace std;
   using namespace say;       //! namespace for 'speaker.h'-verbosity levels
   using namespace fastNLO;   //! namespace for fastNLO constants

   //! --- Set verbosity level
   SetGlobalVerbosity(INFO);

   //! --- Print program purpose
   yell << _CSEPSC << endl;
   info["fnlo-tk-yodaout"] << "Program to read fastNLO tables and write out" << endl;
   info["fnlo-tk-yodaout"] << "QCD cross sections in YODA format for use with Rivet" << endl;
   info["fnlo-tk-yodaout"] << "(If compiled without YODA support only text printout is given)" << endl;
   yell << _SSEPSC << endl;
   info["fnlo-tk-yodaout"] << "For more explanations type:" << endl;
   info["fnlo-tk-yodaout"] << "./fnlo-tk-yodaout -h" << endl;
   yell << _CSEPSC << endl;

   //! --- Parse commmand line
   char buffer[1024];
   yell << "" << endl;
   yell << _CSEPSC << endl;
   shout["fnlo-tk-yodaout"] << "fastNLO YODA Writer" << endl;
   yell << _SSEPSC << endl;
   string tablename;
   if (argc <= 1) {
      error["fnlo-tk-yodaout"] << "No fastNLO table specified!" << endl;
      shout["fnlo-tk-yodaout"] << "For an explanation of command line arguments type:" << endl;
      shout["fnlo-tk-yodaout"] << "./fnlo-tk-yodaout -h" << endl;
      yell << _CSEPSC << endl;
      exit(1);
   } else {
      tablename = (const char*) argv[1];
      //! --- Usage info
      if (tablename == "-h") {
         yell << " #" << endl;
         info["fnlo-tk-yodaout"] << "This program evaluates a fastNLO table and" << endl;
         info["fnlo-tk-yodaout"] << "prints out cross sections with either scale or" << endl;
         info["fnlo-tk-yodaout"] << "PDF uncertainties in YODA format for use with Rivet." << endl;
         info["fnlo-tk-yodaout"] << "For this to work, the scenario description must contain" << endl;
         info["fnlo-tk-yodaout"] << "the Rivet ID in the form 'RIVET_ID=EXP_YYYY_INSPIREID/Dii-xjj-ykk'," << endl;
         info["fnlo-tk-yodaout"] << "where 'ii', 'jj', and 'kk' indicate the first histogram covered by" << endl;
         info["fnlo-tk-yodaout"] << "this table and the capital letter indicates the plot counter to increase." << endl;
         info["fnlo-tk-yodaout"] << "In case the Rivet ID is missing, it can be added e.g. by using 'fnlo-tk-modify'." << endl;
         man << "" << endl;
         man << "Usage: ./fnlo-tk-yodaout <fastNLOtable.tab> [PDF] [uncertainty]" << endl;
         man << "       Arguments: <> mandatory; [] optional." << endl;
         man << "<fastNLOtable.tab>: Table input file, e.g. fnl2342b.tab" << endl;
         man << "[PDF]: PDF set, def. = CT10nlo" << endl;
         man << "   For LHAPDF5: Specify set names WITH filename extension, e.g. \".LHgrid\"." << endl;
         man << "   For LHAPDF6: Specify set names WITHOUT filename extension." << endl;
         man << "   If the PDF set still is not found, then:" << endl;
         man << "   - Check, whether the LHAPDF environment variable is set correctly." << endl;
         man << "   - Specify the PDF set including the absolute path." << endl;
         man << "   - Download the desired PDF set from the LHAPDF web site." << endl;
         man << "[uncertainty]: Uncertainty to show, def. = none" << endl;
         man << "   Alternatives: NN (none, but correct MC sampling average value --> NNPDF PDFs)" << endl;
         man << "                 2P (symmetric 2-point scale factor variation)" << endl;
         man << "                 6P (asymmetric 6-point scale factor variation)" << endl;
         man << "                 HS (symmetric Hessian PDF uncertainty --> ABM, (G)JR PDFs)" << endl;
         man << "                 HA (asymmetric Hessian PDF uncertainty)" << endl;
         man << "                 HP (pairwise asymmetric Hessian PDF uncertainty --> CTEQ|MSTW PDFs)" << endl;
         man << "                 HC (pairwise asymmetric Hessian PDF uncertainty rescaled to CL68 --> CTEQ PDFs)" << endl;
         man << "                 MC (MC sampling PDF uncertainty --> NNPDF PDFs)" << endl;
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
         man << "                 L6 (LHAPDF6 PDF uncertainty --> LHAPDF6 PDFs)" << endl;
#endif
         man << "                 AS (a_s(M_Z) variation uncertainty with GRV evolution)" << endl;
         man << "[order]: Fixed-order precision to use, def. = NLO" << endl;
         man << "   Alternatives: LO, NLO_only, NNLO, NNLO_only (if available)" << endl;
         man << "[norm]: Normalize if applicable, def. = no." << endl;
         man << "   Alternatives: \"yes\" or \"norm\"" << endl;
         man << "[flexscale]: Central scale choice for flex-scale tables." << endl;
         man << "   Default:      \"scale1\",  i.e. mur=muf=scale1," << endl;
         man << "   Alternatives: \"scale2\",  i.e. mur=muf=scale2," << endl;
         man << "                 \"scale12\", i.e. mur=scale1, muf=scale2," << endl;
         man << "                 \"scale21\", i.e. mur=scale2, muf=scale1." << endl;
         man << "[np]: Apply nonperturbative corrections if available, def. = no." << endl;
         man << "   Alternatives: \"yes\" or \"np\"" << endl;
         yell << " #" << endl;
         man << "Use \"_\" to skip changing a default argument." << endl;
         yell << " #" << endl;
         yell << _CSEPSC << endl;
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
   EAsUncertaintyStyle    eAsUnc    = kAsNone;
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
         shout["fnlo-tk-yodaout"] << "Showing symmetric Hessian PDF uncertainty (--> ABM, (G)JR PDFs)." << endl;
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
      } else if ( chunc == "AS" ) {
         eAsUnc = kAsGRV;
         shout["fnlo-tk-yodaout"] << "Showing a_s(M_Z) uncertainty with GRV evolution." << endl;
      } else {
         error["fnlo-tk-yodaout"] << "Illegal choice of uncertainty, " << chunc << ", aborted!" << endl;
         exit(1);
      }
   }
   //! --- Fixed-order choice
   ESMOrder eOrder = kNextToLeading;
   string chord = "NLO";
   bool bonly = false;
   if (argc > 4) {
      chord = (const char*) argv[4];
   }
   if (argc <= 4 || chord == "_") {
      shout["fnlo-tk-yodaout"] << "No request given for fixed-order precision, using NLO." << endl;
   } else {
      if ( chord == "LO" ) {
         eOrder = kLeading;
         shout["fnlo-tk-yodaout"] << "Deriving LO cross sections for comparison." << endl;
      } else if ( chord == "NLO" ) {
         eOrder = kNextToLeading;
         shout["fnlo-tk-yodaout"] << "Deriving NLO cross sections for comparison." << endl;
      } else if ( chord == "NLO_only" ) {
         eOrder = kNextToLeading;
         bonly = true;
         shout["fnlo-tk-yodaout"] << "Experimental: Deriving NLO contribution only for comparison." << endl;
      } else if ( chord == "NNLO" ) {
         eOrder = kNextToNextToLeading;
         shout["fnlo-tk-yodaout"] << "Deriving NNLO cross sections for comparison." << endl;
      } else if ( chord == "NNLO_only" ) {
         eOrder = kNextToNextToLeading;
         shout["fnlo-tk-yodaout"] << "Experimental: Deriving NNLO contribution only for comparison." << endl;
      } else {
         error["fnlo-tk-yodaout"] << "Illegal choice of fixed-order precision, " << chord << ", aborted!" << endl;
         exit(1);
      }
   }

   //! --- Normalization
   string chnorm = "no";
   if (argc > 5) {
      chnorm = (const char*) argv[5];
   }
   if (argc <= 5 || chnorm == "_") {
      shout["fnlo-tk-yodaout"] << "Preparing unnormalized cross sections," << endl;
   } else {
      shout["fnlo-tk-yodaout"] << "Normalizing cross sections. " << endl;
   }

   //--- Scale choice (flex-scale tables only; ignored for fix-scale tables)
   string chflex = "scale1";
   if (argc > 6) {
      chflex = (const char*) argv[6];
   }
   if (argc <= 6 || chflex == "_") {
      shout["fnlo-tk-yodaout"] << "Using default mur=muf=scale 1." << endl;
   } else {
      shout["fnlo-tk-yodaout"] << "Using scale definition "+chflex << endl;
   }

   //! --- Nonperturbative correction
   string chnp = "no";
   if (argc > 7) {
      chnp = (const char*) argv[7];
   }
   if (argc <= 7 || chnp == "_") {
      shout["fnlo-tk-yodaout"] << "Do not apply nonperturbative corrections." << endl;
   } else {
      shout["fnlo-tk-yodaout"] << "Apply nonperturbative corrections if available." << endl;
   }

   //! ---  Too many arguments
   if (argc > 8) {
      error["fnlo-tk-yodaout"] << "Too many arguments, aborting!" << endl;
      exit(1);
   }
   yell << _CSEPSC << endl;

   //! --- fastNLO initialisation, read & evaluate table
   //! Initialise a fastNLO instance with interface to LHAPDF
   //! Note: This also initializes the cross section to the LO/NLO one!
   fastNLOLHAPDF* fnlo = NULL;
   if ( chunc != "AS" ) {
      fnlo = new fastNLOLHAPDF(tablename,PDFFile,0);
   } else {
      fnlo = new fastNLOAlphas(tablename,PDFFile,0);
   }

   //! Print essential table information
   fnlo->PrintContributionSummary(0);

   //! Check on existence of LO (Id = -1 if not existing)
   int ilo   = fnlo->ContrId(kFixedOrder, kLeading);
   if (ilo < 0) {
      if ( bonly ) {
         warn["fnlo-tk-yodaout"] << "LO not found, output for NLO_only or NNLO_only is experimental ...!" << endl;
      } else {
         error["fnlo-tk-yodaout"] << "LO not found; full cross section can not be derived, aborted!" << endl;
         exit(1);
      }
   } else {
      info["fnlo-tk-yodaout"] << "The LO contribution has Id: " << ilo << endl;
      if ( (eOrder > kLeading) & bonly ) {
         fnlo->SetContributionON(kFixedOrder, ilo, false);
      } else {
         fnlo->SetContributionON(kFixedOrder, ilo, true);
      }
   }
   //! Check on existence of NLO (Id = -1 if not existing)
   int inlo  = fnlo->ContrId(kFixedOrder, kNextToLeading);
   if (inlo < 0) {
      info["fnlo-tk-yodaout"] << "No NLO contribution found!" << endl;
      if ( eOrder == kNextToLeading || ((eOrder > kNextToLeading) & !bonly) ) {
         error["fnlo-tk-yodaout"] << "Requested NLO not found, aborted!" << endl;
         exit(1);
      }
   } else {
      info["fnlo-tk-yodaout"] << "The NLO contribution has Id: " << inlo << endl;
      if ( eOrder == kNextToLeading || ((eOrder > kNextToLeading) & !bonly) ) {
         fnlo->SetContributionON(kFixedOrder, inlo, true);
      } else {
         fnlo->SetContributionON(kFixedOrder, inlo, false);
      }
   }
   //! Check on existence of NNLO (Id = -1 if not existing)
   int innlo = fnlo->ContrId(kFixedOrder, kNextToNextToLeading);
   if (innlo < 0) {
      info["fnlo-tk-yodaout"] << "No NNLO contribution found!" << endl;
      if ( eOrder == kNextToNextToLeading || ((eOrder > kNextToNextToLeading) & !bonly) ) {
         error["fnlo-tk-yodaout"] << "Requested NNLO not found, aborted!" << endl;
         exit(1);
      }
   } else {
      info["fnlo-tk-yodaout"] << "The NNLO contribution has Id: " << innlo << endl;
      if ( eOrder == kNextToNextToLeading || ((eOrder > kNextToNextToLeading) & !bonly) ) {
         fnlo->SetContributionON(kFixedOrder, innlo, true);
      } else {
         fnlo->SetContributionON(kFixedOrder, innlo, false);
      }
   }
   //! Check on existence of non-perturbative corrections from LO MC
   int inpc1 = fnlo->ContrId(kNonPerturbativeCorrection, kLeading);
   if (inpc1 > -1 && (chnp == "yes" || chnp == "np" ) ) {
      info["fnlo-tk-yodaout"] << "Found non-perturbative correction factors. Switch on." << endl;
      bool SetOn = fnlo->SetContributionON(kNonPerturbativeCorrection, inpc1, true);
      if (!SetOn) {
         error["fnlo-tk-yodaout"] << "NPC1 not found, nothing to be done!" << endl;
         error["fnlo-tk-yodaout"] << "This should have been caught before!" << endl;
         exit(1);
      }
   }

   //! Normalize?
   bool lNorm = false;
   if ( chnorm == "yes" || chnorm == "norm" ) {
      if ( fnlo->IsNorm() ) {
         lNorm = true;
      } else {
         error["fnlo-read"] << "Normalization requested but not defined for this table, aborted!" << endl;
         exit(1);
      }
   }

   //! --- Determine dimensioning of observable bins in table
   const int NDim = fnlo->GetNumDiffBin();
   if (NDim < 1 || NDim > 2) {
      error["fnlo-tk-yodaout"] << "Found " << NDim << "-dimensional observable binning in table." << endl;
      error["fnlo-tk-yodaout"] << "Only up to two dimensions currently possible with YODA/Rivet, aborted!" << endl;
      exit(1);
   }

   //! --- Get all required info from table
   //! Get binning
   vector < pair < double, double > > bins = fnlo->GetObsBinsBounds(NDim-1);

   //! For flex-scale tables:
   //! Possibility to redefine primary scale Q for mu_r and mu_f from the up to two stored scales
   //! Default choice is the first scale via enum 'kScale1'
   if (fnlo->GetIsFlexibleScaleTable()) {
      if ( chflex == "scale1" ) {
         fnlo->SetMuFFunctionalForm(kScale1);
         fnlo->SetMuRFunctionalForm(kScale1);
         info["fnlo-tk-yodaout"] << "The average scale reported in this example as mu1 is derived "
                                 << "from only the first scale of this flexible-scale table." << endl
                                 << "                        Please check how this table was filled!" << endl;
      } else if ( chflex == "scale2" ) {
         fnlo->SetMuFFunctionalForm(kScale2);
         fnlo->SetMuRFunctionalForm(kScale2);
         info["fnlo-tk-yodaout"] << "The average scale reported in this example as mu2 is derived "
                                 << "from only the second scale of this flexible-scale table." << endl
                                 << "                        Please check how this table was filled!" << endl;
      } else if ( chflex == "scale12" ) {
         fnlo->SetMuFFunctionalForm(kScale2);
         fnlo->SetMuRFunctionalForm(kScale1);
         info["fnlo-tk-yodaout"] << "The average scale reported in this example as mu1 is derived "
                                 << "from only the first scale of this flexible-scale table." << endl
                                 << "                        Please check how this table was filled!" << endl;
      } else if ( chflex == "scale21" ) {
         fnlo->SetMuFFunctionalForm(kScale1);
         fnlo->SetMuRFunctionalForm(kScale2);
         info["fnlo-tk-yodaout"] << "The average scale reported in this example as mu2 is derived "
                                 << "from only the second scale of this flexible-scale table." << endl
                                 << "                        Please check how this table was filled!" << endl;
      } else {
         error["fnlo-tk-yodaout"] << "Unknown scale choice " << chflex << ", aborted!" << endl;
      }
   }

   //! Re-calculate cross sections to potentially include the above-selected non-perturbative factors
   fnlo->CalcCrossSection();
   //! Get cross sections
   vector < double > xs = fnlo->GetCrossSection(lNorm);
   //! If required get uncertainties (only for additive perturbative contributions)
   XsUncertainty XsUnc;
   string LineName;
   if ( chunc == "2P" || chunc == "6P" ) {
      XsUnc = fnlo->GetScaleUncertainty(eScaleUnc, lNorm);
      snprintf(buffer, sizeof(buffer), " # Relative Scale Uncertainties (%s)",chunc.c_str());
      LineName += "_dxscl";
   } else if ( chunc == "AS" ) {
      XsUnc = fnlo->GetAsUncertainty(eAsUnc, lNorm);
      snprintf(buffer, sizeof(buffer), " # Relative a_s(M_Z) Uncertainties (%s)",chunc.c_str());
      LineName += "_dxa_s";
   } else if ( chunc != "none" ) {
      XsUnc = fnlo->GetPDFUncertainty(ePDFUnc, lNorm);
      snprintf(buffer, sizeof(buffer), " # Relative PDF Uncertainties (%s)",chunc.c_str());
      LineName += "_dxpdf";
   }

   if ( XsUnc.xs.size() ) {
      yell << _CSEPSC << endl;
      shout << "fnlo-tk-yodaout: Evaluating uncertainties" << endl;
      yell << _CSEPSC << endl;
      yell << _DSEPSC << endl;
      yell <<  buffer << endl;
      yell << _SSEPSC << endl;
      shout << "bin      cross_section           lower_uncertainty       upper_uncertainty" << endl;
      yell << _SSEPSC << endl;
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
      yell << _SSEPSC << endl;
   }

   //! --- Get RivetID
   //!     For 2-dimensions determine running number in Rivet plot name by spotting the capital letter in "RIVET_ID=" in the fnlo table
   //!     For inverted order, the Id starts with "-" after "/".
   size_t capital_pos  = 0;
   int    invert_order = 1;
   string RivetId = fnlo->GetRivetId();
   if (RivetId.empty()) {
      warn["fnlo-tk-yodaout"] << "No Rivet ID found in fastNLO Table, no YODA formatted output possible, exiting!" << endl;
      exit(0);
   }
   if ( NDim == 2 ) {
      size_t i0 = RivetId.find("/");
      /// Check for inversion of histogram order
      if ( RivetId.substr(i0,2) == "/-" ) {
         invert_order = -1;
         RivetId.erase(i0+1,1);
      }
      for (size_t i = i0; i < RivetId.size(); i++) {
         /// Identify capital letter in "d01-x01-y01" type strings, memorize position and lower it
         if (isupper(RivetId[i])) {
            RivetId[i] = tolower(RivetId[i]);
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
   //   string FileName = chord + "_" + PDFName + "_" + chunc + "_" + chnorm + "_" + chnp;
   string FileName = PDFName;
   if ( chord  != "_" ) {FileName = chord + "_" + PDFName;}
   if ( chunc  != "_" ) {FileName += "_" + chunc;}
   if ( chnorm != "_" ) {FileName += "_" + chnorm;}
   if ( chnp   != "_" ) {FileName += "_" + chnp;}
   LineName = chord + "_" + PDFName + LineName;

   //! --- YODA analysis object creation and storage
#ifdef WITH_YODA
   YODA::Writer & writer = YODA::WriterYODA::create();                          //! Create the writer for the yoda file
   vector< YODA::AnalysisObject * > aos;                                   //! Vector of pointers to each of multiple analysis objects
#endif
   size_t counter = atoi(RivetId.substr(capital_pos +1, 2).c_str());

   //! --- Initialize dimension bin and continuous observable bin counter
   unsigned int NDimBins[NDim];
   NDimBins[0] = fnlo->GetNDim0Bins();
   unsigned int iobs = 0;
   //! Vector of 2D scatter plots
#ifdef WITH_YODA
   vector<YODA::Scatter2D> plots;
#endif

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
#ifdef WITH_YODA
      /// Pointer in order not to be deleted after we exit the loop, so we can then save them into the yoda file
      YODA::Scatter2D * plot = new YODA::Scatter2D(x,y,exminus,explus,eyminus,eyplus,"/" + RivetId,LineName);
      /// Insert the plot pointer into the vector of analysis object pointers
      aos.push_back(plot);
#endif
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
         NDimBins[1] = fnlo->GetNDim1Bins(j);
         for (unsigned int k = 0; k<NDimBins[1]; k++) {
            x.push_back((bins[iobs].second + bins[iobs].first)/2.0);
            explus.push_back((bins[iobs].second - bins[iobs].first)/2.0);
            exminus.push_back((bins[iobs].second - bins[iobs].first)/2.0);
            y.push_back(xs[iobs]);
            eyplus.push_back(dxsu[iobs]);
            eyminus.push_back(abs(dxsl[iobs]));
            iobs++;
         }
         /// Derive histogram counter
         size_t ihist = (invert_order > 0) ? (counter + j) : (counter - j);
         if ( ihist == 0 || ihist > 99 ) {
            error["fnlo-tk-yodaout"] << "Rivet histogram counter out of range, aborted!" << endl;
            error["fnlo-tk-yodaout"] << "ihist = " << ihist << endl;
            exit(1);
         }
         /// Convert size_t into string for naming
         stringstream histno;
         histno << ihist;
         /// Replace counter part in RivetId by histno
         RivetId.replace(capital_pos +3 - histno.str().size(), histno.str().size(), histno.str());
#ifdef WITH_YODA
         /// Pointer in order not to be deleted after we exit the loop, so we can then save the plots into the yoda file
         YODA::Scatter2D * plot = new YODA::Scatter2D(x,y,exminus,explus,eyminus,eyplus,"/" + RivetId,LineName);
         /// Insert the plot pointer into the vector of analysis object pointers
         aos.push_back(plot);
#endif
      }
   }

   //! --- Output
   //! Save histograms into the yoda file
#ifdef WITH_YODA
   writer.write( FileName + ".yoda", aos );
#endif
   yell << "" << endl;
   yell << _CSEPSC << endl;
   yell << " #" << endl;
#ifdef WITH_YODA
   shout << FileName + ".yoda was successfully produced" << endl;
#else
   shout << "Compiled without YODA support ==> No YODA file produced!" << endl;
#endif
   yell << " #" << endl;
   yell << _CSEPSC << endl << endl;
}
