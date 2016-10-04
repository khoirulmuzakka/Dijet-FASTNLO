///********************************************************************
///
///     fastNLO_toolkit: fnlo-tk-rootout
///     Program to read fastNLO tables and write out
///     QCD cross sections into ROOT histograms
///
///     K. Rabbertz
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
//! Includes for filling ROOT histograms
//! Usable only when configured with '--with-root=/path/to/root' option
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
//! End of ROOT part


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
   info["fnlo-tk-rootout"] << "Program to read fastNLO tables and write out" << endl;
   info["fnlo-tk-rootout"] << "QCD cross sections into ROOT histograms" << endl;
   yell << _SSEPSC << endl;
   info["fnlo-tk-rootout"] << "For more explanations type:" << endl;
   info["fnlo-tk-rootout"] << "./fnlo-tk-rootout -h" << endl;
   yell << _CSEPSC << endl;

   //! --- Parse commmand line
   char buffer[1024];
   char titlel[1024];
   char titleu[1024];
   yell << "" << endl;
   yell << _CSEPSC << endl;
   shout["fnlo-tk-rootout"] << "fastNLO ROOT Writer" << endl;
   yell << _SSEPSC << endl;
   string tablename;
   if (argc <= 1) {
      error["fnlo-tk-rootout"] << "No fastNLO table specified!" << endl;
      shout["fnlo-tk-rootout"] << "For an explanation of command line arguments type:" << endl;
      shout["fnlo-tk-rootout"] << "./fnlo-tk-rootout -h" << endl;
      yell << _CSEPSC << endl;
      exit(1);
   } else {
      tablename = (const char*) argv[1];
      //! --- Usage info
      if (tablename == "-h") {
         yell << " #" << endl;
         info["fnlo-tk-rootout"] << "This program evaluates a fastNLO table and" << endl;
         info["fnlo-tk-rootout"] << "writes histograms with cross sections and scale or" << endl;
         info["fnlo-tk-rootout"] << "PDF uncertainties into ROOT." << endl;
         info["fnlo-tk-rootout"] << "" << endl;
         info["fnlo-tk-rootout"] << "TODO: Provide more info on histogram numbering/labelling ..." << endl;
         man << "" << endl;
         man << "Usage: ./fnlo-tk-rootout <fastNLOtable.tab> [PDF] [uncertainty]" << endl;
         man << "       Arguments: <> mandatory; [] optional." << endl;
         man << "<fastNLOtable.tab>: Table input file, e.g. fnl2342b.tab" << endl;
         man << "[PDF]: PDF set, def. = series of CT14nlo, MMHT2014nlo68cl, NNPDF30_nlo_as_0118, PDF4LHC15_nlo_mc" << endl;
         man << "   For LHAPDF5: Specify set names WITH filename extension, e.g. \".LHgrid\"." << endl;
         man << "   For LHAPDF6: Specify set names WITHOUT filename extension." << endl;
         man << "   If the PDF set still is not found, then:" << endl;
         man << "   - Check, whether the LHAPDF environment variable is set correctly." << endl;
         man << "   - Specify the PDF set including the absolute path." << endl;
         man << "   - Download the desired PDF set from the LHAPDF web site." << endl;
         man << "[PDF uncertainty]: Uncertainty to show, def. = none" << endl;
         man << "   Alternatives: NN (none, but correct MC sampling average value --> NNPDF PDFs)" << endl;
         man << "                 HS (symmetric Hessian PDF uncertainty --> ABM, (G)JR PDFs)" << endl;
         man << "                 HA (asymmetric Hessian PDF uncertainty)" << endl;
         man << "                 HP (pairwise asymmetric Hessian PDF uncertainty --> CTEQ|MSTW PDFs)" << endl;
         man << "                 HC (pairwise asymmetric Hessian PDF uncertainty rescaled to CL68 --> CTEQ PDFs)" << endl;
         man << "                 MC (MC sampling PDF uncertainty --> NNPDF PDFs)" << endl;
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
         man << "                 L6 (LHAPDF6 PDF uncertainty --> LHAPDF6 PDFs)" << endl;
#endif
         man << "[order]: Fixed-order precision to use, def. = series up to highest fixed-order available" << endl;
         man << "   Alternatives: LO, NLO, NNLO (if available)" << endl;
         man << "[norm]: Normalize if applicable, def. = no." << endl;
         man << "   Alternatives: \"yes\" or \"norm\"" << endl;
         yell << " #" << endl;
         man << "Use \"_\" to skip changing a default argument." << endl;
         yell << " #" << endl;
         yell << _CSEPSC << endl;
         return 0;
      } else {
         shout["fnlo-tk-rootout"] << "Evaluating table: "  <<  tablename << endl;
      }
   }
   //! --- PDF choice
   string PDFFile;
   if (argc > 2) {
      PDFFile = (const char*) argv[2];
   }
   if (argc <= 2 || PDFFile == "_") {
      PDFFile = "X";
      shout["fnlo-tk-rootout"] << "No PDF set given, taking standard series of PDF sets instead!" << endl;
   } else {
      shout["fnlo-tk-rootout"] << "Using PDF set   : " << PDFFile << endl;
   }
   //! --- Uncertainty choice
   EScaleUncertaintyStyle eScaleUnc = kScaleNone;
   EPDFUncertaintyStyle   ePDFUnc   = kPDFNone;
   string chunc;
   if (argc > 3) {
      chunc = (const char*) argv[3];
   }
   if (argc <= 3 || chunc == "_") {
      chunc = "none";
      shout["fnlo-tk-rootout"] << "No request given for PDF uncertainty, none evaluated." << endl;
   } else {
      if ( chunc == "NN" ) {
         shout["fnlo-tk-rootout"] << "No PDF uncertainty, but correct MC sampling average value as needed for NNPDF." << endl;
      } else if ( chunc == "HS" ) {
         ePDFUnc = kHessianSymmetric;
         shout["fnlo-tk-rootout"] << "Showing symmetric Hessian PDF uncertainty (--> ABM, (G)JR PDFs)." << endl;
      } else if ( chunc == "HA" ) {
         ePDFUnc = kHessianAsymmetric;
         shout["fnlo-tk-rootout"] << "Showing asymmetric Hessian PDF uncertainty." << endl;
      } else if ( chunc == "HP" ) {
         ePDFUnc = kHessianAsymmetricMax;
         shout["fnlo-tk-rootout"] << "Showing pairwise asymmetric Hessian PDF uncertainty (--> CTEQ|MSTW PDFs)." << endl;
      } else if ( chunc == "HC" ) {
         ePDFUnc = kHessianCTEQCL68;
         shout["fnlo-tk-rootout"] << "Showing pairwise asymmetric Hessian PDF uncertainty rescaled to CL68 (--> CTEQ PDFs)." << endl;
      } else if ( chunc == "MC" ) {
         ePDFUnc = kMCSampling;
         shout["fnlo-tk-rootout"] << "Showing MC sampling PDF uncertainty (--> NNPDF PDFs)." << endl;
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
      } else if ( chunc == "L6" ) {
         ePDFUnc = kLHAPDF6;
         shout["fnlo-tk-rootout"] << "Showing LHAPDF6 PDF uncertainty (--> LHAPDF6 PDFs)." << endl;
#endif
      } else {
         error["fnlo-tk-rootout"] << "Illegal choice of uncertainty, " << chunc << ", aborted!" << endl;
         exit(1);
      }
   }
   //! --- Fixed-order choice
   ESMOrder eOrder = kLeading;
   string chord;
   if (argc > 4) {
      chord = (const char*) argv[4];
   }
   if (argc <= 4 || chord == "_") {
      chord = "ALL";
      shout["fnlo-tk-rootout"] << "No request given for fixed-order precision, using default series of orders." << endl;
   } else {
      if ( chord == "LO" ) {
         eOrder = kLeading;
         shout["fnlo-tk-rootout"] << "Deriving LO cross sections for comparison." << endl;
      } else if ( chord == "NLO" ) {
         eOrder = kNextToLeading;
         shout["fnlo-tk-rootout"] << "Deriving NLO cross sections for comparison." << endl;
      } else if ( chord == "NNLO" ) {
         eOrder = kNextToNextToLeading;
         shout["fnlo-tk-rootout"] << "Deriving NNLO cross sections for comparison." << endl;
      } else {
         error["fnlo-tk-rootout"] << "Illegal choice of fixed-order precision, " << chord << ", aborted!" << endl;
         exit(1);
      }
   }
   //! --- Normalization
   string chnorm;
   if (argc > 5) {
      chnorm = (const char*) argv[5];
   }
   if (argc <= 5 || chnorm == "_") {
      chnorm = "no";
      shout["fnlo-tk-rootout"] << "Preparing unnormalized cross sections," << endl;
   } else {
      shout["fnlo-tk-rootout"] << "Normalizing cross sections. " << endl;
   }
   //! ---  Too many arguments
   if (argc > 6) {
      error["fnlo-tk-rootout"] << "Too many arguments, aborting!" << endl;
      exit(1);
   }
   yell << _CSEPSC << endl;

   //! --- Prepare loop over PDF sets
   const int nsets = 4;
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   string StandardPDFSets[nsets] = {"CT14nlo", "MMHT2014nlo68cl", "NNPDF30_nlo_as_0118", "PDF4LHC15_nlo_mc"};
#else
   string StandardPDFSets[nsets] = {"CT10nlo.LHgrid", "MSTW2008nlo68cl.LHgrid", "NNPDF21_100.LHgrid", "NNPDF10_100.lHpdf"};
#endif
   EPDFUncertaintyStyle ePDFUncs[nsets] = {kHessianCTEQCL68, kHessianAsymmetricMax, kMCSampling, kMCSampling};

   vector < string > PDFFiles;
   vector < EPDFUncertaintyStyle > PDFUncStyles;
   if ( PDFFile == "X" ) {
      for ( int i=0; i<nsets; i++ ) {
         PDFFiles.push_back(StandardPDFSets[i]);
         PDFUncStyles.push_back(ePDFUncs[i]);
      }
   } else {
      PDFFiles.push_back(PDFFile);
      PDFUncStyles.push_back(ePDFUnc);
   }

   //! --- fastNLO initialisation, attach table
   fastNLOTable table = fastNLOTable(tablename);
   //! Print essential table information
   table.PrintContributionSummary(0);

   //! Initialise a fastNLO reader instance with interface to LHAPDF
   //! Note: This also initializes the cross section to the LO/NLO one!
   fastNLOLHAPDF fnlo(table,PDFFiles[0],0);

   //! Check on existence of LO (Id = -1 if not existing)
   int ilo   = fnlo.ContrId(kFixedOrder, kLeading);
   if (ilo < 0) {
      error["fnlo-tk-rootout"] << "LO not found, aborted!" << endl;
      exit(1);
   } else {
      info["fnlo-tk-rootout"] << "The LO contribution has Id: " << ilo << endl;
      fnlo.SetContributionON(kFixedOrder, ilo, true);
   }
   unsigned int nOrder = 1;
   //! Check on existence of NLO (Id = -1 if not existing)
   int inlo  = fnlo.ContrId(kFixedOrder, kNextToLeading);
   if (inlo < 0) {
      info["fnlo-tk-rootout"] << "No NLO contribution found!" << endl;
      if ( eOrder >= kNextToLeading ) {
         error["fnlo-tk-rootout"] << "Requested NLO not found, aborted!" << endl;
         exit(1);
      }
   } else {
      info["fnlo-tk-rootout"] << "The NLO contribution has Id: " << inlo << endl;
      if ( eOrder >= kNextToLeading ) {
         fnlo.SetContributionON(kFixedOrder, inlo, true);
      } else {
         fnlo.SetContributionON(kFixedOrder, inlo, false);
      }
      nOrder = 2;
   }
   //! Check on existence of NNLO (Id = -1 if not existing)
   int innlo = fnlo.ContrId(kFixedOrder, kNextToNextToLeading);
   if (innlo < 0) {
      info["fnlo-tk-rootout"] << "No NNLO contribution found!" << endl;
      if ( eOrder >= kNextToNextToLeading ) {
         error["fnlo-tk-rootout"] << "Requested NNLO not found, aborted!" << endl;
         exit(1);
      }
   } else {
      info["fnlo-tk-rootout"] << "The NNLO contribution has Id: " << innlo << endl;
      if ( eOrder >= kNextToNextToLeading ) {
         fnlo.SetContributionON(kFixedOrder, innlo, true);
      } else {
         fnlo.SetContributionON(kFixedOrder, innlo, false);
      }
      nOrder = 3;
   }
   // //! Check on existence of non-perturbative corrections from LO MC
   // int inpc1 = fnlo.ContrId(kNonPerturbativeCorrection, kLeading);
   // if (inpc1 > -1) {
   //    info["fnlo-read"] << "Found non-perturbative correction factors. Switch on." << endl;
   //    bool SetOn = fnlo.SetContributionON(kNonPerturbativeCorrection, inpc1, true);
   //    if (!SetOn) {
   //       error["fnlo-tk-rootout"] << "NPC1 not found, nothing to be done!" << endl;
   //       error["fnlo-tk-rootout"] << "This should have been caught before!" << endl;
   //       exit(1);
   //    }
   // }

   //! Normalize?
   bool lNorm = false;
   if ( chnorm == "yes" || chnorm == "norm" ) {
      if ( fnlo.IsNorm() ) {
         lNorm = true;
      } else {
         error["fnlo-read"] << "Normalization requested but not defined for this table, aborted!" << endl;
         exit(1);
      }
   }

   //! --- Get all required info from table
   //! Determine number and dimensioning of observable bins in table
   //   unsigned int NObsBin = fnlo.GetNObsBin();
   const int NDim = fnlo.GetNumDiffBin();
   if (NDim < 1 || NDim > 2) {
      error["fnlo-tk-rootout"] << "Found " << NDim << "-dimensional observable binning in table." << endl;
      error["fnlo-tk-rootout"] << "Only up to two dimensions currently possible, aborted!" << endl;
      exit(1);
   }

   //! Get binning
   vector < pair < double, double > > bins = fnlo.GetObsBinsBounds(NDim-1);

   //! For flex-scale tables:
   //! Possibility to redefine primary scale Q for mu_r and mu_f from the up to two stored scales
   //! Default choice is the first scale via enum 'kScale1'
   if (fnlo.GetIsFlexibleScaleTable()) {
      fnlo.SetMuFFunctionalForm(kScale1);
      fnlo.SetMuRFunctionalForm(kScale1);
      //      fnlo.SetMuFFunctionalForm(kProd);
      //      fnlo.SetMuRFunctionalForm(kProd);
      warn["fnlo-read"] << "The average scale reported in this example as mu1 is derived "
                        << "from only the first scale of this flexible-scale table." << endl
                        << "                        Please check how this table was filled!" << endl;
   }

   //! --- Create ROOT file
   string BaseName     = tablename.substr(0, max(0,(int)tablename.size() -4));
   if ( PDFFile != "X" )    BaseName = BaseName + "_" + PDFFile;
   if ( chunc   != "none" ) BaseName = BaseName + "_" + chunc;
   if ( chnorm  != "no" )   BaseName = BaseName + "_norm";
   string RootFileName = BaseName + ".root";
   TFile *rootfile = new TFile(RootFileName.c_str(),"NEW");
   //  Initialize histogram counter
   unsigned int nHist = 0;
   const unsigned int nMult = 3;

   //! --- Now loop over PDF sets
   for (unsigned int iPDF=0; iPDF<PDFFiles.size(); iPDF++) {

      //! Initialize table with requested PDF set
      fastNLOLHAPDF fnlo(table,PDFFiles[iPDF],0);

      //! Do multiple fixed-order levels unless desired differently
      unsigned int iOrdMin = 0;
      unsigned int iOrdMax = nOrder;
      if ( chord == "ALL" ) {
         iOrdMin = 0;
         // Full result is stored already in iOrder=0
         iOrdMax = nOrder-1;
      } else if ( chord != "ALL" ) {
         iOrdMin = (unsigned int)eOrder+1;
         iOrdMax = min((unsigned int)eOrder+1,nOrder);
      }
      for (unsigned int iOrder = iOrdMin; iOrder <= iOrdMax; iOrder++) {

         //! Starting default: All fixed-order levels on
         if (iOrder == 0 ) {
         } else if (iOrder == 1 ) {
            if ( innlo > -1 ) fnlo.SetContributionON(kFixedOrder, innlo, false);
            if (  inlo > -1 ) fnlo.SetContributionON(kFixedOrder, inlo, false);
         } else if (iOrder == 2) {
            if ( innlo > -1 ) fnlo.SetContributionON(kFixedOrder, innlo, false);
            if (  inlo > -1 ) fnlo.SetContributionON(kFixedOrder, inlo, true);
         } else if (iOrder == 3) {
            if ( innlo > -1 ) fnlo.SetContributionON(kFixedOrder, innlo, true);
            if (  inlo > -1 ) fnlo.SetContributionON(kFixedOrder, inlo, true);
         } else {
            error["fnlo-tk-rootout"] << "Orders beyond " << nOrder << " are not implemented yet, aborted!" << endl;
            exit(1);
         }
         string sOrder = _OrdName[kFixedOrder][nOrder-1];
         if (iOrder > 0) sOrder = _OrdName[kFixedOrder][iOrder-1];

         //! Re-calculate cross sections for new settings
         fnlo.CalcCrossSection();

         //! Do PDF and scale uncertainties
         unsigned int iOffs[3] = {1,6,8};
         for (unsigned int iUnc = 0; iUnc<3; iUnc++) {

            //! Get cross section & uncertainties (only for additive perturbative contributions)
            XsUncertainty XsUnc;

            //! PDF first
            if ( iUnc==0 ) {
               // Back up zeroth member result when no proper uncertainty choice was made (only approx. correct for NNPDF)
               vector < double > xstmp = fnlo.GetCrossSection(lNorm);
               XsUnc = fnlo.GetPDFUncertainty(PDFUncStyles[iPDF], lNorm);
               if ( PDFUncStyles[iPDF] == kPDFNone ) XsUnc.xs = xstmp;
               snprintf(buffer, sizeof(buffer), " # Relative PDF Uncertainties (%s %s)",sOrder.c_str(),PDFFiles[iPDF].c_str());
               snprintf(titlel, sizeof(titlel), "-dsigma_%s/sigma",PDFFiles[iPDF].c_str());
               snprintf(titleu, sizeof(titleu), "+dsigma_%s/sigma",PDFFiles[iPDF].c_str());
            }
            //! 2P scale uncertainties
            else if ( iUnc==1 ) {
               eScaleUnc = kSymmetricTwoPoint;
               XsUnc = fnlo.GetScaleUncertainty(eScaleUnc, lNorm);
               snprintf(buffer, sizeof(buffer), " # 2P Relative Scale Uncertainties (%s %s)",sOrder.c_str(),PDFFiles[iPDF].c_str());
               snprintf(titlel, sizeof(titlel), "-dsigma_2P/sigma");
               snprintf(titleu, sizeof(titleu), "+dsigma_2P/sigma");
            }
            //! 6P scale uncertainties
            else if ( iUnc==2 ) {
               eScaleUnc = kAsymmetricSixPoint;
               XsUnc = fnlo.GetScaleUncertainty(eScaleUnc, lNorm);
               snprintf(buffer, sizeof(buffer), " # 6P Relative Scale Uncertainties (%s %s)",sOrder.c_str(),PDFFiles[iPDF].c_str());
               snprintf(titlel, sizeof(titlel), "-dsigma_6P/sigma");
               snprintf(titleu, sizeof(titleu), "+dsigma_6P/sigma");
            }

            //! Print out of results
            if ( XsUnc.xs.size() ) {
               yell << _CSEPSC << endl;
               shout << "fnlo-tk-rootout: Evaluating uncertainties" << endl;
               yell << _CSEPSC << endl;
               yell << _DSEPSC << endl;
               yell <<  buffer << endl;
               yell << _SSEPSC << endl;
               shout << "bin      cross section           lower uncertainty       upper uncertainty" << endl;
               yell << _SSEPSC << endl;
               for ( unsigned int iobs=0;iobs<XsUnc.xs.size();iobs++ ) {
                  printf("%5.i      %#18.11E      %#18.11E      %#18.11E\n",iobs+1,XsUnc.xs[iobs],XsUnc.dxsl[iobs],XsUnc.dxsu[iobs]);
               }
               yell << _SSEPSC << endl;
            }

            //! --- Initialize dimension bin and continuous observable bin counter
            unsigned int NDimBins[NDim];
            NDimBins[0] = fnlo.GetNDim0Bins();
            unsigned int nLoops      = 0;
            unsigned int nHistBins   = 0;
            //! --- 1D table
            if (NDim == 1) {
               nLoops  = NDim;
            }
            //! --- 2D table
            else if (NDim == 2) {
               nLoops  = NDimBins[0];
            }
            //! --- 3D table not yet implemented
            else {
               error["fnlo-tk-rootout"] << "Found " << NDim << "-dimensional observable binning in table." << endl;
               error["fnlo-tk-rootout"] << "Only up to two dimensions currently possible, aborted!" << endl;
               exit(1);
            }
            //! --- Book ROOT histos
            TH1D *histo[nMult*nLoops];
            //  Initialize observable bin counter
            unsigned int iobs  = 0;

            //! Loop over bins in outer (1st) dimension (might be only one!)
            for (unsigned int j=0; j<nLoops; j++) {

               //! --- 1D table
               if (NDim == 1) {
                  // No. of bins in differential distribution equals total no. of observable bins
                  nHistBins = NDimBins[0];
               }
               //! --- 2D table
               else if (NDim == 2) {
                  // No. of bins in differential distribution equals no. of observable bins in 2nd dimension
                  NDimBins[1] = fnlo.GetNDim1Bins(j);
                  nHistBins = NDimBins[1];
               }
               //! --- 3D table not yet implemented
               else {
                  error["fnlo-tk-rootout"] << "Found " << NDim << "-dimensional observable binning in table." << endl;
                  error["fnlo-tk-rootout"] << "Only up to two dimensions currently possible, aborted!" << endl;
                  exit(1);
               }

               // Arrays to fill ROOT histos
               // N+1 bin borders
               double xbins[nHistBins+1];
               // N+2 ROOT histo bins
               double ycont[nHistBins+2];
               double dylow[nHistBins+2];
               double dyupp[nHistBins+2];
               // Underflow bin
               ycont[0] = 0.;
               dylow[0] = 0.;
               dyupp[0] = 0.;
               for (unsigned int k = 0; k<nHistBins; k++) {
                  xbins[k]   = bins[iobs].first;
                  ycont[k+1] = XsUnc.xs[iobs];
                  dylow[k+1] = XsUnc.dxsl[iobs];
                  dyupp[k+1] = XsUnc.dxsu[iobs];
                  iobs++;
               }
               // Right edge
               xbins[nHistBins] = bins[iobs-1].second;
               // Overflow bin
               ycont[nHistBins+1] = 0.;
               dylow[nHistBins+1] = 0.;
               dyupp[nHistBins+1] = 0.;

               // fastNLO numbering for ROOT histos
               int iScale = 1;
               int iProc  = 0;
               int iDim0  = j+1;
               int iOff   = 0;
               int iHist = iOrder * 1000000 + iScale * 100000 + iProc * 10000 + iDim0 * 100 + iPDF * 10 + iOff;
               char histno[9];

               if ( iUnc==0 ) {
                  // Cross section
                  sprintf(histno,"h%07i",iHist);
                  histo[j+iUnc] = new TH1D(histno,BaseName.c_str(),nHistBins,xbins);
                  histo[j+iUnc]->SetContent(ycont);
                  histo[j+iUnc]->GetXaxis()->SetTitle(fnlo.GetDimLabels()[NDim-1].c_str());
                  histo[j+iUnc]->GetYaxis()->SetTitle(fnlo.GetXSDescr().c_str());
                  nHist++;
               }

               // Lower uncertainty
               iOff  = iOffs[iUnc];
               iHist = iOrder * 1000000 + iScale * 100000 + iProc * 10000 + iDim0 * 100 + iPDF * 10 + iOff;
               sprintf(histno,"h%07i",iHist);
               histo[j+iUnc] = new TH1D(histno,BaseName.c_str(),nHistBins,xbins);
               histo[j+iUnc]->SetContent(dylow);
               histo[j+iUnc]->GetXaxis()->SetTitle(fnlo.GetDimLabels()[NDim-1].c_str());
               histo[j+iUnc]->GetYaxis()->SetTitle(titlel);
               nHist++;

               // Upper uncertainty
               iOff  = iOffs[iUnc]+1;
               iHist = iOrder * 1000000 + iScale * 100000 + iProc * 10000 + iDim0 * 100 + iPDF * 10 + iOff;
               sprintf(histno,"h%07i",iHist);
               histo[j+iUnc] = new TH1D(histno,BaseName.c_str(),nHistBins,xbins);
               histo[j+iUnc]->SetContent(dyupp);
               histo[j+iUnc]->GetXaxis()->SetTitle(fnlo.GetDimLabels()[NDim-1].c_str());
               histo[j+iUnc]->GetYaxis()->SetTitle(titleu);
               nHist++;
            }
         }
      }
   }

   //! --- Output
   //! Save histograms into the ROOT file
   rootfile->cd();
   rootfile->Write();
   rootfile->Close();
   yell << "" << endl;
   yell << _CSEPSC << endl;
   yell << " #" << endl;
   shout << RootFileName + " with " << nHist << " histograms was successfully produced" << endl;
   yell << " #" << endl;
   yell << _CSEPSC << endl << endl;
}
