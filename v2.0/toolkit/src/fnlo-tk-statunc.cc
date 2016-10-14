///********************************************************************
///
///     fastNLO_reader: fnlo-tk-statunc
///     Program to read a sample of fastNLO tables and write out
///     averaged QCD cross sections with statistical uncertainties
///
///     For more explanations type:
///     ./fnlo-tk-statunc -h
///
///     K. Rabbertz
///
///********************************************************************

#include <config.h>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
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
   info["fnlo-tk-statunc"] << "Program to read a sample of fastNLO tables and write out" << endl;
   info["fnlo-tk-statunc"] << "averaged QCD cross sections with statistical uncertainties" << endl;
   yell << _SSEPSC << endl;
   info["fnlo-tk-statunc"] << "For more explanations type:" << endl;
   info["fnlo-tk-statunc"] << "./fnlo-tk-statunc -h" << endl;
   yell << _CSEPSC << endl;

   //! --- Parse commmand line
   char buffer[1024]; // TODO: Use PATH_MAX instead?
   yell << "" << endl;
   yell << _CSEPSC << endl;
   shout["fnlo-tk-statunc"] << "fastNLO Statistical Table Sample Evaluator"<<endl;
   yell << _SSEPSC << endl;
   string tablebase;
   if (argc <= 1) {
      error["fnlo-tk-statunc"] << "No fastNLO sample specified!" << endl;
      shout["fnlo-tk-statunc"] << "For an explanation of command line arguments type:" << endl;
      shout["fnlo-tk-statunc"] << "./fnlo-tk-statunc -h" << endl;
      yell << _CSEPSC << endl;
      exit(1);
   } else {
      tablebase = (const char*) argv[1];
      //! --- Usage info
      if (tablebase == "-h") {
         yell << " #" << endl;
         info["fnlo-tk-statunc"] << "This program evaluates a sample of equivalent, but" << endl;
         info["fnlo-tk-statunc"] << "statistically independent fastNLO tables and" << endl;
         info["fnlo-tk-statunc"] << "provides their averaged cross sections together with" << endl;
         info["fnlo-tk-statunc"] << "statistical uncertainties." << endl;
         info["fnlo-tk-statunc"] << "In addition, some statistical sample characteristics are" << endl;
         info["fnlo-tk-statunc"] << "printed for each observable bin. The tables files with the" << endl;
         info["fnlo-tk-statunc"] << "largest downwards and upwards fluctuations from the average are indicated." << endl;
         info["fnlo-tk-statunc"] << "All files with any deviation larger than 10 sigma are listed" << endl;
         info["fnlo-tk-statunc"] << "in an extra kill file for removal by the user." << endl;
         info["fnlo-tk-statunc"] << "If optional YODA support is configured and proper RIVET_ID's are" << endl;
         info["fnlo-tk-statunc"] << "available in the fastNLO tables, the result is additionally provided" << endl;
         info["fnlo-tk-statunc"] << "in YODA format." << endl;
         man << "" << endl;
         man << "Usage: ./fnlo-tk-statunc <fastNLOsample> [PDF] [order] [nmin] [nmax]" << endl;
         man << "       Arguments: <> mandatory; [] optional." << endl;
         man << "<fastNLOsample>: Basename of table input files, e.g. fnl2452_I1082936_v23_flex-hhc-born-2jet," << endl;
         man << "                 that will be complemented by '_nnnn.tab'" << endl;
         man << "[PDF]: PDF set, def. = CT10nlo" << endl;
         man << "   For LHAPDF5: Specify set names WITH filename extension, e.g. \".LHgrid\"." << endl;
         man << "   For LHAPDF6: Specify set names WITHOUT filename extension." << endl;
         man << "   If the PDF set still is not found, then:" << endl;
         man << "   - Check, whether the LHAPDF environment variable is set correctly." << endl;
         man << "   - Specify the PDF set including the absolute path." << endl;
         man << "   - Download the desired PDF set from the LHAPDF web site." << endl;
         man << "[order]: Fixed-order precision to use, def. = NLO" << endl;
         man << "   Alternatives: LO, NNLO (if available)" << endl;
         man << "[nmin]: Smallest table number nnnn to start with, def. = 0000." << endl;
         man << "[nmax]: Largest  table number nnnn to end with, def. = 1000." << endl;
         yell << " #" << endl;
         man << "Use \"_\" to skip changing a default argument." << endl;
         yell << " #" << endl;
         yell << _CSEPSC << endl;
         return 0;
      } else {
         shout["fnlo-tk-statunc"] << "Evaluating table sample: "  <<  tablebase << endl;
      }
   }
   //---  PDF set
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
      shout["fnlo-tk-statunc"] << "No PDF set given, taking " << PDFFile << " instead!" << endl;
   } else {
      shout["fnlo-tk-statunc"] << "Using PDF set   : " << PDFFile << endl;
   }
   //! --- Fixed-order choice
   ESMOrder eOrder = kNextToLeading;
   string chord = "NLO";
   if (argc > 3) {
      chord = (const char*) argv[3];
   }
   if (argc <= 3 || chord == "_") {
      shout["fnlo-tk-statunc"] << "No request given for fixed-order precision, using NLO." << endl;
   } else {
      if ( chord == "LO" ) {
         eOrder = kLeading;
         shout["fnlo-tk-statunc"] << "Deriving LO cross sections for comparison." << endl;
      } else if ( chord == "NLO" ) {
         eOrder = kNextToLeading;
         shout["fnlo-tk-statunc"] << "Deriving NLO cross sections for comparison." << endl;
      } else if ( chord == "NNLO" ) {
         eOrder = kNextToNextToLeading;
         shout["fnlo-tk-statunc"] << "Deriving NNLO cross sections for comparison." << endl;
      } else {
         error["fnlo-tk-statunc"] << "Illegal choice of fixed-order precision, " << chord << ", aborted!" << endl;
         exit(1);
      }
   }
   //! --- Start counter
   string chmin = "0000";
   if (argc > 4) {
      chmin = (const char*) argv[4];
   }
   if (argc <= 4 || chmin == "_") {
      chmin = "0000";
      shout["fnlo-tk-statunc"] << "No request given for minimal counter number, using default of: " << chmin << endl;
   }
   int nmin = atoi(chmin.c_str());
   if ( nmin < 0 || nmin > 9999 ) {
      error["fnlo-tk-statunc"] << "Illegal minimal counter number for looping over table sample, " << nmin << ", aborted!" << endl;
      exit(1);
   }
   //! --- End counter
   string chmax = "1000";
   if (argc > 5) {
      chmax = (const char*) argv[4];
   }
   if (argc <= 5 || chmax == "_") {
      chmax = "1000";
      shout["fnlo-tk-statunc"] << "No request given for maximal counter number, using default of: " << chmax << endl;
   }
   int nmax = atoi(chmax.c_str());
   if ( nmax < 0 || nmax > 9999 ) {
      error["fnlo-tk-statunc"] << "Illegal maximal counter number for looping over table sample, " << nmax << ", aborted!" << endl;
      exit(1);
   } else if ( nmax < nmin ) {
      error["fnlo-tk-statunc"] << "Maximal counter number smaller than minimal one, nmin = " << nmin << ", nmax = " << nmax << ", aborted!" << endl;
      exit(1);
   }
   yell << _CSEPSC << endl;

   //! --- Reset verbosity level to warning only from here on
   SetGlobalVerbosity(WARNING);

   //! --- Loop over selected table sample
   //! Initialise fastNLO instances with interface to LHAPDF
   //! Note: This also initializes the cross section to the LO/NLO one!
   unsigned int nfound = 0;
   unsigned int nfail  = 0;
   vector < int > iTabMin;
   vector < int > iTabMax;
   vector < double > xsMin;
   vector < double > xsMax;
   // Use xs0 as assumed mean and add up shifted means ,(xsi - xs0), and
   // shifted means squared, (xsi - xs0)^2. Shifting permits numerically
   // more precise variance calculations.
   vector < double > xs0;
   vector < double > xs;
   vector < double > dxs;
   vector < double > dxstmp3;
   vector < double > dxstmp4;
   string firsttable;

   for ( int itab=nmin; itab<=nmax; itab++) {
      char buftmp[5];
      snprintf(buftmp, sizeof(buftmp), "%04d", itab);
      string tablename = tablebase + "_" + buftmp + ".tab";
      if ( access(tablename.c_str(), R_OK) == 0 ) {
         nfound++;
         fastNLOLHAPDF fnlo(tablename,PDFFile,0);
         //! Print essential table information, but only for very first table
         //! and do some initialization
         if ( nfound == 1 ) {
            firsttable = tablename;
            fnlo.PrintContributionSummary(0);
            xs0 = fnlo.GetCrossSection();
            for (unsigned int i=0; i<xs0.size(); i++) {
               iTabMin.push_back(-1);
               iTabMax.push_back(-1);
               xsMin.push_back(+DBL_MAX);
               xsMax.push_back(-DBL_MAX);
               xs.push_back(0);
               dxs.push_back(0);
            }
         }
         if ( nfound == 1 || (nfound-1) % 10 == 0 ) {shout << "Analyzing table " << tablename << " ..." << endl;}
         //! Check on existence of LO (Id = -1 if not existing)
         int ilo   = fnlo.ContrId(kFixedOrder, kLeading);
         if (ilo < 0) {
            error["fnlo-tk-statunc"] << "LO not found, aborted!" << endl;
            exit(1);
         } else {
            if ( nfound == 1 ) {info["fnlo-tk-statunc"] << "The LO contribution has Id: " << ilo << endl;}
            fnlo.SetContributionON(kFixedOrder, ilo, true);
         }
         //! Check on existence of NLO (Id = -1 if not existing)
         int inlo  = fnlo.ContrId(kFixedOrder, kNextToLeading);
         if (inlo < 0) {
            info["fnlo-tk-statunc"] << "No NLO contribution found!" << endl;
            if ( eOrder >= kNextToLeading ) {
               error["fnlo-tk-statunc"] << "Requested NLO not found, aborted!" << endl;
               exit(1);
            }
         } else {
            if ( nfound == 1 ) {info["fnlo-tk-statunc"] << "The NLO contribution has Id: " << inlo << endl;}
            if ( eOrder >= kNextToLeading ) {
               fnlo.SetContributionON(kFixedOrder, inlo, true);
            } else {
               fnlo.SetContributionON(kFixedOrder, inlo, false);
            }
         }
         //! Check on existence of NNLO (Id = -1 if not existing)
         int innlo = fnlo.ContrId(kFixedOrder, kNextToNextToLeading);
         if (innlo < 0) {
            if ( nfound == 1 ) {info["fnlo-tk-statunc"] << "No NNLO contribution found!" << endl;}
            if ( eOrder >= kNextToNextToLeading ) {
               error["fnlo-tk-statunc"] << "Requested NNLO not found, aborted!" << endl;
               exit(1);
            }
         } else {
            info["fnlo-tk-statunc"] << "The NNLO contribution has Id: " << innlo << endl;
            if ( eOrder >= kNextToNextToLeading ) {
               fnlo.SetContributionON(kFixedOrder, innlo, true);
            } else {
               fnlo.SetContributionON(kFixedOrder, innlo, false);
            }
         }

         //! For flex-scale tables:
         //! Possibility to redefine primary scale Q for mu_r and mu_f from the up to two stored scales
         //! Default choice is the first scale via enum 'kScale1'
         if (fnlo.GetIsFlexibleScaleTable()) {
            fnlo.SetMuFFunctionalForm(kScale1);
            fnlo.SetMuRFunctionalForm(kScale1);
            //      fnlo.SetMuFFunctionalForm(kProd);
            //      fnlo.SetMuRFunctionalForm(kProd);
            if ( nfound == 1 ) {
               warn["fnlo-tk-statunc"] << "The average scale reported in this example as mu1 is derived from " << endl;
               shout << "                     only the first scale of this flexible-scale table." << endl;
               shout << "                     Please check how this table was filled!" << endl;
            }
         }

         //! Re-calculate cross sections for selected fixed order
         fnlo.CalcCrossSection();
         //! Get cross sections
         vector < double > xsTab = fnlo.GetCrossSection();
         for (unsigned int i=0; i<xsTab.size(); i++) {
            if (xsTab[i] < xsMin[i]) {
               iTabMin[i] = itab;
               xsMin[i]   = xsTab[i];
            }
            if (xsTab[i] > xsMax[i]) {
               iTabMax[i] = itab;
               xsMax[i]   = xsTab[i];
            }
            double diff = xsTab[i]-xs0[i];
            xs[i]  += diff;
            dxs[i] += diff*diff;
         }
      } else {
         nfail++;
         //! State that missing tables are ignored at first occurrence
         if ( nfail == 1 ) {warn["fnlo-tk-statunc"] << "Table file " << tablename << " not found, skipped!" << endl;}
      }
   }
   if ( nfound < 2 ) {
      error["fnlo-tk-statunc"] << "Not enough tables found, aborted! nfound = " << nfound << endl;
      exit(1);
   }

   //! Write out statistical fluctuation info
   yell << _CSEPSC << endl;
   shout << "Special info on statistical fluctuations:" << endl;
   yell << _SSEPSC << endl;
   printf(" # No. of filenames skipped:  %5i\n", nfail);
   printf(" # No. of filenames analyzed: %5i\n", nfound);
   yell << _CSEPSC << endl;

   yell << " =========" << _DSEPS << _DSEP40 << endl;
   printf(" #IBIN #IMIN #IMAX         <s>           s_min           s_max       ds/<s>/%% ds_min/<s>/%% ds_max/<s>/%%     ds_min/ds    ds_max/ds\n");
   yell << " ---------" << _SSEPS << _SSEP40 << endl;
   for (unsigned int i=0; i<xs0.size(); i++) {
      double sumw  = xs[i];
      double sumw2 = dxs[i];
      double var   = (sumw2 - pow(sumw,2)/nfound) / (nfound-1);
      // Factor to get agreement with criterium in Fortran code
      // Which version is correct?
      double fac   = sqrt(nfound/(nfound-1.0));
      // Undo shift of mean
      xs[i]  = xs[i]/nfound + xs0[i];
      // Mittlerer Fehler des Mittelwerts
      dxs[i] = +sqrt( var / nfound );
      dxstmp3.push_back((xsMin[i]-xs[i])/sqrt(var)*fac);
      dxstmp4.push_back((xsMax[i]-xs[i])/sqrt(var)*fac);
      printf("%6i%6i%6i%16.5e%16.5e%16.5e%10.3f   %10.3f   %10.3f   %10.3f   %10.3f   \n",
             i+1,iTabMin[i],iTabMax[i],xs[i],xsMin[i],xsMax[i],dxs[i]/xs[i]*100,
             (xsMin[i]/xs[i]-1)*100,(xsMax[i]/xs[i]-1)*100,dxstmp3[i],dxstmp4[i]);
   }

   //! Write out list of files (full path) with largest fluctuations
   string killfile = tablebase + "_kill.txt";
   for ( int itab=nmin; itab<=nmax; itab++) {
      for (unsigned int i=0; i<xs0.size(); i++) {
         if ( (itab == iTabMin[i] && dxstmp3[i] < -10 ) ||
              (itab == iTabMax[i] && dxstmp4[i] > +10 ) ) {
            char *cwd = getcwd( buffer, 1024 );
            char buftmp[5];
            snprintf(buftmp, sizeof(buftmp), "%04d", itab);
            // Stringify first object so everything is +-concatenated to a string and not a char ...
            string tablename = string(cwd) + "/" + tablebase + "_" + buftmp + ".tab";
            ofstream outfile;
            outfile.open(killfile.c_str(),ios::app);
            outfile << tablename << endl;
            outfile.close();
            break; // --> Next table. Need to remove same table only once ...
         }
      }
   }

   // //! Write out cross sections with statistical uncertainty
   // yell << " " << _DSEPS << endl;
   // yell << " Statistical uncertainties of the " << chord << " cross section for the primary scale choice.\n";
   // yell << " " << _SSEPS << endl;
   // yell << "  bin       cross section           rel. average error of the mean\n";
   // yell << " " << _SSEPS << endl;
   // for (unsigned int i=0; i<xs0.size(); i++) {
   //    printf("%5i      %18.11e      %18.11e\n",i+1,xs[i],dxs[i]/xs[i]);
   // }

   //! If YODA is linked, write out cross sections with statistical uncertainty in YODA format as well
   fastNLOLHAPDF fnlo(firsttable,PDFFile,0);

   //! --- Determine dimensioning of observable bins in table
   const int NDim = fnlo.GetNumDiffBin();
   if (NDim < 1 || NDim > 2) {
      error["fnlo-tk-statunc"] << "Found " << NDim << "-dimensional observable binning in table." << endl;
      error["fnlo-tk-statunc"] << "Only up to two dimensions currently possible with YODA/Rivet, aborted!" << endl;
      exit(1);
   }

   //! --- Get all required info from table
   //! Get binning
   vector < pair < double, double > > bins = fnlo.GetObsBinsBounds(NDim-1);

   yell  << _CSEPSC << endl;
   shout << "fnlo-tk-statunc: Evaluating uncertainties" << endl;
   yell  << _CSEPSC << endl;
   yell  << _DSEPSC << endl;
   shout << "Relative Statistical Uncertainties" << endl;
   yell  << _SSEPSC << endl;
   shout << "bin      cross section           lower uncertainty       upper uncertainty" << endl;
   yell  << _SSEPSC << endl;

   vector < double > dxsu;
   vector < double > dxsl;
   for ( unsigned int iobs=0;iobs<xs.size();iobs++ ) {
      printf("%5.i      %#18.11E      %#18.11E      %#18.11E\n",iobs+1,xs[iobs],-dxs[iobs]/xs[iobs],dxs[iobs]/xs[iobs]);
      dxsu.push_back(dxs[iobs]);
      dxsl.push_back(dxs[iobs]);
   }
   yell << _SSEPSC << endl;

   //! Without YODA we can stop here
#ifndef WITH_YODA
   return 0;
#else
   //! --- Get RivetID
   //!     For 2-dimensions determine running number in Rivet plot name by spotting the capital letter in "RIVET_ID=" in the fnlo table
   size_t capital_pos = 0;
   string RivetId = fnlo.GetRivetId();
   if (RivetId.empty()) {
      warn["fnlo-tk-statunc"] << "No Rivet ID found in fastNLO Table, stopping here without YODA output!" << endl;
      exit(0);
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
         error["fnlo-tk-statunc"] << "Rivet ID found in fastNLO table does not indicate the 2-dim. histogram counter, aborted." << endl;
         exit(1);
      }
   }

   //! --- Naming the file (and the legend line!) according to calculation order and uncertainty choice
   string PDFName  = PDFFile.substr(0, min(11,(int)PDFFile.size()) );
   string FileName = chord + "_" + PDFName + "_stat";
   string LineName = FileName + "_dxst";

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
      //      RivetId.replace( capital_pos +3 - plotno.str().size(), plotno.str().size(), plotno.str());   // Next plot name
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
   yell  << "" << endl;
   yell  << _CSEPSC << endl;
   yell << " #" << endl;
   shout << FileName + ".yoda was successfully produced" << endl;
   yell << " #" << endl;
   yell  << _CSEPSC << endl << endl;
#endif
}
