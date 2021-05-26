// Author: Daniel Britzger
// DESY, 20/04/2012

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO_reader_2.1.0                                                //
//  D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch           //
//                                                                      //
//  The projects web page can be found at:                              //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//  If you use this code, please cite:                                  //
//    T. Kluge, K. Rabbertz and M. Wobisch, hep-ph/0609285              //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch,        //
//       arXiv:1109.1310                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//
//  fastNLOAlphas
//  This class inherits the PDF interface from
//  fastNLOLHAPDF, while the alpha_s evolution
//  is superseeded by the Alphas.h class.
//
//////////////////////////////////////////////////////////////////////////

#include <cstring>
#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include "fastnlotk/fastNLOLHAPDF.h"
#include "fastnlotk/fastNLOQCDNUMAS.h"

using namespace std;



//______________________________________________________________________________
//
fastNLOQCDNUMAS::fastNLOQCDNUMAS(std::string name) : fastNLOLHAPDF(name) {
   // Without PDF info use PDG values as default
   SetPDGValues();
   // Print out values for checking
   //   PrintParmValues();
};
fastNLOQCDNUMAS::fastNLOQCDNUMAS(std::string name, std::string LHAPDFFile, int PDFMem) : fastNLOLHAPDF(name,LHAPDFFile,PDFMem) {
   // Set initial values via LHAPDF6 info system
   SetLHAPDFValues(LHAPDFFile, PDFMem);
   // Print out values for checking
   //   PrintParmValues();
};



// Getters
double fastNLOQCDNUMAS::GetQMass(int pdgid) const {
   if (pdgid < 1 || pdgid > 6 ) {
      logger.error["fastNLOQCDNUMAS::GetQMass"]<<"PDG code out of quark index range 1-6! Aborted.\n";
      exit(1);
   }
   return QMass[pdgid];
}
double fastNLOQCDNUMAS::GetMz() const {
   return fMz;
}
std::string fastNLOQCDNUMAS::GetNScheme() const {
   return fnScheme;
}
int fastNLOQCDNUMAS::GetNFlavor(int nflavor) const {
   return nflavor;
}
int fastNLOQCDNUMAS::GetNLoop() const {
   return fnLoop;
}
double fastNLOQCDNUMAS::GetAlphasMz() const {
   return fAlphasMz;
};



// Setters
void fastNLOQCDNUMAS::SetQMass(int pdgid, double qmass) {
   QMass[pdgid] = qmass;
}
void fastNLOQCDNUMAS::SetMz(double Mz) {
   fMz = Mz;
}
void fastNLOQCDNUMAS::SetNFlavor(int  nflavor) {
   fnFlavor = nflavor;
}
void fastNLOQCDNUMAS::SetNLoop(int nloop) {
   if ( nloop < 1 || nloop > 3 ) {
      logger.error["fastNLOQCDNUMAS::SetNLoop"] << "Illegal no. of loops nloop = " << nloop <<
         ", aborted! Only 1, 2, or 3 are allowed with QCDNUM." << endl;
      exit(11);
   }
   fnLoop = nloop;
}
void fastNLOQCDNUMAS::SetAlphasMz(double AlphasMz) {
   fAlphasMz    = AlphasMz;
}



// Combined Setters
void fastNLOQCDNUMAS::SetPDGValues() {
   // Initialize with PDG values
   QMass[0]  = PDG_MD;
   QMass[1]  = PDG_MU;
   QMass[2]  = PDG_MS;
   QMass[3]  = PDG_MC;
   QMass[4]  = PDG_MB;
   QMass[5]  = PDG_MT;
   fMz       = PDG_MZ;
   // Variable flavor number scheme
   fnFlavor = 0;
   // 2-loop alpha_s evolution
   fnLoop = 2;
   fAlphasMz = PDG_ASMZ;
}

void fastNLOQCDNUMAS::SetLHAPDFValues(std::string LHAPDFFile, int PDFMem) {
   // AlphaS_MZ can vary among PDF members, so we really need the PDF member info from LHAPDF
   const LHAPDF::PDFInfo PDFMemInfo(LHAPDFFile, PDFMem);
   QMass[0] = PDFMemInfo.get_entry_as<double>("MDown");
   QMass[1] = PDFMemInfo.get_entry_as<double>("MUp");
   QMass[2] = PDFMemInfo.get_entry_as<double>("MStrange");
   QMass[3] = PDFMemInfo.get_entry_as<double>("MCharm");
   QMass[4] = PDFMemInfo.get_entry_as<double>("MBottom");
   QMass[5] = PDFMemInfo.get_entry_as<double>("MTop");
   fMz      = PDFMemInfo.get_entry_as<double>("MZ");
   fnScheme = PDFMemInfo.get_entry_as<std::string>("FlavorScheme");
   if ( PDFMemInfo.has_key("AlphaS_NumFlavors") ) {
      fnFlavor = PDFMemInfo.get_entry_as<int>("AlphaS_NumFlavors");
   } else {
      fnFlavor = PDFMemInfo.get_entry_as<int>("NumFlavors");
   }
   // Variable flavor numbers are usually set via Nf = 0 in evolution code.
   // Ensure that fnFlavor is maximum Nf for variable flavor number scheme by
   // setting quark masses to 10^10.
   if ( fnFlavor != 0 && fnFlavor < 3 ) {
      logger.error["fastNLOQCDNUMAS::SetLHAPDFValues"] << "Less than 3 flavors is not supported! Aborted." << endl;
      exit(11);
   }
   if ( fnScheme == "variable" && fnFlavor < 6 ) {
      QMass[5] = 1.E10;
      if ( fnFlavor < 5 ) QMass[4] = 1.E10;
      if ( fnFlavor < 4 ) QMass[3] = 1.E10;
      fnFlavor = 0;
   }
   if ( PDFMemInfo.has_key("AlphaS_OrderQCD") ) {
      fnLoop = PDFMemInfo.get_entry_as<int>("AlphaS_OrderQCD") + 1;
   } else {
      fnLoop = PDFMemInfo.get_entry_as<int>("OrderQCD") + 1;
   }
   if ( fnLoop > 3 ) {
      logger.error["fastNLOQCDNUMAS::SetLHAPDFValues"] << "More than 3 loops is not supported! Aborted." << endl;
      exit(11);
   }
   fAlphasMz = PDFMemInfo.get_entry_as<double>("AlphaS_MZ");
}



// Printers
void fastNLOQCDNUMAS::PrintParmValues() {
   for ( int i = 0; i<6; i++ ) {
      cout << "fQMass[" << i << "] = " << QMass[i] << endl;
   }
   cout << "fMz       = " << fMz << endl;
   cout << "fnScheme  = " << fnScheme << endl;
   cout << "fnFlavor  = " << fnFlavor << endl;
   cout << "fnLoop    = " << fnLoop << endl;
   cout << "fAlphasMz = " << fAlphasMz << endl;
}



// Initialisation
void fastNLOQCDNUMAS::InitEvolveAlphas() {
   // Ensure reasonable values are set
   // TODO Really neccessary?
   char filename[] = " ";
   int len_filename = strlen(filename);
   int lun = 6;
   qcinit_(&lun, filename, len_filename);

   // LHAPDF LO=0 while QCDNUM LO=1
   int iord = fnLoop;

   // TODO Set correct Array in q2. maybe fnloreader. getQScale...
   double qarr[2] = {1.0, 1000000};
   double wgt[2] =  {1.0, 1.0};
   // Length of array
   int n= 2;
   // Number of grid points
   int nqin = 140;
   // Real number generated grid points
   int nqout = 0;
   // Create Q2 Grid
   gqmake_(qarr, wgt, &n, &nqin, &nqout);
   setord_(&iord);
   double r2 = fMz * fMz;
   setalf_(&fAlphasMz, &r2);
   // Get Indices of Flavor Thresholds (currently just the Q mass)
   double Q2Mass[6];
   for (int i = 0; i < 6; i++)
      Q2Mass[i] = QMass[i]*QMass[i];

   int iqc = iqfrmq_(&Q2Mass[3]) ;
   int iqb = iqfrmq_(&Q2Mass[4]);
   int iqt = iqfrmq_(&Q2Mass[5]);

   //cout << iqc << " " << iqb << " " << iqt << endl;
   //When fNFlavor = 0 VFNS if >0 then FFNS
   //iqc,b,t are neglected if fnflavor =0
   setcbt_(&fnFlavor, &iqc, &iqb, &iqt);
}

// Evolution
double fastNLOQCDNUMAS::EvolveAlphas(double Q) const {
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   double mu2 = Q*Q;
   int ierr = 9876;
   // Number of really used flavors
   int nf = 9; // KR: Why 9 ???!!!
   double as = asfunc_(&mu2, &nf , &ierr);
   if (ierr > 0)
      logger.error["EvolveAlphas"]<<"Alphas evolution failed. ierr = "<<ierr<<", Q = "<<Q<<endl;
   return as;
}

// Calculation
void fastNLOQCDNUMAS::CalcCrossSection() {
   InitEvolveAlphas();
   fastNLOLHAPDF::CalcCrossSection();
}
