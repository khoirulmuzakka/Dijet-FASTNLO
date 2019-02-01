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

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
#include "fastnlotk/fastNLOLHAPDF.h"
#include "fastnlotk/fastNLOCRunDec.h"

using namespace std;



//______________________________________________________________________________
//
fastNLOCRunDec::fastNLOCRunDec(string name) : fastNLOLHAPDF(name) {
   // Without PDF info use PDG values as default
   InitCRunDec();
   SetPDGValues();
   // Print out values for checking
   //   PrintParmValues();
};
fastNLOCRunDec::fastNLOCRunDec(string name, string LHAPDFFile, int PDFMem) : fastNLOLHAPDF(name, LHAPDFFile, PDFMem) {
   // Set initial values via LHAPDF6 info system
   InitCRunDec();
   SetLHAPDFValues(LHAPDFFile, PDFMem);
   // Print out values for checking
   //   PrintParmValues();
};



// Getters
double fastNLOCRunDec::GetQMass(int pdgid) const {
   if (pdgid < 1 || pdgid > 6 ) {
      logger.error["fastNLOCRunDec::GetQMass"]<<"PDG code out of quark index range 1-6! Aborted.\n";
      exit(1);
   }
   return QMass[pdgid];
}
double fastNLOCRunDec::GetMz() const {
   return fMz;
}
std::string fastNLOCRunDec::GetNScheme() const {
   return fnScheme;
}
int fastNLOCRunDec::GetNFlavor() const {
   return fnFlavor;
}
int fastNLOCRunDec::GetNLoop() const {
   return fnLoop;
}
double fastNLOCRunDec::GetAlphasMz() const {
   return fAlphasMz;
};



// Setters
void fastNLOCRunDec::SetQMass(int pdgid, double qmass) {
   QMass[pdgid] = qmass;
}
void fastNLOCRunDec::SetMz(double Mz) {
   fMz = Mz;
}
void fastNLOCRunDec::SetNFlavor(int nflavor) {
   fnFlavor = nflavor;
}
void fastNLOCRunDec::SetNLoop(int nloop) {
   if ( nloop < 1 || nloop > 4 ) {
      logger.error["fastNLOHoppet::SetNLoop"] << "Illegal no. of loops nloop = " << nloop <<
         ", aborted! Only 1, 2, 3, or 4 are allowed with RUNDEC." << endl;
      exit(11);
   }
   fnLoop = nloop;
}
void fastNLOCRunDec::SetAlphasMz(double AlphasMz) {
   fAlphasMz    = AlphasMz;
}



// Combined Setters
void fastNLOCRunDec::SetPDGValues() {
   // Initialize with PDG values
   QMass[0] = PDG_MD;
   QMass[1] = PDG_MU;
   QMass[2] = PDG_MS;
   QMass[3] = PDG_MC;
   QMass[4] = PDG_MB;
   QMass[5] = PDG_MT;
   fMz      = PDG_MZ;
   // Variable flavor number scheme
   fnFlavor = 0;
   // 2-loop alpha_s evolution
   fnLoop = 2;
   fAlphasMz = PDG_ASMZ;
}

void fastNLOCRunDec::SetLHAPDFValues(std::string LHAPDFFile, int PDFMem) {
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
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
      logger.error["fastNLOCRunDec::SetLHAPDFValues"] << "Less than 3 flavors is not supported! Aborted." << endl;
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
   fAlphasMz = PDFMemInfo.get_entry_as<double>("AlphaS_MZ");
#else
   // TODO: Remove old LHAPDF5 stuff
   cerr << "LHAPDF5 not supported anymore! Please update." << endl;
   exit(1);
   // fAlphasMz = LHAPDF::alphasPDF(fiPDFMember,fMz);
   // fnLoop = LHAPDF::getOrderAlphaS(fiPDFMember) + 1;
   // for (int i = 0; i < 6; i++)
   //    QMass[i] = LHAPDF::getQMass(fiPDFMember,i+1);
   // fAlphasMz = LHAPDF::alphasPDF(fMz);
   // fnLoop = LHAPDF::getOrderAlphaS() + 1;
   // for (int i = 0; i < 6; i++)
   //    QMass[i] = LHAPDF::getQMass(i+1);
#endif
}



// Printers
void fastNLOCRunDec::PrintParmValues() {
   for ( int i = 0; i<6; i++ ) {
      cout << "fQMass[" << i << "] = " << QMass[i] << endl;
   }
   cout << "fMz       = " << fMz << endl;
   cout << "fnScheme  = " << fnScheme << endl;
   cout << "fnFlavor  = " << fnFlavor << endl;
   cout << "fnLoop    = " << fnLoop << endl;
   cout << "fAlphasMz = " << fAlphasMz << endl;
}


//______________________________________________________________________________
void fastNLOCRunDec::InitCRunDec() {

   crundec = new CRunDec();

}



//______________________________________________________________________________
double fastNLOCRunDec::EvolveAlphas(double Q) const {
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   //FFNS
   if (fnFlavor != 0) {
      return crundec->AlphasExact(fAlphasMz, fMz, Q, fnFlavor, fnLoop);
   }

   //TODO: Replace with better code
   //VFNS
   //Always evolving from Mz as starting scale
   //Crossing the mt threshold
   if (Q > QMass[5]) {
      crundec->nfMmu[0].nf = 6;
      crundec->nfMmu[0].Mth = QMass[5];
      crundec->nfMmu[0].muth = QMass[5];
      return crundec->AlL2AlH(fAlphasMz, fMz, crundec->nfMmu, Q, fnLoop);
   }
   //Not Crossing any threshold
   else if ( Q > QMass[4]) {
      return crundec->AlphasExact(fAlphasMz, fMz, Q, 5, fnLoop);
   }
   //Crossing mb threshold
   else if (Q > QMass[3]) {
      crundec->nfMmu[0].nf = 5;
      crundec->nfMmu[0].Mth = QMass[4];
      crundec->nfMmu[0].muth = QMass[4];
      return crundec->AlH2AlL(fAlphasMz, fMz, crundec->nfMmu, Q, fnLoop);
   }
   //Crossing mc and mb threshold
   else {
      crundec->nfMmu[0].nf = 5;
      crundec->nfMmu[0].Mth = QMass[4];
      crundec->nfMmu[0].muth = QMass[4];
      crundec->nfMmu[1].nf = 4;
      crundec->nfMmu[1].Mth = QMass[3];
      crundec->nfMmu[1].muth = QMass[3];
      return crundec->AlH2AlL(fAlphasMz, fMz, crundec->nfMmu, Q, fnLoop);
   }
}
