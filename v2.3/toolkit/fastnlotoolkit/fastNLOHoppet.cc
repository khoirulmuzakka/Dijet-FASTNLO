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
#include "fastnlotk/fastNLOHoppet.h"
#include "fastnlotk/HoppetInterface.h"
//#include "fastnlotk/speaker.h"

using namespace std;



//______________________________________________________________________________
//
fastNLOHoppet::fastNLOHoppet(std::string name, std::string LHAPDFFile, int PDFMem) : fastNLOLHAPDF(name, LHAPDFFile, PDFMem) {
   // Without PDF info use PDG values as default
   //   SetPDGValues();
   //   PrintParmValues();
   // Set initial values via LHAPDF6 info system
   SetLHAPDFValues(LHAPDFFile, PDFMem);
   // Print out values for checking
   //   PrintParmValues();
   fastNLOHoppet::InitPDF();
};



// Getters
double fastNLOHoppet::GetQMass(int pdgid) const {
   if (pdgid < 1 || pdgid > 6 ) {
      logger.error["fastNLOHoppet::GetQMass"]<<"PDG code out of quark index range 1-6! Aborted.\n";
      exit(1);
   }
   return HoppetInterface::QMass[pdgid];
}
double fastNLOHoppet::GetMz() const {
   return HoppetInterface::fMz;
}
std::string fastNLOHoppet::GetNScheme() const {
   return HoppetInterface::fnScheme;
}
int fastNLOHoppet::GetNFlavor() const {
   return HoppetInterface::fnFlavor;
}
int fastNLOHoppet::GetNLoop() const {
   return HoppetInterface::fnLoop;
}
double fastNLOHoppet::GetAlphasMz() const {
   return HoppetInterface::fAlphasMz;
};



// Setters
void fastNLOHoppet::SetQMass(int pdgid, double qmass) {
   HoppetInterface::QMass[pdgid] = qmass;
   HoppetInterface::InitHoppet(*this);
}
void fastNLOHoppet::SetMz(double Mz) {
   HoppetInterface::fMz = Mz;
   HoppetInterface::InitHoppet(*this);
}
void fastNLOHoppet::SetNFlavor(int nflavor) {
   HoppetInterface::fnFlavor = nflavor;
   HoppetInterface::InitHoppet(*this);
}
void fastNLOHoppet::SetNLoop(int nloop) {
   if ( nloop < 1 || nloop > 3 ) {
      logger.error["fastNLOHoppet::SetNLoop"] << "Illegal no. of loops nloop = " << nloop <<
         ", aborted! Only 1, 2, or 3 are allowed with HOPPET." << endl;
      exit(11);
   }
   HoppetInterface::fnLoop = nloop;
   HoppetInterface::InitHoppet(*this);
}
void fastNLOHoppet::SetAlphasMz(double AlphasMz) {
   HoppetInterface::fAlphasMz = AlphasMz;
   HoppetInterface::InitHoppet(*this);
}



// Combined Setters
void fastNLOHoppet::SetPDGValues() {
   // Initialize with PDG values
   HoppetInterface::QMass[0]  = PDG_MD;
   HoppetInterface::QMass[1]  = PDG_MU;
   HoppetInterface::QMass[2]  = PDG_MS;
   HoppetInterface::QMass[3]  = PDG_MC;
   HoppetInterface::QMass[4]  = PDG_MB;
   HoppetInterface::QMass[5]  = PDG_MT;
   HoppetInterface::fMz       = PDG_MZ;
   // Variable flavor number scheme
   HoppetInterface::fnFlavor = 0;
   // 2-loop alpha_s evolution
   HoppetInterface::fnLoop = 2;
   HoppetInterface::fAlphasMz = PDG_ASMZ;
   HoppetInterface::InitHoppet(*this);
}

void fastNLOHoppet::SetLHAPDFValues(std::string LHAPDFFile, int PDFMem) {
   // AlphaS_MZ can vary among PDF members, so we really need the PDF member info from LHAPDF
   const LHAPDF::PDFInfo PDFMemInfo(LHAPDFFile, PDFMem);
   HoppetInterface::QMass[0] = PDFMemInfo.get_entry_as<double>("MDown");
   HoppetInterface::QMass[1] = PDFMemInfo.get_entry_as<double>("MUp");
   HoppetInterface::QMass[2] = PDFMemInfo.get_entry_as<double>("MStrange");
   HoppetInterface::QMass[3] = PDFMemInfo.get_entry_as<double>("MCharm");
   HoppetInterface::QMass[4] = PDFMemInfo.get_entry_as<double>("MBottom");
   HoppetInterface::QMass[5] = PDFMemInfo.get_entry_as<double>("MTop");
   HoppetInterface::fMz      = PDFMemInfo.get_entry_as<double>("MZ");
   HoppetInterface::fnScheme = PDFMemInfo.get_entry_as<std::string>("FlavorScheme");
   if ( PDFMemInfo.has_key("AlphaS_NumFlavors") ) {
      HoppetInterface::fnFlavor = PDFMemInfo.get_entry_as<int>("AlphaS_NumFlavors");
   } else {
      HoppetInterface::fnFlavor = PDFMemInfo.get_entry_as<int>("NumFlavors");
   }
   // Variable flavor numbers are usually set via Nf = 0 in evolution code.
   // Ensure that fnFlavor is maximum Nf for variable flavor number scheme by
   // setting quark masses to 10^10.
   if ( HoppetInterface::fnFlavor != 0 && HoppetInterface::fnFlavor < 3 ) {
      logger.error["fastNLOHoppet::SetLHAPDFValues"] << "Less than 3 flavors is not supported! Aborted." << endl;
      exit(11);
   }
   if ( HoppetInterface::fnScheme == "variable" && HoppetInterface::fnFlavor < 6 ) {
      HoppetInterface::QMass[5] = 1.E10;
      if ( HoppetInterface::fnFlavor < 5 ) HoppetInterface::QMass[4] = 1.E10;
      if ( HoppetInterface::fnFlavor < 4 ) HoppetInterface::QMass[3] = 1.E10;
      HoppetInterface::fnFlavor = 0;
   }
   if ( PDFMemInfo.has_key("AlphaS_OrderQCD") ) {
      HoppetInterface::fnLoop = PDFMemInfo.get_entry_as<int>("AlphaS_OrderQCD") + 1;
   } else {
      HoppetInterface::fnLoop = PDFMemInfo.get_entry_as<int>("OrderQCD") + 1;
   }
   if ( HoppetInterface::fnLoop > 3 ) {
      logger.error["fastNLOHoppet::SetLHAPDFValues"] << "More than 3 loops is not supported! Aborted." << endl;
      exit(11);
   }
   HoppetInterface::fAlphasMz = PDFMemInfo.get_entry_as<double>("AlphaS_MZ");
   HoppetInterface::InitHoppet(*this);
}



// Printers
void fastNLOHoppet::PrintParmValues() {
   for ( int i = 0; i<6; i++ ) {
      cout << "fQMass[" << i << "] = " << HoppetInterface::QMass[i] << endl;
   }
   cout << "fMz       = " << HoppetInterface::fMz << endl;
   cout << "fnScheme  = " << HoppetInterface::fnScheme << endl;
   cout << "fnFlavor  = " << HoppetInterface::fnFlavor << endl;
   cout << "fnLoop    = " << HoppetInterface::fnLoop << endl;
   cout << "fAlphasMz = " << HoppetInterface::fAlphasMz << endl;
}



// Initialisation
bool fastNLOHoppet::InitPDF() {
   bool init = fastNLOLHAPDF::InitPDF();
   HoppetInterface::InitHoppet(*this);
   return init;
}

// Evolution
double fastNLOHoppet::EvolveAlphas(double Q ) const {
   return HoppetInterface::EvolveAlphas(Q);
}

std::vector<double> fastNLOHoppet::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   return HoppetInterface::GetXFX(xp, muf);
}
