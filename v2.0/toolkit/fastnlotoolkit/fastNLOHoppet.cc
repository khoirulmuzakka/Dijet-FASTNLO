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
//lhasub
//////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
//#include "fastnlotk/fastNLOReader.h"
//#include "fastnlotk/speaker.h"
#include "fastnlotk/fastNLOLHAPDF.h"
#include "fastnlotk/fastNLOHoppet.h"
#include "hoppet_v1.h"

using namespace std;



fastNLOHoppet::fastNLOHoppet(string name) : fastNLOLHAPDF(name) {
    //Set some meaningful values
    SetPDGValues();
};
fastNLOHoppet::fastNLOHoppet(string name, string LHAPDFFile, int PDFSet = 0) :
    fastNLOLHAPDF(name,LHAPDFFile,PDFSet),
    fAlphasMz(0.1184) {
        //Set some meaningful values
        SetPDGValues();
    };
// Getters
double fastNLOHoppet::GetMz() const {
    return fMz;
}
double fastNLOHoppet::GetQMass(int pdgid) const {
    return QMass[pdgid];
}
int fastNLOHoppet::GetNFlavor() const {
    return fnFlavor;
}
int fastNLOHoppet::GetNLoop() const {
    return fnLoop;
}
double fastNLOHoppet::GetAlphasMz() const {
    return fAlphasMz;
};


void fastNLOHoppet::SetPDGValues() {
   // Initialize with PDG values
   QMass[0]  = PDG_MD;
   QMass[1]  = PDG_MU;
   QMass[2]  = PDG_MS;
   QMass[3]  = PDG_MC;
   QMass[4]  = PDG_MB;
   QMass[5]  = PDG_MT;
   fMz       = PDG_MZ;
   fAlphasMz = PDG_ASMZ;
   //Variable flavor number scheme
   fnFlavor = -1;
   //2-loop alpha_s evolution
   fnLoop = 2;
   InitHoppet();
}

void fastNLOHoppet::SetLHAPDFValues() {
   //Be sure LHAPDF is initialized when reading the properties
   if (fchksum == 0 || fchksum != CalcChecksum(1.)) {
      InitPDF();
   }
   //How to read LHAPDF Mz???
   fMz = PDG_MZ;
   fAlphasMz = LHAPDF::alphasPDF(fMz);
   fnLoop = LHAPDF::getOrderAlphaS() + 1;
   fnFlavor = LHAPDF::getNf();
   for (int i = 0; i < 6; i++)
      QMass[i] = LHAPDF::getQMass(i+1);
   InitHoppet();
}

void fastNLOHoppet::InitHoppet() {

   //Define Grid for alphaS and PDF evolution
   double ymax = 12.0;
   double dy = 0.1;
   int order = -6;
   double dlnlnQ = dy/4.0;
   double Qmin = 1.0;
   double Qmax = 28000;
   hoppetStartExtended( ymax, dy, Qmin, Qmax, dlnlnQ, fnLoop, order, factscheme_MSbar);

   //If fnFlavor smaller than 1 use VFNS (NNPDF reports nf=-1)
   if (fnFlavor >= 1)
      hoppetsetffn_(fnFlavor);
   else
      hoppetsetpolemassvfn_(QMass[3], QMass[4], QMass[5]);
   //Init PDF/As evolution once with dummy LHAsub
   //Evolved PDFs with Hoppet will be nonsense
   hoppetevolve_(fAlphasMz, fMz, fnLoop, 1.0, &LHAsub, 2.00001);
}

void fastNLOHoppet::SetMz(double Mz) {
   fMz = Mz;
   InitHoppet();
}

void fastNLOHoppet::SetNFlavor(int nflavor) {
   fnFlavor = nflavor;
   InitHoppet();
}

void fastNLOHoppet::SetNLoop(int  nloop) {
   fnLoop = nloop;
   InitHoppet();
}

void fastNLOHoppet::SetQMass(int pdgid, double qmass) {
   QMass[pdgid] = qmass;
   InitHoppet();
}

void fastNLOHoppet::SetAlphasMz(double AlphasMz , bool ReCalcCrossSection) {

   fAlphasMz    = AlphasMz;             // new alpha_s value
   InitHoppet();
}

void fastNLOHoppet::LHAsub(const double & x, const double & Q, double * pdf) {
   //
   //Provides PDF for Hoppet
   for (int i=0; i<13; i++)
   {
   pdf[i] = LHAPDF::xfx(x, Q, i-6);
   }
}

double fastNLOHoppet::EvolveAlphas(double Q ) const {
   return hoppetalphas_(Q);
}

vector<double> fastNLOHoppet::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   vector<double> xfx(13);
   //Hoppet PDF Evolution
   hoppeteval_(xp, muf, &xfx[0]);
   return xfx;
}
