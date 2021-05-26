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
//  FastNLOAlhpas
//  This class inherits the PDF interface from
//  FastNLOLHAPDF, while the alpha_s evolution
//  is superseeded by the Alphas.h class.
//lhasub
//////////////////////////////////////////////////////////////////////////


#ifndef FASTNLOHOPPET
#define FASTNLOHOPPET

#include "FastNLOReader.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
#include "speaker.h"
#include "FastNLOLHAPDF.h"
#include "hoppet_v1.h"

using namespace std;


class FastNLOHoppet : public FastNLOLHAPDF {

   public:
      FastNLOHoppet(string name) : FastNLOLHAPDF(name) {
         //Set some meaningful values
         SetPDGValues();
      };
      FastNLOHoppet(string name, string LHAPDFFile, int PDFSet = 0) : 
         FastNLOLHAPDF(name,LHAPDFFile,PDFSet), 
         fAlphasMz(0.1184) {
            //Set some meaningful values
            SetPDGValues();
         };

      // ---- Alphas vars ---- //
      void InitHoppet();
      // Setters
      void SetMz(double Mz);
      void SetNFlavor(int nflavor);
      void SetNLoop(int nloop);
      void SetQMass(int pdgid, double qmass);
      void SetAlphasMz(double AlphasMz , bool ReCalcCrossSection = false);
      void SetLHAPDFValues();
      void SetPDGValues();
      // Getters
      double GetMz() const {
         return fMz;
      }
      double GetQMass(int pdgid) const {
         return QMass[pdgid];
      }
      int GetNFlavor() const {
         return fnFlavor;
      }
      int GetNLoop() const {
         return fnLoop;
      }
      double GetAlphasMz() const {
         return fAlphasMz;
      };



   protected:

      // inherited functions
      double EvolveAlphas(double Q) const ;
      //bool InitPDF();
      vector<double> GetXFX(double xp, double muf) const ;
      static void LHAsub(const double&, const double&, double*);
      // ---- Alphas vars ---- //
      double fAlphasMz;
      double fMz;
      int fnFlavor;
      int fnLoop;
      double QMass[6];
      // ___ //

};


void FastNLOHoppet::SetPDGValues() {
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

void FastNLOHoppet::SetLHAPDFValues() {
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

void FastNLOHoppet::InitHoppet() {

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

void FastNLOHoppet::SetMz(double Mz) {
   fMz = Mz;
   InitHoppet();
}

void FastNLOHoppet::SetNFlavor(int nflavor) {
   fnFlavor = nflavor;
   InitHoppet();
}

void FastNLOHoppet::SetNLoop(int  nloop) {
   fnLoop = nloop;
   InitHoppet();
}

void FastNLOHoppet::SetQMass(int pdgid, double qmass) {
   QMass[pdgid] = qmass;
   InitHoppet();
}

void FastNLOHoppet::SetAlphasMz(double AlphasMz , bool ReCalcCrossSection) {

   fAlphasMz    = AlphasMz;             // new alpha_s value
   InitHoppet();
}

void FastNLOHoppet::LHAsub(const double & x, const double & Q, double * pdf) {
   //
   //Provides PDF for Hoppet
   //for (int i=0; i<13; i++)
   //{
   //pdf[i] = LHAPDF::xfx(x, Q, i-6);
   //}
}

double FastNLOHoppet::EvolveAlphas(double Q ) const {
   return hoppetalphas_(Q);
}

vector<double> FastNLOHoppet::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   vector<double> xfx(13);
   //Hoppet PDF Evolution
   //hoppeteval_(xp, muf, &xfx[0]);
   xfx = LHAPDF::xfx(xp, muf);
   return xfx;
}

//______________________________________________________________________________

#endif
