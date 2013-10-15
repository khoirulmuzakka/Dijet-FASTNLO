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
//
//////////////////////////////////////////////////////////////////////////


#ifndef FASTNLOCRUNDEC
#define FASTNLOCRUNDEC

#include "FastNLOReader.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
#include "speaker.h"
#include "FastNLOLHAPDF.h"
//#include "Alphas.h"
#include "CRunDec.h"

using namespace std;

//should not be defined here...
namespace PDG {
   static const double Mz = 91.1876;
   static const double AlphasMz = 0.1184;
   static const double QMass[6] = {0.0023, 0.0048, 0.0035, 1.275, 4.18, 173.03};
}

class FastNLOCRunDec : public FastNLOLHAPDF {

   public:
      FastNLOCRunDec(string name) : FastNLOLHAPDF(name) {
         InitCRunDec();
      };
      FastNLOCRunDec(string name, string LHAPDFFile, int PDFSet = 0) : FastNLOLHAPDF(name,LHAPDFFile,PDFSet), fAlphasMz(0.1184) {
         InitCRunDec();
      };

      // ---- Alphas vars ---- //
      // Setters
      void SetMz(double Mz);
      void SetNFlavor(double nflavor);
      void SetNLoop(double nloop);
      void SetQMass(int pdgid, double qmass);
      void SetAlphasMz(double AlphasMz , bool ReCalcCrossSection = false);
      void SetLHAPDFValues();
      // Getters
      double GetMz(double Mz) const {
         return fMz;
      }
      double GetQMass(int pdgid) const {
         return QMass[pdgid];
      }
      double GetNFlavor(double nflavor) const {
         return fnFlavor;
      }
      double GetNLoop(double nloop) const {
         return fnLoop;
      }
      double GetAlphasMz() const {
         return fAlphasMz;
      };



   protected:

   // inherited functions
   double EvolveAlphas(double Q) const ;

   // ---- Alphas vars ---- //
   CRunDec *crundec;
   void InitCRunDec();
   double fAlphasMz;
   double fMz;
   double fnFlavor;
   double fnLoop;
   double QMass[6];


};


//______________________________________________________________________________
void FastNLOCRunDec::InitCRunDec() {
   crundec = new CRunDec();
   fMz = PDG::Mz;
   fAlphasMz = PDG::AlphasMz;
   //Variable Flavors
   fnFlavor = 0;
   fnLoop = 2;
   for (int i = 0; i < 6; i++)
      QMass[i] = PDG::QMass[i];
}

void FastNLOCRunDec::SetLHAPDFValues() {
   FillPDFCache();
   fAlphasMz = LHAPDF::alphasPDF(fMz);
   fnLoop = LHAPDF::getOrderAlphaS() + 1;
   for (int i = 0; i < 6; i++)
      QMass[i] = LHAPDF::getQMass(i+1);
  
}

void FastNLOCRunDec::SetMz(double Mz) {
   fMz = Mz;
}

void FastNLOCRunDec::SetNFlavor(double nflavor) {
   fnFlavor = nflavor;
}

void FastNLOCRunDec::SetNLoop(double nloop) {
   fnLoop = nloop;
}

void FastNLOCRunDec::SetQMass(int pdgid, double qmass) {
   QMass[pdgid] = qmass;
}

void FastNLOCRunDec::SetAlphasMz(double AlphasMz , bool ReCalcCrossSection) {
   debug["SetAlphasMz"]<<"Setting alpha_s(Mz)="<<AlphasMz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the alpha_s value at M_Z
   //
   fAlphasMz    = AlphasMz;             // new alpha_s value
   //if (ReCalcCrossSection) CalcCrossSection();

}

//______________________________________________________________________________
double FastNLOCRunDec::EvolveAlphas(double Q) const {
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'lphasMz' is not used here!
   //
   //return crundec.AlphasExact(fAlphasMz, Mz, Q, nflavor, nloop);
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



//______________________________________________________________________________


#endif
