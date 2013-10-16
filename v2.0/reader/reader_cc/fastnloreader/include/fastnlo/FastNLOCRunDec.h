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

class FastNLOCRunDec : public FastNLOLHAPDF {

   public:
      FastNLOCRunDec(string name) : FastNLOLHAPDF(name) {
         InitCRunDecPDG();
      };
      FastNLOCRunDec(string name, string LHAPDFFile, int PDFSet = 0) : FastNLOLHAPDF(name,LHAPDFFile,PDFSet), fAlphasMz(0.1184) {
         InitCRunDecPDG();
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
   void InitCRunDecPDG();
   double fAlphasMz;
   double fMz;
   double fnFlavor;
   double fnLoop;
   double QMass[6];


};


//______________________________________________________________________________
void FastNLOCRunDec::InitCRunDecPDG() {
   crundec = new CRunDec();
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
   fnFlavor = 0;
   //2-loop alpha_s evolution
   fnLoop = 2;
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
