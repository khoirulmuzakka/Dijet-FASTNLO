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


#ifndef FASTNLOALPHAS
#define FASTNLOALPHAS

#include "FastNLOReader.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
#include "speaker.h"
#include "FastNLOLHAPDF.h"
#include "Alphas.h"

using namespace std;


class FastNLOAlphas : public FastNLOLHAPDF {

public:
   FastNLOAlphas(string name) : FastNLOLHAPDF(name) {
      ;
   };
   FastNLOAlphas(string name, string LHAPDFFile, int PDFSet = 0) : FastNLOLHAPDF(name,LHAPDFFile,PDFSet), fAlphasMz(0.1184) {
      ;
   };

   // ---- Alphas vars ---- //
   void SetAlphasMz(double AlphasMz , bool ReCalcCrossSection = false);
   double GetAlphasMz() const {
      return fAlphasMz;
   };
   void SetGRVtoPDG2012_2loop();

   // ---- getters and setters alphas variables ---- //
   void   SetMz(double Mz , bool ReCalcCrossSection = false);
   double GetMz() const {
      return fMz;
   };
   void   SetNf(int nf , bool ReCalcCrossSection = false);
   int    GetNf() const {
      return fNf;
   };
   void   SetNloop(int nloop , bool ReCalcCrossSection = false);
   int    GetNloop() const {
      return fNloop;
   };


protected:
   // inherited functions
   double EvolveAlphas(double Q) const ;

   // ---- Alphas vars ---- //
   double fAlphasMz;
   double fMz;
   double fNf;
   double fNloop;
};


//______________________________________________________________________________


void FastNLOAlphas::SetAlphasMz(double AlphasMz , bool ReCalcCrossSection) {
   debug["SetAlphasMz"]<<"Setting alpha_s(Mz)="<<AlphasMz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the alpha_s value at M_Z
   //
   fAlphasMz    = AlphasMz;             // new alpha_s value
   if (ReCalcCrossSection) CalcCrossSection();
}

void FastNLOAlphas::SetMz(double Mz , bool ReCalcCrossSection) {
   debug["SetMz"]<<"Setting MZ-"<<Mz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the Z-Boson mass
   //
   fMz    = Mz;             // new alpha_s value
   if (ReCalcCrossSection) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOAlphas::SetNloop(int nloop, bool ReCalcCrossSection) {
   debug["SetNloop"]<<"Setting n-loop="<<nloop<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set n loop calculation
   //
   fNloop    = nloop;             // new alpha_s value
   if (ReCalcCrossSection) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOAlphas::SetNf(int Nf , bool ReCalcCrossSection) {
   debug["SetNf"]<<"Setting number of flavors to "<<Nf<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the number of flavors
   //
   fNf    = Nf;             // new alpha_s value
   //fcrundec.SetConstants(fNf);
   if (ReCalcCrossSection) CalcCrossSection();
}

//______________________________________________________________________________


double FastNLOAlphas::EvolveAlphas(double Q) const {
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'alphasMz' is not used here!
   //
   return Alphas::CalcAlphasMu(Q , fAlphasMz);
}


//______________________________________________________________________________


void FastNLOAlphas::SetGRVtoPDG2012_2loop() {
   info["SetGrVtoPDF2012"]<<"Resetting to GRV Alphas::Alphas evolution."<<endl;
   Alphas::SetMz(91.1876); // PDG 2012
   Alphas::SetNf(5);
   Alphas::SetNLoop(2);
   Alphas::SetFlavorMatchingOn(false);
   if (info.GetSpeak()) {
      info<<"Calling Alphas::PrintInfo()."<<endl;
      info<<"Alpha_s(Mz) value is taken from FastNLOAlphas, instead of Alphas::Alphas."<<endl;
      Alphas::PrintInfo();
   }
}


//______________________________________________________________________________


#endif
