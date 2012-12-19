// Author: Daniel Britzger
// DESY, 16/12/2012

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
//  FastNLOCJpdf
//
//////////////////////////////////////////////////////////////////////////

#ifndef FASTNLOCJPDF
#define FASTNLOCJPDF

#include "FastNLOReader.h"
#include <iostream>
#include "CRunDec.h"

using namespace std;

extern "C" {
   double cjpdf_(int *Iptn, double *x, double *Q);
   void setcj_(int *set);
}


class FastNLOCJpdf : public FastNLOReader {
private:
public:
   FastNLOCJpdf(string name) ;
   FastNLOCJpdf(string name, int iset) ;

   // ----- Printout ---- //
   void PrintRunDecValues();                    // Print values, which are passed to CRunDec for alpha_s evolution

   // ---- getters and setters CRunDec variables ---- //
   void   SetAlphasMz(double AlphasMz , bool ReCalcCrossSection = false);
   double GetAlphasMz() const {
      return fAlphasMz;
   };
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

   // ---- CJpdf ---- //
   void   SetISet(int iset);
   int    GetISet() {
      return fIset;
   }

protected:
   // inherited functions
   double EvolveAlphas(double Q) const ;

   // ---- CRunDec ---- //
   // declare variables as static, that all instances use same alpha_s evolution
   double fAlphasMz;
   static double fMz;
   static int fNf;
   static int fNloop;
   void InitReasonableRunDecValues();
   static CRunDec fcrundec;

   // ---- CJpdf ---- //
   bool InitPDF();
   vector<double> GetXFX(double xp, double muf) const ;
   int fIset;

};


//double FastNLOCJpdf::fAlphasMz=0.1180000654;
double FastNLOCJpdf::fMz=0.1180000654;
int FastNLOCJpdf::fNf = 5;
int FastNLOCJpdf::fNloop=2;
CRunDec FastNLOCJpdf::fcrundec=CRunDec(FastNLOCJpdf::fNf);


//______________________________________________________________________________


FastNLOCJpdf::FastNLOCJpdf(string name) : FastNLOReader(name) , fAlphasMz(0.1180000654) : fIset(100) {
   cout<<" WARNING! This code is not tested properly!"<<endl;
   InitReasonableRunDecValues();
}


//______________________________________________________________________________


FastNLOCJpdf::FastNLOCJpdf(string name, string LHAPDFFile, int PDFSet) : FastNLOLHAPDF(name,LHAPDFFile,PDFSet) , fAlphasMz(0.1180000654) :fIset(100)  {
   // do cross sections calculation, since everything is ready
   cout<<" WARNING! This code is not tested properly!"<<endl;
   InitReasonableRunDecValues();
   CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCJpdf::SetAlphasMz(double AlphasMz , bool ReCalcCrossSection) {
   debug["SetAlphasMz"]<<"Setting alpha_s(Mz)="<<AlphasMz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the alpha_s value at M_Z
   //
   fAlphasMz    = AlphasMz;             // new alpha_s value
   if (ReCalcCrossSection) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCJpdf::SetMz(double Mz , bool ReCalcCrossSection) {
   debug["SetMz"]<<"Setting MZ-"<<Mz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the Z-Boson mass
   //
   fMz    = Mz;             // new alpha_s value
   if (ReCalcCrossSection) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCJpdf::SetNloop(int nloop, bool ReCalcCrossSection) {
   debug["SetNloop"]<<"Setting n-loop="<<nloop<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set n loop calculation
   //
   fNloop    = nloop;             // new alpha_s value
   if (ReCalcCrossSection) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCJpdf::SetNf(int Nf , bool ReCalcCrossSection) {
   debug["SetNf"]<<"Setting number of flavors to "<<Nf<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the number of flavors
   //
   fNf    = Nf;             // new alpha_s value
   //fcrundec.SetConstants(fNf);
   if (ReCalcCrossSection) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCJpdf::InitReasonableRunDecValues() {
   fAlphasMz = 0.11840000000042; // PDG 2012 + epsilon(THE ANSWER ...) to avoid uninitialized a_s cache when explicitly setting the PDG2012 value
   fMz = 91.1876;
   SetNf(5);
   fNloop=2;
   if (info.GetSpeak()) {
      info["InitReasonableRunDecValues"]<<"Printing initialized CRunDecValues."<<endl;
      PrintRunDecValues();
   }
}


//______________________________________________________________________________

void FastNLOCJpdf::PrintRunDecValues() {
   static const string csep41("#########################################");
   cout<<csep41<<csep41<<endl;
   cout<<"CRunDec Values: Alphas(Mz)="<<fAlphasMz
       <<"\tMZ="<<fMz
       <<"\tn-flavors="<<fNf
       <<"\tn-loop="<<fNloop<<endl;
   cout<<csep41<<csep41<<endl;
}


//______________________________________________________________________________

double FastNLOCJpdf::EvolveAlphas(double Q) const {
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   //
   // Implementation of Alpha_s evolution as function of Mu_r.
   //
   // alpha_s evolution as given by the CRunDec program by
   // Barbara Schmidt, Matthias Steinhauser;
   // see also description in
   //   K.~G.~Chetyrkin, J.~H.~Kuhn and M.~Steinhauser,
   //   ``RunDec: A Mathematica package for running and decoupling of the strong
   //   coupling and quark masses,'' Comput.\ Phys.\ Commun.\  {\bf 133} (2000) 43
   //   [arXiv:hep-ph/0004189].

   double as_crundec = fcrundec.AlphasExact(fAlphasMz, fMz,Q, fNf, fNloop);
   return as_crundec;
}


//______________________________________________________________________________
//__________________________  CJpdf part starts here ___________________________
//______________________________________________________________________________

void FastNLOCJpdf::SetISet(int iset) {
   int fIset = iset;
   InitPDF();
}

bool FastNLOCJpdf::InitPDF() {
   // --- fastNLO user:
   //  Initalize PDF parameters if necessary
   //
   // It might be necessary that the PDF grid is recalculated/generated.
   setcj_(&fIset);
   return true;
}


//______________________________________________________________________________



vector<double> FastNLOCJpdf::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  x-proton and factorisation scale.
   //
   vector<double> xfx(13);
   for (int k = -5 ; k<=5 ; k++) {
      xfx[k+6] = xp * cjpdf_(&k,&xp,&muf);
   }
   return xfx;
}


//______________________________________________________________________________







#endif
