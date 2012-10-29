// Author: Daniel Britzger
// DESY, 11/08/2012

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
//  FastNLOCRunDec
//  This class inherits the PDF interface from
//  FastNLOLHAPDF, while the alpha_s evolution
//  is superseeded by the CRunDec evolution code.
//  All paramters of the evolution function are 
//  handled as static memeber, in order to have a
//  consistent evolution for all instances.
//  
//  The value of alpha_s at Mz however has to be
//  individually specified for each FastNLOCRunDec
//   class.
//
//////////////////////////////////////////////////////////////////////////

#ifndef FASTNLOCRUNDEC
#define FASTNLOCRUNDEC

#include "FastNLOReader.h"
#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include "FastNLOLHAPDF.h"
#include "CRunDec.h"

using namespace std;


class FastNLOCRunDec : public FastNLOLHAPDF {
private:
public:
   FastNLOCRunDec(string name) ;
   FastNLOCRunDec(string name, string LHAPDFFile, int PDFSet = 0) ;

   // ----- Printout ---- //
   void PrintRunDecValues();			// Print values, which are passed to CRunDec for alpha_s evolution
   
   // ---- getters and setters CRunDec variables ---- //
   void   SetAlphasMz( double AlphasMz , bool ReCalcCrossSection = false);
   double GetAlphasMz() const { return fAlphasMz; };
   void   SetMz( double Mz , bool ReCalcCrossSection = false );
   double GetMz() const { return fMz; };
   void   SetNf( int nf , bool ReCalcCrossSection = false );
   int    GetNf() const { return fNf; };
   void   SetNloop( int nloop , bool ReCalcCrossSection = false );
   int    GetNloop() const {return fNloop;};

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

};


//double FastNLOCRunDec::fAlphasMz=0.1180000654;
double FastNLOCRunDec::fMz=0.1180000654;
int FastNLOCRunDec::fNf = 5;
int FastNLOCRunDec::fNloop=2;
CRunDec FastNLOCRunDec::fcrundec=CRunDec(FastNLOCRunDec::fNf);


//______________________________________________________________________________


FastNLOCRunDec::FastNLOCRunDec(string name) : FastNLOLHAPDF(name) , fAlphasMz(0.1180000654) {
   info["FastNLOLHAPDF"]<<"Please initialize a PDF file using SetLHAPDFFilename( PDFFile ) and a PDF set using SetLHAPDFSet(int PDFSet)"<<std::endl;
   InitReasonableRunDecValues();
}


//______________________________________________________________________________


FastNLOCRunDec::FastNLOCRunDec(string name, string LHAPDFFile, int PDFSet) : FastNLOLHAPDF(name,LHAPDFFile,PDFSet) , fAlphasMz(0.1180000654)  {
   // do cross sections calculation, since everything is ready
   InitReasonableRunDecValues();
   CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCRunDec::SetAlphasMz( double AlphasMz , bool ReCalcCrossSection ){
   debug["SetAlphasMz"]<<"Setting alpha_s(Mz)="<<AlphasMz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the alpha_s value at M_Z
   //
   fAlphasMz    = AlphasMz;             // new alpha_s value
   if ( ReCalcCrossSection ) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCRunDec::SetMz( double Mz , bool ReCalcCrossSection ){
   debug["SetMz"]<<"Setting MZ-"<<Mz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the Z-Boson mass
   //
   fMz    = Mz;             // new alpha_s value
   if ( ReCalcCrossSection ) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCRunDec::SetNloop( int nloop, bool ReCalcCrossSection ){
   debug["SetNloop"]<<"Setting n-loop="<<nloop<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set n loop calculation
   //
   fNloop    = nloop;             // new alpha_s value
   if ( ReCalcCrossSection ) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCRunDec::SetNf( int Nf , bool ReCalcCrossSection ){
   debug["SetNf"]<<"Setting number of flavors to "<<Nf<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the number of flavors
   //
   fNf    = Nf;             // new alpha_s value
   //fcrundec.SetConstants(fNf);
   if ( ReCalcCrossSection ) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCRunDec::InitReasonableRunDecValues(){
   fAlphasMz = 0.11840000000042; // PDG 2012 + epsilon(THE ANSWER ...) to avoid uninitialized a_s cache when explicitly setting the PDG2012 value  
   fMz = 91.1876;
   SetNf(5);
   fNloop=2;
   if ( info.GetSpeak() ) {
      info["InitReasonableRunDecValues"]<<"Printing initialized CRunDecValues."<<endl;
      PrintRunDecValues();
   }
}

   
//______________________________________________________________________________

void FastNLOCRunDec::PrintRunDecValues(){
   static const string csep41("#########################################");
   cout<<csep41<<csep41<<endl;
   cout<<"CRunDec Values: Alphas(Mz)="<<fAlphasMz
	<<"\tMZ="<<fMz
	<<"\tn-flavors="<<fNf
	<<"\tn-loop="<<fNloop<<endl;
   cout<<csep41<<csep41<<endl;
}


//______________________________________________________________________________

double FastNLOCRunDec::EvolveAlphas(double Q) const {
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

//    // - do some initial print out
//    const string csep41("#########################################");
//    const string cseps = csep41 + csep41;
//    static bool first = true;
//    if ( first ) {
//       first = false;
//       cout << endl << " " << cseps << endl;
//       printf(" # ALPHAS-CRUNDEC: First call:\n");
//       cout << " " << cseps << endl;
//       //    printf(" # ALPHAS-CRUNDEC: PI              = %#18.15f\n",twopi/2.);
//       printf(" # ALPHAS-CRUNDEC: M_Z/GeV         = %#9.6f\n",Alphas::GetMz());
//       printf(" # ALPHAS-CRUNDEC: a_s(M_Z)        = %#9.6f\n",Alphas::GetAlphasMz());
//       printf(" # APLHAS-CRUNDEC: a_s loop        = %2i\n",Alphas::GetNLoop());
//       //    printf(" # APLHAS-CRUNDEC: flavor-matching = %s\n",(bFlavorMatching?"true":"false"));
//       printf(" # APLHAS-CRUNDEC: nf (M_Z)        = %2d\n",Alphas::GetNf());
//       cout << " " << cseps << endl;
//    }

   double as_crundec = fcrundec.AlphasExact(fAlphasMz, fMz,Q, fNf, fNloop);
   return as_crundec; 
}


//______________________________________________________________________________


#endif
