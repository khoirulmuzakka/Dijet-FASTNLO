// Author: Daniel Britzger
// DESY, 02/04/2012

#ifndef FASTNLODIFFTESTDUMMY
#define FASTNLODIFFTESTDUMMY


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  FastNLODiffUSER                                                     //
//                                                                      //
//  FastNLODiffReader is a standalone code for reading                  //
//  diffractive FastNLO tables of version 2.0 for DIS processes         //
//                                                                      //
//  FastNLO is developed by                                             //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch         //
//    (publication in preparation)                                      //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <string>
#include <cstdio>
#include <vector>
#include "FastNLODiffReader.h"

using namespace std;

class FastNLODiffTestDummy : public FastNLODiffReader {

public:
   
   FastNLODiffTestDummy(string filename);
   ~FastNLODiffTestDummy(void){;};
  
protected:
   
   // inherited functions
   double EvolveAlphas(double Q) const ;
   bool InitPDF();
   vector<double> GetDiffXFX(double xpom, double zpom, double muf) const ;

};



//______________________________________________________________________________




FastNLODiffTestDummy::FastNLODiffTestDummy(string filename) : FastNLODiffReader(filename)
{
}


//______________________________________________________________________________


double FastNLODiffTestDummy::EvolveAlphas(double Q) const {
   // --- fastNLO user:
   // Implementation of Alpha_s evolution as function of the
   // factorization scale [and alphas(Mz)].
   //

   static const int NF=4; // from h12006B_wrapper.h
   static const double b0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
   static const  double b1 =  (51. - 19./3.*NF);

   //      double lmd = 0.399; // according to matthias
   static const double lmd = 0.3395; // according to matthias
   double t = log(Q/lmd);
   double asMz = 1.0/(b0*t);

   return asMz*(1.0-b1/b0*asMz*log(2.0*t)) *TWOPI ;
}


//______________________________________________________________________________


bool FastNLODiffTestDummy::InitPDF(){
   // --- fastNLO user:
   //  Initalize PDF parameters if necessary
   //
   // nothing todo!
   return true;
}


//______________________________________________________________________________



vector<double> FastNLODiffTestDummy::GetDiffXFX(double xpom, double zpom, double muf) const {
   //
   //  GetDiffXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  xpom, zpom and factorisation scale.
   //
   vector < double > xfx(13);

   double guess = 1/xpom/(-2.*log10(zpom));
   for(int i = 2; i< 11; i++ )
      xfx[i]=guess;
   return xfx;
}


//______________________________________________________________________________


#endif
