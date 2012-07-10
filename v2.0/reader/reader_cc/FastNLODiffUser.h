// Author: Daniel Britzger
// DESY, 02/04/2012

#ifndef FASTNLODIFFUSER
#define FASTNLODIFFUSER


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

class FastNLODiffUser : public FastNLODiffReader {

public:
   
   FastNLODiffUser(string filename);
   ~FastNLODiffUser(void){;};
  
protected:
   
   // inherited functions
   double EvolveAlphas(double Q, double alphasMz ) const ;
   void InitPDF();
   vector<double> GetDiffXFX(double xpom, double zpom, double muf) const ;

};



//______________________________________________________________________________




FastNLODiffUser::FastNLODiffUser(string filename) : FastNLODiffReader(filename)
{

}


//______________________________________________________________________________


double FastNLODiffUser::EvolveAlphas(double Q, double alphasMz ) const {
   // --- fastNLO user:
   // Implementation of Alpha_s evolution as function of the
   // factorization scale [and alphas(Mz)].
   //

   double as = 0;
   return as;
}


//______________________________________________________________________________


void FastNLODiffUser::InitPDF(){
   // --- fastNLO user:
   //  Initalize PDF parameters if necessary
   //
}


//______________________________________________________________________________



vector<double> FastNLODiffUser::GetDiffXFX(double xpom, double zpom, double muf) const {
   //
   //  GetDiffXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  xpom, zpom and factorisation scale.
   //
   vector < double > a(13);
   // ---- fastNLO user: Interface your pdf here!
   return a;
}


//______________________________________________________________________________


#endif