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

#ifndef FASTNLOUSER
#define FASTNLOUSER

#include "FastNLOReader.h"

using namespace std;


class FastNLOUser : public FastNLOReader {

public:
   FastNLOUser(string name);

protected:
   // inherited functions
   double EvolveAlphas(double Q) const ;
   bool InitPDF();
   vector<double> GetXFX(double xp, double muf) const ;

};



//______________________________________________________________________________


FastNLOUser::FastNLOUser(string name) : FastNLOReader(name) {
   // --- fastNLO user: if you have your own alpha_s routing in FastNLOUser::EvolveAlphas(double,double)
   //     it is convenient to automatically interface it here.
   //     Otherwise the FastNLO alpha_s evolution code Alphas.cc is used, or
   //     you have to call SetAlphasEvolution.
   //     It might be also convenient to make your scale choices here!
}



//______________________________________________________________________________


double FastNLOUser::EvolveAlphas(double Q) const {
   // --- fastNLO user:
   // Implementation of Alpha_s evolution as function of the
   // factorization scale [and alphas(Mz)].
   //
   return 0;
}


//______________________________________________________________________________


bool FastNLOUser::InitPDF() {
   // --- fastNLO user:
   //  Initalize PDF parameters if necessary
   //  return true, if successful initialization
   //  return false, if PDF initialization failed
   //
   return true;
}


//______________________________________________________________________________



vector<double> FastNLOUser::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  x-proton and factorisation scale.

   vector <double > xfx(13);
   // fastNLO user:
   //   include some function here to fill the parton density array xfx
   //   xfx[0]=tbar, ... ,  xfx[6]=gluon, ... , xfx[12]=t
   //
   return xfx;
}


//______________________________________________________________________________


#endif
