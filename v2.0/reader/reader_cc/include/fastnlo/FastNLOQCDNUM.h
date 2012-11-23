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

#ifndef FASTNLOQCDNUM
#define FASTNLOQCDNUM

#include "FastNLOReader.h"

using namespace std;


extern "C"{
   void fpdfxq_(int *iset, const double *x, const double *q2, double *pdfs,int *ichk);
   void evolution_();
   double asfunc_( double* r2, int* nf  , int* ierr);
}


class FastNLOQCDNUM : public FastNLOReader {

public:
   FastNLOQCDNUM(string name);

protected:
   // inherited functions
   double EvolveAlphas(double Q ) const ;
   void InitPDF();
   vector<double> GetXFX(double xp, double muf) const ;
  
};



//______________________________________________________________________________


FastNLOQCDNUM::FastNLOQCDNUM(string name) : FastNLOReader(name) {
   // --- fastNLO user: if you have your own alpha_s routing in FastNLOQCDNUM::EvolveAlphas(double,double)
   //     it is convenient to automatically interface it here.
   //     Otherwise the FastNLO alpha_s evolution code Alphas.cc is used, or
   //     you have to call SetAlphasEvolution.
   //     It might be also convenient to make your scale choices here!
   FillAlphasCache();
}



//______________________________________________________________________________


double FastNLOQCDNUM::EvolveAlphas(double Q) const {
   // --- fastNLO user: 
   // Implementation of Alpha_s evolution as function of the
   // factorization scale [and alphas(Mz)].
   //
   // This function does not make use alphasMz but is using
   // the nominal values as defined in QCDNUM
   //
   double mu2 = Q*Q;
   int ierr = 9876;
   int nf = 9;
   double as = asfunc_( &mu2, &nf  , &ierr);
   if ( ierr > 0 )
      error["EvolveAlphas"]<<"Alphas evolution failed. ierr = "<<ierr<<", Q = "<<Q<<endl;
   return as;
}


//______________________________________________________________________________


void FastNLOQCDNUM::InitPDF(){
   // --- fastNLO user: 
   //  Initalize PDF parameters if necessary
   //
   // It might be necessary that the PDF grid is recalculated/generated.
   evolution_();
}


//______________________________________________________________________________



vector<double> FastNLOQCDNUM::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  x-proton and factorisation scale.
   //
   int iqnset = 1;
   int iqnchk = 0;
   double muf2 = muf*muf;
   vector < double > xfx(13);
   fpdfxq_(&iqnset, &xp, &muf2, &xfx[0], &iqnchk);
   return xfx;
}


//______________________________________________________________________________


#endif
