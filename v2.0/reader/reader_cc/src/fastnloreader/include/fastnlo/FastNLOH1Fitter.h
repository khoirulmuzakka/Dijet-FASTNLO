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

#ifndef FASTNLOH1FITTER
#define FASTNLOH1FITTER

#include "FastNLOReader.h"
#include "get_pdfs.h"

using namespace std;


class FastNLOH1Fitter : public FastNLOReader {

public:
   FastNLOH1Fitter(string name);

protected:
   // inherited functions
   double EvolveAlphas(double Q) const ;
   void InitPDF();
   vector<double> GetXFX(double xp, double muf) const ;
};



//______________________________________________________________________________


FastNLOH1Fitter::FastNLOH1Fitter(string name) : FastNLOReader(name) {
   // --- fastNLO user: if you have your own alpha_s routing in FastNLOH1Fitter::EvolveAlphas(double,double)
   //     it is convenient to automatically interface it here.
   //     Otherwise the FastNLO alpha_s evolution code Alphas.cc is used, or
   //     you have to call SetAlphasEvolution.
   //     It might be also convenient to make your scale choices here!

   SetAlphasEvolution(kExternAs);
}



//______________________________________________________________________________


double FastNLOH1Fitter::EvolveAlphas(double Q) const {
   // --- fastNLO user:
   // Implementation of Alpha_s evolution as function of the
   // factorization scale [and alphas(Mz)].
   //
   // This function does not make use alphasMz but is using
   // the nominal values as defined in QCDNUM
   //

   double mu2 = Q*Q;
   return HF_GET_ALPHAS_WRAP( &mu2 );
}


//______________________________________________________________________________


void FastNLOH1Fitter::InitPDF(){
   // --- fastNLO user:
   //  Initalize PDF parameters if necessary
   //

   // nothing todo
}


//______________________________________________________________________________



vector<double> FastNLOH1Fitter::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  x-proton and factorisation scale.
   //
   //! return  pdf grid 'xfx'
   vector < double > xfx(13);
   double muf2       = muf*muf;
   HF_GET_PDFS_WRAP(&xp, &muf2, &xfx[0]);
   return xfx;
}


//______________________________________________________________________________


#endif
