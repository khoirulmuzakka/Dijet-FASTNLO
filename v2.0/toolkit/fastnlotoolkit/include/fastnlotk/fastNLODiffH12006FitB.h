// Author: Daniel Britzger
// DESY, 02/04/2012

#ifndef fASTNLODIFFH12006FITB
#define fASTNLODIFFH12006FITB


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  FastNLODiffUSER                                                     //
//                                                                      //
//  FastNLODiffReader is a standalone code for reading                  //
//  diffractive FastNLO tables of version 2.0 for DIS processes         //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <string>
#include <cstdio>
#include <vector>
#include "fastnlotk/fastNLODiffReader.h"

using namespace std;

class fastNLODiffH12006FitB : public fastNLODiffReader {

public:

   fastNLODiffH12006FitB(string filename);
   ~fastNLODiffH12006FitB(void) {
      ;
   };

protected:

   // inherited functions
   double EvolveAlphas(double Q) const ;
   bool InitPDF();
   vector<double> GetDiffXFX(double xpom, double zpom, double muf) const ;

};



//______________________________________________________________________________


extern "C" {
   void diffpdf_(double* xpom, double*  zpom, double*  Q2, double *pdfs);
}



//______________________________________________________________________________




fastNLODiffH12006FitB::fastNLODiffH12006FitB(string filename) : fastNLODiffReader(filename) {
}


//______________________________________________________________________________


double fastNLODiffH12006FitB::EvolveAlphas(double Q) const {
   // --- fastNLO user:
   // Implementation of Alpha_s evolution as function of the
   // factorization scale [and alphas(Mz)].
   //

   static const int NF=4; // from h12006B_wrapper.h
   static const double b0 = (11. - 2./3.*NF);  // The beta coefficients of the QCD beta function
   static const  double b1 = (51. - 19./3.*NF);

   //      double lmd = 0.399; // according to matthias
   static const double lmd = 0.3395; // according to matthias
   double t = log(Q/lmd);
   double asMz = 1.0/(b0*t);

   return asMz*(1.0-b1/b0*asMz*log(2.0*t)) *TWOPI ;
}


//______________________________________________________________________________


bool fastNLODiffH12006FitB::InitPDF() {
   // --- fastNLO user:
   //  Initalize PDF parameters if necessary
   //
   // nothing todo!
   return true;
}


//______________________________________________________________________________



vector<double> fastNLODiffH12006FitB::GetDiffXFX(double xpom, double zpom, double muf) const {
   //
   //  GetDiffXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  xpom, zpom and factorisation scale.
   //
   vector < double > xfx(13);
   diffpdf_(&xpom,&zpom,&muf,&xfx[0]);
   //debug<<"xpom="<<xpom<<"\tzpom="<<zpom<<"\tmuf="<<muf<<"\tgluon = "<<xfx[6]<<endl;
   return xfx;
}


//______________________________________________________________________________


#endif
