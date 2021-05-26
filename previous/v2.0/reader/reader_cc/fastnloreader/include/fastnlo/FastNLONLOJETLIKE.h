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

#ifndef FASTNLONLOJETLIKE
#define FASTNLONLOJETLIKE

#include "FastNLOReader.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
#include "speaker.h"

using namespace std;


class FastNLONLOJETLIKE : public FastNLOReader {
public:
   FastNLONLOJETLIKE(string name);
protected:
   // inherited functions
   double EvolveAlphas(double Q) const ;
   void InitPDF();
   vector<double> GetXFX(double xp, double muf) const ;
};



//______________________________________________________________________________


FastNLONLOJETLIKE::FastNLONLOJETLIKE(string name) : FastNLOReader(name) {
   LHAPDF::setVerbosity(LHAPDF::SILENT);
   LHAPDF::initPDFSet("cteq6m.LHpdf");
   LHAPDF::initPDF(0);
   // do cross sections calculation, since everything is yet ready
   FillAlphasCache();
   FillPDFCache();
   CalcCrossSection();
}



//______________________________________________________________________________


double FastNLONLOJETLIKE::EvolveAlphas(double Q) const {
   //
   // Implementation of Alpha_s evolution as function of Mu_r.
   // this is the evolution, which is used by nlojet++ and cteq6m.
   // Be aware of the Mz-value from 2001 of 91.70.
   //
   // Values used in nlojet++ and cteq6m:
   //   alphasmz   = 0.1179
   //   double b0  = 1.2202
   //   double b1  = 0.4897
   //   double Mz  = 91.70
   //
   // as evolution by lhpdf.c
   // please cite the nlojet++ references.
   //
   // the original parameters from the cteq6 pdf
   //   // Alpha QCD //
   //   1, 1, 0, -1, 0.1179, 91.70, 1.3, 4.5, 180.0,
   //   double BETA0 =  (11. - 2./3.*NF); // The beta coefficients of the QCD beta function
   //   double BETA1 =  (51. - 19./3.*NF);

   static const double alphasMZ = 0.1179;
   static const double b0  = 1.2202;
   static const double b1  = 0.4897;
   static const double Mz     = 91.70;
   // double b2 = 0.1913;
   //double Mz   = 91.187;

   double L = log(Q/Mz);
   L = (b0 + alphasMZ*b1)*L;

   return alphasMZ/(1.0 + alphasMZ*L);
}


//______________________________________________________________________________


void FastNLONLOJETLIKE::InitPDF() {
   //
   //  Initalize some necessary LHAPDF parameters
   //
}


//______________________________________________________________________________



vector<double> FastNLONLOJETLIKE::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   return LHAPDF::xfx(xp,muf);
}


#endif
