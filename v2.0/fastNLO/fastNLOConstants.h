#ifndef __fnloconstants__
#define __fnloconstants__

#include <string>
#include <vector>

#include "fastNLOVersion.cc"
#ifndef _FNLO_VERSION
#define _FNLO_VERSION    "2.1.0"
#define _FNLO_SVNREV     "xxxx"
#define _FNLO_AUTHORS    "D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch"
#define _FNLO_WEBPAGE    "http://projects.hepforge.org/fastnlo"
#define _FNLO_AUTHORSv14 "T. Kluge, K. Rabbertz, M. Wobisch"
#define _FNLO_QUOTEv14   "hep-ph/0609285"
#define _FNLO_AUTHORSv2  "D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch"
#define _FNLO_QUOTEv2    "arXiv:1109.1310"
#endif

namespace fastNLO {
  
   // ---- typedefs ---- //
   typedef std::vector<double > v1d;
   typedef std::vector<std::vector<double > > v2d;
   typedef std::vector<std::vector<std::vector<double > > > v3d;
   typedef std::vector<std::vector<std::vector<std::vector<double > > > > v4d;
   typedef std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > v5d;
   typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > > v6d;
   typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double > > > > > > > v7d;
 
   // ---- constants ---- //
   const int tabversion   = 20000;
   const int tablemagicno = 1234567890;
   const double TWOPI = 6.28318530717958647692528;


   // ---- enumerators ---- //
   enum EMuX {
      kMuR                        = 0,    // renormalization scale
      kMuF                        = 1     // factorization scale
   };

   enum EScaleFunctionalForm {
      kScale1                   = 0,    // e.g. mu^2 = Q^2
      kScale2                   = 1,    // e.g. mu^2 = pt^2
      kQuadraticSum             = 2,    // e.g. mu^2 = ( Q^2 + pt^2 )
      kQuadraticMean            = 3,    // e.g. mu^2 = ( Q^2 + pt^2 ) / 2
      kQuadraticSumOver4        = 4,    // e.g. mu^2 = ( Q^2 + pt^2 ) / 4
      kLinearMean               = 5,    // e.g. mu^2 = (( Q + pt ) / 2 )^2
      kLinearSum                = 6,    // e.g. mu^2 = (( Q + pt ))^2
      kScaleMax                 = 7,    // e.g. mu^2 = max( Q^2, pt^2)
      kScaleMin                 = 8,    // e.g. mu^2 = min( Q^2, pt^2)
      kExpProd2                 = 9,    // e.g. mu^2 = (scale1 * exp(0.3 * scale2)) ^2
      kExtern                   = 10    // define an external function for your scale
   };

   enum ESMCalculation {
      kFixedOrder               = 0,    // Fixed order calculation (pQCD)
      kThresholdCorrection      = 1,    // Threshold corrections
      kElectroWeakCorrection    = 2,    // Electroweak corrections
      kNonPerturbativeCorrection= 3     // Non-perturbative corrections|Hadronisation corrections
   };

   enum ESMOrder {
      kLeading                  = 0,    // LO,   1-loop, LO MC
      kNextToLeading            = 1,    // NLO,  2-loop, NLO MC
      kNextToNextToLeading      = 2     // NNLO, 3-loop, NNLO MC
   };

   enum EUnits {
      kAbsoluteUnits            = 0,    // calculate the cross section in barn for each publicated bin
      kPublicationUnits         = 1     // calculate the cross section in units as given in the according publication
   };


   // ---- some names for nice output ---- //
   const std::string _ContrName[20] = {
      "Fixed order calculation", "Threshold corrections", "Electroweak corrections", "Non-perturbative corrections",
      "Undefined", "Undefined", "Undefined", "Undefined", "Undefined", "Undefined", "Undefined",
      "Quark compositeness", "ADD-LED", "TeV 1-ED", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"
   };
   const std::string _OrdName[4][4] = { 
      { "LO",     "NLO",    "NNLO"   , "N3LO"    },
      { "1-loop", "2-loop", "3-loop" , "4-loop"  },
      { "Undef" , "Undef" , "Undef"  , "Undef"   },
      { "LO MC" , "NLO MC", "NNLO MC", "N3LO MC" }
   };
   const std::string _fNSDep[6] = {"v2.0","v2.0","v2.0","v2.1","v2.2","v2.2"};


   // ---- Some shapes for nice output ---- //
   const std::string _CSEP41("#########################################");
   const std::string _DSEP41("=========================================");
   const std::string _SSEP41("-----------------------------------------");
   const std::string _CSEP = _CSEP41 + _CSEP41 + _CSEP41 + _CSEP41;
   const std::string _DSEP = _DSEP41 + _DSEP41 + _DSEP41 + _DSEP41;
   const std::string _SSEP = _SSEP41 + _SSEP41 + _SSEP41 + _SSEP41;

}


#endif
