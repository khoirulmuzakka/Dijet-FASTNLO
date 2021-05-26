// Author: Daniel Britzger
// DESY, 20/04/2012

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO_toolkit                                                     //
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
//  fastNLOAlphas
//  This class inherits the PDF interface from
//  fastNLOLHAPDF, while the alpha_s evolution
//  is superseeded by the Alphas.h class.
//
//////////////////////////////////////////////////////////////////////////

#ifndef FASTNLOCRUNDEC
#define FASTNLOCRUNDEC

#include "fastNLOLHAPDF.h"
#include "CRunDec.h"



class fastNLOCRunDec : public fastNLOLHAPDF {

 public:
   // Parameter initialisation with PDG values
   fastNLOCRunDec(std::string name);
   // Parameter initialisation with LHAPDF set values
   fastNLOCRunDec(std::string name, std::string LHAPDFFile, int PDFMem = 0);

   // Getters
   double GetQMass(int pdgid) const;
   double GetMz() const;
   std::string GetNScheme() const;
   int GetNFlavor() const;
   int GetNLoop() const;
   double GetAlphasMz() const;

   // Setters
   void SetQMass(int pdgid, double qmass);
   void SetMz(double Mz);
   void SetNFlavor(int nflavor);
   void SetNLoop(int nloop);
   void SetAlphasMz(double AlphasMz);
   void SetPDGValues();
   void SetLHAPDFValues(std::string LHAPDFFile, int PDFMem = 0);

   // Printers
   void PrintParmValues();

 protected:
   // Inherited functions
   double EvolveAlphas(double Q) const ;

   // ---- Alphas vars ---- //
   void InitCRunDec();
   CRunDec *crundec;
   double QMass[6];
   double fMz;
   std::string fnScheme;
   int fnFlavor;
   int fnLoop;
   double fAlphasMz;
};

#endif
