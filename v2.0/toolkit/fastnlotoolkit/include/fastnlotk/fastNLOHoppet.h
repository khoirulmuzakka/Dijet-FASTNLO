// Author: Daniel Britzger
// DESY, 20/04/2012

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO_toolkit                                                     //
//  D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch           //
//                                                                      //
//  The projects web page can be found at:                              //
//    https://fastnlo.hepforge.org                                      //
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

#ifndef FASTNLOHOPPET
#define FASTNLOHOPPET

#include "fastNLOLHAPDF.h"



class fastNLOHoppet : public fastNLOLHAPDF {

 public:
   // Only allow constructor with information on LHAPDF set and member
   // as needed for HOPPET initialisation.
   fastNLOHoppet(std::string name) = delete;
   fastNLOHoppet(std::string name, std::string LHAPDFFile, int PDFSet);

   // Getters
   double GetQMass(int pdgid) const;
   int GetNFlavor() const;
   int GetNLoop() const;
   double GetMz() const;
   double GetAlphasMz() const;

   // Setters
   virtual bool InitPDF();
   void SetQMass(int pdgid, double qmass);
   void SetNFlavor(int nflavor);
   void SetNLoop(int nloop);
   void SetMz(double Mz);
   void SetAlphasMz(double AlphasMz, bool ReCalcCrossSection = false);
   void SetPDGValues();
   void SetLHAPDFValues(std::string LHAPDFFile);

   // Printers
   void PrintParmValues();

   protected:
   // Inherited functions
   virtual double EvolveAlphas(double Q) const;
   virtual std::vector<double> GetXFX(double xp, double muf) const;
};

#endif
