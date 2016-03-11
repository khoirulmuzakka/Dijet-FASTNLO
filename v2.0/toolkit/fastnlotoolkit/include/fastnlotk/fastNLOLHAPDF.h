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

#ifndef FASTNLOLHAPDF
#define FASTNLOLHAPDF

#include "fastNLOReader.h"
#include "fastNLOConstants.h"
#include <LHAPDF/LHAPDF.h>
#include <cmath>


class fastNLOLHAPDF : public fastNLOReader {

private:
public:
   fastNLOLHAPDF(std::string name);
   fastNLOLHAPDF(const fastNLOTable&);
   ~fastNLOLHAPDF();
   fastNLOLHAPDF(std::string name, std::string LHAPDFfile, int PDFSet = 0);
   fastNLOLHAPDF(const fastNLOTable&, std::string LHAPDFfile, int PDFSet = 0);

   // Initializer. Necessary for some alternative evolutions.
   virtual void InitEvolveAlphas();
   // Pseudo-Setters. Don´t work with LHAPDF, but print warning instead.
   virtual void SetMz(double Mz);
   virtual void SetNFlavor(int nflavor);
   virtual void SetNLoop(int nloop);
   virtual void SetAlphasMz(double AlphasMz , bool ReCalcCrossSection = false);
   virtual void SetQMass(int pdgid, double mq);
   // Setters
   void SetLHAPDFFilename(std::string filename);
   void SetLHAPDFMember(int set);
   // Getters
   std::string GetLHAPDFFilename() const {return fLHAPDFFilename;}
   int GetIPDFMember() const;
   int GetNPDFMembers() const;
   int GetNPDFMaxMember() const;
   void PrintPDFInformation() const ;
   virtual double GetQMass(int pdgid) const;
   int GetNLoop() const;
   int GetNFlavor() const;
   double GetAlphasMz() const;
   //! Return struct with vectors containing the cross section values and the selected scale uncertainty
   XsUncertainty GetPDFUncertainty( const fastNLO::EPDFUncertaintyStyle ePDFUnc );
   XsUncertainty GetPDFUncertainty( const fastNLO::EPDFUncertaintyStyle ePDFUnc, bool lNorm );
   std::vector<std::vector<double> > GetPDFUncertaintyVec(fastNLO::EPDFUncertaintyStyle);
   // Deprecated: Replaced by struct as return object: Return vector of pairs with all cross section values first and pairs of PDF uncertainties second
   //   vector < pair < double, pair <double, double> > > GetPDFUncertainty(const EPDFUncertaintyStyle ePDFUnc);
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   std::vector<LHAPDF::PDFUncertainty>  GetPDFUncertaintyLHAPDF(double cl=100*erf(1/sqrt(2)), bool alternative=false); //!< return PDF uncertainty, formulae taken from LHAPDF6
   std::vector<double> CalcPDFUncertaintyMinus(const std::vector<LHAPDF::PDFUncertainty>& ) const; //!<get vector<double> for PDF-minus uncertainty. Uncertainties are POSITIVE!
   std::vector<double> CalcPDFUncertaintyPlus(const std::vector<LHAPDF::PDFUncertainty>& ) const; //!<get vector<double> for PDF-up uncertainty
   std::vector<double> CalcPDFUncertaintyRelMinus(const std::vector<LHAPDF::PDFUncertainty>& ) const; //!<get vector<double> for relative PDF-minus uncertainty. Uncertainties are NEGATIVE!
   std::vector<double> CalcPDFUncertaintyRelPlus(const std::vector<LHAPDF::PDFUncertainty>& ) const; //!<get vector<double> for relative PDF-up uncertainty
   std::vector<double> CalcPDFUncertaintySymm(const std::vector<LHAPDF::PDFUncertainty>& ) const; //!< get vector<double> for symmetrized PDF uncertainty
   std::vector<double> CalcPDFUncertaintyCentral(const std::vector<LHAPDF::PDFUncertainty>& ) const; //!< get vector<double> for 'new' central value
#endif

protected:
   // inherited functions
   double EvolveAlphas(double Q) const ;
   virtual bool InitPDF();
   std::vector<double> GetXFX(double xp, double muf) const ;

   // ---- LHAPDF vars ---- //
   std::string fLHAPDFFilename;
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   LHAPDF::PDFSet* PDFSet;
   LHAPDF::PDF* PDF;
   #endif
   int fnPDFs;
   int fiPDFMember;

   double fchksum;


};

#endif
