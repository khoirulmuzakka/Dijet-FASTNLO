// Author: Daniel Britzger
// DESY, 02/04/2012

#ifndef fASTNLODIFFREADER
#define fASTNLODIFFREADER


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLODiffReader                                                   //
//                                                                      //
//  fastNLODiffReader is a standalone code for reading                  //
//  diffractive FastNLO tables of version 2.0 for DIS processes         //
//                                                                      //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <string>
#include <vector>
#include "fastnlotk/fastNLOReader.h"


class fastNLODiffReader : public fastNLOReader {

public:

   fastNLODiffReader(std::string filename);
   virtual ~fastNLODiffReader(void) {};

   void SetXPomSlicing(int nStep, double* xpom, double* dxpom);
   void SetXPomLogSlicing(int nStep, double xpommin, double xpommax);
   void SetXPomLinSlicing(int nStep, double xpommin, double xpommax);
   void SetXPomExpSlicing(int nStep, double xpommin, double xpommax);
   void SetZRange(double zmin , double zmax) {
      fzmin = zmin ;
      fzmax = zmax;
   };
   double GetZRangeMin() { return fzmin; };
   double GetZRangeMax() { return fzmax; };
   void SetPrintCrossSectionVsxIPzIP(bool print=true) {fPrintxIPzIP=print;}; //!< Print these values during integration.
   std::vector < std::map< std::pair<double, double>, double > > GetXSection_vs_xIPzIP() { 
      return fXSection_vs_xIPzIP;}//!< Get cross section for each zIP node and xIP slice.

   std::vector < double > GetCrossSection();
   void CalcCrossSection();
   std::vector < double > GetDiffCrossSection();
   void FillPDFCache(bool ReCalcCrossSection = false);
   std::vector < double > GetReferenceCrossSection();

   // ---- Print outs must be overwritten ---- //
   void PrintCrossSectionsWithReference();

	std::map<double,  std::vector < std::map< double, double > > > Get3DCrossSection() { return xsQ2; }

	void SetProtonE(double E) { fProtonE = E; }

protected:

   double fxpom;
   double fzmin;
   double fzmax;

   std::vector < double > fxPoms;
   std::vector < double > fdxPoms;
   bool fPrintxIPzIP = false;
   std::vector < std::map< std::pair<double, double>, double > > fXSection_vs_xIPzIP; //! Cross section vs. x ( XSection_vsX1[bin][<{x,z},xs>] )

   // inherited functions
   virtual double EvolveAlphas(double Q) const = 0;
   virtual bool InitPDF() = 0;
   std::vector<double> GetXFX(double xp, double muf) const;
   virtual std::vector<double> GetDiffXFX(double xpom, double zpom, double muf) const = 0;


	std::map<double,  std::vector < std::map< double, double > > > xsQ2;
	double fProtonE;


};

#endif

