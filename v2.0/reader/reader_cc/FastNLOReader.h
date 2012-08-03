// Author: Daniel Britzger
// DESY, 23/07/2011

#ifndef FASTNLOREADER
#define FASTNLOREADER


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  FastNLOReader                                                       //
//                                                                      //
//  FastNLOReader is a standalone code for reading                      //
//  FastNLO tables of version 2.0 for DIS processes                     //
//  It is also optimized for an integration into                        //
//  the H1Fitter project.                                               //
//                                                                      //
//  FastNLO is developed by                                             //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch         //
//    (publication in preparation)                                      //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <string>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <vector>
#include "FastNLOBlockB.h"

using namespace std;

class speaker {
public:
   speaker(std::string prefix="",bool quiet=false,bool err=false) : 
      weg(0) , quiet(quiet){
      weg.clear(std::ios::badbit);
      pref=prefix;
      errs=cerr; };
   speaker(const speaker& spk):weg(0){;}; 
   std::ostream& operator() (std::string fct) {
      if (quiet) return weg;
      *this<<"In "<<fct<<". ";
      if(errs) return std::cerr;
      else return std::cout;
   }
   std::ostream& operator[] (std::string fct) {
      if (quiet) return weg;
      return *this<<"["<<fct<<"] ";
      if(errs) return std::cerr;
      else return std::cout;
      //return *this;
   }
   speaker& operator+ (std::string fct){return this->prefix(fct);}
   speaker& prefix(std::string fct) {
      if (!quiet){
	 if (errs) std::cerr<<fct;
	 else std::cout<<fct;;
      }
      return *this;
   }
   template<typename T> std::ostream& operator<< (T arg) {
      if (quiet) return weg;
      else {	
	 if (errs) return std::cerr<<pref<<arg;
	 else std::cout<<pref<<arg;
      }
   };
   template<typename T> std::ostream& operator>> (T arg) {
      return print(arg);
   };
   std::ostream& print(string mes) { 
      if (quiet) return weg;
      else {	
	 if (errs) return std::cerr<<mes;
	 else  return std::cout<<mes;
      }
   }
   void DoSpeak(bool loud){quiet=!loud;};
   bool GetSpeak() const {return !quiet;};
   void SetPrefix(std::string prefix){pref=prefix;};
   std::string GetPrefix(std::string prefix) const {return pref;};
private:
   std::ostream weg;
   bool quiet;
   std::string pref;
   bool errs;
};


class FastNLOReader {

public:

  typedef double(*mu_func)(double,double);

  enum EMuX {
    kMuR			= 0,	// renormalization scale
    kMuF			= 1	// factorization scale
  };
  
  enum EScaleFunctionalForm {
    kScale1			= 0,	// e.g. mu^2 = Q^2 
    kScale2			= 1,	// e.g. mu^2 = pt^2 
    kQuadraticSum		= 2,	// e.g. mu^2 = ( Q^2 + pt^2 )
    kQuadraticMean		= 3,	// e.g. mu^2 = ( Q^2 + pt^2 ) / 2 
    kQuadraticSumOver4		= 4,	// e.g. mu^2 = ( Q^2 + pt^2 ) / 4 
    kLinearMean			= 5,	// e.g. mu^2 = (( Q + pt ) / 2 )^2
    kLinearSum			= 6,	// e.g. mu^2 = (( Q + pt ))^2
    kScaleMax			= 7,	// e.g. mu^2 = max( Q^2, pt^2)
    kScaleMin			= 8,	// e.g. mu^2 = min( Q^2, pt^2) 
    kExpProd2                   = 9,    // e.g. mu^2 = (scale1 * exp(0.3 * scale2)) ^2
    kExtern			= 10	// define an external function for your scale
  };

  enum EAlphasEvolution {
    kGRV			= 0,
    kNLOJET			= 1,      
    kCTEQpdf			= 2,      
    kExternAs			= 4,	// use alphas-evolution as define in FastNLOInterface.cc
    kFixed			= 7     // Always gives back alpha_s(Mz) for testing.
  };

  enum EUnits {
    kAbsoluteUnits		= 0,	// calculate the cross section in barn for each publicated bin
    kPublicationUnits		= 1	// calculate the cross section in units as given in the according publication
  };

   enum Verbosity {
      DevelopDebug		= -100,	// All output, inclusive debug informations
      Debug			= -1,	// All output, also text-output with explanations
      Info			= 0,	// Info, Warning and error messages
      Warning			= 1,	// Warning and error messages
      Error			= 2,	// error messages
      Silent			= 3	// No output (not recommended)
   };

  // Corresponds to IContrFlag1 in v2 table definition
  enum ESMCalculation {
    kFixedOrder		        = 0,	// Fixed order calculation (pQCD)
    kThresholdCorrection	= 1,	// Threshold corrections
    kElectroWeakCorrection	= 2,	// Electroweak corrections
    kNonPerturbativeCorrection	= 3	// Non-perturbative corrections|Hadronisation corrections
  };
  
  // Corresponds to IContrFlag2 in v2 table definition
  enum ESMOrder {
    kLeading		        = 0,	// LO,   1-loop, LO MC
    kNextToLeading        	= 1,	// NLO,  2-loop, NLO MC
    kNextToNextToLeading	= 2	// NNLO, 3-loop, NNLO MC
  };
  
  static const double TWOPI = 6.28318530717958647692528;
  static const double TWOPISQR = 39.47841760435743447533796;

protected:

  static const int tablemagicno	= 1234567890;
  string ffilename;
  int fScalevar;
  double fScaleFacMuR;
  double fScaleFacMuF;
  EAlphasEvolution	fAlphasEvolution;
  EScaleFunctionalForm fMuRFunc;
  EScaleFunctionalForm fMuFFunc;
  EUnits		fUnits;
  mu_func Fct_MuR;				// Function, if you define your functional form for your scale external
  mu_func Fct_MuF;				// Function, if you define your functional form for your scale external
  vector < vector < bool > > bUseSMCalc;		// switch calculations ON/OFF
  vector < vector < bool > > bUseNewPhys;		// switch calculations ON/OFF

  // ---- alpha_s vars ---- //
  double fAlphasMz;

  // ---- Block A1 ---- //
  int Itabversion;
  string ScenName;
  int Ncontrib;
  int Nmult;
  int Ndata;
  int NuserString;
  int NuserInt;
  int NuserFloat;
  int Imachine;

  // ---- Block A2 ---- //
  int Ipublunits;
  vector < int > bla;
  vector <string> ScDescript;
  double Ecms;
  int ILOord;
  int NObsBin;
  int NDim;
  vector <int> RapIndex;
  vector <string> DimLabel;
  vector <int> IDiffBin;
  vector < vector <double> > LoBin;
  vector < vector <double> > UpBin;
  vector <double> BinSize;
  int INormFlag;
  string DenomTable;
  vector <int> IDivLoPointer;
  vector <int> IDivUpPointer;

  // ---- Block B ---- //
  FastNLOBlockB* BlockB_Data;
  FastNLOBlockB* BlockB_LO_Ref;
  FastNLOBlockB* BlockB_NLO_Ref;
  vector < vector < FastNLOBlockB* > > BBlocksSMCalc;	// BlockB's for SM corrections
  vector < vector < FastNLOBlockB* > > BBlocksNewPhys;	// BlockB's for New physics corrections

  // ---- Cross sections ---- //
  vector < double > XSection_LO;
  vector < double > XSection;
  vector < double > kFactor;

  // ----  reference tables ---- //
  vector < double > XSectionRef;
  vector < double > XSectionRefMixed;
  vector < double > XSectionRef_s1;
  vector < double > XSectionRef_s2;

 
public:

  FastNLOReader(string filename);
  FastNLOReader(const FastNLOReader& fnlo);
  virtual ~FastNLOReader(void);

   static void SetVerbosity(FastNLOReader::Verbosity verbosity);
  void SetFilename(string filename) ;
  void InitScalevariation();
  void SetAlphasMz( double AlphasMz , bool ReCalcCrossSection = false );
  void SetAlphasEvolution( EAlphasEvolution AlphasEvolution );
  void SetUnits( EUnits Unit );
  //void SetCalculationOrder( ECalculationOrder order ){ fOrder = order;};
  void SetContributionON( ESMCalculation eCalc , unsigned int Id , bool SetOn = true);	// Set contribution On/Off. Look for Id of this contribution during initialization.
  int ContrId( const ESMCalculation eCalc, const ESMOrder eOrder ) const;
  void SetGRVtoPDG2011_2loop(bool print);

  // ---- setters for scales of MuVar tables ---- //
  void SetMuRFunctionalForm( EScaleFunctionalForm func);// Set the functional form of Mu_R
  void SetMuFFunctionalForm( EScaleFunctionalForm func , bool ReFillCache = true);// Set the functional form of Mu_F
  void SetFunctionalForm( EScaleFunctionalForm func , FastNLOReader::EMuX kMuX);// Set functional form of MuX
  bool SetScaleFactorsMuRMuF( double xmur, double xmuf, bool ReFillCache = true);// Set scale factors for MuR and MuF
  void SetExternalFuncForMuR( mu_func);						// Set external function for scale calculation (optional)
  void SetExternalFuncForMuF( mu_func , bool ReFillCache = true);		// Set external function for scale calculation (optional)


  // ---- Pdf interface ---- //
  void FillPDFCache( bool ReCalcCrossSection = false );					// Prepare for recalculation of cross section with 'new'/updated pdf.


  // ---- Do the cross section calculation ---- //
  void CalcCrossSection();


  // ---- Getters for results---- //
  vector < double > GetCrossSection();
  vector < double > GetReferenceCrossSection();
  vector < double > GetKFactors();


  // ---- Getters for FastNLOReader member variables ---- //
  EScaleFunctionalForm GetMuRFunctionalForm() const { return fMuRFunc; };
  EScaleFunctionalForm GetMuFFunctionalForm() const { return fMuFFunc; };
  EAlphasEvolution GetAlphasEvolution() const { return fAlphasEvolution; };
  EUnits GetUnits() const{ return fUnits; };
  mu_func GetExternalFuncForMuR(){ return Fct_MuR; };
  mu_func GetExternalFuncForMuF(){ return Fct_MuF; };
  double GetAlphasMz() const { return fAlphasMz; };
  double GetScaleFactorMuR() const { return fScaleFacMuR;};
  double GetScaleFactorMuF() const { return fScaleFacMuF;};
  int GetScaleVariation() const { return fScalevar; };

  // ---- Getters for FastNLO table constants ---- //
  int GetNcontrib() const { return Ncontrib; };
  int GetIExpUnit() const { return Ipublunits; };			// exponent of xs units (like -12 for pb)
  string GetScenarioName() const { return ScenName; };			// Get Scenario/Table name
  vector < string > GetScenarioDescription() const { return ScDescript; };	// Get Description of scenario
  double GetCMSEnergy() const { return Ecms; };				// Get center of mass energy
  int GetILOord() const { return ILOord; };				// Get number of alpha_s in leading order (1 or 2 usually)
  int GetNObsBins() const { return NObsBin; };				// Get number of measured bins
  int GetNDiffBin() const { return NDim; };				// Get number of differential measurement. 1: single differential; 2: double differential
  vector < int > GetRapidityIndex() const { return RapIndex;};		// Get rapidity indices
  vector < string > GetDimensionLabel() const { return DimLabel;};	// Get label for measured dimensions
  vector < int > GetIDiffBin() const { return IDiffBin;};		// Get number of differential bins
  vector < vector < double > > GetLowBinEdge() const { return LoBin; };	// Get Lower Bin edge [ObsBin][DiffBin]
  vector < vector < double > > GetUpBinEdge() const { return UpBin; };	// Get Upper Bin edge [ObsBin][DiffBin]
  vector < double > GetBinSize() const { return BinSize; };		// Get Binsize = BinSizeDim1 < * BinSizeDim2 >
  int IsNormalized() const { return INormFlag; };			// Normalized Cross sections?
  string GetScaleDescription(int scalen=0) const { return BBlocksSMCalc[0][0]->ScaleDescript[0][scalen]; };		// Description of renormalization and facorization scale choice
  int GetNScaleVariations() const;					// Get number of available scale variations
  vector < double > GetScaleFactors() const;				// Get list of available scale factors
   bool GetIsFlexibleScaleTable() const { return BBlocksSMCalc[0][0]->NScaleDep == 3; } // Get, if this table is a 'flexible scale' table or not.


  // ---- Print outs ---- //
  void PrintTableInfo(const int iprint = 0) const;
  void PrintFastNLOTableConstants(const int iprint = 2) const;
  void PrintCrossSections() const ;
  void PrintCrossSectionsDefault(vector<double> kthc = vector<double>() ) const ;
  void PrintCrossSectionsWithReference();
  void PrintCrossSectionsData() const;
  void PrintFastNLODemo();


   static speaker debug;
   static speaker error;
   static speaker warn;
   static speaker info;
   static speaker text;
      

protected:

  void Init() ;
  void ReadTable();
  void InitMembers();
  void StripWhitespace(string* s);

  void ReadBlockA1(istream *table);
  void ReadBlockA2(istream *table);
  void ReadBlockB(istream *table);

  void PrintBlockA1() const;
  void PrintBlockA2() const;

  void PrintScaleSettings(EMuX kMuX=kMuR);
  void FillBlockBPDFLCsDISv20( FastNLOBlockB* B );
  void FillBlockBPDFLCsDISv21( FastNLOBlockB* B );
  void FillBlockBPDFLCsHHCv20( FastNLOBlockB* B );
  void FillBlockBPDFLCsHHCv21( FastNLOBlockB* B );
  void CalcAposterioriScaleVariation();
  vector<double> CalcPDFLinearCombDIS(vector<double> pdfx1, int NSubproc );
  vector<double> CalcPDFLinearCombHHC(vector<double> pdfx1, vector<double> pdfx2, int NSubproc );
  void FillAlphasCacheInBlockBv20( FastNLOBlockB* B );
  void FillAlphasCacheInBlockBv21( FastNLOBlockB* B );
  double CalcAlphas(double Q);

  void CalcReferenceCrossSection();

  double CalcAlphasNLOJET(double Q, double alphasMz);
  double CalcAlphasCTEQpdf(double Q, double alphasMz);
  
  double CalcMu(FastNLOReader::EMuX kMuX, double scale1 , double scale2 , double scalefactor);
  double FuncMixedOver1 ( double scale1 , double scale2 ) ;
  double FuncMixedOver2 ( double scale1 , double scale2 ) ;
  double FuncMixedOver4 ( double scale1 , double scale2 ) ;
  double FuncLinearMean ( double scale1 , double scale2 ) ;
  double FuncLinearSum ( double scale1 , double scale2 ) ;
  double FuncMax ( double scale1 , double scale2 ) ;
  double FuncMin ( double scale1 , double scale2 ) ;
  double FuncExpProd2 ( double scale1 , double scale2 ) ;

  void CalcCrossSectionv21(FastNLOBlockB* B , bool IsLO = false );
  void CalcCrossSectionv20(FastNLOBlockB* B , bool IsLO = false);
   
   FastNLOBlockB* B_NLO() { return BBlocksSMCalc[0][1]; };
   FastNLOBlockB* B_LO() { return BBlocksSMCalc[0][0]; };
   FastNLOBlockB* B_ThC(int n=0) { 
      if ( BBlocksSMCalc[kThresholdCorrection].empty() ) return NULL;
      else return BBlocksSMCalc[kThresholdCorrection][n]; };

   // virtual functions for the user interface
   virtual void InitPDF() = 0;
   virtual vector<double> GetXFX(double x, double muf ) const = 0;
   virtual double EvolveAlphas(double Q, double alphasMz) const = 0;

   // ---- setters for scale variation in v2.0 tables  ---- //
   double SetScaleVariation( int scalevar , bool ReFillCache = true);// Choose the MuF scale variation table
   // ---- alphas cache ---- //
   void FillAlphasCache();								// prepare for recalculation of cross section with new alpha_s value.

  // ---- human readable strings ---- //
  static const string fContrName[20];
  static const string fOrdName[4][4];
  static const string fNSDep[4];

protected:
  static int WelcomeOnce;
   
};


#endif
