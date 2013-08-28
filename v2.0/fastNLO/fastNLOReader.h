#ifndef __fastNLOReader__
#define __fastNLOReader__

#include "fastNLOTable.h"
#include "fastNLOPDFLinearCombinations.h"

using namespace std;

class fastNLOReader : public fastNLOTable , public fastNLOPDFLinearCombinations {
   //
   // fastNLOReader. 
   //

public:
   typedef double(*mu_func)(double,double);

   fastNLOReader(string filename);
   virtual ~fastNLOReader();
   void SetFilename(string filename) ;
   void InitScalevariation();
   void SetUnits(fastNLO::EUnits Unit);
   void SetContributionON(fastNLO::ESMCalculation eCalc , unsigned int Id , bool SetOn = true);  // Set contribution On/Off. Look for Id of this contribution during initialization.
   int ContrId(const fastNLO::ESMCalculation eCalc, const fastNLO::ESMOrder eOrder) const;

   // ---- setters for scales of MuVar tables ---- //
   void SetMuRFunctionalForm(fastNLO::EScaleFunctionalForm func); // Set the functional form of Mu_R
   void SetMuFFunctionalForm(fastNLO::EScaleFunctionalForm func);  // Set the functional form of Mu_F
   void SetFunctionalForm(fastNLO::EScaleFunctionalForm func , fastNLO::EMuX kMuX); // Set functional form of MuX
   bool SetScaleFactorsMuRMuF(double xmur, double xmuf);			// Set scale factors for MuR and MuF
   void SetExternalFuncForMuR(mu_func);                                         // Set external function for scale calculation (optional)
   void SetExternalFuncForMuF(mu_func);						// Set external function for scale calculation (optional)


   // ---- Pdf interface ---- //
   void FillPDFCache(double chksum=0.);                                          // Prepare for recalculation of cross section with 'new'/updated pdf.

   // ---- alphas cache ---- //
   void FillAlphasCache();                                                              // prepare for recalculation of cross section with new alpha_s value.

   // ---- Do the cross section calculation ---- //
   void CalcCrossSection();


   // ---- Getters for results---- //
   vector < double > GetCrossSection();
   vector < double > GetReferenceCrossSection();
   vector < double > GetKFactors();
   vector < double > GetQScales(int irelord); // Order (power of alpha_s) rel. to LO: 0 --> LO, 1 --> NLO

   // ---- Getters for FastNLOReader member variables ---- //
   fastNLO::EScaleFunctionalForm GetMuRFunctionalForm() const {
      return fMuRFunc;
   };
   fastNLO::EScaleFunctionalForm GetMuFFunctionalForm() const {
      return fMuFFunc;
   };
   fastNLO::EUnits GetUnits() const {
      return fUnits;
   };
   mu_func GetExternalFuncForMuR() {
      return Fct_MuR;
   };
   mu_func GetExternalFuncForMuF() {
      return Fct_MuF;
   };
   double GetScaleFactorMuR() const {
      return fScaleFacMuR;
   };
   double GetScaleFactorMuF() const {
      return fScaleFacMuF;
   };
   int GetScaleVariation() const {
      return fScalevar;
   };

   // ---- Getters for FastNLO table constants ---- //
   //    int GetIExpUnit() const {						// exponent of xs units (like -12 for pb)
   //       return Ipublunits;
   //    };                       
   //    string GetScaleDescription(int scalen=0) const {			// Description of renormalization and facorization scale choice
   //       return BBlocksSMCalc[0][0]->ScaleDescript[0][scalen];
   //    }; 
   
   int GetNScaleVariations() const;						// Get number of available scale variations
   vector < double > GetScaleFactors() const;					// Get list of available scale factors
   bool GetIsFlexibleScaleTable(fastNLOCoeffAddBase* c=NULL) const {		// Get, if this table is a 'flexible scale' table or not.
      if ( c ) return  c->GetNScaleDep() >= 3;
      else return BBlocksSMCalc[0][0]->GetNScaleDep() >= 3;
   }

   // ---- Print outs ---- //
   void PrintTableInfo(const int iprint = 0) const;				//  Print basic info about fastNLO table and its contributions
   void PrintFastNLOTableConstants(const int iprint = 2) const;			//  Print (technical) constants of fastNLO table (use iprint) for level of details.
   void PrintCrossSections() const;						//  Print cross sections (optimized for double-differential tables)
   void PrintCrossSectionsDefault(vector<double> kthc = vector<double>()) const;//  Print cross sections in the same format as in the fortran version.
   void PrintCrossSectionsWithReference();					
   void PrintCrossSectionsData() const;						//  Print data table. (if available)

   void RunFastNLODemo();							//  Run an example of fastNLO for educational purposes, i.e. calculate and print cross sections for several scale variations
 



   // ---- Test virtual functions for reasonable values. ---- //
   bool TestXFX();								// Test if XFX reasonable values
   bool TestAlphas();								// Test if EvolvaAlphas returns a reasonable value


protected:
   fastNLOReader();
   void Init() ;
   //void ReadTable();
   void StripWhitespace(string* s);

   void PrintScaleSettings(EMuX kMuX=kMuR);
   void FillBlockBPDFLCsDISv20(fastNLOCoeffAddPert* B);
   void FillBlockBPDFLCsDISv21(fastNLOCoeffAddFlex* B);
   void FillBlockBPDFLCsHHCv20(fastNLOCoeffAddPert* B);
   void FillBlockBPDFLCsHHCv21(fastNLOCoeffAddFlex* B);
   void CalcAposterioriScaleVariation();
   void FillAlphasCacheInBlockBv20(fastNLOCoeffAddPert* B);
   void FillAlphasCacheInBlockBv21(fastNLOCoeffAddFlex* B);
   double CalcAlphas(double Q);
   double CalcReferenceAlphas();
   double CalcNewPDFChecksum();
   double CalcChecksum(double mu);
   bool PrepareCache();

   void CalcReferenceCrossSection();

   double CalcMu(fastNLO::EMuX kMuX, double scale1 , double scale2 , double scalefactor);
   double FuncMixedOver1(double scale1 , double scale2) ;
   double FuncMixedOver2(double scale1 , double scale2) ;
   double FuncMixedOver4(double scale1 , double scale2) ;
   double FuncLinearMean(double scale1 , double scale2) ;
   double FuncLinearSum(double scale1 , double scale2) ;
   double FuncMax(double scale1 , double scale2) ;
   double FuncMin(double scale1 , double scale2) ;
   double FuncExpProd2(double scale1 , double scale2) ;

   void CalcCrossSectionv21(fastNLOCoeffAddFlex* B , bool IsLO = false);
   void CalcCrossSectionv20(fastNLOCoeffAddPert* B , bool IsLO = false);

   fastNLOCoeffAddBase* B_NLO() const {
      return (fastNLOCoeffAddBase*) BBlocksSMCalc[fastNLO::kFixedOrder][fastNLO::kNextToLeading];
   };
   fastNLOCoeffAddBase* B_LO() const {
      return (fastNLOCoeffAddBase*) BBlocksSMCalc[fastNLO::kFixedOrder][fastNLO::kLeading];
   };
   fastNLOCoeffBase* B_ThC(int n=0) {
      if (BBlocksSMCalc[fastNLO::kThresholdCorrection].empty()) return NULL;
      else return BBlocksSMCalc[fastNLO::kThresholdCorrection][n];
   };

   // virtual functions for the user interface
   virtual bool InitPDF() = 0;
   virtual vector<double> GetXFX(double x, double muf) const = 0;
   virtual double EvolveAlphas(double Q) const = 0;

   // ---- setters for scale variation in v2.0 tables  ---- //
   double SetScaleVariation(int scalevar , bool FirstCall=false); // Choose the MuF scale variation table

   // ---- human readable strings ---- //
   static const string fContrName[20];
   static const string fOrdName[4][4];
   static const string fNSDep[6];


protected:
   string ffilename;
   int fScalevar;
   double fScaleFacMuR;
   double fScaleFacMuF;
   fastNLO::EScaleFunctionalForm fMuRFunc;
   fastNLO::EScaleFunctionalForm fMuFFunc;
   fastNLO::EUnits               fUnits;
   bool fPDFSuccess;
   double fPDFCached;
   double fAlphasCached;
   mu_func Fct_MuR;                              // Function, if you define your functional form for your scale external
   mu_func Fct_MuF;                              // Function, if you define your functional form for your scale external
   vector < vector < bool > > bUseSMCalc;                // switch calculations ON/OFF
   vector < vector < bool > > bUseNewPhys;               // switch calculations ON/OFF


   // ---- Block B ---- //
   fastNLOCoeffData* fCoeffData;
   fastNLOCoeffAddBase* Coeff_LO_Ref;
   fastNLOCoeffAddBase* Coeff_NLO_Ref;
//    vector< vector < fastNLOCoeffAddBase* > > fCoAdd;
//    vector< vector < fastNLOCoeffMult* > > fCoMult;
   vector < vector < fastNLOCoeffBase* > > BBlocksSMCalc;   // BlockB's for SM corrections
   vector < vector < fastNLOCoeffBase* > > BBlocksNewPhys;  // BlockB's for New physics corrections

   //    FastNLOBlockB* BlockB_Data;
   //    FastNLOBlockB* BlockB_LO_Ref;
   //    FastNLOBlockB* BlockB_NLO_Ref;
   //    vector < vector < FastNLOBlockB* > > BBlocksSMCalc;   // BlockB's for SM corrections
   //    vector < vector < FastNLOBlockB* > > BBlocksNewPhys;  // BlockB's for New physics corrections

   // ---- Cross sections ---- //
   vector < double > XSection_LO;
   vector < double > XSection;
   vector < double > kFactor;
   vector < double > QScale_LO;
   vector < double > QScale;

   // ----  reference tables ---- //
   vector < double > XSectionRef;
   vector < double > XSectionRefMixed;
   vector < double > XSectionRef_s1;
   vector < double > XSectionRef_s2;



//    fastNLOCoeffData* fCoData;
//    vector<fastNLOCoeffMult* > fCoMult
//   vector<fastNLOCoeffAddEvalBase* > fCoAdd

};
#endif
