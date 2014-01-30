// Daniel Britzger
// DESY, 08.08.2013
#ifndef __fastNLOTable__
#define __fastNLOTable__

#include <math.h>
#include <string>
#include <fstream>
#include <vector>
#include <utility>

#include "fastNLOConstants.h"
#include "fastNLOBase.h"
#include "fastNLOCoeffBase.h"

#include "fastNLOCoeffData.h"
#include "fastNLOCoeffMult.h"
#include "fastNLOCoeffAddFix.h"
#include "fastNLOCoeffAddFlex.h"

using namespace std;

class fastNLOTable : public fastNLOBase {

 public:
   fastNLOTable();
   fastNLOTable(string filename);
   ~fastNLOTable();
   fastNLOTable(const fastNLOTable&);

   virtual void Print() const;
   void PrintScenario() const;
   void PrintFastNLOTableConstants(const int iprint = 2) const;                                 //  Print (technical) constants of fastNLO table (use iprint) for level of details.
   void PrintTableInfo(const int iprint = 0) const;                                             //  Print basic info about fastNLO table and its contributions

   int ReadTable();
   void WriteTable();
   void WriteTable(string filename);
   bool IsCompatible(const fastNLOTable& other) const;

   int GetNObsBin() const {return NObsBin;}

   double GetLoBin(int bin, int dimension) const {return Bin[bin][dimension].first;}            // Get lower bin boundary
   double GetUpBin(int bin, int dimension) const {return Bin[bin][dimension].second;}           // Get upper bin boundary

   vector < double > GetBinSize() const {return BinSize;};					// Get Binsize = BinSizeDim1 < * BinSizeDim2 >
   double GetBinSize(int bin) const {return BinSize[bin];};					// Get Binsize = BinSizeDim1 < * BinSizeDim2 >

   void SetNumDiffBin(int iDiff ) { NDim=iDiff; DimLabel.resize(NDim); IDiffBin.resize(NDim);}  // Set dimension of calculation. (Singledifferential, double-differntial, etc...)
   int GetNumDiffBin() const { return NDim; }                                                   // Get dimension of calculation. (Singledifferential, double-differntial, etc...)

   int GetIDiffBin(int bin) const { return IDiffBin[bin]; }					// Get if dimension is 'truly differential' or bin-integrated (divided by bin-width or not)

   void SetDimLabel( string label, int iDim , bool IsDiff = true );
   string GetDimLabel( int iDim  ) const {return DimLabel[iDim];};				// Get label (name) of observable in dimension iDim
   vector<string > GetDimLabels() const {return DimLabel;};					// Get label (name) of all observables

   void SetIpublunits(int unit){Ipublunits = unit;}
   int GetIpublunits() const {return Ipublunits;}

   double GetEcms() const {return Ecms;}
   void SetEcms(double E) {Ecms = E;}

   int GetLoOrder() const {return ILOord;}
   void SetLoOrder(int LOOrd);

   bool IsNorm() const { return INormFlag == 0 ? false : true;}
   string GetDenomTable() const {return DenomTable;}

   fastNLOCoeffBase* GetCoeffTable(int no) const;

   // useful functions
   //    void InitBinning( const int nBins1 , double* bingrid1 , const int* nBins2 = NULL , vector<double*> bingrid2 = vector<double*>() , double binwidth3 = 0 );
   //    void InitBinningKR( const int nBins1 , const double* bingrid1 , const int* nBins2 = NULL , vector< vector<double> > bingrid2 = vector< vector<double> >() , const double bwfactor = 0. );
   int GetBinNumber(double val1 , double val2 = -42. ) const ;                                  // calculate bin number (iObsBin)

   // handle coefficient tables
   //int WriteCoeffTable(int no);
   //int WriteCoeffTable(int no,ofstream* outstream );
   //int WriteCoeffTableDividebyN(int no);
   void DeleteAllCoeffTable();
   //int CreateCoeffBase(int no);
   int CreateCoeffTable(int no,fastNLOCoeffBase *newcoeff);
   void AddTable(const fastNLOTable& rhs);


private:
   bool cmp(const double x1, const double x2) const;
   bool cmp(const vector<double>& x1, const vector<double >& x2) const;
   bool cmp(const vector<vector<double> >& x1,const vector<vector<double > >& x2) const;
   bool cmp(const vector<vector<pair<double,double> > >&  x1,const vector<vector<pair<double,double> > >& x2) const;


protected:
   int WriteScenario(ostream *table);
   int ReadScenario(istream *table);
   int ReadCoeffTables(istream *table);
   fastNLOCoeffBase* ReadRestOfCoeffTable(const fastNLOCoeffBase& cB, istream *table);

   vector < fastNLOCoeffBase* > fCoeff;
   //fastNLOCoeffData* fData;

   double Ecms;
   int ILOord;
   int Ipublunits;
   vector <string> ScDescript;

   int NObsBin;
   int NDim;
   vector <string> DimLabel;
   vector <int> IDiffBin;
   vector < vector <pair<double,double> > > Bin; // every bin has a lower and upper bin boundary and belongs to a 'dimension'. If a truely differential measurment, then upper bin boundary is equal lower one
   vector <double> BinSize;

   vector <int> RapIndex;               //KR: Added possibility to store and read start of new rapidity bin in nobs

   // contributions for normalization
   int INormFlag;
   string DenomTable;
   vector <int> IDivLoPointer;
   vector <int> IDivUpPointer;

};
#endif
