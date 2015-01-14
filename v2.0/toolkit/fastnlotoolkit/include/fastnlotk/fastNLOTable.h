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
   void PrintFastNLOTableConstants(const int iprint = 0) const;                                 //  Print (technical) constants of fastNLO table (use iprint) for level of details.
   void PrintTableInfo(const int iprint = 0) const;                                             //  Print basic info about fastNLO table and its contributions

   virtual void ReadTable();
   virtual void WriteTable();
   virtual void WriteTable(string filename);
   bool IsCompatible(const fastNLOTable& other) const;

   std::string GetRivetId() const;

   // getters for 'ObsBin': ObsBin runs from 0-NObsBin
   unsigned int GetNObsBin() const {return NObsBin;}                                                     // Number of observable bins
   std::vector < std::pair < double, double > > GetObsBinBoundaries(int iObsBin) const { return Bin[iObsBin];}    // Get binning for a given ObsBin

   // getters for dimension specific bins, here called Dim<I>Bins
   unsigned int GetIDim0Bin(unsigned int iobs) const;                                           // Bin number in first dimension for given ObsBin
   unsigned int GetNDim0Bins() const;                                                           // Bin number +1 in first dimension for last ObsBin
                                                                                                // ==> Number of bins in first dimension
   unsigned int GetIDim1Bin(unsigned int iobs) const;                                           // Bin number in second dimension for given ObsBin
   unsigned int GetNDim1Bins(unsigned int iDim0Bin) const;                                      // Number of bins in second dimension for given bin in first dimension
   unsigned int GetIDim2Bin(unsigned int iobs) const;                                           // Bin number in third dimension for given ObsBin
   unsigned int GetNDim2Bins(unsigned int iDim0Bin, unsigned int iDim1Bin) const;               // Number of bins in third dimension for a second and f

   unsigned int GetIDimBin(unsigned int iobs, unsigned int iDim) const;                         // Bin number in third dimension for given ObsBin

   int GetODim0Bin(double var0) const;
   int GetODim1Bin(double var0, double var1) const;
   int GetODim2Bin(double var0, double var1, double var2) const;
   std::vector < std::pair < double, double > > GetBinBoundaries(int iDim0Bin, int iDim1Bin = -1, int iDim2Bin = -1);

   int GetObsBinNumber( const vector < double >& vobs ) const ;                                 // Calculate observable bin number (iObsBin)
   int GetObsBinNumber( double obs0 ) const ;
   int GetObsBinNumber( double obs0, double obs1 ) const ;
   int GetObsBinNumber( double obs0, double obs1, double obs2 ) const ;

   std::pair < double, double > GetObsDimBin(unsigned int iobs, unsigned int iDim) const        // Get binning for given observable and dimension
   {return Bin[iobs][iDim];}
   std::vector < std::pair < double, double > > GetDimBinBoundaries(unsigned int iDim) const;            // Get all bins for given dimension

   std::vector < std::pair < double, double > > GetDim0BinBoundaries() const;                            // Get binning of 1st dimension
   std::vector < std::pair < double, double > > GetDim1BinBoundaries(unsigned int iDim0Bin) const;        // Get binning of 2nd dimension for bin iDim0Bin in 1st dimension
   std::vector < std::pair < double, double > > GetDim2BinBoundaries(unsigned int iDim0Bin, unsigned int iDim1Bin) const;        // Get binning of 2nd dimension for bin iDim0Bin in 1st dimension

   std::vector < double > GetLoBin(unsigned int iDim) const;                                    // Get lower bin boundaries in dim. iDim for all obs. bins
   double GetLoBin(int bin, unsigned int iDim) const {return Bin[bin][iDim].first;}             // Get lower bin boundary in dim. iDim for one obs. bin
   double GetLoBinMin(unsigned int iDim) const;                                                 // Get lowest bin boundary in dim. iDim of all obs. bins
   std::vector < double > GetUpBin(unsigned int iDim) const;                                    // Get upper bin boundaries in dim. iDim for all obs. bins
   double GetUpBin(int bin, unsigned int iDim) const {return Bin[bin][iDim].second;}            // Get upper bin boundary in dim. iDim for one obs. bin
   double GetUpBinMax(unsigned int iDim) const;                                                 // Get uppermost bin boundary in dim. iDim of all obs. bins

   vector < double > GetBinSize() const {return BinSize;};                                      // Get Binsize = BinSizeDim1 < * BinSizeDim2 >
   double GetBinSize(int bin) const {return BinSize[bin];};                                     // Get Binsize = BinSizeDim1 < * BinSizeDim2 >

   void SetNumDiffBin(int iDiff ) { NDim=iDiff; DimLabel.resize(NDim); IDiffBin.resize(NDim);}  // Set dimension of calculation. (Singledifferential, double-differntial, etc...)
   unsigned int GetNumDiffBin() const { return NDim; }                                                   // Get dimension of calculation. (Singledifferential, double-differntial, etc...)

   int GetIDiffBin(int bin) const { return IDiffBin[bin]; }                                     // Get if dimension is 'truly differential' or bin-integrated (divided by bin-width or not)

   void SetDimLabel( string label, unsigned int iDim , bool IsDiff = true );
   string GetDimLabel( int iDim  ) const {return DimLabel[iDim];};                              // Get label (name) of observable in dimension iDim
   vector<string > GetDimLabels() const {return DimLabel;};                                     // Get label (name) of all observables

   void SetIpublunits(int unit){Ipublunits = unit;}
   int GetIpublunits() const {return Ipublunits;}

   double GetEcms() const {return Ecms;}
   void SetEcms(double E) {Ecms = E;}

   int GetLoOrder() const {return ILOord;}
   void SetLoOrder(int LOOrd);

   bool IsNorm() const { return INormFlag == 0 ? false : true;}
   string GetDenomTable() const {return DenomTable;}

   fastNLOCoeffBase* GetCoeffTable(int no) const;
   fastNLOCoeffData* GetDataTable() const;                                                      //!< returns pointer to data table if available, else returns NULL pointer
   fastNLOCoeffAddBase* GetReferenceTable(ESMOrder eOrder) const;                               //!< returns pointer to reference table if available, else returns NULL pointer

   // useful functions

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
   void WriteScenario(ostream& table);
   void ReadScenario(istream& table);
   void ReadCoeffTables(istream& table);
   fastNLOCoeffBase* ReadRestOfCoeffTable(const fastNLOCoeffBase& cB, istream& table);

   vector < fastNLOCoeffBase* > fCoeff;
   //fastNLOCoeffData* fData;

   double Ecms;
   int ILOord;
   int Ipublunits;
   vector <string> ScDescript;

   // CKR: Changed to unsigned int
   unsigned int NObsBin;
   unsigned int NDim;

   vector <string> DimLabel;
   vector <int> IDiffBin;
   vector < vector <pair<double,double> > > Bin; // every bin has a lower and upper bin boundary and belongs to a 'dimension'. If a truely differential measurment, then upper bin boundary is equal lower one
   vector <double> BinSize;

   // contributions for normalization
   int INormFlag;
   string DenomTable;
   vector <int> IDivLoPointer;
   vector <int> IDivUpPointer;

};
#endif
