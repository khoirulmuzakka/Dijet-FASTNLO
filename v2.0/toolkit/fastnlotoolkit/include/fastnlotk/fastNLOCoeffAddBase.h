#ifndef __fastNLOCoeffAddBase__
#define __fastNLOCoeffAddBase__


#include "fastNLOCoeffBase.h"
#include "fastNLOConstants.h"

using namespace std;

class fastNLOCoeffAddBase : public fastNLOCoeffBase {

   friend class fastNLOTable;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddBase(int NObsBin);
   fastNLOCoeffAddBase(const fastNLOCoeffBase& base);
   virtual ~fastNLOCoeffAddBase(){;}
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   int Read(istream *table);
   virtual void Write(ostream *table);
   virtual void Add(const fastNLOCoeffAddBase& other);
   virtual void Print() const;

   void ResizeTable( v2d& v, int dim0 , int* dim1GetNxmaxFromDimI );
   void ResizeTable( v3d& v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2 );
   void ResizeTable( v4d& v, int dim0 , int dim1, int* dim2GetNxmaxFromDimI, int dim3 );
   void ResizeTable( v5d& v, int dim0 , int dim1, int dim2, int* dim3GetNxmaxFromDimI, int dim4 );
   void ResizeTable( v5d& v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2, int dim3, int dim4 );
   void ResizeTable( v7d& v, int dim0 , int dim1, int dim2, int dim3, int dim4, int* dim5GetNxmaxFromDimI , int dim6 );


   int GetIRef() const {return IRef;}
   long long int GetNevt(int NObsBin, int NSubproc) const {
      if ( Nevt > 0 ) return Nevt;
      else {cout<<"Todo. Preparation for v2.2."<<endl; return Nevt;}
   }
   int GetNxmax(int Obsbin) const ;
   int GetXIndex(int Obsbin,int x1bin,int x2bin =0) const ;
   int GetNSubproc() const { return NSubproc;}
   int GetIScaleDep() const { return IScaleDep;}
   int GetNPDF() const {return NPDFPDG.size();}
   int GetPDFPDG(int iPDF) const {return NPDFPDG[iPDF];}
   int GetNPDFDim() const {return NPDFDim;}
   int GetIPDFdef1() const { return IPDFdef1; }
   int GetIPDFdef2() const { return IPDFdef2; }
   int GetIPDFdef3() const { return IPDFdef3; }
   int GetNpow() const {return Npow;}
   //vector<string > GetScaleDescript(int iScale=0) const { return ScaleDescript[iScale]; };
   string GetScaleDescription(int iScale=0) const { return ScaleDescript[0][iScale]; };		// getter for scale description of scale iScale
   vector<vector<string > > GetScaleDescr() const { return ScaleDescript; }
   int GetNxtot1(int iBin ) const { return XNode1[iBin].size(); }
   int GetNxtot2(int iBin ) const { return XNode2[iBin].size(); }

   double GetXNode1(int iObsBin, int iNode) const { return XNode1[iObsBin][iNode]; } 
   double GetXNode2(int iObsBin, int iNode) const { return XNode2[iObsBin][iNode]; } 

   bool IsReference() const {return IRef>0;};
   bool IsCompatible(const fastNLOCoeffAddBase& other) const;

protected:
   fastNLOCoeffAddBase();
   int ReadCoeffAddBase(istream *table);
   int GetScaledimfromvar(int scalevar) const;

   int IRef;
   int IScaleDep;
   unsigned long long int Nevt;
   int Npow;
   vector < int > NPDFPDG;
   int NPDFDim;
   vector < int > NFFPDG;
   int NFFDim;
   int NSubproc;
   int IPDFdef1;
   int IPDFdef2;
   int IPDFdef3;
   // Missing: linear PDF combinations for IPDFdef1=0
   vector < double > Hxlim1;
   v2d XNode1;
   vector < double > Hxlim2;
   v2d XNode2;
   vector < int > Nztot;
   vector < double > Hzlim;
   v2d ZNode;
   int NScaleDim;
   vector < int > Iscale;									// not used
   vector < vector < string > > ScaleDescript;

};

#endif
