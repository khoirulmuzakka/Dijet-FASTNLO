#ifndef __fnloBlockB__
#define __fnloBlockB__

#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <iostream>

#include "fnloconstants.h"
#include "fnloBlockA1.h"
#include "fnloBlockA2.h"

using namespace std;

class fnloTableUser;

class fnloBlockB {
public:
   fnloBlockB(){;}
   fnloBlockB(fnloBlockA1 *blocka1, fnloBlockA2 *blocka2) : BlockA1(blocka1) ,  BlockA2(blocka2)   {;}
   int Read(istream *table);
   int Write(ostream *table, int option = 0);
   int Copy(fnloBlockB* other);

   void ResizeTable( vector<double >* v, int dim0 );
   void ResizeTable( vector<vector<double > >*  v, int dim0 , int dim1 );
   void ResizeTable( vector<vector<double > >*  v, int dim0 , int* dim1GetNxmaxFromDimI );
   void ResizeTable( vector<vector<vector<double > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2 );
   void ResizeTable( vector<vector<vector<double > > >* v, int dim0 , int dim1, int dim2 );
   void ResizeTable( vector<vector<vector<vector<double > > > >* v, int dim0 , int dim1, int dim2, int dim3 );
   void ResizeTable( vector<vector<vector<vector<double > > > >* v, int dim0 , int dim1, int* dim2GetNxmaxFromDimI, int dim3 );
   void ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4 );
   void ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int dim1, int dim2, int* dim3GetNxmaxFromDimI, int dim4 );
   void ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2, int dim3, int dim4 );
   void ResizeTable( vector<vector<vector<vector<vector<vector<double > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 );
   void ResizeTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 );
   void ResizeTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int* dim5GetNxmaxFromDimI , int dim6 );

   int ReadTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, istream *table );
   int ReadTable( vector<vector<vector<vector<vector<vector<double > > > > > >* v, istream *table );
   int ReadTable( vector<vector<vector<vector<vector<double > > > > >* v, istream *table );
   int ReadTable( vector<vector<vector<vector<double > > > >* v, istream *table );
   int ReadTable( vector<vector<vector<double > > >* v, istream *table );
   int ReadTable( vector<vector<double > >* v, istream *table );
   int ReadTable( vector<double >* v, istream *table );

   int ReadFlexibleVector( vector<double >* v, istream *table , bool nProcLast = false );
   int ReadFlexibleVector( vector<vector<double > >* v, istream *table , bool nProcLast = false );
   int ReadFlexibleVector( vector<vector<vector<double > > >* v, istream *table , bool nProcLast = false );
   int ReadFlexibleVector( vector<vector<vector<vector<double > > > >* v, istream *table , bool nProcLast = false );
   int ReadFlexibleVector( vector<vector<vector<vector<vector<double > > > > >* v, istream *table , bool nProcLast = false );
   int ReadFlexibleVector( vector<vector<vector<vector<vector<vector<double > > > > > >* v, istream *table , bool nProcLast = false );
   int ReadFlexibleVector( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, istream *table , bool nProcLast = false );

   int WriteTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, ostream *table , bool DivByNevt=false , int Nevt=1 );
   int WriteTable( vector<vector<vector<vector<vector<vector<double > > > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( vector<vector<vector<vector<vector<double > > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( vector<vector<vector<vector<double > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( vector<vector<vector<double > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( vector<vector<double > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( vector<double >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );

   int WriteFlexibleTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, ostream *table , bool DivByNevt=false , int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<vector<vector<vector<vector<vector<double > > > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<vector<vector<vector<vector<double > > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<vector<vector<vector<double > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<vector<vector<double > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<vector<double > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<double >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );

   void AddTableToAnotherTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* vSum, vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* vAdd, double w1 = 0, double w2 = 0 );
   void AddTableToAnotherTable( vector<vector<vector<vector<vector<vector<double > > > > > >* vSum, vector<vector<vector<vector<vector<vector<double > > > > > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( vector<vector<vector<vector<vector<double > > > > >* vSum, vector<vector<vector<vector<vector<double > > > > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( vector<vector<vector<vector<double > > > >* vSum, vector<vector<vector<vector<double > > > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( vector<vector<vector<double > > >* vSum, vector<vector<vector<double > > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( vector<vector<double > >* vSum, vector<vector<double > >* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( vector<double >* vSum, vector<double >* vAdd, double w1 = 1, double w2 = 1 );

   void Print();

   bool IsCompatible(fnloBlockB* other);
   int GetIRef(){return IRef;}
   int GetIDataFlag(){return IDataFlag;}
   int GetIAddMultFlag(){return IAddMultFlag;}
   int GetIContrFlag1(){return IContrFlag1;}
   int GetIContrFlag2(){return IContrFlag2;}
   int GetIContrFlag3(){return IContrFlag3;}
   int GetNpow(){return Npow;}
   long long int GetNevt(){return Nevt;}
   int GetNxmax(int Obsbin);
   int GetXIndex(int Obsbin,int x1bin,int x2bin =0);
   void Add(fnloBlockB* other);
   bool IsLO(){return IContrFlag1==1 && IContrFlag2==1;}
   bool IsNLO(){return IContrFlag1==1 && IContrFlag2==2;}
   bool IsReference(){return IRef>0;};
   double GetAlphas(double mur, double asmz);
   int GetTotalScalevars();
   int GetTotalScalenodes();

   // from here on "user" code
   void ResetXsection();
   void FillPDFCache(int scalevar, void (fnloTableUser::*GetPdfs)(double x, double muf ,vector<double> &xfx),fnloTableUser *tableptr);
   void FillPDFCache(int scalevar, void (fnloTableUser::*GetPdfs)(double x, double muf,vector<double> &xfx),
                     void (fnloTableUser::*GetPdfs2)(double x, double muf,vector<double> &xfx), fnloTableUser *tableptr);
   void CalcXsection(double asmz,  int scalevar = 0, double rescale=1.);
   double GetXsection(int bin){return Xsection[bin];}
   double GetSmallestX(int ObsBin);
   double GetSmallestX2(int ObsBin);

 private:
   void StripWhitespace(string &str);
   void CalcPDFLinearComb(vector<double> pdfx1, vector<double> pdfx2, vector<double> *result);
   int GetScaledimfromvar(int scalevar);

 public:
   static const int DividebyNevt = 1;

   fnloBlockA1 *BlockA1;
   fnloBlockA2 *BlockA2;
   int IXsectUnits;
   int IDataFlag;
   int IAddMultFlag;
   int IContrFlag1;
   int IContrFlag2;
   int IContrFlag3;	// @MW, @KR: IContrFlag3 was replaced by NScaleDep by DB and is now without any usage
   int NScaleDep;
   // obsolete int NContrDescr;
   vector < string > CtrbDescript;
   // obsolete   int NCodeDescr;
   vector < string > CodeDescript;
   int Nuncorrel;
   vector < string > UncDescr;
   int Ncorrel;
   vector < string > CorDescr;
   vector < double > Xcenter;
   vector < double > Value;
   vector < vector < double > > UncorLo;
   vector < vector < double > > UncorHi;
   vector < vector < double > > CorrLo;
   vector < vector < double > > CorrHi;
   int NErrMatrix;
   vector < vector < double > > matrixelement;
   vector < double > fact;
   int IRef;
   int IScaleDep;
   unsigned long long int Nevt;
   int Npow;
   int NPDF;
   vector < int > NPDFPDG;
   int NPDFDim;
   int NFragFunc;
   vector < int > NFFPDG;
   int NFFDim;
   int NSubproc;
   int IPDFdef1;
   int IPDFdef2;
   int IPDFdef3;
   // Missing: linear PDF combinations for IPDFdef1=0
   vector < int > Nxtot1;
   vector < double > Hxlim1;
   vector < vector < double > > XNode1;
   vector < int >  Nxtot2;
   vector < double > Hxlim2;
   vector < vector < double > > XNode2;
   vector < int > Nztot;
   vector < double > Hzlim;
   vector < vector < double > > ZNode;
   int NScales;
   int NScaleDim;
   vector < int > Iscale;
   // obsolete vector < int > NscaleDescript;
   vector < vector < string > > ScaleDescript;
   vector < int > Nscalevar;
   vector < int > Nscalenode;
   vector < vector < double > > ScaleFac;
   vector < vector < vector < vector < double > > > > ScaleNode;
   vector < vector < vector < vector < double > > > > HScaleNode;

   // DB: todo: those variables should end up in a new class fnloBlockBXS : public fnloBlockB
   vector < vector < vector < vector < vector < double > > > > > SigmaTilde; 
   vector < vector < vector < vector < double > > > > PdfLc; 
   vector < double > Xsection; 

   // --------------------------- mu_f, mu_r variaton scheme --------------------------- //
   // ---- members to write to disc ---- //
   // SigmaTilde [NObsBins] ['n' x-nodes] [n s1-Nodes] [n s2-Nodes] [nsubproc]
   vector < vector < vector < vector < vector < double > > > > > SigmaTildeMuIndep; 
   vector < vector < vector < vector < vector < double > > > > > SigmaTildeMuFDep; 
   vector < vector < vector < vector < vector < double > > > > > SigmaTildeMuRDep; 
   // SigmaRef [NObsBins] [nsubproc]
   vector < vector < double > > SigmaRefMixed; 
   vector < vector < double > > SigmaRef_s1; 
   vector < vector < double > > SigmaRef_s2; 
  
   int NscalenodeScale1;
   int NscalenodeScale2;
   // ScaleNodeXY [ObsBin] [NscalenodeScaleX]  
   vector < vector < double > > ScaleNode1;
   vector < vector < double > > ScaleNode2;
   vector < vector < double > > HScaleNode1;
   vector < vector < double > > HScaleNode2;
  
   // ---- stuff for reading the table ---- //
   //    // @MW, @KR: if we provide a separate reader: we don't need variables here for reading the table.
   //    vector < double > XsectionRef_s1; 
   //    vector < double > XsectionRef_s2; 
   //    vector < double > XsectionRefMixed; 
   //    vector < double > XsectionRefMuf1_MuRMixed; 
   //    vector < double > XsectionMuVar; 
   //    // PdfLcMuVar [ObsBins] [NxNodes] [NQNodes] [NPtNodes] [nsubproc]
   //    vector < vector < vector < vector < vector <double > > > > > PdfLcMuVar; 
   //    // AlphasTwoPi [ObsBins] [N s1-Nodes] [N s2-Nodes]
   //    vector < vector < vector <double > > > AlphasTwoPi; 
  
};
#endif
