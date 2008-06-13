#ifndef __fnloBlockB__
#define __fnloBlockB__

#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>


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
   void Add(fnloBlockB* other);
   bool IsLO(){return IContrFlag1==1 && IContrFlag2==1;}
   bool IsNLO(){return IContrFlag1==1 && IContrFlag2==2;}
   bool IsReference(){return IRef>0;};
   double GetAlphas(double mur, double asmz);

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

 public:
   static const int DividebyNevt = 1;

   fnloBlockA1 *BlockA1;
   fnloBlockA2 *BlockA2;
   int IXsectUnits;
   int IDataFlag;
   int IAddMultFlag;
   int IContrFlag1;
   int IContrFlag2;
   int IContrFlag3;
   int NContrDescr;
   vector < string > CtrbDescript;
   int NCodeDescr;
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
   vector < int > NscaleDescript;
   vector < vector < string > > ScaleDescript;
   vector < int > Nscalevar;
   vector < int > Nscalenode;
   vector < vector < double > > ScaleFac;
   vector < vector < vector < vector < double > > > > ScaleNode;

   // the follwoing is valid only for one scale dimension
   vector < vector < vector < vector < vector < double > > > > > SigmaTilde; 
   vector < vector < vector < vector < double > > > > PdfLc; 
   vector < double > Xsection; 

   static const double TWOPI = 6.28318530717958647692528;
   static const double TWOPISQR = 39.47841760435743447533796;
   static const int NF = 5;
   static const double MZ = 91.1882;

};
#endif
