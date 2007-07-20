#ifndef __fnloBlockB__
#define __fnloBlockB__

#include <string>
#include <fstream>
#include <math.h>
#include <vector>

#include "fnloconstants.h"
#include "fnloBlockA1.h"
#include "fnloBlockA2.h"

using namespace std;

class fnloBlockB {
 public:
   fnloBlockB(fnloBlockA1 *blocka1, fnloBlockA2 *blocka2) : BlockA1(blocka1) ,  BlockA2(blocka2)   {;}
   int Read(istream *table);
   int Write(ostream *table);
   bool IsCompatible(fnloBlockB* other);
   int GetIDataFlag(){return IDataFlag;}
   int GetIAddMultFlag(){return IAddMultFlag;}
   int GetIContrFlag1(){return IContrFlag1;}
   int GetIContrFlag2(){return IContrFlag2;}
   int GetIContrFlag3(){return IContrFlag3;}
   int GetNpow(){return Npow;}
   long long int GetNevt(){return Nevt;}
   void Add(fnloBlockB* other);

 private:
   void StripWhitespace(string &str);
 protected:
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
   vector < int > Nxtot;
   vector < vector < double > > XNode1;
   vector < int >  Nxtot2;
   vector < vector < double > > XNode2;
   vector < int > Nztot;
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
   vector < vector < vector < vector < vector < vector < double > > > > > > SigmaTilde;

};
#endif
