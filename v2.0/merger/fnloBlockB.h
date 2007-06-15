#ifndef __fnloBlockB__
#define __fnloBlockB__

#include <string>
#include <fstream>
#include <math.h>
#include <vector>

#include "fnloconstants.h"

using namespace std;

class fnloBlockB {
 public:
   int Read(istream *table);
   int Write(ostream *table);
   bool IsCompatible(fnloBlockB* other);
 protected:
   int IXsectUnits;
   int IdataFlag;
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
   long int Nevt;
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
   int Nxtot;
   vector < vector < double > > XNode1;
   int Nxtot2;
   vector < vector < double > > XNode2;
   int Nztot;
   vector < vector < double > > ZNode;
   int NScales;
   int NScaleDim;
   vector < int > Iscale;
   vector < string > NscaleDescript;
   vector < vector < string > > ScaleDescript;
   vector < int > Nscalevar;
   vector < int > Nscalenode;
   vector < vector < double > > ScaleFac;
   vector < vector < vector < vector < double > > > > ScaleNode;
   vector < vector < vector < vector < vector < vector < double > > > > > > SigmaTilde;

};
#endif
