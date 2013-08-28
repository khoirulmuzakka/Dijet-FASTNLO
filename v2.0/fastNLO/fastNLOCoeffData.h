#ifndef __fastNLOCoeffData__
#define __fastNLOCoeffData__

#include "fastNLOCoeffBase.h"
#include "fastNLOConstants.h"

using namespace std;

class fastNLOCoeffData : public fastNLOCoeffBase {

   friend class fastNLOTable;

public:
   fastNLOCoeffData(int NObsBin);
   fastNLOCoeffData(const fastNLOCoeffBase&);
   int  Read(istream *table);
   void ReadRest(istream *table);
   virtual int Write(ostream *table, int option = 0);
   virtual int Copy(fastNLOCoeffData* other);
   virtual void Print() const;
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);

protected:
   fastNLOCoeffData();
   int ReadCoeffData(istream *table);

   int Nuncorrel;
   vector<string > UncDescr;
   int Ncorrel;
   vector<string > CorDescr;
   vector<double > Xcenter;
   vector<double > Value;
   v2d UncorLo;
   v2d UncorHi;
   v2d CorrLo;
   v2d CorrHi;
   int NErrMatrix;
   v2d matrixelement;

};

#endif
