#ifndef __fastNLOCoeffMult__
#define __fastNLOCoeffMult__

#include "fastNLOCoeffBase.h"
#include "fastNLOConstants.h"

using namespace std;

class fastNLOCoeffMult : public fastNLOCoeffBase {

   friend class fastNLOTable;
   friend class fastNLOCreate;

public:
   fastNLOCoeffMult();
   fastNLOCoeffMult(int NObsBin);
   fastNLOCoeffMult(const fastNLOCoeffBase&);
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   int Read(istream *table);
   void ReadRest(istream *table);
   virtual int Write(ostream *table, int option = 0);
   virtual int Copy(fastNLOCoeffMult* other);
   virtual void Print() const;
   
   double GetMultFactor(int iObs) const { return fact[iObs]; }
   vector<double > GetMultFactor() const { return fact; }
   vector<string> GetUncDescription() const { return UncDescr; }
   vector<string> GetCorDescription() const { return CorDescr; }
   v2d GetUncorLo() const { return UncorHi; };
   v2d GetUncorHi() const { return UncorLo; };
   v2d GetCorrLo()  const { return CorrLo; };
   v2d GetCorrHi()  const { return CorrHi; };

protected:
   int ReadCoeffMult(istream *table);

   int Nuncorrel;
   vector < string > UncDescr;
   int Ncorrel;
   vector < string > CorDescr;
   v2d UncorLo;
   v2d UncorHi;
   v2d CorrLo;
   v2d CorrHi;
   v1d fact;

};

#endif
