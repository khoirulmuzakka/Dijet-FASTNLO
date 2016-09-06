#ifndef __fastNLOCoeffData__
#define __fastNLOCoeffData__

#include "fastNLOCoeffBase.h"
#include "fastNLOConstants.h"


class fastNLOCoeffData : public fastNLOCoeffBase {

   friend class fastNLOTable;

public:
   fastNLOCoeffData(int NObsBin);
   fastNLOCoeffData(const fastNLOCoeffBase&);
   virtual ~fastNLOCoeffData(){;};
   virtual fastNLOCoeffBase* Clone() const;                                     //!< returns 'new' copy of this instance.
   virtual void Read(std::istream& table);
   virtual void Write(std::ostream& table);
   virtual void Print(int iprint) const;
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);

   // Erase observable bin; iObsIdx is the C++ array index to be removed and
   // not the observable bin no. running from 1 to fNObsBins
   virtual void EraseBin(unsigned int iObsIdx);

protected:
   fastNLOCoeffData();
   void ReadCoeffData(std::istream& table);
   void ReadRest(std::istream& table);

   int Nuncorrel;
   std::vector<std::string > UncDescr;
   int Ncorrel;
   std::vector<std::string > CorDescr;
   std::vector<double > Xcenter;
   std::vector<double > Value;
   fastNLO::v2d UncorLo;
   fastNLO::v2d UncorHi;
   fastNLO::v2d CorrLo;
   fastNLO::v2d CorrHi;
   int NErrMatrix;
   fastNLO::v2d matrixelement;

};

#endif
