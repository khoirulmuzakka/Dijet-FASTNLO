#ifndef __fastNLOCoeffUnc__
#define __fastNLOCoeffUnc__

#include "fastNLOCoeffBase.h"
#include "fastNLOConstants.h"


class fastNLOCoeffUnc : public fastNLOCoeffBase {

   friend class fastNLOTable;
   friend class fastNLOCreate;

public:
   fastNLOCoeffUnc() = delete;
   fastNLOCoeffUnc(int NObsBin);
   explicit fastNLOCoeffUnc(const fastNLOCoeffBase&);
   virtual ~fastNLOCoeffUnc(){;};
   virtual fastNLOCoeffUnc* Clone() const;   //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   virtual void Read(std::istream& table);
   virtual void Write(std::ostream& table, int ITabVersionWrite);
   virtual void Print(int iprint) const;

   double GetRelUncLo(int iObs) const  { return RelUncLo[iObs]; }
   double GetRelUncHi(int iObs) const  { return RelUncHi[iObs]; }
   double GetRelUncSym(int iObs) const { return RelUncSym[iObs]; }
   int GetNUncDescr() const {return NUncDescr;}
   void SetNUncDescr(int n){NUncDescr = n;}
   std::vector<std::string> GetUncDescription() const { return UncDescr; }
   std::vector<double > GetRelUncLo()  const { return RelUncLo; };
   std::vector<double > GetRelUncHi()  const { return RelUncHi; };
   std::vector<double > GetRelUncSym() const { return RelUncSym; };

   // Erase observable bin; iObsIdx is the C++ array index to be removed and
   // not the observable bin no. running from 1 to fNObsBins
   virtual void EraseBin(unsigned int iObsIdx);
   virtual void MultiplyBin(unsigned int iObsIdx, double fact);
   bool IsCatenable(const fastNLOCoeffUnc& other) const;
   // Catenate observable to table
   virtual void CatBin(const fastNLOCoeffUnc& other, unsigned int iObsIdx);



protected:
   void ReadCoeffUnc(std::istream& table);
   void ReadRest(std::istream& table);

   int NUncDescr;
   std::vector < std::string > UncDescr;
   unsigned int AsymUnc;
   unsigned int CorrUnc;
   fastNLO::v1d RelUncLo;
   fastNLO::v1d RelUncHi;
   fastNLO::v1d RelUncSym;

};

#endif
