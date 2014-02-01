#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffMult.h"

using namespace std;
using namespace fastNLO;


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffMult::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet)  {
  if ( c->GetIDataFlag()==0 && c->GetIAddMultFlag()==1 ) return true;
   else if ( c->GetIAddMultFlag()!=1 ) {
      if ( !quiet)
	 say::error["fastNLOCoeffMult::CheckCoeffConstants"]
	    <<"This is not an additive contribution (IAddMultFlag="<<c->GetIAddMultFlag()<<", but must be equals 1)."<<endl;
      return false;
   }
   else if ( c->GetIDataFlag()!=0) {
      say::error["fastNLOCoeffMult::CheckCoeffConstants"]
	 <<"This seems to be an additive contribution, but is also labeled as data table. (IAddMultFlag="<<c->GetIAddMultFlag()<<", IDataFlag="<<c->GetIDataFlag()<<endl;
      exit(1);
      return false;
   }
   else return true;
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffMult::fastNLOCoeffMult(){
   SetClassName("fastNLOCoeffMult");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffMult::fastNLOCoeffMult(int NObsBin) : fastNLOCoeffBase(NObsBin){
   SetClassName("fastNLOCoeffMult");
   fNObsBins = NObsBin;
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffMult::fastNLOCoeffMult(const fastNLOCoeffBase& base) : fastNLOCoeffBase(base) {
   SetClassName("fastNLOCoeffMult");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffBase* fastNLOCoeffMult::Clone() const {
   //! Use has to take care to delete this object later
   return static_cast<fastNLOCoeffBase*>(new fastNLOCoeffMult(*this));
}


///________________________________________________________________________________________________________________ //
int fastNLOCoeffMult::Read(istream *table){
   fastNLOCoeffBase::ReadBase(table);
   ReadRest(table);
   return 0;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffMult::ReadRest(istream *table){
   CheckCoeffConstants(this);
   ReadCoeffMult(table);
   EndReadCoeff(table);
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffMult::ReadCoeffMult(istream *table){
   char buffer[5257];
   *table >> Nuncorrel;
   UncDescr.resize(Nuncorrel);
   table->getline(buffer,5256);
   for(int i=0;i<Nuncorrel;i++){
      table->getline(buffer,5256);
      UncDescr[i] = buffer;
      //         StripWhitespace(UncDescr[i]);
   }
   *table >> Ncorrel;
   CorDescr.resize(Ncorrel);
   table->getline(buffer,5256);
   for(int i=0;i<Ncorrel;i++){
      table->getline(buffer,5256);
      CorDescr[i] = buffer;
      //         StripWhitespace(CorDescr[i]);
   }
   fact.resize(fNObsBins);
   UncorLo.resize(fNObsBins);
   UncorHi.resize(fNObsBins);
   CorrLo.resize(fNObsBins);
   CorrHi.resize(fNObsBins);
   for(int i=0;i<fNObsBins;i++){
      *table >> fact[i];
      UncorLo[i].resize(Nuncorrel);
      UncorHi[i].resize(Nuncorrel);
      for(int j=0;j<Nuncorrel;j++){
	 *table >> UncorLo[i][j];
	 *table >> UncorHi[i][j];
      }
      CorrLo[i].resize(Ncorrel);
      CorrHi[i].resize(Ncorrel);
      for(int j=0;j<Ncorrel;j++){
	 *table >> CorrLo[i][j];
	 *table >> CorrHi[i][j];
      }
   }
   // end of IAddMultFlag==1
   return 0;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffMult::Write(ostream *table) {
   fastNLOCoeffBase::Write(table);
   CheckCoeffConstants(this);
   *table << Nuncorrel << endl;
   for(int i=0;i<Nuncorrel;i++){
      *table << UncDescr[i]  << endl;
   }
   *table << Ncorrel << endl;
   for(int i=0;i<Ncorrel;i++){
      *table << CorDescr[i]  << endl;
   }
   for(int i=0;i<fNObsBins;i++){
      *table << fact[i] << endl;
      for(int j=0;j<Nuncorrel;j++){
	 *table << UncorLo[i][j] << endl;
	 *table << UncorHi[i][j] << endl;
      }
      for(int j=0;j<Ncorrel;j++){
	 *table << CorrLo[i][j] << endl;
	 *table << CorrHi[i][j] << endl;
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffMult::Print() const {
  printf(" **************** FastNLO Table: CoeffMult ****************\n");
  fastNLOCoeffBase::Print();
  if(IAddMultFlag==1){
    printf(" B   fastNLOCoeffMult::Print(). Printing of multiplicative factors not yet implemented (IAddMultFlag==1).\n");
  }
  printf(" *******************************************************\n");
}


//________________________________________________________________________________________________________________ //
