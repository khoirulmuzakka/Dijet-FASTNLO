#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffMult.h"
#include "fastnlotk/fastNLOTools.h"

using namespace std;
using namespace fastNLO;


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffMult::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet)  {
   if ( c->GetIAddMultFlag()==0 && c->GetIDataFlag()==0 ) {
      // Additive contribution
      return false;
   } else if ( c->GetIAddMultFlag()==1 && c->GetIDataFlag()==0 ) {
      // Multiplicative contribution
      return true;
   } else if ( c->GetIAddMultFlag()==0 && c->GetIDataFlag()==1 ) {
      // Data contribution
      return false;
   } else {
      // Unknown contribution
      say::error["fastNLOCoeffMult::CheckCoeffConstants"]
         << "Unknown contribution type, aborting! "
         << "IAddMultFlag = " << c->GetIAddMultFlag()
         << ", IDataFlag ="   << c->GetIDataFlag() <<endl;
      exit(1);
   }
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
void fastNLOCoeffMult::Read(istream& table){
   fastNLOCoeffBase::ReadBase(table);
   ReadRest(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffMult::ReadRest(istream& table){
   CheckCoeffConstants(this);
   ReadCoeffMult(table);
   EndReadCoeff(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffMult::ReadCoeffMult(istream& table){
   char buffer[5257];
   table >> Nuncorrel;
   UncDescr.resize(Nuncorrel);
   table.getline(buffer,5256);
   for(int i=0;i<Nuncorrel;i++){
      table.getline(buffer,5256);
      UncDescr[i] = buffer;
      //         StripWhitespace(UncDescr[i]);
   }
   table >> Ncorrel;
   CorDescr.resize(Ncorrel);
   table.getline(buffer,5256);
   for(int i=0;i<Ncorrel;i++){
      table.getline(buffer,5256);
      CorDescr[i] = buffer;
      //         StripWhitespace(CorDescr[i]);
   }
   fact.resize(fNObsBins);
   UncorLo.resize(fNObsBins);
   UncorHi.resize(fNObsBins);
   CorrLo.resize(fNObsBins);
   CorrHi.resize(fNObsBins);
   for(int i=0;i<fNObsBins;i++){
      table >> fact[i];
      UncorLo[i].resize(Nuncorrel);
      UncorHi[i].resize(Nuncorrel);
      for(int j=0;j<Nuncorrel;j++){
         table >> UncorLo[i][j];
         table >> UncorHi[i][j];
      }
      CorrLo[i].resize(Ncorrel);
      CorrHi[i].resize(Ncorrel);
      for(int j=0;j<Ncorrel;j++){
         table >> CorrLo[i][j];
         table >> CorrHi[i][j];
      }
   }
   // end of IAddMultFlag==1
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffMult::Write(ostream& table) {
   fastNLOCoeffBase::Write(table);
   CheckCoeffConstants(this);
   table << Nuncorrel << endl;
   for(int i=0;i<Nuncorrel;i++){
      table << UncDescr[i]  << endl;
   }
   table << Ncorrel << endl;
   for(int i=0;i<Ncorrel;i++){
      table << CorDescr[i]  << endl;
   }
   for(int i=0;i<fNObsBins;i++){
      table << fact[i] << endl;
      for(int j=0;j<Nuncorrel;j++){
         table << UncorLo[i][j] << endl;
         table << UncorHi[i][j] << endl;
      }
      for(int j=0;j<Ncorrel;j++){
         table << CorrLo[i][j] << endl;
         table << CorrHi[i][j] << endl;
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffMult::Print(int iprint) const {
   if ( !(iprint < 0) ) {
      fastNLOCoeffBase::Print(iprint);
      cout << fastNLO::_DSEP20C << " fastNLO Table: CoeffMult " << fastNLO::_DSEP20 << endl;
   } else {
      cout << endl << fastNLO::_CSEP20C << " fastNLO Table: CoeffMult " << fastNLO::_CSEP20 << endl;
   }
   double minfact = *min_element(fact.begin(),fact.end());
   double maxfact = *max_element(fact.begin(),fact.end());
   printf(" # Minimal correction factor (fact[])  %f\n",minfact);
   printf(" # Maximal correction factor (fact[])  %f\n",maxfact);
   printf(" # No. of uncorr. unc. (Nuncorrel)     %d\n",Nuncorrel);
   if (Nuncorrel > 0) {fastNLOTools::PrintVector(UncDescr,"Uncorr. uncertainties (UncDescr)","#");}
   printf(" # No. of corr. unc. (Ncorrel)         %d\n",Ncorrel);
   if (Ncorrel > 0) {fastNLOTools::PrintVector(CorDescr,"Corr. uncertainties (CorDescr)","#");}
   if ( abs(iprint) > 1 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 1) " << fastNLO::_SSEP20 << endl;
      fastNLOTools::PrintVector(fact,"Correction factors (fact)","#    ");
   }
   if ( abs(iprint) > 2 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 2) " << fastNLO::_SSEP20 << endl;
      for (int i=0; i<fNObsBins; i++) {
         // Print only for first and last observable bin
         if (i==0 || i==fNObsBins-1) {
            printf(" #       Observable bin no. %d\n",i+1);
            if (Nuncorrel > 0) {
               fastNLOTools::PrintVector(UncorLo[i],"Lower uncorr. uncertainties (UncorLo)","#      ");
               fastNLOTools::PrintVector(UncorHi[i],"Upper uncorr. uncertainties (UncorHi)","#      ");
            }
            if (Ncorrel > 0) {
               fastNLOTools::PrintVector(CorrLo[i],"Lower corr. uncertainties (CorrLo)","#      ");
               fastNLOTools::PrintVector(CorrHi[i],"Upper corr. uncertainties (CorrHi)","#      ");
            }
         }
      }
   }
   cout << fastNLO::_CSEPSC << endl;
}


//________________________________________________________________________________________________________________ //
// void fastNLOCoeffMult::Cat(const fastNLOCoeffMult& other){
//    //! Concatenate bins of another coefficient table to this table
//    const fastNLOCoeffMult& cother = (const fastNLOCoeffMult&)other;
//    for ( int iObs=0; iObs<cother.GetNObsBin(); iObs++ ) {
//       CatBin(cother,iObs);
//    }
// }



// Erase observable bin
void fastNLOCoeffMult::EraseBin(unsigned int iObsIdx) {
   info["fastNLOCoeffMult::EraseBin"]<<"Erasing table entries in CoeffMult for bin index " << iObsIdx << endl;
   fact.erase(fact.begin()+iObsIdx);
   UncorLo.erase(UncorLo.begin()+iObsIdx);
   UncorHi.erase(UncorHi.begin()+iObsIdx);
   CorrLo.erase(CorrLo.begin()+iObsIdx);
   CorrHi.erase(CorrHi.begin()+iObsIdx);
   fastNLOCoeffBase::EraseBin(iObsIdx);
}

// Catenate observable bin
void fastNLOCoeffMult::CatBin(const fastNLOCoeffMult& other, unsigned int iObsIdx) {
   debug["fastNLOCoeffMult::CatBin"]<<"Catenating observable bin in CoeffMult corresponding to bin index " << iObsIdx << endl;
   if ( fact.size() == 0 ) {
      say::error["CatBin"]<<"fact size cannot be zero for a multiplicative table. Aborted!" << endl;
      exit(1);
   }
   unsigned int nold = fact.size();
   fact.resize(nold+1);
   fact[nold] = other.fact[iObsIdx];
   UncorLo.resize(nold+1);
   UncorLo[nold] = other.UncorLo[iObsIdx];
   UncorHi.resize(nold+1);
   UncorHi[nold] = other.UncorHi[iObsIdx];
   CorrLo.resize(nold+1);
   CorrLo[nold] = other.CorrLo[iObsIdx];
   CorrHi.resize(nold+1);
   CorrHi[nold] = other.CorrHi[iObsIdx];
   fastNLOCoeffBase::CatBin(other, iObsIdx);
}

//________________________________________________________________________________________________________________ //
bool fastNLOCoeffMult::IsCatenableContribution(const fastNLOCoeffMult& other) const {
   //! Check for compatibility of catenating observable bins
   if ( ! ((fastNLOCoeffBase*)this)->IsCatenableContribution(other)) return false;
   if( Nuncorrel != other.GetNuncorrel() ){
      debug["IsCatenableContribution"]<<"Nuncorrel != other.GetNuncorrel()"<<endl;
      return false;
   }
   if( Ncorrel != other.GetNcorrel() ){
      debug["IsCatenableContribution"]<<"Ncorrel != other.GetNcorrel()"<<endl;
      return false;
   }
   info["IsCatenableContribution"]<<"Both multiplicable contributions are catenable"<<endl;
   return true;
}
