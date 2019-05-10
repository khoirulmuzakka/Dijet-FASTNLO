#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffUnc.h"
#include "fastnlotk/fastNLOTools.h"

using namespace std;
using namespace fastNLO;


//________________________________________________________________________________________________________________ //
fastNLOCoeffUnc::fastNLOCoeffUnc(int NObsBin)
   : fastNLOCoeffBase(NObsBin), NUncDescr(), UncDescr(), AsymUnc(), CorrUnc(),
     RelUncLo(), RelUncHi(), RelUncSym() {
   SetClassName("fastNLOCoeffUnc");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffUnc::fastNLOCoeffUnc(const fastNLOCoeffBase& base) : fastNLOCoeffBase(base) {
   SetClassName("fastNLOCoeffUnc");
}


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffUnc::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet)  {
   if ( c->GetIAddMultFlag()==0 && c->GetIDataFlag()==0 ) {
      // Additive contribution
      return false;
   } else if ( c->GetIAddMultFlag()==1 && c->GetIDataFlag()==0 ) {
      // Multiplicative contribution
      return false;
   } else if ( c->GetIAddMultFlag()==0 && c->GetIDataFlag()==1 ) {
      // Data contribution
      return false;
   } else if ( c->GetIAddMultFlag()==2 && c->GetIDataFlag()==0 ) {
      // Uncertainty contribution
      return true;
   } else {
      // Unknown contribution
      say::error["fastNLOCoeffUnc::CheckCoeffConstants"]
         << "Unknown contribution type, aborting! "
         << "IAddMultFlag = " << c->GetIAddMultFlag()
         << ", IDataFlag ="   << c->GetIDataFlag() <<endl;
      exit(1);
   }
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffUnc* fastNLOCoeffUnc::Clone() const {
   //! User has to take care to delete this object later
   return new fastNLOCoeffUnc(*this);
}


///________________________________________________________________________________________________________________ //
void fastNLOCoeffUnc::Read(istream& table){
   debug["Read"]<<"Start reading base content of uncertainty table contribution ..."<<endl;
   fastNLOCoeffBase::ReadBase(table);
   ReadRest(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffUnc::ReadRest(istream& table){
   debug["Read"]<<"Start reading rest of uncertainty table contribution ..."<<endl;
   CheckCoeffConstants(this);
   ReadCoeffUnc(table);
   EndReadCoeff(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffUnc::ReadCoeffUnc(istream& table){
   // char buffer[5257];
   // table >> NUncDescr;
   // cout << "nuncdescr = " << NUncDescr << endl;
   // for (int i=0; i<NUncDescr; i++) {
   //    table.getline(buffer,5256);
   //    UncDescr.push_back(buffer);
   //    cout << "descr = " << buffer << endl;
   // }
   table >> AsymUnc;
   cout << "asymunc = " << AsymUnc << endl;
   table >> CorrUnc;
   cout << "corrunc = " << CorrUnc << endl;
   // int itest;
   // table >> itest;
   // cout << "itest = " << itest << endl;
   // table >> itest;
   // cout << "itest = " << itest << endl;
   // table >> itest;
   // cout << "itest = " << itest << endl;
   // table >> itest;
   // cout << "itest = " << itest << endl;
   // table >> itest;
   // cout << "itest = " << itest << endl;
   RelUncLo.resize(fNObsBins);
   RelUncHi.resize(fNObsBins);
   RelUncSym.resize(fNObsBins);
   for(int i=0;i<fNObsBins;i++){
      if ( AsymUnc > 0 ) {
         table >> RelUncLo[i];
         table >> RelUncHi[i];
         RelUncSym[i] = 0.5*(RelUncHi[i]-RelUncLo[i]);
         cout << "AS: RelUncLo["<<i<<"] = "<<RelUncLo[i]<< endl;
         cout << "AS: RelUncHi["<<i<<"] = "<<RelUncHi[i]<< endl;
         cout << "AS: RelUncSym["<<i<<"] = "<<RelUncSym[i]<< endl;
      } else {
         table >> RelUncSym[i];
         RelUncLo[i] = -RelUncSym[i];
         RelUncHi[i] = +RelUncSym[i];
         cout << "SY: RelUncLo["<<i<<"] = "<<RelUncLo[i]<< endl;
         cout << "SY: RelUncHi["<<i<<"] = "<<RelUncHi[i]<< endl;
         cout << "SY: RelUncSym["<<i<<"] = "<<RelUncSym[i]<< endl;
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffUnc::Write(ostream& table, int itabversion) {
   fastNLOCoeffBase::Write(table,itabversion);
   CheckCoeffConstants(this);
   table << NUncDescr << sep;
   for (unsigned int i=0; i<UncDescr.size(); i++) {
      table << UncDescr[i] << sep;
   }
   table << AsymUnc << sep;
   table << CorrUnc << sep;
   for(int i=0;i<fNObsBins;i++){
      if ( AsymUnc > 0 ) {
         table << RelUncLo[i] << sep;
         table << RelUncHi[i] << sep;
      } else {
         table << RelUncSym[i] << sep;
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffUnc::Print(int iprint) const {
   if ( !(iprint < 0) ) {
      fastNLOCoeffBase::Print(iprint);
      cout << fastNLO::_DSEP20C << " fastNLO Table: CoeffUnc " << fastNLO::_DSEP20 << endl;
   } else {
      cout << endl << fastNLO::_CSEP20C << " fastNLO Table: CoeffUnc " << fastNLO::_CSEP20 << endl;
   }
   printf(" # No. of description lines (NUncDescr)     %d\n", NUncDescr);
   if (NUncDescr > 0) {fastNLOTools::PrintVector(UncDescr,"Uncertainty description (UncDescr)","#");}
   printf(" # Asymmetric uncertainty (AsymUnc)     %d\n", AsymUnc);
   printf(" # Correlated uncertainty (CorrUnc)     %d\n", CorrUnc);
   //
   // double minfact = *min_element(fact.begin(),fact.end());
   // double maxfact = *max_element(fact.begin(),fact.end());
   // printf(" # Minimal correction factor (fact[])  %f\n",minfact);
   // printf(" # Maximal correction factor (fact[])  %f\n",maxfact);
   // if ( std::abs(iprint) > 1 ) {
   //    cout << fastNLO::_SSEP20C << " Extended information (iprint > 1) " << fastNLO::_SSEP20 << endl;
   //    fastNLOTools::PrintVector(fact,"Correction factors (fact)","#    ");
   // }
   // if ( std::abs(iprint) > 2 ) {
   //    cout << fastNLO::_SSEP20C << " Extended information (iprint > 2) " << fastNLO::_SSEP20 << endl;
   //    for (int i=0; i<fNObsBins; i++) {
   //       // Print only for first and last observable bin
   //       if (i==0 || i==fNObsBins-1) {
   //          printf(" #       Observable bin no. %d\n",i+1);
   //          if (Nuncorrel > 0) {
   //             fastNLOTools::PrintVector(UncorLo[i],"Lower uncorr. uncertainties (UncorLo)","#      ");
   //             fastNLOTools::PrintVector(UncorHi[i],"Upper uncorr. uncertainties (UncorHi)","#      ");
   //          }
   //          if (Ncorrel > 0) {
   //             fastNLOTools::PrintVector(CorrLo[i],"Lower corr. uncertainties (CorrLo)","#      ");
   //             fastNLOTools::PrintVector(CorrHi[i],"Upper corr. uncertainties (CorrHi)","#      ");
   //          }
   //       }
   //    }
   // }
   cout << fastNLO::_CSEPSC << endl;
}


//________________________________________________________________________________________________________________ //
// Erase observable bin
void fastNLOCoeffUnc::EraseBin(unsigned int iObsIdx) {
   debug["fastNLOCoeffUnc::EraseBin"]<<"Erasing table entries in CoeffUnc for bin index " << iObsIdx << endl;
   if ( RelUncSym.size() == 0 ) {
      say::error["EraseBin"]<<"All uncertainty bins deleted already. Aborted!" << endl;
      exit(1);
   }
   if ( RelUncLo.size() != 0 )  RelUncLo.erase(RelUncLo.begin()+iObsIdx);
   if ( RelUncHi.size() != 0 )  RelUncHi.erase(RelUncHi.begin()+iObsIdx);
   if ( RelUncSym.size() != 0 ) RelUncSym.erase(RelUncSym.begin()+iObsIdx);
   fastNLOCoeffBase::EraseBin(iObsIdx);
}

// Catenate observable bin
void fastNLOCoeffUnc::CatBin(const fastNLOCoeffUnc& other, unsigned int iObsIdx) {
   debug["fastNLOCoeffUnc::CatBin"]<<"Catenating observable bin in CoeffUnc corresponding to bin index " << iObsIdx << endl;
   if ( RelUncSym.size() == 0 ) {
      say::error["CatBin"]<<"Initial uncertainty table is empty. Aborted!" << endl;
      exit(1);
   }
   unsigned int nold = RelUncSym.size();
   if ( RelUncLo.size() != 0 ) {
      RelUncLo.resize(nold+1);
      RelUncLo[nold] = other.RelUncLo[iObsIdx];
   }
   if ( RelUncHi.size() != 0 ) {
      RelUncHi.resize(nold+1);
      RelUncHi[nold] = other.RelUncHi[iObsIdx];
   }
   if ( RelUncSym.size() != 0 ) {
      RelUncSym.resize(nold+1);
      RelUncSym[nold] = other.RelUncSym[iObsIdx];
   }
   fastNLOCoeffBase::CatBin(other, iObsIdx);
}

// Multiply observable bin
void fastNLOCoeffUnc::MultiplyBin(unsigned int iObsIdx, double nfact) {
   debug["fastNLOCoeffUnc::MultiplyBin"]<<"Multiplying table entries. Nothing to be done in CoeffUnc." << endl;
}

//________________________________________________________________________________________________________________ //
bool fastNLOCoeffUnc::IsCatenable(const fastNLOCoeffUnc& other) const {
   //! Check for compatibility of catenating observable bins
   if ( ! ((fastNLOCoeffBase*)this)->IsCatenable(other)) return false;
   info["IsCatenable"]<<"Uncertainty contributions are catenable"<<endl;
   return true;
}
