#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffBase.h"
#include "fastnlotk/fastNLOConstants.h"
#include "fastnlotk/fastNLOTools.h"

using namespace std;

//________________________________________________________________________________________________________________ //
fastNLOCoeffBase::fastNLOCoeffBase() : PrimalScream("fastNLOCoeffBase"){
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffBase::fastNLOCoeffBase(int NObsBin) : PrimalScream("fastNLOCoeffBase") {
   fNObsBins = NObsBin;
}

//________________________________________________________________________________________________________________ //
fastNLOCoeffBase* fastNLOCoeffBase::Clone() const {
   //! User has to take care to delete this object later
   return new fastNLOCoeffBase(*this);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::Read(istream& table){
   // basic read function.
   // reads in only 'base'-variables.
   debug["Read"]<<endl;
   ReadBase(table);
   EndReadCoeff(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::ReadBase(istream& table){
   debug["ReadBase"]<<endl;
   table.peek();
   if (table.eof()){
      //printf("fastNLOCoeffBase::Read: Cannot read from file.\n");
      error["ReadBase"]<<"Cannot read from file."<<endl;
   }

   if (!fastNLOTools::ReadMagicNo(table)) {
      say::error["ReadBase"]<<"Did not find initial magic number, aborting!"<<endl;
      say::error["ReadBase"]<<"Please check compatibility of tables and program version!"<<endl;
      exit(1);
   }

   table >> IXsectUnits;
   table >> IDataFlag;
   table >> IAddMultFlag;
   table >> IContrFlag1;
   table >> IContrFlag2;
   table >> NScaleDep;
   int NContrDescr;
   table >> NContrDescr;
   //   printf("  *  fastNLOCoeffBase::Read().  IDataFlag: %d, IAddMultFlag: %d, IContrFlag1: %d, IContrFlag2: %d,, NScaleDep: %d\n",IDataFlag,IAddMultFlag,IContrFlag1,IContrFlag2,NScaleDep );
   CtrbDescript.resize(NContrDescr);
   char buffer[5257];
   table.getline(buffer,5256);
   for(int i=0;i<NContrDescr;i++){
      table.getline(buffer,256);
      CtrbDescript[i] = buffer;
      //      StripWhitespace(CtrbDescript[i]);
   }
   int NCodeDescr;
   table >> NCodeDescr;
   CodeDescript.resize(NCodeDescr);
   table.getline(buffer,256);
   for(int i=0;i<NCodeDescr;i++){
      table.getline(buffer,256);
      CodeDescript[i] = buffer;
      //      StripWhitespace(CodeDescript[i]);
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::EndReadCoeff(istream& table){
   if (!fastNLOTools::ReadMagicNo(table)) {
      say::error["ReadBase"]<<"Did not find final magic number, aborting!"<<endl;
      say::error["ReadBase"]<<"Please check compatibility of tables and program version!"<<endl;
      say::error["ReadBase"]<<"This might also be provoked by lines with unexpected non-numeric content like 'inf' or 'nan'!"<<endl;
      exit(1);
   }
   fastNLOTools::PutBackMagicNo(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::Write(ostream& table) {
   say::debug["Write"]<<"Writing fastNLOCoeffBase."<<endl;
   table << fastNLO::tablemagicno << endl;
   table << IXsectUnits << endl;
   table << IDataFlag << endl;
   table << IAddMultFlag << endl;
   table << IContrFlag1 << endl;
   table << IContrFlag2 << endl;
   table << NScaleDep << endl;
   table << CtrbDescript.size() << endl;
   //printf("  *  fastNLOCoeffBase::Write().  IDataFlag: %d, IAddMultFlag: %d, IContrFlag1: %d, IContrFlag2: %d, NScaleDep: %d\n",
   //IDataFlag,IAddMultFlag,IContrFlag1,IContrFlag2,NScaleDep);
   for(unsigned int i=0;i<CtrbDescript.size();i++){
      table << CtrbDescript[i] << endl;
   }
   table << CodeDescript.size() << endl;
   for(unsigned int i=0;i<CodeDescript.size();i++){
      table << CodeDescript[i] << endl;
   }
}


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffBase::IsCompatible(const fastNLOCoeffBase& other) const {
   if( fNObsBins != other.GetNObsBin() ){
      warn["IsCompatible"]<<"fNObsBins != other.GetNObsBin()"<<endl;
      return false;
   }
   if( IXsectUnits != other.GetIXsectUnits() ){
      warn["IsCompatible"]<<"IXsectUnits != other.GetIXsectUnits()"<<endl;
      return false;
   }
   if( IDataFlag != other.GetIDataFlag() ){
      debug["IsCompatible"]<<"IDataFlag != other.GetIDataFlag()"<<endl;
      return false;
   }
   if( IAddMultFlag != other.GetIAddMultFlag() ){
      debug["IsCompatible"]<<"IAddMultFlag != other.GetIAddMultFlag()"<<endl;
      return false;
   }
   if( IContrFlag1 != other.GetIContrFlag1() ){
      debug["IsCompatible"]<<"IContrFlag1 != other.GetIContrFlag1()"<<endl;
      return false;
   }
   if( IContrFlag2 != other.GetIContrFlag2() ){
      debug["IsCompatible"]<<"IContrFlag2 != other.GetIContrFlag2()"<<endl;
      return false;
   }
   if( NScaleDep != other.GetNScaleDep() ){
      debug["IsCompatible"]<<"NScaleDep != other.GetNScaleDep()"<<endl;
      if ( (NScaleDep==5 && other.GetNScaleDep()==6) || (NScaleDep==6 && other.GetNScaleDep()==5) ) {
         debug["IsCompatible"]<<"One table with NScale=5 and one with NScaleDep=6"<<endl;
         // continue;
      }
      else {
         warn["IsCompatible"]<<"Incompatible NScaleDep found!()"<<endl;
         return false;
      }
   }
   debug["IsCompatible"]<<"Both tables are compatible"<<endl;
   // check descripts here ?!
   //bool potentialcompatible = true;
   //vector < string > CtrbDescript;
   //vector < string > CodeDescript;
   return true;
}


//________________________________________________________________________________________________________________ //


void fastNLOCoeffBase::SetCoeffAddDefaults(){
  SetIDataFlag(0);
  SetIAddMultFlag(0);
  SetIContrFlag1(1);
  SetIContrFlag2(100); // specify if LO or NLO
  SetNScaleDep(0);
  SetIXsectUnits(12);
};

//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::Print(int iprint) const {
   if ( !(iprint < 0) ) {
      cout << fastNLO::_DSEP20C << " fastNLO Table: CoeffBase " << fastNLO::_DSEP20 << endl;
   } else {
      cout << endl << fastNLO::_CSEP20C << " fastNLO Table: CoeffBase " << fastNLO::_CSEP20 << endl;
   }
   fastNLOTools::PrintVector(CtrbDescript,"Contribution description (CtrbDescript)","#");
   fastNLOTools::PrintVector(CodeDescript,"Code description (CodeDescript)","#");
   if ( abs(iprint) > 0 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 0) " << fastNLO::_SSEP20 << endl;
      printf(" #   IXsectUnits                       %d\n",IXsectUnits);
      printf(" #   IDataFlag                         %d\n",IDataFlag);
      printf(" #   IAddMultFlag                      %d\n",IAddMultFlag);
      printf(" #   IContrFlag1                       %d\n",IContrFlag1);
      printf(" #   IContrFlag2                       %d\n",IContrFlag2);
      printf(" #   NScaleDep                         %d\n",NScaleDep);
   }
   if ( iprint < 0 ) {
      cout << fastNLO::_CSEPSC << endl;
   } else {
      //      cout << fastNLO::_DSEPSC << endl;
   }
}


//________________________________________________________________________________________________________________ //

// Erase observable bin
void fastNLOCoeffBase::EraseBin(unsigned int iObsIdx) {
   info["fastNLOCoeffBase::EraseBin"]<<"Erasing table entries in CoeffBase for bin index " << iObsIdx << endl;
   SetNObsBin(GetNObsBin()-1);
}

// Multiply observable bin
void fastNLOCoeffBase::MultiplyBin(unsigned int iObsIdx, double nfact) {
   info["fastNLOCoeffBase::MultiplyBin"]<<"Multiplying table entries. Nothing to be done in CoeffBase." << endl;
}
