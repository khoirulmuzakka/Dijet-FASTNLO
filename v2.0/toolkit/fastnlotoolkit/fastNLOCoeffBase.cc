#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffBase.h"
#include "fastnlotk/fastNLOTools.h"
#include "fastnlotk/fastNLOGeneratorConstants.h"

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
   //! Use has to take care to delete this object later
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

   fastNLOTools::ReadMagicNo(table);

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
   //debug["EndReadCoeff"]<<endl;
   fastNLOTools::ReadMagicNo(table);
   fastNLOTools::PutBackMagicNo(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::Write(ostream& table) {
   say::debug["Write"]<<"Writing fastNLOCoeffBase."<<endl;
   table << tablemagicno << endl;
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
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   if( IXsectUnits != other.GetIXsectUnits() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   if( IDataFlag != other.GetIDataFlag() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   if( IAddMultFlag != other.GetIAddMultFlag() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   if( IContrFlag1 != other.GetIContrFlag1() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   if( IContrFlag2 != other.GetIContrFlag2() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   if( NScaleDep != other.GetNScaleDep() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   // chcekc descripts here ?!
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


void fastNLOCoeffBase::Print() const {
  printf("\n **************** FastNLO Table: CoeffBase ****************\n");
  printf(" B   Scenario::GetNObsBin()        %d\n",fNObsBins);
  printf(" B   IXsectUnits                   %d\n",IXsectUnits);
  printf(" B   IDataFlag                     %d\n",IDataFlag);
  printf(" B   IAddMultFlag                  %d\n",IAddMultFlag);
  printf(" B   IContrFlag1                   %d\n",IContrFlag1);
  printf(" B   IContrFlag2                   %d\n",IContrFlag2);
  printf(" B   NScaleDep                     %d\n",NScaleDep);
  fastNLOTools::PrintVector(CtrbDescript,"CtrbDescript","B");
  fastNLOTools::PrintVector(CodeDescript,"CodeDescript","B");
  printf(" *******************************************************\n");

}


//________________________________________________________________________________________________________________ //
