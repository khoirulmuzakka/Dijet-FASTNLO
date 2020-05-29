#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffBase.h"
#include "fastnlotk/fastNLOConstants.h"
#include "fastnlotk/fastNLOTools.h"

using namespace std;
using namespace fastNLO;

//________________________________________________________________________________________________________________ //
fastNLOCoeffBase::fastNLOCoeffBase(int NObsBin)
   : PrimalScream("fastNLOCoeffBase"), fNObsBins(NObsBin), IXsectUnits(),
     IDataFlag(), IAddMultFlag(), IContrFlag1(), IContrFlag2(), NScaleDep(),
     CtrbDescript(), CodeDescript() {}

//________________________________________________________________________________________________________________ //
fastNLOCoeffBase* fastNLOCoeffBase::Clone() const {
   //! User has to take care to delete this object later
   return new fastNLOCoeffBase(*this);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::Read(istream& table, int ITabVersionRead){
   // basic read function.
   // reads in only 'base'-variables.
   debug["ReadCoeffBase::Read"]<<"Start reading table ..."<<endl;
   ReadBase(table, ITabVersionRead);
   ReadCoeffInfoBlocks(table, ITabVersionRead);
   EndReadCoeff(table, ITabVersionRead);
   debug["ReadCoeffBase::Read"]<<"Finished reading table."<<endl;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::ReadBase(istream& table, int ITabVersionRead){
   debug["ReadCoeffBase::ReadBase"]<<"Start reading base coefficient table ..."<<endl;
   fastNLOTools::ReadMagicNo(table);
   table >> IXsectUnits;
   table >> IDataFlag;
   table >> IAddMultFlag;
   table >> IContrFlag1;
   table >> IContrFlag2;
   table >> NScaleDep;
   fastNLOTools::ReadFlexibleVector(CtrbDescript,table);
   fastNLOTools::ReadFlexibleVector(CodeDescript,table);
   debug["ReadCoeffBase::ReadBase"]<<"Finished reading base coefficient table."<<endl;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::EndReadCoeff(istream& table, int ITabVersionRead){
   debug["fastNLOCoeffBase::EndReadCoeff"]<<"Should have reached end of coefficient table for table version "<<ITabVersionRead<<endl;
   fastNLOTools::ReadMagicNo(table);
   fastNLOTools::PutBackMagicNo(table);
   debug["fastNLOCoeffBase::EndReadCoeff"]<<"Finished reading coefficient table for table version "<<ITabVersionRead<<endl;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::Write(ostream& table, int ITabVersionWrite) {
   say::debug["Write"]<<"Writing fastNLOCoeffBase for table version " << ITabVersionWrite << "." << endl;
   table << fastNLO::tablemagicno << sep;
   // if ( itabversion >= 24000 ) table << fastNLO::tabversion << sep;
   // if ( itabversion >= 24000 ) table << "fastNLO_CoeffAddBase" << sep;
   // if ( itabversion >= 24000 ) table << 0 << sep; // v2.4, but yet unused
   table << IXsectUnits << sep;
   table << IDataFlag << sep;
   table << IAddMultFlag << sep;
   table << IContrFlag1 << sep;
   table << IContrFlag2 << sep;
   table << NScaleDep << sep;
   fastNLOTools::WriteFlexibleVector(CtrbDescript,table);
   fastNLOTools::WriteFlexibleVector(CodeDescript,table);
   // if ( itabversion >= 24000 ) table << 0 << sep; // v2.4, but yet unuse
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
bool fastNLOCoeffBase::IsCatenable(const fastNLOCoeffBase& other) const {
   if( IXsectUnits != other.GetIXsectUnits() ){
      debug["IsCatenable"]<<"IXsectUnits != other.GetIXsectUnits(). Skipped."<<endl;
      return false;
   }
   if( IDataFlag != other.GetIDataFlag() ){
      debug["IsCatenable"]<<"IDataFlag != other.GetIDataFlag(). Skipped."<<endl;
      return false;
   }
   if( IAddMultFlag != other.GetIAddMultFlag() ){
      debug["IsCatenable"]<<"IAddMultFlag != other.GetIAddMultFlag(). Skipped."<<endl;
      return false;
   }
   if( IContrFlag1 != other.GetIContrFlag1() ){
      debug["IsCatenable"]<<"IContrFlag1 != other.GetIContrFlag1(). Skipped."<<endl;
      return false;
   }
   if( IContrFlag2 != other.GetIContrFlag2() ){
      debug["IsCatenable"]<<"IContrFlag2 != other.GetIContrFlag2(). Skipped."<<endl;
      return false;
   }
   if( NScaleDep != other.GetNScaleDep() ){
      debug["IsCatenable"]<<"NScaleDep != other.GetNScaleDep(). Skipped."<<endl;
      return false;
   }
   info["IsCatenable"]<<"Base parameters of contribution allow catenation"<<endl;
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
   if ( std::abs(iprint) > 0 ) {
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
   debug["fastNLOCoeffBase::EraseBin"]<<"Erasing table entries in CoeffBase for bin index " << iObsIdx << endl;
   SetNObsBin(GetNObsBin()-1);
}

// Catenate observable bin
void fastNLOCoeffBase::CatBin(const fastNLOCoeffBase& other, unsigned int iObsIdx) {
   debug["fastNLOCoeffBase::CatBin"]<<"Catenating observable bin in CoeffBase corresponding to bin index " << iObsIdx << endl;
   SetNObsBin(GetNObsBin()+1);
}

// Multiply observable bin
void fastNLOCoeffBase::MultiplyBin(unsigned int iObsIdx, double nfact) {
   debug["fastNLOCoeffBase::MultiplyBin"]<<"Multiplying table entries. Nothing to be done in CoeffBase." << endl;
}

//
//________________________________________________________________________________________________________________
// Added to include CoeffInfoBlocks
bool fastNLOCoeffBase::HasCoeffInfoBlock(int fICoeffInfoBlockFlag1) {
   bool result = false;
   for (int i=0; i<NCoeffInfoBlocks; i++) {
      if ( ICoeffInfoBlockFlag1[i] == fICoeffInfoBlockFlag1 ) result = true;
   }
   return result;
}

bool fastNLOCoeffBase::HasCoeffInfoBlock(int fICoeffInfoBlockFlag1, int fICoeffInfoBlockFlag2) {
   bool result = false;
   for (int i=0; i<NCoeffInfoBlocks; i++) {
      if ( ICoeffInfoBlockFlag1[i] == fICoeffInfoBlockFlag1 && ICoeffInfoBlockFlag2[i] == fICoeffInfoBlockFlag2 ) result = true;
   }
   return result;
}

int fastNLOCoeffBase::GetCoeffInfoBlockIndex(int fICoeffInfoBlockFlag1) {
   for (int i=0; i<NCoeffInfoBlocks; i++) {
      if ( ICoeffInfoBlockFlag1[i] == fICoeffInfoBlockFlag1 ) return i;
   }
   return -1;
}

int fastNLOCoeffBase::GetCoeffInfoBlockIndex(int fICoeffInfoBlockFlag1, int fICoeffInfoBlockFlag2) {
   for (int i=0; i<NCoeffInfoBlocks; i++) {
      if ( ICoeffInfoBlockFlag1[i] == fICoeffInfoBlockFlag1 && ICoeffInfoBlockFlag2[i] == fICoeffInfoBlockFlag2 ) return i;
   }
   return -1;
}

void fastNLOCoeffBase::ReadCoeffInfoBlocks(istream& table, int ITabVersionRead) {
   if (ITabVersionRead < 25000) {
      debug["fastNLOCoeffBase::ReadCoeffInfoBlocks"]<<"No additional info blocks allowed for table versions < 25000"<<endl;
   } else {
      table >> NCoeffInfoBlocks;
      debug["fastNLOCoeffBase::ReadCoeffInfoBlocks"]<<"Found "<<NCoeffInfoBlocks<<" additional info blocks for coefficient table."<<endl;
      for (int i=0; i<NCoeffInfoBlocks; i++) {
         int iflag;
         table >> iflag;
         ICoeffInfoBlockFlag1.push_back(iflag);
         if (ICoeffInfoBlockFlag1[i] == 0) {
            debug["fastNLOCoeffBase::ReadCoeffInfoBlocks"]<<"Found info block of type ICoeffInfoBlockFlag1 = "<<ICoeffInfoBlockFlag1[i]<<endl;
         } else {
            error["fastNLOCoeffBase::ReadCoeffInfoBlocks"]<<"Found info block of unknown type ICoeffInfoBlockFlag1 = "<<ICoeffInfoBlockFlag1[i]<<endl;
            exit(111);
         }
         table >> iflag;
         ICoeffInfoBlockFlag2.push_back(iflag);
         if (ICoeffInfoBlockFlag2[i] == 0) {
            debug["fastNLOCoeffBase::ReadCoeffInfoBlocks"]<<"Found info block of type ICoeffInfoBlockFlag2 = "<<ICoeffInfoBlockFlag2[i]<<endl;
         } else {
            error["fastNLOCoeffBase::ReadCoeffInfoBlocks"]<<"Found info block of unknown type ICoeffInfoBlockFlag2 = "<<ICoeffInfoBlockFlag2[i]<<endl;
            exit(222);
         }
         std::vector < std::string > Description;
         fastNLOTools::ReadFlexibleVector(Description,table);
         CoeffInfoBlockDescript.push_back(Description);
         if ( ICoeffInfoBlockFlag1[i] == 0 ) { // Entry per ObsBin
            std::vector < double > Content;
            int nlines = fastNLOTools::ReadFlexibleVector(Content,table);
            if ( nlines-1 != fNObsBins ) {
               error["fastNLOCoeffBase::ReadCoeffInfoBlocks"]<<"Found info block of type ICoeffInfoBlockFlag1 = "<<ICoeffInfoBlockFlag1[i]<<" , but # of content lines = "<<nlines-1<<" differs from fNObsBins = "<<fNObsBins<<"! Aborted."<<endl;
               exit(223);
            } else {
               debug["fastNLOCoeffBase::ReadCoeffInfoBlocks"]<<"Read "<<nlines-1<<" lines into InfoBlock content vector."<<endl;
            }
            CoeffInfoBlockContent.push_back(Content);
         }
      }
   }
}

void fastNLOCoeffBase::WriteCoeffInfoBlocks(ostream& table, int ITabVersionWrite) {
   if (ITabVersionWrite < 25000) {
      debug["fastNLOCoeffBase::WriteCoeffInfoBlocks"]<<"No additional info blocks allowed for table versions < 25000"<<endl;
   } else {
      // Test with zero InfoBlocks
      debug["fastNLOCoeffBase::WriteCoeffInfoBlocks"]<<"Writing additional line "<<NCoeffInfoBlocks<<endl;
      table << NCoeffInfoBlocks << sep;
      for (int i=0; i<NCoeffInfoBlocks; i++) {
         table << ICoeffInfoBlockFlag1[i];
         table << ICoeffInfoBlockFlag2[i];
         table << NCoeffInfoBlockDescr[i];
         for (int j=0; j<NCoeffInfoBlockDescr[i];j++) {
            table << CoeffInfoBlockDescript[i][j];
         }
         int nlines = fastNLOTools::WriteFlexibleVector(CoeffInfoBlockContent[i],table);
      }
   }
}
