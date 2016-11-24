#include "config.h"
#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <unistd.h>
#include <set>
#include "fastnlotk/fastNLOTable.h"
#include "fastnlotk/fastNLOTools.h"


using namespace std;

// ___________________________________________________________________________________________________
bool fastNLOTable::fWelcomeOnce = false;


// ___________________________________________________________________________________________________
//fastNLOTable::fastNLOTable() : PrimalScream("fastNLOTable") {
fastNLOTable::fastNLOTable() : fPrecision(8), logger("fastNLOTable") {
   if (!fWelcomeOnce) PrintWelcomeMessage();
}


// ___________________________________________________________________________________________________
fastNLOTable::fastNLOTable(string name)  : ffilename(name), fPrecision(8), logger("fastNLOTable")  {
   //logger.SetClassName("fastNLOTable");
   if (!fWelcomeOnce) PrintWelcomeMessage();
   ReadTable();
}


// ___________________________________________________________________________________________________
fastNLOTable::~fastNLOTable(){
   // delete fCoeff tables...
   DeleteAllCoeffTable();
}

// ___________________________________________________________________________________________________
fastNLOTable::fastNLOTable(const fastNLOTable& other)
   : ffilename(other.ffilename), fPrecision(other.fPrecision),
     Itabversion(other.Itabversion), ScenName(other.ScenName),
     //Ncontrib(other.Ncontrib),
     //Nmult(other.Nmult),
     //Ndata(other.Ndata),
     logger("fastNLOTable"),
    fCoeff(other.fCoeff.size()),
    Ecms(other.Ecms), ILOord(other.ILOord), Ipublunits(other.Ipublunits),
    ScDescript(other.ScDescript), NObsBin(other.NObsBin), NDim(other.NDim),
    DimLabel(other.DimLabel), IDiffBin(other.IDiffBin), Bin(other.Bin),
    BinSize(other.BinSize), INormFlag(other.INormFlag),
    DenomTable(other.DenomTable), IDivLoPointer(other.IDivLoPointer),
    IDivUpPointer(other.IDivUpPointer)
{
   //! Copy constructor
   logger.SetClassName("fastNLOTable");
   for (size_t i = 0; i < other.fCoeff.size(); ++i) {
      fCoeff[i] = other.fCoeff[i]->Clone();
   }
}


// ___________________________________________________________________________________________________
void fastNLOTable::DeleteAllCoeffTable(){
   for (size_t i = 0; i < fCoeff.size(); ++i) {
      delete fCoeff[i];
   }
   fCoeff.clear();
}


// ___________________________________________________________________________________________________
void fastNLOTable::ReadTable(){
   //! Read file
   zstr::ifstream* strm = OpenFileRead();
   // read header
   logger.debug["ReadTabl"]<<"Reading header."<<endl;
   int nCoeff = ReadHeader(*strm);
   // read scenario
   logger.debug["ReadTabl"]<<"Reading scenario."<<endl;
   ReadScenario(*strm);
   // read b-blocks
   logger.debug["ReadTabl"]<<"Reading coeff tables."<<endl;
   ReadCoeffTables(*strm, nCoeff);
   // close stream
   logger.debug["ReadTabl"]<<"Reading done closing files."<<endl;
   CloseFileRead(*strm);
}


//______________________________________________________________________________
int fastNLOTable::ReadHeader(istream& table) {
   //!<
   //!< Read table header (formely named BlockA1 and BlockA2)
   //!< return number of contributions to follow
   //!<
   table.peek();
   if (table.eof()) {
      logger.error["ReadHeader"]<<"Cannot read from stream."<<endl;
   }

   if (!fastNLOTools::ReadMagicNo(table)) {
#ifndef HAVE_LIBZ
      // check if filename ends with .gz
      std::string ending = ".gz";
      if (ffilename.length() >= ending.length() && ffilename.compare(ffilename.length() - ending.length(), ending.length(), ending) == 0) {
         logger.error["ReadScenario"]<<"Input file has a .gz file extension but zlib support is not enabled!"<<endl;
      }
#endif
      logger.error["ReadHeader"]<<"Did not find initial magic number, aborting!"<<endl;
      logger.error["ReadHeader"]<<"Please check compatibility of tables and program version!"<<endl;
      exit(1);
   }
   table >> Itabversion;
   table >> ScenName;
   std::string test;
   std::getline(table,test);
   if ( test != "" )  {
      logger.warn["ReadHeader"]<<"Scenario name is not allowed to contain white spaces!!"<<endl;
   }
   // check if ScenName contains spaces
   int unused, Ncontrib,Ndata;
   table >> Ncontrib;
   table >> unused; //Nmult;
   table >> Ndata;
   table >> unused; //NuserString;
   if ( unused>0) logger.warn["ReadHeader"]<<"Reading 'userInt' (NuserString="<<unused<<"). This is not usual."<<endl;
   table >> unused; //NuserInt;
   if ( unused>0) logger.warn["ReadHeader"]<<"Reading 'userInt' (NuserInt="<<unused<<"). This is not usual."<<endl;
   for (int i = 0 ; i<unused ; i++) {
      table >> unused; //IUserLines;
      // future code if 'user-blocks' are used ...
      // int NUserFlag;
      // string NUserBlockDescr;
      // table >> NUserFlag;
      // table >> NUserBlockDescr;;
      // if ( known-user-block ) { read-known-userblock... }
      // else { // skip meaningful reading
      //    for ( int i = 2 ; i<NuserInt ; i++ ) {
      //       double devnull;
      //       table >> devnull;
      //    }
      // }
      // ...sofar skip reading
      for (int i = 0 ; i<unused ; i++) {
         table >> unused;//devnull;
      }
   }
   table >> unused;//NuserFloat;
   if ( unused>0) logger.warn["ReadHeader"]<<"Reading 'userFloat' (NuserInt="<<unused<<"). This is not usual."<<endl;
   table >> unused;//Imachine;
   if (!fastNLOTools::ReadMagicNo(table)) {
      Print(1);
      logger.error["ReadHeader"]<<"Did not find final magic number, aborting!"<<endl;
      logger.error["ReadHeader"]<<"Please check compatibility of tables and program version!"<<endl;
      exit(1);
   }
   fastNLOTools::PutBackMagicNo(table);
   return Ncontrib+Ndata;
}


// ___________________________________________________________________________________________________
void fastNLOTable::ReadCoeffTables(istream& table, int nCoeff){
   //!< read nCoeff Coefficient tables (additive, multiplicative and data)
   for (int i=0; i<nCoeff; i++) {
      fastNLOCoeffBase cTemp(NObsBin);
      cTemp.ReadBase(table);
      fastNLOCoeffBase* cN = ReadRestOfCoeffTable(cTemp, table);
      CreateCoeffTable(i, cN);
   }
}


// ___________________________________________________________________________________________________
fastNLOCoeffBase* fastNLOTable::ReadRestOfCoeffTable(const fastNLOCoeffBase& cB, istream& table){
   // take coeffbase and identify type of contribution.
   //  - create instance of correct full coefficient table
   //  - read in 'rest' of coeff table

   // identify coeff-table:
   bool quiet = true;
   if ( fastNLOCoeffData::CheckCoeffConstants(&cB,quiet) ) {
      logger.debug["ReadRestOfCoeffTable"]<<"Found data table. Now reading in."<<endl;
      fastNLOCoeffData* cN = new fastNLOCoeffData(cB);
      cN->ReadRest(table);
      return cN;
   } else if ( fastNLOCoeffMult::CheckCoeffConstants(&cB,quiet) ) {
      logger.debug["ReadRestOfCoeffTable"]<<"Found multiplicative contribution. Now reading in."<<endl;
      fastNLOCoeffMult* cN = new fastNLOCoeffMult(cB);
      cN->ReadRest(table);
      return cN;
   } else if ( fastNLOCoeffAddFix::CheckCoeffConstants(&cB,quiet) ) {
      logger.debug["ReadRestOfCoeffTable"]<<"Found additive fixed order contribution (v2.0). Now reading in."<<endl;
      fastNLOCoeffAddFix* cN = new fastNLOCoeffAddFix(cB);
      cN->ReadRest(table);
      return cN;
   } else if ( fastNLOCoeffAddFlex::CheckCoeffConstants(&cB,quiet) ) {
      logger.debug["ReadRestOfCoeffTable"]<<"Found additive flexible scale contribution. Now reading in."<<endl;
      fastNLOCoeffAddFlex* cN = new fastNLOCoeffAddFlex(cB,ILOord);
      cN->ReadRest(table);
      return cN;
   } else {
      logger.error["ReadRestOfCoeffTable"]<<"Could not identify coefficient table. Print and exiting ... "<<endl;
      cB.Print(5);
      exit(1);
   }
   return NULL;
}


// ___________________________________________________________________________________________________
void fastNLOTable::WriteTable() {
   //!<
   //!< WriteTable(). writes the full FastNLO table to
   //!< the previously defined ffilename on disk.
   //!< Write fastNLO table to file 'ffilename' (member)
   std::string extension = ".gz";
   bool compress = false;
   if ((ffilename.length() >= extension.length()) and 
         (ffilename.compare(ffilename.length() - extension.length(), extension.length(), extension) == 0)) 
   {
      logger.info["WriteTable"]<<"Filename ends with .gz, therefore enable compression." << endl;
      compress = true; 
   }

   logger.info["WriteTable"]<<"Writing fastNLO table with " << GetNcontrib() << " theory contributions to file: " << ffilename << endl;
   std::ostream* table = OpenFileWrite(compress);
   logger.debug["WriteTable"]<<"Writing table header to file ..."<<endl;
   WriteHeader(*table);
   logger.debug["WriteTable"]<<"Writing scenario to file ..."<<endl;
   WriteScenario(*table);
   for(int i=0;i<GetNcontrib()+GetNdata();i++){
      logger.debug["WriteTable"]<<"Writing coefficient table #"<<i<<endl;
      GetCoeffTable(i)->Write(*table);
   }
   CloseFileWrite(*table);
}


// ___________________________________________________________________________________________________
void fastNLOTable::WriteTable(string filename) {
   //! Write fastNLO table to file 'filename'
   string tempfilename = ffilename;
   SetFilename(filename);
   WriteTable();
   SetFilename(tempfilename);
}



//______________________________________________________________________________
void fastNLOTable::WriteHeader(ostream& table) {
   table << fastNLO::tablemagicno << endl;
   table << Itabversion << endl;
   if ( ScenName.find(" ")!=string::npos )  {
      logger.warn["WriteHeader"]<<"Scenario name is not allowed to contain white spaces!!"<<endl;
      ScenName = ScenName.substr(0,ScenName.find(" "));
      logger.warn["WriteHeader"]<<"Write ScenarioName: "<<ScenName<<endl;
   }
   table << ScenName << endl;
   // table << Ncontrib << endl;
   // table << Nmult << endl;
   // table << Ndata << endl;
   table << GetNcontrib() << endl;
   table << GetNmult() << endl;
   table << GetNdata() << endl;
   table << 0 << endl; //NUserString
   table << 0 << endl; //NuserInt
   table << 0 << endl; //NuserFloat
   table << 0 << endl; //Imachine
}


// ___________________________________________________________________________________________________
void fastNLOTable::ReadScenario(istream& table){
   table.peek();
   if (table.eof()){
      logger.warn["ReadScenario"]<<"Cannot read from file."<<endl;
   }

   if (!fastNLOTools::ReadMagicNo(table)) {
      logger.error["ReadScenario"]<<"Did not find initial magic number, aborting!"<<endl;
      logger.error["ReadScenario"]<<"Please check compatibility of tables and program version!"<<endl;
      exit(1);
   }

   table >> Ipublunits;
   size_t NScDescript = 0;
   table >> NScDescript;
   ScDescript.resize(NScDescript);
   char buffer[257];
   table.getline(buffer,256);
   for(size_t i=0;i<NScDescript;i++){
      table.getline(buffer,256);
      ScDescript[i] = buffer;
      //      StripWhitespace(ScDescript[i]);
   }

   table >> Ecms;
   table >> ILOord;
   table >> NObsBin;
   table >> NDim;
   DimLabel.resize(NDim);
   table.getline(buffer,256);
   for(int i=NDim-1;i>=0;i--){
      table.getline(buffer,256);
      DimLabel[i] = buffer;
   }

   IDiffBin.resize(NDim);
   for(int i=NDim-1;i>=0;i--){
      table >>  IDiffBin[i];
   }
   Bin.resize(NObsBin);
   for(unsigned int i=0;i<NObsBin;i++){
      Bin[i].resize(NDim);
      for(int j=NDim-1;j>=0;j--){
         table >> Bin[i][j].first;
         if (IDiffBin[j]==0 || IDiffBin[j]==2) {
            table >> Bin[i][j].second;
         } else {
            // For point-wise differential, IDiffBin = 1, set UpBin equal to LoBin
            Bin[i][j].second = Bin[i][j].first;
         }
      }
   }

   BinSize.resize(NObsBin);
   for(unsigned int i=0;i<NObsBin;i++){
      table >> BinSize[i];
   }

   table >> INormFlag;
   if( INormFlag < 0 ){
      table >> DenomTable;
   }
   if( INormFlag != 0 ){
      IDivLoPointer.resize(NObsBin);
      IDivUpPointer.resize(NObsBin);
      for(unsigned int i=0;i<NObsBin;i++){
         table >> IDivLoPointer[i];
         table >> IDivUpPointer[i];
      }
   }

   if (!fastNLOTools::ReadMagicNo(table)) {
      logger.error["ReadScenario"]<<"Did not find final magic number, aborting!"<<endl;
      logger.error["ReadScenario"]<<"Please check compatibility of tables and program version!"<<endl;
      exit(1);
   }
   fastNLOTools::PutBackMagicNo(table);
}


// ___________________________________________________________________________________________________
void fastNLOTable::WriteScenario(ostream& table){
   table << fastNLO::tablemagicno << endl;
   table << Ipublunits << endl;
   size_t NScDescript = ScDescript.size();
   table << NScDescript << endl;
   for(size_t i=0;i<NScDescript;i++){
      table << ScDescript[i] << endl;
   }
   table << Ecms << endl;
   table << ILOord << endl;
   logger.debug["WriteScenario"]<<"Writing NObsBin to be "<<NObsBin<<endl;
   table << NObsBin << endl;
   table << NDim << endl;
   for(int i=NDim-1;i>=0;i--){
      table << DimLabel[i] << endl;
   }
   for(int i=NDim-1;i>=0;i--){
      table << IDiffBin[i] << endl;
   }
   logger.debug["WriteScenario"]<<"Bin border size is "<<Bin.size()<<endl;
   for(unsigned int i=0;i<NObsBin;i++){
      for(int j=NDim-1;j>=0;j--){
         table <<  Bin[i][j].first  << endl;
         if(IDiffBin[j]==0 || IDiffBin[j]==2) table <<  Bin[i][j].second  << endl;
      }
   }
   for(unsigned int i=0;i<NObsBin;i++){
     table << BinSize[i]  << endl;
   }

   table << INormFlag << endl;
   if( INormFlag < 0 ){
      table << DenomTable << endl;
   }
   if( INormFlag != 0 ){
      for(unsigned int i=0;i<NObsBin;i++){
         table << IDivLoPointer[i] << endl;
         table << IDivUpPointer[i] << endl;
      }
   }
}


// ___________________________________________________________________________________________________
bool fastNLOTable::IsCompatible(const fastNLOTable& other) const {
   if ( !IsCompatibleHeader(other) ) return false;
   if ( !IsCompatibleScenario(other) ) return false;
   logger.info["IsCompatible"]<<"Tables seem to be compatible for merging/appending. Continuing."<<endl;
   return true;
}


// ___________________________________________________________________________________________________
bool fastNLOTable::IsCompatibleScenario(const fastNLOTable& other) const {
   bool potentialcompatible = true;
   if(Ipublunits != other.Ipublunits){
      logger.warn["IsCompatibleScenario"]<<"Differing cross section units found: "<<Ipublunits<<" and "<<other.Ipublunits<<endl;
      return false;
   }
   if(ScDescript != other.ScDescript){
      logger.warn["IsCompatibleScenario"]<<"Differing scenario description found."<<endl;
      potentialcompatible = false;
   }
   if(!cmp(Ecms,other.Ecms)){
      logger.warn["IsCompatibleScenario"]<<"Differing center-of-mass energy found: "<<Ecms<<" and "<<other.Ecms<<endl;
      return false;
   }
   if(ILOord != other.ILOord){
      logger.warn["IsCompatibleScenario"]<<"Differing ILOord found: "<<ILOord<<" and "<<other.GetLoOrder()<<endl;
      return false;
   }
   if(NObsBin != other.NObsBin){
      logger.warn["IsCompatibleScenario"]<<"Differing NObsBin found: "<<NObsBin<<" and "<<other.NObsBin<<endl;
      return false;
   }
   if(NDim != other.NDim){
      logger.warn["IsCompatibleScenario"]<<"Differing NDim found: "<<NDim<<" and "<<other.NDim<<endl;
      return false;
   }
   if(DimLabel != other.DimLabel){
      logger.warn["IsCompatibleScenario"]<<"Differing label of observables found."<<endl;
      potentialcompatible = false;
   }
   if(IDiffBin != other.IDiffBin){
      logger.warn["IsCompatibleScenario"]<<"Differing IDiffBin found."<<endl;
      return false;
   }
   if(!cmp(Bin,other.Bin)){
      logger.warn["IsCompatibleScenario"]<<"Differing Bin boundaries found."<<endl;
      return false;
   }
   if(!cmp(BinSize,other.BinSize)){
      logger.warn["IsCompatibleScenario"]<<"Differing bin sizes found."<<endl;
      return false;
   }
   if(INormFlag != other.INormFlag){
      logger.warn["IsCompatibleScenario"]<<"Differing INormFlag found: "<<INormFlag<<" and "<<other.INormFlag<<endl;
      return false;
   }
   if(INormFlag<0){
      if(DenomTable != other.DenomTable){
         logger.warn["IsCompatibleScenario"]<<"Differing DenomTable found."<<endl;
         return false;
      }
   }
   if(INormFlag!=0){
      for(unsigned int i=0;i<NObsBin;i++){
         if(IDivLoPointer[i] != other.IDivLoPointer[i]){
            logger.warn["IsCompatibleScenario"]<<"Differing IDivLoPointer["<<i<<"] found"<<endl;
            return false;
         }
         if(IDivUpPointer[i] != other.IDivUpPointer[i]){
            logger.warn["IsCompatibleScenario"]<<"Differing IDivUpPointer["<<i<<"] found."<<endl;
            return false;
         }
      }
   }
   if ( !potentialcompatible ) logger.warn["IsCompatibleScenario"]<<"Some labels have differing values, but relevant variables seem to be compatible. Continuing."<<endl;
   return true;
}


// ___________________________________________________________________________________________________
bool fastNLOTable::IsCatenable(const fastNLOTable& other) const {
   if ( !IsCatenableHeader(other) ) return false;
   if ( !IsCatenableScenario(other) ) return false;

   const bool quiet = true;
   const int nc = other.GetNcontrib() + other.GetNdata();
   // loop over all contributions from 'other'-table to check catenability
   int matches[nc];
   for ( int ic=0 ; ic<nc; ic++ ) {
      matches[ic] = 0;
      // check against all contributions from 'this'-table to check catenability
      for (unsigned int j = 0; j<fCoeff.size() ; j++) {
         fastNLOCoeffBase* cother = (fastNLOCoeffBase*)other.GetCoeffTable(ic);
         // data?
         if ( fastNLOCoeffData::CheckCoeffConstants(cother,quiet) ) {
            fastNLOCoeffData* cdatthis  = (fastNLOCoeffData*)fCoeff[j];
            fastNLOCoeffData* cdatother = (fastNLOCoeffData*)other.GetCoeffTable(ic);
            if ( cdatthis->IsCatenable(*cdatother) ) {
               matches[ic]++;
               continue;
            }
         }
         // multiplicative?
         else if ( fastNLOCoeffMult::CheckCoeffConstants(cother,quiet) ) {
            fastNLOCoeffMult* cmultthis  = (fastNLOCoeffMult*)fCoeff[j];
            fastNLOCoeffMult* cmultother = (fastNLOCoeffMult*)other.GetCoeffTable(ic);
            if ( cmultthis->IsCatenable(*cmultother) ) {
               matches[ic]++;
               continue;
            }
         }
         // additive?
         else if ( fastNLOCoeffAddBase::CheckCoeffConstants(cother,quiet) ) {
            if ( fastNLOCoeffAddFix::CheckCoeffConstants(cother,quiet) ) {
               fastNLOCoeffAddFix* cfixthis  = (fastNLOCoeffAddFix*)fCoeff[j];
               fastNLOCoeffAddFix* cfixother = (fastNLOCoeffAddFix*)other.GetCoeffTable(ic);
               if ( cfixthis->IsCatenable(*cfixother) ) {
                  matches[ic]++;
                  continue;
               }
            } else if ( fastNLOCoeffAddFlex::CheckCoeffConstants(cother,quiet) ) {
               fastNLOCoeffAddFlex* cflexthis  = (fastNLOCoeffAddFlex*)fCoeff[j];
               fastNLOCoeffAddFlex* cflexother = (fastNLOCoeffAddFlex*)other.GetCoeffTable(ic);
               if ( cflexthis->IsCatenable(*cflexother) ) {
                  matches[ic]++;
                  continue;
               }
            }
         } else {
            logger.error["IsCatenable"] << "Unknown contribution found. Aborted!" <<endl;
            exit(1);
         }
      }
   }

   // check match count; all numbers must be unity
   bool catenable = true;
   for ( int ic=0 ; ic<nc; ic++ ) {
      if ( matches[ic] != 1 ) {
         catenable = false;
         logger.warn["IsCatenable"] << "Table contributions do not match. Catenation of observable bins not possible!" <<endl;
         break;
      }
   }
   if ( catenable ) {
      logger.info["IsCatenable"]<<"Table contributions seem to be compatible for catenating observable bins. Continuing."<<endl;
   }
   return catenable;
}


// ___________________________________________________________________________________________________
bool fastNLOTable::IsCatenableScenario(const fastNLOTable& other) const {
   bool potentialcatenable = true;
   if(Ipublunits != other.Ipublunits){
      logger.warn["IsCatenableScenario"]<<"Differing cross section units found: "<<Ipublunits<<" and "<<other.Ipublunits<<endl;
      return false;
   }
   if(ScDescript != other.ScDescript){
      logger.warn["IsCatenableScenario"]<<"Differing scenario description found. Only the first one is kept."<<endl;
      potentialcatenable = false;
   }
   if(!cmp(Ecms,other.Ecms)){
      logger.warn["IsCatenableScenario"]<<"Differing center-of-mass energy found: "<<Ecms<<" and "<<other.Ecms<<endl;
      return false;
   }
   if(ILOord != other.ILOord){
      logger.warn["IsCatenableScenario"]<<"Differing ILOord found: "<<ILOord<<" and "<<other.GetLoOrder()<<endl;
      return false;
   }
   if(NDim != other.NDim){
      logger.warn["IsCatenableScenario"]<<"Differing NDim found: "<<NDim<<" and "<<other.NDim<<endl;
      return false;
   }
   if(DimLabel != other.DimLabel){
      logger.warn["IsCatenableScenario"]<<"Differing label of observables found."<<endl;
      potentialcatenable = false;
   }
   if(IDiffBin != other.IDiffBin){
      logger.warn["IsCatenableScenario"]<<"Differing IDiffBin found."<<endl;
      return false;
   }
   if(INormFlag != other.INormFlag){
      logger.warn["IsCatenableScenario"]<<"Differing INormFlag found: "<<INormFlag<<" and "<<other.INormFlag<<endl;
      return false;
   }
   if(INormFlag<0){
      if(DenomTable != other.DenomTable){
         logger.warn["IsCatenableScenario"]<<"Differing DenomTable found."<<endl;
         return false;
      }
   }
   if ( !potentialcatenable ) logger.warn["IsCatenableScenario"]<<"Some labels have differing values, but relevant variables seem to be catenable. Continuing."<<endl;
   return true;
}


// ___________________________________________________________________________________________________
void fastNLOTable::AddTable(const fastNLOTable& other) {
   // Add another table to this table.
   // Either increase statistics of existing fixed-order contribution or
   // add further contributions (or both, if many tables are merged)
   //
   if ( !IsCompatible(other) ) {
      logger.error["AddTable"]<<"Tried to add/merge incompatible tables. Aborted!"<<endl;
      exit(1);
   }

   // These are counters for the newly read other table ...,
   const int ntot = other.GetNcontrib() + other.GetNdata();
   // but we need also to bookkeep this for the current table!
   int newnc = fCoeff.size();
   int newnd = 0;
   bool quiet = true;
   for ( unsigned int jc=0; jc<fCoeff.size(); jc++) {
      fastNLOCoeffBase* cthis = (fastNLOCoeffBase*)fCoeff[jc];
      if ( fastNLOCoeffData::CheckCoeffConstants(cthis,quiet) ) {
         newnc--;
         newnd++;
      }
   }
   // Loop over all contributions from 'other'-table
   for ( int ic=0; ic<ntot; ic++ ) {
      logger.info["fastNLOTable::AddTable"]<<"Adding contribution no. " << ic << endl;
      bool wasAdded = false;
      // Find matching contribution from 'this'-table
      for ( unsigned int jc=0; jc<fCoeff.size(); jc++) {
         fastNLOCoeffBase* cother = (fastNLOCoeffBase*)other.GetCoeffTable(ic);
         // Identify type of other coeff table
         // Additive fixed-order?
         if ( fastNLOCoeffAddBase::CheckCoeffConstants(cother,quiet) ) {
            if ( fastNLOCoeffAddFix::CheckCoeffConstants(cother,quiet) ) {
               fastNLOCoeffAddFix* clhs = (fastNLOCoeffAddFix*)fCoeff[jc];
               fastNLOCoeffAddFix* crhs = (fastNLOCoeffAddFix*)other.GetCoeffTable(ic);
               if ( clhs->IsCompatible(*crhs) ) {
                  logger.info["AddTable"]<<"Found matching fix-scale additive contribution." << endl;
                  logger.debug["AddTable"]<<"Summing contribution "<<ic<<" to fCoeff #"<<jc<<endl;
                  clhs->Add(*crhs);
                  wasAdded = true;
               }
            }
            else if ( fastNLOCoeffAddFlex::CheckCoeffConstants(cother,quiet) ) {
               fastNLOCoeffAddFlex* clhs = (fastNLOCoeffAddFlex*)fCoeff[jc];
               fastNLOCoeffAddFlex* crhs = (fastNLOCoeffAddFlex*)other.GetCoeffTable(ic);
               if ( clhs->IsCompatible(*crhs) ) {
                  logger.info["AddTable"]<<"Found matching flex-scale additive contribution." << endl;
                  logger.debug["AddTable"]<<"Summing contribution "<<ic<<" to fCoeff #"<<jc<<endl;
                  clhs->Add(*crhs);
                  wasAdded = true;
               }
            }
         }
         // Multiplicative?
         else if ( fastNLOCoeffMult::CheckCoeffConstants(cother,quiet) ) {
            fastNLOCoeffMult* clhs = (fastNLOCoeffMult*)fCoeff[jc];
            fastNLOCoeffMult* crhs = (fastNLOCoeffMult*)other.GetCoeffTable(ic);
            if ( clhs->IsCompatible(*crhs) ) {
               logger.error["AddTable"]<<"Found matching multiplicative contribution. This is not allowed. Aborted!" << endl;
               wasAdded = true;
               exit(1);
            }
         }
         // Data?
         else if ( fastNLOCoeffData::CheckCoeffConstants(cother,quiet) ) {
            fastNLOCoeffData* clhs = (fastNLOCoeffData*)fCoeff[jc];
            fastNLOCoeffData* crhs = (fastNLOCoeffData*)other.GetCoeffTable(ic);
            if ( clhs->IsCompatible(*crhs) ) {
               logger.error["AddTable"]<<"Found matching data contribution. This is not allowed. Aborted!" << endl;
               wasAdded = true;
               exit(1);
            }
         }
         // Unknown
         else {
            logger.error["AddTable"]<<"Could not identify contribution. Print and abort!" << endl;
            cother->Print(-1);
            exit(1);
         }
      }
      // Couldn't find a corresponding contribution,
      // so add this contribution as new.
      if ( !wasAdded ) {
         logger.info["AddTable"]<<"Adding new contribution to table."<<endl;
         fastNLOCoeffBase* add = other.GetCoeffTable(ic);
         if ( fastNLOCoeffAddFix::CheckCoeffConstants(add,quiet) ) {
            add = new fastNLOCoeffAddFix((fastNLOCoeffAddFix&)*add);
            // Adjust new theory contribution counter
            newnc++;
         } else if ( fastNLOCoeffAddFlex::CheckCoeffConstants(add,quiet) ) {
            add = new fastNLOCoeffAddFlex((fastNLOCoeffAddFlex&)*add);
            newnc++;
         } else if ( fastNLOCoeffMult::CheckCoeffConstants(add,quiet) ) {
            add = new fastNLOCoeffMult((fastNLOCoeffMult&)*add);
            newnc++;
         } else if ( fastNLOCoeffData::CheckCoeffConstants(add,quiet) ) {
            add = new fastNLOCoeffData((fastNLOCoeffData&)*add);
            // Adjust new data counter
            newnd++;
         }
         CreateCoeffTable(fCoeff.size(),add);
      }
   }
   // Check # of coefficients
   if ( (int)fCoeff.size() != newnc + newnd ) {
      logger.error["AddTable"]<<"Sorry, I'm confused about the no. of contributions. Aborted!" << endl;
      logger.error["AddTable"]<<"newnc = " << newnc << ", newnd = " << newnd << ", fCoeff.size() = " << fCoeff.size() << endl;
      exit(1);
   }
   // Set nc and nd for current table to be written out eventually
   // SetNcontrib(newnc);
   // SetNdata(newnd);
}


// ___________________________________________________________________________________________________
void fastNLOTable::CatenateTable(const fastNLOTable& other) {
   // Catenate another table to this table.
   // All contributions must be identically defined.
   //
   static unsigned int table_count = 0;
   if ( !IsCatenable(other) ) {
      logger.error["CatenateTable"]<<"Tried to catenate incompatible tables. Aborted!"<<endl;
      exit(1);
   } else {
      table_count++;
   }
   for ( unsigned int iObs=0; iObs<other.GetNObsBin(); iObs++ ) {
      this->CatBinToTable(other,iObs,table_count);
   }
}


// ___________________________________________________________________________________________________
int fastNLOTable::CreateCoeffTable(int no, fastNLOCoeffBase *newblockb) {
   // Attention: Proper adaptation of Ncontrib and Ndata, which are set each time a table is read,
   // to the current value of the table in memory must be done in the calling routine!
   logger.debug["CreateCoeffTable"]<<"Old: Ncontrib = " << GetNcontrib() << ", Ndata = " << GetNdata() << ", fCoeff.size() = " << fCoeff.size() << endl;
   logger.debug["CreateCoeffTable"]<<"Creating coefficient table no. " << no << ", actual fCoeff.size() is: " << fCoeff.size() << endl;
   if ( (no+1) > (int)fCoeff.size() ) {
      fCoeff.resize(no+1);
      logger.debug["CreateCoeffTable"]<<"Creating new coefficient table no. " << no << endl;
   }
   fCoeff[no] = newblockb;
   //Ncontrib = fCoeff.size();
   return 0;
}


// ___________________________________________________________________________________________________
bool fastNLOTable::cmp(const double x1, const double x2) const {
   double norm = (x1>0.) ? x1 : 1.; // If x1 is 0, do not try to calculate relative deviation, use absolute
   return((fabs(x1-x2)/norm)<1e-7);
}

bool fastNLOTable::cmp(const vector<double>& x1,const vector<double>& x2) const {
   bool result = true;
   for(unsigned int i = 0; i<x1.size() ;i++ ){
      result &= cmp (x1[i],x2[i]);
   }
   return result;
}

bool fastNLOTable::cmp(const vector<vector<double> >& x1, const vector<vector<double> >& x2) const {
   bool result = true;
   for(unsigned int i = 0; i<x1.size() ;i++ ){
      result &=  cmp (x1[i],x2[i]);
   }
   return result;
}

bool fastNLOTable::cmp(const vector<vector<pair<double,double> > >& x1, const vector<vector<pair<double,double> > >& x2) const {
   bool result = true;
   for(unsigned int i = 0; i<x1.size() ;i++ ){
      for(unsigned int j = 0; j<x1[i].size() ;j++ ){
         result = result & (cmp(x1[i][j].first,x2[i][j].first) && cmp(x1[i][j].second,x2[i][j].second));
      }
   }
   return result;
}


// ___________________________________________________________________________________________________
void fastNLOTable::SetLoOrder(int LOOrd){
   ILOord = LOOrd;
}


// ___________________________________________________________________________________________________
void fastNLOTable::SetDimLabel( string label, unsigned int iDim , bool IsDiff ){
   //! Set label for dimension
   //! In this method, we also set IDiffBin.
   //! The IDiffBin flag defines, if this dimension is
   //!    0 (not differential, two bin borders required),
   //!    1 (pointwise differential, one value required), (not yet completely implemented)
   //!    2 (binwise differential, two bin borders required)
   //!    In case 2 the cross section is divided by the corresponding bin width in this dimension.
   //
   // TODO: KR: The IsDiff boolean should be changed into an int to accommodate IDiffBin 0,1,2
   //           possibility?!
   //
   // int iDim: counting starts from 0


   // check validity of call
   if ( NDim < iDim ) {
      logger.error["SetDimLabel"]<<"Sorry, you have only initialized "<<NDim<<" dimensions, but you want to label a dimension with number "<<iDim<<endl;
      exit(1);
   }
   if ( iDim < 1) {
      logger.error["SetDimLabel"]<<"The dimension must be a natural number. iDim="<<iDim<<endl;
      exit(1);
   }

   if ( DimLabel.size() != NDim ){
      logger.error["SetDimLabel"]<<"You have to call SetNumDiffBin with a reasonable number before."<<endl;
      exit(1);
   }

   DimLabel[iDim] = label;
   IDiffBin[iDim] = IsDiff ? 2 : 0 ;
}


// ___________________________________________________________________________________________________
int fastNLOTable::GetNmult() const {
   int ret = 0;
   for ( unsigned int i = 0 ; i<fCoeff.size() ; i++ )
      if ( (fCoeff[i]->GetIDataFlag()==0) && (fCoeff[i]->GetIAddMultFlag()==1) ) ret++;
   return ret;
}


// ___________________________________________________________________________________________________
int fastNLOTable::GetNcontrib() const {
   int ret = 0;
   for ( unsigned int i = 0 ; i<fCoeff.size() ; i++ )
      if ( (fCoeff[i]->GetIDataFlag()==0) /*&& (fCoeff[i]->GetIAddMultFlag()==0)*/ ) ret++;
   return ret;
}


// ___________________________________________________________________________________________________
int fastNLOTable::GetNdata() const {
   int ret = 0;
   for ( unsigned int i = 0 ; i<fCoeff.size() ; i++ )
      if ( (fCoeff[i]->GetIDataFlag()==1) && (fCoeff[i]->GetIAddMultFlag()==0) ) ret++;
   return ret;
}


// ___________________________________________________________________________________________________
fastNLOCoeffBase* fastNLOTable::GetCoeffTable(int no) const {
   if ( no >= (int)fCoeff.size() ){
      logger.warn["GetCoeffTable"]<<"There is no contribution with number "<<no<<" but only "<<fCoeff.size()<<". Returning null pointer."<<endl;
      return NULL;
   }
   else
      return fCoeff[no];
}


// ___________________________________________________________________________________________________
fastNLOCoeffData* fastNLOTable::GetDataTable() const {
   for (unsigned int i= 0; i<fCoeff.size() ; i++ ) {
      fastNLOCoeffBase* c = GetCoeffTable(i);
      if ( fastNLOCoeffData::CheckCoeffConstants(c,true) ) {
         return (fastNLOCoeffData*)c;
      }
   }
   return NULL;
}


// ___________________________________________________________________________________________________
fastNLOCoeffAddBase* fastNLOTable::GetReferenceTable(fastNLO::ESMOrder eOrder) const {
   for (unsigned int i= 0; i<fCoeff.size() ; i++ ){
      fastNLOCoeffBase* c = GetCoeffTable(i);
      if ( fastNLOCoeffAddBase::CheckCoeffConstants(c,true) ) {
         if ( ((fastNLOCoeffAddBase*)c)->IsReference() ) {
            if ( eOrder == fastNLO::kLeading && c->IsLO() )
               return (fastNLOCoeffAddBase*)c;
            else if ( eOrder == fastNLO::kNextToLeading && c->IsNLO() )
               return (fastNLOCoeffAddBase*)c;
            else if ( eOrder == fastNLO::kNextToNextToLeading && c->IsNNLO() )
               return (fastNLOCoeffAddBase*)c;
         }
      }
   }
   return NULL;
}


// ___________________________________________________________________________________________________
// Getters for binning structure
// ___________________________________________________________________________________________________


// ___________________________________________________________________________________________________
// Return lower bin bound for obs. bin iObs in dim. iDim
double fastNLOTable::GetObsBinLoBound(unsigned int iObs, unsigned int iDim) const {
   if ( ! (iObs < NObsBin) ) {
      logger.error["GetObsBinLoBound"]<<"Observable bin iObs " << iObs << " out of range, NObsBin = " << NObsBin << ", aborted!" << endl;
      exit(1);
   }
   if ( ! (iDim < NDim) ) {
      logger.error["GetObsBinLoBound"]<<"Dimension iDim " << iDim << " out of range, NDim = " << NDim << ", aborted!" << endl;
      exit(1);
   }
   return Bin[iObs][iDim].first;
}


// ___________________________________________________________________________________________________
// Return upper bin bound for obs. bin iObs in dim. iDim
double fastNLOTable::GetObsBinUpBound(unsigned int iObs, unsigned int iDim) const {
   if ( ! (iObs < NObsBin) ) {
      logger.error["GetObsBinUpBound"]<<"Observable bin iObs " << iObs << " out of range, NObsBin = " << NObsBin << ", aborted!" << endl;
      exit(1);
   }
   if ( ! (iDim < NDim) ) {
      logger.error["GetObsBinUpBound"]<<"Dimension iDim " << iDim << " out of range, NDim = " << NDim << ", aborted!" << endl;
      exit(1);
   }
   return Bin[iObs][iDim].second;
}


// ___________________________________________________________________________________________________
// Return vector of lower bin bounds in dim. iDim for all obs. bins
vector < double > fastNLOTable::GetObsBinsLoBounds(unsigned int iDim) const {
   if ( ! (iDim < NDim) ) {
      logger.error["GetObsBinsLoBounds"]<<"Dimension iDim " << iDim << " out of range, NDim = " << NDim << ", aborted!" << endl;
      exit(1);
   }
   // Get lower bin edge of all observable bins for dimension 'iDim'
   vector < double > LoBin;
   for (size_t i = 0; i < Bin.size(); ++i) {
      LoBin.push_back(Bin[i][iDim].first);
   }
   return LoBin;
}


// ___________________________________________________________________________________________________
// Return vector of upper bin bounds in dim. iDim for all obs. bins
vector < double > fastNLOTable::GetObsBinsUpBounds(unsigned int iDim) const {
   if ( ! (iDim < NDim) ) {
      logger.error["GetObsBinsUpBounds"]<<"Dimension iDim " << iDim << " out of range, NDim = " << NDim << ", aborted!" << endl;
      exit(1);
   }
   // Get upper bin edge of all observable bins for dimension 'iDim'
   vector < double > UpBin;
   for (size_t i = 0; i < Bin.size(); ++i) {
      UpBin.push_back(Bin[i][iDim].second);
   }
   return UpBin;
}


// ___________________________________________________________________________________________________
// Return minimum value of all lower bin bounds for dim. iDim
double fastNLOTable::GetObsBinsLoBoundsMin(unsigned int iDim) const {
   if ( ! (iDim < NDim) ) {
      logger.error["GetObsBinsLoBoundsMin"]<<"Dimension iDim " << iDim << " out of range, NDim = " << NDim << ", aborted!" << endl;
      exit(1);
   }
   //! Get lowest bin edge of all observable bins for dimension 'iDim'
   double LoBinMin = DBL_MAX;
   for (size_t i = 0; i < Bin.size(); ++i) {
      logger.debug["GetObsBinsLoBoundsMin"]<<"iDim = " << iDim << ", i = " << i << ", Bin[i][iDim].first = " << Bin[i][iDim].first << ", LoBinMin = " << LoBinMin << endl;
      LoBinMin = ( (Bin[i][iDim].first < LoBinMin) ? Bin[i][iDim].first : LoBinMin );
   }
   logger.debug["GetObsBinsLoBoundsMin"]<<"Minimum found for dimension " << iDim << " is: " << LoBinMin << endl;
   return LoBinMin;
}


// ___________________________________________________________________________________________________
// Return maximum value of all upper bin bounds for dim. iDim
double fastNLOTable::GetObsBinsUpBoundsMax(unsigned int iDim) const {
   if ( ! (iDim < NDim) ) {
      logger.error["GetObsBinsUpBoundsMax"]<<"Dimension iDim " << iDim << " out of range, NDim = " << NDim << ", aborted!" << endl;
      exit(1);
   }
   //! Get uppermost bin edge of all observable bins for dimension 'iDim'
   double UpBinMax = -DBL_MAX;
   for (size_t i = 0; i < Bin.size(); ++i) {
      logger.debug["GetObsBinsUpBoundsMax"]<<"iDim = " << iDim << ", i = " << i << ", Bin[i][iDim].second = " << Bin[i][iDim].second << ", UpBinMax = " << UpBinMax << endl;
      UpBinMax = ( (Bin[i][iDim].second > UpBinMax) ? Bin[i][iDim].second : UpBinMax );
   }
   logger.debug["GetObsBinsUpBoundsMax"]<<"Maximum found for dimension " << iDim << " is: " << UpBinMax << endl;
   return UpBinMax;
}


// ___________________________________________________________________________________________________
// Return vector of pairs with lower and upper bin bounds in dim. iDim for all obs. bins
vector < pair < double, double > > fastNLOTable::GetObsBinsBounds(unsigned int iDim) const {
   if ( ! (iDim < NDim) ) {
      logger.error["GetObsBinsBounds"]<<"Dimension iDim " << iDim << " out of range, NDim = " << NDim << ", aborted!" << endl;
      exit(1);
   }
   // Get bin edges of all observable bins for dimension 'iDim'
   vector < pair < double, double > > Bounds;
   for (size_t i = 0; i < Bin.size(); ++i) {
      Bounds.push_back( Bin[i][iDim] );
   }
   return Bounds;
}


// ___________________________________________________________________________________________________
// Return vector of pairs with unique bin bounds of 1st dim.
vector < pair < double, double > > fastNLOTable::GetDim0BinBounds() const {
   vector< pair<double, double > > Bins = GetObsBinsBounds(0);
   set< pair< double,double> > set (Bins.begin(), Bins.end());
   Bins.assign(set.begin(),set.end());
   return Bins;
}


// ___________________________________________________________________________________________________
// Return vector of pairs with unique bin bounds of 2nd dim. for 'iDim0Bin' of 1st dim.
vector < pair < double, double > > fastNLOTable::GetDim1BinBounds(unsigned int iDim0Bin) const {
   vector< pair<double, double > > Bins;
   if ( NDim < 2 ) {
      logger.error["fastNLOTable::GetDim1BinBounds"] << "No second dimension available, aborted!" << endl;
      exit(1);
   }
   pair< double, double> bin0 = GetDim0BinBounds()[iDim0Bin];
   for (size_t iobs = 0; iobs < Bin.size(); iobs++) {
      if (Bin[iobs][0] == bin0) {
         Bins.push_back(Bin[iobs][1]);
      }
   }
   set< pair< double,double> > set (Bins.begin(), Bins.end());
   Bins.assign(set.begin(),set.end());
   return Bins;
}


// ___________________________________________________________________________________________________
// Return vector of pairs with unique bin bounds of 3rd dim. for 'iDim0Bin' and 'iDim1Bin' of 1st two dim.
vector < pair < double, double > > fastNLOTable::GetDim2BinBounds(unsigned int iDim0Bin, unsigned int iDim1Bin) const {
   vector< pair<double, double > > Bins;
   if ( NDim < 3 ) {
      logger.error["fastNLOTable::GetDim2BinBounds"] << "No third dimension available, aborted!" << endl;
      exit(1);
   }
   pair< double, double> bin0 = GetDim0BinBounds()[iDim0Bin];
   pair< double, double> bin1 = GetDim1BinBounds(iDim0Bin)[iDim1Bin];
   for (size_t iobs = 0; iobs < Bin.size(); iobs++) {
      if (Bin[iobs][0] == bin0 && Bin[iobs][1] == bin1) {
         Bins.push_back(Bin[iobs][2]);
      }
   }
   set< pair< double,double> > set (Bins.begin(), Bins.end());
   Bins.assign(set.begin(),set.end());
   return Bins;
}


// ___________________________________________________________________________________________________
// Return vector of pairs with lower and upper bin bounds for all dimensions for a given obs. bin
vector < pair < double, double > > fastNLOTable::GetObsBinDimBounds(unsigned int iObs) const {
   if ( ! ( iObs < NObsBin ) ) {
      logger.error["GetObsBinDimBounds"]<<"Observable bin iObs " << iObs << " out of range, NObsBin = " << NObsBin << ", aborted!" << endl;
      exit(1);
   }
   // // Get bin edges of all dimensions for observable bin 'iObs'
   // vector < pair < double, double > > Bounds;
   // for (size_t i = 0; i < Bin[iObs].size(); ++i) {
   //    Bounds.push_back( Bin[iObs][i] );
   // }
   // return Bounds;
   return Bin[iObs];
}


// ___________________________________________________________________________________________________
// Return pair with lower and upper bin bounds for given obs. bin and dim. iDim
pair < double, double > fastNLOTable::GetObsBinDimBounds(unsigned int iObs, unsigned int iDim) const {
   if ( ! ( iObs < NObsBin ) ) {
      logger.error["GetObsBinDimBounds"]<<"Observable bin iObs " << iObs << " out of range, NObsBin = " << NObsBin << ", aborted!" << endl;
      exit(1);
   }
   if ( ! ( iDim < NDim ) ) {
      logger.error["GetObsBinDimBounds"]<<"Dimension iDim " << iDim << " out of range, NDim = " << NDim << ", aborted!" << endl;
      exit(1);
   }
   return Bin[iObs][iDim];
}


// ___________________________________________________________________________________________________
// Return bin no. in 1st dim. for obs. bin iObs
unsigned int fastNLOTable::GetIDim0Bin(unsigned int iObs) const {
   //! Returns bin number in first dimension
   //! Valid for up to triple differential binnings
   //  There always must be at least one bin!
   if ( Bin.size() == 0 || Bin[0].size() == 0 ) {
      logger.error["fastNLOTable::GetIDim0Bin"] << "No observable bins defined, aborted!" << endl;
      exit(1);
   }
   if ( ! (iObs < NObsBin) ) {
      logger.error["GetIDim0Bin"]<<"Observable bin iObs " << iObs << " out of range, NObsBin = " << NObsBin << ", aborted!" << endl;
      exit(1);
   }
   unsigned int i0bin = 0;
   double lo0bin = Bin[0][0].first;
   for ( unsigned int i = 0; i<Bin.size(); i++ ) {
      if ( lo0bin < Bin[i][0].first ) {
         lo0bin = Bin[i][0].first;
         i0bin++;
      }
      if ( i == iObs ) {
         return i0bin;
      }
   }
   logger.error["fastNLOTable::GetIDim0Bin"] << "Observable bin not found. This should never happen, aborted!" << endl;
   exit(1);
}


// ___________________________________________________________________________________________________
// Return bin no. in 2nd dim. for obs. bin iObs
unsigned int fastNLOTable::GetIDim1Bin(unsigned int iObs) const {
   //! Returns bin number in second dimension
   //! Valid for up to triple differential binnings
   //  1d binning --> logger.error exit
   if ( NDim < 2 ) {
      logger.error["fastNLOTable::GetIDim1Bin"] << "No second dimension available, aborted!" << endl;
      exit(1);
   }
   //  Otherwise there always must be at least one bin!
   if ( Bin.size() == 0 || Bin[0].size() == 0 ) {
      logger.error["fastNLOTable::GetIDim1Bin"] << "No observable bins defined, aborted!" << endl;
      exit(1);
   }
   if ( ! (iObs < NObsBin) ) {
      logger.error["GetIDim1Bin"]<<"Observable bin iObs " << iObs << " out of range, NObsBin = " << NObsBin << ", aborted!" << endl;
      exit(1);
   }
   unsigned int i0bin = 0;
   unsigned int i1bin = 0;
   double lo0bin = Bin[0][0].first;
   double lo1bin = Bin[0][1].first;
   for ( unsigned int i = 0; i<Bin.size(); i++ ) {
      if ( lo0bin < Bin[i][0].first ) {
         lo0bin = Bin[i][0].first;
         lo1bin = Bin[i][1].first;
         i0bin++;
         i1bin = 0;
      } else if ( lo1bin < Bin[i][1].first ) {
         lo1bin = Bin[i][1].first;
         i1bin++;
      }
      if ( i == iObs ) {
         return i1bin;
      }
   }
   logger.error["fastNLOTable::GetIDim1Bin"] << "Observable bin not found. This should never happen, aborted!" << endl;
   exit(1);
}


// ___________________________________________________________________________________________________
// Return bin no. in 3rd dim. for obs. bin iObs
unsigned int fastNLOTable::GetIDim2Bin(unsigned int iObs) const {
   //! Returns bin number in third dimension
   //! Valid for up to triple differential binnings
   //  1d, 2d binning --> logger.error exit
   if ( NDim < 3 ) {
      logger.error["fastNLOTable::GetIDim2Bin"] << "No third dimension available, aborted!" << endl;
      exit(1);
   }
   //  Otherwise there always must be at least one bin!
   if ( Bin.size() == 0 || Bin[0].size() == 0 ) {
      logger.error["fastNLOTable::GetIDim2Bin"] << "No observable bins defined, aborted!" << endl;
      exit(1);
   }
   if ( iObs >= NObsBin ) {
      logger.error["fastNLOTable::GetIDim2Bin"] << "Observable bin out of range, aborted!" << endl;
      exit(1);
   }
   unsigned int i0bin = 0;
   unsigned int i1bin = 0;
   unsigned int i2bin = 0;
   double lo0bin = Bin[0][0].first;
   double lo1bin = Bin[0][1].first;
   double lo2bin = Bin[0][2].first;
   for ( unsigned int i = 0; i<Bin.size(); i++ ) {
      if ( lo0bin < Bin[i][0].first ) {
         lo0bin = Bin[i][0].first;
         lo1bin = Bin[i][1].first;
         lo2bin = Bin[i][2].first;
         i0bin++;
         i1bin = 0;
         i2bin = 0;
      } else if ( lo1bin < Bin[i][1].first ) {
         lo1bin = Bin[i][1].first;
         lo2bin = Bin[i][2].first;
         i1bin++;
         i2bin = 0;
      } else if ( lo2bin < Bin[i][2].first ) {
         lo2bin = Bin[i][2].first;
         i2bin++;
      }
      if ( i == iObs ) {
         return i2bin;
      }
   }
   logger.error["fastNLOTable::GetIDim2Bin"] << "Observable bin not found. This should never happen, aborted!" << endl;
   exit(1);
}


// ___________________________________________________________________________________________________
// Return no. of bins in 1st dimension
unsigned int fastNLOTable::GetNDim0Bins() const {
   //! Returns number of bins in first dimension
   //! Valid for up to triple differential binnings
   //  There always must be at least one bin!
   return (GetIDim0Bin(NObsBin-1) + 1);
}


// ___________________________________________________________________________________________________
// Return no. of bins in 2nd dimension for given bin in 1st dim.
unsigned int fastNLOTable::GetNDim1Bins(unsigned int iDim0Bin) const {
   //! Returns number of bins in second dimension for iDim0Bin in first dimension
   //! Valid for up to triple differential binnings
   //  1d binning --> logger.error exit
   if ( NDim < 2 ) {
      logger.error["fastNLOTable::GetNDim1Bins"] << "No second dimension available, aborted!" << endl;
      exit(1);
   }
   for ( unsigned int i = 0; i<Bin.size(); i++ ) {
      if ( GetIDim0Bin(i) == iDim0Bin +1) {
         return GetIDim1Bin(i-1) +1;
      } else if ( i == Bin.size() -1 ) {
         return GetIDim1Bin(i) +1;
      }
   }
   logger.error["fastNLOTable::GetNDim1Bins"] << "Observable bin not found. This should never happen, aborted!" << endl;
   exit(1);
}


// ___________________________________________________________________________________________________
// Return no. of bins in 3rd dimension for given bins in 1st and 2nd dim.
unsigned int fastNLOTable::GetNDim2Bins(unsigned int iDim0Bin, unsigned int iDim1Bin) const {
   //! Returns number of bins in third dimension for iDim0Bin in first and iDim1Bin in second dimension
   //! Valid for up to triple differential binnings
   //  1d, 2d binning --> logger.error exit
   if ( NDim < 3 ) {
      logger.error["fastNLOTable::GetNDim2Bins"] << "No third dimension available, aborted!" << endl;
      exit(1);
   }
   for ( unsigned int i = 0; i<Bin.size(); i++ ) {
      if ( GetIDim0Bin(i) == iDim0Bin && GetIDim1Bin(i) == iDim1Bin +1) {
         return GetIDim2Bin(i-1) +1;
      } else if ( GetIDim0Bin(i) == iDim0Bin +1 && GetIDim1Bin(i-1) == iDim1Bin) {
         return GetIDim2Bin(i-1) +1;
      } else if ( i == Bin.size() -1 ) {
         return GetIDim2Bin(i) +1;
      }
   }
   logger.error["fastNLOTable::GetNDim2Bins"] << "Observable bin not found. This should never happen, aborted!" << endl;
   exit(1);
}


// ___________________________________________________________________________________________________
// Return bin no. in 1st dim. for obs0=var0; -1 if outside range
int fastNLOTable::GetODim0Bin(double obs0) const {
   int iDim0Bin = -1;
   for ( unsigned int i=0; i<NObsBin; i++ ) {
      if ( IDiffBin[0] == 1 ) {
         logger.error["fastNLOTable::GetODim0Bin"] << "Point-wise differential not yet implemented, aborted!" << endl;
         exit(1);
      } else {
         if ( Bin[i][0].first <= obs0 && obs0 < Bin[i][0].second ) {
            iDim0Bin = GetIDim0Bin(i);
            break;
         }
      }
   }
   return iDim0Bin;
}


// ___________________________________________________________________________________________________
// Return bin no. in 2nd dim. for obs0=var0,obs1=var1; -1 if outside range
int fastNLOTable::GetODim1Bin(double obs0, double obs1) const {
   int iDim1Bin = -1;
   for ( unsigned int i=0; i<NObsBin; i++ ) {
      if ( IDiffBin[0] == 1 ) {
         logger.error["fastNLOTable::GetODim1Bin"] << "Point-wise differential not yet implemented, aborted!" << endl;
         exit(1);
      } else {
         if ( Bin[i][0].first <= obs0 && obs0 < Bin[i][0].second &&
              Bin[i][1].first <= obs1 && obs1 < Bin[i][1].second ) {
            iDim1Bin = GetIDim1Bin(i);
            break;
         }
      }
   }
   return iDim1Bin;
}


// ___________________________________________________________________________________________________
// Return bin no. in 3rd dim. for obs0=var0,obs1=var1,obs2=var2; -1 if outside range
int fastNLOTable::GetODim2Bin(double obs0, double obs1, double obs2) const {
   int iDim2Bin = -1;
   for ( unsigned int i=0; i<NObsBin; i++ ) {
      if ( IDiffBin[0] == 1 ) {
         logger.error["fastNLOTable::GetODim2Bin"] << "Point-wise differential not yet implemented, aborted!" << endl;
         exit(1);
      } else {
         if ( Bin[i][0].first <= obs0 && obs0 < Bin[i][0].second &&
              Bin[i][1].first <= obs1 && obs1 < Bin[i][1].second &&
              Bin[i][2].first <= obs2 && obs2 < Bin[i][2].second ) {
            iDim2Bin = GetIDim2Bin(i);
            break;
         }
      }
   }
   return iDim2Bin;
}


// ___________________________________________________________________________________________________
// Return observable bin no. for vector of values obs0=var0,obs1=var1,...; -1 if outside range
int fastNLOTable::GetObsBinNumber( const vector<double>& vobs ) const {
   //! Returns first matching observable bin number for vector of observations
   //! (assumes none or exactly one matching bin!)
   //! Returns -1 if outside range
   if ( vobs.size() != NDim ) {
      logger.error["GetObsBinNumber"] << "Number of observable values not equal dimensionality of the binning, aborted" << endl;
      logger.error["GetObsBinNumber"] << "NDim = " << NDim << ", vobs.size() = " << vobs.size() << endl;
      exit(1);
   }
   if ( ! (NDim < 4) ) {
      logger.error["fastNLOTable::GetObsBinNumber"] << "More than 3-dimensional binning not yet implemented, aborted!" << endl;
      exit(1);
   }

   for ( unsigned int i=0; i<NObsBin; i++ ) {
      bool lmatch = true;
      for ( unsigned int j=0; j<NDim; j++ ) {
         if ( IDiffBin[j] == 1 ) { // Point-wise differential
            lmatch = lmatch && ( fabs(Bin[i][j].first - vobs[j]) < DBL_MIN );
         } else {                  // Non- or bin-wise differential
            lmatch = lmatch && ( Bin[i][j].first <= vobs[j] && vobs[j] < Bin[i][j].second );
         }
      }
      if ( lmatch ) {return i;}
   }
   return -1;
}


// ___________________________________________________________________________________________________
// Return observable bin no. for obs0=var0 in 1D binning; -1 if outside range
int fastNLOTable::GetObsBinNumber( double obs0 ) const {
   //! Returns first matching observable bin number for one observation
   //! (assumes none or exactly one matching bin!)
   //! Returns -1 if outside range
   vector < double > vobs(1);
   vobs[0] = obs0;
   return GetObsBinNumber(vobs);
}


// ___________________________________________________________________________________________________
// Return observable bin no. for obs0=var0,obs1=var1 in 2D binning; -1 if outside range
int fastNLOTable::GetObsBinNumber( double obs0, double obs1 ) const {
   //! Returns first matching observable bin number for two observations
   //! (assumes none or exactly one matching bin!)
   //! Returns -1 if outside range
   vector < double > vobs(2);
   vobs[0] = obs0;
   vobs[1] = obs1;
   return GetObsBinNumber(vobs);
}


// ___________________________________________________________________________________________________
// Return observable bin no. for obs0=var0,obs1=var1,obs2=var2 in 3D binning; -1 if outside range
int fastNLOTable::GetObsBinNumber( double obs0, double obs1, double obs2 ) const {
   //! Returns first matching observable bin number for three observations
   //! (assumes none or exactly one matching bin!)
   //! Returns -1 if outside range
   vector < double > vobs(3);
   vobs[0] = obs0;
   vobs[1] = obs1;
   vobs[2] = obs2;
   return GetObsBinNumber(vobs);
}


// ___________________________________________________________________________________________________
// Some other info getters/setters
// ___________________________________________________________________________________________________
string fastNLOTable::GetRivetId() const {
   string identifier("RIVET_ID");
   string found;
   for (size_t i=0; i < ScDescript.size(); ++i) {
      if (ScDescript[i].find(identifier) != string::npos){
         size_t RivetIdx = ScDescript[i].find(identifier);
         size_t RivetValIdx = ScDescript[i].find("=", RivetIdx) + 1;
         size_t RivetValLen = ScDescript[i].find(",", RivetValIdx) - RivetValIdx;
         found = ScDescript[i].substr(RivetValIdx, RivetValLen);
         break;
      }
   }
   return found;
}

string fastNLOTable::GetXSDescr() const {
   string identifier("sigma");
   for (size_t i=0; i < ScDescript.size(); ++i) {
      if (ScDescript[i].find(identifier) != string::npos){
         return ScDescript[i];
      }
   }
   return "Undefined";
}

vector <string> fastNLOTable::GetScDescr() const {
   return ScDescript;
}

void fastNLOTable::SetScDescr(std::vector <std::string> ScDescr) {
   size_t NScDescript = ScDescr.size();
   fastNLOTable::ScDescript.resize(NScDescript);
   for (size_t i=0; i < NScDescript; ++i) {
      fastNLOTable::ScDescript[i] = ScDescr[i];
   }
}

void fastNLOTable::SetNObsBin(int NObs) {
   fastNLOTable::NObsBin = NObs;
}

void fastNLOTable::SetBinSize(std::vector < double > NewBinSize) {
   size_t NewSize = NewBinSize.size();
   fastNLOTable::BinSize.resize(NewSize);
   for (size_t i=0; i < NewSize; ++i) {
      fastNLOTable::BinSize[i] = NewBinSize[i];
   }
}

void fastNLOTable::SetBins(std::vector < std::vector <std::pair<double,double> > > NewBins) {
   size_t NewSize = NewBins.size();
   fastNLOTable::Bin.resize(NewSize);
   for (size_t i=0; i < NewSize; ++i) {
      fastNLOTable::Bin[i] = NewBins[i];
   }
}



// ___________________________________________________________________________________________________
// Info print out functionality
// ___________________________________________________________________________________________________


// ___________________________________________________________________________________________________
void fastNLOTable::Print(int iprint) const {
   // Define different levels of detail for printing out table content
   // The minimum (iprint = 0) just gives basic scenario information
   // including the employed code with references. The additional levels
   // are mostly for debugging purposes.
   // (Partially to be implemented!)
   //
   // iprint = 0: No additional printout
   //          1: Print Block A1 & A2 (A1, A2)
   //          2: Also print values of Block B (B0)
   //          Not implemented yet
   //          3: Also print x nodes of Block B for each contribution (BX)
   //          4: Also print scale nodes of Block B for each contribution (BS)
   //          5: Also print sigma tilde of Block B (not implemented yet)

   //
   // Print table header
   char buffer[1024];
   cout  << endl;
   cout  << fastNLO::_CSEPSC << endl;
   snprintf(buffer, sizeof(buffer), "Information on table header");
   logger.shout << buffer << endl;
   cout  << fastNLO::_SSEPSC << endl;
   PrintHeader(iprint);

   //
   // Print scenario information
   PrintScenario(iprint);

   //
   // Loop over available contributions
   for (unsigned int j = 0 ; j<fCoeff.size() ; j++) {
      fastNLOCoeffBase* c = fCoeff[j];
      char buffer[1024];
      cout  << endl;
      cout  << fastNLO::_CSEPSC << endl;
      snprintf(buffer, sizeof(buffer), "Information on table contribution no. %d: %s",j,c->CtrbDescript[0].data());
      logger.shout << buffer << endl;
      cout  << fastNLO::_SSEPSC << endl;
      // Print information for each contribution
      c->Print(iprint);
   }
}


// ___________________________________________________________________________________________________
void fastNLOTable::PrintScenario(int iprint) const {
   //
   //  Print scenario information for table
   //  iprint: iprint > 0: Print more info ...
   //
   logger.debug["PrintScenario"] << "Printing info on scenario: " << ScenName.data() << endl;
   char buffer[1024];
   cout  << endl;
   cout  << fastNLO::_CSEPSC << endl;
   snprintf(buffer, sizeof(buffer), "Information on fastNLO scenario: %s",ScenName.data());
   logger.shout << buffer << endl;
   cout  << fastNLO::_SSEPSC << endl;
   if ( !(iprint < 0) ) {
      cout << fastNLO::_DSEP20C << " fastNLO Table: Scenario " << fastNLO::_DSEP20 << endl;
   } else {
      cout << endl << fastNLO::_CSEP20C << " fastNLO Table: Scenario " << fastNLO::_CSEP20 << endl;
   }
   fastNLOTools::PrintVector(ScDescript,"Scenario description (ScDescript)","#");
   printf(" #\n");
   printf(" # Publ. x section (10^-Ipublunits b)  %d\n",Ipublunits);
   printf(" # Centre-of-mass energy (Ecms/GeV)    %5.0f\n",Ecms);
   printf(" # Power in a_s of LO process (ILOord) %d\n",ILOord);
   printf(" # No. of observable bins (NObsBin)    %d\n",NObsBin);
   printf(" # Dim. of observable binning (NDim)   %d\n",NDim);
   printf(" #\n");
   fastNLOTools::PrintVector(DimLabel,"Dimension labels (DimLabel)","#");
   fastNLOTools::PrintVector(IDiffBin,"Differential dimension (IDiffBin)","#");
   printf(" #\n");
   if ( abs(iprint) > 1 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 1) " << fastNLO::_SSEP20 << endl;
      for (unsigned int i=0; i<NObsBin; i++) {
         // Print only for first and last observable bin
         if (i==0 || i==NObsBin-1) {
            for (unsigned int j=0; j<NDim; j++) {
               printf(" #   LoBin[%d][%d]                        %7.4f\n", i,j,Bin[i][j].first);
               if ( IDiffBin[j]==2 ) {
                  printf(" #   UpBin[%d][%d]                       %7.4f\n", i,j,Bin[i][j].second);
               }
            }
         }
      }
      for (unsigned int i=0; i<NObsBin; i++) {
         // Print only for first and last observable bin
         if (i==0 || i==NObsBin-1) {
            printf(" #   BinSize[%d]                       %7.4f\n", i,BinSize[i]);
         }
      }
   }
   if( INormFlag != 0 ) {
      printf(" # Normalization flag (INormFlag)      %d\n",INormFlag);
      if ( INormFlag<0 ) {
         printf(" # Normalization table (DenomTable)    %s\n",DenomTable.data());
      }
      if ( abs(iprint) > 1 ) {
         cout << fastNLO::_SSEP20C << " Extended information (iprint > 1) " << fastNLO::_SSEP20 << endl;
         for (unsigned int i=0; i<NObsBin; i++) {
            // Print only for first and last observable bin
            if (i==0 || i==NObsBin-1) {
               printf(" #  IDivLoPointer[%d]               %d\n",i,IDivLoPointer[i]);
               printf(" #  IDivUpPointer[%d]               %d\n",i,IDivUpPointer[i]);
            }
         }
      }
      printf(" #\n");
   }
   printf(" # Total no. of contributions (theory + optional data) in this table: %d\n",(int)fCoeff.size());
   cout << fastNLO::_CSEPSC << endl;
}


//______________________________________________________________________________
void fastNLOTable::PrintContributionSummary(int iprint) const {
   //
   //  Print summary of contributions available in table
   //  iprint: iprint > 0: Print full descriptions (all lines) for each contribution
   //
   logger.debug["PrintContributionSummary"] << "Printing flag iprint = " << iprint << endl;
   char buffer[1024];
   cout  << endl;
   cout  << fastNLO::_CSEPSC << endl;
   logger.shout << "Overview on contribution types and numbers contained in table: " << ffilename << endl;
   cout  << fastNLO::_SSEPSC << endl;
   snprintf(buffer, sizeof(buffer), "Total number of contributions: %2i", (int)fCoeff.size());
   logger.shout << buffer << endl;

   int iccount[21] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   int ictype = 0;
   string coeffname;
   for (unsigned int j = 0 ; j<fCoeff.size() ; j++) {
      fastNLOCoeffBase* c = fCoeff[j];
      if ( fastNLOCoeffData::CheckCoeffConstants(c,true) ) {
         ictype = 20;
         coeffname = "Data";
      } else {
         ictype = fCoeff[j]->GetIContrFlag1()-1;
         coeffname = fastNLO::_ContrName[ictype];
      }
      iccount[ictype]++;
      snprintf(buffer, sizeof(buffer), "  No.: %d, type: %-30.30s, Id: %d, order: %-20.20s, by: %s",
               j+1,
               coeffname.c_str(),
               iccount[ictype]-1,
               c->GetContributionDescription()[0].c_str(),
               c->GetCodeDescription()[0].c_str());
      logger.shout << buffer << endl;

      if (iprint > 0) {
         for (unsigned int k = 0 ; k<c->GetCodeDescription().size(); k++) {
            snprintf(buffer, sizeof(buffer), "          %s",c->GetCodeDescription()[k].c_str());
            logger.shout << buffer << endl;
         }
      }
   }
   int nc = 0;
   for ( int icnt=0; icnt<21; icnt++ ) {
      nc += iccount[icnt];
   }
   if ( iccount[2]+iccount[3] != GetNmult() ) {
      logger.warn["PrintContributionSummary"] << "Multiplicative contribution not correctly advertised in table header." << endl;
      logger.warn["PrintContributionSummary"] << "Nmult = " << GetNmult() <<
         " should equal " << iccount[2]+iccount[3] << " instead. Continue anyway, since not actually used." << endl;
   }
   if ( iccount[20] > 1 ) {
      logger.error["PrintContributionSummary"] << "Maximally one data contribution allowed per table," << endl;
      logger.error["PrintContributionSummary"] << "but found " << iccount[20] << "! Aborted!" << endl;
      exit(1);
   }
   if ( iccount[20] != GetNdata() ) {
      if ( nc == GetNcontrib()+GetNdata() ) {
         logger.warn["PrintContributionSummary"] << "Data contribution not correctly advertised in table header." << endl;
         logger.warn["PrintContributionSummary"] << "Ncontrib = " << GetNcontrib() << " and Ndata = " << GetNdata() <<
            " should equal " << nc-1 << " and " << iccount[20] << " instead. Continue anyway, since only sum is actually used." << endl;
      } else {
         logger.error["PrintContributionSummary"] << "Inconsistent number of contributions found!" << endl;
         logger.error["PrintContributionSummary"] << "Ncontrib = " << GetNcontrib() << " and Ndata = " << GetNdata() <<
            " should be " << nc-1 << " and " << iccount[20] << " instead. Aborted!" << endl;
         exit(1);
      }
   }
   cout << fastNLO::_CSEPSC << endl;
}


//_DEPRECATED___________________________________________________________________
void fastNLOTable::PrintFastNLOTableConstants(const int iprint) const {
   logger.error["PrintFastNLOTableConstants"]<<"This function is deprecated, aborted!"<<endl;
   logger.error["PrintFastNLOTableConstants"]<<"Please use Print instead."<<endl;
}


//_DEPRECATED___________________________________________________________________
void fastNLOTable::PrintTableInfo(const int iprint) const {
   logger.error["PrintTableInfo"]<<"This function is deprecated, aborted!"<<endl;
   logger.error["PrintTableInfo"]<<"Please use PrintContributionSummary instead."<<endl;
}


// DO NOT USE ANYTHING BELOW! DOES NOT WORK YET!
// ___________________________________________________________________________________________________
// DO NOT USE! DOES NOT WORK YET!
// unsigned int fastNLOTable::GetIDimBin(unsigned int iObsBin, unsigned int iDim) const {
//    //! Returns bin number in dimension iDim
//    //  iDim larger than table dimensions --> logger.error exit
//    logger.error["fastNLOTable::GetIDimBin"] << "DO NOT USE! DOES NOT WORK YET!" << endl;
//    const unsigned int idiff = GetNumDiffBin();
//    if ( ! (iDim < idiff) ) {
//       logger.error["fastNLOTable::GetIDimBin"] << "Requested dimension iDim not available, aborted!" << endl;
//       exit(1);
//    }
//    //  Otherwise there always must be at least one bin!
//    if ( Bin.size() == 0 || Bin[0].size() == 0 ) {
//       logger.error["fastNLOTable::GetIDimBin"] << "No observable bins defined, aborted!" << endl;
//       exit(1);
//    }
//    if ( iObsBin >= NObsBin ) {
//       logger.error["fastNLOTable::GetIDimBin"] << "Observable bin out of range, aborted!" << endl;
//       exit(1);
//    }
//    vector < unsigned int > ibin(NDim);
//    vector < double > lobin(NDim);
//    for ( unsigned int jdim = 0; jdim<=iDim; jdim++ ) {
//       ibin.push_back(0);
//       lobin.push_back(Bin[0][jdim].first);
//    }
//    for ( unsigned int iobs = 0; iobs<Bin.size(); iobs++ ) {
//       if (iDim == 2) {
//          if ( lobin[0] < Bin[iobs][0].first ) {
//             for (unsigned int k=0; k<=iDim; k++) {
//                lobin[k] = Bin[iobs][k].first;
//             }
//             ibin[0]++;
//             for (unsigned int k=1; k<=iDim; k++) {
//                ibin[k] = 0;
//             }
//          } else if ( lobin[1] < Bin[iobs][1].first ) {
//             for (unsigned int k=1; k<=iDim; k++) {
//                lobin[k] = Bin[iobs][k].first;
//             }
//             ibin[1]++;
//             for (unsigned int k=2; k<=iDim; k++) {
//                ibin[k] = 0;
//             }
//          } else if ( lobin[2] < Bin[iobs][2].first ) {
//             for (unsigned int k=2; k<=iDim; k++) {
//                lobin[k] = Bin[iobs][k].first;
//             }
//             ibin[2]++;
//          }
//       } else if (iDim == 1) {
//          if ( lobin[0] < Bin[iobs][0].first ) {
//             for (unsigned int k=0; k<=iDim; k++) {
//                lobin[k] = Bin[iobs][k].first;
//             }
//             ibin[0]++;
//             for (unsigned int k=1; k<=iDim; k++) {
//                ibin[k] = 0;
//             }
//          } else if ( lobin[1] < Bin[iobs][1].first ) {
//             for (unsigned int k=1; k<=iDim; k++) {
//                lobin[k] = Bin[iobs][k].first;
//             }
//             ibin[1]++;
//          }
//       } else if (iDim == 0 ) {
//          if ( lobin[0] < Bin[iobs][0].first ) {
//             for (unsigned int k=0; k<=iDim; k++) {
//                lobin[k] = Bin[iobs][k].first;
//             }
//             ibin[0]++;
//          }
//       }
//       if ( iobs == iObsBin ) {
//          return ibin[iDim];
//       }
//    }
//    logger.error["fastNLOTable::GetIDimBin"] << "Observable bin not found. This should never happen, aborted!" << endl;
//    exit(1);
// }
//
//
// ___________________________________________________________________________________________________
// DO NOT USE! DOES NOT WORK YET!
// vector < pair < double, double > > fastNLOTable::GetBinBoundaries(int iDim0Bin, int iDim1Bin, int iDim2Bin) {
//    //!
//    //! Get bin boundaries for first, second, and third dimension
//    //!    Assuming for instance following 2-dimensional binning scheme:
//    //!
//    //!    iDim0Bin  ________________________________
//    //!       0      |___|___|___|_______|__|__|____|
//    //!  D    1      |____|____|____|_____|____|____|
//    //!  I    2      |__|__|___|__|__|___|__|___|___|
//    //!  M    3      |______|_______|_________|_____|
//    //!       4      |__________|_______|_________|_|
//    //!  0    5      |______________|______|___|____|
//    //!       6      |_______|_______|______|_______|
//    //!                            DIM 1
//    //!  iDim1Bin may be different for each iDim0Bin
//    //!
//    //! usage e.g.:
//    //! int LowerBoundary = GetBinBoundaries(ibin)[dim].first;
//    //! int UpperBoundary = GetBinBoundaries(ibin)[dim].second;
//    //! 'dim' must be smaller than number of parameters passed to GetBinBoundaries
//    //! usage e.g.:
//    //! int LoYBin  = GetBinBoundaries(2)[0].first;
//    //! int UpPtBin = GetBinBoundaries(2,5)[1].second;
//    //! int LoYBin  = GetBinBoundaries(2,5)[0].second;
//    logger.error["GetBinBoundaries"] << "DO NOT USE! DOES NOT WORK YET!" << endl;
//    vector<pair<double,double> > BinRet(NDim);
//    const int idiff = GetNumDiffBin();
//
//    if ( idiff==1 ) {
//       if ( iDim0Bin<0 || iDim0Bin >= (int)NObsBin ) {
//          logger.warn["GetBinBoundaries"]<<"0th dimension does only have "<<NObsBin<<" but bin "<<iDim0Bin<<" was requested."<<endl;
//          return BinRet;
//       }
//       return Bin[iDim0Bin];
//    }
//    else if ( idiff==2) {
//       unsigned int nDim1 = GetNDim1Bins(iDim0Bin);
//       unsigned int nDim0 = GetNDim0Bins();
//       // sanity
//       if ( iDim0Bin < 0 || iDim0Bin >= (int)nDim0 ) {
//          logger.warn["GetBinBoundaries"]<<"Dimension 2 does only have "<<nDim0<<" but bin "<<iDim0Bin<<" was requested."<<endl;
//       }
//       if ( iDim1Bin < 0 || iDim1Bin >= (int)nDim1 ) {
//          logger.warn["GetBinBoundaries"]<<"Dimension 1 does only have "<<nDim1<<" but bin "<<iDim1Bin<<" was requested."<<endl;
//       }
//       int iObs = 0;
//       for ( int i0 = 0 ; i0<iDim0Bin ; i0++ ) {
//          iObs += GetNDim1Bins(i0);
//       }
//       iObs += iDim1Bin;
//       return Bin[iObs];
//    }
//    else if ( idiff == 3 ) {
//       unsigned int nDim0 = GetNDim0Bins();
//       unsigned int nDim1 = GetNDim1Bins(iDim0Bin);
//       unsigned int nDim2 = GetNDim2Bins(iDim0Bin,iDim1Bin);
//       // sanity
//       if ( iDim0Bin < 0 || iDim0Bin >= (int)nDim0 ) {
//          logger.warn["GetBinBoundaries"]<<"Dimension 0 does only have "<<nDim0<<" but bin "<<iDim0Bin<<" was requested."<<endl;
//       }
//       if ( iDim1Bin < 0 || iDim1Bin >= (int)nDim1 ) {
//          logger.warn["GetBinBoundaries"]<<"Dimension 1 does only have "<<nDim1<<" but bin "<<iDim1Bin<<" was requested."<<endl;
//       }
//       if ( iDim2Bin < 0 || iDim2Bin >= (int)nDim2 ) {
//          logger.warn["GetBinBoundaries"]<<"Dimension 2 does only have "<<nDim2<<" but bin "<<iDim2Bin<<" was requested."<<endl;
//       }
//       logger.error["GetBinBoundaries"]<<"todo. Further code no yet implemented."<<endl;
//
//    }
//    else {
//       logger.error["GetBinBoundaries"]<<"Higher than triple-differential binnings are not implemented."<<endl;
//    }
//    return BinRet;
//
// }


// Erase observable bin; iObsIdx is the C++ array index to be removed and
// not the observable bin no. running from 1 to NObsBin
template<typename T> void fastNLOTable::EraseBin(vector<T>& v, unsigned int idx) {
   if ( v.empty() ) {
      logger.warn["fastNLOTable::EraseBin"]<<"Empty vector, nothing to erase!" << endl;
   } else if ( idx < v.size() ) {
      logger.info["fastNLOTable::EraseBin"]<<"Erasing vector index no. " << idx << endl;
      v.erase(v.begin()+idx);
   } else {
      logger.error["fastNLOTable::EraseBin"]<<"Bin no. larger than vector size, aborted!" << endl;
      exit(1);
   }
}

void fastNLOTable::EraseBinFromTable(unsigned int iObsIdx) {
   logger.info["fastNLOTable::EraseBinFromTable"]<<"Erasing from table the observable index no. " << iObsIdx << endl;
   // Changes to table header block A2
   EraseBin(fastNLOTable::Bin,iObsIdx);
   EraseBin(fastNLOTable::BinSize,iObsIdx);
   if ( fastNLOTable::INormFlag != 0 ) {
      EraseBin(fastNLOTable::IDivLoPointer,iObsIdx);
      EraseBin(fastNLOTable::IDivUpPointer,iObsIdx);
   }
   // Changes to table contributions block B
   cout << "Ncontrib = " << GetNcontrib() << ", Ndata = " << GetNdata() << endl;
   for ( int ic = 0; ic<GetNcontrib()+GetNdata(); ic++ ) {
      logger.info["fastNLOTable::EraseBinFromTable"]<<"Erasing the observable index no. " << iObsIdx << " from contribution no. " << ic << endl;
      fastNLOCoeffAddBase* ctmp = (fastNLOCoeffAddBase*)fCoeff[ic];

      // Identify type of coeff-table
      bool quiet = true;
      if ( fastNLOCoeffData::CheckCoeffConstants(ctmp,quiet) ) {
         logger.info["EraseBinFromTable"]<<"Found data contribution. Now erasing index no. " << iObsIdx << endl;
         fastNLOCoeffData* cdata = (fastNLOCoeffData*)fCoeff[ic];
         cdata->EraseBin(iObsIdx);
      } else if ( fastNLOCoeffMult::CheckCoeffConstants(ctmp,quiet) ) {
         logger.info["EraseBinFromTable"]<<"Found multiplicative contribution. Now erasing index no. " << iObsIdx << endl;
         fastNLOCoeffMult* cmult = (fastNLOCoeffMult*)fCoeff[ic];
         cmult->EraseBin(iObsIdx);
      } else if ( fastNLOCoeffAddFix::CheckCoeffConstants(ctmp,quiet) ) {
         logger.info["EraseBinFromTable"]<<"Found additive fix-table contribution. Now erasing index no. " << iObsIdx << endl;
         fastNLOCoeffAddFix* cfix = (fastNLOCoeffAddFix*)fCoeff[ic];
         cfix->EraseBin(iObsIdx);
      } else if ( fastNLOCoeffAddFlex::CheckCoeffConstants(ctmp,quiet) ) {
         logger.info["EraseBinFromTable"]<<"Found additive flex-table contribution. Now erasing index no. " << iObsIdx << endl;
         fastNLOCoeffAddFlex* cflex = (fastNLOCoeffAddFlex*)fCoeff[ic];
         cflex->EraseBin(iObsIdx);
      } else {
         logger.error["EraseBinFromTable"]<<"Could not identify contribution. Print and abort!" << endl;
         ctmp->Print(-1);
         exit(1);
      }
   }
   // Reduce no. of observable bins
   SetNObsBin(GetNObsBin()-1);
}



// Multiply observable bin; iObsIdx is the C++ array index to be multiplied and
// not the observable bin no. running from 1 to NObsBin
template<typename T> void fastNLOTable::MultiplyBin(vector<T>& v, unsigned int idx, double fact) {
   if ( v.empty() ) {
      logger.warn["fastNLOTable::MultiplyBin"]<<"Empty vector, nothing to multiply!" << endl;
   } else if ( idx < v.size() ) {
      logger.info["fastNLOTable::MultiplyBin"]<<"Multiplying vector index no. " << idx << endl;
      v[idx] *= fact;
   } else {
      logger.error["fastNLOTable::MultiplyBin"]<<"Bin no. larger than vector size, aborted!" << endl;
      exit(1);
   }
}

void fastNLOTable::MultiplyBinSize(unsigned int iObsIdx, double fact) {
   logger.debug["fastNLOTable::MultiplyBinSize"]<<"Multiplying the bin size of the observable index no. " << iObsIdx << " by " << fact << endl;
   MultiplyBin(fastNLOTable::BinSize,iObsIdx,fact);
}

void fastNLOTable::MultiplyBinInTable(unsigned int iObsIdx, double fact) {
   logger.debug["fastNLOTable::MultiplyBinInTable"]<<"Multiplying the observable index no. " << iObsIdx << endl;
   // Changes to table header block A2
   // Changes to table contributions block B
   for ( int ic = 0; ic<GetNcontrib()+GetNdata(); ic++ ) {
      logger.debug["fastNLOTable::MultiplyBinInTable"]<<"Multiplying the observable index no. " << iObsIdx << " from contribution no. " << ic << endl;
      fastNLOCoeffAddBase* ctmp = (fastNLOCoeffAddBase*)fCoeff[ic];

      // Identify type of coeff-table
      bool quiet = true;
      if ( fastNLOCoeffData::CheckCoeffConstants(ctmp,quiet) ) {
         logger.debug["MultiplyBinInTable"]<<"Found data contribution. Skipped! Index no. " << iObsIdx << endl;
         fastNLOCoeffData* cdata = (fastNLOCoeffData*)fCoeff[ic];
         cdata->MultiplyBin(iObsIdx,fact);
      } else if ( fastNLOCoeffMult::CheckCoeffConstants(ctmp,quiet) ) {
         logger.debug["MultiplyBinInTable"]<<"Found multiplicative contribution. Skipped! Index no. " << iObsIdx << endl;
         fastNLOCoeffMult* cmult = (fastNLOCoeffMult*)fCoeff[ic];
         cmult->MultiplyBin(iObsIdx,fact);
      } else if ( fastNLOCoeffAddFix::CheckCoeffConstants(ctmp,quiet) ) {
         logger.debug["MultiplyBinInTable"]<<"Found additive fix-table contribution. Now multiplying index no. " << iObsIdx << endl;
         fastNLOCoeffAddFix* cfix = (fastNLOCoeffAddFix*)fCoeff[ic];
         cfix->MultiplyBin(iObsIdx,fact);
      } else if ( fastNLOCoeffAddFlex::CheckCoeffConstants(ctmp,quiet) ) {
         logger.debug["MultiplyBinInTable"]<<"Found additive flex-table contribution. Now multiplying index no. " << iObsIdx << endl;
         fastNLOCoeffAddFlex* cflex = (fastNLOCoeffAddFlex*)fCoeff[ic];
         cflex->MultiplyBin(iObsIdx,fact);
      } else {
         logger.error["MultiplyBinInTable"]<<"Could not identify contribution. Print and abort!" << endl;
         ctmp->Print(-1);
         exit(1);
      }
   }
}

void fastNLOTable::CatBinToTable(const fastNLOTable& other, unsigned int iObsIdx, unsigned int table_count) {
   logger.info["fastNLOTable::CatBinToTable"]<<"Catenating the observable bin index no. " << iObsIdx << " from other table to this." << endl;
   // Changes to table header block A2
   CatBin(other,iObsIdx,table_count);
   // Changes to table contributions block B
   // Loop over all contributions from 'other'-table
   for ( int ic=0; ic<other.GetNcontrib()+other.GetNdata(); ic++ ) {
      logger.info["fastNLOTable::CatBinToTable"]<<"Catenating the observable index no. " << iObsIdx << " from contribution no. " << ic << endl;
      // Find matching contribution from 'this'-table
      for ( unsigned int jc=0; jc<fCoeff.size(); jc++) {
         bool quiet = true;
         //         fastNLOCoeffBase* cthis  = (fastNLOCoeffBase*)fCoeff[jc];
         fastNLOCoeffBase* cother = (fastNLOCoeffBase*)other.GetCoeffTable(ic);
         // Identify type of other coeff table
         // Additive fixed-order?
         if ( fastNLOCoeffAddBase::CheckCoeffConstants(cother,quiet) ) {
            if ( fastNLOCoeffAddFix::CheckCoeffConstants(cother,quiet) ) {
               fastNLOCoeffAddFix* clhs = (fastNLOCoeffAddFix*)fCoeff[jc];
               fastNLOCoeffAddFix* crhs = (fastNLOCoeffAddFix*)other.GetCoeffTable(ic);
               if ( clhs->IsCatenable(*crhs) ) {
                  logger.info["CatBinToTable"]<<"Found fix-scale additive contribution. Now catenating index no. " << iObsIdx << endl;
                  clhs->CatBin(*crhs,iObsIdx);
                  continue;
               }
            }
            else if ( fastNLOCoeffAddFlex::CheckCoeffConstants(cother,quiet) ) {
               fastNLOCoeffAddFlex* clhs = (fastNLOCoeffAddFlex*)fCoeff[jc];
               fastNLOCoeffAddFlex* crhs = (fastNLOCoeffAddFlex*)other.GetCoeffTable(ic);
               if ( clhs->IsCatenable(*crhs) ) {
                  logger.info["CatBinToTable"]<<"Found flex-scale additive contribution. Now catenating index no. " << iObsIdx << endl;
                  clhs->CatBin(*crhs,iObsIdx);
                  continue;
               }
            }
         }
         // Multiplicative?
         else if ( fastNLOCoeffMult::CheckCoeffConstants(cother,quiet) ) {
            fastNLOCoeffMult* clhs = (fastNLOCoeffMult*)fCoeff[jc];
            fastNLOCoeffMult* crhs = (fastNLOCoeffMult*)other.GetCoeffTable(ic);
            if ( clhs->IsCatenable(*crhs) ) {
               logger.info["CatBinToTable"]<<"Found multiplicative contribution. Now catenating index no. " << iObsIdx << endl;
               clhs->CatBin(*crhs,iObsIdx);
               continue;
            }
         }
         // Data?
         else if ( fastNLOCoeffData::CheckCoeffConstants(cother,quiet) ) {
            fastNLOCoeffData* clhs = (fastNLOCoeffData*)fCoeff[jc];
            fastNLOCoeffData* crhs = (fastNLOCoeffData*)other.GetCoeffTable(ic);
            if ( clhs->IsCatenable(*crhs) ) {
               logger.info["CatBinToTable"]<<"Found data contribution. Now catenating index no. " << iObsIdx << endl;
               clhs->CatBin(*crhs,iObsIdx);
               continue;
            }
         }
         // Unknown
         else {
            logger.error["CatBinToTable"]<<"Could not identify contribution. Print and abort!" << endl;
            cother->Print(-1);
            exit(1);
         }
      }
   }
   // Increase no. of observable bins
   SetNObsBin(GetNObsBin()+1);
}

// Catenate observable bin
void fastNLOTable::CatBin(const fastNLOTable& other, unsigned int iObsIdx, unsigned int table_count) {
   logger.debug["fastNLOTable::CatBin"]<<"Catenating observable bin in scenario header corresponding to bin index " << iObsIdx << endl;
   if ( Bin.size() == 0 ) {
      say::error["CatBin"]<<"Bin size cannot be zero for a fastNLO table. Aborted!" << endl;
      exit(1);
   }
   static unsigned int noff = 0;
   static unsigned int ntab = 0;
   unsigned int nold = Bin.size();
   if ( ntab != table_count ) {
      ntab = table_count;
      noff = nold;
   }
   Bin.resize(nold+1);
   Bin[nold] = other.Bin[iObsIdx];
   BinSize.resize(nold+1);
   BinSize[nold] = other.BinSize[iObsIdx];
   if ( fastNLOTable::INormFlag != 0 ) {
      IDivLoPointer.resize(nold+1);
      IDivUpPointer.resize(nold+1);
      if ( fastNLOTable::INormFlag == 2 ) {
         IDivLoPointer[nold] = noff + other.IDivLoPointer[iObsIdx];
         IDivUpPointer[nold] = noff + other.IDivUpPointer[iObsIdx];
      } else {
         say::error["CatBin"]<<"Table catenation not yet implemented for INormFlag = " << fastNLOTable::INormFlag << ". Aborted!" << endl;
         exit(1);
      }
   }
}

//
// functions previously included in fastNLOBase
//

//______________________________________________________________________________
zstr::ifstream* fastNLOTable::OpenFileRead() {
   //! Open file-stream for reading table
   // does file exist?
   if (access(ffilename.c_str(), R_OK) != 0) {
      logger.error["OpenFileRead"]<<"File does not exist! Was looking for: "<<ffilename<<". Exiting."<<endl;
      exit(1);
   }
   zstr::ifstream* strm = new zstr::ifstream(ffilename.c_str(),ios::in);
   return strm;
}


//______________________________________________________________________________
void fastNLOTable::CloseFileRead(zstr::ifstream& strm) {
   //! Close file-stream
   // strm.close();
   delete &strm;
}


//______________________________________________________________________________
std::ostream* fastNLOTable::OpenFileWrite(bool compress) {
   //! open ostream for writing tables
   //! do overwrite existing table
   if (access(ffilename.c_str(), F_OK) == 0) {
      logger.info["OpenFileWrite"]<<"Overwriting the already existing table file: " << ffilename << endl;
   }
   std::ostream* stream = NULL;
   if (compress) {
      stream = new zstr::ofstream(ffilename.c_str(),ios::out);
   } else {
      stream = new std::ofstream(ffilename.c_str(),ios::out);
   }
   if (!stream->good()) {
      logger.error["OpenFileWrite"]<<"Cannot open file '"<<ffilename<<"' for writing. Aborting."<<endl;
      exit(2);
   }
   stream->precision(fPrecision);
   return stream;
}


//______________________________________________________________________________
void fastNLOTable::CloseFileWrite(std::ostream& table) {
   //! close stream and delete object;
   table << fastNLO::tablemagicno << endl;
   table << fastNLO::tablemagicno << endl;
   // table.close();
   delete &table;
}


//______________________________________________________________________________
bool fastNLOTable::IsCompatibleHeader(const fastNLOTable& other) const {
   if (Itabversion!= other.GetItabversion()) {
      logger.warn["IsCompatibleHeader"]<<"Differing versions of table format: "<<Itabversion<<" and "<< other.GetItabversion()<<endl;
      return false;
   }
   if (GetNdata() + other.GetNdata() > 1) {
      logger.warn["IsCompatibleHeader"]<<"Two tables containing both experimental data are incompatible"<<endl;
      return false;
   }
   if (ScenName!= other.GetScenName()) {
      logger.warn["IsCompatibleHeader"]<<"Differing names of scenarios: "<<ScenName.c_str()<<" and "<<other.ScenName.c_str()<<endl;
      // continue...
   }
   return true;
}


//______________________________________________________________________________
bool fastNLOTable::IsCatenableHeader(const fastNLOTable& other) const {
   if ( trunc(Itabversion/10000) != trunc(other.GetItabversion()/10000)) {
      logger.error["IsCatenableHeader"]<<"Differing versions of table format: "<<Itabversion<<" and "<<other.GetItabversion()<<endl;
      return false;
   } else if ( Itabversion != other.GetItabversion() ) {
      logger.warn["IsCatenableHeader"]<<"Differing sub-versions of table format: "<<Itabversion<<" and "<<other.GetItabversion()<<endl;
      logger.warn["IsCatenableHeader"]<<"Please check your result carefully!"<<endl;
   }
   if (GetNcontrib() != other.GetNcontrib()) {
      logger.warn["IsCatenableHeader"]<<"Differing number of contributions: "<<GetNcontrib()<<" and "<<other.GetNcontrib()<<endl;
      return false;
   }
   if (GetNmult() != other.GetNmult()) {
      logger.warn["IsCatenableHeader"]<<"Differing number of multiplicative contributions: "<<GetNmult()<<" and "<<other.GetNmult()<<endl;
      return false;
   }
   if (GetNdata() != other.GetNdata()) {
      logger.warn["IsCatenableHeader"]<<"Differing number of data contributions: "<<GetNdata()<<" and "<<other.GetNdata()<<endl;
      return false;
   }
   return true;
}

//______________________________________________________________________________
void fastNLOTable::PrintHeader(int iprint) const {
   if ( !(iprint < 0) ) {
      cout << fastNLO::_DSEP20C << " fastNLO Table: Header " << fastNLO::_DSEP20 << endl;
   } else {
      cout << endl << fastNLO::_CSEP20C << " fastNLO Table: Header " << fastNLO::_CSEP20 << endl;
   }
   printf(" # Table version (Itabversion)         %d\n",Itabversion);
   printf(" # Scenario name (ScenName)            %s\n",ScenName.data());
   printf(" # Theory contributions (Ncontrib)     %d\n",GetNcontrib());
   printf(" # Data contribution 0/1 (Ndata)       %d\n",GetNdata());
   if ( abs(iprint) > 0 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 0) " << fastNLO::_SSEP20 << endl;
      printf(" #   Separator (tablemagicno)            %d\n",fastNLO::tablemagicno);
      printf(" #   Unused (Nmult)                      %d\n",GetNmult());
   }
   cout << fastNLO::_CSEPSC << endl;
}


//______________________________________________________________________________
void fastNLOTable::PrintWelcomeMessage() {

   char fnlo[100];
   sprintf(fnlo,"%c[%d;%dmfast%c[%d;%dmNLO\033[0m",27,0,31,27,0,34);
   char subproject[100]      = FNLO_SUBPROJECT;
   char package_version[100] = FNLO_VERSION;
   char svnrev[100]          = FNLO_SVNREV;
   char authors[500]         = FNLO_AUTHORS;
   char webpage[500]         = FNLO_WEBPAGE;
   char authorsv14[200]      = FNLO_AUTHORSv14;
   char quotev14[200]        = FNLO_QUOTEv14;
   char authorsv2[200]       = FNLO_AUTHORSv2;
   char quotev2[200]         = FNLO_QUOTEv2;
   char years[100]           = FNLO_YEARS;

   cout  << endl;
   cout  << fastNLO::_CSEPSC << endl;
   speaker &shout = logger.shout;
   cout << " #" << endl;
   shout << fnlo << "_" << subproject << endl;
   shout << "Version " << package_version << "_" << svnrev << endl;
   cout << " #" << endl;
   shout << "C++ program and toolkit to read and create fastNLO v2 tables and" << endl;
   shout << "derive QCD cross sections using PDFs, e.g. from LHAPDF" << endl;
   cout << " #" << endl;
   cout  << fastNLO::_SSEPSC << endl;
   cout << " #" << endl;
   shout << "Copyright © " << years << " " << fnlo << " Collaboration" << endl;
   shout << authors << endl;
   cout << " #" << endl;
   shout << "This program is free software: you can redistribute it and/or modify" << endl;
   shout << "it under the terms of the GNU General Public License as published by" << endl;
   shout << "the Free Software Foundation, either version 3 of the License, or" << endl;
   shout << "(at your option) any later version." << endl;
   cout << " #" << endl;
   shout << "This program is distributed in the hope that it will be useful," << endl;
   shout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
   shout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the" << endl;
   shout << "GNU General Public License for more details." << endl;
   cout << " #" << endl;
   shout << "You should have received a copy of the GNU General Public License" << endl;
   shout << "along with this program. If not, see <http://www.gnu.org/licenses/>." << endl;
   cout << " #" << endl;
   cout  << fastNLO::_SSEPSC << endl;
   cout << " #" << endl;
   shout << "The projects web page can be found at:" << endl;
   shout << "  " << webpage << endl;
   cout << " #" << endl;
   shout << "If you use this code, please cite:" << endl;
   shout << "  " << authorsv14 << ", " << quotev14 << endl;
   shout << "  " << authorsv2 << ", " << quotev2 << endl;
   cout << " #" << endl;
   cout  << fastNLO::_CSEPSC << endl;
   fWelcomeOnce = true;
}
