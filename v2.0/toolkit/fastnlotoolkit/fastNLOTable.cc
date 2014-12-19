#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <set>
#include "fastnlotk/fastNLOTable.h"
#include "fastnlotk/fastNLOTools.h"


using namespace std;

// ___________________________________________________________________________________________________
//fastNLOTable::fastNLOTable() : PrimalScream("fastNLOTable") {
fastNLOTable::fastNLOTable(){
   SetClassName("fastNLOTable");
}


// ___________________________________________________________________________________________________
fastNLOTable::fastNLOTable(string name) : fastNLOBase(name) {
   SetClassName("fastNLOTable");
   ReadTable();
}


// ___________________________________________________________________________________________________
fastNLOTable::~fastNLOTable(){
   // delete fCoeff tables...
   DeleteAllCoeffTable();
}

// ___________________________________________________________________________________________________
fastNLOTable::fastNLOTable(const fastNLOTable& tab)
   : fastNLOBase(tab), fCoeff(tab.fCoeff.size()),
    Ecms(tab.Ecms), ILOord(tab.ILOord), Ipublunits(tab.Ipublunits),
    ScDescript(tab.ScDescript), NObsBin(tab.NObsBin), NDim(tab.NDim),
    DimLabel(tab.DimLabel), IDiffBin(tab.IDiffBin), Bin(tab.Bin),
    BinSize(tab.BinSize), INormFlag(tab.INormFlag),
    DenomTable(tab.DenomTable), IDivLoPointer(tab.IDivLoPointer),
    IDivUpPointer(tab.IDivUpPointer)
{
   //! Copy constructor
   SetClassName("fastNLOTable");
   for (std::size_t i = 0; i < tab.fCoeff.size(); ++i) {
      fCoeff[i] = tab.fCoeff[i]->Clone();
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
   ifstream* strm = OpenFileRead();
   // read header
   ReadHeader(*strm);
   // read scenario
   ReadScenario(*strm);
   // read b-blocks
   ReadCoeffTables(*strm);
   // close stream
   CloseFileRead(*strm);
}


// ___________________________________________________________________________________________________
void fastNLOTable::ReadCoeffTables(istream& table){
   int nblocks = GetNcontrib()+GetNdata();
   for(int i=0;i<nblocks;i++){
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
      debug["ReadRestOfCoeffTable"]<<"Found data table. Now reading in."<<endl;
      fastNLOCoeffData* cN = new fastNLOCoeffData(cB);
      cN->ReadRest(table);
      return cN;
   }
   else if ( fastNLOCoeffMult::CheckCoeffConstants(&cB,quiet) ) {
      debug["ReadRestOfCoeffTable"]<<"Found multiplicative contribution. Now reading in."<<endl;
      fastNLOCoeffMult* cN = new fastNLOCoeffMult(cB);
      cN->ReadRest(table);
      return cN;
   }
   else if ( fastNLOCoeffAddFix::CheckCoeffConstants(&cB,quiet) ) {
      debug["ReadRestOfCoeffTable"]<<"Found additive fixed order contribution (v2.0). Now reading in."<<endl;
      fastNLOCoeffAddFix* cN = new fastNLOCoeffAddFix(cB);
      cN->ReadRest(table);
      return cN;
   }
   else if ( fastNLOCoeffAddFlex::CheckCoeffConstants(&cB,quiet) ) {
      debug["ReadRestOfCoeffTable"]<<"Found additive flexible scale contribution. Now reading in."<<endl;
      fastNLOCoeffAddFlex* cN = new fastNLOCoeffAddFlex(cB,ILOord);
      cN->ReadRest(table);
      return cN;
   }
   else {
      error["ReadRestOfCoeffTable"]<<"Could not identify coefficient table. Print and exiting... "<<endl;
      cB.Print();
      exit(1);
   }
   return NULL;
}


// ___________________________________________________________________________________________________
void fastNLOTable::WriteTable() {
   //! Write fastNLO table to file 'ffilename' (member)
   info["WriteTable"]<<"Writing fastNLO table to file: "<<ffilename<<endl;
   ofstream* table = OpenFileWrite();
   WriteHeader(*table);
   WriteScenario(*table);
   for(int i=0;i<GetNcontrib();i++){
      debug["WriteTable"]<<"Writing coefficient table #"<<i<<endl;
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


// ___________________________________________________________________________________________________
void fastNLOTable::ReadScenario(istream& table){
   table.peek();
   if (table.eof()){
      warn["ReadScenario"]<<"Cannot read from file."<<endl;
   }

   fastNLOTools::ReadMagicNo(table);

   table >> Ipublunits;
   int  NScDescript = 0;
   table >> NScDescript;
   ScDescript.resize(NScDescript);
   char buffer[257];
   table.getline(buffer,256);
   for(int i=0;i<NScDescript;i++){
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
//    LoBin.resize(NObsBin);
//    UpBin.resize(NObsBin);
   for(unsigned int i=0;i<NObsBin;i++){
      Bin[i].resize(NDim);
//       LoBin[i].resize(NDim);
//       UpBin[i].resize(NDim);
      for(int j=NDim-1;j>=0;j--){
         //table >>  LoBin[i][j];
         //if(IDiffBin[j]==2) table >>  UpBin[i][j];
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
      // maxime pre-v2.0 conversion
      //    if ( NDim == 1 ){
      //         double binsize = 1;
      //         if ( IDiffBin[0] == 2 ) binsize *=  UpBin[i][0] - LoBin[i][0];
      //         printf(" binszie bin %d  = %7.4f\n",i,binsize);
      //         BinSize[i] = binsize;

      //    }
      //    else if ( NDim == 2 || NDim == 3 ){
      //       // warning: the variables are exchanged here!
      //       // what is bound[0] corresponds to bingrid2[nBins][nBins2]
      //       // what is bound[1] corresponds to bingrid1[nBins]

      //       double binsize = 1;
      //       // warning: the variables are exchanged here!
      //       // what is DimLabel[0] corresponds to bingrid2[nBins][nBins2]
      //       // what is DimLabel[1] corresponds to bingrid1[nBins]
      //       printf("UpBin[.][0] =  %7.4f, LoBin[.][0] =  %7.4f , UpBin[.][1] =  %7.4f  LoBin[.][1] =  %7.4f\n",
      //             UpBin[i][0],LoBin[i][0],UpBin[i][1],LoBin[i][1]);
      //       if ( IDiffBin[0] == 2 ) binsize *= UpBin[i][0] - LoBin[i][0];
      //       if ( IDiffBin[1] == 2 ) binsize *= UpBin[i][1] - LoBin[i][1];
      //       printf(" binszie 2Dim bin %d  = %7.4f\n",i,binsize);
      //       BinSize[i] = binsize;
      //    }
   }

   table >> INormFlag;
   if(INormFlag>1){
      table >> DenomTable;
   }
   if(INormFlag>0){
      IDivLoPointer.resize(NObsBin);
      IDivUpPointer.resize(NObsBin);
      for(unsigned int i=0;i<NObsBin;i++){
         table >> IDivLoPointer[i];
         table >> IDivUpPointer[i];
      }
   }

   fastNLOTools::ReadMagicNo(table);
   fastNLOTools::PutBackMagicNo(table);
}


// ___________________________________________________________________________________________________
void fastNLOTable::WriteScenario(ostream& table){
   table << tablemagicno << endl;
   table << Ipublunits << endl;
   int NScDescript =  ScDescript.size();
   table << NScDescript << endl;
   for(int i=0;i<NScDescript;i++){
      table << ScDescript[i] << endl;
   }
   table << Ecms << endl;
   table << ILOord << endl;
   table << NObsBin << endl;
   table << NDim << endl;
   for(int i=NDim-1;i>=0;i--){
      table << DimLabel[i] << endl;
   }
   for(int i=NDim-1;i>=0;i--){
      table << IDiffBin[i] << endl;
   }
   for(unsigned int i=0;i<NObsBin;i++){
      for(int j=NDim-1;j>=0;j--){
         table <<  Bin[i][j].first  << endl;
         //         if(IDiffBin[j]==2) table <<  UpBin[i][j]  << endl;
         if(IDiffBin[j]==0 || IDiffBin[j]==2) table <<  Bin[i][j].second  << endl;
      }
   }
   for(unsigned int i=0;i<NObsBin;i++){
     table << BinSize[i]  << endl;
   }

   table << INormFlag << endl;
   if(INormFlag>1){
      table << DenomTable << endl;
   }
   if(INormFlag>0){
      for(unsigned int i=0;i<NObsBin;i++){
         table << IDivLoPointer[i] << endl;
         table << IDivUpPointer[i] << endl;
      }
   }
}


// ___________________________________________________________________________________________________
bool fastNLOTable::IsCompatible(const fastNLOTable& other) const {
   if ( !IsCompatibleHeader(other) ) return false;
   bool potentialcompatible = true;
   if(Ipublunits != other.Ipublunits){
      warn["IsCompatible"]<<"Differing cross section units found: "<<Ipublunits<<" and "<<other.Ipublunits<<endl;
      return false;
   }
   if(ScDescript != other.ScDescript){
      warn["IsCompatible"]<<"Differing scale description found."<<endl;
      potentialcompatible = false;
   }
   if(!cmp(Ecms,other.Ecms)){
      warn["IsCompatible"]<<"Differing center-of-mass energy found: "<<Ecms<<" and "<<other.Ecms<<endl;
      return false;
   }
   if(ILOord != other.ILOord){
      warn["IsCompatible"]<<"Differing ILOord found: "<<ILOord<<" and "<<other.GetLoOrder()<<endl;
      return false;
   }
   if(NObsBin != other.NObsBin){
      warn["IsCompatible"]<<"Differing NObsBin found: "<<NObsBin<<" and "<<other.NObsBin<<endl;
      return false;
   }
   if(NDim != other.NDim){
      warn["IsCompatible"]<<"Differing NDim found: "<<NDim<<" and "<<other.NDim<<endl;
      return false;
   }
   if(DimLabel != other.DimLabel){
      warn["IsCompatible"]<<"Differing label of observables found."<<endl;
      potentialcompatible = false;
   }
   if(IDiffBin != other.IDiffBin){
      warn["IsCompatible"]<<"Differing IDiffBin found."<<endl;
      return false;
   }
   //    if(!cmp(LoBin,other.LoBin)){
   if(!cmp(Bin,other.Bin)){
      warn["IsCompatible"]<<"Differing Bin boundaries found."<<endl;
      return false;
   }
   //    if(!cmp(UpBin,other.UpBin)){
   //       warn["IsCompatible"]<<"Differing UpBin found."<<endl;
   //       return false;
   //    }
   if(!cmp(BinSize,other.BinSize)){
      warn["IsCompatible"]<<"Differing bin sizes found."<<endl;
      return false;
   }
   if(INormFlag != other.INormFlag){
      warn["IsCompatible"]<<"Differing INormFlag found: "<<INormFlag<<" and "<<other.INormFlag<<endl;
      return false;
   }
   if(INormFlag>1){
      if(DenomTable != other.DenomTable){
         warn["IsCompatible"]<<"Differing DenomTable found."<<endl;
         return false;
      }
   }
   if(INormFlag>0){
      for(unsigned int i=0;i<NObsBin;i++){
         if(IDivLoPointer[i] != other.IDivLoPointer[i]){
            warn["IsCompatible"]<<"Differing IDivLoPointer["<<i<<"] found"<<endl;
            return false;
         }
         if(IDivUpPointer[i] != other.IDivUpPointer[i]){
            warn["IsCompatible"]<<"Differing IDivUpPointer["<<i<<"] found."<<endl;
            return false;
         }
      }
   }
   if ( !potentialcompatible ) warn["IsCompatible"]<<"Some labels have differing values, but relevant variables seem to be compatible. Continuing."<<endl;
   return true;

}


// ___________________________________________________________________________________________________
void fastNLOTable::AddTable(const fastNLOTable& other){
   // add another table to this table.
   // Add either further contributions (higher-orders, data, non-pert corr, etc...)
   // or increase statistics of fixed-order calc.
   //
   if ( !IsCompatible(other) ) {
      warn["AddTable"]<<"Table not compatible with this table. Ignoring command."<<endl;
      return;
   }

   // loop over all contributions from 'other'-table
   const bool quiet = true;
   const int nc = other.GetNcontrib() + other.GetNdata();
   for ( int ic=0 ; ic<nc; ic++ ) {
      bool wasAdded = false;

      // is additive?
      if ( other.GetCoeffTable(ic)->GetIAddMultFlag()==0) {
         fastNLOCoeffAddBase* cadd = (fastNLOCoeffAddBase*)other.GetCoeffTable(ic);

         // find compatible contribution, or add
         for (unsigned int j = 0 ; j<fCoeff.size() ; j++) {
            fastNLOCoeffAddBase* lhs = (fastNLOCoeffAddBase*)fCoeff[j];
            if ( lhs->IsCompatible(*cadd) ) { // found compatible table
               if ( wasAdded )
                  error["AddTable"]<<"This contribution was already added. It seems that there is one contribution twice in the table."<<endl;
               else {
                  debug["AddTable"]<<"Summing contribution "<<ic<<" to fCoeff #"<<j<<endl;
                  if ( fastNLOCoeffAddFlex::CheckCoeffConstants(lhs,quiet) )
                     lhs->Add(*cadd);
                  else if ( fastNLOCoeffAddFix::CheckCoeffConstants(lhs,quiet) )
                     lhs->Add(*cadd);
                  wasAdded = true;
               }
            }
         }
      }
      else {
         // check if this data or 'mult' contribution already exists (which should not happen)
         cout<<"todo. Check if data table already exists!."<<endl;
      }

      // couldn't find a corresponding contribution.
      // add this contribution as new contrib.
      if ( !wasAdded ) {
         info["AddTable"]<<"Adding new contribution to table."<<endl;
         fastNLOCoeffBase* add = other.GetCoeffTable(ic);
         if ( fastNLOCoeffData::CheckCoeffConstants(add,quiet) ) {
            add = new fastNLOCoeffData((fastNLOCoeffData&)*add);
            Ndata++;
         }
         else if ( fastNLOCoeffMult::CheckCoeffConstants(add,quiet) )
            add = new fastNLOCoeffMult((fastNLOCoeffMult&)*add);
         else if ( fastNLOCoeffAddFix::CheckCoeffConstants(add,quiet) )
            add = new fastNLOCoeffAddFix((fastNLOCoeffAddFix&)*add);
         else if ( fastNLOCoeffAddFlex::CheckCoeffConstants(add,quiet) )
            add = new fastNLOCoeffAddFlex((fastNLOCoeffAddFlex&)*add);
         CreateCoeffTable(fCoeff.size(),add);
         ///     Ndata++, and ncontrib++, n
      }
   }
}


// // ___________________________________________________________________________________________________
// int fastNLOTable::CreateCoeffBase(int no){
//    //fastNLOCoefficients* blockb = new fastNLOCoefficients(NObsBin,ILOord);
//    fastNLOCoeffBase* blockb = new fastNLOCoeffBase(NObsBin);
//    return CreateCoeffTable(no,blockb);
// }


// ___________________________________________________________________________________________________
//int fastNLOTable::CreateCoeffTable(int no,fastNLOCoefficients *newblockb){
int fastNLOTable::CreateCoeffTable(int no,fastNLOCoeffBase *newblockb){
   if((no+1)>(int)fCoeff.size())
      fCoeff.resize(no+1);
   fCoeff[no] = newblockb;
   //Ncontrib++; // member of fastNLOBase
   Ncontrib = fCoeff.size();
   return 0;
}


// ___________________________________________________________________________________________________
bool fastNLOTable::cmp(const double x1, const double x2) const {
   double norm;
   if (x1>0.){
      norm = x1;
   }else{
      norm = 1.; // If x1 is 0, do not try to calculate relative deviation, use absolute
   }
   return((fabs(x1-x2)/norm)<1e-7);
}

bool fastNLOTable::cmp(const vector<double>& x1,const vector<double>& x2) const {
   bool result = true;
   for(unsigned int i = 0; i<x1.size() ;i++ ){
      result = result & cmp (x1[i],x2[i]);
   }
   return result;
}

bool fastNLOTable::cmp(const vector<vector<double> >& x1, const vector<vector<double> >& x2) const {
   bool result = true;
   for(unsigned int i = 0; i<x1.size() ;i++ ){
      result = result & cmp (x1[i],x2[i]);
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
   //    for(unsigned int i = 0; i<fCoeff.size() ;i++ ){
   //       fCoeff[i]->fILOord = LOOrd; // fCoeff should not need this member!
   //    }
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
      error["SetDimLabel"]<<"Sorry, you have only initialized "<<NDim<<" dimensions, but you want to label a dimension with number "<<iDim<<endl;
      exit(1);
   }
   if ( iDim < 1) {
      error["SetDimLabel"]<<"The dimension must be a natural number. iDim="<<iDim<<endl;
      exit(1);
   }

   if ( DimLabel.size() != NDim ){
      error["SetDimLabel"]<<"You have to call SetNumDiffBin with a reasonable number before."<<endl;
      exit(1);
   }

   DimLabel[iDim] = label;
   IDiffBin[iDim] = IsDiff ? 2 : 0 ;
}


// ___________________________________________________________________________________________________
//fastNLOCoefficients* fastNLOTable::GetCoeffTable(int no) const {
fastNLOCoeffBase* fastNLOTable::GetCoeffTable(int no) const {
   if ( no >= (int)fCoeff.size() ){
      warn["GetCoeffTable"]<<"There is no contribution with number "<<no<<" but only "<<fCoeff.size()<<". Returning null pointer."<<endl;
      return NULL;
   }
   else
      return fCoeff[no];
}


// ___________________________________________________________________________________________________
fastNLOCoeffData* fastNLOTable::GetDataTable() const {
   for (unsigned int i= 0; i<fCoeff.size() ; i++ ){
      fastNLOCoeffBase* c = GetCoeffTable(i);
      if ( fastNLOCoeffData::CheckCoeffConstants(c,true) ) {
         return (fastNLOCoeffData*)c;
      }
   }
   return NULL;
}


// ___________________________________________________________________________________________________
fastNLOCoeffAddBase* fastNLOTable::GetReferenceTable(ESMOrder eOrder) const {
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
void fastNLOTable::Print() const {
   fastNLOBase::Print();
   PrintScenario();
}


// ___________________________________________________________________________________________________
void fastNLOTable::PrintScenario() const {
  printf("\n **************** FastNLO Table: Scenario ****************\n\n");
  printf("    Ipublunits                    %d\n",Ipublunits);
  for(unsigned int i=0;i<ScDescript.size();i++){
    printf("    ScDescript[%d]                 %s\n",i,ScDescript[i].data());
  }
  printf("    Ecms                          %7.4f\n",Ecms);
  printf("    ILOord                        %d\n",ILOord);
  printf("    NDim                          %d\n",NDim);
  for(unsigned int i=0;i<NDim;i++){
    printf("     - DimLabel[%d]                %s\n",i,DimLabel[i].data());
  }
  for(unsigned int i=0;i<NDim;i++){
    printf("     - IDiffBin[%d]               %d\n",i,IDiffBin[i]);
  }
  printf("    NObsBin                       %d\n",NObsBin);
  for(unsigned int i=0;i<NObsBin;i++){
    for(unsigned int j=0;j<NDim;j++){
      printf("     -  - LoBin[%d][%d]             %7.4f\n", i,j,Bin[i][j].first);
      if(IDiffBin[j]==2)
        printf("     -  - UpBin[%d][%d]             %7.4f\n", i,j,Bin[i][j].second);
    }
   }
   for(unsigned int i=0;i<NObsBin;i++){
     printf("     - BinSize[%d]                %7.4f\n", i,BinSize[i]);
   }
   printf("    INormFlag                     %d\n",INormFlag);

   if(INormFlag>1){
     printf("    DenomTable                    %s\n",DenomTable.data());
   }
   if(INormFlag>0){
      for(unsigned int i=0;i<NObsBin;i++){
        printf("     - IDivLoPointer[%d]               %d\n",i,IDivLoPointer[i]);
        printf("     - IDivUpPointer[%d]               %d\n",i,IDivUpPointer[i]);
      }
   }
   printf("\n ********************************************************\n\n");

}


//______________________________________________________________________________
void fastNLOTable::PrintFastNLOTableConstants(const int iprint) const {
   //if ( debug.GetSpeak() ) iprint=10000;
   //
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
   // Print basic scenario information (always)
   //
   char buffer[1024];
   cout  << endl;
   cout  << _CSEPSC << endl;
   snprintf(buffer, sizeof(buffer), "Information on fastNLO scenario: %s",ScenName.data());
   shout << buffer << endl;
   cout  << _SSEPSC << endl;
   snprintf(buffer, sizeof(buffer), "Description:");
   shout << buffer << endl;
   for (unsigned int i=0; i<ScDescript.size(); i++) {
      printf(" #   %s\n",ScDescript[i].data());
   }
   printf(" #\n");
   printf(" # Centre-of-mass energy Ecms: % -#10.4g GeV\n",Ecms);
   printf(" #\n");
   printf(" # Tot. no. of observable bins: %3i in %1i dimensions:\n",NObsBin,NDim);
   printf(" #\n");
   printf(" # No. of contributions: %1i\n",Ncontrib);

   //
   // Print basic contribution information (always)
   //
   for (unsigned int j = 0 ; j<fCoeff.size() ; j++) {
      fastNLOCoeffBase* c = fCoeff[j];
      if ( iprint == 0 ) {
         printf(" # Contribution %1i:\n",j+1);
         for (unsigned int i=0; i<c->CtrbDescript.size(); i++) {
            printf(" #   %s\n",c->CtrbDescript[i].data());
         }
         if ( fastNLOCoeffAddBase::CheckCoeffConstants(c,true) ) {
            fastNLOCoeffAddBase* cA = (fastNLOCoeffAddBase*)c;
            double nevt = cA->GetNevt();
            printf(" #   No. of events: %#17.0F\n",nevt);
         }
         printf(" #   provided by:\n");
         for (unsigned int i=0; i<c->CodeDescript.size(); i++) {
            printf(" #   %s\n",c->CodeDescript[i].data());
         }

         if ( fastNLOCoeffAddFix::CheckCoeffConstants(c,true) ) {
            fastNLOCoeffAddFix* cFix = (fastNLOCoeffAddFix*)c;
            int NScaleDim = cFix->GetNScaleDim();
            printf(" #   Scale dimensions: %1i\n",NScaleDim);
            for (int i=0; i<NScaleDim; i++) {
               for (unsigned int j=0; j<cFix->ScaleDescript[i].size(); j++) {
                  printf(" #     Scale description for dimension %1i:          %s\n",i+1,cFix->ScaleDescript[i][j].data());
               }
               printf(" #     Number of scale variations for dimension %1i: %1i\n",NScaleDim,cFix->GetNScalevar());
               printf(" #     Available scale settings for dimension %1i:\n",NScaleDim);
               for (int k=0; k<cFix->GetNScalevar(); k++) { // fastNLOReader has method: 'GetNScaleVariations()', which returns nr of scale variations for all active tables!
                  printf(" #       Scale factor number %1i:                   % #10.4f\n",k+1,cFix->GetScaleFactor(k));
               }
               printf(" #     Number of scale nodes for dimension %1i:      %1i\n",NScaleDim,cFix->GetNScaleNode());
            }
         } else if ( fastNLOCoeffAddFlex::CheckCoeffConstants(c,true) ) {
            fastNLOCoeffAddFlex* cFlex = (fastNLOCoeffAddFlex*)c;
            int NScaleDim = cFlex->GetNScaleDim();
            if (! (NScaleDim == 1)) {
               error["PrintFastNLOTableConstants"] << "Flex-scale tables must have scale dimensions of one, aborting! "
                                                   << "NScaleDim = " << NScaleDim <<endl;
               exit(1);
            }
            //            printf(" #   Scale dimensions: %1i\n",NScaleDim);
            for (int i=0; i<NScaleDim; i++) {
               //               unsigned int NFlexScales = cFlex->ScaleDescript[i].size();
               //               for (unsigned int j=0; j<NFlexScales; j++) {
               //                  printf(" #   Scale description for flexible scale %1i:          %s\n",j+1,cFlex->ScaleDescript[i][j].data());
               //               }
               printf(" #   Scale description for flexible scale %1i:          %s\n",1,cFlex->ScaleDescript[i][0].data());
               printf(" #     Number of scale nodes in first observable bin:      %2i\n",cFlex->GetNScaleNode1(0));
               printf(" #     Number of scale nodes in last  observable bin:      %2i\n",cFlex->GetNScaleNode1(NObsBin-1));
               printf(" #   Scale description for flexible scale %1i:          %s\n",2,cFlex->ScaleDescript[i][1].data());
               printf(" #     Number of scale nodes in first observable bin:      %2i\n",cFlex->GetNScaleNode2(0));
               printf(" #     Number of scale nodes in last  observable bin:      %2i\n",cFlex->GetNScaleNode2(NObsBin-1));
            }
         } else {
            // Anything else to write for multiplicative or data contributions?
         }
      } else {
         fCoeff[j]->Print();
      }
   }
   if (iprint > 0) {
      Print();
   }
   shout << "" << endl;
   cout  << _CSEPSC << endl;
}


//______________________________________________________________________________
void fastNLOTable::PrintTableInfo(const int iprint) const {
   debug["PrintTableInfo"] << "Printing flag iprint = " << iprint << endl;
   //
   //  Print basic info about fastNLO table and its contributions
   //   - iprint: iprint > 0: print also contribution descriptions
   //

   char buffer[1024];
   cout  << endl;
   cout  << _CSEPSC << endl;
   shout << "Overview on contribution types and numbers contained in table:" << endl;
   cout  << _SSEPSC << endl;
   snprintf(buffer, sizeof(buffer), "Number of contributions: %2i", Ncontrib);
   shout << buffer << endl;

   int iccount[21] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
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
      // KR: How can I access ContrId from here ???
      iccount[ictype]++;
      shout << "  No.: " << j+1 << ", type: " << coeffname <<", Id: " << iccount[ictype]
            << ", order: " << c->GetContributionDescription()[0]
            << ", by: " << c->GetCodeDescription()[0] << endl;
      if (iprint > 0) {
         for (unsigned int k = 0 ; k<c->GetCodeDescription().size(); k++) {
            snprintf(buffer, sizeof(buffer), "\t\t%s",c->GetCodeDescription()[k].c_str());
            shout << buffer << endl;
         }
      }
   }
   cout << _CSEPSC << endl;

}


// ___________________________________________________________________________________________________
// Collect all binning structure stuff
// ___________________________________________________________________________________________________

// ___________________________________________________________________________________________________
vector < pair <double, double > > fastNLOTable::GetDimBinBoundaries(unsigned int iDim) const {
   //! Get all bins for given dimension 'iDim'
   std::vector< std::pair<double, double > > Bins;
   for (size_t i = 0; i < Bin.size(); ++i)
      Bins.push_back(Bin[i][iDim]);
   return Bins;
}


// ___________________________________________________________________________________________________
vector < pair < double, double > > fastNLOTable::GetDim0BinBoundaries() const {
   //! Get binning of first dimension
   std::vector< std::pair<double, double > > Bins = GetDimBinBoundaries(0);
   std::set< pair< double,double> > set (Bins.begin(), Bins.end());
   Bins.assign(set.begin(),set.end());
   return Bins;
}


// ___________________________________________________________________________________________________
vector < pair < double, double > > fastNLOTable::GetDim1BinBoundaries(unsigned int iDim0Bin) const {
   //! Get binning of second dimension for bin 'iDim0Bin' of first dimension
   std::vector< std::pair<double, double > > Bins;
   const int idiff = GetNumDiffBin();
   if ( idiff < 2 ) {
      error["fastNLOTable::GetDim1Bins"] << "No second dimension available, aborted!" << endl;
      exit(1);
   }
   pair< double, double> bin0 = GetDim0BinBoundaries()[iDim0Bin];
   for (size_t iobs = 0; iobs < Bin.size(); iobs++) {
      if (Bin[iobs][0] == bin0) {
         Bins.push_back(Bin[iobs][1]);
      }
   }
   std::set< pair< double,double> > set (Bins.begin(), Bins.end());
   Bins.assign(set.begin(),set.end());
   return Bins;
}


// ___________________________________________________________________________________________________
vector < pair < double, double > > fastNLOTable::GetDim2BinBoundaries(unsigned int iDim0Bin, unsigned int iDim1Bin) const {
   //! Get binning of third dimension for bins 'iDim0Bin' and 'iDim1Bin' of first two dimensions
   std::vector< std::pair<double, double > > Bins;
   const int idiff = GetNumDiffBin();
   if ( idiff < 3 ) {
      error["fastNLOTable::GetDim2Bins"] << "No third dimension available, aborted!" << endl;
      exit(1);
   }
   pair< double, double> bin0 = GetDim0BinBoundaries()[iDim0Bin];
   pair< double, double> bin1 = GetDim1BinBoundaries(iDim0Bin)[iDim1Bin];
   for (size_t iobs = 0; iobs < Bin.size(); iobs++) {
      if (Bin[iobs][0] == bin0 && Bin[iobs][1] == bin1) {
         Bins.push_back(Bin[iobs][2]);
      }
   }
   std::set< pair< double,double> > set (Bins.begin(), Bins.end());
   Bins.assign(set.begin(),set.end());
   return Bins;
}


// ___________________________________________________________________________________________________
vector < double > fastNLOTable::GetLoBin( int dimension) const {
   //! Get lower bin edge of all observable bins for dimension 'dimension'
   vector < double > LoBin;
   for (size_t i = 0; i < Bin.size(); ++i)
      LoBin.push_back(Bin[i][dimension].first);
   return LoBin;
}


// ___________________________________________________________________________________________________
vector < double > fastNLOTable::GetUpBin( int dimension) const {
   //! Get upper bin edge of all observable bins for dimension 'dimension'
   vector < double > UpBin;
   for (size_t i = 0; i < Bin.size(); ++i)
      UpBin.push_back(Bin[i][dimension].second);
   return UpBin;
}


// ___________________________________________________________________________________________________
unsigned int fastNLOTable::GetIDim0Bin(unsigned int iObsBin) const {
   //! Returns bin number in first dimension
   //! Valid for up to triple differential binnings
   //  There always must be at least one bin!
   if ( Bin.size() == 0 || Bin[0].size() == 0 ) {
      error["fastNLOTable::GetIDim0Bin"] << "No observable bins defined, aborted!" << endl;
      exit(1);
   }
   if ( iObsBin >= NObsBin ) {
      error["fastNLOTable::GetIDim0Bin"] << "Observable bin out of range, aborted!" << endl;
      exit(1);
   }
   unsigned int i0bin = 0;
   double lo0bin = Bin[0][0].first;
   for ( unsigned int iobs = 0; iobs<Bin.size(); iobs++ ) {
      if ( lo0bin < Bin[iobs][0].first ) {
         lo0bin = Bin[iobs][0].first;
         i0bin++;
      }
      if ( iobs == iObsBin ) {
         return i0bin;
      }
   }
   error["fastNLOTable::GetIDim0Bin"] << "Observable bin not found. This should never happen, aborted!" << endl;
   exit(1);
}


// ___________________________________________________________________________________________________
unsigned int fastNLOTable::GetNDim0Bins() const {
   //! Returns number of bins in first dimension
   //! Valid for up to triple differential binnings
   //  There always must be at least one bin!
   return (GetIDim0Bin(NObsBin-1) + 1);
}


// ___________________________________________________________________________________________________
unsigned int fastNLOTable::GetIDim1Bin(unsigned int iObsBin) const {
   //! Returns bin number in second dimension
   //! Valid for up to triple differential binnings
   //  1d binning --> error exit
   const int idiff = GetNumDiffBin();
   if ( idiff < 2 ) {
      error["fastNLOTable::GetIDim1Bin"] << "No second dimension available, aborted!" << endl;
      exit(1);
   }
   //  Otherwise there always must be at least one bin!
   if ( Bin.size() == 0 || Bin[0].size() == 0 ) {
      error["fastNLOTable::GetIDim1Bin"] << "No observable bins defined, aborted!" << endl;
      exit(1);
   }
   if ( iObsBin >= NObsBin ) {
      error["fastNLOTable::GetIDim1Bin"] << "Observable bin out of range, aborted!" << endl;
      exit(1);
   }
   unsigned int i0bin = 0;
   unsigned int i1bin = 0;
   double lo0bin = Bin[0][0].first;
   double lo1bin = Bin[0][1].first;
   for ( unsigned int iobs = 0; iobs<Bin.size(); iobs++ ) {
      if ( lo0bin < Bin[iobs][0].first ) {
         lo0bin = Bin[iobs][0].first;
         lo1bin = Bin[iobs][1].first;
         i0bin++;
         i1bin = 0;
      } else if ( lo1bin < Bin[iobs][1].first ) {
         lo1bin = Bin[iobs][1].first;
         i1bin++;
      }
      if ( iobs == iObsBin ) {
         return i1bin;
      }
   }
   error["fastNLOTable::GetIDim1Bin"] << "Observable bin not found. This should never happen, aborted!" << endl;
   exit(1);
}


// ___________________________________________________________________________________________________
unsigned int fastNLOTable::GetNDim1Bins(unsigned int iDim0Bin) const {
   //! Returns number of bins in second dimension for iDim0Bin in first dimension
   //! Valid for up to triple differential binnings
   //  1d binning --> error exit
   const int idiff = GetNumDiffBin();
   if ( idiff < 2 ) {
      error["fastNLOTable::GetNDim1Bins"] << "No second dimension available, aborted!" << endl;
      exit(1);
   }
   for ( unsigned int iobs = 0; iobs<Bin.size(); iobs++ ) {
      if ( GetIDim0Bin(iobs) == iDim0Bin +1) {
         return GetIDim1Bin(iobs-1) +1;
      } else if ( iobs == Bin.size() -1 ) {
         return GetIDim1Bin(iobs) +1;
      }
   }
   error["fastNLOTable::GetNDim1Bins"] << "Observable bin not found. This should never happen, aborted!" << endl;
   exit(1);
}


// ___________________________________________________________________________________________________
unsigned int fastNLOTable::GetIDim2Bin(unsigned int iObsBin) const {
   //! Returns bin number in third dimension
   //! Valid for up to triple differential binnings
   //  1d, 2d binning --> error exit
   const int idiff = GetNumDiffBin();
   if ( idiff < 3 ) {
      error["fastNLOTable::GetIDim2Bin"] << "No third dimension available, aborted!" << endl;
      exit(1);
   }
   //  Otherwise there always must be at least one bin!
   if ( Bin.size() == 0 || Bin[0].size() == 0 ) {
      error["fastNLOTable::GetIDim2Bin"] << "No observable bins defined, aborted!" << endl;
      exit(1);
   }
   if ( iObsBin >= NObsBin ) {
      error["fastNLOTable::GetIDim2Bin"] << "Observable bin out of range, aborted!" << endl;
      exit(1);
   }
   unsigned int i0bin = 0;
   unsigned int i1bin = 0;
   unsigned int i2bin = 0;
   double lo0bin = Bin[0][0].first;
   double lo1bin = Bin[0][1].first;
   double lo2bin = Bin[0][2].first;
   for ( unsigned int iobs = 0; iobs<Bin.size(); iobs++ ) {
      if ( lo0bin < Bin[iobs][0].first ) {
         lo0bin = Bin[iobs][0].first;
         lo1bin = Bin[iobs][1].first;
         lo2bin = Bin[iobs][2].first;
         i0bin++;
         i1bin = 0;
         i2bin = 0;
      } else if ( lo1bin < Bin[iobs][1].first ) {
         lo1bin = Bin[iobs][1].first;
         lo2bin = Bin[iobs][2].first;
         i1bin++;
         i2bin = 0;
      } else if ( lo2bin < Bin[iobs][2].first ) {
         lo2bin = Bin[iobs][2].first;
         i2bin++;
      }
      if ( iobs == iObsBin ) {
         return i2bin;
      }
   }
   error["fastNLOTable::GetIDim2Bin"] << "Observable bin not found. This should never happen, aborted!" << endl;
   exit(1);
}


// ___________________________________________________________________________________________________
unsigned int fastNLOTable::GetNDim2Bins(unsigned int iDim0Bin, unsigned int iDim1Bin) const {
   //! Returns number of bins in third dimension for iDim0Bin in first and iDim1Bin in second dimension
   //! Valid for up to triple differential binnings
   //  1d, 2d binning --> error exit
   const int idiff = GetNumDiffBin();
   if ( idiff < 3 ) {
      error["fastNLOTable::GetNDim2Bins"] << "No third dimension available, aborted!" << endl;
      exit(1);
   }
   for ( unsigned int iobs = 0; iobs<Bin.size(); iobs++ ) {
      if ( GetIDim0Bin(iobs) == iDim0Bin && GetIDim1Bin(iobs) == iDim1Bin +1) {
         return GetIDim2Bin(iobs-1) +1;
      } else if ( GetIDim0Bin(iobs) == iDim0Bin +1 && GetIDim1Bin(iobs-1) == iDim1Bin) {
         return GetIDim2Bin(iobs-1) +1;
      } else if ( iobs == Bin.size() -1 ) {
         return GetIDim2Bin(iobs) +1;
      }
   }
   error["fastNLOTable::GetNDim2Bins"] << "Observable bin not found. This should never happen, aborted!" << endl;
   exit(1);
}


// ___________________________________________________________________________________________________
// DOES NOT WORK YET!
unsigned int fastNLOTable::GetIDimBin(unsigned int iObsBin, unsigned int iDim) const {
   //! Returns bin number in dimension iDim
   //  iDim larger than table dimensions --> error exit
   const unsigned int idiff = GetNumDiffBin();
   if ( ! (iDim < idiff) ) {
      error["fastNLOTable::GetIDimBin"] << "Requested dimension iDim not available, aborted!" << endl;
      exit(1);
   }
   //  Otherwise there always must be at least one bin!
   if ( Bin.size() == 0 || Bin[0].size() == 0 ) {
      error["fastNLOTable::GetIDimBin"] << "No observable bins defined, aborted!" << endl;
      exit(1);
   }
   if ( iObsBin >= NObsBin ) {
      error["fastNLOTable::GetIDimBin"] << "Observable bin out of range, aborted!" << endl;
      exit(1);
   }
   vector < unsigned int > ibin(NDim);
   vector < double > lobin(NDim);
   for ( unsigned int jdim = 0; jdim<=iDim; jdim++ ) {
      ibin.push_back(0);
      lobin.push_back(Bin[0][jdim].first);
   }
   for ( unsigned int iobs = 0; iobs<Bin.size(); iobs++ ) {
      if (iDim == 2) {
         if ( lobin[0] < Bin[iobs][0].first ) {
            for (unsigned int k=0; k<=iDim; k++) {
               lobin[k] = Bin[iobs][k].first;
            }
            ibin[0]++;
            for (unsigned int k=1; k<=iDim; k++) {
               ibin[k] = 0;
            }
         } else if ( lobin[1] < Bin[iobs][1].first ) {
            for (unsigned int k=1; k<=iDim; k++) {
               lobin[k] = Bin[iobs][k].first;
            }
            ibin[1]++;
            for (unsigned int k=2; k<=iDim; k++) {
               ibin[k] = 0;
            }
         } else if ( lobin[2] < Bin[iobs][2].first ) {
            for (unsigned int k=2; k<=iDim; k++) {
               lobin[k] = Bin[iobs][k].first;
            }
            ibin[2]++;
         }
      } else if (iDim == 1) {
         if ( lobin[0] < Bin[iobs][0].first ) {
            for (unsigned int k=0; k<=iDim; k++) {
               lobin[k] = Bin[iobs][k].first;
            }
            ibin[0]++;
            for (unsigned int k=1; k<=iDim; k++) {
               ibin[k] = 0;
            }
         } else if ( lobin[1] < Bin[iobs][1].first ) {
            for (unsigned int k=1; k<=iDim; k++) {
               lobin[k] = Bin[iobs][k].first;
            }
            ibin[1]++;
         }
      } else if (iDim == 0 ) {
         if ( lobin[0] < Bin[iobs][0].first ) {
            for (unsigned int k=0; k<=iDim; k++) {
               lobin[k] = Bin[iobs][k].first;
            }
            ibin[0]++;
         }
      }
      if ( iobs == iObsBin ) {
         return ibin[iDim];
      }
   }
   error["fastNLOTable::GetIDimBin"] << "Observable bin not found. This should never happen, aborted!" << endl;
   exit(1);
}


// ___________________________________________________________________________________________________
int fastNLOTable::GetODim0Bin(double obs0) const {
   //! Returns bin number in first dimension
   //! Returns -1 if outside range

   int iDim0Bin = -1;
   for ( unsigned int i=0; i<NObsBin; i++ ) {
      if ( IDiffBin[0] == 1 ) {
         error["fastNLOTable::GetODim0Bin"] << "Point-wise differential not yet implemented, aborted!" << endl;
         exit(1);
      } else {
         if ( Bin[i][0].first < obs0 && obs0 <= Bin[i][0].second ) {
            iDim0Bin = GetIDim0Bin(i);
            break;
         }
      }
   }
   return iDim0Bin;
}


// ___________________________________________________________________________________________________
int fastNLOTable::GetODim1Bin(double obs0, double obs1) const {
   //! Returns bin number in second dimension
   //! Returns -1 if outside range

   int iDim1Bin = -1;
   for ( unsigned int i=0; i<NObsBin; i++ ) {
      if ( IDiffBin[0] == 1 ) {
         error["fastNLOTable::GetODim1Bin"] << "Point-wise differential not yet implemented, aborted!" << endl;
         exit(1);
      } else {
         if ( Bin[i][0].first < obs0 && obs0 <= Bin[i][0].second &&
              Bin[i][1].first < obs1 && obs1 <= Bin[i][1].second ) {
            iDim1Bin = GetIDim1Bin(i);
            break;
         }
      }
   }
   return iDim1Bin;
}


// ___________________________________________________________________________________________________
int fastNLOTable::GetODim2Bin(double obs0, double obs1, double obs2) const {
   //! Returns bin number in third dimension
   //! Returns -1 if outside range

   int iDim2Bin = -1;
   for ( unsigned int i=0; i<NObsBin; i++ ) {
      if ( IDiffBin[0] == 1 ) {
         error["fastNLOTable::GetODim2Bin"] << "Point-wise differential not yet implemented, aborted!" << endl;
         exit(1);
      } else {
         if ( Bin[i][0].first < obs0 && obs0 <= Bin[i][0].second &&
              Bin[i][1].first < obs1 && obs1 <= Bin[i][1].second &&
              Bin[i][2].first < obs2 && obs2 <= Bin[i][2].second ) {
            iDim2Bin = GetIDim2Bin(i);
            break;
         }
      }
   }
   return iDim2Bin;
}


// ___________________________________________________________________________________________________
std::vector < std::pair < double, double > > fastNLOTable::GetBinBoundaries(int iDim0Bin, int iDim1Bin, int iDim2Bin) {
   //!
   //! Get bin boundaries for first, second, and third dimension
   //!    Assuming for instance following 2-dimensional binning scheme:
   //!
   //!    iDim0Bin  ________________________________
   //!       0      |___|___|___|_______|__|__|____|
   //!  D    1      |____|____|____|_____|____|____|
   //!  I    2      |__|__|___|__|__|___|__|___|___|
   //!  M    3      |______|_______|_________|_____|
   //!       4      |__________|_______|_________|_|
   //!  0    5      |______________|______|___|____|
   //!       6      |_______|_______|______|_______|
   //!                            DIM 1
   //!  iDim1Bin may be different for each iDim0Bin
   //!
   //! usage e.g.:
   //! int LowerBoundary = GetBinBoundaries(ibin)[dim].first;
   //! int UpperBoundary = GetBinBoundaries(ibin)[dim].second;
   //! 'dim' must be smaller than number of parameters passed to GetBinBoundaries
   //! usage e.g.:
   //! int LoYBin  = GetBinBoundaries(2)[0].first;
   //! int UpPtBin = GetBinBoundaries(2,5)[1].second;
   //! int LoYBin  = GetBinBoundaries(2,5)[0].second;

   cout<<"todo.GetBinBoundaries() ."<<endl;
   vector<pair<double,double> > BinRet(NDim);
   const int idiff = GetNumDiffBin();

   if ( idiff==1 ) {
      if ( iDim0Bin<0 || iDim0Bin >= (int)NObsBin ) {
         warn["GetBinBoundaries"]<<"0th dimension does only have "<<NObsBin<<" but bin "<<iDim0Bin<<" was requested."<<endl;
         return BinRet;
      }
      return Bin[iDim0Bin];
   }
   else if ( idiff==2) {
      unsigned int nDim1 = GetNDim1Bins(iDim0Bin);
      unsigned int nDim0 = GetNDim0Bins();
      // sanity
      if ( iDim0Bin < 0 || iDim0Bin >= (int)nDim0 ) {
         warn["GetBinBoundaries"]<<"Dimension 2 does only have "<<nDim0<<" but bin "<<iDim0Bin<<" was requested."<<endl;
      }
      if ( iDim1Bin < 0 || iDim1Bin >= (int)nDim1 ) {
         warn["GetBinBoundaries"]<<"Dimension 1 does only have "<<nDim1<<" but bin "<<iDim1Bin<<" was requested."<<endl;
      }
      int iObs = 0;
      for ( int i0 = 0 ; i0<iDim0Bin ; i0++ ) {
         iObs += GetNDim1Bins(i0);
      }
      iObs += iDim1Bin;
      return Bin[iObs];
   }
   else if ( idiff == 3 ) {
      unsigned int nDim0 = GetNDim0Bins();
      unsigned int nDim1 = GetNDim1Bins(iDim0Bin);
      unsigned int nDim2 = GetNDim2Bins(iDim0Bin,iDim1Bin);
      // sanity
      if ( iDim0Bin < 0 || iDim0Bin >= (int)nDim0 ) {
         warn["GetBinBoundaries"]<<"Dimension 0 does only have "<<nDim0<<" but bin "<<iDim0Bin<<" was requested."<<endl;
      }
      if ( iDim1Bin < 0 || iDim1Bin >= (int)nDim1 ) {
         warn["GetBinBoundaries"]<<"Dimension 1 does only have "<<nDim1<<" but bin "<<iDim1Bin<<" was requested."<<endl;
      }
      if ( iDim2Bin < 0 || iDim2Bin >= (int)nDim2 ) {
         warn["GetBinBoundaries"]<<"Dimension 2 does only have "<<nDim2<<" but bin "<<iDim2Bin<<" was requested."<<endl;
      }
      error["GetBinBoundaries"]<<"todo. Further code no yet implemented."<<endl;

   }
   else {
      error["GetBinBoundaries"]<<"Higher than triple-differential binnings are not implemented."<<endl;
   }
   return BinRet;

}


// ___________________________________________________________________________________________________
int fastNLOTable::GetObsBinNumber( const vector<double>& vobs ) const {
   //! Returns first matching observable bin number for vector of observations
   //! (assumes none or exactly one matching bin!)
   //! Returns -1 if outside range

   if ( vobs.size() != NDim ) {
      error["GetObsBinNumber"] << "Number of observable values not equal dimensionality of the binning, aborted" << endl;
      error["GetObsBinNumber"] << "NDim = " << NDim << ", vobs.size() = " << vobs.size() << endl;
      exit(1);
   }
   if ( NDim > 3 ) {
      error["fastNLOTable::GetObsBinNumber"] << "More than 3-dimensional binning not yet implemented, aborted!" << endl;
      exit(1);
   }

   for ( unsigned int i=0; i<NObsBin; i++ ) {
      bool lmatch = true;
      for ( unsigned int j=0; j<NDim; j++ ) {
         if ( IDiffBin[j] == 1 ) { // Point-wise differential
            lmatch = lmatch && ( fabs(Bin[i][j].first - vobs[j]) < DBL_MIN );
         } else {                  // Non- or bin-wise differential
            lmatch = lmatch && ( Bin[i][j].first < vobs[j] && vobs[j] <= Bin[i][j].second );
         }
      }
      if ( lmatch ) {return i;}
   }
   return -1;
}
// ___________________________________________________________________________________________________
int fastNLOTable::GetObsBinNumber( double obs0 ) const {
   //! Returns first matching observable bin number for one observation
   //! (assumes none or exactly one matching bin!)
   //! Returns -1 if outside range
   vector < double > vobs(1);
   vobs[0] = obs0;
   return GetObsBinNumber(vobs);
}
// ___________________________________________________________________________________________________
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
// Some info getters
// ___________________________________________________________________________________________________

// ___________________________________________________________________________________________________
std::string fastNLOTable::GetRivetId() const {
   std::string identifier("RIVET_ID");
   std::string found;
   for (size_t i=0; i < ScDescript.size(); ++i) {
      if (ScDescript[i].find(identifier) != std::string::npos){
         size_t RivetIdx = ScDescript[i].find(identifier);
         size_t RivetValIdx = ScDescript[i].find("=", RivetIdx) + 1;
         size_t RivetValLen = ScDescript[i].find(",", RivetValIdx) - RivetValIdx;
         found = ScDescript[i].substr(RivetValIdx, RivetValLen);
         break;
      }
   }
   return found;
}


// ___________________________________________________________________________________________________
// Deprecated
// ___________________________________________________________________________________________________
