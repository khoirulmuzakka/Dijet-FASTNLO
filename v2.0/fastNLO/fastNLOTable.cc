#include <cstdlib>
#include "fastNLOTable.h"


using namespace std;
using namespace fastNLO;

//fastNLOTable::fastNLOTable() : PrimalScream("fastNLOTable") {
fastNLOTable::fastNLOTable(){
   SetClassName("fastNLOTable");
}

fastNLOTable::fastNLOTable(string name) : fastNLOBase(name) {
   SetClassName("fastNLOTable");
   ReadTable();
}

fastNLOTable::~fastNLOTable(){
   // delete fCoeff tables...
}

int fastNLOTable::ReadTable(){
   // open file and read header
   fastNLOBase::ReadTable();
   // read scenario
   ReadScenario(ifilestream);
   // read b-blocks
   ReadCoeffTables(ifilestream);
   return 0;
}

int fastNLOTable::ReadCoeffTables(istream* table){
   int nblocks = GetNcontrib()+GetNdata();
   for(int i=0;i<nblocks;i++){
      fastNLOCoeffBase cTemp(NObsBin);
      cTemp.ReadBase(table);
      fastNLOCoeffBase* cN = ReadRestOfCoeffTable(cTemp, table);
      CreateCoeffTable(i, cN);
   }
   return 0;
}

fastNLOCoeffBase* fastNLOTable::ReadRestOfCoeffTable(const fastNLOCoeffBase& cB, istream *table){
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


int fastNLOTable::WriteTable() {
   info["WriteTable"]<<"Writing fastNLO table to file: "<<ffilename<<endl;
   OpenFileRewrite();
   fastNLOBase::WriteHeader(ofilestream);
   //    for (int k=0;k<GetNcontrib();k++){
   //       GetBlockB(k)->Nevt = (long long int)nevents;
   //    }
   WriteScenario(ofilestream);
   for(int i=0;i<GetNcontrib();i++){
      info["WriteTable"]<<"Writing coefficient table #"<<i<<endl;
      WriteCoeffTableDividebyN(i);
   }
   CloseFileWrite();
   return 0;
}

int fastNLOTable::WriteTable(string filename) {
   string tempfilename = filename;
   SetFilename(filename);
   int ret = WriteTable();
   SetFilename(tempfilename);
   return ret;
}

int fastNLOTable::ReadScenario(istream *table){
   table->peek();
   if (table->eof()){
      warn["ReadScenario"]<<"Cannot read from file."<<endl;
      return(2);
   }

   ReadMagicNo(table);

   *table >> Ipublunits;
   int  NScDescript = 0;
   *table >> NScDescript;
   ScDescript.resize(NScDescript);
   char buffer[257];
   table->getline(buffer,256);
   for(int i=0;i<NScDescript;i++){
      table->getline(buffer,256);
      ScDescript[i] = buffer;
      //      StripWhitespace(ScDescript[i]);
   }

   *table >> Ecms;
   *table >> ILOord;
   *table >> NObsBin;
   *table >> NDim;
   DimLabel.resize(NDim);
   table->getline(buffer,256);
   for(int i=0;i<NDim;i++){
      table->getline(buffer,256);
      DimLabel[i] = buffer;
      //      StripWhitespace(DimLabel[i]);
   }

   IDiffBin.resize(NDim);
   for(int i=0;i<NDim;i++){
      *table >>  IDiffBin[i];
   }
   Bin.resize(NObsBin);
//    LoBin.resize(NObsBin);
//    UpBin.resize(NObsBin);
   //KR: Set rapidity index also when reading a table
   RapIndex.push_back(0);
   //   int irap = 0;
   for(int i=0;i<NObsBin;i++){
      Bin[i].resize(NDim);
//       LoBin[i].resize(NDim);
//       UpBin[i].resize(NDim);
      for(int j=0;j<NDim;j++){
         //*table >>  LoBin[i][j];
         //if(IDiffBin[j]==2) *table >>  UpBin[i][j];
	 *table >>  Bin[i][j].first;
         if(IDiffBin[j]==2) *table >>  Bin[i][j].second;
      }
      //      cout << "iobs1: " << i << ", LoBin i: " << LoBin[i][1] << endl;
      if ( i > 0 ) {
	 //if ( LoBin[i][1] != LoBin[i-1][1] ) {
	if ( Bin[i][1].first != Bin[i-1][1].first ) {
	  //	  cout << "iobs2: " << i << ", LoBin i-1: " << LoBin[i-1][1] << ", LoBin i: " << LoBin[i][1] << endl;
	  RapIndex.push_back(i);
	  //	  irap++;
	  //	  cout << "irap: " << irap << ", RapIndex: " << RapIndex[irap] << endl;
	} 
      }
   }

   BinSize.resize(NObsBin);
   for(int i=0;i<NObsBin;i++){
      *table >> BinSize[i];
      // maxime pre-v2.0 conversion
      //    if ( NDim == 1 ){
      // 	 double binsize = 1;
      // 	 if ( IDiffBin[0] == 2 ) binsize *=  UpBin[i][0] - LoBin[i][0];
      // 	 printf(" binszie bin %d  = %7.4f\n",i,binsize);
      // 	 BinSize[i] = binsize;
	 
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
      // 	     UpBin[i][0],LoBin[i][0],UpBin[i][1],LoBin[i][1]); 
      //       if ( IDiffBin[0] == 2 ) binsize *= UpBin[i][0] - LoBin[i][0];
      //       if ( IDiffBin[1] == 2 ) binsize *= UpBin[i][1] - LoBin[i][1];
      //       printf(" binszie 2Dim bin %d  = %7.4f\n",i,binsize);
      //       BinSize[i] = binsize;
      //    }
   }

   *table >> INormFlag;
   if(INormFlag>1){
      *table >> DenomTable;
   }
   if(INormFlag>0){
      IDivLoPointer.resize(NObsBin);
      IDivUpPointer.resize(NObsBin);
      for(int i=0;i<NObsBin;i++){
         *table >> IDivLoPointer[i];
         *table >> IDivUpPointer[i];
      }
   }

   ReadMagicNo(table);
   PutBackMagicNo(table);
   return 0;
}

int fastNLOTable::WriteScenario(ostream *table){
   *table << tablemagicno << endl;
   *table << Ipublunits << endl;
   int NScDescript =  ScDescript.size();
   *table << NScDescript << endl;
   for(int i=0;i<NScDescript;i++){
      *table << ScDescript[i] << endl;
   }
   *table << Ecms << endl;
   *table << ILOord << endl;
   *table << NObsBin << endl;
   *table << NDim << endl;
   for(int i=0;i<NDim;i++){
      *table << DimLabel[i] << endl;
   }
   for(int i=0;i<NDim;i++){
      *table << IDiffBin[i] << endl;
   }
   for(int i=0;i<NObsBin;i++){
      for(int j=0;j<NDim;j++){
         *table <<  Bin[i][j].first  << endl;
	 //         if(IDiffBin[j]==2) *table <<  UpBin[i][j]  << endl;
	 if(IDiffBin[j]==2) *table <<  Bin[i][j].second  << endl;
      }
   }
   for(int i=0;i<NObsBin;i++){
     *table << BinSize[i]  << endl;
   }

   *table << INormFlag << endl;
   if(INormFlag>1){
      *table << DenomTable << endl;
   }
   if(INormFlag>0){
      for(int i=0;i<NObsBin;i++){
         *table << IDivLoPointer[i] << endl;
         *table << IDivUpPointer[i] << endl;
      }
   }
   return 0;
}


bool fastNLOTable::IsCompatible(fastNLOTable* other) const {
   if(Ipublunits != other->Ipublunits){
      warn["IsCompatible"]<<"Differing cross section units found: "<<Ipublunits<<" and "<<other->Ipublunits<<endl;
      return false;
   }
   if(ScDescript != other->ScDescript){
      warn["IsCompatible"]<<"Differing ScDescript found."<<endl;
      return false;
   }  
   if(!cmp(Ecms,other->Ecms)){
      warn["IsCompatible"]<<"Differing Ecms found: "<<Ecms<<" and "<<other->Ecms<<endl;
      return false;
   }  
   if(ILOord != other->ILOord){
      warn["IsCompatible"]<<"Differing ILOord found: "<<ILOord<<" and "<<other->GetLoOrder()<<endl;
      return false;
   }  
   if(NObsBin != other->NObsBin){
      warn["IsCompatible"]<<"Differing NObsBin found: "<<NObsBin<<" and "<<other->NObsBin<<endl;
      return false;
   }  
   if(NDim != other->NDim){
      warn["IsCompatible"]<<"Differing NDim found: "<<NDim<<" and "<<other->NDim<<endl;
      return false;
   }  
   if(DimLabel != other->DimLabel){
      warn["IsCompatible"]<<"Differing DimLabel found."<<endl;
      return false;
   }  
   if(IDiffBin != other->IDiffBin){
      warn["IsCompatible"]<<"Differing IDiffBin found."<<endl;
      return false;
   }  
//    if(!cmp(LoBin,other->LoBin)){
   if(!cmp(Bin,other->Bin)){
      warn["IsCompatible"]<<"Differing Bin boundaries found."<<endl;
      return false;
   }
   //    if(!cmp(UpBin,other->UpBin)){
   //       warn["IsCompatible"]<<"Differing UpBin found."<<endl;
   //       return false;
   //    }
   if(!cmp(BinSize,other->BinSize)){
      warn["IsCompatible"]<<"Differing BinSize found."<<endl;
      return false;
   }
   if(INormFlag != other->INormFlag){
      warn["IsCompatible"]<<"Differing INormFlag found: "<<INormFlag<<" and "<<other->INormFlag<<endl;
      return false;
   }  
   if(INormFlag>1){
      if(DenomTable != other->DenomTable){
         warn["IsCompatible"]<<"Differing DenomTable found."<<endl;
         return false;
      }
   }
   if(INormFlag>0){
      for(int i=0;i<NObsBin;i++){
         if(IDivLoPointer[i] != other->IDivLoPointer[i]){
            warn["IsCompatible"]<<"Differing IDivLoPointer["<<i<<"] found"<<endl;
            return false;
         }
         if(IDivUpPointer[i] != other->IDivUpPointer[i]){
            warn["IsCompatible"]<<"Differing IDivUpPointer["<<i<<"] found."<<endl;
            return false;
         }
      }
   }  
   
   return true;
};


// ___________________________________________________________________________________________________
int fastNLOTable::CreateCoeffBase(int no){
   //fastNLOCoefficients* blockb = new fastNLOCoefficients(NObsBin,ILOord);
   fastNLOCoeffBase* blockb = new fastNLOCoeffBase(NObsBin);
   return CreateCoeffTable(no,blockb);
}
//int fastNLOTable::CreateCoeffTable(int no,fastNLOCoefficients *newblockb){
int fastNLOTable::CreateCoeffTable(int no,fastNLOCoeffBase *newblockb){
   if((no+1)>(int)fCoeff.size())
      fCoeff.resize(no+1);
   fCoeff[no] = newblockb;
   //Ncontrib++; // member of fastNLOBase
   Ncontrib = fCoeff.size();
   return 0;
}
int fastNLOTable::WriteCoeffTable(int no){
   if((no)<(int)fCoeff.size()){
      return fCoeff[no]->Write(ofilestream);
   }else{
      error["WriteCoeffTable"]<<"Table no. "<<no<<" does not exist, only up to "<<fCoeff.size()<<". Stopping."<<endl;
      exit(2);
   }
}
int fastNLOTable::WriteCoeffTable(int no,ofstream* outstream ){
   if((no)<(int)fCoeff.size()){
      return fCoeff[no]->Write(outstream);
   }else{
      error["WriteCoeffTable"]<<"Table no. "<<no<<" does not exist, only up to "<<fCoeff.size()<<". Stopping."<<endl;
      exit(2);
   }
}
int fastNLOTable::WriteCoeffTableDividebyN(int no){
   if((no)<(int)fCoeff.size()){
      //return fCoeff[no]->Write(ofilestream,fastNLOCoefficients::DividebyNevt);
      fCoeff[no]->Print();
      return fCoeff[no]->Write(ofilestream,fastNLOCoeffBase::DividebyNevt);
   }else{
      error["WriteCoeffTable"]<<"Table no. "<<no<<" does not exist, only up to "<<fCoeff.size()<<". Stopping."<<endl;
      exit(2);
   }
}
void fastNLOTable::DeleteAllCoeffTable(){
   for (size_t i = 0; i < fCoeff.size(); ++i)
      delete fCoeff[i];
   fCoeff.clear();
}


// ___________________________________________________________________________________________________
bool fastNLOTable::cmp(const double x1, const double x2) const{
   double norm;
   if (x1>0.){
      norm = x1;
   }else{
      norm = 1.; // If x1 is 0, do not try to calculate relative deviation, use absolute
   }
   return((fabs(x1-x2)/norm)<1e-7);
}

bool fastNLOTable::cmp(const vector < double > x1,const vector < double > x2) const{
   bool result = true;
   for(unsigned int i = 0; i<x1.size() ;i++ ){
      result = result & cmp (x1[i],x2[i]);
   }
   return result;
}

bool fastNLOTable::cmp(vector < vector < double > > x1,  vector < vector < double > > x2) const{
   bool result = true;
   for(unsigned int i = 0; i<x1.size() ;i++ ){
      result = result & cmp (x1[i],x2[i]);
   }
   return result;
}

bool fastNLOTable::cmp(vector < vector < pair<double,double > > > x1,  vector < vector < pair<double,double > > > x2) const{
   bool result = true;
   for(unsigned int i = 0; i<x1.size() ;i++ ){
      for(unsigned int j = 0; j<x1[i].size() ;j++ ){
	 result = result & (cmp(x1[i][j].first,x2[i][j].first) && cmp(x1[i][j].second,x2[i][j].second));
      }
   }
   return result;
}


void fastNLOTable::SetLoOrder(int LOOrd){
   ILOord = LOOrd;
   //    for(unsigned int i = 0; i<fCoeff.size() ;i++ ){
   //       fCoeff[i]->fILOord = LOOrd; // fCoeff should not need this member!
   //    }
}


void fastNLOTable::SetDimLabel( string label, int iDim , bool IsDiff ){
   // Set label for dimension
   // 
   // In this method, we also set IDiffBin.
   // IDiffBin defines, if this dimension is a ('truely') differential (=1) oder 
   // binned distribution (=2).
   // Since we assume here (in ::SetDimLabel and ::InitBinning) that we only
   // use binned distributions, we use IDiffBin to identify, if the publication
   // was divided by this bin(-width) (=2) or not.
   // At the end of ::InitBinning() we have to set IDiffBin then always to 2
   // to identify this dimension to be a 'binned' dimension.
   //


   // check validity of call
   if ( NDim < iDim ) {
      error["SetDimLabel"]<<"Sorry, you have only initialized "<<NDim<<" dimensions, but you want to label a dimension with number "<<iDim<<endl;
      exit(1);
   }
   if ( iDim < 1) {
      error["SetDimLabel"]<<"The dimension must be a natural number. iDim="<<iDim<<endl;
      exit(1);
   }
   
   if ( (int)DimLabel.size() != NDim ){
      error["SetDimLabel"]<<"You have to call SetNumDiffBin with a reasonable number before."<<endl;
      exit(1);
   }

   DimLabel[iDim-1] = label;
   IDiffBin[iDim-1] = IsDiff ? 2 : 0 ;
}


//fastNLOCoefficients* fastNLOTable::GetCoeffTable(int no) const {
fastNLOCoeffBase* fastNLOTable::GetCoeffTable(int no) const {
   if ( no >= (int)fCoeff.size() ){
      warn["GetCoeffTable"]<<"There is no contribution with number "<<no<<" but only "<<fCoeff.size()<<". Returning null pointer."<<endl;
      return NULL;
   }
   else
      return fCoeff[no];
}


int fastNLOTable::GetBinNumber( double val1 , double val2 ) const {
   // Get Bin number of this event if you use a single or double differential binning
   // return -1 if no bin was found

   //
   //  calculate the bin number as define in Scenario::LoBin and Scenario::UpBin
   //  initialized by Scenario::InitBinning.
   //
   //  returns the bin number, that has to be passse to FillEvent()
   //  return -1 if values are out of bin-ranges
   //

   //    if ( val2==-42 && NDim!=1){
   //       printf("fastNLOTable::GetBinNumber(%6.3f,%6.3f). Error. A single differential table only has one variable.\n",var1,var2);exit(1);
   //    }
   //    if ( val2!==-42 && NDim!=2){
   //       printf("fastNLOTable::GetBinNumber(%6.3f,%6.3f). Error. A double differential table only needs two variables.\n",var1,var2);exit(1);
   //    }
   
   if ( (val2 == -42 && NDim != 1 ) || ( val2 != -42 && NDim != 2 ) ) {
      error["GetBinNumber"]<<"This function can calculate the bin number only for single and double differential binnings (NDim = "<<NDim<<" )."<<endl;
      exit(1);
   }
   int obsbin = -1;

   static const double eps = 1.e-8;
   if ( NDim == 2 ) {
      if ( IDiffBin[0] != 1 || IDiffBin[1] != 1 ) { // bin integrated calculation
 	 for(int j = 0; j < NObsBin; j++) {
	    if ( val1 >= Bin[j][0].first  && val1 <  Bin[j][0].second &&
		 val2 >= Bin[j][1].first  && val2 <  Bin[j][1].second) {
	       obsbin=j;
	       break;
	    }
	 }
      }
      else {  // truly differential calculation
	 for(int j = 0; j < NObsBin; j++) {
	    if ( fabs(val1 - Bin[j][0].first ) < eps  && fabs(val2 - Bin[j][1].first) < eps ) {
	       obsbin=j;
	       break;
	    }
	 }
      }
   }
   else if ( NDim == 1 ) {
      if ( IDiffBin[0] != 1 ) { // bin integrated calculation
	 for(int j = 0; j < NObsBin; j++) {
	    if ( val1 >= Bin[j][0].first  && val1 <  Bin[j][0].second ){
	       obsbin = j;
	       break;
	    }
	 }
      }
      else {  // truly differential calculation
	 for(int j = 0; j < NObsBin; j++) {
	    if ( fabs(val1 - Bin[j][0].first ) < eps ) {
	       obsbin=j;
	       break;
	    }
	 }
      }
   }
   else {
      error["GetBinNumber"]<<"("<<val1<<","<<val2<<"). Error. Only single, or double differential tables are supported. Scenario::NDim = "<<NDim<<"."<<endl;
      exit(1);
   }
   //cout<<"--- calcbin: val1="<<val1<<", val2="<<val2<<", obsbin="<<obsbin<<endl;
   return obsbin;

}


/*
void fastNLOTable::InitBinning( const int nBins1 , double* bingrid1 , const int* nBins2  , vector<double*> bingrid2  , double binwidth3 ){

   // ------------------------------------------------------------------- //
   //
   //  InitBinning. This method initalizes and (partly) checks the binning
   //  that is used in the scenario.
   //
   //   We must know NDim before calling InitBinning.
   //	NDim tells us, in how many dimensions/variables the measurement was performed
   //   this method only supports two dimensional measurements and a third dimension for 
   //   a pseudo-dimensional binning (if only one bin was measured e.g. in the pseudorapidity).
   //   Still, fastNLO could support higher dimensional binnings.
   // 
   //  input.
   //     nBins1	number of bins in 1st dimension
   //     bingrid1	binning in 1st dimension
   //     nBins2	number of bins of second dimension for each 1st-dimension variable
   //     bingrid	binning in 2nd dimension for each 1st dimension bin
   //     binwidth3	binwidth for a 3rd dimension. If the publ. cross sections are e.g. divided by the eta-range.
   //			   if this is dependent on 1st or 2nd dimension binning, this method has to be updated.
   //                   or you can use binwidth3 as a scalling factor if your binning is e.g. in TeV, but you want to have pb/GeV
   //
   //  output.
   //     no output
   //
   //  initalized values
   //     LoBin, UpBin, BinSize, NObsBin
   //
   //
   //  logic for binwidth
   //     using IDiffBin you can specify, if the publ cross section was divided by this bin.
   //     if IDiffBin==2, then you divide the 'final' fastNLO-cross section by this number too -> we multiply the binwidth by this 
   //        dimensional width.
   //     if IDiffBin==1. then the cross section is not divided by this dimensional-binning width. However, we store
   //        the bingrid since the 'ObsBin' is binned in this binning.
   //        MENTION: the UpBin is the NOT stored in the table!
   //
   // ------------------------------------------------------------------- //


   if ( NDim == 2 && ( nBins2==NULL || bingrid2.empty()) ) printf("fastNLOTable::InitBinning. Error. NDim=2, but you have defined only one dimensional bin grid\n");
   if ( NDim == 3 ) printf("fastNLOTable::InitBinning. Error. NDim=3, you must define the binwidth of dimension 3!\n");
   //if ( NDim != 1 && NDim != 2 )  printf("fastNLOTable::InitBinning. Error. NDim must be 1 or 2. three is not yet fully implemented.\n");

   vector <double> bound(NDim);
   int nbins = 0;   // --- count total No. bins

   if ( (int)DimLabel.size() != NDim ) printf("Error. you do not have the same number of DimLabel than NDim.\n");
   if ( (int)IDiffBin.size() != NDim ) printf("Error. you do not have the same number of IDiffBin than NDim.\n");
 
   if ( NDim == 1 ){
      for(int i=0;i<nBins1;i++){
	 nbins++;
	 bound[0] = bingrid1[i];
         LoBin.push_back(bound);
         bound[0] = bingrid1[i+1];
         UpBin.push_back(bound);

	 double binsize = IDiffBin[0] == 2 ? (UpBin.back())[0] - (LoBin.back())[0] : 1;
         BinSize.push_back(binsize);
      }
      // here we always assume, that all dimensions are
      // 'binned' dimensions (and not 'differential'). We were using IDiffBin
      // to tag, if the publication was divided by this binwidth or not, so we have 
      // to set NOW IDiffBin = 2
      IDiffBin[0] = 2 ;
   }
   else if ( NDim == 2 || NDim == 3 ){
      for(int i=0;i<nBins1;i++){
	 for(int j=0;j<nBins2[i];j++){
	    nbins ++;
	    // warning: the variables are exchanged here!
	    // what is bound[0] corresponds to bingrid2[nBins][nBins2]
	    // what is bound[1] corresponds to bingrid1[nBins]
	    bound[0] = bingrid2[i][j];
	    bound[1] = bingrid1[i];
	    //if ( NDim == 3 ) bound[2] = 0;
	    if ( NDim == 3 ) bound[2] = binwidth3;
	    LoBin.push_back(bound);
	    bound[0] = bingrid2[i][j+1];
	    bound[1] = bingrid1[i+1];
	    UpBin.push_back(bound);
	    //if ( NDim == 3 ) bound[2] = binwidth3;
	    if ( NDim == 3 ) bound[2] = 0;
	    //if ( binwidth3 != 0 ) bound[2] = binwidth3;
	    
	    double binsize = 1;
	
	    // warning: the variables are exchanged here!
	    // what is DimLabel[0] corresponds to bingrid2[nBins][nBins2]
	    // what is DimLabel[1] corresponds to bingrid1[nBins]
	    if ( IDiffBin[0] == 2 ) binsize *= bingrid2[i][j+1] - bingrid2[i][j];
	    if ( IDiffBin[1] == 2 ) binsize *= bingrid1[i+1] - bingrid1[i];
	    if ( NDim==3 ) { 
	       //if (IDiffBin[2] == 2 ) 
		  binsize *= binwidth3;
	    }
	    else if ( binwidth3 != 0 && NDim != 3 ) {
	       binsize *= binwidth3;
	    }
	    BinSize.push_back(binsize);
	    
	 }
      }
      // here we always assume, that all dimensions are
      // 'binned' dimensions (and not 'differential'). We were using IDiffBin
      // to tag, if the publication was divided by this binwidth or not, so we have 
      // to set NOW IDiffBin = 2
      // The 'third' dimension in this method however, is NOT a binned distribution
      IDiffBin[0] = 2 ;
      IDiffBin[1] = 2 ;
      if ( NDim==3 )
	 IDiffBin[2] = 1 ;
   }
   else error["InitBinning"]<<"Unknown NDim."<<endl;

   info["InitBinning"]<<"Tot. No. observable bins = "<<nbins<<endl;

   NObsBin = nbins;
 
   INormFlag = 0;    // --- fastNLO user: default=0 - set =1 if observable is 
			 //     to be normalized by own integral (in 1st dimension)
			 //     see documentation for details and for other options
   
}


void fastNLOTable::InitBinningKR( const int nBins1, const double* bingrid1, const int* nBins2, vector< vector<double> > bingrid2, const double bwfactor ){
  
  // ------------------------------------------------------------------- //
  //
  //  InitBinning. This method initalizes and (partly) checks the binning
  //  that is used in the scenario.
  //
  //   We must know NDim before calling InitBinning.
  //	NDim tells us, in how many dimensions/variables the measurement was performed
  //   this method only supports two dimensional measurements and a third dimension for 
  //   a pseudo-dimensional binning (if only one bin was measured e.g. in the pseudorapidity).
  //   Still, fastNLO could support higher dimensional binnings.
  // 
  //  input.
  //     nBins1	number of bins in 1st dimension
  //     bingrid1	binning in 1st dimension
  //     nBins2	number of bins of second dimension for each 1st-dimension variable
  //     bingrid	binning in 2nd dimension for each 1st dimension bin
  //     bwfactor	additional factor to take into account e.g. factors of 2 for binning in abs. rapidity
  //
  //  output.
  //     no output
  //
  //  initalized values
  //     LoBin, UpBin, BinSize, NObsBin
  //
  //
  //  logic for binwidth
  //     using IDiffBin you can specify, if the publ cross section was divided by this bin.
  //     if IDiffBin==2, then you divide the 'final' fastNLO-cross section by this number too -> we multiply the binwidth by this 
  //        dimensional width.
  //     if IDiffBin==1. then the cross section is not divided by this dimensional-binning width. However, we store
  //        the bingrid since the 'ObsBin' is binned in this binning.
  //        MENTION: the UpBin is the NOT stored in the table! (DB. maybe this could be changed.)
  //
  // ------------------------------------------------------------------- //


  if ( NDim == 2 && ( nBins2==NULL || bingrid2.empty()) ) printf("fastNLOTable::InitBinning. Error. NDim=2, but you have defined only one dimensional bin grid\n");
  if ( NDim != 1 && NDim != 2 )  printf("fastNLOTable::InitBinning. Error. NDim must be 1 or 2. three is not yet fully implemented.\n");

  vector <double> bound(NDim);
  double binsize = 0.;
  int nbins = 0;   // --- count total No. bins

  if ( (int)DimLabel.size() != NDim ) printf("Error. you do not have the same number of DimLabel than NDim.\n");
  if ( (int)IDiffBin.size() != NDim ) printf("Error. you do not have the same number of IDiffBin than NDim.\n");
 
  if ( NDim == 1 ){
    for(int i=0;i<nBins1;i++){
      nbins++;
      bound[0] = bingrid1[i];
      LoBin.push_back(bound);
      bound[0] = bingrid1[i+1];
      UpBin.push_back(bound);

      if ( bwfactor > 0. ) {
	binsize = bwfactor;
      } else {
	binsize = 1;
      }
      if ( IDiffBin[0] == 2 ) binsize *=  ((UpBin.back())[0] - (LoBin.back())[0]);
      BinSize.push_back(binsize);
    }
  } else if ( NDim == 2 ){
    for(int i=0;i<nBins1;i++){
      for(int j=0;j<nBins2[i];j++){
	nbins ++;
	// warning: the variables are exchanged here!
	// what is bound[0] corresponds to bingrid2[nBins][nBins2]
	// what is bound[1] corresponds to bingrid1[nBins]
	bound[0] = bingrid2[i][j];
	bound[1] = bingrid1[i];
	LoBin.push_back(bound);
	bound[0] = bingrid2[i][j+1];
	bound[1] = bingrid1[i+1];
	UpBin.push_back(bound);
	    
	if ( bwfactor > 0. ) {
	  binsize = bwfactor;
	} else {
	  binsize = 1;
	}
	// warning: the variables are exchanged here!
	// what is DimLabel[0] corresponds to bingrid2[nBins][nBins2]
	// what is DimLabel[1] corresponds to bingrid1[nBins]
	if ( IDiffBin[0] == 2 ) binsize *= bingrid2[i][j+1] - bingrid2[i][j];
	if ( IDiffBin[1] == 2 ) binsize *= bingrid1[i+1] - bingrid1[i];
	BinSize.push_back(binsize);
      }
    }
  }
  else error["InitBinningKR"]<<"Error. unknown NDim."<<endl;;
  
  info["InitBinningKR"]<<" tot. No. observable bins = "<<nbins<<endl;

  NObsBin = nbins;
 
  INormFlag = 0;    // --- fastNLO user: default=0 - set =1 if observable is 
  //     to be normalized by own integral (in 1st dimension)
  //     see documentation for details and for other options
   
}
*/

void fastNLOTable::Print() const {
   fastNLOBase::Print();
   PrintScenario();
}

void fastNLOTable::PrintScenario() const {
  printf("\n **************** FastNLO Table: Scenario ****************\n\n");
  printf("    Ipublunits                    %d\n",Ipublunits);
  for(unsigned int i=0;i<ScDescript.size();i++){
    printf("    ScDescript[%d]                 %s\n",i,ScDescript[i].data());
  }
  printf("    Ecms                          %7.4f\n",Ecms);
  printf("    ILOord                        %d\n",ILOord);
  printf("    NDim                          %d\n",NDim);
  for(int i=0;i<NDim;i++){
    printf("     - DimLabel[%d]                %s\n",i,DimLabel[i].data());
  }
  for(int i=0;i<NDim;i++){
    printf("     - IDiffBin[%d]               %d\n",i,IDiffBin[i]);
  }
  printf("    NObsBin                       %d\n",NObsBin);
  for(int i=0;i<NObsBin;i++){
    for(int j=0;j<NDim;j++){
      printf("     -  - LoBin[%d][%d]             %7.4f\n", i,j,Bin[i][j].first);
      if(IDiffBin[j]==2) 
	printf("     -  - UpBin[%d][%d]             %7.4f\n", i,j,Bin[i][j].second);
    }
   }
   for(int i=0;i<NObsBin;i++){
     printf("     - BinSize[%d]                %7.4f\n", i,BinSize[i]);
   }
   printf("    INormFlag                     %d\n",INormFlag);
   
   if(INormFlag>1){
     printf("    DenomTable                    %s\n",DenomTable.data());
   }
   if(INormFlag>0){
      for(int i=0;i<NObsBin;i++){
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

   //---  Initialization for nice printing
   const string CSEPS = "##################################################################################\n";
   const string LSEPS = "#---------------------------------------------------------------------------------\n";

   //
   // Print basic scenario information (always)
   //
   printf("\n");
   printf(" %s",CSEPS.c_str());
   printf(" # Information on fastNLO scenario: %s\n",ScenName.data());
   printf(" %s",LSEPS.c_str());
   printf(" # Description:\n");
   for (unsigned int i=0; i<ScDescript.size(); i++) {
      printf(" #   %s\n",ScDescript[i].data());
   }
   printf(" #\n");
   printf(" # Centre-of-mass energy Ecms: % -#10.4g GeV\n",Ecms);
   printf(" #\n");
   printf(" # Tot. no. of observable bins: %3i in %1i dimensions:\n",NObsBin,NDim);
   printf(" #\n");
   printf(" # No. of contributions: %1i\n",Ncontrib);
   for (unsigned int j = 0 ; j<fCoeff.size() ; j++) {
      if ( iprint == 0 ) fCoeff[j]->fastNLOCoeffBase::Print();
      else               fCoeff[j]->Print();
   }
   if (iprint > 0) {
      Print();
   }
   printf(" #\n");
   printf(" %s",CSEPS.c_str());
}


//______________________________________________________________________________
void fastNLOTable::PrintTableInfo(const int iprint) const {
   //
   //  Print basic info about fastNLO table and its contributions
   //   - iprint: iprint > 0: print also contribution descriptions
   //

   //---  Initialization for nice printing
   const string CSEPS = "##################################################################################\n";
   const string LSEPS = "#---------------------------------------------------------------------------------\n";
   printf("\n");
   printf(" %s",CSEPS.c_str());
   printf(" # Overview on contribution types and numbers contained in table:\n");
   printf(" %s",LSEPS.c_str());
   printf(" # Number of contributions: %2i\n",Ncontrib);


   for (unsigned int j = 0 ; j<fCoeff.size() ; j++) {
      string coeffname = fastNLO::_ContrName[fCoeff[j]->GetIContrFlag1()-1];
      fastNLOCoeffBase* c = fCoeff[j];
      cout << " # "<< "  No.: " << j << ", type: " << coeffname <<", Id: " << c->GetIContrFlag1()-1 
	   << ", order: " << c->GetContributionDescription()[0] 
	   << ", by: " << c->GetCodeDescription()[0] << endl;
      if (iprint > 0) {
	 for (unsigned int k = 0 ; k<c->GetCodeDescription().size(); k++) {
	    printf(" # \t\t%s\n",c->GetCodeDescription()[k].c_str());
	 }
      }
   }
   printf(" %s",CSEPS.c_str());

}
