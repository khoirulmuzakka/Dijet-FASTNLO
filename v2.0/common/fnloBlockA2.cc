#include "fnloBlockA2.h"
//#include <iostream>

int fnloBlockA2::Read(istream *table){
   table->peek();
   if (table->eof()){
      printf("fnloBlockA2::Read: Cannot read from file.\n");
      return(2);
   }

   int key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockA2::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
      return 1;
   };

   *table >> Ipublunits;
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
   LoBin.resize(NObsBin);
   UpBin.resize(NObsBin);
   //KR: Set rapidity index also when reading a table
   RapIndex.push_back(0);
   //   int irap = 0;
   for(int i=0;i<NObsBin;i++){
      LoBin[i].resize(NDim);
      UpBin[i].resize(NDim);
      for(int j=0;j<NDim;j++){
         *table >>  LoBin[i][j];
         if(IDiffBin[j]==2) *table >>  UpBin[i][j];
      }
      //      cout << "iobs1: " << i << ", LoBin i: " << LoBin[i][1] << endl;
      if ( i > 0 ) {
	if ( LoBin[i][1] != LoBin[i-1][1] ) {
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

   key=0;
   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockA2::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
      return 1;
   };
   // Put magic number back
   for(int i=0;i<(int)(log10((double)key)+1);i++){
      table->unget();
   }
   return 0;
}

int fnloBlockA2::Write(ostream *table){
   *table << tablemagicno << endl;
   *table << Ipublunits << endl;
   NScDescript	= ScDescript.size();
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
         *table <<  LoBin[i][j]  << endl;
         if(IDiffBin[j]==2) *table <<  UpBin[i][j]  << endl;
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


bool fnloBlockA2::IsCompatible(fnloBlockA2* other){
   if(Ipublunits != other->Ipublunits){
      printf("fnloBlockA2::IsCompatible: Differing cross section units found: %d and %d\n",Ipublunits,other->Ipublunits);
      return false;
   }
   if(NScDescript != other->NScDescript){
      printf("fnloBlockA2::IsCompatible: Differing NScDescript found: %d and %d\n",NScDescript,other->NScDescript);
      return false;
   }  
   if(ScDescript != other->ScDescript){
      printf("fnloBlockA2::IsCompatible: Differing ScDescript found.\n");
      return false;
   }  
   if(!cmp(Ecms,other->Ecms)){
      printf("fnloBlockA2::IsCompatible: Differing Ecms found: %f and %f.\n",Ecms,other->Ecms);
      return false;
   }  
   if(ILOord != other->ILOord){
      printf("fnloBlockA2::IsCompatible: Differing ILOord found: %d and %d\n",ILOord,other->ILOord);
      return false;
   }  
   if(NObsBin != other->NObsBin){
      printf("fnloBlockA2::IsCompatible: Differing NObsBin found: %d and %d\n",NObsBin,other->NObsBin);
      return false;
   }  
   if(NDim != other->NDim){
      printf("fnloBlockA2::IsCompatible: Differing NDim found: %d and %d\n",NDim,other->NDim);
      return false;
   }  
   if(DimLabel != other->DimLabel){
      printf("fnloBlockA2::IsCompatible: Differing DimLabel found.\n");
      return false;
   }  
   if(IDiffBin != other->IDiffBin){
      printf("fnloBlockA2::IsCompatible: Differing IDiffBin found.\n");
      return false;
   }  
   if(!cmp(LoBin,other->LoBin)){
      printf("fnloBlockA2::IsCompatible: Differing LoBin found.\n");
      return false;
   }
   if(!cmp(UpBin,other->UpBin)){
      printf("fnloBlockA2::IsCompatible: Differing UpBin found.\n");
      return false;
   }
   if(!cmp(BinSize,other->BinSize)){
      printf("fnloBlockA2::IsCompatible: Differing BinSize found.\n");
      return false;
   }
   if(INormFlag != other->INormFlag){
      printf("fnloBlockA2::IsCompatible: Differing INormFlag found: %d and %d\n",INormFlag,other->INormFlag);
      return false;
   }  
   if(INormFlag>1){
      if(DenomTable != other->DenomTable){
         printf("fnloBlockA2::IsCompatible: Differing DenomTable found.\n");
         return false;
      }
   }
   if(INormFlag>0){
      for(int i=0;i<NObsBin;i++){
         if(IDivLoPointer[i] != other->IDivLoPointer[i]){
            printf("fnloBlockA2::IsCompatible: Differing IDivLoPointer[%d] found.\n",i);
            return false;
         }
         if(IDivUpPointer[i] != other->IDivUpPointer[i]){
            printf("fnloBlockA2::IsCompatible: Differing IDivUpPointer[%d] found.\n",i);
            return false;
         }
      }
   }  
   
   return true;
};

void fnloBlockA2::StripWhitespace(string &str) const{
   for(string::iterator achar = str.end(); achar>str.begin();achar--) {
      if (*achar==0x20 || *achar==0x00){
         str.erase(achar);
      }else{
         break;
      }
   }
}

bool fnloBlockA2::cmp(const double x1, const double x2) const{
   double norm;
   if (x1>0.){
      norm = x1;
   }else{
      norm = 1.; // If x1 is 0, do not try to calculate relative deviation, use absolute
   }
   return((fabs(x1-x2)/norm)<1e-7);
}

bool fnloBlockA2::cmp(const vector < double > x1,const vector < double > x2) const{
   bool result = true;
   for(int i = 0; i<x1.size() ;i++ ){
      result = result & cmp (x1[i],x2[i]);
   }
   return result;
}

bool fnloBlockA2::cmp(vector < vector < double > > x1,  vector < vector < double > > x2) const{
   bool result = true;
   for(int i = 0; i<x1.size() ;i++ ){
      result = result & cmp (x1[i],x2[i]);
   }
   return result;
}



void fnloBlockA2::InitBinning( const int nBins1 , double* bingrid1 , const int* nBins2  , vector<double*> bingrid2  , double binwidth3 ){

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
   //                   or you can use binwidth3 if your binning is e.g. in TeV, but you want to have pb/GeV
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
   //     if IDiffBin==2, then you divide the 'final' fnlo-cross section by this number too -> we multiply the binwidth by this 
   //        dimensional width.
   //     if IDiffBin==1. then the cross section is not divided by this dimensional-binning width. However, we store
   //        the bingrid since the 'ObsBin' is binned in this binning.
   //        MENTION: the UpBin is the NOT stored in the table! (DB. maybe this could be changed.)
   //
   // ------------------------------------------------------------------- //


   if ( NDim == 2 && ( nBins2==NULL || bingrid2.empty()) ) printf("fnloBlockA2::InitBinning. Error. NDim=2, but you have defined only one dimensional bin grid\n");
   if ( NDim == 3 ) printf("fnloBlockA2::InitBinning. Error. NDim=3, you must define the binwidth of dimension 3!\n");
   //if ( NDim != 1 && NDim != 2 )  printf("fnloBlockA2::InitBinning. Error. NDim must be 1 or 2. three is not yet fully implemented.\n");

   vector <double> bound(NDim);
   double binsize = 0.;
   int nbins = 0;   // --- count total No. bins

   if ( DimLabel.size() != NDim ) printf("Error. you do not have the same number of DimLabel than NDim.\n");
   if ( IDiffBin.size() != NDim ) printf("Error. you do not have the same number of IDiffBin than NDim.\n");
 
   if ( NDim == 1 ){
      for(int i=0;i<nBins1;i++){
	 nbins++;
	 bound[0] = bingrid1[i];
         LoBin.push_back(bound);
         bound[0] = bingrid1[i+1];
         UpBin.push_back(bound);

	 binsize = 1;
	 if ( IDiffBin[0] == 2 ) binsize *=  ((UpBin.back())[0] - (LoBin.back())[0]);
         BinSize.push_back(binsize);
      }
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
	    if ( NDim == 3 ) bound[2] = 0;
	    LoBin.push_back(bound);
	    bound[0] = bingrid2[i][j+1];
	    bound[1] = bingrid1[i+1];
	    UpBin.push_back(bound);
	    if ( NDim == 3 ) bound[2] = binwidth3;
	    //if ( binwidth3 != 0 ) bound[2] = binwidth3;
	    
	    binsize = 1;
	    // warning: the variables are exchanged here!
	    // what is DimLabel[0] corresponds to bingrid2[nBins][nBins2]
	    // what is DimLabel[1] corresponds to bingrid1[nBins]
	    if ( IDiffBin[0] == 2 ) binsize *= bingrid2[i][j+1] - bingrid2[i][j];
	    if ( IDiffBin[1] == 2 ) binsize *= bingrid1[i+1] - bingrid1[i];
	    if ( NDim==3 ) { 
	       if (IDiffBin[2] == 2 ) binsize *= binwidth3;
	    }
	    else if ( binwidth3 != 0 && NDim != 3 ) { // go from e.g. from GeV to TeV
	       binsize *= 1000;
	    }

	    BinSize.push_back(binsize);
	 }
      }
   }
   else printf("fnloBlockA2::InitBinning. Error. unknown NDim.\n");

   printf(" tot. No. observable bins = %d\n",nbins);

   NObsBin = nbins;
 

   INormFlag = 0;    // --- fastNLO user: default=0 - set =1 if observable is 
			 //     to be normalized by own integral (in 1st dimension)
			 //     see documentation for details and for other options
   
}


void fnloBlockA2::Print(){
  printf("\n **************** FastNLO Table: BlockA1 ****************\n\n");
  printf(" A2  Ipublunits                    %d\n",Ipublunits);
  //   NScDescript = NScDescript.size();
  //   printf(" A2  NScDescript                   %d\n",NScDescript);
  for(int i=0;i<ScDescript.size();i++){
    printf(" A2  ScDescript[%d]                 %s\n",i,ScDescript[i].data());
  }
  printf(" A2  Ecms                          %7.4f\n",Ecms);
  printf(" A2  ILOord                        %d\n",ILOord);
  printf(" A2  NDim                          %d\n",NDim);
  for(int i=0;i<NDim;i++){
    printf(" A2   - DimLabel[%d]                %s\n",i,DimLabel[i].data());
  }
  for(int i=0;i<NDim;i++){
    printf(" A2   - IDiffBin[%d]               %d\n",i,IDiffBin[i]);
  }
  printf(" A2  NObsBin                       %d\n",NObsBin);
  for(int i=0;i<NObsBin;i++){
    for(int j=0;j<NDim;j++){
      printf(" A2   -  - LoBin[%d][%d]             %7.4f\n", i,j,LoBin[i][j]);
      if(IDiffBin[j]==2) 
	printf(" A2   -  - UpBin[%d][%d]             %7.4f\n", i,j,UpBin[i][j]);
    }
   }
   for(int i=0;i<NObsBin;i++){
     printf(" A2   - BinSize[%d]                %7.4f\n", i,BinSize[i]);
   }
   printf(" A2  INormFlag                     %d\n",INormFlag);
   
   if(INormFlag>1){
     printf(" A2  DenomTable                    %s\n",DenomTable.data());
   }
   if(INormFlag>0){
      for(int i=0;i<NObsBin;i++){
        printf(" A2   - IDivLoPointer[%d]               %d\n",i,IDivLoPointer[i]);
	printf(" A2   - IDivUpPointer[%d]               %d\n",i,IDivUpPointer[i]);
      }
   }
   printf("\n ********************************************************\n\n");
  
}
