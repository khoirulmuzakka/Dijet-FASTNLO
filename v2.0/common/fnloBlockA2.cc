#include "fnloBlockA2.h"


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
   for(int i=0;i<NObsBin;i++){
      LoBin[i].resize(NDim);
      UpBin[i].resize(NDim);
      for(int j=0;j<NDim;j++){
         *table >>  LoBin[i][j];
         if(IDiffBin[j]==2) *table >>  UpBin[i][j];
      }
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
