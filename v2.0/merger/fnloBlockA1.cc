#include "fnloBlockA1.h"


int fnloBlockA1::Read(istream *table){
   table->peek();
   if (table->eof()){
      printf("fnloBlockA1::Read: Cannot read from file.\n");
      return(2);
   }

   int key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockA1::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
      return 1;
   };
   *table >> Itabversion;
   *table >> ScenName;
   *table >> Ncontrib;
   *table >> Nmult;
   *table >> Ndata;
   key=0;
   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockA1::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
      return 1;
   };
   // Put magic number back
   for(int i=0;i<(int)(log10((double)key)+1);i++){
      table->unget();
   }
   return 0;
}

int fnloBlockA1::Write(ostream *table){
   *table << tablemagicno << endl;
   *table << Itabversion << endl;
   *table << ScenName << endl;
   *table << Ncontrib << endl;
   *table << Nmult << endl;
   *table << Ndata << endl;
   return 0;
}



bool fnloBlockA1::IsCompatible(fnloBlockA1* other){
   if(Itabversion!= other->GetItabversion()){
      printf("fnloBlockA1::IsCompatible: Differing Versions of table format: %d and %d\n",Itabversion,other->GetItabversion());
      return false;
   }
   if(ScenName!= other->GetScenName()){
      printf("fnloBlockA1::IsCompatible: Differing Names of Scenarios: %s and %s\n",ScenName.c_str(),other->ScenName.c_str());
      return false;
   }
   if(Ndata + other->GetNdata() > 1){
      printf("fnloBlockA1::IsCompatible: Two tables containing both experimental data are incompatible\n");
      return false;
   }
   
   return true;
};
