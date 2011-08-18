#include "fnloBlockA1.h"
#include <iostream>

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
   *table >> NuserString;
   *table >> NuserInt;
   *table >> NuserFloat;
   *table >> Imachine;
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
   *table << NuserString << endl;
   *table << NuserInt << endl;
   *table << NuserFloat << endl;
   *table << Imachine << endl;
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


void fnloBlockA1::Print(){
  printf("\n **************** FastNLO Table: BlockA1 ****************\n\n");
  printf(" A1  tablemagicno                  %d\n",tablemagicno);
  printf(" A1  Itabversion                   %d\n",Itabversion);
  printf(" A1  ScenName                       %s\n",ScenName.data());
  printf(" A1  Ncontrib                      %d\n",Ncontrib);
  printf(" A1  Nmult                         %d\n",Nmult);
  printf(" A1  Ndata                         %d\n",Ndata);
  printf(" A1  NuserString                   %d\n",NuserString);
  printf(" A1  NuserInt                      %d\n",NuserInt);
  printf(" A1  NuserFloat                    %d\n",NuserFloat);
  printf(" A1  Imachine                      %d\n",Imachine);
  printf("\n ********************************************************\n\n");
}
