#include "fnloBlockB.h"


int fnloBlockB::Read(istream *table){
   table->peek();
   if (table->eof()){
      printf("fnloBlockB::Read: Cannot read from file.\n");
      return(2);
   }

   int key = 0;
   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockB::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
      return 1;
   };







   *table >> key;
   if(key != tablemagicno){
      printf("fnloBlockB::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
      return 1;
   };
   // Put magic number back
   for(int i=0;i<(int)(log10((double)key)+1);i++){
      table->unget();
   }
   return 0;
}

int fnloBlockB::Write(ostream *table){
   *table << tablemagicno << endl;
   return 0;
}



bool fnloBlockB::IsCompatible(fnloBlockB* other){
   return true;
};
