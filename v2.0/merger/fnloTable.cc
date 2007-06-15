#include "fnloTable.h"

int fnloTable::OpenFileRead(){
   ifilestream = new ifstream(filename.c_str(),ios::in);
}
int fnloTable::OpenFileWrite(){
   if (access(filename.c_str(), F_OK) == 0){
      printf("fnloTable::OpenFileWrite: File for writing the table exists: %s.\nPlease remove it.\n",filename.c_str());
      exit(2);
   }
   ofilestream = new ofstream(filename.c_str(),ios::out);
   if(!ofilestream->good()){
       printf("fnloTable::OpenFileWrite: Cannot open %s for writing. Stopping.\n",filename.c_str());
       exit(2);
    }

}

fnloTable::~fnloTable(){
   if(ifilestream) delete ifilestream;
   if(ofilestream) delete ofilestream;
}
