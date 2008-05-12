#include "fnloTable.h"

#ifdef __CINT__
   ClassImp(fnloTable)
#endif

int fnloTable::OpenFileRead(){
      //   if (ifilestream) delete ifilestream;
   ifilestream = new ifstream(filename.c_str(),ios::in);
}

void fnloTable::RewindRead(){
   ifilestream->close();
   ifilestream->open(filename.c_str(),ios::in);
}

void fnloTable::SkipBlockA1A2(){
   char buffer[257];
   int key;
   int count = 0;
  while(!ifilestream->eof()){
      ifilestream->getline(buffer,256);
      sscanf(buffer,"%d",&key);
      if(key == tablemagicno) count++;
      if(count == 3){
         // Put magic number back
         ifilestream->unget();
         for(int i=0;i<(int)(log10((double)key)+1);i++){
            ifilestream->unget();
         }
         break;
      }
   }
}

ofstream *fnloTable::OpenFileWrite(){
   if (access(filename.c_str(), F_OK) == 0){
      printf("fnloTable::OpenFileWrite: File for writing the table exists: %s.\nPlease remove it.\n",filename.c_str());
      exit(2);
   }
   ofilestream = new ofstream(filename.c_str(),ios::out);
   if(!ofilestream->good()){
       printf("fnloTable::OpenFileWrite: Cannot open %s for writing. Stopping.\n",filename.c_str());
       exit(2);
    }
   ofilestream->precision(8);
   return ofilestream;
}

ofstream *fnloTable::OpenFileRewrite(){
   ofilestream = new ofstream(filename.c_str(),ios::out);
   if(!ofilestream->good()){
       printf("fnloTable::OpenFileWrite: Cannot open %s for writing. Stopping.\n",filename.c_str());
       exit(2);
    }
   ofilestream->precision(8);
   return ofilestream;
}

void fnloTable::CloseFileWrite(){
   *ofilestream << tablemagicno << endl;
   *ofilestream << tablemagicno << endl;
   ofilestream->close();
}

int fnloTable::ReadBlockB(int no){
   if( CreateBlockB(no) ){
      printf("fnloTable::ReadBlockB: Cannot create block B #%d. Already existing blocks: %d. Stopping.\n",no,BlockB.size());
   }
   return BlockB[no]->Read(ifilestream);
}

int fnloTable::CreateBlockB(int no){
   fnloBlockB* blockb = new fnloBlockB(&BlockA1,&BlockA2);
   if((no+1)>BlockB.size()){
      BlockB.resize(no+1);
   }
   BlockB[no] = blockb;
   return 0;
}

int fnloTable::CreateBlockB(int no,fnloBlockB *newblockb){
   if((no+1)>BlockB.size()){
      BlockB.resize(no+1);
   }
   BlockB[no] = newblockb;
   return 0;
}

int fnloTable::WriteBlockB(int no){
   fnloBlockB* blockb;
   if((no)<BlockB.size()){
      blockb = BlockB[no];
   }else{
      printf("fnloTable::WriteBlockB: Table no. %d does not exist, only up to %d. Stopping.\n",no,BlockB.size());
      exit(2);
      
   }
   return blockb->Write(ofilestream);
}

int fnloTable::WriteBlockBDividebyN(int no){
   fnloBlockB* blockb;
   if((no)<BlockB.size()){
      blockb = BlockB[no];
   }else{
      printf("fnloTable::WriteBlockB: Table no. %d does not exist, only up to %d. Stopping.\n",no,BlockB.size());
      exit(2);
      
   }
   return blockb->Write(ofilestream,fnloBlockB::DividebyNevt);
}


int fnloTable::WriteBlockB(int no,ofstream* outstream ){
   fnloBlockB* blockb;
   if((no)<BlockB.size()){
      blockb = BlockB[no];
   }else{
      printf("fnloTable::WriteBlockB: Table no. %d does not exist, only up to %d. Stopping.\n",no,BlockB.size());
      exit(2);
      
   }
   return blockb->Write(outstream);
}



fnloTable::~fnloTable(){
   if(ifilestream) delete ifilestream;
   if(ofilestream) delete ofilestream;
}
