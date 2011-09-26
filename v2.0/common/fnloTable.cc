// KR: Add include because of header clean-up in gcc-4.3
#include <cstdlib>

#include "fnloTable.h"

int fnloTable::ReadTable(){
   OpenFileRead();
   ReadBlockA1();
   ReadBlockA2();
   fnloBlockA1 *blocka1 = GetBlockA1();  
   int nblocks = blocka1->GetNcontrib()+blocka1->GetNdata();
   for(int i=0;i<nblocks;i++){
      ReadBlockB(i);
      fnloBlockB *B = GetBlockB(i);
      if(B->IContrFlag1==1 && B->IContrFlag2==1 && B->IRef==0) BlockIndexLO=i;
      if(B->IContrFlag1==1 && B->IContrFlag2==2 && B->IRef==0) BlockIndexNLO=i;
      if(B->IContrFlag1==1 && B->IContrFlag2==1 && B->IRef==1) BlockIndexLORef=i;
      if(B->IContrFlag1==1 && B->IContrFlag2==2 && B->IRef==1) BlockIndexNLORef=i;
   }
   return 0;
}

int fnloTable::OpenFileRead(){
   ifilestream = new ifstream(filename.c_str(),ios::in);
   return 0;
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

void fnloTable::WriteTable ( long long int nevents ){
   for (int k=0;k<this->GetBlockA1()->GetNcontrib();k++){
      this->GetBlockB(k)->Nevt = (long long int)nevents;
   }
   this->OpenFileRewrite();
   this->WriteBlockA1();
   this->WriteBlockA2();
   for(int i=0;i< this->GetBlockA1()->GetNcontrib();i++){
      this->WriteBlockBDividebyN(i);
   }
   this->CloseFileWrite();
}

int fnloTable::GetBinNumber( double var1, double var2 ){
   //
   //  calculate the bin number as define in A2::LoBin and A2::UpBin
   //  initialized by A2::InitBinning.
   //
   //  returns the bin number, that has to be passse to FillEvent()
   //  return -1 if values are out of bin-ranges
   // 

   fnloBlockA2* A2 =  this->GetBlockA2();
   if ( var2==0 && A2->NDim!=1){
      printf("fnloTable::GetBinNumber(%6.3f,%6.3f). Error. A single differential table only has one variable.\n",var1,var2);exit(1);
   }
   if ( var2!=0 && A2->NDim!=2){
      printf("fnloTable::GetBinNumber(%6.3f,%6.3f). Error. A double differential table only needs two variables.\n",var1,var2);exit(1);
   }
   
   int bin = -1;
   if ( A2->NDim == 2) {
	for(int j = 0; j < A2->GetNObsBin(); j++) {
	   if ( var1 >= A2->LoBin[j][0]  && var1 <  A2->UpBin[j][0] &&
		var2 >= A2->LoBin[j][1]  && var2 <  A2->UpBin[j][1]) {
	      bin=j;
	      break;
	   }
	}
   }
   else if ( A2->NDim == 1) {
      for(int j = 0; j < A2->GetNObsBin(); j++) {
	 if ( var1 >= A2->LoBin[j][0]  && var1 <  A2->UpBin[j][0] ){
	    bin=j;
	    break;
	 }
      }
   }
   else {
      printf("fnloTable::GetBinNumber(%6.3f,%6.3f). Error. Only single, or double differential tables are supported. A2::NDim = %d.\n",var1,var2,A2->NDim);exit(1);
   }

   return bin;
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
      printf("fnloTable::ReadBlockB: Cannot create block B #%d. Already existing blocks: %zu. Stopping.\n",no,BlockB.size());
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
      printf("fnloTable::WriteBlockB: Table no. %d does not exist, only up to %zu. Stopping.\n",no,BlockB.size());
      exit(2);
      
   }
   return blockb->Write(ofilestream);
}

int fnloTable::WriteBlockBDividebyN(int no){
   fnloBlockB* blockb;
   if((no)<BlockB.size()){
      blockb = BlockB[no];
   }else{
      printf("fnloTable::WriteBlockB: Table no. %d does not exist, only up to %zu. Stopping.\n",no,BlockB.size());
      exit(2);
      
   }
   return blockb->Write(ofilestream,fnloBlockB::DividebyNevt);
}


int fnloTable::WriteBlockB(int no,ofstream* outstream ){
   fnloBlockB* blockb;
   if((no)<BlockB.size()){
      blockb = BlockB[no];
   }else{
      printf("fnloTable::WriteBlockB: Table no. %d does not exist, only up to %zu. Stopping.\n",no,BlockB.size());
      exit(2);
      
   }
   return blockb->Write(outstream);
}


void fnloTable::DeleteAllBlockB()
{
	for (size_t i = 0; i < BlockB.size(); ++i)
		delete BlockB[i];
	BlockB.clear();
}

fnloTable::~fnloTable(){
   if(ifilestream) delete ifilestream;
   if(ofilestream) delete ofilestream;
}


int fnloTable::GetScale1index(int scalevar){
   fnloBlockB* blockb;
   blockb = BlockB[0];
   int scalei1 = scalevar;
   if(blockb->NScaleDim==2){
      scalei1 = scalevar / blockb->Nscalevar[1];
   }
   if(blockb->NScaleDim>2){
      printf("fnloTable::GetScale1index: more than two scale dimensions not yet implemented. Stopping.\n");
      exit(2);
   }
   return scalei1;
}

int fnloTable::GetScale2index(int scalevar){
   fnloBlockB* blockb;
   blockb = BlockB[0];
   int scalei2 = scalevar;
   if(blockb->NScaleDim==2){
      scalei2 = scalevar % blockb->Nscalevar[1];
   }
   if(blockb->NScaleDim>2){
      printf("fnloTable::GetScale1index: more than two scale dimensions not yet implemented. Stopping.\n");
      exit(2);
   }
   return scalei2;
}
