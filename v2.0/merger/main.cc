#include "stdio.h"
#include <vector>

#include "fnloTable.h"

using namespace std;

int main(int argc,void** argv)
{
  // check parameters
  if(argc<3){
    printf("Usage: fnlo-merge file1.txt [filex.txt]+ result.txt\n");
    return(1);
  }
  // loop over arguments
  int nfiles = argc-1;
  vector <fnloTable*> table_list;
  for(int i=0;i<nfiles-1;i++){
     const char* path=(char *)argv[i+1];
     // File there?
     if (access(path, R_OK) == 0){
        fnloTable* thistable = new fnloTable(path);
        table_list.push_back(thistable);
        thistable->OpenFileRead();
        if(thistable->ReadBlockA1() != 0){
           printf("Invalid format in %s, skipping.\n",path);      
           table_list.pop_back();
           delete thistable;
           break;
        }
        if(thistable->ReadBlockA2() != 0){
           printf("Invalid format in %s, skipping.\n",path);      
           table_list.pop_back();
           delete thistable;
           break;
        }
     }else{
        printf("Cannot access %s, skipping.\n",path);      
     }
  }
  
  int ntables = table_list.size();
  printf("Found %d table(s).\n",ntables);
  if(ntables<1) exit(1);

  if(ntables>1){
     // check for compatibility by looking at block A1 and A2
     vector <fnloTable*>::iterator table;

     fnloBlockA1 *blocka1 = table_list[0]->GetBlockA1();
     fnloBlockA2 *blocka2 = table_list[0]->GetBlockA2();
     for( table=table_list.begin()+1; table!=table_list.end();  table++){
        if(! blocka1->IsCompatible((*table)->GetBlockA1())){
           printf("main.cc: Non compatible A1 blocks found. Files: %s and %s. Stopping.\n",
                  table_list[0]->GetFilename().c_str(),(*table)->GetFilename().c_str());
           exit(1);
        }
        if(! blocka2->IsCompatible((*table)->GetBlockA2())){
           printf("main.cc: Non compatible A2 blocks found. Files: %s and %s. Stopping.\n",
                  table_list[0]->GetFilename().c_str(),(*table)->GetFilename().c_str());
           exit(1);
        }
     }

     // extract the contribution types from the tables
  }

  //Write result
  fnloTable* result = table_list[0];
  const char* path=(char *)argv[nfiles];
  printf("Write merged results to file %s.\n",path);
  result->SetFilename(path);
  result->OpenFileWrite();
  result->WriteBlockA1();
  result->WriteBlockA2();
  result->CloseFileWrite();
  printf("Done.\n");

  return 0;
}
