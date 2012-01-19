// KR: Add include because of header clean-up in gcc-4.3
#include <cstdlib>

#include "stdio.h"
#include <vector>


#include "fnloTable.h"
#include "entry.h"

using namespace std;

// KR: Change void** to char** for gcc-4.3
// int main(int argc,void** argv)
int main(int argc, char** argv)
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
        if(thistable->GetBlockA1()->GetItabversion() != 20000){
           printf("This Merger likes only tableformat V20.000, skipping %s.\n",path);      
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
  printf("Found %d table file(s).\n",ntables);
  if(ntables<1) exit(1);

  Entry *oneentry;
  vector <Entry*> entries;
  Contrib contribution;
  tableptr tablepointer;

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

  // extract the contribution types from all the block Bs of the tablefiles
  for( table=table_list.begin(); table!=table_list.end();  table++){
     fnloBlockA1 *blocka1 = (*table)->GetBlockA1();  
     int nblocks = blocka1->GetNcontrib()+blocka1->GetNdata();
     (*table)->RewindRead();
     (*table)->SkipBlockA1A2();
     for(int i=0;i<nblocks;i++){
        if((*table)->ReadBlockB(i)>0){
           printf("main.cc: Bad format in a B block found, file %s. Stopping.\n",(*table)->GetFilename().c_str());
           exit(1);
        }
        tablepointer.table = table;
        tablepointer.blockbno = i;
        fnloBlockB *blockb = (*table)->GetBlockB(i);
        contribution.IRef = blockb->GetIRef();
        contribution.IDataFlag = blockb->GetIDataFlag();
        contribution.IAddMultFlag = blockb->GetIAddMultFlag();
        contribution.IContrFlag1 = blockb->GetIContrFlag1();
        contribution.IContrFlag2 = blockb->GetIContrFlag2();
        contribution.NScaleDep = blockb->GetNScaleDep();
        contribution.Npow = blockb->GetNpow();
        bool newentry = true;
        for(int j=0;j<entries.size();j++){
           if(entries[j]->contribution==contribution){
              entries[j]->tables.push_back(tablepointer);
              newentry = false;
              break;
           }
        }
        if(newentry){
           oneentry = new Entry();
           oneentry->contribution = contribution;
           oneentry->tables.clear();
           oneentry->tables.push_back(tablepointer);
           entries.push_back(oneentry);
        }
     }
  }

  printf("No of types of contributions: %zu\n",entries.size());
  int Ncontrib = 0;
  int Nmult = 0;
  int Ndata = 0;
  for(int i=0;i<entries.size();i++){
     Nmult += entries[i]->contribution.IAddMultFlag;
     Ndata  += entries[i]->contribution.IDataFlag;
     Ncontrib += entries[i]->contribution.IAddMultFlag + (entries[i]->contribution.IContrFlag1>0?1:0);
  }

  //Write result
  fnloTable* result = table_list[0];
  const char* path=(char *)argv[nfiles];
  printf("Write merged results to file %s.\n",path);
  result->SetFilename(path);
  ofstream *outstream = result->OpenFileWrite();
  result->GetBlockA1()->SetNcontrib(Ncontrib);
  result->GetBlockA1()->SetNmult(Nmult);
  result->GetBlockA1()->SetNdata(Ndata);
  result->WriteBlockA1();
  result->WriteBlockA2();
  int nblocks = entries.size();
  for(int i=0;i<nblocks;i++){
     printf(" %zu file(s) containing",entries[i]->tables.size());
     printf(" %s %s",entries[i]->contribution.GetName1().c_str(),entries[i]->contribution.GetName2(table_list[0]->GetBlockA2()->GetILOord()).c_str());
     vector <fnloTable*>::iterator table = entries[i]->tables[0].table;
     // Addition 
     for(int j=1;j<entries[i]->tables.size();j++){
        fnloBlockB *otherblock = (*(entries[i]->tables[j].table))->GetBlockB(entries[i]->tables[j].blockbno);
        (*table)->GetBlockB(entries[i]->tables[0].blockbno)->Add(otherblock);
     }
     (*table)->WriteBlockB(entries[i]->tables[0].blockbno,outstream);
     long long int nevents = (*table)->GetBlockB(entries[i]->tables[0].blockbno)->GetNevt();
     printf("  (%#4.2g events).\n",(double)nevents);
     
  }
  result->CloseFileWrite();
  printf("Done.\n");

  return 0;
}
