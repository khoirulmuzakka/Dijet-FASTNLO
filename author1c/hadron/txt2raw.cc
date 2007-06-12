#include <stdio.h>
#include <unistd.h>             
#include <iostream>
#include <fstream>

#define WRITE(n) outtable.write(reinterpret_cast<char *>(&n),sizeof(n))

using namespace std;

int main(int argc,void** argv)
{
  // check parameters
  if(argc!=3){
    printf("Usage: txt2raw table.txt table.raw \n");
    return(1);
  }
  const char* inpath=(char *)argv[1];
  // File there?
  if (access(inpath, R_OK) != 0){
     printf("Cannot access %s.\n",inpath);      
     return(2);
  }
  fstream intable(inpath,ios::in); // open file

  const char* outpath=(char *)argv[2];
  // check if result file is writable
  if (access(outpath, F_OK) == 0){
     printf("File for writing the table exists: %s.\n Please remove it first.\n",outpath);      
     return(2);
  }  
  fstream outtable(outpath,ios::out|ios::binary); // open file for writing

  char buffer[256];
  double doublebuffer;
  int intbuffer;
  
  while(!intable.eof()){
     intable.getline(buffer,256);
     //is it a double?
     if (strstr(buffer,".")){
        if(sscanf(buffer,"%lg",&doublebuffer)){
           WRITE(doublebuffer) ;// it is a double!
        }else{
           printf("Warning: '%s' is treated as a string.\n",buffer); // strange, contains a dot but is no double
           outtable << buffer  <<endl;
        }
     }else{ // next test for integer?
        if(sscanf(buffer,"%d",&intbuffer)){
           WRITE(intbuffer) ;// it is a int!
        }else{
           outtable << buffer  <<endl; // last resort: a string
        }
     }
     
  }
  
  return(0);
}
