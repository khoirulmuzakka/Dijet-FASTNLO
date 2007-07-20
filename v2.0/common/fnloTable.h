#ifndef __fnloTable__
#define __fnloTable__

#include <fstream>

#include "fnloBlockA1.h"
#include "fnloBlockA2.h"
#include "fnloBlockB.h"

class fnloTable{
 public:
   fnloTable(string name){filename = name;}
   int OpenFileRead();
   void RewindRead();
   void SkipBlockA1A2();
   ofstream *OpenFileWrite();
   void CloseFileWrite();
   int ReadBlockA1(){return BlockA1.Read(ifilestream);}
   int WriteBlockA1(){return BlockA1.Write(ofilestream);}
   fnloBlockA1* GetBlockA1(){return &BlockA1;}
   int ReadBlockA2(){return BlockA2.Read(ifilestream);}
   int WriteBlockA2(){return BlockA2.Write(ofilestream);}
   fnloBlockA2* GetBlockA2(){return &BlockA2;}
   int CreateBlockB(int no);
   int ReadBlockB(int no);
   int WriteBlockB(int no);
   int WriteBlockB(int no,ofstream* outstream );
   fnloBlockB* GetBlockB(int no){return BlockB[no];}
   string GetFilename(){return filename;}
   void SetFilename(string name){filename=name;}
   
   ~fnloTable();
 private:
   string filename;
   ifstream *ifilestream;
   ofstream *ofilestream;
   fnloBlockA1 BlockA1;
   fnloBlockA2 BlockA2;
   vector < fnloBlockB* > BlockB;
};
#endif
