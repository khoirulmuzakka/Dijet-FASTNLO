#ifndef __fnloTable__
#define __fnloTable__

#include <fstream>

#include "fnloBlockA1.h"
#include "fnloBlockA2.h"
#include "fnloBlockB.h"

class fnloTable{
public:
   fnloTable(){filename = "";BlockIndexLO=BlockIndexNLO=BlockIndexLORef=BlockIndexNLORef=-1;}
   fnloTable(string name){filename = name;BlockIndexLO=BlockIndexNLO=BlockIndexLORef=BlockIndexNLORef=-1;}
   int ReadTable();
   int OpenFileRead();
   void RewindRead();
   void SkipBlockA1A2();
   ofstream *OpenFileWrite();
   ofstream *OpenFileRewrite();
   void CloseFileWrite();
   int ReadBlockA1(){return BlockA1.Read(ifilestream);}
   int WriteBlockA1(){return BlockA1.Write(ofilestream);}
   fnloBlockA1* GetBlockA1(){return &BlockA1;}
   int ReadBlockA2(){return BlockA2.Read(ifilestream);}
   int WriteBlockA2(){return BlockA2.Write(ofilestream);}
   fnloBlockA2* GetBlockA2(){return &BlockA2;}
   int CreateBlockB(int no);
   int CreateBlockB(int no,fnloBlockB *newblock);
   int ReadBlockB(int no);
   int WriteBlockB(int no);
   int WriteBlockB(int no,ofstream* outstream );
   int WriteBlockBDividebyN(int no);   
   fnloBlockB* GetBlockB(int no){return BlockB[no];}
   string GetFilename(){return filename;}
   void SetFilename(string name){filename=name;}

   int GetNObsBin(){return BlockA2.GetNObsBin();}
   int GetNScaleVar(int dimension){return BlockB[0]->Nscalevar[dimension];}
   double GetScaleFac(int dimension,int scalevar){return BlockB[0]->ScaleFac[dimension][scalevar];}
   string GetScaleDescript(int dimension,int line){return BlockB[0]->ScaleDescript[dimension][line];}
   double GetLoBin(int bin, int dimension){return BlockA2.LoBin[bin][dimension];}
   double GetUpBin(int bin, int dimension){return BlockA2.UpBin[bin][dimension];}

   ~fnloTable();
private:
   string filename;
   ifstream *ifilestream;
   ofstream *ofilestream;
   fnloBlockA1 BlockA1;
   fnloBlockA2 BlockA2;
   vector < fnloBlockB* > BlockB;
protected:
   int BlockIndexLO;
   int BlockIndexNLO;
   int BlockIndexLORef;
   int BlockIndexNLORef;

};
#endif
