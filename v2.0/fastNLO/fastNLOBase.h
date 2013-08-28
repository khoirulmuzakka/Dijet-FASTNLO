#ifndef __fastNLOBase__
#define __fastNLOBase__

#include <fstream>
#include <iostream>
#include <istream>
#include <string>
#include "speaker.h"

using namespace std;

class fastNLOBase : public PrimalScream {

public:
   fastNLOBase();
   fastNLOBase(string name);
   ~fastNLOBase();

   // i/o
   int ReadTable();
   int ReadHeader(istream *table);
   
   int OpenFileRead();
   void WriteTable();
   int WriteHeader(ostream *table);
   void RewindRead();
   void SkipBlockA1A2();
   ofstream *OpenFileWrite();
   ofstream *OpenFileRewrite();
   void CloseFileWrite();
   void CloseStream();
   //    int WriteBlockB(int no);
//    int WriteBlockB(int no,ofstream* outstream );
//    int WriteBlockBDividebyN(int no);
//    void DeleteAllBlockB(); // FIX: Correct way would be to fix object ownership
//    int CreateBlockB(int no);
//    int CreateBlockB(int no,fastNLOBlockB *newblock);


   // useful getters
   //fastNLOBlockB* GetBlockB(int no){return BlockB[no];}
   //    int GetNObsBin(){return fScenario.GetNObsBin();}
   //    int GetNScaleVar(int dimension){return BlockB[0]->Nscalevar[dimension];}
   //    int GetNScaleVar(){return BlockB[0]->GetTotalScalevars();}
   //    int GetIscale(int scale){return BlockB[0]->Iscale[scale];}
   //    double GetScaleFac(int dimension,int scalevar){return BlockB[0]->ScaleFac[dimension][scalevar];}
   //    string GetScaleDescript(int dimension,int line){return BlockB[0]->ScaleDescript[dimension][line];}
   //    double GetLoBin(int bin, int dimension){return fScenario.LoBin[bin][dimension];}
   //    double GetUpBin(int bin, int dimension){return fScenario.UpBin[bin][dimension];}

   virtual void Print() const;
   
   // header
   void PrintHeader() const;
   void SetHeaderDefaults();
   void ResetHeader();
   void SetContributionHeader();
   bool IsCompatibleHeader(fastNLOBase* other) const;

   // getter/setters
   string GetFilename() const {return ffilename;}
   void   SetFilename(string name){ffilename=name;}

   int  GetItabversion() const {return Itabversion;}
   void SetItabversion(int version){Itabversion = version;}

   string GetScenName() const {return ScenName;}
   void   SetScenName(string name){ScenName = name;}

   int  GetNmult() const {return Nmult;}
   void SetNmult(int n){Nmult = n;}

   int  GetNuserString() const {return NuserString;}
   void SetNuserString(int n){NuserString = n;}

   int  GetNuserFloat() const {return NuserFloat;}
   void SetNuserFloat(int n){NuserFloat = n;}

   int  GetNcontrib() const {return Ncontrib;}
   void SetNcontrib(int n){Ncontrib = n;}

   int  GetNdata() const {return Ndata;}
   void SetNdata(int n){Ndata = n;}

   int  GetNuserInt() const {return NuserInt;}
   void SetNuserInt(int n){NuserInt = n;}

   int  GetImachine() const {return Imachine;}
   void SetImachine(int n){Imachine = n;}

   int  GetOutputPrecision() const {return fPrecision;}
   void SetOutputPrecision(int precision) {fPrecision = precision;}


protected:
   void StripWhitespace(string &str) const;
   void PutBackMagicNo(istream* table);
   bool ReadMagicNo(istream *table);
   void PrintWelcomeMessage();

   string ffilename;
   ifstream *ifilestream;
   ofstream *ofilestream;
   int fPrecision;
   // header
   int Itabversion;
   string ScenName;
   int Ncontrib;
   int Nmult;
   int Ndata;
   int NuserString;
   int NuserInt;
   int NuserFloat;
   int Imachine;

   static bool fWelcomeOnce;

};
#endif
