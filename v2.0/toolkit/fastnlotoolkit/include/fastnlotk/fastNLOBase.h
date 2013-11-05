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
   virtual ~fastNLOBase();

   // i/o
   int ReadTable();							// read table
   int ReadHeader(istream *table);					// read header of table (BlockA1)
   
   int OpenFileRead();							// open stream
   void WriteTable();							// write full table to disk
   int WriteHeader(ostream *table);					// write hader using ostream
   void RewindRead();							
   void SkipBlockA1A2();
   ofstream *OpenFileWrite();
   ofstream *OpenFileRewrite();
   void CloseFileWrite();
   void CloseStream();

   virtual void Print() const;
   
   // header
   void PrintHeader() const;						// Print header variables (BlockA1) to screen
   void SetHeaderDefaults();						// Set some default values 
   void ResetHeader();							// Reset variables to default values
   void SetContributionHeader();					// 	
   bool IsCompatibleHeader(fastNLOBase* other) const;			// Compare header with header of another table

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
   void PutBackMagicNo(istream* table);					// Reset magic number, such that it can be recognized by other reading routines
   bool ReadMagicNo(istream *table);					// read and crosscheck magic number
   void PrintWelcomeMessage();						// Say hello to fastNLO user

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