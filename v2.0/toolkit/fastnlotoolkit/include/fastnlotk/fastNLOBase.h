#ifndef __fastNLOBase__
#define __fastNLOBase__

#include <fstream>
#include <ostream>
#include <istream>
#include <string>
#include "speaker.h"


class fastNLOBase {

public:
   fastNLOBase();
   fastNLOBase(std::string name);
   fastNLOBase(const fastNLOBase&);
   virtual ~fastNLOBase();

   // i/o
   virtual void ReadTable();                                                    //!< read table
   virtual void WriteTable();                                                   //!< write full table to disk
   virtual void Print() const;
   // header
   void PrintHeader() const;                                            //!< Print header variables (BlockA1) to screen
   void SetHeaderDefaults();                                            //!< Set some default values
   void ResetHeader();                                                  //!< Reset variables to default values
   void SetContributionHeader();                                        //
   bool IsCompatibleHeader(const fastNLOBase& other) const;             //!< Compare header with header of another table

   // getter/setters
   std::string GetFilename() const {return ffilename;}
   void   SetFilename(std::string name){ffilename=name;}

   int  GetItabversion() const {return Itabversion;}
   void SetItabversion(int version){Itabversion = version;}

   std::string GetScenName() const {return ScenName;}
   void   SetScenName(std::string name){ScenName = name;}

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
   void PrintWelcomeMessage();                                          //!< Say hello to fastNLO user
   std::ofstream* OpenFileWrite();                                           //!< open std::ofstream for writing tables to ffilename
   std::ifstream* OpenFileRead();                                            //!< open std::ifstream for reading table
   //std::ofstream *OpenFileRewrite();
   void WriteHeader(std::ostream& table);                                    //!< write (or cout) hader using std::ostream
   void ReadHeader(std::istream& table);                                     //!< read header of table (BlockA1)
   void CloseFileWrite(std::ofstream& table);
   void CloseFileRead(std::ifstream& table);
   //void CloseStream();

   std::string ffilename;
   //std::ifstream *ifilestream;
   //std::ofstream *ofilestream;
   int fPrecision;
   // header
   int Itabversion;
   std::string ScenName;
   int Ncontrib;
   int Nmult;
   int Ndata;
   int NuserString;
   int NuserInt;
   int NuserFloat;
   int Imachine;

   PrimalScream logger;
   static bool fWelcomeOnce;

};
#endif
