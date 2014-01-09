#include <cstdlib>
#include <unistd.h>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include "fastnlotk/fastNLOConstants.h"
#include "fastnlotk/fastNLOBase.h"

using namespace std;
using namespace fastNLO;


//______________________________________________________________________________
bool fastNLOBase::fWelcomeOnce = false;


//______________________________________________________________________________
fastNLOBase::fastNLOBase() : PrimalScream("fastNLOBase") ,  ifilestream(0), ofilestream(0), fPrecision(8) {
   if (!fWelcomeOnce) PrintWelcomeMessage();
}


//______________________________________________________________________________
fastNLOBase::fastNLOBase(string name) : PrimalScream("fastNLOBase") , ffilename(name), ifilestream(0), ofilestream(0), fPrecision(8)  {
   if (!fWelcomeOnce) PrintWelcomeMessage();
}


//______________________________________________________________________________
fastNLOBase::~fastNLOBase(){
   if(ifilestream) delete ifilestream;
   if(ofilestream) delete ofilestream;
}


//______________________________________________________________________________
int fastNLOBase::ReadTable(){
   // does file exist?
   if (access(ffilename.c_str(), R_OK) == 0) {
      error["ReadTable"]<<"File does not exist! Was looking for: "<<ffilename<<endl;
      return 1;
   }
   // open file
   ifilestream = new ifstream(ffilename.c_str(),ios::in);
   // read header
   ReadHeader(ifilestream);
   return 0;
}


//______________________________________________________________________________
bool fastNLOBase::ReadMagicNo(istream *table) {
   int key = 0;
   *table >> key;
   if(key != tablemagicno){
      error["ReadMagicNo"]<<"Found "<<key<<" instead of "<<tablemagicno<<"."<<endl;
      return false;
   };
   return true;
}


//______________________________________________________________________________
void fastNLOBase::PutBackMagicNo(istream* table){
   // Put magic number back
   for(int i=0;i<(int)(log10((double)tablemagicno)+1);i++){
      table->unget();
   }
}


//______________________________________________________________________________
void fastNLOBase::StripWhitespace(string &str) const{
   for(string::iterator achar = str.end(); achar>str.begin();achar--) {
      if (*achar==0x20 || *achar==0x00){
         str.erase(achar);
      }else{
         break;
      }
   }
}


//______________________________________________________________________________
int fastNLOBase::ReadHeader(istream *table){
   table->peek();
   if (table->eof()){
      error["ReadHeader"]<<"Cannot read from file."<<endl;
      return(2);
   }

   ReadMagicNo(table);
   *table >> Itabversion;
   *table >> ScenName;
   *table >> Ncontrib;
   *table >> Nmult;
   *table >> Ndata;
   *table >> NuserString;
   *table >> NuserInt;
   *table >> NuserFloat;
   *table >> Imachine;

   ReadMagicNo(table);
   PutBackMagicNo(table);
   return 0;
}


//______________________________________________________________________________
void fastNLOBase::RewindRead(){
   ifilestream->close();
   ifilestream->open(ffilename.c_str(),ios::in);
}


//______________________________________________________________________________
void fastNLOBase::SkipBlockA1A2(){
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


//______________________________________________________________________________
void fastNLOBase::WriteTable (){
   //
   // WriteTable(). writes the full FastNLO table to
   // the previously defined ffilename on disk.
   //
   // this function should be overwritten by
   // fastNLOTable::WriteTable();
   //
   OpenFileRewrite();
   WriteHeader(ofilestream);
   CloseFileWrite();
}


//______________________________________________________________________________
int fastNLOBase::WriteHeader(ostream *table){
   *table << tablemagicno << endl;
   *table << Itabversion << endl;
   *table << ScenName << endl;
   *table << Ncontrib << endl;
   *table << Nmult << endl;
   *table << Ndata << endl;
   *table << NuserString << endl;
   *table << NuserInt << endl;
   *table << NuserFloat << endl;
   *table << Imachine << endl;
   return 0;
}


//______________________________________________________________________________
ofstream *fastNLOBase::OpenFileWrite(){
   // do not overwrite existing table
   if (access(ffilename.c_str(), F_OK) == 0){
      printf("fastNLOBase::OpenFileWrite: File for writing the table exists: %s.\nPlease remove it.\n",ffilename.c_str());
      exit(2);
   }
   return OpenFileRewrite();
}


//______________________________________________________________________________
ofstream *fastNLOBase::OpenFileRewrite(){
   ofilestream = new ofstream(ffilename.c_str(),ios::out);
   if(!ofilestream->good()){
       printf("fastNLOBase::OpenFileWrite: Cannot open %s for writing. Stopping.\n",ffilename.c_str());
       exit(2);
    }
   ofilestream->precision(fPrecision);
   //ofilestream->precision(18);
   return ofilestream;
}


//______________________________________________________________________________
void fastNLOBase::CloseFileWrite(){
   *ofilestream << tablemagicno << endl;
   *ofilestream << tablemagicno << endl;
   CloseStream();
}


//______________________________________________________________________________
void fastNLOBase::CloseStream(){
   ofilestream->close();
}



//______________________________________________________________________________
bool fastNLOBase::IsCompatibleHeader(const fastNLOBase& other) const {
   if(Itabversion!= other.GetItabversion()){
      warn["IsCompatibleHeader"]<<"Differing versions of table format: "<<Itabversion<<" and "<< other.GetItabversion()<<endl;
      return false;
   }
   if(Ndata + other.GetNdata() > 1){
      warn["IsCompatibleHeader"]<<"Two tables containing both experimental data are incompatible"<<endl;
      return false;
   }
   if(ScenName!= other.GetScenName()){
      warn["IsCompatibleHeader"]<<"Differing names of scenarios: "<<ScenName.c_str()<<" and "<<other.ScenName.c_str()<<endl;
      // continue...
   }
   return true;
}


//______________________________________________________________________________
void fastNLOBase::SetHeaderDefaults(){
   // TableMagicNo and ITabVersion are defined as constant in fastNLOConstants.h
   SetScenName("tns2000");
   SetContributionHeader();
}


//______________________________________________________________________________
void fastNLOBase::SetContributionHeader(){
   SetNcontrib(1);
   SetNmult(0);
   SetNdata(0);
   SetNuserString(0);
   SetNuserInt(0);
   SetNuserFloat(0);
   SetImachine(0);
}


//______________________________________________________________________________
void fastNLOBase::ResetHeader(){
   debug["ResetHeader"]<<endl;
   SetNcontrib(0);
   SetNmult(0);
   SetNdata(0);
   SetNuserString(0);
   SetNuserInt(0);
   SetNuserFloat(0);
   SetImachine(0);
}


//______________________________________________________________________________
void fastNLOBase::Print() const {
   PrintHeader();
}


//______________________________________________________________________________
void fastNLOBase::PrintHeader() const {
   printf("\n **************** FastNLO Table Header ******************\n\n");
   printf("   tablemagicno                  %d\n",tablemagicno);
   printf("   Itabversion                   %d\n",Itabversion);
   printf("   ScenName                      %s\n",ScenName.data());
   printf("   Ncontrib                      %d\n",Ncontrib);
   printf("   Nmult                         %d\n",Nmult);
   printf("   Ndata                         %d\n",Ndata);
   printf("   NuserString                   %d\n",NuserString);
   printf("   NuserInt                      %d\n",NuserInt);
   printf("   NuserFloat                    %d\n",NuserFloat);
   printf("   Imachine                      %d\n",Imachine);
   printf("\n ********************************************************\n\n");
}


//______________________________________________________________________________
void fastNLOBase::PrintWelcomeMessage() {
   //---  Initialization for nice printing
   const string CSEPS = " ##################################################################################\n";
   const string LSEPS = " #---------------------------------------------------------------------------------\n";

   char fnlo[100];
   //sprintf(fnlo,"27[0;31mfast27[0;34mNLO\033[0m",27,0,31,27,0,34);
   sprintf(fnlo,"%c[%d;%dmfast%c[%d;%dmNLO\033[0m",27,0,31,27,0,34);
   char package_version[100]    = FNLO_VERSION;
   char svnrev[100]             = FNLO_SVNREV;
   char authors[500]            = FNLO_AUTHORS;
   char webpage[500]    = FNLO_WEBPAGE;
   char authorsv14[200] = FNLO_AUTHORSv14;
   char quotev14[200]   = FNLO_QUOTEv14;
   char authorsv2[200]  = FNLO_AUTHORS;
   char quotev2[200]    = FNLO_QUOTEv2;

   shout>>"\n";
   shout>>""<<CSEPS;
   shout<<"\n";
   shout<<" "<<fnlo<<endl;
   shout<<" Version "<<package_version<<"_"<<svnrev<<endl;
   shout<<"\n";
   shout<<" C++ program and toolkit to read and create fastNLO v2 tables and"<<endl;
   shout<<" derive QCD cross sections using PDFs, e.g. from LHAPDF"<<endl;
   shout<<"\n";
   shout>>""<<LSEPS;
   shout<<"\n";
   shout<<" Copyright © 2011,2012,2013 "<<fnlo<<" Collaboration"<<endl;
   shout<<" "<<authors<<endl;
   shout<<"\n";
   shout>>" # This program is free software: you can redistribute it and/or modify"<<endl;
   shout>>" # it under the terms of the GNU General Public License as published by"<<endl;
   shout>>" # the Free Software Foundation, either version 3 of the License, or"<<endl;
   shout>>" # (at your option) any later version."<<endl;
   shout>>" #\n";
   shout>>" # This program is distributed in the hope that it will be useful,"<<endl;
   shout>>" # but WITHOUT ANY WARRANTY; without even the implied warranty of"<<endl;
   shout>>" # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the"<<endl;
   shout>>" # GNU General Public License for more details."<<endl;
   shout>>" #\n";
   shout>>" # You should have received a copy of the GNU General Public License"<<endl;
   shout>>" # along with this program. If not, see <http://www.gnu.org/licenses/>."<<endl;
   shout>>" #\n";
   shout>>""<<LSEPS;
   shout>>" #\n";
   shout<<" The projects web page can be found at:"<<endl;
   shout<<"   "<<webpage<<endl;
   shout<<"\n";
   shout<<" If you use this code, please cite:"<<endl;
   shout<<"   "<<authorsv14<<", "<<quotev14<<endl;
   shout<<"   "<<authorsv2<<", "<<quotev2<<endl;
   shout<<"\n";
   shout>>""<<CSEPS;
   fWelcomeOnce = true;
}

