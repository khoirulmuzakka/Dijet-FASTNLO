#ifndef __fnloBlockA1__
#define __fnloBlockA1__

#include <string>
#include <fstream>
#include <math.h>

#include "fnloConstants.h"

using namespace std;

class fnloBlockA1 {
 public:
   fnloBlockA1(){Itabversion=tabversion;}
   int Read(istream *table);
   int Write(ostream *table);
   bool IsCompatible(fnloBlockA1* other);
   int GetItabversion(){return Itabversion;}
   string GetScenName(){return ScenName;}
   int GetNcontrib(){return Ncontrib;}
   int GetNmult(){return Nmult;}
   int GetNdata(){return Ndata;}
   int GetNuserString(){return NuserString;}
   int GetNuserInt(){return NuserInt;}
   int GetNuserFloat(){return NuserFloat;}
   int GetImachine(){return Imachine;}
   void SetScenName(string name){ScenName = name;}
   void SetNcontrib(int n){Ncontrib = n;}
   void SetNmult(int n){Nmult = n;}
   void SetNdata(int n){Ndata = n;}
   void SetNuserString(int n){NuserString = n;}
   void SetNuserInt(int n){NuserInt = n;}
   void SetNuserFloat(int n){NuserFloat = n;}
   void SetImachine(int n){Imachine = n;}
   void Print();

 protected:
   int Itabversion;
   string ScenName;
   int Ncontrib;
   int Nmult;
   int Ndata;
   int NuserString;
   int NuserInt;
   int NuserFloat;
   int Imachine;
};
#endif
