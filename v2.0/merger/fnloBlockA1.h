#ifndef __fnloBlockA1__
#define __fnloBlockA1__

#include <string>
#include <fstream>
#include <math.h>

#include "fnloconstants.h"

using namespace std;

class fnloBlockA1 {
 public:
   int Read(istream *table);
   int Write(ostream *table);
   bool IsCompatible(fnloBlockA1* other);
   int GetItabversion(){return Itabversion;}
   string GetScenName(){return ScenName;}
   int GetNcontrib(){return Ncontrib;}
   int GetNmult(){return Nmult;}
   int GetNdata(){return Ndata;}
   void SetNcontrib(int n){Ncontrib = n;}
   void SetNmult(int n){Nmult = n;}
   void SetNdata(int n){Ndata = n;}
 protected:
   int Itabversion;
   string ScenName;
   int Ncontrib;
   int Nmult;
   int Ndata;
};
#endif
