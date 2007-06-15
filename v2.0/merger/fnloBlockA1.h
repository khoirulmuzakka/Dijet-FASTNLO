#ifndef __fnloBlockA1__
#define __fnloBlockA1__

#include <string>
#include <fstream>

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
 protected:
   int Itabversion;
   string ScenName;
   int Ncontrib;
   int Nmult;
   int Ndata;
};
#endif
