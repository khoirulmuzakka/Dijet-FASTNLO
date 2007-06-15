#ifndef __fnloBlockA2__
#define __fnloBlockA2__

#include <math.h>
#include <string>
#include <fstream>
#include <vector>

#include "fnloconstants.h"

using namespace std;

class fnloBlockA2 {
 public:
   int Read(istream *table);
   int Write(ostream *table);
   bool IsCompatible(fnloBlockA2* other);
 private:
   void StripWhitespace(string &str);
 protected:
   int Ipublunits;
   int NScDescript;
   vector <string> ScDescript;
   double Ecms;
   int ILOord;
   int NObsBin;
   int NDim;
   vector <string> DimLabel;
   vector <int> IDiffBin;
   vector < vector <double> > LoBin;
   vector < vector <double> > UpBin;
   int INormFlag;
   string DenomTable;
   vector <int> IDivPointer;
};
#endif
