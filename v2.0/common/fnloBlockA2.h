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
   int GetNObsBin(){return NObsBin;}
   int GetILOord(){return ILOord;}
   void SetIpublunits(int unit){Ipublunits = unit;}
 private:
   void StripWhitespace(string &str) const;
   bool cmp(const double x1, const double x2) const;
   bool cmp(const vector < double > x1, const vector < double > x2) const;
   bool cmp(const vector < vector < double > > x1,const vector < vector < double > > x2) const;
   
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
