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
   
 public:
   int Ipublunits;
   int NScDescript;
   vector <string> ScDescript;
   double Ecms;
   int ILOord;
   int NObsBin;
   int NDim;
   //KR: Added possibility to store and read start of new rapidity bin in nobs
   vector <int> RapIndex;
   vector <string> DimLabel;
   vector <int> IDiffBin;
   vector < vector <double> > LoBin;
   vector < vector <double> > UpBin;
   int INormFlag;
   string DenomTable;
   vector <int> IDivLoPointer;
   vector <int> IDivUpPointer;
};
#endif
