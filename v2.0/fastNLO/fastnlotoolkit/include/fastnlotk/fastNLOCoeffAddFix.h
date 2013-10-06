#ifndef __fastNLOCoeffAddFix__
#define __fastNLOCoeffAddFix__

#include "fastNLOCoeffAddBase.h"
#include "fastNLOConstants.h"

using namespace std;

class fastNLOCoeffAddFix : public fastNLOCoeffAddBase {

   friend class fastNLOReader;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddFix(int NObsBin);
   fastNLOCoeffAddFix(const fastNLOCoeffBase&);
   virtual ~fastNLOCoeffAddFix(){;}
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   int Read(istream *table);
   void ReadRest(istream *table);
   virtual int Write(ostream *table, int option = 0);
   virtual int Copy(fastNLOCoeffAddFix* other);
   void Add(fastNLOCoeffAddFix* other);
   virtual void Print() const;
   
   int GetTotalScalevars() const ;
   int GetTotalScalenodes() const ;
   int GetNScalevar() const { return Nscalevar[0];}
   v1d GetAvailableScaleFactors() const { return ScaleFac[0]; }
   int GetNScaleNode() const { return Nscalenode[0]; }
   double GetScaleFactor(int iVar) const { 
      if ( iVar >= (int)ScaleFac[0].size() ) error["GetScaleFactor"]<<"Scalevariation no. "<<iVar<<" not available. There are only "<<GetNScalevar()<<" available in this table."<<endl;
      return ScaleFac[0][iVar];
   }

   double GetSigmaTilde(int iObs, int iSvar, int ix, int is, int iN ) const { return SigmaTilde[iObs][iSvar][ix][is][iN];}
   double GetScaleNode(int iObs, int iSvar, int iNode ) const { return ScaleNode[iObs][0][iSvar][iNode];}

protected:
   fastNLOCoeffAddFix();
   int ReadCoeffAddFix(istream *table);

   vector < int > Nscalevar;
   vector < int > Nscalenode;
   v2d  ScaleFac;
   v4d ScaleNode;
   v5d SigmaTilde;

public:
   v2d AlphasTwoPi_v20;
   v4d PdfLc;
};

#endif
