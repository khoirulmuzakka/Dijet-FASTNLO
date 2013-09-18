#ifndef __fastNLOCoeffAddFlex__
#define __fastNLOCoeffAddFlex__

#include "fastNLOCoeffAddBase.h"
#include "fastNLOConstants.h"

using namespace std;

class fastNLOCoeffAddFlex : public fastNLOCoeffAddBase {

   friend class fastNLOTable;
   friend class fastNLOReader;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddFlex(int NObsBin, int iLOord);
   fastNLOCoeffAddFlex(const fastNLOCoeffBase& base , int iLOord);
   virtual ~fastNLOCoeffAddFlex(){;}
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false) ;
   int Read(istream *table);
   void ReadRest(istream *table);
   virtual int Write(ostream *table, int option = 0);
   virtual int Copy(fastNLOCoeffAddFlex* other);
   virtual void Print() const;
   void Add(fastNLOCoeffAddFlex* other);

   template<typename T>  int ReadFlexibleVector(vector<T>* v, istream* table, bool nProcLast=false);
   int ReadFlexibleVector( vector<double >* v, istream *table , bool nProcLast = false );

   unsigned int GetNScaleNode1(int iObsBin) const { return ScaleNode1[iObsBin].size(); };
   unsigned int GetNScaleNode2(int iObsBin) const { return ScaleNode2[iObsBin].size(); };
   double GetScaleNode1(int iObsBin, int iNode) const { return ScaleNode1[iObsBin][iNode]; };
   double GetScaleNode2(int iObsBin, int iNode) const { return ScaleNode2[iObsBin][iNode]; };

protected:

   fastNLOCoeffAddFlex();
   int ReadCoeffAddFlex(istream *table);
   
   int fILOord;   // obtained from Scenario
   
   // SigmaTilde [NObsBins] ['n' x-nodes] [n s1-Nodes] [n s2-Nodes] [nsubproc]
   v5d SigmaTildeMuIndep; 
   v5d SigmaTildeMuFDep; 
   v5d SigmaTildeMuRDep; 
   v5d SigmaTildeMuRRDep; 
   v5d SigmaTildeMuFFDep; 
   v5d SigmaTildeMuRFDep; 
   // SigmaRef [NObsBins] [nsubproc]
   v2d SigmaRefMixed; 
   v2d SigmaRef_s1; 
   v2d SigmaRef_s2; 
   //int NscalenodeScale1;
   //int NscalenodeScale2;
   // ScaleNodeXY [ObsBin] [NscalenodeScaleX]  
   v2d ScaleNode1;
   v2d ScaleNode2;
   
public:
   v3d AlphasTwoPi;
   v5d PdfLcMuVar;

};


template<typename T>
int fastNLOCoeffAddFlex::ReadFlexibleVector(vector<T>* v, istream* table, bool nProcLast){
   int nn = 0;
   int size = 0;
   *table >> size; nn++;
   v->resize(size);
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += ReadFlexibleVector(&(v->at(i0)),table,nProcLast);
   }
   return nn;
};

#endif
