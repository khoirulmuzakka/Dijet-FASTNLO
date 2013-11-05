#ifndef __fastNLOCoeffBase__
#define __fastNLOCoeffBase__

#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <iostream>

#include "fastNLOConstants.h"
#include "speaker.h"

using namespace std;
using namespace fastNLO;

class fastNLOCoeffBase : public PrimalScream {
   friend class fastNLOTable;


public:
   fastNLOCoeffBase(int NObsBin);
   virtual ~fastNLOCoeffBase(){;};
   //fastNLOCoeffBase(const fastNLOCoeffBase& coeff);
   int Read(istream *table);
   int ReadBase(istream *table);
   virtual int Write(ostream *table, int option = 0);
   virtual int Copy(fastNLOCoeffBase* other);
   //void Add(fastNLOCoeffBase* other);
   int EndReadCoeff(istream *table);
   virtual void Print() const;

   void SetNlojetDefaults();
   void SetNlojetDescr();

   void StripWhitespace(string& s) const;

   // there are nicer  options in c++11
   void ResizeTable( v1d* v, int dim0 );
   void ResizeTable( v2d*  v, int dim0 , int dim1 );
   void ResizeTable( v3d* v, int dim0 , int dim1, int dim2 );
   void ResizeTable( v4d* v, int dim0 , int dim1, int dim2, int dim3 );
   void ResizeTable( v5d* v, int dim0 , int dim1, int dim2, int dim3, int dim4 );
   void ResizeTable( v6d* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 );
   void ResizeTable( v7d* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 );

   template<typename T> void ResizeFlexibleVector(vector<T>* v, vector<T>* nom);
   void ResizeFlexibleVector(vector<double >* v, vector<double >*nom) { v->resize(nom->size());}

   template<typename T> int ReadTable( vector<T>* v, istream *table );
   int ReadTable( vector<double>* v, istream *table );
   template<typename T> int WriteFlexibleTable( vector<T>* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   int WriteFlexibleTable( vector<double >* v, ostream *table , bool DivByNevt=false, int Nevt=1 , bool nProcLast = false );
   template<typename T> void AddTableToAnotherTable( vector<T>* vSum, vector<T>* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( vector<double >* vSum, vector<double >* vAdd, double w1 = 1, double w2 = 1 );

   template<typename T> int WriteTable( vector<T>* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   int WriteTable( vector<double >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   // replacement of functions with template not entirely tested
   //    int WriteTable( vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, ostream *table , bool DivByNevt=false , int Nevt=1 );
   //    int WriteTable( vector<vector<vector<vector<vector<vector<double > > > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   //    int WriteTable( vector<vector<vector<vector<vector<double > > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   //    int WriteTable( vector<vector<vector<vector<double > > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   //    int WriteTable( vector<vector<vector<double > > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   //    int WriteTable( vector<vector<double > >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );
   //    int WriteTable( vector<double >* v, ostream *table , bool DivByNevt=false, int Nevt=1 );


   int GetIDataFlag() const {return IDataFlag;}
   int GetIAddMultFlag() const {return IAddMultFlag;}
   int GetIContrFlag1() const {return IContrFlag1;}
   int GetIContrFlag2() const {return IContrFlag2;}
   int GetNScaleDep() const {return NScaleDep;}
   vector<string > GetContributionDescription() const { return CtrbDescript; }
   vector<string > GetCodeDescription() const { return CodeDescript; }

   void SetIXsectUnits(int n){IXsectUnits = n;}
   void SetIDataFlag(int n){IDataFlag = n;}
   void SetIAddMultFlag(int n){IAddMultFlag = n;}
   void SetIContrFlag1(int n){IContrFlag1 = n;} 
   void SetIContrFlag2(int n){IContrFlag2 = n;} 
   void SetNScaleDep(int n){NScaleDep = n;}
   
   bool IsLO() const {return IContrFlag1==1 && IContrFlag2==1;}
   bool IsNLO() const {return IContrFlag1==1 && IContrFlag2==2;}

   static const int DividebyNevt = 1; // shitty definition of a global constant

protected:
   fastNLOCoeffBase();

   int fNObsBins; // obtained from Scenario

   int IXsectUnits;
   int IDataFlag;
   int IAddMultFlag;
   int IContrFlag1;
   int IContrFlag2;
   int NScaleDep;
   vector < string > CtrbDescript;
   vector < string > CodeDescript;

};



template<typename T>
void fastNLOCoeffBase::ResizeFlexibleVector(vector<T>* v, vector<T>* nom) {
   v->resize(nom->size());
   for (unsigned int i = 0 ; i<v->size() ; i++) {
      ResizeFlexibleVector(&((*v)[i]),&((*nom)[i]));
   }
};

template<typename T> 
int fastNLOCoeffBase::WriteTable( vector<T>* v, ostream *table , bool DivByNevt, int Nevt ){
   int nn = 0;
   for(unsigned int i=0;i<v->size();i++){
      nn += WriteTable( &(*v)[i] , table , DivByNevt , Nevt );
   }   
   return nn;
};


template<typename T> 
int fastNLOCoeffBase::ReadTable( vector<T>* v, istream *table ){
   int nn = 0;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn+= ReadTable(&(*v)[i0],table);
   }
   return nn;
};

template<typename T> 
int fastNLOCoeffBase::WriteFlexibleTable( vector<T>* v, ostream *table , bool DivByNevt, int Nevt , bool nProcLast ){
   int nn = 1;
   *table << v->size() << endl;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn += WriteFlexibleTable( &(v->at(i0)) , table , DivByNevt, Nevt , nProcLast );
   }
   return nn;
};
 
template<typename T> 
void fastNLOCoeffBase::AddTableToAnotherTable( vector<T>* vSum, vector<T>* vAdd, double w1, double w2 ){
   if ( vSum->size() != vAdd->size() ) {error["AddTableToAnotherTable"]<<"Cannot add tables with different size. s1="<<vSum->size()<<", s2="<<vAdd->size()<<endl; return;}
   for ( unsigned int i = 0 ; i<vSum->size() ; i++ ){
      AddTableToAnotherTable( &(vSum->at(i)), &(vAdd->at(i)), w1 , w2  );
   }
}


#endif