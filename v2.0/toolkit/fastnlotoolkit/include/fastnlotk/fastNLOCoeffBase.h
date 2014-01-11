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
   virtual void Write(ostream *table);
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

   template<typename T> int ReadTable( vector<T>* v, istream *table , unsigned long long Nevt = 1);
   int ReadTable( vector<double>* v, istream *table , unsigned long long Nevt = 1);

   template<typename T> void AddTableToAnotherTable( vector<T>* vSum, vector<T>* vAdd, double w1 = 1, double w2 = 1 );
   void AddTableToAnotherTable( vector<double >* vSum, vector<double >* vAdd, double w1 = 1, double w2 = 1 );

   template<typename T> void AddTableToAnotherTable( vector<T>& vSum, const vector<T>& vAdd, double w1 = 1, double w2 = 1 ) const;
   void AddTableToAnotherTable( vector<double >& vSum, const vector<double >& vAdd, double w1 = 1, double w2 = 1 ) const;

   template<typename T> int WriteFlexibleTable( const vector<T>& v, ostream *table , bool nProcLast=true, double Nevt=1 ) const ;
   int WriteFlexibleTable( const vector<double >& v, ostream *table , bool nProcLast=true, double Nevt=1 ) const ;

   template<typename T> int WriteTable( const vector<T>& v, ostream *table , unsigned long long Nevt=1 ) const;
   int WriteTable( const vector<double >& v, ostream *table , unsigned long long Nevt=1 ) const;

   int GetIDataFlag() const {return IDataFlag;}
   void SetIDataFlag(int n){IDataFlag = n;}

   int GetIAddMultFlag() const {return IAddMultFlag;}
   void SetIAddMultFlag(int n){IAddMultFlag = n;}

   int GetIContrFlag1() const {return IContrFlag1;}
   void SetIContrFlag1(int n){IContrFlag1 = n;}

   int GetIContrFlag2() const {return IContrFlag2;}
   void SetIContrFlag2(int n){IContrFlag2 = n;}

   int GetNScaleDep() const {return NScaleDep;}
   void SetNScaleDep(int n){NScaleDep = n;}

   int GetIXsectUnits() const { return IXsectUnits;}
   void SetIXsectUnits(int n){IXsectUnits = n;}

   int GetNObsBin() const { return fNObsBins;}

   bool GetIsFlexibleScale() const { return (NScaleDep>=3) && (IAddMultFlag==0); }

   vector<string > GetContributionDescription() const { return CtrbDescript; }
   vector<string > GetCodeDescription() const { return CodeDescript; }


   bool IsLO() const {return IContrFlag1==1 && IContrFlag2==1;}
   bool IsNLO() const {return IContrFlag1==1 && IContrFlag2==2;}

   bool IsCompatible(const fastNLOCoeffBase& other) const;
   //bool operator==(const fastNLOCoeffBase& other) const { return IsCompatible(other); }


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
int fastNLOCoeffBase::ReadTable( vector<T>* v, istream *table , unsigned long long Nevt){
   int nn = 0;
   for(unsigned int i0=0;i0<v->size();i0++){
      nn+= ReadTable(&(*v)[i0],table, Nevt);
   }
   return nn;
};


template<typename T>
int fastNLOCoeffBase::WriteTable( const vector<T>& v, ostream *table , unsigned long long Nevt) const {
   int nn = 0;
   for(unsigned int i=0;i<v.size();i++){
      nn += WriteTable( v[i] , table , Nevt );
   }
   return nn;
};


template<typename T>
int fastNLOCoeffBase::WriteFlexibleTable( const vector<T>& v, ostream *table , bool nProcLast, double Nevt ) const {
   if ( Nevt== 0 ) return -1000;
   int nn = 1;
   *table << v.size() << endl;
   for(unsigned int i0=0;i0<v.size();i0++){
      nn += WriteFlexibleTable( v[i0] , table , nProcLast, Nevt );
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

template<typename T>
void fastNLOCoeffBase::AddTableToAnotherTable( vector<T>& vSum, const vector<T>& vAdd, double w1, double w2 ) const{
   if ( vSum.size() != vAdd.size() ) {
      error["AddTableToAnotherTable"]<<"Cannot add tables with different size. s1="<<vSum.size()<<", s2="<<vAdd.size()<<endl;
      return;
   }
   for ( unsigned int i = 0 ; i<vSum.size() ; i++ ){
      AddTableToAnotherTable( vSum[i], vAdd[i], w1 , w2  );
   }
}


#endif
