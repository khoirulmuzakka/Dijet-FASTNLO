#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffBase.h"

using namespace std;
using namespace fastNLO;

//________________________________________________________________________________________________________________ //
fastNLOCoeffBase::fastNLOCoeffBase() : PrimalScream("fastNLOCoeffBase"){
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffBase::fastNLOCoeffBase(int NObsBin) : PrimalScream("fastNLOCoeffBase") {
   fNObsBins = NObsBin;
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::Read(istream* table){
   // basic read function.
   // reads in only 'base'-variables.
   debug["Read"]<<endl;
   ReadBase(table);
   EndReadCoeff(table);
   return 0;
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::ReadBase(istream* table){
   debug["ReadBase"]<<endl;
   table->peek();
   if (table->eof()){
      //printf("fastNLOCoeffBase::Read: Cannot read from file.\n");
      error["ReadBase"]<<"Cannot read from file."<<endl;
      return(2);
   }

   int key = 0;
   *table >> key;
   if(key != tablemagicno){
      //printf("fastNLOCoeffBase::Read: At beginning of block found %d instead of %d.\n",key,tablemagicno);
      warn["ReadBase"]<<"At beginning of block found "<<key<<" instead of "<<tablemagicno<<endl;
      return 1;
   };

   *table >> IXsectUnits;
   *table >> IDataFlag;
   *table >> IAddMultFlag;
   *table >> IContrFlag1;
   *table >> IContrFlag2;
   *table >> NScaleDep;
   int NContrDescr;
   *table >> NContrDescr;
   //   printf("  *  fastNLOCoeffBase::Read().  IDataFlag: %d, IAddMultFlag: %d, IContrFlag1: %d, IContrFlag2: %d,, NScaleDep: %d\n",IDataFlag,IAddMultFlag,IContrFlag1,IContrFlag2,NScaleDep );
   CtrbDescript.resize(NContrDescr);
   char buffer[5257];
   table->getline(buffer,5256);
   for(int i=0;i<NContrDescr;i++){
      table->getline(buffer,256);
      CtrbDescript[i] = buffer;
      //      StripWhitespace(CtrbDescript[i]);
   }
   int NCodeDescr;
   *table >> NCodeDescr;
   CodeDescript.resize(NCodeDescr);
   table->getline(buffer,256);
   for(int i=0;i<NCodeDescr;i++){
      table->getline(buffer,256);
      CodeDescript[i] = buffer;
      //      StripWhitespace(CodeDescript[i]);
   }
   return 0;
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::EndReadCoeff(istream* table){
   //debug["EndReadCoeff"]<<endl;
   int key = 0;
   *table >> key;
   if(key != tablemagicno){
      //printf("fastNLOCoeffBase::Read: At end of block found %d instead of %d.\n",key,tablemagicno);
      //printf("                  You might have 'nan' in your table.\n");
      error["EndReadCoeff"]<<"At end of block found "<<key<<" instead of "<<tablemagicno<<endl;
      return 1;
   };
   // Put magic number back
   for(int i=0;i<(int)(log10((double)key)+1);i++){
      table->unget();
   }
   return 0;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::Write(ostream* table,double) {
   debug["Write"]<<"Writing fastNLOCoeffBase."<<endl;
   *table << tablemagicno << endl;
   *table << IXsectUnits << endl;
   *table << IDataFlag << endl;
   *table << IAddMultFlag << endl;
   *table << IContrFlag1 << endl;
   *table << IContrFlag2 << endl;
   *table << NScaleDep << endl;
   *table << CtrbDescript.size() << endl;
   //printf("  *  fastNLOCoeffBase::Write().  IDataFlag: %d, IAddMultFlag: %d, IContrFlag1: %d, IContrFlag2: %d, NScaleDep: %d\n",
   //IDataFlag,IAddMultFlag,IContrFlag1,IContrFlag2,NScaleDep);
   for(unsigned int i=0;i<CtrbDescript.size();i++){
      *table << CtrbDescript[i] << endl;
   }
   *table << CodeDescript.size() << endl;
   for(unsigned int i=0;i<CodeDescript.size();i++){
      *table << CodeDescript[i] << endl;
   }
}


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffBase::IsCompatible(const fastNLOCoeffBase& other) const {
   if( fNObsBins != other.GetNObsBin() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   if( IXsectUnits != other.GetIXsectUnits() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   if( IDataFlag != other.GetIDataFlag() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   if( IAddMultFlag != other.GetIAddMultFlag() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   if( IContrFlag1 != other.GetIContrFlag1() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   if( IContrFlag2 != other.GetIContrFlag2() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   if( NScaleDep != other.GetNScaleDep() ){
      //warn["IsCompatible"]<<"."<<endl;
      return false;
   }
   // chcekc descripts here ?!
   //bool potentialcompatible = true;
   //vector < string > CtrbDescript;
   //vector < string > CodeDescript;
   return true;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::StripWhitespace(string &str) const {
   for(string::iterator achar = str.end(); achar>str.begin();achar--) {
      if (*achar==0x20 || *achar==0x00){
         str.erase(achar);
      }else{
         break;
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::ResizeTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2, dim3, dim4, dim5, dim6 );
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::ResizeTable( vector<vector<vector<vector<vector<vector<double > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2, dim3, dim4, dim5 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}




//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2, dim3, dim4 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::ResizeTable( vector<vector<vector<vector<double > > > >* v, int dim0 , int dim1, int dim2, int dim3 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2, dim3 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::ResizeTable( vector<vector<vector<double > > >* v, int dim0 , int dim1, int dim2 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1, dim2 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::ResizeTable( vector<vector<double > >*  v, int dim0 , int dim1 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
    for ( int i= 0 ; i<dim0 ; i++){
      ResizeTable( &(v->at(i)) , dim1 );
    }
  } else {
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::ResizeTable( vector<double >* v, int dim0 ){
  if ( dim0 > 0 ){
    v->resize(dim0);
  }
  else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}



//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::ReadTable(vector<double >* v, istream *table ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    *table >> v->at(i0);
    nn++;
  }
  return nn;
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::WriteTable( const vector<double >& v, ostream *table , double Nevt ) const {
   if( Nevt == 0) return -1000;
   int nn = 0;
   for(unsigned int i0=0;i0<v.size();i0++){
      *table << v[i0] / Nevt << endl;
      nn++;
   }
   return nn;
}



//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::WriteFlexibleTable(const vector<double >& v, ostream *table, bool nProcLast, double Nevt) const {
   int nn = 1;
   if ( Nevt == 0 ) return -1000;
   if ( !nProcLast )*table << v.size() << endl;
   for(unsigned int i0=0;i0<v.size();i0++){
      *table << v[i0] / Nevt << endl;
      nn++;
   }
   return nn;
}

//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::AddTableToAnotherTable( vector<double >* vSum, vector<double >* vAdd, double w1, double w2){
   if ( vSum->size() != vAdd->size() ) {cout<<"Error in fastNLOCoeffBase::AddTableToAnotherTable. Cannot add tables with different size. [v1] s1="<<vSum->size()<<", s2="<<vAdd->size()<<endl; return;}
  for ( unsigned int i = 0 ; i<vSum->size() ; i++ ){
    (*vSum)[i] =  w1*(*vSum)[i] + w2*(*vAdd)[i];
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffBase::AddTableToAnotherTable( vector<double >& vSum, const vector<double >& vAdd, double w1, double w2) const {
   if ( vSum.size() != vAdd.size() ) {
      error["AddTableToAnotherTable"]<<"Cannot add tables with different size. [v1] s1="<<vSum.size()<<", s2="<<vAdd.size()<<endl;
   }
   else {
      for ( unsigned int i = 0 ; i<vSum.size() ; i++ ){
	 vSum[i] =  w1*vSum[i] + w2*vAdd[i];
      }
   }
}



//________________________________________________________________________________________________________________ //


void fastNLOCoeffBase::SetNlojetDefaults(){
  SetIDataFlag(0);
  SetIAddMultFlag(0);
  SetIContrFlag1(1);
  SetIContrFlag2(100); // specify if LO or NLO
  SetNScaleDep(0);
  SetIXsectUnits(12);
  SetNlojetDescr();
};

void fastNLOCoeffBase::SetNlojetDescr(){
   CodeDescript.push_back("NLOJet++_4.1.3");
   CodeDescript.push_back("Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),");
   CodeDescript.push_back("Z. Nagy, Phys. Rev. D68, 094002 (2003).");
 }


void fastNLOCoeffBase::Print() const {
  printf("\n **************** FastNLO Table: CoeffBase ****************\n");
  printf(" B   Scenario::GetNObsBin()        %d\n",fNObsBins);
  printf(" B   IXsectUnits                   %d\n",IXsectUnits);
  printf(" B   IDataFlag                     %d\n",IDataFlag);
  printf(" B   IAddMultFlag                  %d\n",IAddMultFlag);
  printf(" B   IContrFlag1                   %d\n",IContrFlag1);
  printf(" B   IContrFlag2                   %d\n",IContrFlag2);
  //  printf(" B   IContrFlag3 (always 0)        %d\n",IContrFlag3);
  printf(" B   NScaleDep                     %d\n",NScaleDep);
  for(unsigned int i=0;i<CtrbDescript.size();i++){
    printf(" B   CtrbDescript[%d]               %s\n",i,CtrbDescript[i].data());
  }
  //printf(" B   NCodeDescr                    %d\n",NCodeDescr);
  for(unsigned int i=0;i<CodeDescript.size();i++){
    printf(" B   CodeDescript[%d]               %s\n",i,CodeDescript[i].data());
  }
  printf(" *******************************************************\n");

}


//________________________________________________________________________________________________________________ //
