#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffBase.h"

using namespace std;
using namespace fastNLO;

fastNLOCoeffBase::fastNLOCoeffBase() : PrimalScream("fastNLOCoeffBase"){
}

fastNLOCoeffBase::fastNLOCoeffBase(int NObsBin) : PrimalScream("fastNLOCoeffBase") {
   fNObsBins = NObsBin;
}

int fastNLOCoeffBase::Read(istream *table){
   // basic read function.
   // reads in only 'base'-variables.
   debug["Read"]<<endl;
   ReadBase(table);
   EndReadCoeff(table);
   return 0;
}


int fastNLOCoeffBase::ReadBase(istream *table){
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


int fastNLOCoeffBase::EndReadCoeff(istream *table){
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


int fastNLOCoeffBase::Write(ostream *table, int option){
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

   return 0;
}

int fastNLOCoeffBase::Copy(fastNLOCoeffBase* other){
   debug["Copy"]<<endl;
   streambuf* streambuf = new stringbuf(ios_base::in | ios_base::out);
   iostream* buffer = new iostream(streambuf);
   other->Write(buffer);
   *buffer << tablemagicno << endl;
   this->Read(buffer);
   delete buffer;
   delete streambuf;

   return(0);
}



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

/*
//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::WriteTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, ostream *table , bool DivByNevt , int Nevt  ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
	    for(unsigned int i5=0;i5<v->at(i0)[i1][i2][i3][i4].size();i5++){
	      for(unsigned int i6=0;i6<v->at(i0)[i1][i2][i3][i4][i5].size();i6++){
		if( DivByNevt && Nevt>0){
 		  *table << v->at(i0)[i1][i2][i3][i4][i5][i6] / Nevt << endl;
 		}else{
 		  *table << v->at(i0)[i1][i2][i3][i4][i5][i6] << endl;
		}
		nn++;
	      }
	    }
	  }
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::WriteTable(vector<vector<vector<vector<vector<vector<double > > > > > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
	    for(unsigned int i5=0;i5<v->at(i0)[i1][i2][i3][i4].size();i5++){
	      if( DivByNevt && Nevt>0){
 		*table << v->at(i0)[i1][i2][i3][i4][i5] / Nevt << endl;
	      }else{
		*table << v->at(i0)[i1][i2][i3][i4][i5] << endl;
	      }
	      nn++;
	    }
	  }
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::WriteTable(vector<vector<vector<vector<vector<double > > > > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  for(unsigned int i4=0;i4<v->at(i0)[i1][i2][i3].size();i4++){
	    if( DivByNevt && Nevt>0){
 	      *table << v->at(i0)[i1][i2][i3][i4] / Nevt << endl;
		}else{
	      *table << v->at(i0)[i1][i2][i3][i4] << endl;
	    }
	    nn++;
	  }
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::WriteTable(vector<vector<vector<vector<double > > > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
 	for(unsigned int i3=0;i3<v->at(i0)[i1][i2].size();i3++){
	  if( DivByNevt && Nevt>0){
 	    *table << v->at(i0)[i1][i2][i3] / Nevt << endl;
	  }else{
	    *table << v->at(i0)[i1][i2][i3] << endl;
	  }
	  nn++;
	}
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::WriteTable(vector<vector<vector<double > > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      for(unsigned int i2=0;i2<v->at(i0)[i1].size();i2++){
	if( DivByNevt && Nevt>0){
	  *table << v->at(i0)[i1][i2] / Nevt << endl;
	}else{
	  *table << v->at(i0)[i1][i2] << endl;
	}
	nn++;
      }
    }
  }
  return nn;
}

//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::WriteTable(vector<vector<double > >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    for(unsigned int i1=0;i1<v->at(i0).size();i1++){
      if( DivByNevt && Nevt>0){
	*table << v->at(i0)[i1] / Nevt << endl;
      }else{
	*table << v->at(i0)[i1] << endl;
      }
      nn++;
    }
  }
  return nn;
}

*/

//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::WriteTable(vector<double >* v, ostream *table , bool DivByNevt , int Nevt ){
  int nn = 0;
  for(unsigned int i0=0;i0<v->size();i0++){
    if( DivByNevt && Nevt>0){
      *table << v->at(i0) / Nevt << endl;
    }else{
      *table << v->at(i0) << endl;
    }
    nn++;
  }
  return nn;
}



//________________________________________________________________________________________________________________ //
int fastNLOCoeffBase::WriteFlexibleTable(vector<double >* v, ostream *table , bool DivByNevt , int Nevt , bool nProcLast ){
   int nn = 1;
   if ( !nProcLast )*table << v->size() << endl;
   for(unsigned int i0=0;i0<v->size();i0++){
      if( DivByNevt && Nevt>0)	*table << v->at(i0) / Nevt << endl;
      else			*table << v->at(i0) << endl;
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
