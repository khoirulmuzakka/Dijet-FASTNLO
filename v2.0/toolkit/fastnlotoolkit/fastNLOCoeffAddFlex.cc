#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOTools.h"
#include "fastnlotk/fastNLOCoeffAddFlex.h"

using namespace std;
using namespace fastNLO;


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffAddFlex::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet ) {
   bool ret = fastNLOCoeffAddBase::CheckCoeffConstants(c,quiet);
   if ( ret &&  c->GetNScaleDep() >= 3) return true;
   else if ( c->GetNScaleDep() < 3 ) {
      if ( !quiet )
	 say::error["CheckCoeffConstants"]<<"This is not a flexible scale table. NScaleDep must be >= 3 but is NScaleDep="
					  <<c->GetNScaleDep()<<endl;
      return false;
   }
   else return false;
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFlex::fastNLOCoeffAddFlex(){
   SetClassName("fastNLOCoeffAddFlex");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFlex::fastNLOCoeffAddFlex(int NObsBin, int iLOord) : fastNLOCoeffAddBase(NObsBin){
   SetClassName("fastNLOCoeffAddFlex");
   fILOord = iLOord; // only necessary for fixing NScaleDep 3 -> 4,5
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFlex::fastNLOCoeffAddFlex(const fastNLOCoeffBase& base , int iLOord ) : fastNLOCoeffAddBase(base)  {
   SetClassName("fastNLOCoeffAddFlex");
   fILOord = iLOord;
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddFlex::Read(istream *table){
   fastNLOCoeffBase::ReadBase(table);
   ReadRest(table);
   return 0;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::ReadRest(istream *table){
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::ReadCoeffAddBase(table);
   ReadCoeffAddFlex(table);
   EndReadCoeff(table);
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddFlex::ReadCoeffAddFlex(istream *table){
   CheckCoeffConstants(this);

   //  ---- order of reading... ---- //
   //    - nscalenode q2
   //    - scalenode Q
   //    - nscalenode pt
   //    - scalenode pt
   //    - simgatilde mu indep
   //    - simgatilde mu_f dep
   //    - simgatilde mu_r dep
   //    - sigmarefmixed
   //    - sigmaref scale 1
   //    - sigmaref scale 2
   // ------------------------------ //
   int nn3 = 0;

   nn3 += ReadFlexibleVector  ( &ScaleNode1 , table );
   nn3 += ReadFlexibleVector  ( &ScaleNode2 , table );
   //NscalenodeScale1 = ScaleNode1[0].size();
   //NscalenodeScale2 = ScaleNode2[0].size();

   nn3 += ReadFlexibleVector  ( &SigmaTildeMuIndep , table , true , Nevt );
   //if ( NScaleDep==3 || fScen->ILOord!=Npow || NScaleDep==5 ){
   if ( NScaleDep==3 || NScaleDep>=5 ){
      nn3 += ReadFlexibleVector  ( &SigmaTildeMuFDep , table , true , Nevt );
      nn3 += ReadFlexibleVector  ( &SigmaTildeMuRDep , table , true , Nevt );
      if ( NScaleDep>=6 ){
	 nn3 += ReadFlexibleVector  ( &SigmaTildeMuRRDep , table , true , Nevt );
	 nn3 += ReadFlexibleVector  ( &SigmaTildeMuFFDep , table , true , Nevt );
	 nn3 += ReadFlexibleVector  ( &SigmaTildeMuRFDep , table , true , Nevt );
      }
   }
   // fixing old convention
   if ( NScaleDep == 3 ) {
      info["ReadCoeffAddFlex"]<<"This is a table with a deprecated convention (NScaleDep=3). Fixing it."<<endl;
      if (Npow!=fILOord) NScaleDep = 5;
      else NScaleDep = 3;
   }
   nn3 += ReadFlexibleVector  ( &SigmaRefMixed , table , true , Nevt );
   nn3 += ReadFlexibleVector  ( &SigmaRef_s1 , table , true , Nevt );
   nn3 += ReadFlexibleVector  ( &SigmaRef_s2 , table , true , Nevt );
   printf("  *  fastNLOCoeffAddFlex::Read(). Read %d lines of flexible-scale tables.\n",nn3);

   // init table for evaluation
   ResizeFlexibleVector( PdfLcMuVar , SigmaTildeMuIndep );
   AlphasTwoPi.resize(ScaleNode1.size());
   for (unsigned int i=0; i<AlphasTwoPi.size() ; i++) {
      AlphasTwoPi[i].resize(ScaleNode1[i].size());
      for (unsigned int j=0; j<AlphasTwoPi[i].size() ; j++) {
	 AlphasTwoPi[i][j].resize(ScaleNode2[i].size());
      }
   }

   return 0;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Write(ostream *table) {
   CheckCoeffConstants(this);
   // update to latest version
    if ( NScaleDep==3 ) {
       if ( Npow==fILOord) {
 	 info["Write"]<<" * Increase NScaleDep from 3 to 4, because LO!"<<endl;
 	 NScaleDep=4;
       }
       else if ( Npow==fILOord+1 ) {
 	 info["Write"]<<" * Increase NScaleDep from 3 to 5 because NLO!"<<endl;
 	 NScaleDep=5;
       }
       else if ( Npow==fILOord+2 ) {
 	 info["Write"]<<" * Increase NScaleDep from 3 to 6 because NNLO!"<<endl;
 	 NScaleDep=6;
       }
    }
   fastNLOCoeffAddBase::Write(table);

   int nn3 = 0;
   nn3 += WriteFlexibleVector( ScaleNode1 , table );
   nn3 += WriteFlexibleVector( ScaleNode2 , table );

   nn3 += WriteFlexibleVector( SigmaTildeMuIndep, table , NSubproc , Nevt);
   if ( NScaleDep==3 || NScaleDep>=5) {
      //cout<<"Write NLO FlexTable. NScaleDep="<<NScaleDep<<"\tNpow="<<Npow<<"\tfScen->ILOord="<<fScen->ILOord<<endl;
      nn3 += WriteFlexibleVector( SigmaTildeMuFDep , table , NSubproc, Nevt);
      nn3 += WriteFlexibleVector( SigmaTildeMuRDep , table , NSubproc, Nevt);
      if ( NScaleDep>=6) {
	 nn3 += WriteFlexibleVector( SigmaTildeMuRRDep , table , NSubproc, Nevt);
	 nn3 += WriteFlexibleVector( SigmaTildeMuFFDep , table , NSubproc, Nevt);
	 nn3 += WriteFlexibleVector( SigmaTildeMuRFDep , table , NSubproc, Nevt);
      }
   }
   if ( SigmaRefMixed.empty() ) fastNLO::ResizeVector(SigmaRefMixed,fNObsBins,NSubproc);
   if ( SigmaRef_s1.empty() )   fastNLO::ResizeVector(SigmaRef_s1,fNObsBins,NSubproc);
   if ( SigmaRef_s2.empty() )   fastNLO::ResizeVector(SigmaRef_s2,fNObsBins,NSubproc);
   nn3 += WriteFlexibleVector( SigmaRefMixed	, table , NSubproc, Nevt);
   nn3 += WriteFlexibleVector( SigmaRef_s1	, table , NSubproc, Nevt);
   nn3 += WriteFlexibleVector( SigmaRef_s2	, table , NSubproc, Nevt);

   /*
   nn3 += WriteFlexibleTable( &SigmaTildeMuIndep, table , (bool)(option & DividebyNevt) , Nevt , true );

   //if ( NScaleDep==3 || Npow!=fScen->ILOord || NScaleDep==5) {
   if ( NScaleDep==3 || NScaleDep>=5) {
      //cout<<"Write NLO FlexTable. NScaleDep="<<NScaleDep<<"\tNpow="<<Npow<<"\tfScen->ILOord="<<fScen->ILOord<<endl;
      nn3 += WriteFlexibleTable( &SigmaTildeMuFDep , table , (bool)(option & DividebyNevt) , Nevt , true );
      nn3 += WriteFlexibleTable( &SigmaTildeMuRDep , table , (bool)(option & DividebyNevt) , Nevt , true );
      if ( NScaleDep>=6) {
	 nn3 += WriteFlexibleTable( &SigmaTildeMuRRDep , table , (bool)(option & DividebyNevt) , Nevt , true );
	 nn3 += WriteFlexibleTable( &SigmaTildeMuFFDep , table , (bool)(option & DividebyNevt) , Nevt , true );
	 nn3 += WriteFlexibleTable( &SigmaTildeMuRFDep , table , (bool)(option & DividebyNevt) , Nevt , true );
      }
   }
   if ( SigmaRefMixed.empty() ) fastNLOCoeffBase::ResizeTable(&SigmaRefMixed,fNObsBins,NSubproc);
   if ( SigmaRef_s1.empty() )   fastNLOCoeffBase::ResizeTable(&SigmaRef_s1,fNObsBins,NSubproc);
   if ( SigmaRef_s2.empty() )   fastNLOCoeffBase::ResizeTable(&SigmaRef_s2,fNObsBins,NSubproc);
   nn3 += WriteFlexibleTable( &SigmaRefMixed	, table , (bool)(option & DividebyNevt) , Nevt , true );
   nn3 += WriteFlexibleTable( &SigmaRef_s1	, table , (bool)(option & DividebyNevt) , Nevt , true );
   nn3 += WriteFlexibleTable( &SigmaRef_s2	, table , (bool)(option & DividebyNevt) , Nevt , true );
   */
   printf("  *  fastNLOCoeffAddFlex::Write(). Wrote %d lines of v2.1 Tables.\n",nn3);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Add(const fastNLOCoeffAddFlex& other){
   cout<<"Adding! fastNLOCoeffAddFlex::Add"<<endl;
   CheckCoeffConstants(this);
   Nevt += other.Nevt;
   
   /*
   double w1 = (double)Nevt / (Nevt+other.Nevt);
   double w2 = (double)other.Nevt / (Nevt+other.Nevt);
   AddTableToAnotherTable( SigmaTildeMuIndep , other.SigmaTildeMuIndep ,w1 , w2 );
   if ( NScaleDep==3 || NScaleDep>=5 ) {
      AddTableToAnotherTable( SigmaTildeMuFDep , other.SigmaTildeMuFDep ,w1 , w2 );
      AddTableToAnotherTable( SigmaTildeMuRDep , other.SigmaTildeMuRDep ,w1 , w2 );
      if ( NScaleDep>=6 ) {
	 AddTableToAnotherTable( SigmaTildeMuRRDep , other.SigmaTildeMuRRDep ,w1 , w2 );
	 AddTableToAnotherTable( SigmaTildeMuFFDep , other.SigmaTildeMuFFDep ,w1 , w2 );
	 AddTableToAnotherTable( SigmaTildeMuRFDep , other.SigmaTildeMuRFDep ,w1 , w2 );
      }
   }
   AddTableToAnotherTable( SigmaRefMixed , other.SigmaRefMixed ,w1 , w2 );
   AddTableToAnotherTable( SigmaRef_s1 , other.SigmaRef_s1 ,w1 , w2 );
   AddTableToAnotherTable( SigmaRef_s2 , other.SigmaRef_s2 ,w1 , w2 );
   */

   AddVectors( SigmaTildeMuIndep , other.SigmaTildeMuIndep );
   if ( NScaleDep==3 || NScaleDep>=5 ) {
      AddVectors( SigmaTildeMuFDep , other.SigmaTildeMuFDep );
      AddVectors( SigmaTildeMuRDep , other.SigmaTildeMuRDep );
      if ( NScaleDep>=6 ) {
	 AddVectors( SigmaTildeMuRRDep , other.SigmaTildeMuRRDep );
	 AddVectors( SigmaTildeMuFFDep , other.SigmaTildeMuFFDep );
	 AddVectors( SigmaTildeMuRFDep , other.SigmaTildeMuRFDep );
      }
   }
   AddVectors( SigmaRefMixed , other.SigmaRefMixed );
   AddVectors( SigmaRef_s1 , other.SigmaRef_s1 );
   AddVectors( SigmaRef_s2 , other.SigmaRef_s2 );
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddFlex::ReadFlexibleVector(vector<double >* v, istream *table , bool nProcLast , unsigned long long nevts ){
   int nn = 0;
   if ( !nProcLast ) {
      int size = 0;
      *table >> size; nn++;
      v->resize(size);
   }
   else {
      v->resize(NSubproc);
   }
   for(unsigned int i0=0;i0<v->size();i0++){
      *table >> v->at(i0);
      v->at(i0) *= nevts;
      nn++;
   }
   return nn;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Print() const {
   fastNLOCoeffAddBase::Print();
   printf(" **************** FastNLO Table: fastNLOCoeffAddFlex ****************\n");
   printf(" B   NscalenodeScale1              %d\n",ScaleNode1[0].size());
   printf(" B   NscalenodeScale2              %d\n",ScaleNode2[0].size());
   printf(" *******************************************************\n");
}


//________________________________________________________________________________________________________________ //
