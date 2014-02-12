#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOTools.h"
#include "fastnlotk/fastNLOCoeffAddFlex.h"

using namespace std;


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
fastNLOCoeffBase* fastNLOCoeffAddFlex::Clone() const {
   //! Use has to take care to delete this object later
   return static_cast<fastNLOCoeffBase*>(new fastNLOCoeffAddFlex(*this));
}


///________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Read(istream& table){
   fastNLOCoeffBase::ReadBase(table);
   ReadRest(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::ReadRest(istream& table){
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::ReadCoeffAddBase(table);
   ReadCoeffAddFlex(table);
   EndReadCoeff(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::ReadCoeffAddFlex(istream& table){
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

   nn3 += fastNLOTools::ReadFlexibleVector  ( ScaleNode1 , table );
   nn3 += fastNLOTools::ReadFlexibleVector  ( ScaleNode2 , table );
   //NscalenodeScale1 = ScaleNode1[0].size();
   //NscalenodeScale2 = ScaleNode2[0].size();

   nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaTildeMuIndep , table , NSubproc , Nevt );
   //if ( NScaleDep==3 || fScen->ILOord!=Npow || NScaleDep==5 ){
   if ( NScaleDep==3 || NScaleDep>=5 ){
      nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaTildeMuFDep , table , NSubproc , Nevt );
      nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaTildeMuRDep , table , NSubproc , Nevt );
      if ( NScaleDep>=6 ){
         nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaTildeMuRRDep , table , NSubproc , Nevt );
         nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaTildeMuFFDep , table , NSubproc , Nevt );
         nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaTildeMuRFDep , table , NSubproc , Nevt );
      }
   }
   // fixing old convention
   if ( NScaleDep == 3 ) {
      info["ReadCoeffAddFlex"]<<"This is a table with a deprecated convention (NScaleDep=3). Fixing it."<<endl;
      if (Npow!=fILOord) NScaleDep = 5;
      else NScaleDep = 3;
   }
   nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaRefMixed , table , NSubproc , Nevt );
   nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaRef_s1 , table , NSubproc , Nevt );
   nn3 += fastNLOTools::ReadFlexibleVector  ( SigmaRef_s2 , table , NSubproc , Nevt );
   debug["ReadCoeffAddFlex"]<<"Read "<<nn3<<" lines of flexible-scale coefficients."<<endl;

   // init table for evaluation
   fastNLOTools::ResizeFlexibleVector( PdfLcMuVar , SigmaTildeMuIndep );
   AlphasTwoPi.resize(ScaleNode1.size());
   for (unsigned int i=0; i<AlphasTwoPi.size() ; i++) {
      AlphasTwoPi[i].resize(ScaleNode1[i].size());
      for (unsigned int j=0; j<AlphasTwoPi[i].size() ; j++) {
         AlphasTwoPi[i][j].resize(ScaleNode2[i].size());
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Write(ostream& table) {
   CheckCoeffConstants(this);
   // update to latest version
   if ( NScaleDep==3 ) {
      if ( Npow==fILOord) {
	 debug["Write"]<<" * Increase NScaleDep from 3 to 4, because LO!"<<endl;
	 NScaleDep=4;
      }
      else if ( Npow==fILOord+1 ) {
	 debug["Write"]<<" * Increase NScaleDep from 3 to 5 because NLO!"<<endl;
	 NScaleDep=5;
      }
      else if ( Npow==fILOord+2 ) {
	 debug["Write"]<<" * Increase NScaleDep from 3 to 6 because NNLO!"<<endl;
	  NScaleDep=6;
      }
   }
   fastNLOCoeffAddBase::Write(table);

   int nn3 = 0;
   nn3 += fastNLOTools::WriteFlexibleVector( ScaleNode1 , table );
   nn3 += fastNLOTools::WriteFlexibleVector( ScaleNode2 , table );

   nn3 += fastNLOTools::WriteFlexibleVector( SigmaTildeMuIndep, table , NSubproc , Nevt);
   if ( NScaleDep==3 || NScaleDep>=5) {
      //cout<<"Write NLO FlexTable. NScaleDep="<<NScaleDep<<"\tNpow="<<Npow<<"\tfScen->ILOord="<<fScen->ILOord<<endl;
      nn3 += fastNLOTools::WriteFlexibleVector( SigmaTildeMuFDep , table , NSubproc, Nevt);
      nn3 += fastNLOTools::WriteFlexibleVector( SigmaTildeMuRDep , table , NSubproc, Nevt);
      if ( NScaleDep>=6) {
         nn3 += fastNLOTools::WriteFlexibleVector( SigmaTildeMuRRDep , table , NSubproc, Nevt);
         nn3 += fastNLOTools::WriteFlexibleVector( SigmaTildeMuFFDep , table , NSubproc, Nevt);
         nn3 += fastNLOTools::WriteFlexibleVector( SigmaTildeMuRFDep , table , NSubproc, Nevt);
      }
   }
   if ( SigmaRefMixed.empty() ) fastNLOTools::ResizeVector(SigmaRefMixed,fNObsBins,NSubproc);
   if ( SigmaRef_s1.empty() )   fastNLOTools::ResizeVector(SigmaRef_s1,fNObsBins,NSubproc);
   if ( SigmaRef_s2.empty() )   fastNLOTools::ResizeVector(SigmaRef_s2,fNObsBins,NSubproc);
   nn3 += fastNLOTools::WriteFlexibleVector( SigmaRefMixed      , table , NSubproc, Nevt);
   nn3 += fastNLOTools::WriteFlexibleVector( SigmaRef_s1        , table , NSubproc, Nevt);
   nn3 += fastNLOTools::WriteFlexibleVector( SigmaRef_s2        , table , NSubproc, Nevt);

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
   nn3 += WriteFlexibleTable( &SigmaRefMixed    , table , (bool)(option & DividebyNevt) , Nevt , true );
   nn3 += WriteFlexibleTable( &SigmaRef_s1      , table , (bool)(option & DividebyNevt) , Nevt , true );
   nn3 += WriteFlexibleTable( &SigmaRef_s2      , table , (bool)(option & DividebyNevt) , Nevt , true );
   */
   //printf("  *  fastNLOCoeffAddFlex::Write(). Wrote %d lines of v2.1 Tables.\n",nn3);
   debug["Write"]<<"Wrote "<<nn3<<" lines of v2.1 Tables."<<endl;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Add(const fastNLOCoeffAddBase& other){
   bool ok = CheckCoeffConstants(this);
   if ( !ok ) {
      error["Add"]<<"Incompatible table."<<endl;
   }
   const fastNLOCoeffAddFlex& othflex = (const fastNLOCoeffAddFlex&) other;
   Nevt += othflex.Nevt;
   fastNLOTools::AddVectors( SigmaTildeMuIndep , othflex.SigmaTildeMuIndep );
   if ( NScaleDep==3 || NScaleDep>=5 ) {
      fastNLOTools::AddVectors( SigmaTildeMuFDep , othflex.SigmaTildeMuFDep );
      fastNLOTools::AddVectors( SigmaTildeMuRDep , othflex.SigmaTildeMuRDep );
      if ( NScaleDep>=6 ) {
         fastNLOTools::AddVectors( SigmaTildeMuRRDep , othflex.SigmaTildeMuRRDep );
         fastNLOTools::AddVectors( SigmaTildeMuFFDep , othflex.SigmaTildeMuFFDep );
         fastNLOTools::AddVectors( SigmaTildeMuRFDep , othflex.SigmaTildeMuRFDep );
      }
   }
   fastNLOTools::AddVectors( SigmaRefMixed , othflex.SigmaRefMixed );
   fastNLOTools::AddVectors( SigmaRef_s1 , othflex.SigmaRef_s1 );
   fastNLOTools::AddVectors( SigmaRef_s2 , othflex.SigmaRef_s2 );
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Clear() {
   //! Set all elements of sigma tilde to zero
   fastNLOCoeffAddBase::Clear();
   fastNLOTools::ClearVector(SigmaTildeMuIndep);
   fastNLOTools::ClearVector(SigmaTildeMuFDep); 
   fastNLOTools::ClearVector(SigmaTildeMuRDep); 
   fastNLOTools::ClearVector(SigmaTildeMuRRDep); 
   fastNLOTools::ClearVector(SigmaTildeMuFFDep); 
   fastNLOTools::ClearVector(SigmaTildeMuRFDep); 
   fastNLOTools::ClearVector(SigmaRefMixed);
   fastNLOTools::ClearVector(SigmaRef_s1); 
   fastNLOTools::ClearVector(SigmaRef_s2); 
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Print() const {
   fastNLOCoeffAddBase::Print();
   printf(" **************** FastNLO Table: fastNLOCoeffAddFlex ****************\n");
   printf(" B   NscalenodeScale1              %lu\n",ScaleNode1[0].size());
   printf(" B   NscalenodeScale2              %lu\n",ScaleNode2[0].size());
   printf(" *******************************************************\n");
}


//________________________________________________________________________________________________________________ //
