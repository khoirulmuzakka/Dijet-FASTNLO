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
   //! User has to take care to delete this object later
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
      }
      if ( NScaleDep>=7 ){
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
         NScaleDep=5;
         if ( !fastNLOTools::IsEmptyVector(SigmaTildeMuRRDep) ) {
            NScaleDep=6;
            debug["Write"]<<" * Increase NScaleDep from 3 to 6 because NLO and log^2(mur) terms!"<<endl;
         }
         else
            debug["Write"]<<" * Increase NScaleDep from 3 to 5 because NLO!"<<endl;

      }
      else if ( Npow==fILOord+2 ) {
         debug["Write"]<<" * Increase NScaleDep from 3 to 6 because NNLO!"<<endl;
         NScaleDep=7;
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
      }
      if ( NScaleDep>=7) {
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
      if (( NScaleDep>=6 || !SigmaTildeMuRRDep.empty())  // both tables contain log^2 contributions (default case)
          && (othflex.NScaleDep>=6 || !othflex.SigmaTildeMuRRDep.empty()) ) {
         fastNLOTools::AddVectors( SigmaTildeMuRRDep , othflex.SigmaTildeMuRRDep );
      }
      else if ( NScaleDep==6 && othflex.NScaleDep==5 ) { // this tables contains log^2 contributions, but the other does not
         // nothing todo.
      }
      else if ( NScaleDep==5 && othflex.NScaleDep==6 ) { // this tables does not contain log^2 contributions, but the other does !
         SigmaTildeMuRRDep = othflex.SigmaTildeMuRRDep;
         NScaleDep = 6;
      }
      if ( NScaleDep>=7 || !SigmaTildeMuFFDep.empty()) {
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
bool  fastNLOCoeffAddFlex::IsCompatible(const fastNLOCoeffAddFlex& other) const {
   //! Check for compatibility for merging/adding of two contributions
   if ( ! ((fastNLOCoeffAddBase*)this)->IsCompatible(other)) return false;
   for ( int i=0 ; i<fNObsBins ; i++ ){
      if ( GetNScaleNode1(i) != other.GetNScaleNode1(i) ) {
         say::warn["fastNLOCoeffAddFlex::IsCompatible"]<<"Incompatible number of scale nodes found."<<endl;
         return false;
      }
      if ( GetNScaleNode2(i) != other.GetNScaleNode2(i) ) {
         say::warn["fastNLOCoeffAddFlex::IsCompatible"]<<"Incompatible number of scale nodes found."<<endl;
         return false;
      }
      for ( unsigned int is1 = 0 ; is1<GetNScaleNode1(i) ; is1++ ) {
         if ( GetScaleNode1(i,is1) != other.GetScaleNode1(i,is1) ) {
            say::warn["fastNLOCoeffAddFlex::IsCompatible"]<<"Incompatible scale1 node found."<<endl;
            return false;
         }
      }
      for ( unsigned int is2 = 0 ; is2<GetNScaleNode2(i) ; is2++ ) {
         if ( GetScaleNode2(i,is2) != other.GetScaleNode2(i,is2) ) {
            say::warn["fastNLOCoeffAddFlex::IsCompatible"]<<"Incompatible scale2 node found."<<endl;
            return false;
         }
      }
   }
   return true;
}


//________________________________________________________________________________________________________________ //
bool  fastNLOCoeffAddFlex::IsCatenable(const fastNLOCoeffAddFlex& other) const {
   //! Check for compatibility of catenating observable bins
   if ( ! ((fastNLOCoeffAddBase*)this)->IsCatenable(other)) return false;
   // TODO: More checks
   info["IsCatenable"]<<"Flex-scale contributions are catenable"<<endl;
   return true;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::NormalizeCoefficients(){
   //!< Set number of events to 1 and normalize coefficients accordingly.
   //! This means, that the information about the
   //! number of events is essentially lost
   MultiplyCoefficientsByConstant(1./Nevt);
   Nevt = 1;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::MultiplyCoefficientsByConstant(double fact) {
   for (unsigned int i=0; i<SigmaTildeMuIndep.size(); i++) {
      MultiplyBin(i,fact);
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFlex::Print(int iprint) const {
   if ( !(iprint < 0) ) {
      fastNLOCoeffAddBase::Print(iprint);
      cout << fastNLO::_DSEP20C << " fastNLO Table: CoeffAddFlex " << fastNLO::_DSEP20 << endl;
   } else {
      cout << endl << fastNLO::_CSEP20C << " fastNLO Table: CoeffAddFlex " << fastNLO::_CSEP20 << endl;
   }
   printf(" # No. of scale1 nodes (Nscalenode1)   %d\n",(int)ScaleNode1[0].size());
   printf(" # No. of scale2 nodes (Nscalenode2)   %d\n",(int)ScaleNode2[0].size());
   if ( abs(iprint) > 0 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 0) " << fastNLO::_SSEP20 << endl;
      char buffer[1024];
      for (int i=0; i<fNObsBins; i++) {
         // Print only for first and last observable bin
         if (i==0 || i==fNObsBins-1) {
            printf(" #   Observable bin no. %d\n",i+1);
            snprintf(buffer, sizeof(buffer), "Scale nodes 1 (ScaleNode1[%d][])",i);
            fastNLOTools::PrintVector(GetScaleNodes1(i),buffer,"#  ");
            if ( ScaleNode2.size() != 0 ) {
               snprintf(buffer, sizeof(buffer), "Scale nodes 2 (ScaleNode2[%d][])",i);
               fastNLOTools::PrintVector(GetScaleNodes2(i),buffer,"#  ");
            }
         }
      }
   }
   cout << fastNLO::_CSEPSC << endl;
}


//________________________________________________________________________________________________________________ //

// Erase observable bin
void fastNLOCoeffAddFlex::EraseBin(unsigned int iObsIdx) {
   debug["fastNLOCoeffAddFlex::EraseBin"]<<"Erasing table entries in CoeffAddFlex for bin index " << iObsIdx << endl;
   if ( ScaleNode1.size() != 0 ) {ScaleNode1.erase(ScaleNode1.begin()+iObsIdx);}
   if ( ScaleNode2.size() != 0 ) {ScaleNode2.erase(ScaleNode2.begin()+iObsIdx);}
   if ( SigmaTildeMuIndep.size() != 0 ) {SigmaTildeMuIndep.erase(SigmaTildeMuIndep.begin()+iObsIdx);}
   if ( SigmaTildeMuFDep.size()  != 0 ) {SigmaTildeMuFDep.erase(SigmaTildeMuFDep.begin()+iObsIdx);}
   if ( SigmaTildeMuRDep.size()  != 0 ) {SigmaTildeMuRDep.erase(SigmaTildeMuRDep.begin()+iObsIdx);}
   if ( SigmaTildeMuRRDep.size() != 0 ) {SigmaTildeMuRRDep.erase(SigmaTildeMuRRDep.begin()+iObsIdx);}
   if ( SigmaTildeMuFFDep.size() != 0 ) {SigmaTildeMuFFDep.erase(SigmaTildeMuFFDep.begin()+iObsIdx);}
   if ( SigmaTildeMuRFDep.size() != 0 ) {SigmaTildeMuRFDep.erase(SigmaTildeMuRFDep.begin()+iObsIdx);}
   fastNLOCoeffAddBase::EraseBin(iObsIdx);
}

// Catenate observable bin
void fastNLOCoeffAddFlex::CatBin(const fastNLOCoeffAddFlex& other, unsigned int iObsIdx) {
   debug["fastNLOCoeffAddFlex::CatBin"]<<"Catenating observable bin in CoeffAddFlex corresponding to bin index " << iObsIdx << endl;
   if ( ScaleNode1.size() == 0 ) {
      say::error["CatBin"]<<"Scale node1 size cannot be zero for a flexible-scale table. Aborted!" << endl;
      exit(1);
   }
   unsigned int nold = ScaleNode1.size();
   ScaleNode1.resize(nold+1);
   ScaleNode1[nold] = other.ScaleNode1[iObsIdx];
   if ( ScaleNode2.size() != 0 ) {
      ScaleNode2.resize(nold+1);
      ScaleNode2[nold] = other.ScaleNode2[iObsIdx];
   }
   if ( SigmaTildeMuIndep.size() != 0 ) {
      SigmaTildeMuIndep.resize(nold+1);
      SigmaTildeMuIndep[nold] = other.SigmaTildeMuIndep[iObsIdx];
   }
   if ( SigmaTildeMuFDep.size() != 0 ) {
      SigmaTildeMuFDep.resize(nold+1);
      SigmaTildeMuFDep[nold] = other.SigmaTildeMuFDep[iObsIdx];
   }
   if ( SigmaTildeMuRDep.size() != 0 ) {
      SigmaTildeMuRDep.resize(nold+1);
      SigmaTildeMuRDep[nold] = other.SigmaTildeMuRDep[iObsIdx];
   }
   if ( SigmaTildeMuRRDep.size() != 0 ) {
      SigmaTildeMuRRDep.resize(nold+1);
      SigmaTildeMuRRDep[nold] = other.SigmaTildeMuRRDep[iObsIdx];
   }
   if ( SigmaTildeMuFFDep.size() != 0 ) {
      SigmaTildeMuFFDep.resize(nold+1);
      SigmaTildeMuFFDep[nold] = other.SigmaTildeMuFFDep[iObsIdx];
   }
   if ( SigmaTildeMuRFDep.size() != 0 ) {
      SigmaTildeMuRFDep.resize(nold+1);
      SigmaTildeMuRFDep[nold] = other.SigmaTildeMuRFDep[iObsIdx];
   }
   fastNLOCoeffAddBase::CatBin(other, iObsIdx);
}

// Multiply observable bin
void fastNLOCoeffAddFlex::MultiplyBin(unsigned int iObsIdx, double fact) {
   debug["fastNLOCoeffAddFlex::MultiplyBin"]<<"Multiplying table entries in CoeffAddFlex for bin index " << iObsIdx << " by factor " << fact << endl;
   const bool WithMuR  = !SigmaTildeMuRDep.empty();
   const bool WithMuF  = !SigmaTildeMuFDep.empty();
   const bool WithMuRR = !SigmaTildeMuRRDep.empty();
   const bool WithMuFF = !SigmaTildeMuFFDep.empty();
   const bool WithMuRF = !SigmaTildeMuRFDep.empty();
   int nxmax = GetNxmax(iObsIdx);
   for (unsigned int jS1=0; jS1<GetNScaleNode1(iObsIdx); jS1++) {
      for (unsigned int kS2=0; kS2<GetNScaleNode2(iObsIdx); kS2++) {
         for (int x=0; x<nxmax; x++) {
            for (int n=0; n<GetNSubproc(); n++) {
               SigmaTildeMuIndep[iObsIdx][x][jS1][kS2][n] *= fact;
               if (WithMuR)  {SigmaTildeMuRDep[iObsIdx][x][jS1][kS2][n]  *= fact;}
               if (WithMuF)  {SigmaTildeMuFDep[iObsIdx][x][jS1][kS2][n]  *= fact;}
               if (WithMuRR) {SigmaTildeMuRRDep[iObsIdx][x][jS1][kS2][n] *= fact;}
               if (WithMuFF) {SigmaTildeMuFFDep[iObsIdx][x][jS1][kS2][n] *= fact;}
               if (WithMuRF) {SigmaTildeMuRFDep[iObsIdx][x][jS1][kS2][n] *= fact;}
            }
         }
      }
   }
   fastNLOCoeffAddBase::MultiplyBin(iObsIdx, fact);
}
