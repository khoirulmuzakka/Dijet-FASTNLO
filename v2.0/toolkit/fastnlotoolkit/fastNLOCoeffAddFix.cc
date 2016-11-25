#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffAddFix.h"
#include "fastnlotk/fastNLOTools.h"
#include "fastnlotk/speaker.h"

using namespace std;
using namespace fastNLO;


//________________________________________________________________________________________________________________ //
bool fastNLOCoeffAddFix::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet) {
   bool ret = fastNLOCoeffAddBase::CheckCoeffConstants(c,quiet);
   if ( ret && c->GetNScaleDep() == 0 ) return true;
   else if ( c->GetNScaleDep() >= 3 ) {
      if ( !quiet)
         say::error["fastNLOCoeffAddFix::CheckCoeffConstants"]
            <<"This is not a fixed order v2.0  table. NScaleDep must be equal 0 but is NScaleDep="
            <<c->GetNScaleDep()<<endl;
      return false;
   }
   else return false;
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFix::fastNLOCoeffAddFix(int NObsBin) : fastNLOCoeffAddBase(NObsBin) {
   SetClassName("fastNLOCoeffAddFix");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffAddFix::fastNLOCoeffAddFix(const fastNLOCoeffBase& base) : fastNLOCoeffAddBase(base) {
   SetClassName("fastNLOCoeffAddFix");
}


//________________________________________________________________________________________________________________ //
fastNLOCoeffBase* fastNLOCoeffAddFix::Clone() const {
   //! Use has to take care to delete this object later
   return static_cast<fastNLOCoeffBase*>(new fastNLOCoeffAddFix(*this));
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Read(istream& table){
   fastNLOCoeffBase::ReadBase(table);
   ReadRest(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ReadRest(istream& table){
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::ReadCoeffAddBase(table);
   ReadCoeffAddFix(table);
   EndReadCoeff(table);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ReadCoeffAddFix(istream& table){
   CheckCoeffConstants(this);

   Nscalevar.resize(NScaleDim);
   vector<int> Nscalenode(NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      table >> Nscalevar[i];
      table >> Nscalenode[i];
   }
   //    printf("  *  fastNLOCoeffAddFix::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d,  NScaleDim %d  \n",
   //    fNObsBins, Nscalevar[0] , Nscalenode[0] , NScaleDim );
   ScaleFac.resize(NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      ScaleFac[i].resize(Nscalevar[i]);
      for(int j=0;j<Nscalevar[i];j++){
         table >> ScaleFac[i][j];
      }
   }
   //printf("  *  fastNLOCoeffAddFix::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d, ScaleFac[0][0] %d,  NScaleDim %d  \n",
   //fNObsBins, Nscalevar[0] , Nscalenode[0] , ScaleFac[0][0], NScaleDim );
   fastNLOTools::ResizeVector( ScaleNode , fNObsBins, 1 , Nscalevar[0] , Nscalenode[0] ); // should work, since NScaleDim==1, but is not yet tested for 100%
   int nsn = fastNLOTools::ReadVector( ScaleNode , table );
   debug["fastNLOCoeffAddFix::Read()"]<<"Read "<<nsn<<" lines of ScaleNode."<<endl;

   ResizeSigmaTilde();
   ResizePdfLC();
   ResizePdfSplLC();
   int nst = fastNLOTools::ReadVector( SigmaTilde , table , Nevt);
   debug["fastNLOCoeffAddFix::Read()"]<<"Read "<<nsn+nst<<" lines of fastNLO v2 tables."<<endl;

   // prepare members for evaluation
   fastNLOTools::ResizeVector(AlphasTwoPi_v20 , fNObsBins, GetTotalScalenodes());
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ResizeSigmaTilde(){
   //! resize SigmaTilde and ensure that all entries are empty
   SigmaTilde.resize(fNObsBins);
   for( int i=0 ; i<fNObsBins ; i++ ){
      int nxmax = GetNxmax(i);
      SigmaTilde[i].resize(GetTotalScalevars());
      for( int k=0 ; k<GetTotalScalevars() ; k++ ){
         SigmaTilde[i][k].resize(GetTotalScalenodes());
         for( int l=0 ; l<GetTotalScalenodes() ; l++ ){
            //ResizeVector(SigmaTilde[i][k][l],nxmax,NSubproc);
            SigmaTilde[i][k][l].resize(nxmax);
            for( int m=0 ; m<nxmax ; m++ ){
               SigmaTilde[i][k][l][m].resize(NSubproc);
               for( int n=0 ; n<NSubproc ; n++ ){
                  SigmaTilde[i][k][l][m][n] = 0.;
               }
            }
         }
      }
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ResizePdfLC(){
   //! resize PdfLC
   PdfLc.resize(fNObsBins);
   for( int i=0 ; i<fNObsBins ; i++ ){
      int nxmax = GetNxmax(i);
      int totalscalenodes = GetTotalScalenodes();
      PdfLc[i].resize(totalscalenodes);
      for( int l=0 ; l<totalscalenodes ; l++ ){
         fastNLOTools::ResizeVector(PdfLc[i][l],nxmax,NSubproc);
      }
   }
}

//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ResizePdfSplLC(){
   //! resize PdfSplLC
   PdfSplLc1.resize(fNObsBins);
   PdfSplLc2.resize(fNObsBins);
   for( int i=0 ; i<fNObsBins ; i++ ){
      int nxmax = GetNxmax(i);
      int totalscalenodes = GetTotalScalenodes();
      PdfSplLc1[i].resize(totalscalenodes);
      PdfSplLc2[i].resize(totalscalenodes);
      for( int l=0 ; l<totalscalenodes ; l++ ){
         fastNLOTools::ResizeVector(PdfSplLc1[i][l],nxmax,NSubproc);
         fastNLOTools::ResizeVector(PdfSplLc2[i][l],nxmax,NSubproc);
      }
   }
}



//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Write(ostream& table){
   //! Write coefficient table to disk (ostream)
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::Write(table);

   for(int i=0;i<NScaleDim;i++){
      table << Nscalevar[i] << sep;
      table << GetNScaleNode() << sep;
   }
   fastNLOTools::WriteVector( ScaleFac , table );
   int nsn = fastNLOTools::WriteVector( ScaleNode , table );
   int nst = fastNLOTools::WriteVector( SigmaTilde , table , Nevt);
   info["Write"]<<"Wrote "<<nst+nsn<<" lines into fastNLO table."<<endl;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Add(const fastNLOCoeffAddBase& other){
   //! Add another coefficient table to this table
   bool ok = CheckCoeffConstants(this);
   if ( !ok ) {
      error["Add"]<<"Cannot add tables."<<endl;
      return;
   }
   const fastNLOCoeffAddFix& othfix = (const fastNLOCoeffAddFix&)other;
   Nevt += othfix.Nevt;
   fastNLOTools::AddVectors( SigmaTilde , othfix.SigmaTilde);
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddFix::GetTotalScalevars() const {
   //! Get nuber of scale-variations
   int totalscalevars=1;
   for(int scaledim=0;scaledim<NScaleDim;scaledim++){
      totalscalevars *= Nscalevar[scaledim];
   }
   return totalscalevars;
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddFix::GetTotalScalenodes() const {
   //! Get number of scale nodes
   if ( !ScaleNode.empty() ) return ScaleNode[0][0][0].size();
   else return 0;
   //    int totalscalenodes=1;
   //    for(int scaledim=0;scaledim<NScaleDim;scaledim++){
   //       if ( !ScaleNode.empty() )
   //  totalscalenodes *= ScaleNode[0][0][0].size();
   //       //Nscalenode[scaledim];
   //    }
   //    return totalscalenodes;
}


//________________________________________________________________________________________________________________ //
bool  fastNLOCoeffAddFix::IsCompatible(const fastNLOCoeffAddFix& other) const {
   //! Check for compatibility for merging/adding of two contributions
   if ( ! ((fastNLOCoeffAddBase*)this)->IsCompatible(other)) return false;
   if ( GetNScaleNode() != other.GetNScaleNode() ) {
      say::warn["IsCompatible"]<<"Incompatible number of scale nodes found."<<endl;
      return false;
   }
   if ( GetNScalevar() != other.GetNScalevar() ) {
      say::warn["IsCompatible"]<<"Incompatible number of scale variations found."<<endl;
      return false;
   }
   if ( GetAvailableScaleFactors()[GetNScalevar()-1] != other.GetAvailableScaleFactors()[GetNScalevar()-1] ) {
      say::warn["IsCompatible"]<<"Incompatible scale variations found."<<endl;
      return false;
   }
   for ( int i=0 ; i<fNObsBins ; i++ ){
      for ( int is=0 ; is<GetNScaleNode() ; is++ ){
         if ( GetScaleNode(i,0,is) != other.GetScaleNode(i,0,is) ) {
            say::warn["IsCompatible"]<<"Incompatible scale node found."<<endl;
            return false;
         }
      }
   }
   return true;
}


//________________________________________________________________________________________________________________ //
bool  fastNLOCoeffAddFix::IsCatenable(const fastNLOCoeffAddFix& other) const {
   //! Check for compatibility of catenating observable bins
   if ( ! ((fastNLOCoeffAddBase*)this)->IsCatenable(other)) return false;
   if ( GetNScaleNode() != other.GetNScaleNode() ) {
      debug["IsCatenable"]<<"Incompatible number of scale nodes found. Skipped."<<endl;
      return false;
   }
   if ( GetNScalevar() != other.GetNScalevar() ) {
      debug["IsCatenable"]<<"Incompatible number of scale variations found. Skipped."<<endl;
      return false;
   }
   if ( GetAvailableScaleFactors()[GetNScalevar()-1] != other.GetAvailableScaleFactors()[GetNScalevar()-1] ) {
      debug["IsCatenable"]<<"Incompatible scale variations found. Skipped."<<endl;
      return false;
   }
   info["IsCatenable"]<<"Fix-scale contributions are catenable"<<endl;
   return true;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Clear() {
   //! Set all elelments of SigmaTilde to zero.
   fastNLOCoeffAddBase::Clear();
   fastNLOTools::ClearVector(SigmaTilde);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::NormalizeCoefficients(){
   //!< Set number of events to 1 and normalize coefficients accordingly.
   //! This means, that the information about the
   //! number of events is essentially lost
   MultiplyCoefficientsByConstant(1./Nevt);
   Nevt = 1;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::MultiplyCoefficientsByConstant(double fact) {
   for (unsigned int i=0; i<SigmaTilde.size(); i++) {
      MultiplyBin(i,fact);
   }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Print(int iprint) const {
   if ( !(iprint < 0) ) {
      fastNLOCoeffAddBase::Print(iprint);
      cout << fastNLO::_DSEP20C << " fastNLO Table: CoeffAddFix " << fastNLO::_DSEP20 << endl;
   } else {
      cout << endl << fastNLO::_CSEP20C << " fastNLO Table: CoeffAddFix " << fastNLO::_CSEP20 << endl;
   }
   for (int i=0; i<NScaleDim; i++) {
      printf(" # No. of scale variations (Nscalevar) %d\n",GetNScalevar());
      fastNLOTools::PrintVector(GetAvailableScaleFactors(),"Scale factors (ScaleFac[0][])","#");
      printf(" # No. of scale nodes (Nscalenode)     %d\n",GetNScaleNode());
   }
   if ( abs(iprint) > 0 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 0) " << fastNLO::_SSEP20 << endl;
      char buffer[1024];
      for (int i=0; i<fNObsBins; i++) {
         // Print only for first and last observable bin
         if (i==0 || i==fNObsBins-1) {
            printf(" #   Observable bin no. %d\n",i+1);
            snprintf(buffer, sizeof(buffer), "Scale nodes (ScaleNode[%d][0][0][])",i);
            fastNLOTools::PrintVector(GetScaleNodes(i,0),buffer,"#  ");
         }
      }
   }
   if ( abs(iprint) > 1 ) {
      cout << fastNLO::_SSEP20C << " Extended information (iprint > 1) " << fastNLO::_SSEP20 << endl;
      printf(" #   Printing of SigmaTilde not yet implemented.\n");
   }
   cout << fastNLO::_CSEPSC << endl;
}


//________________________________________________________________________________________________________________ //

// Erase observable bin
void fastNLOCoeffAddFix::EraseBin(unsigned int iObsIdx) {
   debug["fastNLOCoeffAddFix::EraseBin"]<<"Erasing table entries in CoeffAddFix for bin index " << iObsIdx << endl;
   if ( ScaleNode.size() == 0 ) {
      say::error["EraseBin"]<<"All fix-scale bins deleted already. Aborted!" << endl;
      exit(1);
   }
   if ( ScaleNode.size() != 0 ) ScaleNode.erase(ScaleNode.begin()+iObsIdx);
   if ( SigmaTilde.size() != 0 ) SigmaTilde.erase(SigmaTilde.begin()+iObsIdx);
   fastNLOCoeffAddBase::EraseBin(iObsIdx);
}

// Catenate observable bin
void fastNLOCoeffAddFix::CatBin(const fastNLOCoeffAddFix& other, unsigned int iObsIdx) {
   debug["fastNLOCoeffAddFix::CatBin"]<<"Catenating observable bin in CoeffAddFix corresponding to bin index " << iObsIdx << endl;
   if ( ScaleNode.size() == 0 ) {
      say::error["CatBin"]<<"Initial fix-scale table is empty. Aborted!" << endl;
      exit(1);
   }
   unsigned int nold = ScaleNode.size();
   if ( ScaleNode.size() != 0 ) {
      ScaleNode.resize(nold+1);
      ScaleNode[nold] = other.ScaleNode[iObsIdx];
   }
   if ( SigmaTilde.size() != 0 ) {
      SigmaTilde.resize(nold+1);
      SigmaTilde[nold] = other.SigmaTilde[iObsIdx];
   }
   fastNLOCoeffAddBase::CatBin(other, iObsIdx);
}

// Multiply observable bin
void fastNLOCoeffAddFix::MultiplyBin(unsigned int iObsIdx, double fact) {
   debug["fastNLOCoeffAddFix::MultiplyBin"]<<"Multiplying table entries in CoeffAddFix for bin index " << iObsIdx << " by factor " << fact << endl;
   for (unsigned int s=0 ; s<SigmaTilde[iObsIdx].size() ; s++) {
      for (unsigned int x=0 ; x<SigmaTilde[iObsIdx][s].size() ; x++) {
         for (unsigned int l=0 ; l<SigmaTilde[iObsIdx][s][x].size() ; l++) {
            for (unsigned int m=0 ; m<SigmaTilde[iObsIdx][s][x][m].size() ; m++) {
               SigmaTilde[iObsIdx][s][x][l][m] *= fact;
            }
         }
      }
   }
   fastNLOCoeffAddBase::MultiplyBin(iObsIdx, fact);
}
