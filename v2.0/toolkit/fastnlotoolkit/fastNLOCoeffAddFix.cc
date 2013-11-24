#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffAddFix.h"

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
fastNLOCoeffAddFix::fastNLOCoeffAddFix(){
   SetClassName("fastNLOCoeffAddFix");
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
int fastNLOCoeffAddFix::Read(istream *table){
   fastNLOCoeffBase::ReadBase(table);
   ReadRest(table);
   return 0;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::ReadRest(istream *table){
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::ReadCoeffAddBase(table);
   ReadCoeffAddFix(table);
   EndReadCoeff(table);
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddFix::ReadCoeffAddFix(istream *table){
   CheckCoeffConstants(this);

   Nscalevar.resize(NScaleDim);
   vector<int> Nscalenode(NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      *table >> Nscalevar[i];
      *table >> Nscalenode[i];
   }
   // 	 printf("  *  fastNLOCoeffAddFix::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d,  NScaleDim %d  \n",
   // 	 fNObsBins, Nscalevar[0] , Nscalenode[0] , NScaleDim );
   ScaleFac.resize(NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      ScaleFac[i].resize(Nscalevar[i]);
      for(int j=0;j<Nscalevar[i];j++){
	 *table >> ScaleFac[i][j];
      }
   }
   //printf("  *  fastNLOCoeffAddFix::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d, ScaleFac[0][0] %d,  NScaleDim %d  \n",
   //fNObsBins, Nscalevar[0] , Nscalenode[0] , ScaleFac[0][0], NScaleDim );
   fastNLOCoeffBase::ResizeTable( &ScaleNode , fNObsBins, 1 , Nscalevar[0] , Nscalenode[0] ); // should work, since NScaleDim==1, but is not yet tested for 100%
   int nsn = ReadTable  ( &ScaleNode , table );
   //printf("  *  fastNLOCoeffAddFix::Read(). Read %d lines of ScaleNode.\n",nsn);

   int XmaxFromI[1] = {0};
   //printf(" &SigmaTilde  %i  %i  %i  *%i  %i\n", fNObsBins, GetTotalScalevars(), GetTotalScalenodes(), XmaxFromI[0], NSubproc);
   fastNLOCoeffAddBase::ResizeTable( &SigmaTilde , fNObsBins, GetTotalScalevars(), GetTotalScalenodes(), XmaxFromI, NSubproc );
   int nst = ReadTable  ( &SigmaTilde , table );
   //printf("  *  fastNLOCoeffAddFix::Read(). Read %d lines of SigmaTilde.\n",nst);
   //printf("  *  fastNLOCoeffAddFix::Read(). Read %d lines of fastNLO v2.0 tables.\n",nst+nsn);
   info["Read"]<<"Read "<<nst+nsn<<" lines of fastNLO v2.0 tables."<<endl;

   // prepare members for evaluation
   fastNLOCoeffAddBase::ResizeTable(&PdfLc , fNObsBins, GetTotalScalenodes(), XmaxFromI, NSubproc);
   fastNLOCoeffBase::ResizeTable(&AlphasTwoPi_v20 , fNObsBins, GetTotalScalenodes());

   return 0;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Write(ostream* table,double Nevt){
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::Write(table);

   for(int i=0;i<NScaleDim;i++){
      *table << Nscalevar[i] << endl;
      *table << GetNScaleNode() << endl;
   }
   for(int i=0;i<NScaleDim;i++){
      for(int j=0;j<Nscalevar[i];j++){
	 *table << ScaleFac[i][j] << endl;
      }
   }

   int nsn = WriteTable( ScaleNode  , table , Nevt);
   //printf("  *  fastNLOCoeffAddFix::Write(). Wrote %d lines of ScaleNode.\n",nsn);
   //int nst = WriteTable( &SigmaTilde , table , Nevt );
   int nst = WriteTable( SigmaTilde , table , Nevt);
   //printf("  *  fastNLOCoeffAddFix::Write(). Wrote %d lines of SigmaTilde.\n",nst);
   printf("  *  fastNLOCoeffAddFix::Write(). Wrote %d lines of FASTNLO v2.0 tables.\n",nst+nsn);
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Add(const fastNLOCoeffAddFix& other){
   CheckCoeffConstants(this);
   double w1 = (double)Nevt / (Nevt+other.Nevt);
   double w2 = (double)other.Nevt / (Nevt+other.Nevt);
   Nevt += other.Nevt;
   AddTableToAnotherTable( SigmaTilde , other.SigmaTilde ,w1 , w2 );
}



//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddFix::GetTotalScalevars() const {
   int totalscalevars=1;
   for(int scaledim=0;scaledim<NScaleDim;scaledim++){
      totalscalevars *= Nscalevar[scaledim];
   }
   return totalscalevars;
}


//________________________________________________________________________________________________________________ //
int fastNLOCoeffAddFix::GetTotalScalenodes() const {
   if ( !ScaleNode.empty() ) return ScaleNode[0][0][0].size();
   else return 0;
   //    int totalscalenodes=1;
   //    for(int scaledim=0;scaledim<NScaleDim;scaledim++){
   //       if ( !ScaleNode.empty() )
   // 	 totalscalenodes *= ScaleNode[0][0][0].size();
   //       //Nscalenode[scaledim];
   //    }
   //    return totalscalenodes;
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddFix::Print() const {
   fastNLOCoeffAddBase::Print();
   printf(" **************** FastNLO Table: fastNLOCoeffAddFix ****************\n");
   for(int i=0;i<NScaleDim;i++){
      printf(" B    - Nscalenode[%d]              %d\n",i,GetNScaleNode());
      printf(" B    - Nscalevar[%d]               %d\n",i,Nscalevar[i]);
      for(int j=0;j<Nscalevar[i];j++){
	 printf(" B    -  - ScaleFac[%d][%d]          %6.4f\n",i,j,ScaleFac[i][j]);
      }
   }
   printf(" B   No printing of ScaleNode implemented yet.\n");
   printf(" B   No printing of SigmaTilde implemented yet.\n");
   printf(" *******************************************************\n");
}


//________________________________________________________________________________________________________________ //
