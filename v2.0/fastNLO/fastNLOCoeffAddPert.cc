#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastNLOCoeffAddPert.h"

using namespace std;
using namespace fastNLO;

bool fastNLOCoeffAddPert::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet) {
   bool ret = fastNLOCoeffAddBase::CheckCoeffConstants(c,quiet);
   if ( ret && c->GetNScaleDep() == 0 ) return true;
   else if ( c->GetNScaleDep() >= 3 ) {
      if ( !quiet) 
	 say::error["fastNLOCoeffAddPert::CheckCoeffConstants"]
	    <<"This is not a fixed order v2.0  table. NScaleDep must be equal 0 but is NScaleDep="
	    <<c->GetNScaleDep()<<endl;
      return false;
   }
   else return false;
}

fastNLOCoeffAddPert::fastNLOCoeffAddPert(){
   SetClassName("fastNLOCoeffAddPert");
}

fastNLOCoeffAddPert::fastNLOCoeffAddPert(int NObsBin) : fastNLOCoeffAddBase(NObsBin) {
   SetClassName("fastNLOCoeffAddPert");
}

fastNLOCoeffAddPert::fastNLOCoeffAddPert(const fastNLOCoeffBase& base) : fastNLOCoeffAddBase(base) {
   SetClassName("fastNLOCoeffAddPert");
}


int fastNLOCoeffAddPert::Read(istream *table){
   fastNLOCoeffBase::ReadBase(table);
   ReadRest(table);
   return 0;
}

void fastNLOCoeffAddPert::ReadRest(istream *table){
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::ReadCoeffAddBase(table);
   ReadCoeffAddPert(table);
   EndReadCoeff(table);
}



int fastNLOCoeffAddPert::ReadCoeffAddPert(istream *table){
   CheckCoeffConstants(this);

   Nscalevar.resize(NScaleDim);
   Nscalenode.resize(NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      *table >> Nscalevar[i];
      *table >> Nscalenode[i];
   }
   // 	 printf("  *  fastNLOCoeffAddPert::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d,  NScaleDim %d  \n",
   // 	 fNObsBins, Nscalevar[0] , Nscalenode[0] , NScaleDim );
   ScaleFac.resize(NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      ScaleFac[i].resize(Nscalevar[i]);
      for(int j=0;j<Nscalevar[i];j++){
	 *table >> ScaleFac[i][j];
      }
   }
   //printf("  *  fastNLOCoeffAddPert::Read().bins %d, NScalevar[0] %d, Nscalenode[0] %d, ScaleFac[0][0] %d,  NScaleDim %d  \n",
   //fNObsBins, Nscalevar[0] , Nscalenode[0] , ScaleFac[0][0], NScaleDim );
   fastNLOCoeffBase::ResizeTable( &ScaleNode , fNObsBins, 1 , Nscalevar[0] , Nscalenode[0] ); // should work, since NScaleDim==1, but is not yet tested for 100%
   int nsn = ReadTable  ( &ScaleNode , table );
   //printf("  *  fastNLOCoeffAddPert::Read(). Read %d lines of ScaleNode.\n",nsn);
	 
   int XmaxFromI[1] = {0};
   //printf(" &SigmaTilde  %i  %i  %i  *%i  %i\n", fNObsBins, GetTotalScalevars(), GetTotalScalenodes(), XmaxFromI[0], NSubproc);
   fastNLOCoeffAddBase::ResizeTable( &SigmaTilde , fNObsBins, GetTotalScalevars(), GetTotalScalenodes(), XmaxFromI, NSubproc );
   int nst = ReadTable  ( &SigmaTilde , table );
   //printf("  *  fastNLOCoeffAddPert::Read(). Read %d lines of SigmaTilde.\n",nst);
   //printf("  *  fastNLOCoeffAddPert::Read(). Read %d lines of fastNLO v2.0 tables.\n",nst+nsn);
   info["Read"]<<"Read "<<nst+nsn<<" lines of fastNLO v2.0 tables."<<endl;

   // prepare members for evaluation
   fastNLOCoeffAddBase::ResizeTable(&PdfLc , fNObsBins, GetTotalScalenodes(), XmaxFromI, NSubproc);
   fastNLOCoeffBase::ResizeTable(&AlphasTwoPi_v20 , fNObsBins, GetTotalScalenodes());  
	 
   return 0;
}


int fastNLOCoeffAddPert::Write(ostream *table, int option){
   CheckCoeffConstants(this);
   fastNLOCoeffAddBase::Write(table,option);

   for(int i=0;i<NScaleDim;i++){
      *table << Nscalevar[i] << endl;
      *table << Nscalenode[i] << endl;
   }
   for(int i=0;i<NScaleDim;i++){
      for(int j=0;j<Nscalevar[i];j++){
	 *table << ScaleFac[i][j] << endl;
      }
   }
   int nsn = WriteTable( &ScaleNode  , table );
   //printf("  *  fastNLOCoeffAddPert::Write(). Wrote %d lines of ScaleNode.\n",nsn);
   int nst = WriteTable( &SigmaTilde , table , (bool)(option & DividebyNevt) , Nevt );
   //printf("  *  fastNLOCoeffAddPert::Write(). Wrote %d lines of SigmaTilde.\n",nst);
   printf("  *  fastNLOCoeffAddPert::Write(). Wrote %d lines of FASTNLO v2.0 tables.\n",nst+nsn);
   return 0;
}

int fastNLOCoeffAddPert::Copy(fastNLOCoeffAddPert* other){
   streambuf* streambuf = new stringbuf(ios_base::in | ios_base::out); 
   iostream* buffer = new iostream(streambuf);
   other->Write(buffer);
   *buffer << tablemagicno << endl;
   this->Read(buffer);
   delete buffer;
   delete streambuf;

   return(0);
}

void fastNLOCoeffAddPert::Add(fastNLOCoeffAddPert* other){
   double w1 = (double)Nevt / (Nevt+other->Nevt);
   double w2 = (double)other->Nevt / (Nevt+other->Nevt);
   Nevt += other->Nevt;
   CheckCoeffConstants(this);
   AddTableToAnotherTable( &SigmaTilde , &(other->SigmaTilde) ,w1 , w2 );
}



int fastNLOCoeffAddPert::GetTotalScalevars() const {
   int totalscalevars=1;
   for(int scaledim=0;scaledim<NScaleDim;scaledim++){
      totalscalevars *= Nscalevar[scaledim];
   }
   return totalscalevars;
}

int fastNLOCoeffAddPert::GetTotalScalenodes() const {
   int totalscalenodes=1;
   for(int scaledim=0;scaledim<NScaleDim;scaledim++){
      totalscalenodes *= Nscalenode[scaledim];
   }
   return totalscalenodes;
}


void fastNLOCoeffAddPert::Print() const {
   fastNLOCoeffAddBase::Print();
   printf(" **************** FastNLO Table: fastNLOCoeffAddPert ****************\n");
   for(int i=0;i<NScaleDim;i++){
      printf(" B    - Nscalenode[%d]              %d\n",i,Nscalenode[i]);
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
