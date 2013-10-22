#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastnlotk/fastNLOCoeffAddBase.h"

using namespace std;
using namespace fastNLO;

bool fastNLOCoeffAddBase::CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet) {
   //   if(!(!(IDataFlag==1) && !(IAddMultFlag==1))){
   if ( c->GetIDataFlag()==0 && c->GetIAddMultFlag()==0 ) return true;
   else if ( c->GetIDataFlag()==1 || c->GetIAddMultFlag()==1 ) {
      if ( !quiet)say::error["fastNLOCoeffAddBase::CheckCoeffConstants"]
	 <<"This must be a table with multiplicative perturbative coefficients. IDataFlag="
	 <<c->GetIDataFlag()<<", IAddMultFlag="<<c->GetIAddMultFlag()
	 <<", but none is allowed to be 1."<<endl;
      return false;
   }
   else return false;
}

fastNLOCoeffAddBase::fastNLOCoeffAddBase(){
}

fastNLOCoeffAddBase::fastNLOCoeffAddBase(int NObsBin) : fastNLOCoeffBase(NObsBin) {
   NScaleDep = 0;
}

fastNLOCoeffAddBase::fastNLOCoeffAddBase(const fastNLOCoeffBase& base ) : fastNLOCoeffBase(base) {
}



int fastNLOCoeffAddBase::Read(istream *table){
   fastNLOCoeffBase::ReadBase(table);
   CheckCoeffConstants(this);
   ReadCoeffAddBase(table);
   EndReadCoeff(table);
   return 0;
}


int fastNLOCoeffAddBase::ReadCoeffAddBase(istream *table){
   CheckCoeffConstants(this);
   char buffer[5257];
   *table >> IRef;
   *table >> IScaleDep;
   *table >> Nevt;
   *table >> Npow;
   int NPDF;
   *table >> NPDF;
   if(NPDF>0){
      NPDFPDG.resize(NPDF);
      for(int i=0;i<NPDF;i++){
	 *table >>  NPDFPDG[i];
      }
   }
   *table >> NPDFDim;
   int NFragFunc;
   *table >> NFragFunc;
   if(NFragFunc>0){
      NFFPDG.resize(NFragFunc);
      for(int i=0;i<NFragFunc;i++){
	 *table >>  NFFPDG[i];
      }
   }
   *table >> NFFDim;
   *table >> NSubproc;
   *table >> IPDFdef1;
   *table >> IPDFdef2;
   *table >> IPDFdef3;
   //printf("  *  fastNLOCoeffAddBase::Read(). IRef : %d, IScaleDep: %d, Nevt: %d, Npow: %d, NPDF: %d, NPDFDim: %d\n", IRef ,IScaleDep  ,Nevt  , Npow ,NPDF , NPDFDim  );

   if(IPDFdef1==0){
      for(int i=0;i<NSubproc;i++){
	 // Missing: linear PDF combinations for IPDFdef1=0
	 if(NPDF==1){
	 }else{
	    if(NPDF==2){
	    }
	 }
      }
   }
   //Nxtot1.resize(fNObsBins);
   XNode1.resize(fNObsBins);
   for(int i=0;i<fNObsBins;i++){
      int xtot;
      *table >> xtot;
      //*table >> Nxtot1[i];
      //XNode1[i].resize(Nxtot1[i]);
      XNode1[i].resize(xtot);
      for(int j=0;j<xtot;j++){
	 *table >> XNode1[i][j];
      }
   }
   if(NPDFDim==2){
      //Nxtot2.resize(fNObsBins);
      XNode2.resize(fNObsBins);
      for(int i=0;i<fNObsBins;i++){
	 int xtot;
	 *table >> xtot;
	 XNode2[i].resize(xtot);
	 //*table >> Nxtot2[i];
	 //XNode2[i].resize(Nxtot2[i]);
	 for(int j=0;j<xtot;j++){
	    *table >> XNode2[i][j];
	 }
      }
   }
   if(NFragFunc>0){
      Nztot.resize(fNObsBins);
      ZNode.resize(fNObsBins);
      for(int i=0;i<fNObsBins;i++){
	 *table >> Nztot[i];
	 ZNode[i].resize(Nztot[i]);
	 for(int j=0;j<Nztot[i];j++){
	    *table >> ZNode[i][j];
	 }
      }
   }

   int NScales;
   *table >> NScales;
   *table >> NScaleDim;
   Iscale.resize(NScales);
   for(int i=0;i<NScales;i++){
      *table >> Iscale[i];
   }
   int NscaleDescript;
   ScaleDescript.resize(NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      *table >> NscaleDescript;
      ScaleDescript[i].resize(NscaleDescript);
      table->getline(buffer,256);
      for(int j=0;j<NscaleDescript;j++){
	 table->getline(buffer,256);
	 ScaleDescript[i][j] = buffer;
	 //            StripWhitespace(ScaleDescript[i][j]);
      }
   }
   return 0;
}

int fastNLOCoeffAddBase::Write(ostream *table, int option){
   cout<<"fastNLOCoeffAddBase::Write(), calling Base::Write()."<<endl;
   fastNLOCoeffBase::Write(table,option);
   CheckCoeffConstants(this);

   *table << IRef << endl;
   *table << IScaleDep << endl;
   *table << Nevt << endl;
   *table << Npow << endl;
   *table << NPDFPDG.size() << endl;
   for(unsigned int i=0;i<NPDFPDG.size();i++){
      *table <<  NPDFPDG[i] << endl;
   }
   *table << NPDFDim << endl;
   int NFragFunc = NFFPDG.size();
   *table << NFragFunc << endl;
   if(NFragFunc>0){
      for(int i=0;i<NFragFunc;i++){
	 *table <<  NFFPDG[i] << endl;
      }
   }
   *table << NFFDim << endl;
   *table << NSubproc << endl;
   *table << IPDFdef1 << endl;
   *table << IPDFdef2 << endl;
   *table << IPDFdef3 << endl;
   if(IPDFdef1==0){
      for(int i=0;i<NSubproc;i++){
	 // Missing: linear PDF combinations for IPDFdef1=0
	 if(NPDFPDG.size()==1){
	 }else{
	    if(NPDFPDG.size()==2){
	    }
	 }
      }
   }
   for(int i=0;i<fNObsBins;i++){
      *table << XNode1[i].size() << endl;
      for(unsigned int j=0;j<XNode1[i].size();j++){
	 *table << XNode1[i][j] << endl;
      }
   }
   if(NPDFDim==2){
      for(int i=0;i<fNObsBins;i++){
	 *table << XNode2[i].size() << endl;
	 for(unsigned int j=0;j<XNode2[i].size();j++){
	    *table << XNode2[i][j] << endl;
	 }
      }
   }
   cout<<" 10"<<endl;
   if(NFragFunc>0){
      for(int i=0;i<fNObsBins;i++){
	 *table << Nztot[i] << endl;
	 for(int j=0;j<Nztot[i];j++){
	    *table << ZNode[i][j] << endl;
	 }
      }
   }
   int NScales = Iscale.size();
   *table << NScales << endl;
   *table << NScaleDim << endl;
   for(int i=0;i<NScales;i++){
      *table << Iscale[i] << endl;
   }
   for(int i=0;i<NScaleDim;i++){
      *table << ScaleDescript[i].size() << endl;
      for(unsigned int j=0;j<ScaleDescript[i].size();j++){
	 *table << ScaleDescript[i][j] << endl;
      }
   }

   //! v2.1 store NScaleDep here
   //! *table << NScaleDep << endl;
   cout<<"fastNLOCoeffAddBase. Writing coefficients."<<endl;
   return 0;
}

int fastNLOCoeffAddBase::Copy(fastNLOCoeffAddBase* other){
   streambuf* streambuf = new stringbuf(ios_base::in | ios_base::out);
   iostream* buffer = new iostream(streambuf);
   other->Write(buffer);
   *buffer << tablemagicno << endl;
   this->Read(buffer);
   delete buffer;
   delete streambuf;
   return(0);
}

void fastNLOCoeffAddBase::Add(fastNLOCoeffAddBase* other){
   //    double w1 = (double)Nevt / (Nevt+other->Nevt);
   //    double w2 = (double)other->Nevt / (Nevt+other->Nevt);
   Nevt += other->Nevt;
}


int fastNLOCoeffAddBase::GetNxmax(int i) const {
   int nxmax = 0;
   switch (NPDFDim) {
   case 0: nxmax = (int)XNode1[i].size();
      break;
      //   case 1: nxmax = ((int)pow((double)Nxtot1[i],2)+Nxtot1[i])/2;
   case 1: nxmax = ((int)pow((double)XNode1[i].size(),2)+XNode1[i].size())/2;
      break;
   case 2: nxmax = XNode1[i].size()*XNode2[i].size();
      break;
   default: ;
   }
   return nxmax;
};

int fastNLOCoeffAddBase::GetXIndex(int Obsbin,int x1bin,int x2bin) const {
   int xbin = 0;
   switch (NPDFDim) {
   case 0: xbin = x1bin; // linear
      break;
   case 1: xbin = x1bin + (x2bin*(x2bin+1)/2);    // half matrix
      break;
   case 2: xbin = x1bin + x2bin * XNode1[Obsbin].size(); // full matrix
      break;
   default: ;
   }
   return xbin;
};



//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::ResizeTable(vector<vector<vector<vector<vector<vector<vector<double > > > > > > >* v, int dim0 , int dim1, int dim2, int dim3, int dim4, int* dim5GetNxmaxFromDimI , int dim6 ){
  if ( dim0 > 0 ){
    if ( dim5GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
	fastNLOCoeffBase::ResizeTable( &(v->at(i)) , dim1, dim2, dim3, dim4, GetNxmax(i), dim6 );
      }
    }
    else if ( dim5GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }

}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int dim1, int dim2, int* dim3GetNxmaxFromDimI, int dim4 ){
  if ( dim0 > 0 ){
    if ( dim3GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
	 fastNLOCoeffBase::ResizeTable( &(v->at(i)) , dim1, dim2, GetNxmax(i), dim4 );
      }
    }
    else if ( dim3GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::ResizeTable( vector<vector<vector<vector<vector<double > > > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2, int dim3, int dim4 ){
  if ( dim0 > 0 ){
    if ( dim1GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
	fastNLOCoeffBase::ResizeTable( &(v->at(i)) , GetNxmax(i), dim2, dim3, dim4 );
      }
    }
    else if ( dim1GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }

}

//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::ResizeTable( vector<vector<vector<vector<double > > > >* v, int dim0 , int dim1, int* dim2GetNxmaxFromDimI, int dim3 ){
  if ( dim0 > 0 ){
    if ( dim2GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
	fastNLOCoeffBase::ResizeTable( &(v->at(i)) , dim1, GetNxmax(i), dim3 );
      }
    }
    else if ( dim2GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }

}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::ResizeTable( vector<vector<vector<double > > >* v, int dim0 , int* dim1GetNxmaxFromDimI, int dim2 ){
  if ( dim0 > 0 ){
    if ( dim1GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
	fastNLOCoeffBase::ResizeTable( &(v->at(i)) , GetNxmax(i), dim2 );
      }
    }
    else if ( dim1GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }

}


//________________________________________________________________________________________________________________ //
void fastNLOCoeffAddBase::ResizeTable( vector<vector<double > >*  v, int dim0 , int* dim1GetNxmaxFromDimI ){
  if ( dim0 > 0 ){
    if ( dim1GetNxmaxFromDimI[0] == 0 ) {
      v->resize(dim0);
      for ( int i= 0 ; i<dim0 ; i++){
	fastNLOCoeffBase::ResizeTable( &(v->at(i)) , GetNxmax(i) );
      }
    }
    else if ( dim1GetNxmaxFromDimI[0] != 0 ){
      cout << "Error in Resize Table. This is not yet implemented" << endl;
      exit(1);
    }
  } else{
    cout << "Error in Resize Table." << endl;
    exit(1);
  }
}

//________________________________________________________________________________________________________________ //


void fastNLOCoeffAddBase::Print() const {
   fastNLOCoeffBase::Print();
   printf(" **************** FastNLO Table: fastNLOCoeffAddBase ****************\n");
   printf(" B   IRef                          %d\n",IRef);
   printf(" B   IScaleDep                     %d\n",IScaleDep);
   printf(" B   Nevt                          %llu\n",Nevt);
   printf(" B   Npow                          %d\n",Npow);
   printf(" B   NPDF                          %d\n",NPDFPDG.size());
   for(unsigned int i=0;i<NPDFPDG.size();i++){
      printf(" B    - NPDFPDG[%d]                 %d\n",i,NPDFPDG[i]);
   }
   printf(" B   NPDFDim                       %d\n",NPDFDim);
   printf(" B   NFragFunc                     %d\n",NFFPDG.size());
   for(unsigned int i=0;i<NFFPDG.size();i++){
      printf(" B    - NFFPDG[%d]               %d\n",i,NFFPDG[i]);
   }
   printf(" B   NFFDim                        %d\n",NFFDim);
   printf(" B   NSubproc                      %d\n",NSubproc);
   printf(" B   IPDFdef1                      %d\n",IPDFdef1);
   printf(" B   IPDFdef2                      %d\n",IPDFdef2);
   printf(" B   IPDFdef3                      %d\n",IPDFdef3);
   printf(" B   Nxtot1[0-%d]             ",fNObsBins);
   for(int i=0;i<fNObsBins;i++){
      printf("%d ,",XNode1[i].size());
   }
   printf("\n");
   //     for(int i=0;i<fNObsBins;i++){
   //       printf(" B    XNode1[%d]             ",i);
   //       for(int j=0;j<Nxtot1[i];j++){
   // 	printf(" B   %8.4f ,",XNode1[i][j]);
   //       }
   //       printf(" B   \n");
   //     }
   printf(" B   if (NPDFDim==2), you could print xnodes2 here. (NPDFDim = %d)\n",NPDFDim);
   printf(" B   if (NFragFunc>0), you could print xnodes2 here. (NFragFunc = %d)\n",NFFPDG.size());
   printf(" B   NScales                       %d\n",Iscale.size());
   for(int i=0;i<Iscale.size();i++){
      printf(" B    - Iscale[%d]                  %d\n",i,Iscale[i]);
   }
   printf(" B   NScaleDim                     %d\n",NScaleDim);
   for(int i=0;i<NScaleDim;i++){
      //printf(" B    -  NscaleDescript[%d]         %d\n",i,NscaleDescript[i]);
      for(unsigned int j=0;j<ScaleDescript[i].size();j++){
	 printf(" B    -  - ScaleDescript[%d][%d]     %s\n",i,j,ScaleDescript[i][j].data());
      }
   }
   printf(" *******************************************************\n");
}


//________________________________________________________________________________________________________________ //
