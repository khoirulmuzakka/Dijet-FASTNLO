#include <cstdlib>
#include <unistd.h>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include "fastnlotk/fastNLOUserBlock.h"

using namespace std;
using namespace fastNLO;

//______________________________________________________________________________
fastNLOUserBlock::fastNLOUserBlock() : PrimalScream("fastNLOUserBlock") {

}


//______________________________________________________________________________
fastNLOUserBlock::~fastNLOUserBlock(){
}



//______________________________________________________________________________
void fastNLOUserBlock::Read(istream *table){
   //! Read user block from disk
   *table >> fUserFlag;
   int nLinesDescr;
   *table >> nLinesDescr;
   fDescr.resize(nLinesDescr);
   for ( int i = 0 ; i<nLinesDescr ; i++ )
      *table >> fDescr[i];
   *table >> fNLines;
   fContent.resize(fNLines);
   for ( int i = 0 ; i<fNLines ; i++ )
      *table >> fContent[i];
}



//______________________________________________________________________________
void fastNLOUserBlock::Write (ostream *table){
   //! Write UserBlock to disk
   if ( fNlines != 3+fDescr.size()+fContent.size() ) {
      error["Write"]<<"Number of lines in this user block is inconsistent with content."<<endl;
      error<<"\tfNlines="<<fNlines<<",\tCalculated number of lines: "<<3+fDescr.size()+fContent.size()<<endl;
      exit(1);
   }
   *table << fUserFlag << endl;
   *table << fDescr.size() <<endl;
   for ( unsigned int i = 0 ; i < fDescr.size() ; i++ )
      *table << fDescr[i] << endl;
   *table << fNLines << endl;
   *table << fContent.size() << endl;
   for ( unsigned int i = 0 ; i< fContent.size(); i++ )
      *table << fContent[i] << endl;
}



//______________________________________________________________________________
void fastNLOUserBlock::Print() const {
   printf("\n **************** FastNLO UserBlock ****************\n\n");
   printf("    UserFlag                       %d\n",fUserFlag);
   printf("    No. of lines                   %d\n",fNLines);
   printf("\n");
   printf("    Description:\n");
   for(unsigned int i=0;i<fDescr.size();i++)
      printf("    %s:\n",fDescr[i].c_str());
   printf("\n");
   printf("    Content:\n");
   for(unsigned int i=0;i<fContent.size();i++)
      printf("    %s:\n",fContent[i].c_str());
   printf("\n");
}

