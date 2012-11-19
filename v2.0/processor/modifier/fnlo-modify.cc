// author: Daniel Britzger, 19/11/12

// ------------------------------------------------- //
//   fnlo-modify
//   Program to modify fastNLO tables.
//
//
//   This prgram reads in the Steerfile:  SteerModify.str
//   If another Steerfile is wanted, pass it via the command line
//   using the tag 'steerfile'
//      fnlo-modify steerfile=AnotherFile.steer
//
//   See the header of SteerModify.str for more detail
// ------------------------------------------------- //

#include <cstdlib>
#include <vector>
#include <iostream>

#include "fnloTable.h"
#include "read_steer.h"

using namespace std;


template<typename T>
void EraseBin(vector<T>& v, unsigned int bin, unsigned int nobs){
   if ( v.empty() ){
      //cout<<"nothing todo."<<endl;
   }
   else if ( v.size() == nobs ) {
      //cout<<"erase."<<endl;
      v.erase(v.begin()+bin-1);
   } else {
      cout<<"Different number of bins."<<endl;
   }
}



int main(int argc, char** argv)
{
  
   if ( !PARSE(argc,argv) ) READ("SteerModify.str");
   if (BOOL(PrintSteeringCard) ) PRINTALL();
   
   cout<<"Reading table "<<STRING(InTable)<<endl;

   fnloTable table(CHAR(InTable));
   table.ReadTable();
   
   if ( BOOL(PrintInputA1) )  table.GetBlockA1()->Print();
   if ( BOOL(PrintInputA2) )  table.GetBlockA2()->Print();


   // Block A1
   if ( BOOL(Itabversion) ) {  
      cout<<"Modifying table version: from "<<table.GetBlockA1()->GetItabversion()<<" to " << INT(Itabversion)<<endl;
      table.GetBlockA1()->SetItabversion(INT(Itabversion));
   }
   if ( BOOL(ScenName) ) {  
      cout<<"Modifying table version: from "<<table.GetBlockA1()->GetScenName()<<" to " << STRING(ScenName)<<endl;
      table.GetBlockA1()->SetScenName(STRING(ScenName));
   }
   // Block A2
   if ( BOOL(Ipublunits) ) {  
      cout<<"Modifying table version: from '"<<table.GetBlockA2()->GetIpublunits()<<"' to '" << INT(Ipublunits)<<"'."<<endl;
      table.GetBlockA2()->SetIpublunits(INT(Ipublunits));
   }

   if ( !STRING_ARR(ScDescript).empty() ){
      cout<<"Modifying Scenario Description. ";
      if ( BOOL(AttachScDescription) ){
	 cout<<"Attaching "<<STRING_ARR(ScDescript).size()<<" lines to the scenario description"<< endl;
	 for ( unsigned int i = 0 ; i< STRING_ARR(ScDescript).size() ; i++ ){
	    table.GetBlockA2()->ScDescript.push_back(STRING_ARR(ScDescript)[i]);
	 }
	 table.GetBlockA2()->NScDescript = (int)table.GetBlockA2()->ScDescript.size();
      }
      else {
	 cout<<"Setting new scenario description with "<<STRING_ARR(ScDescript).size()<<" lines."<<endl;
	 table.GetBlockA2()->ScDescript  = STRING_ARR(ScDescript);
	 table.GetBlockA2()->NScDescript = (int)STRING_ARR(ScDescript).size();
      }
   }
   
   if ( BOOL(Ecms) ) {  
      cout<<"Modifying table version: from "<<table.GetBlockA2()->Ecms<<" to " << DOUBLE(Ecms)<<endl;
      table.GetBlockA2()->Ecms = DOUBLE(Ecms);
   }  

   if ( BOOL(ILOord) ) {  
      cout<<"Modifying table version: from "<<table.GetBlockA2()->ILOord<<" to " << INT(ILOord)<<endl;
      table.GetBlockA2()->ILOord = INT(ILOord);
   }  

   if ( BOOL(Ecms) ) {  
      cout<<"Modifying table version: from "<<table.GetBlockA2()->Ecms<<" to " << DOUBLE(Ecms)<<endl;
      table.GetBlockA2()->Ecms = DOUBLE(Ecms);
   }  

   if ( !INT_ARR(RemoveBins).empty() ) {  
      cout<<"Removing Bins!"<<endl;
      unsigned int nobs = table.GetNObsBin();
      for ( int i = (int)INT_ARR(RemoveBins).size()-1 ; i>=0 ; i-- ){
	 unsigned int b = INT_ARR(RemoveBins)[i];
	 cout<<"Erasing bin "<<i<<endl;
	 if ( b > table.GetNObsBin() ) {
	    cout<<"Error. Cannot erase bin number "<<b<<". There are only "<<table.GetNObsBin()<<" in the table."<<endl;
	    continue;
	 }
	 // A2
	 EraseBin(table.GetBlockA2()->LoBin,b,nobs);
	 EraseBin(table.GetBlockA2()->UpBin,b,nobs);
	 EraseBin(table.GetBlockA2()->BinSize,b,nobs);
	 EraseBin(table.GetBlockA2()->IDivLoPointer,b,nobs);
	 EraseBin(table.GetBlockA2()->IDivUpPointer,b,nobs);
 	 for ( unsigned int idx = 0; idx<table.GetBlockA1()->GetNcontrib() ; idx++ ){
	    EraseBin(table.GetBlockB(idx)->Xcenter,b,nobs);
	    EraseBin(table.GetBlockB(idx)->Value,b,nobs);
	    EraseBin(table.GetBlockB(idx)->fact,b,nobs);
	    EraseBin(table.GetBlockB(idx)->UncorLo,b,nobs);
	    EraseBin(table.GetBlockB(idx)->UncorHi,b,nobs);
	    EraseBin(table.GetBlockB(idx)->CorrLo,b,nobs);
	    EraseBin(table.GetBlockB(idx)->CorrHi,b,nobs);
	    EraseBin(table.GetBlockB(idx)->Nxtot1,b,nobs);
	    EraseBin(table.GetBlockB(idx)->XNode1,b,nobs);
	    EraseBin(table.GetBlockB(idx)->Nxtot2,b,nobs);
	    EraseBin(table.GetBlockB(idx)->XNode2,b,nobs);
	    EraseBin(table.GetBlockB(idx)->Nztot,b,nobs);
	    EraseBin(table.GetBlockB(idx)->ZNode,b,nobs);
	    EraseBin(table.GetBlockB(idx)->Xsection,b,nobs);

	    EraseBin(table.GetBlockB(idx)->ScaleNode,b,nobs);
	    EraseBin(table.GetBlockB(idx)->SigmaTilde,b,nobs);
	    EraseBin(table.GetBlockB(idx)->PdfLc,b,nobs);

	    EraseBin(table.GetBlockB(idx)->ScaleNode1,b,nobs);
	    EraseBin(table.GetBlockB(idx)->ScaleNode2,b,nobs);
	    // 	    EraseBin(table.GetBlockB(idx)->HScaleNode,b,nobs);
	    // 	    EraseBin(table.GetBlockB(idx)->HScaleNode1,b,nobs);
	    // 	    EraseBin(table.GetBlockB(idx)->HScaleNode2,b,nobs);

	    EraseBin(table.GetBlockB(idx)->SigmaTildeMuIndep,b,nobs);
	    EraseBin(table.GetBlockB(idx)->SigmaTildeMuFDep,b,nobs);
	    EraseBin(table.GetBlockB(idx)->SigmaTildeMuRDep,b,nobs);
	    EraseBin(table.GetBlockB(idx)->SigmaRefMixed,b,nobs);
	    EraseBin(table.GetBlockB(idx)->SigmaRef_s1,b,nobs);
	    EraseBin(table.GetBlockB(idx)->SigmaRef_s2,b,nobs);
 	 }		 
	 nobs--;
      }
      table.GetBlockA2()->NObsBin = nobs;
   }

   if ( BOOL(PrintOutputA1) )  table.GetBlockA1()->Print();
   if ( BOOL(PrintOutputA2) )  table.GetBlockA2()->Print();

   // writing modified table
   table.SetFilename(CHAR(OutTable));
   table.OpenFileWrite();  
   table.WriteBlockA1();
   table.WriteBlockA2();
   for ( unsigned int idx = 0; idx<table.GetBlockA1()->GetNcontrib() ; idx++ )  table.WriteBlockB(idx);
   table.CloseFileWrite();

   // goodbye
   cout<<"Wrote modified table to "<<STRING(OutTable)<<endl;
   
   return 0;
}
