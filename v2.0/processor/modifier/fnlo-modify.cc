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
