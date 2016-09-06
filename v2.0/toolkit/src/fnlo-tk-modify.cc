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

#include "fastnlotk/fastNLOBase.h"
#include "fastnlotk/fastNLOTable.h"
#include "fastnlotk/read_steer.h"
#include "fastnlotk/speaker.h"

//! --- namespaces
using namespace std;
using namespace say;       //! namespace for 'speaker.h'-verbosity levels




int main(int argc, char** argv)
{

   if ( !PARSE(argc,argv) ) READ("SteerModify.str");
   if (BOOL(PrintSteeringCard) ) PRINTALL();

   cout << "Reading table " << STRING(InTable) << endl;

   fastNLOTable table(CHAR(InTable));
   table.ReadTable();

   if ( BOOL(PrintInputA1) ) table.PrintHeader(1);
   if ( BOOL(PrintInputA2) ) table.PrintScenario(1);

   // Block A1
   if ( EXIST(Itabversion) ) {
      info["fnlo-tk-modify"]<<"Modifying table version: from "<< table.GetItabversion() << " to " << INT(Itabversion) << endl;
      table.SetItabversion(INT(Itabversion));
   }
   if ( EXIST(ScenName) ) {
      info["fnlo-tk-modify"]<<"Modifying scenario name: from "<< table.GetScenName() << " to " << STRING(ScenName) << endl;
      table.SetScenName(STRING(ScenName));
   }
   // Block A2
   if ( EXIST(Ipublunits) ) {
      info["fnlo-tk-modify"]<<"Modifying publication units: from "<< table.GetIpublunits() << " to " << INT(Ipublunits) << endl;
      table.SetIpublunits(INT(Ipublunits));
   }

   if ( !STRING_ARR(ScDescript).empty() ){
      vector <string> ScDescr = table.GetScDescr();
      size_t NScSize = ScDescr.size();
      info["fnlo-tk-modify"]<<"Modifying existing scenario description:" << endl;
      for ( size_t i = 0; i < NScSize; i++ ) {
         cout << "Line no. " << i << ": " << ScDescr[i] << endl;
      }
      if ( BOOL(AttachScDescription) ){
         info["fnlo-tk-modify"]<<"Attaching lines:" << endl;
         size_t NewNScSize = NScSize + STRING_ARR(ScDescript).size();
         ScDescr.resize(NewNScSize);
         for ( size_t i = NScSize; i < NewNScSize; i++ ) {
            ScDescr[i] = STRING_ARR(ScDescript)[i-NScSize];
            cout << "Line no. " << i << ": " << ScDescr[i] << endl;
         }
      } else {
         info["fnlo-tk-modify"]<<"Replacing lines with:" << endl;
         size_t NewNScSize = STRING_ARR(ScDescript).size();
         ScDescr.resize(NewNScSize);
         for ( size_t i = 0; i < NewNScSize; i++ ) {
            ScDescr[i] = STRING_ARR(ScDescript)[i];
            cout << "Line no. " << i << ": " << ScDescr[i] << endl;
         }
      }
      table.SetScDescr(ScDescr);
   }

   if ( EXIST(Ecms) ) {
      info["fnlo-tk-modify"]<<"Modifying center-of-mass energy: from "<< table.GetEcms() << " to " << DOUBLE(Ecms) << endl;
      table.SetEcms(DOUBLE(Ecms));
   }

   if ( EXIST(ILOord) ) {
      info["fnlo-tk-modify"]<<"Modifying LO: from "<< table.GetLoOrder() << " to " << INT(ILOord) << endl;
      table.SetLoOrder(INT(ILOord));
   }

   if ( !INT_ARR(RemoveBins).empty() ) {
      info["fnlo-tk-modify"]<<"Removing observable bins from interpolation table!"<<endl;
      unsigned int nobs = table.GetNObsBin();
      info["fnlo-tk-modify"]<<"Existing observable bins: " << nobs <<endl;
      //      vector < double > NewBinSize = table.GetBinSize();
      for ( int i = (int)INT_ARR(RemoveBins).size(); i>0; i-- ) {
         unsigned int iObs = INT_ARR(RemoveBins)[i-1];
         info["fnlo-tk-modify"]<<"Bins left to remove: " << i << endl;
         info["fnlo-tk-modify"]<<"Erasing bin no. " << iObs << endl;
         if ( iObs > table.GetNObsBin() ) {
            warn["fnlo-tk-modify"]<<"Cannot erase bin no. " << iObs << ". There are only " << table.GetNObsBin() << " bins in the table. Ignored!"<<endl;
            continue;
         }
         table.EraseBinFromTable(iObs-1);
         nobs--;
      }
      //      table.SetBinSize(NewBinSize);
      //table.SetNObsBin(nobs);
   }

   // if ( !DOUBLE_ARR(MultFac).empty() ){
   //    vector<double> fac = DOUBLE_ARR(MultFac);
   //    cout<<"Applying multiplicative Factors!"<<endl;
   //    unsigned int nobs = table.GetNObsBin();
   //    if ( nobs!=fac.size() ){
   //       cout<<"ERROR. You need same number of multiplicative factors, than bins in the table."<<endl;
   //       exit(1);
   //    }
   //    // SigmaTilde [NObsBins] ['n' x-nodes] [n s1-Nodes] [n s2-Nodes] [nsubproc]
   //    //    vector < vector < vector < vector < vector < double > > > > > SigmaTildeMuIndep;
   //    //    vector < vector < vector < vector < vector < double > > > > > SigmaTildeMuFDep;
   //    //    vector < vector < vector < vector < vector < double > > > > > SigmaTildeMuRDep;
   //    // SigmaRef [NObsBins] [nsubproc]
   //    //    vector < vector < double > > SigmaRefMixed;
   //    //    vector < vector < double > > SigmaRef_s1;
   //    //    vector < vector < double > > SigmaRef_s2;
   //    //    vector < vector < vector < vector < vector < double > > > > > SigmaTilde;
   //    for ( unsigned int idx = 0; idx<table.GetBlockA1()->GetNcontrib() ; idx++ ){
   //       fnloBlockB* b = table.GetBlockB(idx);
   //       for ( int i = 0 ; i<nobs;i++ ){
   //          if ( !b->SigmaTildeMuIndep.empty() ) {
   //             for ( unsigned int j = 0 ; j < b->SigmaTildeMuIndep[i].size() ; j++ ){
   //                for ( unsigned int k = 0 ; k < b->SigmaTildeMuIndep[i][j].size() ; k++ ){
   //                   for ( unsigned int l = 0 ; l < b->SigmaTildeMuIndep[i][j][k].size() ; l++ ){
   //                      for ( unsigned int m = 0 ; m < b->SigmaTildeMuIndep[i][j][k][l].size() ; m++ ){
   //                         b->SigmaTildeMuIndep[i][j][k][l][m] *= fac[i];
   //                         if ( !b->SigmaTildeMuFDep.empty() )
   //                            b->SigmaTildeMuFDep[i][j][k][l][m] *= fac[i];
   //                         if ( !b->SigmaTildeMuRDep.empty() )
   //                            b->SigmaTildeMuRDep[i][j][k][l][m] *= fac[i];
   //                      }
   //                   }
   //                }
   //             }
   //          }
   //          if ( !b->SigmaTilde.empty() ) {
   //             for ( unsigned int j = 0 ; j < b->SigmaTilde[i].size() ; j++ ){
   //                for ( unsigned int k = 0 ; k < b->SigmaTilde[i][j].size() ; k++ ){
   //                   for ( unsigned int l = 0 ; l < b->SigmaTilde[i][j][k].size() ; l++ ){
   //                      for ( unsigned int m = 0 ; m < b->SigmaTilde[i][j][k][l].size() ; m++ ){
   //                         b->SigmaTilde[i][j][k][l][m] *= fac[i];
   //                      }
   //                   }
   //                }
   //             }
   //          }
   //          if ( !b->SigmaRefMixed.empty() ) {
   //             for ( unsigned int j = 0 ; j < b->SigmaRefMixed[i].size() ; j++ ){
   //                   b->SigmaRefMixed[i][j] *= fac[i];
   //             }
   //          }
   //          if ( !b->SigmaRef_s1.empty() ) {
   //             for ( unsigned int j = 0 ; j < b->SigmaRef_s1[i].size() ; j++ ){
   //                b->SigmaRef_s1[i][j] *= fac[i];
   //             }
   //          }
   //          if ( !b->SigmaRef_s2.empty() ) {
   //             for ( unsigned int j = 0 ; j < b->SigmaRef_s2[i].size() ; j++ ){
   //                b->SigmaRef_s2[i][j] *= fac[i];
   //             }
   //          }
   //       }
   //    }
   // }

   if ( BOOL(PrintOutputA1) ) table.PrintHeader(0);
   if ( BOOL(PrintOutputA2) ) table.PrintScenario(0);

   // writing modified table
   table.SetFilename(CHAR(OutTable));
   table.WriteTable();

   // Finished
   info["fnlo-tk-modify"]<<"Wrote modified table to " << STRING(OutTable) << endl;

   return 0;
}
