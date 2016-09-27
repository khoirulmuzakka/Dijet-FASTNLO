///********************************************************************
///
///     fastNLO_toolkit: fnlo-tk-modify
///     Program to modify fastNLO v2 tables
///
///     The desired changes are set up in a steering file.
///     By default the steering file is named 'SteerModify.str'.
///     If another steering filename should be used, this filename can
///     be passed via the command line using the tag 'steerfile=',
///     e.g. fnlo-tk-modify steerfile=AnotherFileName.str
///
///     See the provided default steering file 'SteerModify.str' for more details,
///     or type:
///     ./fnlo-tk-modify -h
///
///     D. Britzger, K. Rabbertz
///
///********************************************************************

#include <cstdlib>
#include <iostream>
#include <vector>
#include "fastnlotk/fastNLOBase.h"
#include "fastnlotk/fastNLOTable.h"
#include "fastnlotk/read_steer.h"
#include "fastnlotk/speaker.h"

//__________________________________________________________________________________________________________________________________
int main(int argc, char** argv) {

   //! --- namespaces
   using namespace std;
   using namespace say;          //! namespace for 'speaker.h'-verbosity levels
   using namespace fastNLO;      //! namespace for fastNLO constants

   //! --- Set verbosity level
   SetGlobalVerbosity(INFO);

   //! --- Print program purpose
   cout << _CSEPSC << endl;
   info["fnlo-tk-modify"] << "Program to modify fastNLO v2 tables" << endl;
   cout << _SSEPSC << endl;
   info["fnlo-tk-modify"] << "The desired changes are set up in a steering file." << endl;
   info["fnlo-tk-modify"] << "By default the steering file is named 'SteerModify.str'." << endl;
   info["fnlo-tk-modify"] << "If another steering filename should be used, this filename can" << endl;
   info["fnlo-tk-modify"] << "be passed via the command line using the tag 'steerfile='," << endl;
   info["fnlo-tk-modify"] << "e.g. fnlo-tk-modify steerfile=AnotherFileName.str" << endl;
   info["fnlo-tk-modify"] << endl;
   info["fnlo-tk-modify"] << "See the provided default steering file 'SteerModify.str' for more details," << endl;
   info["fnlo-tk-modify"] << "or type:" << endl;
   info["fnlo-tk-modify"] << "./fnlo-tk-modify -h" << endl;
   cout << _CSEPSC << endl;

   //! ---  Parse commmand line
   cout << endl;
   cout << _CSEPSC << endl;
   shout["fnlo-tk-modify"] << "fastNLO Toolkit"<<endl;
   cout << _SSEPSC << endl;
   //! Test for default steering file "SteerModify.str"
   string steername = "SteerModify.str";
   if (argc <= 1) {
      ifstream ffile;
      ffile.open(steername.c_str());
      if (!ffile) {
         error["fnlo-tk-modify"] << "Neither mandatory parameters specified in command line " << endl;
         error["fnlo-tk-modify"] << "nor default steering file 'SteerModify.str' found. Aborted!" << endl;
         shout["fnlo-tk-modify"] << "For more explanations type:" << endl;
         shout["fnlo-tk-modify"] << "./fnlo-tk-modify -h" << endl;
         cout << _CSEPSC << endl;
         exit(1);
      }
   } else {
      steername = (const char*) argv[1];
      //! --- Usage info
      if (steername == "-h") {
         cout << " #" << endl;
         shout << "Usage: ./fnlo-tk-modify [steerfile=SteerFile.str] <InTable=fastNLOtableIn.tab> <OutTable=fastNLOtableOut.tab> [OptArg=option]" << endl;
         shout << "       Specification: <> mandatory; [] optional." << endl;
         shout << "       All desired table modifications like" << endl;
         shout << "       - changing the scenario name" << endl;
         shout << "       - changing or complementing the scenario description" << endl;
         shout << "       - adapting the cross section output units" << endl;
         shout << "       - correcting the power of the LO process" << endl;
         shout << "       - correcting the cms energy" << endl;
         shout << "       - cutting out unused bins or" << endl;
         shout << "       - multiplying bins by a set of factors" << endl;
         shout << "       are assumed to be controlled via steering parameters " << endl;
         shout << "       similar to the ones that may be used in table creation." << endl;
         shout << "       By default it is expected that at least the mandatory ones are" << endl;
         shout << "       specified in the steering file 'SteerModify.str' or" << endl;
         shout << "       are given via command line arguments:" << endl;
         shout << "[steerfile=SteerFile.str]:      Alternative steering filename" << endl;
         shout << "<InTable=fastNLOtableIn.tab>:   Table input filename, if not specified in steering file" << endl;
         shout << "<OutTable=fastNLOtableOut.tab>: Table output filename, if not specified in steering file" << endl;
         shout << "       For more steering options please check the default steering file delivered" << endl;
         shout << "       by the fastNLO Tolkit (usually in $prefix/share/fastnlo_toolkit/steerfiles)." << endl;
         cout << " #" << endl;
         cout  << _CSEPSC << endl;
         return 0;
      } else {
         shout["fnlo-tk-modify"] << "Parsing requested modifications ..." << endl;
      }
   }

   if ( !PARSE(argc,argv) ) {
      if ( ! (EXIST(InTable) && EXIST(OutTable)) ) {
         shout["fnlo-tk-modify"] << "Mandatory parameters not specified in command line," << endl;
         shout["fnlo-tk-modify"] << "trying to read from default steering file 'SteerModify.str'" << endl;
         int retcode = READ("SteerModify.str");
         cout << "retcode = " << retcode << endl;
         if ( retcode != 0 ) {
            error["fnlo-tk-modify"] << "Reading of mandatory parameters from default steering file 'SteerModify.str' unsuccessful. Aborted!" << endl;
            exit(retcode);
         }
         if ( ! (EXIST(InTable) && EXIST(OutTable)) ) {
            error["fnlo-tk-modify"] << "Mandatory parameters also not specified in default steering file 'SteerModify.str'. Aborted!" << endl;
            shout["fnlo-tk-modify"] << "For more explanations type:" << endl;
            shout["fnlo-tk-modify"] << "./fnlo-tk-modify -h" << endl;
            PRINTALL();
            exit(1);
         }
      }
   } else {
      shout["fnlo-tk-modify"] << "Parsing alternative steering file ..." << endl;
   }

   if ( ! (EXIST(InTable) && EXIST(OutTable)) ) {
      error["fnlo-tk-modify"] << "Input and/or output table not specified. Aborted!" << endl;
      PRINTALL();
      exit(1);
   }

   if ( EXIST(PrintSteeringCard) ) {
      if ( BOOL(PrintSteeringCard) ) PRINTALL();
   }

   shout["fnlo-tk-modify"] << "Trying to read input table " << STRING(InTable) << endl;

   fastNLOTable table(CHAR(InTable));
   table.ReadTable();

   if ( EXIST(PrintInputA1) ) {
      if ( BOOL(PrintInputA1) ) table.PrintHeader(1);
   }
   if ( EXIST(PrintInputA2) ) {
      if ( BOOL(PrintInputA2) ) table.PrintScenario(1);
   }

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

   if ( EXIST(BinSizeFactor) ) {
      double fac = DOUBLE(BinSizeFactor);
      info["fnlo-tk-modify"]<<"Multiplying all bin sizes by factor " << fac << "!" << endl;
      unsigned int nobs = table.GetNObsBin();
      for (unsigned int iObs=0; iObs<nobs; iObs++) {
         table.MultiplyBinSize(iObs,fac);
      }
   }

   if ( !DOUBLE_ARR(BinSize).empty() ) {
      vector<double> fac = DOUBLE_ARR(BinSize);
      info["fnlo-tk-modify"]<<"Multiplying bin sizes by provided factors!"<<endl;
      unsigned int nobs = table.GetNObsBin();
      if ( nobs != fac.size() ) {
         error["fnlo-tk-modify"]<<"You need the same number of multiplicative factors, nfact = " << fac.size() <<
            ", than bins in the table, nobsbin = " << nobs << ". Aborted!" << endl;
         exit(1);
      }
      for (unsigned int iObs=0; iObs<nobs; iObs++) {
         table.MultiplyBinSize(iObs,fac[iObs]);
      }
   }

   if ( !DOUBLE_ARR(MultCoeff).empty() ) {
      vector<double> fac = DOUBLE_ARR(MultCoeff);
      info["fnlo-tk-modify"]<<"Multiplying by provided factors all coefficients of additive contributions to observable bins!"<<endl;
      unsigned int nobs = table.GetNObsBin();
      if ( nobs != fac.size() ) {
         error["fnlo-tk-modify"]<<"You need the same number of multiplicative factors, nfact = " << fac.size() <<
            ", than bins in the table, nobsbin = " << nobs << ". Aborted!" << endl;
         exit(1);
      }
      for (unsigned int iObs=0; iObs<nobs; iObs++) {
         table.MultiplyBinInTable(iObs,fac[iObs]);
      }
   }

   if ( !INT_ARR(RemoveBins).empty() ) {
      info["fnlo-tk-modify"]<<"Removing observable bins from interpolation table!"<<endl;
      unsigned int nobs = table.GetNObsBin();
      info["fnlo-tk-modify"]<<"Existing observable bins: " << nobs <<endl;
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
   }

   if ( EXIST(PrintOutputA1) ) {
      if ( BOOL(PrintOutputA1) ) table.PrintHeader(0);
   }

   if ( EXIST(PrintOutputA2) ) {
      if ( BOOL(PrintOutputA2) ) table.PrintScenario(0);
   }

   // writing modified table
   table.SetFilename(CHAR(OutTable));
   table.WriteTable();

   // Finished
   info["fnlo-tk-modify"]<<"Wrote modified table to " << STRING(OutTable) << endl;

   return 0;
}
