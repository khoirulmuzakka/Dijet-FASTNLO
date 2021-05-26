///********************************************************************
///
///     fastNLO_toolkit: fnlo-tk-modify
///     Tool to manipulate a fastNLO table
///
///     For more explanations type:
///     ./fnlo-tk-modify -h
///     and consult the provided default steering file 'SteerModify.str'.
///
///     D. Britzger, K. Rabbertz
///
///********************************************************************

#include <cstdlib>
#include <iostream>
#include <vector>
#include "fastnlotk/fastNLOTable.h"
#include "fastnlotk/fastNLOTools.h"
#include "fastnlotk/read_steer.h"
#include "fastnlotk/speaker.h"

//__________________________________________________________________________________________________________________________________
int main(int argc, char** argv) {

   //! --- namespaces
   using namespace std;
   using namespace say;          //! namespace for 'speaker.h'-verbosity levels
   using namespace fastNLO;      //! namespace for fastNLO constants

   //! --- Set initial verbosity level
   SetGlobalVerbosity(INFO);

   //! --- Parse command line
   string steername = "SteerModify.str"; //! Default steering file "SteerModify.str"
   string intable;
   string outtable;
   string testname;
   if (argc <= 1) {
      yell << "" << endl;
      yell << _CSEPSC << endl;
      shout["fnlo-tk-modify"] << "fastNLO Table Manipulator" << endl;
      yell << _SSEPSC << endl;
      shout["fnlo-tk-modify"] << "Since no argument is given, all settings must be provided via" << endl;
      shout["fnlo-tk-modify"] << "  the default steering file 'SteerModify.str'," << endl;
      shout["fnlo-tk-modify"] << "  which can be consulted for further information as well." << endl;
      shout["fnlo-tk-modify"] << "For more explanations type:" << endl;
      shout["fnlo-tk-modify"] << "./fnlo-tk-modify -h" << endl;
      shout["fnlo-tk-modify"] << "For version number printout type:" << endl;
      shout["fnlo-tk-modify"] << "./fnlo-tk-modify -v" << endl;
      yell << _CSEPSC << endl;
   } else {
      testname = (const char*) argv[1];
      if (testname == "-v") {
         fastNLOTools::PrintFastnloVersion();
         return 0;
      }
      //! --- Print program purpose
      yell << _CSEPSC << endl;
      shout["fnlo-tk-modify"] << "Tool to manipulate a fastNLO table" << endl;
      yell << _SSEPSC << endl;
      shout["fnlo-tk-modify"] << "For more explanations type:" << endl;
      shout["fnlo-tk-modify"] << "./fnlo-tk-modify -h" << endl;
      shout["fnlo-tk-modify"] << "and consult the default steering file 'SteerModify.str'." << endl;
      shout["fnlo-tk-modify"] << "For version number printout type:" << endl;
      shout["fnlo-tk-modify"] << "./fnlo-tk-modify -v" << endl;
      yell << _CSEPSC << endl;
      yell << "" << endl;
      //! --- Usage info
      if (testname == "-h") {
         yell << _CSEPSC << endl;
         shout["fnlo-tk-modify"] << "fastNLO Table Manipulator" << endl;
         yell << _SSEPSC << endl;
         yell << " #" << endl;
         info["fnlo-tk-modify"] << "The purpose of this tool is to allow the user to perform a number of" << endl;
         info["fnlo-tk-modify"] << "modifications on a fastNLO table, e.g. the adaptation of the scenario description." << endl;
         info["fnlo-tk-modify"] << "Most importantly, superfluous observable bins can be removed or additional factors" << endl;
         info["fnlo-tk-modify"] << "can be applied to each bin. Moreover, additional uncertainties like from numerical" << endl;
         info["fnlo-tk-modify"] << "integrations can be added to each perturbative contribution." << endl;
         info["fnlo-tk-modify"] << "Because simple command line arguments or options would not be flexible enough," << endl;
         info["fnlo-tk-modify"] << "it is assumed that the desired changes are set up in a steering file." << endl;
         info["fnlo-tk-modify"] << "By default the steering file is named 'SteerModify.str'." << endl;
         info["fnlo-tk-modify"] << "A template for such a file is contained in the release and" << endl;
         info["fnlo-tk-modify"] << "usually should have been installed under share/fastnlo_toolkit/Modify/." << endl;
         info["fnlo-tk-modify"] << "If another steering filename should be used, this filename can" << endl;
         info["fnlo-tk-modify"] << "be passed via the command line using the tag 'steerfile='," << endl;
         info["fnlo-tk-modify"] << "e.g. fnlo-tk-modify steerfile=AnotherFileName.str" << endl;
         info["fnlo-tk-modify"] << endl;
         man << "" << endl;
         man << "Usage: ./fnlo-tk-modify [steerfile=SteerFile.str] [InTable=fastNLOtableIn.tab.gz] [OutTable=fastNLOtableOut.tab.gz] [OptArg=option]" << endl;
         man << "       Specification: <> mandatory; [] optional." << endl;
         man << "[steerfile=SteerFile.str]:         Alternative steering filename." << endl;
         man << "[InTable=fastNLOtableIn.tab.gz]:   Table input filename, if not specified in steering file." << endl;
         man << "[OutTable=fastNLOtableOut.tab.gz]: Table output filename, if not specified in steering file." << endl;
         man << "   All desired table modifications like" << endl;
         man << "   - changing the scenario name," << endl;
         man << "   - changing or complementing the scenario description," << endl;
         man << "   - cutting out unused bins," << endl;
         man << "   - multiplying bins by a set of factors," << endl;
         man << "   - adding uncertainty information," << endl;
         man << "   - adapting the cross section output units," << endl;
         man << "   - correcting the power of the LO process, or" << endl;
         man << "   - correcting the cms energy" << endl;
         man << "   are assumed to be controlled via steering parameters" << endl;
         man << "   similar to the ones that may be used in table creation." << endl;
         man << "   By default it is expected that at least the mandatory ones are" << endl;
         man << "   either given via command line arguments or are specified" << endl;
         man << "   in a steering file by default named 'SteerModify.str'." << endl;
         man << "   For more steering options please check the default steering file delivered" << endl;
         man << "   by the fastNLO Tolkit (usually in $prefix/share/fastnlo_toolkit/steerfiles)." << endl;
         yell << " #" << endl;
         yell  << _CSEPSC << endl;
         return 0;
      }
      yell  << _CSEPSC << endl;
      info["fnlo-tk-modify"] << "Parsing command line ..." << endl;
      for ( int i = 1; i < argc; i++ ) {
         std::string test = argv[i];
         size_t ipos = test.find("steerfile=");
         if ( ipos == 0 ) {
            ipos = test.find("=");
            steername = test.substr(ipos+1,std::string::npos);
            info["fnlo-tk-modify"] << "Found argument for steerfile: " << steername << endl;
         }
         ipos = test.find("InTable=");
         if ( ipos == 0 ) {
            ipos = test.find("=");
            intable = test.substr(ipos+1,std::string::npos);
            info["fnlo-tk-modify"] << "Found argument for InTable: " << intable << endl;
         }
         ipos = test.find("OutTable=");
         if ( ipos == 0 ) {
            ipos = test.find("=");
            outtable = test.substr(ipos+1,std::string::npos);
            info["fnlo-tk-modify"] << "Found argument for OutTable: " << outtable << endl;
         }
      }
   }

   // Reading all settings from command line; success means provided steering file could be read!
   if ( PARSE(argc,argv) ) {
      info["fnlo-tk-modify"] << "Read all settings from command line and provided steering file." << endl;
      if ( ! (EXIST(InTable) && EXIST(OutTable)) ) {
         error["fnlo-tk-modify"] << "Mandatory parameters InTable and OutTable are not defined," << endl;
         error["fnlo-tk-modify"] << "neither on command line nor in steering file. Aborted!" << endl;
         error["fnlo-tk-modify"] << "For more explanations type:" << endl;
         error["fnlo-tk-modify"] << "./fnlo-tk-modify -h" << endl;
         PRINTALL();
         exit(1);
      }
   } else { // Steering file could not be read, try default one!
      info["fnlo-tk-modify"] << "Provided steering file, if any, could not be read." << endl;
      info["fnlo-tk-modify"] << "Trying default steering file 'SteerModify.str' instead ..." << endl;
      int retcode = READ("SteerModify.str");
      if ( retcode != 0 ) {
         warn["fnlo-tk-modify"] << "No steering file found!" << endl;
         warn["fnlo-tk-modify"] << "All manipulations must have been defined by" << endl;
         warn["fnlo-tk-modify"] << "command line options!" << endl;
      }
      if ( ! (EXIST(InTable) && EXIST(OutTable)) ) {
         error["fnlo-tk-modify"] << "Mandatory parameters InTable and OutTable not defined," << endl;
         error["fnlo-tk-modify"] << "neither on command line nor in any steering file. Aborted!" << endl;
         error["fnlo-tk-modify"] << "For more explanations type:" << endl;
         error["fnlo-tk-modify"] << "./fnlo-tk-modify -h" << endl;
         PRINTALL();
         exit(1);
      }
   }

   //! --- If desired print all steering information
   if ( EXIST(PrintSteeringCard) ) {
      info["fnlo-tk-modify"] << "Print all steering information ..." << endl;
      if ( BOOL(PrintSteeringCard) ) PRINTALL();
   }

   //! --- Reset verbosity level
   if ( EXIST(Verbosity) ) {
      info["fnlo-tk-modify"] << "Resetting verbosity to: " << STRING(Verbosity) << endl;
      SetGlobalVerbosity(toVerbosity()[STRING(Verbosity)]);
   }

   //! --- Print input table information
   info["fnlo-tk-modify"] << "Trying to read input table " << STRING(InTable) << endl;
   fastNLOTable table(CHAR(InTable));
   unsigned int nobs = table.GetNObsBin();
   if ( EXIST(PrintInputA1) ) {
      if ( BOOL(PrintInputA1) ) table.PrintHeader(1);
   }
   if ( EXIST(PrintInputA2) ) {
      if ( BOOL(PrintInputA2) ) table.PrintScenario(1);
   }

   // Block A1: Table header
   if ( EXIST(Itabversion) ) {
      info["fnlo-tk-modify"]<<"Modifying table version: from "<< table.GetITabVersionRead() << " to " << INT(Itabversion) << endl;
      table.SetITabVersionWrite(INT(Itabversion));
   }
   if ( EXIST(ScenName) ) {
      info["fnlo-tk-modify"]<<"Modifying scenario name: from "<< table.GetScenName() << " to " << STRING(ScenName) << endl;
      table.SetScenName(STRING(ScenName));
   }

   // Block A2: Table scenario
   if ( EXIST(Ipublunits) ) {
      info["fnlo-tk-modify"]<<"Modifying publication units: from "<< table.GetIpublunits() << " to " << INT(Ipublunits) << endl;
      table.SetIpublunits(INT(Ipublunits));
   }

   if ( !STRING_ARR(ScDescript).empty() ){
      vector <string> ScDescr = table.GetScDescr();
      size_t NScSize = ScDescr.size();
      info["fnlo-tk-modify"]<<"Modifying existing scenario description:" << endl;
      for ( size_t i = 0; i < NScSize; i++ ) {
         shout << "Line no. " << i << ": " << ScDescr[i] << endl;
      }
      if ( BOOL(AttachScDescription) ){
         info["fnlo-tk-modify"]<<"Attaching to scenario description:" << endl;
         size_t NewNScSize = NScSize + STRING_ARR(ScDescript).size();
         ScDescr.resize(NewNScSize);
         for ( size_t i = NScSize; i < NewNScSize; i++ ) {
            ScDescr[i] = STRING_ARR(ScDescript)[i-NScSize];
            shout << "Line no. " << i << ": " << ScDescr[i] << endl;
         }
      } else {
         info["fnlo-tk-modify"]<<"Replacing scenario description by:" << endl;
         size_t NewNScSize = STRING_ARR(ScDescript).size();
         ScDescr.resize(NewNScSize);
         for ( size_t i = 0; i < NewNScSize; i++ ) {
            ScDescr[i] = STRING_ARR(ScDescript)[i];
            shout << "Line no. " << i << ": " << ScDescr[i] << endl;
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
      for (unsigned int iObs=0; iObs<nobs; iObs++) {
         table.MultiplyBinSize(iObs,fac);
      }
   }

   if ( !DOUBLE_ARR(BinSize).empty() ) {
      vector<double> fac = DOUBLE_ARR(BinSize);
      info["fnlo-tk-modify"]<<"Multiplying bin sizes by provided factors!"<<endl;
      if ( nobs != fac.size() ) {
         error["fnlo-tk-modify"]<<"You need the same number of multiplicative factors, nfact = " << fac.size() <<
            ", than bins in the table, nobsbin = " << nobs << ". Aborted!" << endl;
         exit(1);
      }
      for (unsigned int iObs=0; iObs<nobs; iObs++) {
         table.MultiplyBinSize(iObs,fac[iObs]);
      }
   }

   // Block B's: Scenario contributions
   if ( !DOUBLE_ARR(MultCoeff).empty() ) {
      vector<double> fac = DOUBLE_ARR(MultCoeff);
      info["fnlo-tk-modify"]<<"Multiplying by provided factors all coefficients of additive contributions to observable bins!"<<endl;
      if ( nobs != fac.size() ) {
         error["fnlo-tk-modify"]<<"You need the same number of multiplicative factors, nfact = " << fac.size() <<
            ", than bins in the table, nobsbin = " << nobs << ". Aborted!" << endl;
         exit(1);
      }
      for (unsigned int iObs=0; iObs<nobs; iObs++) {
         table.MultiplyBinInTable(iObs,fac[iObs]);
      }
   }

   //! Erase observable bins from table (Do NOT simultaneously to adding InfoBlocks!)
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

   //! Replace CodeDescription in each perturbative contribution
   if ( !STRING_ARR(CodeDescript).empty() ){
      int Ncontrib = table.GetNcontrib();
      size_t NCodeDescript = STRING_ARR(CodeDescript).size();
      std::vector<std::string> Description;
      info["fnlo-tk-modify"]<<"Replacing code description by:" << endl;
      for ( size_t i = 0; i < NCodeDescript; i++ ) {
         info["fnlo-tk-modify"]<<"Line no. " << i << ": " << STRING_ARR(CodeDescript)[i] << endl;
         Description.push_back(STRING_ARR(CodeDescript)[i]);
      }
      for ( int i = 0; i < Ncontrib; i++ ) {
         fastNLOCoeffBase* c = table.GetCoeffTable(i);
         if ( fastNLOCoeffAddBase::CheckCoeffConstants(c,true) ) {
            c->SetCodeDescription(Description);
         }
      }
   }

   //! Add InfoBlocks with statistical uncertainty from steering file
   //! Default flag1=0: Statistical/numerical uncertainty
   int IBFlag1 = 0;
   //! Default flag2=1: Quadratic addition, alternative: 0: linear addition
   int IBFlag2 = 1;
   //! Default description line
   std::string Default = "Please provide description!";
   if ( EXIST(InfoBlockStatUnc) ) {
      if ( !INT_ARR(RemoveBins).empty() ) {
         info["fnlo-tk-modify"]<<"Do NOT erase bins while adding InfoBlocks or vice versa! Aborted."<<endl;
         exit(25);
      } else {
         info["fnlo-tk-modify"]<<"Adding InfoBlocks to contributions."<<endl;
      }
      static vector<double> dstrel_LO   = DOUBLE_COL(InfoBlockStatUnc,dstrel_LO);
      static vector<double> dstrel_NLO  = DOUBLE_COL(InfoBlockStatUnc,dstrel_NLO);
      static vector<double> dstrel_NNLO = DOUBLE_COL(InfoBlockStatUnc,dstrel_NNLO);
      unsigned int NDescr  = STRING_ARR(InfoBlockDescr).size();
      if ( NDescr > 1 ) {
         error["fnlo-tk-modify"]<<"Only one description line allowed for all blocks, aborted! NDescr = " << NDescr << endl;
         exit(39);
      }
      int Ncontrib = table.GetNcontrib();
      int ic = 0;
      for ( int i = 0; i < Ncontrib; i++ ) {
         fastNLOCoeffBase* c = table.GetCoeffTable(i);
         if ( fastNLOCoeffAddBase::CheckCoeffConstants(c,true) ) {
            int iFlag1 = IBFlag1;
            int iFlag2 = IBFlag2;
            if ( EXIST(InfoBlockFlag1) ) { iFlag1 = INT(InfoBlockFlag1); }
            if ( EXIST(InfoBlockFlag2) ) { iFlag2 = INT(InfoBlockFlag2); }
            if ( c->IsLO() ) {
               info["fnlo-tk-modify"]<<"Found LO contribution " << i << endl;
               std::vector<std::string> Description;
               if ( NDescr > 0 ) {
                  Description.push_back(STRING_ARR(InfoBlockDescr)[0]);
               } else {
                  Description.push_back(Default);
               }
               if ( dstrel_LO.size() == 0 ) {
                  warn["fnlo-tk-modify"]<<"Found LO contribution, but no uncertainties! Nothing added." << endl;
               } else if ( dstrel_LO.size() != nobs ) {
                  error["fnlo-tk-modify"]<<"You need the same number of uncertainties, dstrel_LO = " << dstrel_LO.size() <<
                     ", than bins in the table, nobsbin = " << nobs << ". Aborted!" << endl;
                  exit(1);
               } else {
                  c->AddCoeffInfoBlock(iFlag1,iFlag2,Description,dstrel_LO);
               }
               ic += 1;
            } else if ( c->IsNLO() ) {
               info["fnlo-tk-modify"]<<"Found NLO contribution " << i << endl;
               std::vector <std:: string> Description;
               if ( NDescr > 0 ) {
                  Description.push_back(STRING_ARR(InfoBlockDescr)[0]);
               } else {
                  Description.push_back(Default);
               }
               if ( dstrel_NLO.size() == 0 ) {
                  warn["fnlo-tk-modify"]<<"Found NLO contribution, but no uncertainties! Nothing added." << endl;
               } else if ( dstrel_NLO.size() != nobs ) {
                  error["fnlo-tk-modify"]<<"You need the same number of uncertainties, dstrel_NLO = " << dstrel_NLO.size() <<
                     ", than bins in the table, nobsbin = " << nobs << ". Aborted!" << endl;
                  exit(1);
               } else {
                  c->AddCoeffInfoBlock(iFlag1,iFlag2,Description,dstrel_NLO);
               }
               ic += 1;
            } else if ( c->IsNNLO() ) {
               info["fnlo-tk-modify"]<<"Found NNLO contribution " << i << endl;
               std::vector <std:: string> Description;
               if ( NDescr > 0 ) {
                  Description.push_back(STRING_ARR(InfoBlockDescr)[0]);
               } else {
                  Description.push_back(Default);
               }
               if ( dstrel_NNLO.size() == 0 ) {
                  warn["fnlo-tk-modify"]<<"Found NNLO contribution, but no uncertainties! Nothing added." << endl;
               } else if ( dstrel_NNLO.size() != nobs ) {
                  error["fnlo-tk-modify"]<<"You need the same number of uncertainties, dstrel_NNLO = " << dstrel_NNLO.size() <<
                     ", than bins in the table, nobsbin = " << nobs << ". Aborted!" << endl;
                  exit(1);
               } else {
                  c->AddCoeffInfoBlock(iFlag1,iFlag2,Description,dstrel_NNLO);
               }
               ic += 1;
            } else {
               info["fnlo-tk-modify"]<<"Unknown contribution " << i << endl;
               info["fnlo-tk-modify"]<<"Nothing changed." << endl;
            }
         }
      }
   }

   //! Add InfoBlocks with statistical uncertainty from file (NNLOJET .dat, fnlo-tk-statunc .log, or .txt)
   else if ( !STRING_ARR(InfoBlockFiles).empty() &&
             !STRING_ARR(InfoBlockOrders).empty() ) {
      if ( !INT_ARR(RemoveBins).empty() ) {
         info["fnlo-tk-modify"]<<"Do NOT erase bins while adding InfoBlocks or vice versa! Aborted."<<endl;
         exit(25);
      } else {
         info["fnlo-tk-modify"]<<"Adding InfoBlocks to contributions."<<endl;
      }
      unsigned int NFiles  = STRING_ARR(InfoBlockFiles).size();
      unsigned int NCols   = INT_ARR(InfoBlockFileColumns).size();
      unsigned int NOrders = STRING_ARR(InfoBlockOrders).size();
      unsigned int NDescr  = STRING_ARR(InfoBlockDescr).size();
      if ( NFiles != NOrders ) {
         error["fnlo-tk-modify"]<<"Need one order specification per file, aborted! Found NFiles = " << NFiles << ", and NOrders = " << NOrders <<endl;
         exit(37);
      }
      unsigned int icola = 0;
      unsigned int icolb = 0;
      if ( NCols == 0 ) {
      } else if ( NCols == 1 ) {
         icola = INT_ARR(InfoBlockFileColumns)[0];
      } else if ( NCols == 2 ) {
         icola = INT_ARR(InfoBlockFileColumns)[0];
         icolb = INT_ARR(InfoBlockFileColumns)[1];
      } else {
         error["fnlo-tk-modify"]<<"Up to two column numbers allowed, but found more. Aborted! NCols = " << NCols <<endl;
         exit(38);
      }
      if ( NDescr > 1 ) {
         error["fnlo-tk-modify"]<<"Only one description line allowed for all blocks, aborted! NDescr = " << NDescr << endl;
         exit(39);
      }
      for ( unsigned int i = 0; i < NFiles; i++ ){
         info["fnlo-tk-modify"]<<"InfoBlock file no. " << i << " is: " << STRING_ARR(InfoBlockFiles)[i] << endl;
      }
      for ( unsigned int i = 0; i < NOrders; i++ ){
         info["fnlo-tk-modify"]<<"InfoBlock order no. " << i << " is: " << STRING_ARR(InfoBlockOrders)[i] << endl;
      }
      int Ncontrib = table.GetNcontrib();
      int ic = 0;
      std::string Default = "Please provide description!";
      for ( int i = 0; i < Ncontrib; i++ ) {
         fastNLOCoeffBase* c = table.GetCoeffTable(i);
         if ( fastNLOCoeffAddBase::CheckCoeffConstants(c,true) ) {
            int iFlag1 = IBFlag1;
            int iFlag2 = IBFlag2;
            if ( EXIST(InfoBlockFlag1) ) { iFlag1 = INT(InfoBlockFlag1); }
            if ( EXIST(InfoBlockFlag2) ) { iFlag2 = INT(InfoBlockFlag2); }
            if ( c->IsLO() ) {
               info["fnlo-tk-modify"]<<"Found LO contribution " << i << endl;
               int ilo = 0;
               if ( NOrders > 1 ) {
                  for ( unsigned int j = 0; j < NFiles; j++ ) {
                     if (STRING_ARR(InfoBlockOrders)[j] == "LO") {
                        ilo = j;
                     }
                  }
               }
               std::vector<std::string> Description;
               if ( NDescr > 0 ) {
                  Description.push_back(STRING_ARR(InfoBlockDescr)[0]);
               } else {
                  Description.push_back(Default);
               }
               c->AddCoeffInfoBlock(iFlag1,iFlag2,Description,STRING_ARR(InfoBlockFiles)[ilo],icola,icolb);
               ic += 1;
            } else if ( c->IsNLO() ) {
               info["fnlo-tk-modify"]<<"Found NLO contribution " << i << endl;
               int inlo = 0;
               if ( NOrders > 1 ) {
                  for ( unsigned int j = 0; j < NFiles; j++ ) {
                     if (STRING_ARR(InfoBlockOrders)[j] == "NLO") {
                        inlo = j;
                     }
                  }
               }
               std::vector <std:: string> Description;
               if ( NDescr > 0 ) {
                  Description.push_back(STRING_ARR(InfoBlockDescr)[0]);
               } else {
                  Description.push_back(Default);
               }
               c->AddCoeffInfoBlock(iFlag1,iFlag2,Description,STRING_ARR(InfoBlockFiles)[inlo],icola,icolb);
               ic += 1;
            } else if ( c->IsNNLO() ) {
               info["fnlo-tk-modify"]<<"Found NNLO contribution " << i << endl;
               int innlo = 0;
               if ( NOrders > 1 ) {
                  for ( unsigned int j = 0; j < NFiles; j++ ) {
                     if (STRING_ARR(InfoBlockOrders)[j] == "NNLO") {
                        innlo = j;
                     }
                  }
               }
               std::vector <std:: string> Description;
               if ( NDescr > 0 ) {
                  Description.push_back(STRING_ARR(InfoBlockDescr)[0]);
               } else {
                  Description.push_back(Default);
               }
               c->AddCoeffInfoBlock(iFlag1,iFlag2,Description,STRING_ARR(InfoBlockFiles)[innlo],icola,icolb);
               ic += 1;
            } else {
               info["fnlo-tk-modify"]<<"Unknown contribution " << i << endl;
               info["fnlo-tk-modify"]<<"Nothing changed." << endl;
            }
         }
      }
   }

   //! --- Print output table information
   if ( EXIST(PrintOutputA1) ) {
      if ( BOOL(PrintOutputA1) ) table.PrintHeader(0);
   }

   if ( EXIST(PrintOutputA2) ) {
      if ( BOOL(PrintOutputA2) ) table.PrintScenario(0);
   }

   //! Writing modified table
   table.SetFilename(CHAR(OutTable));
   table.WriteTable();

   // Finished
   info["fnlo-tk-modify"]<<"Wrote modified table to " << STRING(OutTable) << endl;

   return 0;
}
