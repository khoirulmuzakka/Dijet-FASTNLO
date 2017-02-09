///********************************************************************
///
///     fastNLO_toolkit: fnlo-tk-merge2
///     Tool to merge fastNLO tables with different contributions or
///     to combine identical statistically independent contributions
///
///     fnlo-tk-merge2 makes use of additional information on event
///     weights, like sumw2, sumsig2, etc... and allows to chose
///     various options for weighting the tables.
///
///     For more explanations type:
///     ./fnlo-tk-merge2 -h
///
///********************************************************************
// DB, 07.02.17

#include <cstdlib>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <unistd.h>
#include "fastnlotk/fastNLOCoeffAddBase.h"
#include "fastnlotk/fastNLOTable.h"
#include "fastnlotk/speaker.h"

std::map<std::string,std::string> _validoptions{
   {"-f","force.  Force to (overwrite) output table."},
   {"-h","help.   Print this message."},
   {"-p","Plot.   -p filename.ps. Plot some statistics of all input files. Option must be followed by filename."},
   {"-1","Once.   Read files only once, and then keep all in memory at the same time."},
   {"-w","Weight.  -w <option>. Calculate (un)weighted average. Options: GenWgt (default), unweighted, append, median, mean, NumEvt, SumW2, SumSig2, SumSig2BinProc, NumEvtBinProc or SumW2BinProc."}
};

//__________________________________________________________________________________________________________________________________
void PrintHelpMessage() {
   using namespace std;
   using namespace say;          //! namespace for 'speaker.h'-verbosity levels
   using namespace fastNLO;      //! namespace for fastNLO constants
   man << "" << endl;
   man << "Run:  $ ./fnlo-tk-merge2 [options] <InTable_1.tab> [InTable_n.tab] <OutTable.tab>" << endl;
   man << "" << endl;
   man << "  <InTable_1.tab>:   First table input file to be merged" << endl;
   man << "  <OutTable.tab>:    Output filename, to which the merged table is written" << endl;
   man << ""<<endl;
   man << "  This program essentially takes one input table and calls fastNLOTable::MergeTable() for each other table. "<<endl;
   man << " "<<endl;
   for ( auto iop : _validoptions ) {
      man << "  "<<iop.first<<"\t\t"<<iop.second<<endl;
   }
   man << " " <<endl;
}

//__________________________________________________________________________________________________________________________________
std::vector<std::vector<std::vector<double> > > CalculateWeight(const std::vector<fastNLO::WgtStat >& allWgtStat, const std::string& wgtoption) {
   using namespace std;
   using namespace say;
   if ( allWgtStat.empty() ) exit(3);

   //! calculate weight for weighted average.
   unsigned int nTab  = allWgtStat.size();
   unsigned int nProc = allWgtStat[0].WgtObsNumEv.size();
   unsigned int nBin  = allWgtStat[0].WgtObsNumEv[0].size();
   
   vector<vector<vector<double> > > ret(nTab);
   for ( unsigned int iTab = 0 ; iTab<nTab ; iTab++ ) { // loop over all input tables.
      ret[iTab].resize(nProc); // proc
      for ( unsigned int iProc = 0 ; iProc < nProc ; iProc++ ) { // proc
	 ret[iTab][iProc].resize(nBin); // proc
      }
   }
   // struct fastNLO::WgtStat
   // vector<double> allWgtNevt;
   // vector<unsigned long long> allWgtNumEv;
   // vector<double> allWgtSumW2; 
   // vector<double> allSigSumW2; 
   // vector<double> allSigSum; 
   // vector<fastNLO::v2d> allWgtObsSumW2;
   // vector<fastNLO::v2d> allSigObsSumW2;
   // vector<fastNLO::v2d> allSigObsSum;  
   // vector<std::vector < std::vector < unsigned long long > > > allWgtObsNumEv; 

   if ( wgtoption=="median" ) return ret;

   const double wgt00 = -1000;
   vector<vector<double> > WgtPerProc(nTab);
   vector<vector<double> > WgtPerBin(nTab);
   for ( unsigned int iTab = 0 ; iTab<nTab ; iTab++ ) { // loop over all input tables.
      double wgt=wgt00;
      // --- one weight for an entire table
      if ( wgtoption.find("tot") != string::npos ) {
	 if ( wgtoption.find("NumEvt")  != string::npos )      wgt = allWgtStat[iTab].WgtNumEv;
	 else if ( wgtoption.find("SumW2")   != string::npos ) wgt = allWgtStat[iTab].WgtSumW2;
	 else if ( wgtoption.find("SumSig2") != string::npos ) wgt = allWgtStat[iTab].SigSumW2;      
      }
      else if ( wgtoption == "GenWgt" ) wgt = allWgtStat[iTab].WgtNevt;
      else if ( wgtoption == "unweighted" ) wgt = 1;
      else if ( wgtoption == "append" ) wgt = -1;
      if ( wgt == 0 ) {
	 error["CalculateWeight"]<<"The requested weight '"<<wgtoption<<"' is zero and thus not available in the table."<<endl;
	 exit(3);
      }
      else if ( wgt != wgt00 ) {
	 for ( unsigned int iProc = 0 ; iProc < nProc ; iProc++ ) { // proc
	    for ( unsigned int iBin = 0 ; iBin < nBin ; iBin++ ) { // ObsBin
	       ret[iTab][iProc][iBin]  = wgt;
	    }
	 }
	 return ret; // done !
      }
      // --- different weights for a single table
      else {
	 WgtPerProc[iTab].resize(nProc);
	 WgtPerBin [iTab].resize(nBin);
	 for ( unsigned int iProc = 0 ; iProc < nProc ; iProc++ ) { // proc
	    for ( unsigned int iBin = 0 ; iBin < nBin ; iBin++ ) { // ObsBin
	       if      ( wgtoption.find("NumEvt")  != string::npos ) wgt = allWgtStat[iTab].WgtObsNumEv[iProc][iBin] ;
	       else if ( wgtoption.find("SumW2")   != string::npos ) wgt = allWgtStat[iTab].WgtObsSumW2[iProc][iBin] ;
	       else if ( wgtoption.find("SumSig2") != string::npos ) wgt = allWgtStat[iTab].SigObsSumW2[iProc][iBin] ;
	       else {
		  say::error["CalculateWeight"]<<"Unrecognized 'weight' option: "<<wgtoption<<endl; exit(2); 
	       }
	       if ( wgt == 0 ) {
		  error["CalculateWeight"]<<"The requested weight '"<<wgtoption<<"' is zero for proc="<<iProc<<" and bin="<<iBin<<" and thus not available in the table."<<endl;
		  exit(3);
	       }
	       ret       [iTab][iProc][iBin]   = wgt;
	       WgtPerProc[iTab][iProc]        += wgt;
	       WgtPerBin [iTab][iBin]         += wgt;
	    }
	 }
	 if ( wgtoption.find("BinProc") ) return ret; // done!
	 else if ( wgtoption.find("Bin") != string::npos ) {
	    for ( unsigned int iProc = 0 ; iProc < nProc ; iProc++ ) { // proc
	       for ( unsigned int iBin = 0 ; iBin < nBin ; iBin++ ) { // ObsBin
		  ret  [iTab][iProc][iBin] = WgtPerBin[iTab][iBin];// 
	       }
	    }
	 }
	 else if ( wgtoption.find("Proc") != string::npos ) {
	    for ( unsigned int iProc = 0 ; iProc < nProc ; iProc++ ) { // proc
	       for ( unsigned int iBin = 0 ; iBin < nBin ; iBin++ ) { // ObsBin
		  ret  [iTab][iProc][iBin] = WgtPerProc[iTab][iProc];// 
	       }
	    }
	 }
	 else {
	    error["CalculateWeight"]<<"Unrecognized weight option: "<<wgtoption<<endl;
	    exit(3);
	 }
	 return ret;
      }
   }
   return ret;
}

//__________________________________________________________________________________________________________________________________
int main(int argc, char** argv) {

   //! --- namespaces
   using namespace std;
   using namespace say;          //! namespace for 'speaker.h'-verbosity levels
   using namespace fastNLO;      //! namespace for fastNLO constants

   //! --- Set verbosity level
   SetGlobalVerbosity(INFO);

   //! --- Print program purpose
   yell << _CSEPSC << endl;
   info["fnlo-tk-merge2"] << "Tool to merge fastNLO tables with different contributions or" << endl;
   info["fnlo-tk-merge2"] << "to combine identical statistically independent contributions" << endl;
   yell << _SSEPSC << endl;
   info["fnlo-tk-merge2"] << "For more explanations type:" << endl;
   info["fnlo-tk-merge2"] << "./fnlo-tk-merge2 -h" << endl;
   yell << _CSEPSC << endl;

   //! --- Parse commmand line
   yell << "" << endl;
   yell << _CSEPSC << endl;
   shout["fnlo-tk-merge2"] << "fastNLO Table Merger"<<endl;
   yell << _SSEPSC << endl;
   if (argc <= 1) {
      error["fnlo-tk-merge2"] << "No arguments given, but need at least two!" << endl;
      shout["fnlo-tk-merge2"] << "For an explanation of command line arguments use option -h" << endl;
      PrintHelpMessage();
      exit(1);
   }
   // } else {
   //    tablename = (const char*) argv[1];
   //    //! --- Usage info
   //    if (tablename == "-h") {
   //       yell << " #" << endl;
   //       info["fnlo-tk-merge2"] << "The purpose of this tool is to merge fastNLO tables with different contributions or" << endl;
   //       info["fnlo-tk-merge2"] << "to combine identical statistically independent additive contributions to improve" << endl;
   //       info["fnlo-tk-merge2"] << "the statistical precision." << endl;
   //       info["fnlo-tk-merge2"] << "The statistical information of each additive contribution is checked." << endl;
   //       info["fnlo-tk-merge2"] << "An event number of unity indicates that this contribution" << endl;
   //       info["fnlo-tk-merge2"] << "has been combined from multiple contributions losing the" << endl;
   //       info["fnlo-tk-merge2"] << "the event normalisation information that is stored contribution-" << endl;
   //       info["fnlo-tk-merge2"] << "and not bin-wise. Further merging of such tables is not possible." << endl;
   //       yell << " #" << endl;
   //       yell << _CSEPSC << endl;
   //       return 0;
   //    }
   // }

   set<string> setfiles;
   vector<string> files;
   set<string> options;
   string outfile = argv[argc-1];
   string wgtoption = "GenWgt";
   string plotfile = "fnlo-tk-merge2.ps";
   
   for (int iarg=1; iarg<argc-1; iarg++) {
      string sarg = argv[iarg];
      if ( sarg.find(".tab") != string::npos ) {
	 if (access(sarg.c_str(), R_OK) != 0) //! --- File there?
	    warn["fnlo-tk-merge2"]<<"Unable to access file '" << sarg << "', skipped!" << endl;
	 else { // keep it
	    if ( setfiles.count(sarg) ) 
	       warn["fnlo-tk-merge2"]<<"File '"<<sarg<<"' already added once (but duplication is allowed)."<<endl;
	    files.push_back(sarg);
	    setfiles.insert(sarg);
	 }
      }
      else if ( sarg.find("-") == 0 ) {
	 options.insert(sarg);
	 if ( sarg == "-p" ) plotfile=argv[++iarg];
	 if ( sarg == "-w" ) wgtoption=argv[++iarg];
	 if ( _validoptions.count(sarg) == 0 ) { 
	    error["fnlo-tk-merge2"]<<"Invalid option '"<<sarg<<"'."<<endl<<endl;;
	    PrintHelpMessage();
	    exit(1);
	 }
      }
      else { //error
	 error["fnlo-tk-merge2"]<<"Input argument not valid: '"<<sarg<<"'. Only in-/out-filenames or options (-XYZ) allowed."<<endl;
	 exit(1);
      }
   }
   // --- help message
   if ( options.count("-h") || outfile=="-h" ){  PrintHelpMessage(); return 0; }
   // output files
   if ( outfile.find(".tab") == string::npos ) {
      error["fnlo-tk-merge2"]<<"Last argument must be output file (containing '.tab')"<<endl;
      exit(1);
   }
   if (access(outfile.c_str(), R_OK) == 0) {
      if ( options.count("-f") ) {
	 warn["fnlo-tk-merge2"]<<"Output file " << outfile << " exists already. Overwriting it."<<endl;
      }
      else {
	 error["fnlo-tk-merge2"]<<"Output file " << outfile << " exists already!" << endl;
	 shout["fnlo-tk-merge2"]<<"Please remove it first." << endl;
	 exit(1);
      }
   }
   //! --- Check no. of file names
   if (files.empty()) {
      error["fnlo-tk-merge2"] << "No input filenames given, need at least one!" << endl;
      exit(1);
   }
   // i/o done
   // --------------------------------------------------------------------------------------- //

   info["fnlo-tk-merge2"]<<"Using weighting option: "<<wgtoption<<"."<<endl;
   fastNLOTable::EMerge moption = fastNLOTable::kUndefined;
   if      ( wgtoption=="GenWgt" ) moption = fastNLOTable::kMerge ;
   else if ( wgtoption=="append" ) moption = fastNLOTable::kAppend ;
   else if ( wgtoption=="unweighted" ) moption = fastNLOTable::kUnweighted ;
   else if ( wgtoption=="median" ) moption = fastNLOTable::kMedian ;
   else if ( wgtoption=="mean" ) moption = fastNLOTable::kMean ;
   else if ( wgtoption=="NumEvt" ) moption = fastNLOTable::kNumEvent ;
   else if ( wgtoption=="SumW2" ) moption = fastNLOTable::kSumW2 ;
   else if ( wgtoption=="SumSig2" ) moption = fastNLOTable::kSumSig2 ;
   else if ( wgtoption=="NumEvtBinProc" ) moption = fastNLOTable::kNumEventBinProc ;
   else if ( wgtoption=="SumW2BinProc" ) moption = fastNLOTable::kSumW2BinProc ;
   else if ( wgtoption=="SumSig2BinProc" ) moption = fastNLOTable::kSumSig2BinProc ;
   else moption = fastNLOTable::kUndefined;

   if ( moption==fastNLOTable::kUndefined ) {
      error["fnlo-tk-merge2"]<<"Cannot recognize merge option: "<<wgtoption<<endl;
      exit(1);
   }
   if ( moption==fastNLOTable::kMedian || moption == fastNLOTable::kMean ) {
      error["fnlo-tk-merge2"]<<"Option median or median not yet implemented."<<endl;
      exit(1);
   }

   // --------------------------------------------------------------------------------------- //
   // --- calculate statistics for mergeing/normalisations
   // --------------------------------------------------------------------------------------- //
   map<string,fastNLOTable*> alltables;
   map<string,unsigned int> lookup;
   // all variables
   vector<fastNLO::WgtStat > allWgtStat(files.size());
   
   for ( auto path : files ) {
      info["fnlo-tk-merge"]<<"Reading table '" << path << "'" << endl;
      fastNLOTable* tab = new fastNLOTable(path);

      if ( tab->GetItabversion() < 23000 ) {
	 warn["fnlo-tk-merge2"]<<"This program is only compatible with fastNLO table versions > v2.3,000."<<endl;
	 warn["fnlo-tk-merge2"]<<"Please use program fnlo-tk-merge and/or fnlo-tk-append instead."<<endl;
	 //exit(2);
      }

      alltables[path] = tab;
      int tabid = alltables.size()-1;
      lookup[path] = tabid;
      
      //! --- Check statistical information of additive contributions
      if ( options.count("-p") ) {
	 int nc = tab->GetNcontrib() + tab->GetNdata();
	 if ( nc != 1 ) {
	    warn["fnlo-tk-merge2"]<<"Program fnlo-tk-merge2 currently can only handle fastNLO tables with exactly one contributions. Please use program fnlo-tk-merge fnlo-tk-append."<<endl;
	    //exit(2);
	 }
	 for ( int ic=0 ; ic<nc; ic++ ) {
	    bool quiet = true;
	    fastNLOCoeffBase* cnew = (fastNLOCoeffBase*)tab->GetCoeffTable(ic);
	    // Identify type of new coeff table, only additive ones have event numbers
	    if ( fastNLOCoeffAddBase::CheckCoeffConstants(cnew,quiet) ) { // additive?
	       fastNLOCoeffAddBase* cadd = (fastNLOCoeffAddBase*)cnew;
	       if ( cadd->GetNevt() == 1 ) {
		  warn["fnlo-tk-merge"]<<"Contribution #" << ic << " in table " << path << " has event number '1', which is usually invalid."<<endl;
		  // error["fnlo-tk-merge"]<<"Contribution #" << ic << " in table " << path << endl;
		  // error["fnlo-tk-merge"]<<"has no valid number-of-events information and cannot be merged. Aborted!" << endl;
		  // error["fnlo-tk-merge"]<<"Nevt = " << cadd->GetNevt() << endl;
		  //exit(1);
	       }
	       //
	       allWgtStat[tabid] = cadd->GetWgtStat();
	    }
	    else {
	       error["fnlo-tk-merge2"]<<"Program fnlo-tk-merge2 can only deal with additive contributions. Exiting"<<endl;
	       exit(2);
	    }
	 }
	 if ( options.count("-1") ) { 
	    error<<"Option -1 not yet fully implemented and tested."<<endl;
	    //delete tab; 
	    //alltables.erase(path); 
	 } 
      }
   }
   
   // --- calculate weight
   if ( options.count("-p") ) {
      vector<vector<vector<double> > > ProcBinWgt = CalculateWeight(allWgtStat,wgtoption);
   }
   // --------------------------------------------------------------------------------------- //
   //! --- Initialize output table
   fastNLOTable* resultTable = NULL;

   // --------------------------------------------------------------------------------------- //
   //! --- Loop input files and merge them
   // --------------------------------------------------------------------------------------- //
   int nValidTables = 0;
   for ( auto path : files ) {
      fastNLOTable* tab = alltables[path]; // fastNLOTable tab(path);
      if ( tab==NULL ) {
	 warn["fnlo-tk-merge2"]<<"Skipping table "<<path<<", because it was not found in list of valid tables."<<endl;
	 continue;
      }

      // --- Initialising result with first read table
      if ( resultTable==NULL ) {
	 info["fnlo-tk-merge"]<<"Initialising result table '" << outfile << "'" << endl;
	 //resultTable = new fastNLOTable(*tab);
	 resultTable = tab;
	 nValidTables++;
	 continue;
      }

      // --- Adding further tables to result table
      resultTable->MergeTable(*tab, moption);
      nValidTables++;
   }

   info["fnlo-tk-merge2"]<<"Merged "<<nValidTables<<" table file(s)."<<endl;
   if (nValidTables == 0 ) {
      error["fnlo-tk-merge"]<<"Found less than two valid tables, no output possible!"<<endl;
      return 0;
   }

   // --- Write result
   resultTable->SetFilename(outfile);
   info["fnlo-tk-merge"]<<"Write merged results to file "<<resultTable->GetFilename()<<"."<<endl;
   resultTable->WriteTable();

   // --- fine
   return 0;
}
