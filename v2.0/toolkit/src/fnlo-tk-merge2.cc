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
   {"-f","force     Force to (overwrite) output table."},
   {"-h","help      Print this message."},
   {"-o","output    -o <output>. Specify output file name (the last argument is then considered to be an input table)."},
//   {"-p","Plot.   -p <filename>. Plot some statistics of all input files. Option must be followed by filename."},
//   {"-1","Once.   Read files only once, and then keep all in memory at the same time."},
   {"-w","Weight    -w <option>. Calculate (un)weighted average. Option: GenWgt (default), unweighted, append, median, mean, NumEvt, SumW2, SumSig2, SumSig2BinProc, NumEvtBinProc or SumW2BinProc."},
   {"-pre","pre-avg -pre <n> <option>. 2 step mergeing: Build pre-averaged of n tables using weighting procedure <option>."},
};


std::map<std::string,fastNLO::EMerge> _wgtoptions {
   {"GenWgt",fastNLO::kMerge },
   {"append",fastNLO::kAppend },
   {"unweighted",fastNLO::kUnweighted },
   {"median",fastNLO::kMedian },
   {"mean",fastNLO::kMean },
   {"NumEvt",fastNLO::kNumEvent },
   {"SumW2",fastNLO::kSumW2 },
   {"SumSig2",fastNLO::kSumSig2 },
   {"NumEvtBinProc",fastNLO::kNumEventBinProc },
   {"SumW2BinProc",fastNLO::kSumW2BinProc },
   {"SumSig2BinProc",fastNLO::kSumSig2BinProc },
//   {"",fastNLO::kUndefined}
};

   //    fastNLOTable::EMerge moption = fastNLO::kUndefined;
   // if      ( wgtoption=="GenWgt" ) moption = fastNLO::kMerge ;
   // else if ( wgtoption=="append" ) moption = fastNLO::kAppend ;
   // else if ( wgtoption=="unweighted" ) moption = fastNLO::kUnweighted ;
   // else if ( wgtoption=="median" ) moption = fastNLO::kMedian ;
   // else if ( wgtoption=="mean" ) moption = fastNLO::kMean ;
   // else if ( wgtoption=="NumEvt" ) moption = fastNLO::kNumEvent ;
   // else if ( wgtoption=="SumW2" ) moption = fastNLO::kSumW2 ;
   // else if ( wgtoption=="SumSig2" ) moption = fastNLO::kSumSig2 ;
   // else if ( wgtoption=="NumEvtBinProc" ) moption = fastNLO::kNumEventBinProc ;
   // else if ( wgtoption=="SumW2BinProc" ) moption = fastNLO::kSumW2BinProc ;
   // else if ( wgtoption=="SumSig2BinProc" ) moption = fastNLO::kSumSig2BinProc ;
   // else moption = fastNLO::kUndefined;


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
   info["fnlo-tk-merge2"] << "Tool to merge fastNLO tables with different contributions or" << endl;
   info["fnlo-tk-merge2"] << "to combine identical statistically independent contributions" << endl;
   info["fnlo-tk-merge2"] << "For more explanations type:" << endl;
   info["fnlo-tk-merge2"] << "  $  fnlo-tk-merge2 -h" << endl<<endl;

   //! --- Parse commmand line
   if (argc <= 1) {
      error["fnlo-tk-merge2"] << "No arguments given, but need at least two!" << endl;
      PrintHelpMessage();   exit(1);
   }
   set<string> setfiles;
   vector<string> files;
   set<string> options;
   string outfile = argv[argc-1];
   string wgtoption = "GenWgt";
   string plotfile = "fnlo-tk-merge2.ps";
   int pre = 0;
   string preoptin = "NumEvtBinProc"; // NumEvt

   int narg = argc-1;
   for (int iarg=1; iarg<narg; iarg++) {
      string sarg = argv[iarg];
      if ( sarg.find(".tab") != string::npos ) {
	 if (access(sarg.c_str(), R_OK) != 0) //! --- File there?
	    warn["fnlo-tk-merge2"]<<"Unable to access file '" << sarg << "', skipped!" << endl;
	 else { // keep it
	    if ( setfiles.count(sarg) ) {
	       warn["fnlo-tk-merge2"]<<"File '"<<sarg<<"' already added once (but duplication is allowed)."<<endl;
	    }
	    files.push_back(sarg);
	    setfiles.insert(sarg);
	 }
      }
      else if ( sarg.find("-") == 0 ) {
	 if ( options.count(sarg) ) {error["fnlo-tk-merge2"]<<"Duplicate option "<<sarg<<" recognized. Exiting."<<endl;exit(1); }
	 options.insert(sarg);
	 if ( sarg == "-p" ) plotfile=argv[++iarg];
	 if ( sarg == "-w" ) wgtoption=argv[++iarg];
	 if ( sarg == "-o" ) { outfile=argv[++iarg]; narg++; }
	 if ( sarg == "-pre" ) { 
	    pre=atoi(argv[++iarg]); 
	    if ( _wgtoptions.count(argv[iarg+1]) ) preoptin=argv[++iarg]; 
	    else info["fnlo-tk-merge2"]<<"Using pre-average weighting default: "<<preoptin<<endl;
	 }
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

   if ( _wgtoptions.count(wgtoption)==0 ) {
      error["fnlo-tk-merge2"]<<"Cannot recognize merge option: "<<wgtoption<<endl;
      exit(1);
   }
   fastNLO::EMerge moption = _wgtoptions[wgtoption];
   fastNLO::EMerge prewgt  = _wgtoptions[preoptin];// alrady checked



   // --------------------------------------------------------------------------------------- //
   // --- calculate statistics for mergeing/normalisations
   // --------------------------------------------------------------------------------------- //
   // DB: Currently nothing is happening in the following code block
   //     Updates for option '-1' are needed as well...
   map<string,fastNLOTable*> alltables;
   map<string,unsigned int> lookup;
   vector<fastNLO::WgtStat > allWgtStat(files.size());
   for ( auto path : files ) {
      info["fnlo-tk-merge"]<<"Reading table '" << path << "'" << endl;
      fastNLOTable* tab = new fastNLOTable(path);

      if ( tab->GetItabversion() < 23000 ) {
	 warn["fnlo-tk-merge2"]<<"This program is maybe only compatible with fastNLO table versions > v2.3,000."<<endl;
	 warn["fnlo-tk-merge2"]<<"Consider to use program fnlo-tk-merge and/or fnlo-tk-append instead."<<endl;
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
   
   // --------------------------------------------------------------------------------------- //
   // --- calculate statistics
   if ( options.count("-p") ) {
      vector<vector<vector<double> > > ProcBinWgt = CalculateWeight(allWgtStat,wgtoption);
   }

   
   // --------------------------------------------------------------------------------------- //
   // --- calculating pre-averages if requested.
   if ( pre > 1 ) {
      fastNLOTable* keeptab = NULL;
      vector<fastNLOTable*> addtab;
      vector<string> remainingfiles;
      for ( auto path : files ) {
	 fastNLOTable* tab = alltables[path]; // fastNLOTable tab(path);
	 if ( keeptab == NULL ) {
	    remainingfiles.push_back(path);
	    keeptab = tab;
	    continue;
	 }
	 else {
	    addtab.push_back(tab);
	    alltables.erase(path);
	 }
	 if ( int(addtab.size()) == pre-1 || files.back() == path) { // merge
	    keeptab->MergeTables(addtab,prewgt);
	    addtab.clear();
	    keeptab = NULL;
	 }
      }
      files=remainingfiles;
   }

   // --------------------------------------------------------------------------------------- //
   //! --- Initialize output table
   fastNLOTable* resultTable = NULL;

   // --------------------------------------------------------------------------------------- //
   //! --- Loop input files and merge them
   // --------------------------------------------------------------------------------------- //
   int nValidTables = 0;
   vector<fastNLOTable*> allvec;
   for ( auto path : files ) {
      if ( alltables.count(path)==0 && pre<=1 ) {
	 warn["fnlo-tk-merge2"]<<"Skipping table "<<path<<", because it was not found in list of valid tables."<<endl;
	 continue;
      }
      fastNLOTable* tab = alltables[path]; // fastNLOTable tab(path);

      // --- Initialising result with first read table
      if ( resultTable==NULL ) {
	 info["fnlo-tk-merge"]<<"Initialising result table '" << outfile << "'" << endl;
	 //resultTable = new fastNLOTable(*tab);
	 resultTable = tab;
	 nValidTables++;
	 continue;
      }

      // --- Adding further tables to result table
      if ( moption==fastNLO::kMedian || moption == fastNLO::kMean )  
	 allvec.push_back(tab); // prepare for mergeing
      else 
	 resultTable->MergeTable(*tab, moption); // merge
      nValidTables++;
   }

   // merge, if not yet done so
   if ( moption == fastNLO::kMedian || moption == fastNLO::kMean ) {
      resultTable->MergeTables(allvec,moption);
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
