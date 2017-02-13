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

std::string _progname = "fnlo-tk-merge2";

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
   info[_progname] << "Tool to merge fastNLO tables with different contributions or" << endl;
   info[_progname] << "to combine identical statistically independent contributions" << endl;
   info[_progname] << "For more explanations type:" << endl;
   info[_progname] << "  $  fnlo-tk-merge2 -h" << endl<<endl;

   //! --- Parse commmand line
   if (argc <= 1) {
      error[_progname] << "No arguments given, but need at least two!" << endl;
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
	    warn[_progname]<<"Unable to access file '" << sarg << "', skipped!" << endl;
	 else { // keep it
	    if ( setfiles.count(sarg) ) {
	       warn[_progname]<<"File '"<<sarg<<"' already added once (but duplication is allowed)."<<endl;
	    }
	    files.push_back(sarg);
	    setfiles.insert(sarg);
	 }
      }
      else if ( sarg.find("-") == 0 ) {
	 if ( options.count(sarg) ) {error[_progname]<<"Duplicate option "<<sarg<<" recognized. Exiting."<<endl;exit(1); }
	 options.insert(sarg);
	 if ( sarg == "-p" ) plotfile=argv[++iarg];
	 if ( sarg == "-w" ) wgtoption=argv[++iarg];
	 if ( sarg == "-o" ) { outfile=argv[++iarg]; narg++; }
	 if ( sarg == "-pre" ) { 
	    pre=atoi(argv[++iarg]); 
	    if ( _wgtoptions.count(argv[iarg+1]) ) preoptin=argv[++iarg]; 
	    else info[_progname]<<"Using pre-average weighting default: "<<preoptin<<endl;
	 }
	 if ( _validoptions.count(sarg) == 0 ) { 
	    error[_progname]<<"Invalid option '"<<sarg<<"'."<<endl<<endl;;
	    PrintHelpMessage();
	    exit(1);
	 }
      }
      else { //error
	 error[_progname]<<"Input argument not valid: '"<<sarg<<"'. Only in-/out-filenames or options (-XYZ) allowed."<<endl;
	 exit(1);
      }
   }
   // --- help message
   if ( options.count("-h") || outfile=="-h" ){  PrintHelpMessage(); return 0; }
   // output files
   if ( outfile.find(".tab") == string::npos ) {
      error[_progname]<<"Last argument must be output file (containing '.tab')"<<endl;
      exit(1);
   }
   if (access(outfile.c_str(), R_OK) == 0) {
      if ( options.count("-f") ) {
	 warn[_progname]<<"Output file " << outfile << " exists already. Overwriting it."<<endl;
      }
      else {
	 error[_progname]<<"Output file " << outfile << " exists already!" << endl;
	 shout[_progname]<<"Please remove it first." << endl;
	 exit(1);
      }
   }
   //! --- Check no. of file names
   if (files.empty()) {
      error[_progname] << "No input filenames given, need at least one!" << endl;
      exit(1);
   }
   // i/o done
   // --------------------------------------------------------------------------------------- //
   
   info[_progname]<<"Using weighting option: "<<wgtoption<<"."<<endl;

   if ( _wgtoptions.count(wgtoption)==0 ) {
      error[_progname]<<"Cannot recognize merge option: "<<wgtoption<<endl;
      exit(1);
   }
   fastNLO::EMerge moption = _wgtoptions[wgtoption];
   fastNLO::EMerge prewgt  = _wgtoptions[preoptin];// alrady checked



   // --------------------------------------------------------------------------------------- //
   // --- calculate statistics for mergeing/normalisations
   // --------------------------------------------------------------------------------------- //
   // DB: Currently nothing is happening in the following code block
   //     Updates for option '-1' are needed as well...
   //map<string,fastNLOTable*> alltables;
   vector<fastNLOTable*> alltables;
   //vector<fastNLOTable*,string> allpaths;
   //map<string,unsigned int> lookup;
   //vector<fastNLO::WgtStat > allWgtStat(files.size());
   for ( auto path : files ) {
      info[_progname]<<"Reading table '" << path << "'" << endl;
      fastNLOTable* tab = new fastNLOTable(path);

      if ( tab->GetItabversion() < 23000 ) {
	 warn[_progname]<<"This program is maybe only compatible with fastNLO table versions > v2.3,000."<<endl;
	 warn[_progname]<<"Consider to use program fnlo-tk-merge and/or fnlo-tk-append instead."<<endl;
	 //exit(2);
      }

      //alltables[path] = tab;
      alltables.push_back(tab);
      //allpaths[tab] = path;
      // int tabid = alltables.size()-1;
      // lookup[path] = tabid;
      
      //! --- Check statistical information of additive contributions
      if ( options.count("-p") ) {
	 int nc = tab->GetNcontrib() + tab->GetNdata();
	 if ( nc != 1 ) {
	    warn[_progname]<<"Program fnlo-tk-merge2 currently can only handle fastNLO tables with exactly one contributions. Please use program fnlo-tk-merge fnlo-tk-append."<<endl;
	    //exit(2);
	 }
	 for ( int ic=0 ; ic<nc; ic++ ) {
	    bool quiet = true;
	    fastNLOCoeffBase* cnew = (fastNLOCoeffBase*)tab->GetCoeffTable(ic);
	    // Identify type of new coeff table, only additive ones have event numbers
	    if ( fastNLOCoeffAddBase::CheckCoeffConstants(cnew,quiet) ) { // additive?
	       fastNLOCoeffAddBase* cadd = (fastNLOCoeffAddBase*)cnew;
	       if ( cadd->GetNevt() == 1 ) {
		  warn[_progname]<<"Contribution #" << ic << " in table " << path << " has event number '1', which is usually invalid."<<endl;
		  // error[_progname]<<"Contribution #" << ic << " in table " << path << endl;
		  // error[_progname]<<"has no valid number-of-events information and cannot be merged. Aborted!" << endl;
		  // error[_progname]<<"Nevt = " << cadd->GetNevt() << endl;
		  //exit(1);
	       }
	       //
	       //allWgtStat[tabid] = cadd->GetWgtStat();
	    }
	    else {
	       error[_progname]<<"Program fnlo-tk-merge2 can only deal with additive contributions. Exiting"<<endl;
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
   // if ( options.count("-p") ) {
   //    vector<vector<vector<double> > > ProcBinWgt = CalculateWeight(allWgtStat,wgtoption);
   // }

   
   // --------------------------------------------------------------------------------------- //
   // --- calculating pre-averages if requested.
   if ( pre > 1 ) {
      fastNLOTable* keeptab = NULL;
      vector<fastNLOTable*> addtab;
      vector<fastNLOTable*> remainingtabs;
      for ( auto tab : alltables ) {
	 //fastNLOTable* tab = alltables[path]; // fastNLOTable tab(path);
	 if ( keeptab == NULL ) {
	    remainingtabs.push_back(tab);
	    keeptab = tab;
	    continue;
	 }
	 else {
	    addtab.push_back(tab);
	    //alltables.erase(path);
	 }
	 if ( int(addtab.size()) == pre-1 || alltables.back() == tab) { // merge
	    keeptab->MergeTables(addtab,prewgt);
	    addtab.clear();
	    keeptab = NULL;
	 }
      }
      alltables=remainingtabs;
   }

   // --------------------------------------------------------------------------------------- //
   //! --- Initialize output table
   fastNLOTable* resultTable = NULL;

   // --------------------------------------------------------------------------------------- //
   //! --- Loop input files and merge them
   // --------------------------------------------------------------------------------------- //
   int nValidTables = 0;
   vector<fastNLOTable*> allvec;
   for ( auto tab : alltables ) {
      // --- Initialising result with first read table
      if ( resultTable==NULL ) {
	 info[_progname]<<"Initializing result table '" << outfile << "'" << endl;
	 //resultTable = new fastNLOTable(*tab);
	 resultTable = tab;
	 nValidTables++;
	 continue;
      }
      else {
	 allvec.push_back(tab); // prepare for mergeing
	 // // --- Adding further tables to result table
	 // if ( moption==fastNLO::kMedian || moption == fastNLO::kMean )  
	 //    allvec.push_back(tab); // prepare for mergeing
	 // else 
	 //    resultTable->MergeTable(*tab, moption); // merge
	 nValidTables++;
      }
   }

   // ---- merge, if not yet done so
   // if ( moption == fastNLO::kMedian || moption == fastNLO::kMean ) {
      resultTable->MergeTables(allvec,moption);
   // }

   
   info[_progname]<<"Merged "<<nValidTables<<" table file(s)."<<endl;
   if (nValidTables == 0 ) {
      error[_progname]<<"Found less than two valid tables, no output possible!"<<endl;
      return 0;
   }

   // --- Write result
   resultTable->SetFilename(outfile);
   info[_progname]<<"Write merged results to file "<<resultTable->GetFilename()<<"."<<endl;
   resultTable->WriteTable();

   // --- fine
   return 0;
}
