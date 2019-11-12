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
#include <memory>
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
   {"-w","Weight    -w <option>. Calculate (un)weighted average. Option: GenWgt (default), unweighted, median, mean, NumEvt, SumW2, SumSig2, SumSig2BinProc, NumEvtBinProc, SumW2BinProc or User [file-of-weights]."},
   {"-attach","Do not 'merge' tables, but attach one to the other, i.e. result is the sum of all tables."},
   {"-add",   "Do not 'merge' tables, but  add   one to the other, i.e. result is the sum of all tables (formerly called 'append')."},
   {"-pre","pre-avg -pre <n> <option>. 2 step mergeing: Build pre-averaged of n tables using weighting procedure <option>."},
   {"-cutRMS","cur RMS -cutRMS <sigma>. Calculate for each node the mean and RMS of all tables and discard values greater than <sigma>*RMS during mergeing."},
};


std::map<std::string,fastNLO::EMerge> _wgtoptions {
   {"GenWgt",fastNLO::kMerge },
   {"add",fastNLO::kAdd },
   {"append",fastNLO::kAdd }, // deprecated
   {"attach",fastNLO::kAttach }, // deprecated
   {"unweighted",fastNLO::kUnweighted },
   {"median",fastNLO::kMedian },
   {"mean",fastNLO::kMean },
   {"NumEvt",fastNLO::kNumEvent },
   {"SumW2",fastNLO::kSumW2 },
   {"SumSig2",fastNLO::kSumSig2 },
   {"NumEvtBinProc",fastNLO::kNumEventBinProc },
   {"SumW2BinProc",fastNLO::kSumW2BinProc },
   {"SumSig2BinProc",fastNLO::kSumSig2BinProc },
   {"NNLOJET",fastNLO::kSumUserBinProc }, // user weigths are stored in 'SumSig'
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
   man << "Warning: Some mergeing options (e.g. 'median') may require lots of RAM as all numbers must be kept in memory at a time."<<endl;
   man <<"          In these cases it may become benefitial to use an (unbiased) pre-averageing, e.g. with weighting options mean, NumEvt or NumEvtBinProc"<<endl;
   man << " " <<endl;
}

//__________________________________________________________________________________________________________________________________
std::map<std::string,std::vector<double> > ReadNnlojetWgtFile(std::string wgtFile,std::vector<std::string> files);


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

   // --------------------------------------------------------------------------------------- //
   //! --- Parse commmand line
   if (argc <= 1) {
      info[_progname] << "For more explanations type:" << endl;
      info[_progname] << "  $  fnlo-tk-merge2 -h" << endl<<endl;
      error[_progname] << "No arguments given, but need at least two!" << endl;
      PrintHelpMessage();   exit(1);
   }
   set<string> setfiles;
   vector<string> files;
   set<string> options;
   string outfile = argv[argc-1];
   string wgtoption = "GenWgt";
   string wgtFile;
   string plotfile = "fnlo-tk-merge2.ps";
   int pre = 0;
   string preoptin = "NumEvtBinProc"; // NumEvt
   double cutRMS = 0;

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
	 if ( sarg == "-w" ) { 
	    wgtoption=argv[++iarg];
	    if ( wgtoption == "NNLOJET" ) {	       
	       wgtFile=argv[++iarg];
	       info[_progname]<<"NNLOJET specified weights are read from file: "<<wgtFile<<endl;
	    }
	 }
	 if ( sarg == "-o" ) { outfile=argv[++iarg]; narg++; }
	 if ( sarg == "-cutRMS" ) cutRMS=atof(argv[++iarg]); 
	 if ( sarg == "-add" ) { wgtoption = "add"; }
	 if ( sarg == "-attach" ) { wgtoption = "attach"; }
	 if ( sarg == "-append" ) { wgtoption = "add"; }
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

   // --- option 'NNLOJET'
   map<string,vector<double> > mUserWgts; // [filename] [weight-per-obsbin]
   if ( moption==fastNLO::kSumUserBinProc ) {
      if ( pre > 0 ) {
	 error[_progname]<<"User specified weights cannot be combined with option -pre."<<endl;
	 exit(1);
      }
      // read file
      mUserWgts = ReadNnlojetWgtFile(wgtFile,files);
   }



   // --------------------------------------------------------------------------------------- //
   // --- calculate statistics for mergeing/normalisations
   // --------------------------------------------------------------------------------------- //
   // DB: Currently nothing is happening in the following code block
   //     Updates for option '-1' are needed as well...
   //map<string,fastNLOTable*> alltables;
   //vector<fastNLOTable*,string> allpaths;
   //map<string,unsigned int> lookup;
   //vector<fastNLO::WgtStat > allWgtStat(files.size());
   /*
   if ( options.count("-1") == 0 ) { // read in all files now!
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
	 }
      }
   }   
   */

   // --------------------------------------------------------------------------------------- //
   // --- calculate statistics
   // if ( options.count("-p") ) {
   //    vector<vector<vector<double> > > ProcBinWgt = CalculateWeight(allWgtStat,wgtoption);
   // }

   
   // --------------------------------------------------------------------------------------- //
   // --- calculating pre-averages if requested.
   vector<fastNLOTable*> alltables;
   if ( pre > 1 ) {
      fastNLOTable* keeptab = NULL;
      vector<fastNLOTable*> addtab;
      //for ( auto tab : alltables ) {
      for ( const string& path : files ) {
	 //fastNLOTable* tab = alltables[path]; // fastNLOTable tab(path);
	 fastNLOTable* tab = new fastNLOTable(path);
	 if ( keeptab == NULL ) {
	    alltables.push_back(tab);
	    keeptab = tab;
	    continue;
	 }
	 else {
	    addtab.push_back(tab);
	 }
	 if ( int(addtab.size()) == pre-1 /*|| files.back() == path*/ ) { // merge
	    keeptab->MergeTables(addtab,prewgt);
	    for ( fastNLOTable* t : addtab ) delete t;
	    addtab.clear();
	    keeptab = NULL;
	 }
      }
      if ( addtab.size() ) {
	 keeptab->MergeTables(addtab,prewgt);
	 for ( fastNLOTable* t : addtab ) delete t;
      }
   }

   // --------------------------------------------------------------------------------------- //
   // --- Loop input files and merge them
   // --------------------------------------------------------------------------------------- //
   fastNLOTable* resultTable = NULL;
   int nValidTables = 1;
   if ( alltables.empty() ) {
      for ( const string& path : files ) {
	 if ( resultTable==NULL ) {
	    resultTable = new fastNLOTable(path);
	    if ( moption==fastNLO::kSumUserBinProc ) 
	       resultTable->SetUserWeights(mUserWgts[path]);
	 }
	 else {
	    if ( moption == fastNLO::kMedian || moption == fastNLO::kMean || cutRMS!=0 ) 
	       alltables.push_back(new fastNLOTable(path));
	    else {
	       //resultTable->MergeTable(fastNLOTable(path), moption); // merge
	       fastNLOTable tab(path);
	       if ( moption==fastNLO::kSumUserBinProc ) tab.SetUserWeights(mUserWgts[path]);
	       resultTable->MergeTable(tab, moption); // merge
	    }
	    nValidTables++;
	 }
      }
   }
   else {
      nValidTables = alltables.size();
      resultTable = alltables.back();
      alltables.pop_back();
   }
   resultTable->MergeTables(alltables,moption,cutRMS);


   // --------------------------------------------------------------------------------------- //
   // --- Some concluding remarks
   // --------------------------------------------------------------------------------------- //
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



//__________________________________________________________________________________________________________________________________
std::map<std::string,std::vector<double> > ReadNnlojetWgtFile( std::string wgtFile, std::vector<std::string> files) {
   //!< Read ascii file with weigts from NNLOJET combination script
   using namespace std;

   // --- open wgt file
   std::ifstream is (wgtFile.c_str());
   if (!is) {
      std::cout<<"[ReadNnlojetWgtFile] ERROR. Cannot open file!"<<std::endl;
      exit(3);
   }
   cout<<endl;
   cout<<" ---- Reading file with weights ---- "<<endl;

   // --- read wgt file
   std::string line;
   std::string item;
   vector<vector<string> > content;
   while (getline(is, line) ) {
      //cout<<line<<endl;
      std::stringstream ss(line);
      content.push_back(vector<string>());
      while (std::getline(ss, item, ' ')) {
	 //cout<<item<<endl;
	 content.back().push_back(item);
      }
      //for ( auto pp : tokens ) cout<<pp<<endl;
   }
   std::cout<<"Read weights for "<<content.size()-1<< " output files and "<<content[0].size()-1<<" bins."<<endl;

   // --- fill return map
   std::map<std::string,std::vector<double> > wgts;
   set<string> uniquefiles;
   for ( auto ff : files ) {
      if ( uniquefiles.count(ff) ) {
	 cout<<"Info. Duplicate input file detected."<<endl;
	 continue;
      }
      uniquefiles.insert(ff);
      std::string ftag = ff;
      if ( ftag.find(".dat") != string::npos) ftag.resize(ftag.find(".dat"));
      if ( ftag.find(".tab") != string::npos) ftag.resize(ftag.find(".tab"));
      if ( ftag.find("/") != string::npos) ftag=ftag.substr(ftag.find_last_of("/")+1,ftag.size());
      ftag += ".";
      //cout<<"ftag="<<ftag<<endl;
      for ( vector<string>& lit : content ) {
	 string nnfile = lit[0];
	 if ( nnfile.find(ftag) != string::npos ){
	    // found wgts for input table ff
	    if ( wgts.count(ff) ) {
	       cout<<"Warning. Weights already found for input table"
		   <<ff<<" (was searching for '"<<ftag<<"'). Maybe duplicates for input?"<<endl; 
	       exit(4);
	    };
	    for ( string cc : lit ) {
	       if ( wgts.count(ff) == 0 ) wgts[ff]; // first entry is the filename, instantiate new vector
	       else {
		  double dval = strtod(cc.c_str(),NULL);
		  if (dval==0 ) {
		     cout<<"Warning. Weight is zero, using tiny number instead. (input: "
			 <<cc<<") file: "<<ff<<endl;
		     dval=1.e-20;
		  }
		  wgts[ff].push_back(dval);
	       }
	    }
	    std::cout<<"Found "<<wgts[ff].size()<<" weights for input file: "<<ff<<std::endl;
	    if ( wgts[ff].size() != content[0].size()-1 ) {
	       std::cout<<"ERROR. Too little weights for file "<<ff<<endl;
	       exit(3);
	    }
	 }
      }
      if ( wgts.count(ff) == 0) {
	 std::cout<<"Error. No weights for table "<<ff<<" found in wight-file "<<wgtFile<<endl;
	 exit(3);
      }
   }
   // --- weights for all input fastNLO files found?
   if ( wgts.size()!= uniquefiles.size() ){
      std::cout<<"ERROR. Could not find weights for all input fastNLO files. Exiting."<<endl;
      exit(3);
   }

   return wgts;
}
