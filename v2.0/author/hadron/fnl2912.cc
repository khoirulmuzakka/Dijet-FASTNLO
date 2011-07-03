//
// fastNLO v2 author code for fnl2912:
//     CMS LHC 3-jet Mass Scenario, E_cms = 7 TeV
//     for fastjet anti-kT algo with R=0.5 in E-scheme
//
// 
// ============== fastNLO user: ===================================
// To create your own scenario, it is recommended to take 
// this code, make a copy and edit the relevant changes.
// Important:
// Edit only those lines which are labeled as "fastNLO user"
// and refer to the documentation ("fastNLO authorcode in 
// NLOJET++") for a detailed explanation of the parameters 
// and variables.
// If a code fragment is not explicitely labeled as "fastNLO user",
// it is likely that a modification will interfere with
// the fastNLO routines.
//
// This file contains the following routines:
//   inputfunc    (-> user edits)
//   psinput      (-> user edits)
//   initfunc     (don't touch)
//   userfunc     (-> user edits)
//   writetable   (don't touch)
//   end_of_event (don't touch)
//   phys_output  (don't touch)
//   inittable    (-> user edits)
//
// Implementing a new scenario requires to edit:
//  - the jet algorithm ("#include" statement and assignment of "jetclus")
//  - number of jets (inputfunc)
//  - center-of-mass energy (psinput)
//  - compute observable, determine bin No. (userfunc)
//  - declare all variables for table, define bin boundaries (inittable)
//  
// ================================================================
// 
// last modifications
// 2011/06/29 KR Try 3-jet Mass Scenario, implement warm-up ext. vs. int. table, add debug mode
// 2011/06/17 KR Try CMS Dijet Mass Scenario
// 2011/01/13 KR unify jet sizes into one .cc and .h file
// 2010/09/24 MW make user-friendly
// 2010/09/22 MW implement D0 phase space
// 2009/01/15 TK make code V2.0 compatible
//

//------ DON'T TOUCH THIS PART! ------
#include <iostream>
#include <bits/hhc-phasespace.h>
#include <bits/hhc-process.h>
#include <bits/hhc-jetfunc.h>

//----- used namespaces -----
using namespace nlo;
using namespace std;

//----- declaration of the user defined functions -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
void psinput(phasespace_hhc *, double&);
user_base_hhc * userfunc();

//----- array of the symbols symbols -----
extern "C"{
  struct { 
    const char *name;
    void *address;
  } user_defined_functions[] = 
    {
      //   process index: hhc for hadron-hadron --> jets
      {"procindex", (void *) "hhc"},
      //   input function 
      {"inputfunc", (void *) inputfunc},
      //   phase space input function 
      {"psinput", (void *) psinput},
      //   user defined functions
      {"userfunc",  (void *) userfunc},
      //  end of the list
      {0, 0}
    };
}
//------ END OF THE DO-NOT-TOUCH-PART ------

//------ USER DEFINED PART STARTS HERE ------
#include <algorithm>

#include "fj-ak.h"   // fastNLO user: .h file for jet algorithm

#include "pdf-cteq6.h"
#include "pdf-hhc-dummy.h"
#include "fnloTable.h"
#include "fnloBlockBNlojet.h"

class UserHHC : public basic_user_set<user0d_hhc, user1h_hhc, user2h_hhc>
{
public:
  //   init and user function
  void initfunc(unsigned int);
  void userfunc(const event_hhc&, const amplitude_hhc&);
  virtual void end_of_event();  
  virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);

private:
  // --- pdf
  pdf_cteq6 pdf;
  pdf_hhc_dummy dummypdf;
  
  // --- jet algorithm
  fj_ak jetclus;   // fastNLO user: define jet algorithm (consistent with .h file above)
  
  bounded_vector<lorentzvector<double> > pj;    // the jet structure 
   
  // --- fastNLO definitions (not for user)
  fnloTable *table;
  double nevents;        // No of events calculated so far
  unsigned long nwrite;  // No of events after to write out the table
  string tablefilename;  // The table file to write to
  time_t start_time;
   
  // --- fastNLO user FYI: steering of these flags is encoded in NLOJet++ run name (option -n name)
  //     if filename matches "deb", "ref", or "wrm" the respective flag is set to true
  bool doDebug, doReference, doWarmUp;
  bool nlo;

  void inittable();
  void writetable();
};

user_base_hhc * userfunc() {
  return new UserHHC;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
  // --- fastNLO user: select the number of jets for your observable
  //nj = 1U;
  //nj = 2U;    
  nj = 3U;

  // --- number of the up and down type flavours (don't touch)
  nu = 2U;
  nd = 3U;
} 

void psinput(phasespace_hhc *ps, double& s)
{
  // --- fastNLO user: set the total c.m. energy squared in GeV^2
  //s =     40000.; // RHIC               200 GeV
  //s =   3240000.; // TeV Run I         1800 GeV
  //s =   3841600.; // TeV Run II        1960 GeV
  //s =    810000.; // LHC Injection Run  900 GeV
  //s =   5569600.; // LHC Initial Run   2360 GeV
  s =  49000000.; // LHC First Run     7000 GeV
  //s = 100000000.; // LHC Start-up Run 10000 GeV
  //s = 196000000.; // LHC Design Run   14000 GeV

  //   You can use your own phase generator. 
  //   Here we use the default.
  ps = 0;
} 

void UserHHC::initfunc(unsigned int)
{
  // --- Initialize event counters
  nevents = 0;
  // Set some defaults
  if (nwrite==0) nwrite = 5000000;
  start_time = std::time(0);
}

// --- fastNLO user: modify jet selection in userfunc (default = cutting in |y| min, |y| max and pt min)
//     (return value must be true for jets to be UNselected)
struct fNLOSelector {
  fNLOSelector(double ymin, double ymax, double ptmin):
    _ymin (ymin), _ymax (ymax), _ptmin (ptmin){};
  double _ymin, _ymax, _ptmin;
  bool operator() (const lorentzvector<double> &a) {return ! (_ymin <= abs(a.rapidity()) && abs(a.rapidity()) < _ymax && _ptmin <= a.perp());};
};

// --- fastNLO user: modify jet sorting in userfunc (default = descending in jet pt)
struct fNLOSorter {
  bool operator() (const lorentzvector<double> &a, const lorentzvector<double> &b) {return (a.perp() > b.perp());};
};

void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
  
  // --- fastNLO: Don't touch this piece of code!
  fnloBlockA2 *A2 =  table->GetBlockA2();
  double x1 = p[-1].Z()/p[hadron(-1)].Z();
  double x2 = p[0].Z()/p[hadron(0)].Z();
  
  // --- fastNLO user: Set jet size and run the jet algorithm
  double jetsize = 0.5;
  pj = jetclus(p,jetsize);
  unsigned int nj = pj.upper(); 

  // --- give some debug output before selection and sorting
  if ( doDebug ) {
    for (unsigned int i=1; i<=nj; i++) {
      double pti = pj[i].perp();
      double yi  = abs(pj[i].rapidity());
      cout << "before cuts: jet # i, pt, |y|: " << i << ", " << pti << ", " << yi << endl;
    }
  }
  
  // --- check on maximal no. of jets: 4 (should never be more in NLOJet++)
  if (nj > 4) {
    cout << "fastNLO: ERROR! This scenario is not suited for " << nj <<
      " jets. Aborted!" << endl;
    exit(1);
  }

  // --- fastNLO user:
  //     Here is your playground where you compute your observable 
  //     and the bin number ("obsbin") which gets passed to
  //     fastNLO's table filling code.
  //     (all pT and E are in GeV)
  
  // --- declare and initialize phase space cut variables
  // smallest |rapidity| for jets to be considered
  double yjmin  = 0.0;
  // largest |rapidity| for jets to be considered
  double yjmax  = 1.0;
  // lowest pT for jets to be considered
  double ptjmin = 40.;
  
  // --- select jets in y and ptjmin (failing jets are moved to the end of the jet array pj!)
  static fNLOSelector SelJets(yjmin,yjmax,ptjmin);
  // --- count number of selected jets left at this stage
  size_t njet = std::remove_if(pj.begin(), pj.end(), SelJets) - pj.begin();
  
  // --- sort selected n jets at beginning of jet array pj, by default decreasing in pt
  static fNLOSorter SortJets;
  std::sort(pj.begin(), pj.begin() + njet, SortJets);

  // --- give some debug output after selection and sorting
  if ( doDebug ) {
    cout << "# jets before and after phase space cuts: nj, njet = " << nj << ", " << njet << endl;
    cout << "phase space cuts: yjmin, yjmax, ptjmin: " << yjmin << ", " << yjmax << ", " << ptjmin << endl;
    for (unsigned int i=1; i<=njet; i++) {
      double pti = pj[i].perp();
      double yi  = abs(pj[i].rapidity());
      cout << "after cuts and sorting: jet # i, pt, |y|: " << i << ", " << pti << ", " << yi << endl;
    }
  }

  // --- declare and initialize additional cut variables
  // minimal dijet mass for events to be considered
  double m3jmin = 200.;
  // minimum pT fraction of leading jet required for 3rd jet
  double ptrelmin = 0.2;
  
  // 3-jet mass requires at least 3 jets
  if (njet > 2) {
    
    // Derive 3-jet variables
    // 3-jet mass
    lorentzvector<double> pj123 = pj[1] + pj[2] + pj[3]; 
    double m3j = pj123.mag();
    if (m3j < 0.) {cout << "Warning!: Negative Mass" << m3j << endl;}
    
    // Maximal (pseudo-)rapidity
    double y3jmax = max(abs(pj[1].rapidity()),abs(pj[2].rapidity()));
    y3jmax = max(y3jmax,abs(pj[3].rapidity()));

    // Ratio of pt3 to pt1 (ptrel)
    double ptrel = pj[3].perp()/pj[1].perp();
    
    // --- Further 3-jets phase space cuts?
    if ( m3jmin <= m3j && y3jmax < yjmax && ptrelmin <= ptrel ) {
      
      // --- set the renormalization and factorization scale to average 3-jet pT
      double mu = (pj[1].perp() + pj[2].perp() + pj[3].perp()) / 3.0;
      
      // --- identify bin number (dim1,dim2) here (m3j,|y3jmax|)
      int obsbin = -1;
      for (int j = 0; j < A2->GetNObsBin(); j++) {
	if (A2->LoBin[j][0] <= m3j    && m3j    < A2->UpBin[j][0] && 
	    A2->LoBin[j][1] <= y3jmax && y3jmax < A2->UpBin[j][1]) {
	  obsbin = j;
	  break;
	}
      }
      
      // --- fill fastNLO arrays - don't touch this piece of code!
      if (obsbin >= 0) {
	double prefactor = 1./A2->BinSize[obsbin]; // - divide by binwidth
	for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
	  if(table->GetBlockB(k)->GetIRef()>0){
	    ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,pdf,prefactor);
	  }else{
	    ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,dummypdf,prefactor);
	  }
	}
      } // --- end: fill fastNLO array
    } // --- end: final selection
  } // --- end: 3-jet events only
} // --- end: fastNLO user playground

void UserHHC::writetable(){
  table->OpenFileRewrite();
  table->WriteBlockA1();
  table->WriteBlockA2();
  for(int i=0;i< table->GetBlockA1()->GetNcontrib();i++){
    table->WriteBlockBDividebyN(i);
  }
  table->CloseFileWrite();

}

void UserHHC::end_of_event(){
  nevents += 1;
  // --- store table
  if (( (unsigned long)nevents % nwrite)==0){
    time_t hour, min, time = std::time(0) - start_time;
      
    hour = time/3600L;
    time -= hour*3600L;
    min  = time/60L;
    time -= min*60L;
      
    std::cout<<"--->     "
	     <<(hour < 10 ? "0" : "")<<hour
	     <<(min < 10 ? ":0" : ":")<<min
	     <<(time < 10 ? ":0" : ":")<<time<<std::endl;
    printf ("fastNLO: No. events: %.3G writing table ...\n",nevents);
    cout.flush();
    for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
      table->GetBlockB(k)->Nevt = (long long int)nevents;
    }
    writetable();
    printf("fastNLO: Table written.\n");
  }
}

void UserHHC::phys_output(const std::basic_string<char>& __file_name, 
                          unsigned long __save, bool __txt) 
{
  tablefilename.assign(__file_name.c_str());
  tablefilename += ".tab";
   
  // --- determine whether we are running LO or NLO
  const char* const file = __file_name.c_str(); 

  if(strstr(file,"born")!=NULL){
    nlo = false;
    printf("fastNLO: This is a LO run!\n");
  }else{
    if(strstr(file,"nlo")!=NULL){
      nlo = true;
      printf("fastNLO: This is a NLO run!\n");
    }else{
      printf("fastNLO: ERROR! This module can only be run at Born level or at NLO.\n");
      exit(1);
    }
  }

  // --- determine whether this is a debug, reference, or warm-up run
  doDebug = false;
  if (strstr(file,"deb")!=NULL) {
    doDebug = true;
    printf("fastNLO: This is a debug run. Attention, huge output\n");
  }
  doReference = false;
  if (strstr(file,"ref")!=NULL) {
    doReference = true;
    printf("fastNLO: This is a reference run!\n");
  }
  doWarmUp = false;
  if (strstr(file,"wrm")!=NULL) {
    doWarmUp = true;
    printf("fastNLO: This is a warm-up run!\n");
    if ( ! nlo ) {
      printf("fastNLO: WARNING! Warm-up runs are better done at NLO!\n");
    }
  }
  if ( doWarmUp && doReference ) {
    printf("fastNLO: ERROR! Warm-up and reference runs cannot be done simultaneously:\n");
    printf("         doWarmUp = %i, doReference = %i\n",doWarmUp,doReference);
    exit(2);
  }

  nwrite = __save;
  inittable();
}

void UserHHC::inittable(){

  // --- fastNLO user: set the total c.m. energy squared in GeV^2
  //double s =     40000.; // RHIC               200 GeV
  //double s =   3240000.; // TeV Run I         1800 GeV
  //double s =   3841600.; // TeV Run II        1960 GeV
  //double s =    810000.; // LHC Injection Run  900 GeV
  //double s =   5569600.; // LHC Initial Run   2360 GeV
  double s =  49000000.; // LHC First Run     7000 GeV
  //double s = 100000000.; // LHC Start-up Run 10000 GeV
  //double s = 196000000.; // LHC Design Run   14000 GeV

  // --- fastNLO user: decide whether to include a reference table (for 
  //                   precision studies, not for production jobs)
  // KR: Now set via filename string match to "ref"
  //const bool doReference = true;
  //const bool doReference = false;

  // --- set up fastNLO
  table = new fnloTable(tablefilename);

  // --- fastNLO: fill variable for table header block A1
  table->GetBlockA1()->SetScenName("fnl2912");  // - fastNLO user: set scenario name
  table->GetBlockA1()->SetNcontrib(1);
  table->GetBlockA1()->SetNmult(0);
  table->GetBlockA1()->SetNdata(0);
  // KR Add vars for Markus updated header
  table->GetBlockA1()->SetNuserString(0);
  table->GetBlockA1()->SetNuserInt(0);
  table->GetBlockA1()->SetNuserFloat(0);
  table->GetBlockA1()->SetImachine(0);
  // KR Ende
  table->GetBlockA2()->SetIpublunits(12);  // - fastNLO user: set cross section units
  //                 (negative power of ten)

  // --- fastNLO: fill variables for table header block A2
  fnloBlockA2 *A2 =  table->GetBlockA2();

  // --- fastNLO user: up to 20 strings to describe the scenario
  A2->ScDescript.push_back("d2sigma-3-jet_dM3Jd|y_max|_(pb_GeV)");
  A2->ScDescript.push_back("CMS_Collaboration");
  A2->ScDescript.push_back("3-jet_Mass");
  A2->ScDescript.push_back("anti-kT_R=0.5");
  A2->ScDescript.push_back("Test");
  
  A2->NScDescript = A2->ScDescript.size();
  A2->Ecms = sqrt(s);
  A2->ILOord = 3;   // --- fastNLO user: power of LO contr. for process (2 for incl. jets, 3 for 3-jet mass)
  A2->NDim = 2;     // --- fastNLO user: no. of dimensions in which observable is binned
  A2->DimLabel.push_back("M3J_[GeV]");  // --- fastNLO user: label of 1st dimension
  A2->IDiffBin.push_back(2);
  A2->DimLabel.push_back("|y_max|");   // --- fastNLO user: label of 2nd dimension
  A2->IDiffBin.push_back(2);

  vector <double> bound;
  bound.resize(2);

  // --- fastNLO user: bin definitions - here in |y_max| and m3j
  const int nrapbins = 1;
  double rapbins[nrapbins+1] = {0.0, 1.0};
  
  const int nmjjbins[nrapbins] = {10};
  
  cout << "------------------------" << endl;
  cout << "Binning in dimension 1: " << endl;
  for (int i=0; i<nrapbins+1; i++) {
    cout << "i, rapbins: " << i << ", " << rapbins[i] << endl;
  }
  
  vector< vector<double> >mjjbins;
  mjjbins.resize(nrapbins);
  for (int i=0; i<nrapbins; i++) {
    mjjbins[i].resize(nmjjbins[i]+1);
  }
  double mjj0[11] = { 200., 295., 348., 404., 529., 598., 751., 926., 1126., 1356., 1483. };
  for (int j=0; j<nmjjbins[0]+1; j++) { 
    mjjbins[0][j] = mjj0[j];
  }

  cout << "------------------------" << endl;
  cout << "Binning in dimension 2: " << endl;
  for (int i=0; i<nrapbins; i++) {
    for (int j=0; j<nmjjbins[i]+1; j++) {
      cout << "i, j, mjjbins: " << i << ", " << j << ", " << mjjbins[i][j] << endl;
    }
  }
  cout << "------------------------" << endl;

  // --- fastNLO user:
  //     define below the bin width ("binsize") by which
  //     the cross section is divided to obtain the 
  //     (multi-) differential result.
  //     default: divide by bin width in dim 1 and dim 2
  //              INCLUDING a factor of 2 for abs. rapidity |y| !
  
  int nbins = 0;   // --- count total No. bins
  for (int i=0;i<nrapbins;i++){
    for (int j=0;j<nmjjbins[i];j++){
      double binsize = 1.; // --- start each bin with preset value = 1.
      nbins += 1;
      bound[0] = mjjbins[i][j];
      bound[1] = rapbins[i];
      A2->LoBin.push_back(bound);
      bound[0] = mjjbins[i][j+1];
      bound[1] = rapbins[i+1];
      A2->UpBin.push_back(bound);
      binsize = binsize 
	* (mjjbins[i][j+1]-mjjbins[i][j]) // ... e.g. times dpT
	* 2. * (rapbins[i+1]-rapbins[i]); // ... e.g. times d|y|
      A2->BinSize.push_back(binsize);
    }
  }
  printf("fastNLO: Total no. of observable bins = %d\n",nbins+1);

  A2->NObsBin = nbins;
  A2->INormFlag = 0;   // --- fastNLO user: default=0 - set =1 if observable is 
  //     to be normalized by own integral (in 1st dimension)
  //     see documentation for details and for other options

  // --- fastNLO table block B
  fnloBlockB *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
  table->CreateBlockB(0,B);
  B->IXsectUnits = 12;    // --- fastNLO user: set to same value as "SetIpublunits"
  B->IDataFlag = 0;
  B->IAddMultFlag = 0;
  B->IContrFlag1 = 1;
  B->IContrFlag3 = 0;
  B->CodeDescript.push_back("NLOJet++_4.1.3");  // --- fastNLO user: enter NLOJET++ version
  B->NCodeDescr = B->CodeDescript.size();
  B->IRef = 0;
  if (nlo || A2->ILOord > 2) {
    B->NSubproc = 7;
  } else {
    B->NSubproc = 6;
  }
  printf("         This job uses %d subprocesses!\n",B->NSubproc);
  if(nlo){
    B->CtrbDescript.push_back("NLO");
    B->IContrFlag2 = 2;
    B->IScaleDep = 1;
    B->Npow = A2->ILOord+1;
  }else{
    B->CtrbDescript.push_back("LO");      
    B->IContrFlag2 = 1;
    B->IScaleDep = 0;
    B->Npow = A2->ILOord;
  }
  B->NContrDescr = B->CtrbDescript.size();

  B->NPDF = 2;
  B->NPDFPDG.push_back(2212);   // --- fastNLO user: PDG code for 1st hadron
  B->NPDFPDG.push_back(2212);  // --- fastNLO user: PDG code for 2nd hadron
  B->NPDFDim = 1;
  B->NFragFunc = 0;
  B->NFFDim = 0;
  B->IPDFdef1 = 3;
  B->IPDFdef2 = 1;
  if(B->NSubproc == 7) {
    B->IPDFdef3 = 2;
  } else {
    B->IPDFdef3 = 1;
  }
  printf("         Set IPDFdef3 = %d, consistent with %d subprocesses.\n",B->IPDFdef3,B->NSubproc);
  B->XNode1.resize(A2->NObsBin);

  // --- initialize variables for WarmUp run
  // KR: Now set via filename string match to "wrm"
  //B->IWarmUp = 1;     // --- fastNLO user: do the Warm-Up run
  B->IWarmUp = 0;     //                   or do production run(s)
  if ( doWarmUp ) {B->IWarmUp = 1;} 
  // KR: This is caught in an error condition now 
  // - fastNLO user: remember to disable reference-mode in
  //                 Warm-Up run: "doReference = false" (above)
  B->IWarmUpPrint = 100000;
  B->xlo.resize(A2->NObsBin);
  B->scalelo.resize(A2->NObsBin);
  B->scalehi.resize(A2->NObsBin);

  // --- arrays for extreme x and (default) scale values (computed in Warm-Up run)
  double xlim[A2->NObsBin];
  double mulo[A2->NObsBin];
  double muup[A2->NObsBin];

  // --- fastNLO user: before running a new scenario for the first time,
  //     the following block (between "start" and "end") should be 
  //     completely removed. The first run must be a "Warm-Up Run" (IWarmUp=1)
  //     which produces an initialization block as output. This should be copied
  //     and pasted below. These initialization values must be used for all
  //     production jobs (IWarmUp=0) for a given scenario (otherwise the result 
  //     tables can not be merged). 
  //
  if ( doWarmUp ) {
    cout << endl << "fastNLO: Initializing x limits for warm-up run ..." << endl;
    for(int i=0;i<A2->NObsBin;i++){
      xlim[i] = 1.1e-07, mulo[i] = 3.0, muup[i]=9.9e10; // - safe initializations
    }
  } else {
    cout << endl << "fastNLO: Initializing x limits ..." << endl;
    FILE * infile;
    infile = fopen("fastNLO-warmup.dat","r");
    if ( ! infile ) {
      cout << "fastNLO: WARNING! Could not read x limits from file: fastNLO-warmup.dat" << endl;
      cout << "         Trying to find and use x limits included in scenario author code ..." << endl;
      // --------- fastNLO: Warm-Up run results (start)
      // fastNLO user: paste warm-up run results here ...
      // 883400000 contributions (!= events) in warm-up run
      xlim[ 0 ] = 1.063511e-02 , mulo[ 0 ] = 4.472458e+01 , muup[ 0 ] = 3.534750e+02 ;
      xlim[ 1 ] = 1.567526e-02 , mulo[ 1 ] = 6.518156e+01 , muup[ 1 ] = 4.399650e+02 ;
      xlim[ 2 ] = 1.849052e-02 , mulo[ 2 ] = 7.608355e+01 , muup[ 2 ] = 4.996462e+02 ;
      xlim[ 3 ] = 2.150961e-02 , mulo[ 3 ] = 8.834811e+01 , muup[ 3 ] = 6.086999e+02 ;
      xlim[ 4 ] = 2.814612e-02 , mulo[ 4 ] = 1.160387e+02 , muup[ 4 ] = 6.392365e+02 ;
      xlim[ 5 ] = 3.178483e-02 , mulo[ 5 ] = 1.312755e+02 , muup[ 5 ] = 7.066228e+02 ;
      xlim[ 6 ] = 3.993076e-02 , mulo[ 6 ] = 1.642966e+02 , muup[ 6 ] = 7.494784e+02 ;
      xlim[ 7 ] = 4.917051e-02 , mulo[ 7 ] = 2.024171e+02 , muup[ 7 ] = 7.763640e+02 ;
      xlim[ 8 ] = 5.949442e-02 , mulo[ 8 ] = 2.459496e+02 , muup[ 8 ] = 8.187805e+02 ;
      xlim[ 9 ] = 7.195217e-02 , mulo[ 9 ] = 2.964002e+02 , muup[ 9 ] = 8.354523e+02 ;
      // --------- fastNLO: Warm-Up run results (end)
      // Safety check:
      // Count no. of xlim values > 0 (should be equal to NObsBin!)   
      int nxlim = 0;
      for (int i=0; i<A2->NObsBin; i++) {
	if (xlim[i] > 0.) {nxlim++;}
      }
      if (nxlim != A2->NObsBin ) {
	cerr << "fastNLO: ERROR! Could not find proper x limits in scenario code!" << endl;
	cerr << "         Do a warm-up run first." << endl;
	exit(1);
      }
    } else {
      cout << "fastNLO: Reading x limits from file fastNLO-warmup.dat ..." << endl;
      char line[256];
      // Ignore first documentation line
      if ( ! fgets(line,sizeof(line),infile) ) {
	cerr << "fastNLO: ERROR! Reading empty file: fastNLO-warmup.dat" << endl;
	exit(1);
      }
      printf("fastNLO: Read first, documentation line from file: fastNLO-warmup.dat\n");
      printf("%s",line);
      // Now read all limits
      int i = 0;
      while ( fgets(line,sizeof(line),infile) ) {
	int nvar = sscanf(line,"      xlim[ %*u ] = %le , mulo[ %*u ] = %le , muup[ %*u ] = %le ;",
			  &xlim[i],&mulo[i],&muup[i]);
	if (nvar != 3) {
	  cerr << "fastNLO: ERROR! x limits line did not match expected format:" << endl;
	  printf("%s",line);
	  exit(1);
	}
	i++;
      }
      if (i != A2->NObsBin ) {
	cerr << "fastNLO: ERROR! Number of x limits read != NObsBin: i = " << i << ", NObsBin = " << A2->NObsBin << endl;
	exit(1);
      }
    }
  }
  // --- print initialized values 
  cout << endl << "fastNLO: Print initialized x and mu limits:" << endl;
  for (int i=0; i<A2->NObsBin; i++) {
    printf("      xlim[ %u ] = %e , mulo[ %u ] = %e , muup[ %u ] = %e ;\n",
	   i,xlim[i],i,mulo[i],i,muup[i]);
  }
  
  for(int i=0;i<A2->NObsBin;i++){
    int nxtot = 15;
    if (i == ((A2->NObsBin)-1)) nxtot += 1; // Darf's etwas mehr sein?
    B->Nxtot1.push_back(nxtot);
    double hxlim = -sqrt(-log10(xlim[i]));   // use value from Warm-Up run
    //printf("%d %g %g \n",i,pow(10,-pow(hxlim,2)),xlim[i]);
    B->Hxlim1.push_back(hxlim);
    for(int j=0;j<nxtot;j++){
      double hx = hxlim*( 1.- ((double)j)/(double)nxtot);
      B->XNode1[i].push_back(pow(10,-pow(hx,2))); 
    }
  }

  B->NScales = 2;  // two scales: mur and muf
  B->NScaleDim = 1; // one variable used in scales: 3-jet pT average
  B->Iscale.push_back(0);  // mur=mur(pT), pT = index 0 
  B->Iscale.push_back(0);  // muf=muf(pT), pT = index 0 
  B->ScaleDescript.resize(B->NScaleDim);

  B->ScaleDescript[0].push_back("<pT_1,2,3>"); // --- fastNLO user: give name for selected scale
  B->NscaleDescript.push_back(B->ScaleDescript[0].size());
  //B->Nscalenode.push_back(4); // number of scale nodes for pT
  B->Nscalenode.push_back(6); // number of scale nodes for pT

  B->ScaleFac.resize(B->NScaleDim);

  B->ScaleFac[0].push_back(1.0);    // --- fastNLO: central scale (don't change)
  if(nlo){
    B->ScaleFac[0].push_back(2.0);  // --- fastNLO user: add any number of
    B->ScaleFac[0].push_back(0.5);  //             additional scale variations
    B->ScaleFac[0].push_back(0.25); //             as desired.
    //B->ScaleFac[0].push_back(4.0);
    //B->ScaleFac[0].push_back(8.0);
  }
  B->Nscalevar.push_back(B->ScaleFac[0].size());

  const double mu0scale = .25; // In GeV
  B->ScaleNode.resize(A2->NObsBin);
  B->HScaleNode.resize(A2->NObsBin);
  for(int i=0;i<A2->NObsBin;i++){
    B->ScaleNode[i].resize(B->NScaleDim);
    B->HScaleNode[i].resize(B->NScaleDim);
    for(int j=0;j<B->NScaleDim;j++){
      B->ScaleNode[i][j].resize(B->Nscalevar[j]);
      B->HScaleNode[i][j].resize(B->Nscalevar[j]);
      for(int k=0;k<B->Nscalevar[j];k++){
	B->ScaleNode[i][j][k].resize(B->Nscalenode[j]);
	B->HScaleNode[i][j][k].resize(B->Nscalenode[j]);
	if(B->Nscalenode[j]==1){
	  B->ScaleNode[i][j][k][0] = B->ScaleFac[0][k] * (muup[i]+mulo[i])/2.; // assume only one scale dimension
	  B->HScaleNode[i][j][k][0] = log(log((B->ScaleFac[0][k]*(muup[i]+mulo[i])/2.)/mu0scale));
	}else{
	  double llscalelo = log(log((B->ScaleFac[0][k]*mulo[i])/mu0scale));  
	  double llscalehi = log(log((B->ScaleFac[0][k]*muup[i])/mu0scale));  
	  for(int l=0;l<B->Nscalenode[j];l++){
	    B->HScaleNode[i][j][k][l] = llscalelo +
	      double(l)/double(B->Nscalenode[j]-1)*(llscalehi-llscalelo);
	    B->ScaleNode[i][j][k][l] = mu0scale * exp(exp(B->HScaleNode[i][j][k][l]));
	  }
	}
      }            
    }
  }

  B->SigmaTilde.resize(A2->NObsBin);
  for(int i=0;i<A2->NObsBin;i++){
    B->SigmaTilde[i].resize(B->Nscalevar[0]);
    for(int k=0;k<B->Nscalevar[0];k++){
      B->SigmaTilde[i][k].resize(B->Nscalenode[0]);
      for(int l=0;l<B->Nscalenode[0];l++){
	int nxmax = B->GetNxmax(i);
	B->SigmaTilde[i][k][l].resize(nxmax);
	for(int m=0;m<nxmax;m++){
	  B->SigmaTilde[i][k][l][m].resize(B->NSubproc);
	  for(int n=0;n<B->NSubproc;n++){
	    B->SigmaTilde[i][k][l][m][n] = 0.;
	  }
	}
      }            
    }
  }   
   
  // --- reference table
  if(doReference){
    fnloBlockB *refB = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
    if (nlo || A2->ILOord > 2) {
      refB->NSubproc = 7;
    } else {
      refB->NSubproc = 6;
    }
    printf("         This reference job uses %d subprocesses!\n",B->NSubproc);
    table->CreateBlockB(1,refB);
    refB->Copy(table->GetBlockB(0));
    refB->IRef = 1;
    refB->Nscalenode[0] = 1;
    refB->Nxtot1.clear();
    refB->Hxlim1.clear();
    for(int i=0;i<A2->NObsBin;i++){
      refB->Nxtot1.push_back(1);
      refB->Hxlim1.push_back(0.);
      refB->XNode1[i].clear(); 
      refB->XNode1[i].push_back(0.); 
    }
    table->GetBlockA1()->SetNcontrib(2);
  }

}
