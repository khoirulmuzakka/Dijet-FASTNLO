//
// fastNLO v2 author code for fnl2622e:
//     CMS LHC Dijet Chi Scenario, E_cms = 7 TeV
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
  nj = 2U;    
  //nj = 3U;

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
  const double yjmin  = 0.0;
  // largest |rapidity| for jets to be considered
  const double yjmax  = 5.0;
  // lowest pT for jets to be considered
  const double ptjmin = 0.;
  
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

  // Dijet chi requires at least 2 jets
  if (njet > 1) {

    // --- declare and initialize additional cut variables
    // minimal dijet mass for events to be considered
    const double mjjmin = 400.;
    // maximal rapidity of the two leading jets
    const double yjjmax = 2.5;
    // maximal allowed y_boost = 0.5 * |y1 + y2|
    const double yboostmax = 1.11;
    // maximal allowed y_star = 0.5 * |y1 - y2|
    const double ystarmax  = 1.5;
    
    // Derive dijet variables
    // Dijet mass
    lorentzvector<double> pj12 = pj[1] + pj[2]; 
    double mjj = pj12.mag();
    if (mjj < 0.) {cout << "Warning!: Negative Mass" << mjj << endl;}

    // Rapidities of two leading jets
    double y1 = pj[1].rapidity();
    double y2 = pj[2].rapidity();

    // Sum, diff and chi for these rapidities
    double syjj  = abs(y1+y2);
    double dyjj  = abs(y1-y2);
    double chijj = exp(dyjj);

    // --- Further dijet phase space cuts?
    if (syjj < 2.*yboostmax && dyjj < 2.*ystarmax &&
	abs(y1) < yjjmax && abs(y2) < yjjmax &&  
	mjj >= mjjmin ) {
      
      // --- set the renormalization and factorization scale to average dijet pT
      double mu = (pj[1].perp() + pj[2].perp()) / 2.0;
      
      // --- identify bin number (dim1,dim2) here (chijj,mjj)
      int obsbin = -1;
      for (int j = 0; j < A2->GetNObsBin(); j++) {
	if (A2->LoBin[j][0] <= chijj && chijj < A2->UpBin[j][0] && 
	    A2->LoBin[j][1] <= mjj   && mjj   < A2->UpBin[j][1]) {
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
  } // --- end: dijet+ events only
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
    printf("fastNLO: This is a debug run. Attention, huge output!\n");
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
  table->GetBlockA1()->SetScenName("fnl2622e");  // - fastNLO user: set scenario name
  table->GetBlockA1()->SetNcontrib(1);
  table->GetBlockA1()->SetNmult(0);
  table->GetBlockA1()->SetNdata(0);
  // KR Add vars for Markus updated header
  table->GetBlockA1()->SetNuserString(0);
  table->GetBlockA1()->SetNuserInt(0);
  table->GetBlockA1()->SetNuserFloat(0);
  table->GetBlockA1()->SetImachine(0);
  // KR Ende
  table->GetBlockA2()->SetIpublunits(9);  // - fastNLO user: set cross section units
  //                 (negative power of ten)

  // --- fastNLO: fill variables for table header block A2
  fnloBlockA2 *A2 =  table->GetBlockA2();

  // --- fastNLO user: up to 20 strings to describe the scenario
  A2->ScDescript.push_back("d2sigma-dijet_dChidMJJ_(pb_GeV)");
  A2->ScDescript.push_back("CMS_Collaboration");
  A2->ScDescript.push_back("Dijet_Chi");
  A2->ScDescript.push_back("anti-kT_R=0.5");
  A2->ScDescript.push_back("Test");
  
  A2->NScDescript = A2->ScDescript.size();
  A2->Ecms = sqrt(s);
  A2->ILOord = 2;   // --- fastNLO user: power of LO contr. for process (2 for incl. jets, 3 for 3-jet mass)
  A2->NDim = 2;     // --- fastNLO user: no. of dimensions in which observable is binned
  A2->DimLabel.push_back("Chi");  // --- fastNLO user: label of 1st dimension
  A2->IDiffBin.push_back(2);
  A2->DimLabel.push_back("MJJ_[GeV]");   // --- fastNLO user: label of 2nd dimension
  A2->IDiffBin.push_back(2);

  vector <double> bound;
  bound.resize(2);

  // --- fastNLO user: bin definitions - here in chijj and mjj
  const int ndim2bins = 13;
  const double dim2bins[ndim2bins+1] = {400., 600., 900., 1200., 1500.,
					1900., 2300., 2400., 2800., 3000.,
					3200., 4000., 5000., 7000.};
  
  const int ndim1bins[ndim2bins] = {12,12,12,12,12,12,12,12,12,12,12,12,12};
  
  cout << endl << "------------------------" << endl;
  cout << "Binning in dimension 2: " << A2->DimLabel[1] << endl;
  cout << "------------------------" << endl;
  for (int i=0; i<ndim2bins+1; i++) {
    cout << "i, dim2bins: " << i << ", " << dim2bins[i] << endl;
  }
  
  vector< vector<double> >dim1bins;
  dim1bins.resize(ndim2bins);
  for (int i=0; i<ndim2bins; i++) {
    dim1bins[i].resize(ndim1bins[i]+1);
  }
  const double dim1[13] = { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.,
			    12., 14., 16. };
  for (int i=0; i<ndim2bins; i++) {
    for (int j=0; j<ndim1bins[i]+1; j++) { 
      dim1bins[i][j] = dim1[j];
    }
  }
  
  cout << endl << "------------------------" << endl;
  cout << "Binning in dimension 1: " << A2->DimLabel[0] << endl;
  cout << "------------------------" << endl;
  for (int i=0; i<ndim2bins; i++) {
    for (int j=0; j<ndim1bins[i]+1; j++) {
      cout << "i, j, dim1bins: " << i << ", " << j << ", " << dim1bins[i][j] << endl;
    }
  }
  cout << "========================" << endl;
  
  // --- fastNLO user:
  //     define below the bin width ("binsize") by which
  //     the cross section is divided to obtain the 
  //     (multi-) differential result.
  //     default: divide by bin width in dim 1 and dim 2
  //              ATTENTION: Don't forget to include a factor of 2 for abs. rapidity |y| !
  
  int nbins = 0;   // --- count total No. bins
  for (int i=0;i<ndim2bins;i++){
    for (int j=0;j<ndim1bins[i];j++){
      double binsize = 1.; // --- start each bin with preset value = 1.
      nbins += 1;
      bound[0] = dim1bins[i][j];
      bound[1] = dim2bins[i];
      A2->LoBin.push_back(bound);
      bound[0] = dim1bins[i][j+1];
      bound[1] = dim2bins[i+1];
      A2->UpBin.push_back(bound);
      binsize = binsize 
	* (dim1bins[i][j+1]-dim1bins[i][j]) // ... e.g. times dpT
	* (dim2bins[i+1]-dim2bins[i]); // ... e.g. times d|y|
      A2->BinSize.push_back(binsize);
    }
  }
  printf("fastNLO: Total no. of observable bins = %d\n\n",nbins);

  A2->NObsBin = nbins;
  A2->INormFlag = 0;   // --- fastNLO user: default=0 - set =1 if observable is 
  //     to be normalized by own integral (in 1st dimension)
  //     see documentation for details and for other options

  // --- fastNLO table block B
  fnloBlockB *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
  table->CreateBlockB(0,B);
  B->IXsectUnits = 9;    // --- fastNLO user: set to same value as "SetIpublunits"
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
  printf("fastNLO: This job uses %d subprocesses!\n",B->NSubproc);
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
  printf("fastNLO: Set IPDFdef3 = %d, consistent with %d subprocesses.\n",B->IPDFdef3,B->NSubproc);
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
      // 3337300000 contributions (!= events) in warm-up run
      xlim[ 0 ] = 1.838033e-02 , mulo[ 0 ] = 1.862831e+02 , muup[ 0 ] = 3.447568e+02 ;
      xlim[ 1 ] = 1.760530e-02 , mulo[ 1 ] = 1.712016e+02 , muup[ 1 ] = 3.187850e+02 ;
      xlim[ 2 ] = 1.702450e-02 , mulo[ 2 ] = 1.581962e+02 , muup[ 2 ] = 2.849849e+02 ;
      xlim[ 3 ] = 1.636166e-02 , mulo[ 3 ] = 1.473751e+02 , muup[ 3 ] = 2.603587e+02 ;
      xlim[ 4 ] = 1.608440e-02 , mulo[ 4 ] = 1.380772e+02 , muup[ 4 ] = 2.392580e+02 ;
      xlim[ 5 ] = 1.602083e-02 , mulo[ 5 ] = 1.311352e+02 , muup[ 5 ] = 2.225461e+02 ;
      xlim[ 6 ] = 1.567624e-02 , mulo[ 6 ] = 1.243378e+02 , muup[ 6 ] = 2.096778e+02 ;
      xlim[ 7 ] = 1.561177e-02 , mulo[ 7 ] = 1.189285e+02 , muup[ 7 ] = 1.994550e+02 ;
      xlim[ 8 ] = 1.549138e-02 , mulo[ 8 ] = 1.136094e+02 , muup[ 8 ] = 1.900092e+02 ;
      xlim[ 9 ] = 1.523885e-02 , mulo[ 9 ] = 1.055191e+02 , muup[ 9 ] = 1.820832e+02 ;
      xlim[ 10 ] = 1.494064e-02 , mulo[ 10 ] = 9.888522e+01 , muup[ 10 ] = 1.692218e+02 ;
      xlim[ 11 ] = 1.493291e-02 , mulo[ 11 ] = 9.314830e+01 , muup[ 11 ] = 1.578008e+02 ;
      xlim[ 12 ] = 2.784589e-02 , mulo[ 12 ] = 2.792793e+02 , muup[ 12 ] = 5.181139e+02 ;
      xlim[ 13 ] = 2.687849e-02 , mulo[ 13 ] = 2.560397e+02 , muup[ 13 ] = 4.788359e+02 ;
      xlim[ 14 ] = 2.563577e-02 , mulo[ 14 ] = 2.372437e+02 , muup[ 14 ] = 4.306908e+02 ;
      xlim[ 15 ] = 2.501458e-02 , mulo[ 15 ] = 2.208059e+02 , muup[ 15 ] = 3.906922e+02 ;
      xlim[ 16 ] = 2.438827e-02 , mulo[ 16 ] = 2.075924e+02 , muup[ 16 ] = 3.596057e+02 ;
      xlim[ 17 ] = 2.381605e-02 , mulo[ 17 ] = 1.961222e+02 , muup[ 17 ] = 3.339229e+02 ;
      xlim[ 18 ] = 2.412516e-02 , mulo[ 18 ] = 1.866014e+02 , muup[ 18 ] = 3.152671e+02 ;
      xlim[ 19 ] = 2.348781e-02 , mulo[ 19 ] = 1.780562e+02 , muup[ 19 ] = 2.993105e+02 ;
      xlim[ 20 ] = 2.348008e-02 , mulo[ 20 ] = 1.706216e+02 , muup[ 20 ] = 2.849732e+02 ;
      xlim[ 21 ] = 2.301544e-02 , mulo[ 21 ] = 1.581492e+02 , muup[ 21 ] = 2.734097e+02 ;
      xlim[ 22 ] = 2.284484e-02 , mulo[ 22 ] = 1.479939e+02 , muup[ 22 ] = 2.527343e+02 ;
      xlim[ 23 ] = 2.224351e-02 , mulo[ 23 ] = 1.396973e+02 , muup[ 23 ] = 2.372201e+02 ;
      xlim[ 24 ] = 4.196160e-02 , mulo[ 24 ] = 4.188582e+02 , muup[ 24 ] = 6.902142e+02 ;
      xlim[ 25 ] = 4.059052e-02 , mulo[ 25 ] = 3.846306e+02 , muup[ 25 ] = 6.402576e+02 ;
      xlim[ 26 ] = 3.928574e-02 , mulo[ 26 ] = 3.556175e+02 , muup[ 26 ] = 5.728551e+02 ;
      xlim[ 27 ] = 3.872026e-02 , mulo[ 27 ] = 3.317293e+02 , muup[ 27 ] = 5.202949e+02 ;
      xlim[ 28 ] = 3.733201e-02 , mulo[ 28 ] = 3.115621e+02 , muup[ 28 ] = 4.792282e+02 ;
      xlim[ 29 ] = 3.715982e-02 , mulo[ 29 ] = 2.937143e+02 , muup[ 29 ] = 4.450311e+02 ;
      xlim[ 30 ] = 3.645947e-02 , mulo[ 30 ] = 2.796687e+02 , muup[ 30 ] = 4.189089e+02 ;
      xlim[ 31 ] = 3.662981e-02 , mulo[ 31 ] = 2.665700e+02 , muup[ 31 ] = 3.979759e+02 ;
      xlim[ 32 ] = 3.627672e-02 , mulo[ 32 ] = 2.561997e+02 , muup[ 32 ] = 3.808425e+02 ;
      xlim[ 33 ] = 3.526240e-02 , mulo[ 33 ] = 2.372256e+02 , muup[ 33 ] = 3.645720e+02 ;
      xlim[ 34 ] = 3.447388e-02 , mulo[ 34 ] = 2.216282e+02 , muup[ 34 ] = 3.370639e+02 ;
      xlim[ 35 ] = 3.461035e-02 , mulo[ 35 ] = 2.091767e+02 , muup[ 35 ] = 3.167088e+02 ;
      xlim[ 36 ] = 5.645768e-02 , mulo[ 36 ] = 5.582200e+02 , muup[ 36 ] = 8.616126e+02 ;
      xlim[ 37 ] = 5.516051e-02 , mulo[ 37 ] = 5.120625e+02 , muup[ 37 ] = 7.956605e+02 ;
      xlim[ 38 ] = 5.410381e-02 , mulo[ 38 ] = 4.744835e+02 , muup[ 38 ] = 7.180584e+02 ;
      xlim[ 39 ] = 5.298734e-02 , mulo[ 39 ] = 4.415656e+02 , muup[ 39 ] = 6.503359e+02 ;
      xlim[ 40 ] = 5.242960e-02 , mulo[ 40 ] = 4.143210e+02 , muup[ 40 ] = 6.002063e+02 ;
      xlim[ 41 ] = 5.108038e-02 , mulo[ 41 ] = 3.915452e+02 , muup[ 41 ] = 5.567739e+02 ;
      xlim[ 42 ] = 5.030743e-02 , mulo[ 42 ] = 3.723603e+02 , muup[ 42 ] = 5.251538e+02 ;
      xlim[ 43 ] = 4.921029e-02 , mulo[ 43 ] = 3.557029e+02 , muup[ 43 ] = 4.988581e+02 ;
      xlim[ 44 ] = 4.921671e-02 , mulo[ 44 ] = 3.413277e+02 , muup[ 44 ] = 4.757709e+02 ;
      xlim[ 45 ] = 4.842957e-02 , mulo[ 45 ] = 3.166815e+02 , muup[ 45 ] = 4.556294e+02 ;
      xlim[ 46 ] = 4.704039e-02 , mulo[ 46 ] = 2.957585e+02 , muup[ 46 ] = 4.224119e+02 ;
      xlim[ 47 ] = 4.757501e-02 , mulo[ 47 ] = 2.797115e+02 , muup[ 47 ] = 3.947379e+02 ;
      xlim[ 48 ] = 7.056452e-02 , mulo[ 48 ] = 6.981994e+02 , muup[ 48 ] = 1.089944e+03 ;
      xlim[ 49 ] = 6.948753e-02 , mulo[ 49 ] = 6.406503e+02 , muup[ 49 ] = 1.007560e+03 ;
      xlim[ 50 ] = 6.847015e-02 , mulo[ 50 ] = 5.928715e+02 , muup[ 50 ] = 9.060336e+02 ;
      xlim[ 51 ] = 6.750558e-02 , mulo[ 51 ] = 5.513220e+02 , muup[ 51 ] = 8.271063e+02 ;
      xlim[ 52 ] = 6.675446e-02 , mulo[ 52 ] = 5.179406e+02 , muup[ 52 ] = 7.592594e+02 ;
      xlim[ 53 ] = 6.646742e-02 , mulo[ 53 ] = 4.900802e+02 , muup[ 53 ] = 7.070293e+02 ;
      xlim[ 54 ] = 6.621235e-02 , mulo[ 54 ] = 4.648495e+02 , muup[ 54 ] = 6.642236e+02 ;
      xlim[ 55 ] = 6.531106e-02 , mulo[ 55 ] = 4.458514e+02 , muup[ 55 ] = 6.305833e+02 ;
      xlim[ 56 ] = 6.438555e-02 , mulo[ 56 ] = 4.274236e+02 , muup[ 56 ] = 6.027141e+02 ;
      xlim[ 57 ] = 6.425720e-02 , mulo[ 57 ] = 3.959310e+02 , muup[ 57 ] = 5.761456e+02 ;
      xlim[ 58 ] = 6.348101e-02 , mulo[ 58 ] = 3.706904e+02 , muup[ 58 ] = 5.354364e+02 ;
      xlim[ 59 ] = 6.192475e-02 , mulo[ 59 ] = 3.493920e+02 , muup[ 59 ] = 4.997832e+02 ;
      xlim[ 60 ] = 8.939671e-02 , mulo[ 60 ] = 8.852353e+02 , muup[ 60 ] = 1.323014e+03 ;
      xlim[ 61 ] = 8.913427e-02 , mulo[ 61 ] = 8.117494e+02 , muup[ 61 ] = 1.221255e+03 ;
      xlim[ 62 ] = 8.851473e-02 , mulo[ 62 ] = 7.512781e+02 , muup[ 62 ] = 1.096617e+03 ;
      xlim[ 63 ] = 8.886233e-02 , mulo[ 63 ] = 6.991158e+02 , muup[ 63 ] = 9.989677e+02 ;
      xlim[ 64 ] = 8.823463e-02 , mulo[ 64 ] = 6.557360e+02 , muup[ 64 ] = 9.181593e+02 ;
      xlim[ 65 ] = 8.778598e-02 , mulo[ 65 ] = 6.203939e+02 , muup[ 65 ] = 8.517534e+02 ;
      xlim[ 66 ] = 8.730206e-02 , mulo[ 66 ] = 5.903435e+02 , muup[ 66 ] = 8.048767e+02 ;
      xlim[ 67 ] = 8.682946e-02 , mulo[ 67 ] = 5.627329e+02 , muup[ 67 ] = 7.629807e+02 ;
      xlim[ 68 ] = 8.615306e-02 , mulo[ 68 ] = 5.390193e+02 , muup[ 68 ] = 7.297032e+02 ;
      xlim[ 69 ] = 8.561654e-02 , mulo[ 69 ] = 5.010628e+02 , muup[ 69 ] = 6.973662e+02 ;
      xlim[ 70 ] = 8.582495e-02 , mulo[ 70 ] = 4.683080e+02 , muup[ 70 ] = 6.465750e+02 ;
      xlim[ 71 ] = 8.659799e-02 , mulo[ 71 ] = 4.413293e+02 , muup[ 71 ] = 6.069323e+02 ;
      xlim[ 72 ] = 1.080593e-01 , mulo[ 72 ] = 1.069860e+03 , muup[ 72 ] = 1.380683e+03 ;
      xlim[ 73 ] = 1.080010e-01 , mulo[ 73 ] = 9.819989e+02 , muup[ 73 ] = 1.269847e+03 ;
      xlim[ 74 ] = 1.080887e-01 , mulo[ 74 ] = 9.081687e+02 , muup[ 74 ] = 1.141779e+03 ;
      xlim[ 75 ] = 1.080929e-01 , mulo[ 75 ] = 8.457281e+02 , muup[ 75 ] = 1.040627e+03 ;
      xlim[ 76 ] = 1.080099e-01 , mulo[ 76 ] = 7.943921e+02 , muup[ 76 ] = 9.593618e+02 ;
      xlim[ 77 ] = 1.080124e-01 , mulo[ 77 ] = 7.511458e+02 , muup[ 77 ] = 8.898864e+02 ;
      xlim[ 78 ] = 1.083083e-01 , mulo[ 78 ] = 7.125759e+02 , muup[ 78 ] = 8.395757e+02 ;
      xlim[ 79 ] = 1.081744e-01 , mulo[ 79 ] = 6.818009e+02 , muup[ 79 ] = 7.974780e+02 ;
      xlim[ 80 ] = 1.080166e-01 , mulo[ 80 ] = 6.547537e+02 , muup[ 80 ] = 7.602179e+02 ;
      xlim[ 81 ] = 1.081383e-01 , mulo[ 81 ] = 6.053877e+02 , muup[ 81 ] = 7.284519e+02 ;
      xlim[ 82 ] = 1.080570e-01 , mulo[ 82 ] = 5.672063e+02 , muup[ 82 ] = 6.755840e+02 ;
      xlim[ 83 ] = 1.080358e-01 , mulo[ 83 ] = 5.352622e+02 , muup[ 83 ] = 6.337936e+02 ;
      xlim[ 84 ] = 1.175774e-01 , mulo[ 84 ] = 1.115865e+03 , muup[ 84 ] = 1.609001e+03 ;
      xlim[ 85 ] = 1.175749e-01 , mulo[ 85 ] = 1.025059e+03 , muup[ 85 ] = 1.479689e+03 ;
      xlim[ 86 ] = 1.175980e-01 , mulo[ 86 ] = 9.483152e+02 , muup[ 86 ] = 1.330121e+03 ;
      xlim[ 87 ] = 1.176084e-01 , mulo[ 87 ] = 8.822787e+02 , muup[ 87 ] = 1.209659e+03 ;
      xlim[ 88 ] = 1.176275e-01 , mulo[ 88 ] = 8.316384e+02 , muup[ 88 ] = 1.117892e+03 ;
      xlim[ 89 ] = 1.175868e-01 , mulo[ 89 ] = 7.836943e+02 , muup[ 89 ] = 1.036613e+03 ;
      xlim[ 90 ] = 1.175979e-01 , mulo[ 90 ] = 7.449401e+02 , muup[ 90 ] = 9.768398e+02 ;
      xlim[ 91 ] = 1.175875e-01 , mulo[ 91 ] = 7.104760e+02 , muup[ 91 ] = 9.305579e+02 ;
      xlim[ 92 ] = 1.175875e-01 , mulo[ 92 ] = 6.818178e+02 , muup[ 92 ] = 8.872402e+02 ;
      xlim[ 93 ] = 1.176448e-01 , mulo[ 93 ] = 6.314395e+02 , muup[ 93 ] = 8.477351e+02 ;
      xlim[ 94 ] = 1.175995e-01 , mulo[ 94 ] = 5.907095e+02 , muup[ 94 ] = 7.877628e+02 ;
      xlim[ 95 ] = 1.176168e-01 , mulo[ 95 ] = 5.577124e+02 , muup[ 95 ] = 7.368155e+02 ;
      xlim[ 96 ] = 1.600280e-01 , mulo[ 96 ] = 1.301979e+03 , muup[ 96 ] = 1.721581e+03 ;
      xlim[ 97 ] = 1.600164e-01 , mulo[ 97 ] = 1.196334e+03 , muup[ 97 ] = 1.588456e+03 ;
      xlim[ 98 ] = 1.600537e-01 , mulo[ 98 ] = 1.106014e+03 , muup[ 98 ] = 1.427332e+03 ;
      xlim[ 99 ] = 1.600610e-01 , mulo[ 99 ] = 1.031670e+03 , muup[ 99 ] = 1.298634e+03 ;
      xlim[ 100 ] = 1.600377e-01 , mulo[ 100 ] = 9.675515e+02 , muup[ 100 ] = 1.193610e+03 ;
      xlim[ 101 ] = 1.600613e-01 , mulo[ 101 ] = 9.132725e+02 , muup[ 101 ] = 1.112169e+03 ;
      xlim[ 102 ] = 1.600855e-01 , mulo[ 102 ] = 8.682879e+02 , muup[ 102 ] = 1.050535e+03 ;
      xlim[ 103 ] = 1.600316e-01 , mulo[ 103 ] = 8.288235e+02 , muup[ 103 ] = 9.964844e+02 ;
      xlim[ 104 ] = 1.600316e-01 , mulo[ 104 ] = 7.949343e+02 , muup[ 104 ] = 9.513860e+02 ;
      xlim[ 105 ] = 1.600865e-01 , mulo[ 105 ] = 7.373336e+02 , muup[ 105 ] = 9.102386e+02 ;
      xlim[ 106 ] = 1.600558e-01 , mulo[ 106 ] = 6.899330e+02 , muup[ 106 ] = 8.437490e+02 ;
      xlim[ 107 ] = 1.602237e-01 , mulo[ 107 ] = 6.519156e+02 , muup[ 107 ] = 7.906797e+02 ;
      xlim[ 108 ] = 1.837267e-01 , mulo[ 108 ] = 1.395977e+03 , muup[ 108 ] = 1.834314e+03 ;
      xlim[ 109 ] = 1.837017e-01 , mulo[ 109 ] = 1.281249e+03 , muup[ 109 ] = 1.689585e+03 ;
      xlim[ 110 ] = 1.837104e-01 , mulo[ 110 ] = 1.185185e+03 , muup[ 110 ] = 1.522632e+03 ;
      xlim[ 111 ] = 1.837295e-01 , mulo[ 111 ] = 1.104970e+03 , muup[ 111 ] = 1.384471e+03 ;
      xlim[ 112 ] = 1.838069e-01 , mulo[ 112 ] = 1.038134e+03 , muup[ 112 ] = 1.267276e+03 ;
      xlim[ 113 ] = 1.837293e-01 , mulo[ 113 ] = 9.803749e+02 , muup[ 113 ] = 1.184308e+03 ;
      xlim[ 114 ] = 1.838375e-01 , mulo[ 114 ] = 9.317506e+02 , muup[ 114 ] = 1.118677e+03 ;
      xlim[ 115 ] = 1.837979e-01 , mulo[ 115 ] = 8.900752e+02 , muup[ 115 ] = 1.063128e+03 ;
      xlim[ 116 ] = 1.837471e-01 , mulo[ 116 ] = 8.530940e+02 , muup[ 116 ] = 1.013048e+03 ;
      xlim[ 117 ] = 1.836832e-01 , mulo[ 117 ] = 7.900507e+02 , muup[ 117 ] = 9.722671e+02 ;
      xlim[ 118 ] = 1.837471e-01 , mulo[ 118 ] = 7.385310e+02 , muup[ 118 ] = 8.996902e+02 ;
      xlim[ 119 ] = 1.837412e-01 , mulo[ 119 ] = 6.996902e+02 , muup[ 119 ] = 8.414069e+02 ;
      xlim[ 120 ] = 2.090425e-01 , mulo[ 120 ] = 1.487597e+03 , muup[ 120 ] = 2.269231e+03 ;
      xlim[ 121 ] = 2.090242e-01 , mulo[ 121 ] = 1.367516e+03 , muup[ 121 ] = 2.110569e+03 ;
      xlim[ 122 ] = 2.089889e-01 , mulo[ 122 ] = 1.262935e+03 , muup[ 122 ] = 1.883985e+03 ;
      xlim[ 123 ] = 2.090379e-01 , mulo[ 123 ] = 1.175420e+03 , muup[ 123 ] = 1.723509e+03 ;
      xlim[ 124 ] = 2.090604e-01 , mulo[ 124 ] = 1.108157e+03 , muup[ 124 ] = 1.588440e+03 ;
      xlim[ 125 ] = 2.090139e-01 , mulo[ 125 ] = 1.043982e+03 , muup[ 125 ] = 1.481340e+03 ;
      xlim[ 126 ] = 2.090241e-01 , mulo[ 126 ] = 9.929854e+02 , muup[ 126 ] = 1.398105e+03 ;
      xlim[ 127 ] = 2.090635e-01 , mulo[ 127 ] = 9.493163e+02 , muup[ 127 ] = 1.326228e+03 ;
      xlim[ 128 ] = 2.091161e-01 , mulo[ 128 ] = 9.121355e+02 , muup[ 128 ] = 1.267510e+03 ;
      xlim[ 129 ] = 2.090421e-01 , mulo[ 129 ] = 8.426684e+02 , muup[ 129 ] = 1.214176e+03 ;
      xlim[ 130 ] = 2.090139e-01 , mulo[ 130 ] = 7.898450e+02 , muup[ 130 ] = 1.127886e+03 ;
      xlim[ 131 ] = 2.090443e-01 , mulo[ 131 ] = 7.450683e+02 , muup[ 131 ] = 1.053236e+03 ;
      xlim[ 132 ] = 3.265548e-01 , mulo[ 132 ] = 1.860186e+03 , muup[ 132 ] = 2.632010e+03 ;
      xlim[ 133 ] = 3.265913e-01 , mulo[ 133 ] = 1.710122e+03 , muup[ 133 ] = 2.478940e+03 ;
      xlim[ 134 ] = 3.267028e-01 , mulo[ 134 ] = 1.581652e+03 , muup[ 134 ] = 2.288315e+03 ;
      xlim[ 135 ] = 3.266768e-01 , mulo[ 135 ] = 1.473005e+03 , muup[ 135 ] = 2.105574e+03 ;
      xlim[ 136 ] = 3.266602e-01 , mulo[ 136 ] = 1.381755e+03 , muup[ 136 ] = 1.962299e+03 ;
      xlim[ 137 ] = 3.265490e-01 , mulo[ 137 ] = 1.307293e+03 , muup[ 137 ] = 1.847643e+03 ;
      xlim[ 138 ] = 3.265705e-01 , mulo[ 138 ] = 1.240142e+03 , muup[ 138 ] = 1.741549e+03 ;
      xlim[ 139 ] = 3.266435e-01 , mulo[ 139 ] = 1.184665e+03 , muup[ 139 ] = 1.651825e+03 ;
      xlim[ 140 ] = 3.266274e-01 , mulo[ 140 ] = 1.136246e+03 , muup[ 140 ] = 1.576987e+03 ;
      xlim[ 141 ] = 3.266750e-01 , mulo[ 141 ] = 1.053528e+03 , muup[ 141 ] = 1.510050e+03 ;
      xlim[ 142 ] = 3.266563e-01 , mulo[ 142 ] = 9.883950e+02 , muup[ 142 ] = 1.400557e+03 ;
      xlim[ 143 ] = 3.266586e-01 , mulo[ 143 ] = 9.311350e+02 , muup[ 143 ] = 1.322186e+03 ;
      xlim[ 144 ] = 5.102258e-01 , mulo[ 144 ] = 2.323464e+03 , muup[ 144 ] = 3.498976e+03 ;
      xlim[ 145 ] = 5.102258e-01 , mulo[ 145 ] = 2.140177e+03 , muup[ 145 ] = 3.297241e+03 ;
      xlim[ 146 ] = 5.102980e-01 , mulo[ 146 ] = 1.976208e+03 , muup[ 146 ] = 3.025862e+03 ;
      xlim[ 147 ] = 5.102078e-01 , mulo[ 147 ] = 1.842157e+03 , muup[ 147 ] = 2.793854e+03 ;
      xlim[ 148 ] = 5.104406e-01 , mulo[ 148 ] = 1.727088e+03 , muup[ 148 ] = 2.605578e+03 ;
      xlim[ 149 ] = 5.104831e-01 , mulo[ 149 ] = 1.630159e+03 , muup[ 149 ] = 2.444859e+03 ;
      xlim[ 150 ] = 5.103696e-01 , mulo[ 150 ] = 1.554251e+03 , muup[ 150 ] = 2.312008e+03 ;
      xlim[ 151 ] = 5.103140e-01 , mulo[ 151 ] = 1.482670e+03 , muup[ 151 ] = 2.193219e+03 ;
      xlim[ 152 ] = 5.105427e-01 , mulo[ 152 ] = 1.422654e+03 , muup[ 152 ] = 2.094378e+03 ;
      xlim[ 153 ] = 5.103140e-01 , mulo[ 153 ] = 1.313974e+03 , muup[ 153 ] = 2.009439e+03 ;
      xlim[ 154 ] = 5.103357e-01 , mulo[ 154 ] = 1.232926e+03 , muup[ 154 ] = 1.860840e+03 ;
      xlim[ 155 ] = 5.106212e-01 , mulo[ 155 ] = 1.165936e+03 , muup[ 155 ] = 1.741017e+03 ;
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
  B->NScaleDim = 1; // one variable used in scales: dijet pT average
  B->Iscale.push_back(0);  // mur=mur(pT), pT = index 0 
  B->Iscale.push_back(0);  // muf=muf(pT), pT = index 0 
  B->ScaleDescript.resize(B->NScaleDim);

  B->ScaleDescript[0].push_back("<pT_1,2>"); // --- fastNLO user: give name for selected scale
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
