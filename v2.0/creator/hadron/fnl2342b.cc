//
// fastNLO v2 author code for fnl2342b:
//     CMS LHC Inclusive Jets Scenario, E_cms = 7 TeV
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
// 2011/01/26 KR Add CMS Inclusive Jets Scenario
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
   
  bool doReference;
  bool doWarmUp;
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


void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
  fnloBlockA2 *A2 =  table->GetBlockA2();
  double x1 = p[-1].Z()/p[hadron(-1)].Z();
  double x2 = p[0].Z()/p[hadron(0)].Z();

  // --- run the jet algorithm
  double jetsize = 0.5;
  pj = jetclus(p,jetsize);
  unsigned int nj = pj.upper(); 

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
  // lowest pT for jets to be considered
  double ptjmin = 18.;
  // highest (pseudo-)rapidity for jets to be considered
  double yjmax  = 3.0;

  // Analyze inclusive jets in jet loop
  for (unsigned int i = 1; i <= nj; i++) {
    double pt  = pj[i].perp(); 
    double rap = abs(pj[i].rapidity());
    
    // --- jet in phase space?
    if (ptjmin < pt && rap < yjmax) {

      // - set the renormalization and factorization scale to jet pT
      double mu = pt;

      // --- identify bin number (dim1,dim2) e.g. (pT,y)
      int obsbin = -1;
      for (int j = 0; j < A2->GetNObsBin(); j++) {
	if (A2->LoBin[j][0] <= pt  && pt  < A2->UpBin[j][0] && 
	    A2->LoBin[j][1] <= rap && rap < A2->UpBin[j][1]) {
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
    } // --- end: phase space selection
  } // --- end: jet loop
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

  // --- determine whether this is a warm-up or reference job
  doWarmUp = false;
  if (strstr(file,"wrm")!=NULL) {
    doWarmUp = true;
    printf("fastNLO: This is a warm-up run!\n");
    if ( ! nlo ) {
      printf("fastNLO: WARNING! Warm-up runs are better done at NLO!\n");
    }
  }
  doReference = false;
  if (strstr(file,"ref")!=NULL) {
    doReference = true;
    printf("fastNLO: This is a reference run!\n");
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
  table->GetBlockA1()->SetScenName("fnl2342b");  // - fastNLO user: set scenario name
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
  A2->ScDescript.push_back("d2sigma-jet_dpTd|y|_(pb_GeV)");
  A2->ScDescript.push_back("CMS_Collaboration");
  A2->ScDescript.push_back("Inclusive_Jet_pT");
  A2->ScDescript.push_back("anti-kT_R=0.5");
  A2->ScDescript.push_back("arXiv:1106.0208");

  A2->NScDescript = A2->ScDescript.size();
  A2->Ecms = sqrt(s);
  A2->ILOord = 2;   // --- fastNLO user: power of LO contribution for process
  A2->NDim = 2;     // --- fastNLO user: No of dimensions in which observable is binned
  A2->DimLabel.push_back("pT_[GeV]");  // --- fastNLO user: label of 1st dimension
  A2->IDiffBin.push_back(2);
  A2->DimLabel.push_back("|y|");   // --- fastNLO user: label of 2nd dimension
  A2->IDiffBin.push_back(2);

  vector <double> bound;
  bound.resize(2);

  // --- fastNLO user: bin definitions - here in |y| and pT
  const int nrapbins = 6;
  double rapbins[nrapbins+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
  
  const int nptbins[nrapbins] = {34, 33, 32, 29, 26, 22};

  cout << "------------------------" << endl;
  cout << "Binning in dimension 1: " << endl;
  for (int i=0; i<nrapbins+1; i++) {
    cout << "i, rapbins: " << i << ", " << rapbins[i] << endl;
  }

  vector< vector<double> >ptbins;
  ptbins.resize(nrapbins);
  for (int i=0; i<nrapbins; i++) {
    ptbins[i].resize(nptbins[i]+1);
  }
  double ptb0[35] = { 18., 21., 24., 28., 32., 37., 43., 49., 56., 64., 74., 84.,
		      97., 114., 133., 153., 174., 196., 220., 245., 272., 300., 330., 362.,
		      395., 430., 468., 507., 548., 592., 638., 686., 737., 846., 1684. };
  double ptb1[34] = { 18., 21., 24., 28., 32., 37., 43., 49., 56., 64., 74., 84.,
		      97., 114., 133., 153., 174., 196., 220., 245., 272., 300., 330., 362.,
		      395., 430., 468., 507., 548., 592., 638., 686., 790., 1684. };
  double ptb2[33] = { 18., 21., 24., 28., 32., 37., 43., 49., 56., 64., 74., 84.,
		      97., 114., 133., 153., 174., 196., 220., 245., 272., 300., 330., 362.,
		      395., 430., 468., 507., 548., 592., 638., 686., 1410. };
  double ptb3[30] = { 18., 21., 24., 28., 32., 37., 43., 49., 56., 64., 74., 84.,
		      97., 114., 133., 153., 174., 196., 220., 245., 272., 300., 330., 362.,
		      395., 430., 468., 507., 548., 1032. };
  double ptb4[27] = { 18., 21., 24., 28., 32., 37., 43., 49., 56., 64., 74., 84.,
		      97., 114., 133., 153., 174., 196., 220., 245., 272., 300., 330., 362.,
		      395., 430., 737. };
  double ptb5[23] = { 18., 21., 24., 28., 32., 37., 43., 49., 56., 64., 74., 84.,
		      97., 114., 133., 153., 174., 196., 220., 245., 272., 300., 468. };
  for (int j=0; j<nptbins[0]+1; j++) { 
    ptbins[0][j] = ptb0[j];
  }
  for (int j=0; j<nptbins[1]+1; j++) { 
    ptbins[1][j] = ptb1[j];
  }
  for (int j=0; j<nptbins[2]+1; j++) { 
    ptbins[2][j] = ptb2[j];
  }
  for (int j=0; j<nptbins[3]+1; j++) { 
    ptbins[3][j] = ptb3[j];
  }
  for (int j=0; j<nptbins[4]+1; j++) { 
    ptbins[4][j] = ptb4[j];
  }
  for (int j=0; j<nptbins[5]+1; j++) { 
    ptbins[5][j] = ptb5[j];
  }

  cout << "------------------------" << endl;
  cout << "Binning in dimension 2: " << endl;
  for (int i=0; i<nrapbins; i++) {
    for (int j=0; j<nptbins[i]+1; j++) {
      cout << "i, j, ptbins: " << i << ", " << j << ", " << ptbins[i][j] << endl;
    }
  }
  cout << "------------------------" << endl;



  // --- fastNLO user:
  //     define below the bin width ("binsize") by which
  //     the cross section is divided to obtain the 
  //     (multi-) differential result.
  // fnl2342b: divide by bin width in pT and |y|

  int nbins = 0;   // --- count total No. bins
  for (int i=0;i<nrapbins;i++){
    for (int j=0;j<nptbins[i];j++){
      double binsize = 1.;
      nbins += 1;
      bound[0] = ptbins[i][j];
      bound[1] = rapbins[i];
      A2->LoBin.push_back(bound);
      bound[0] = ptbins[i][j+1];
      bound[1] = rapbins[i+1];
      A2->UpBin.push_back(bound);
      //printf(" %d %d  |  %f %f\n",i,j,bound[0],bound[1]);
      binsize = binsize // fnl2342b: Start with preset value = 1 
	* (ptbins[i][j+1]-ptbins[i][j]) // ... times dpT
	* 2. * (rapbins[i+1]-rapbins[i]); // ... times d|y|
      A2->BinSize.push_back(binsize);
    }
  }
  printf("fastNLO: Total no. of observable bins = %d\n",nbins+1);

  A2->NObsBin = nbins;
  A2->INormFlag = 0;   // --- fastNLO user: default=0 - set =1 if observable is 
  //     to be normalized by own integral (in 1st dimension)
  //     see documentation for details and for other options

  // --- fastNLO table block B
  fnloBlockBNlojet *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
  table->CreateBlockB(0,B);
  B->IXsectUnits = 12;    // --- fastNLO user: set to same value as "SetIpublunits"
  B->IDataFlag = 0;
  B->IAddMultFlag = 0;
  B->IContrFlag1 = 1;
  B->IContrFlag3 = 0;
  B->CodeDescript.push_back("NLOJet++_4.1.3");  // --- fastNLO user: enter NLOJET++ version
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
  B->IWarmUpPrint = 10000000;
  //B->IWarmUpPrint = 10000;
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
    for(int i=0;i<A2->NObsBin;i++){
      xlim[i] = 1.1e-07, mulo[i] = 3.0, muup[i]=9.9e10; // - safe initializations
    }
  } else {
    FILE * infile;
    infile = fopen("fastNLO-warmup.dat","r");
    if ( ! infile ) {
      cerr << "fastNLO: ERROR! Could not read x limits from file: fastNLO-warmup.dat" << endl;
      exit(1);
    }
    char line[256];
    // Ignore first documentation line
    if ( ! fgets(line,sizeof(line),infile) ) {
      cerr << "fastNLO: ERROR! Reading empty file: fastNLO-warmup.dat" << endl;
      exit(1);
    }
    printf("fastNLO: Reading x limits from file: fastNLO-warmup.dat\n");
    printf("%s",line);
    // Now read and print out all limits
    int i = 0;
    while ( fgets(line,sizeof(line),infile) ) {
      sscanf(line,"%lf, %lf, %lf;",&xlim[i],&mulo[i],&muup[i]);
      printf(" xlim[%d]=%7.5f, mulo[%d]=%8.3f, muup[%d]=%8.3f;\n",i,xlim[i],i,mulo[i],i,muup[i]);
      i++;
    }
    if (i != A2->NObsBin ) {
      cerr << "fastNLO: ERROR! Number of x limits read != NObsBin: i = " << i << ", NObsBin = " << A2->NObsBin << endl;
      exit(1);
    }
  }

  // --------- fastNLO: Warm-Up run results (start)
  // if ( ! doWarmUp ) {
  //  // 16160000000 contributions (!= events) in warm-up run
  // 1.570869e-03, 1.800000e+01, 2.100000e+01;
  // 1.842185e-03, 2.100000e+01, 2.400000e+01;
  // 2.102049e-03, 2.400000e+01, 2.800000e+01;
  // 2.453217e-03, 2.800000e+01, 3.200000e+01;
  // 2.807311e-03, 3.200000e+01, 3.700000e+01;
  // 3.247694e-03, 3.700000e+01, 4.300000e+01;
  // 3.786345e-03, 4.300000e+01, 4.900000e+01;
  // 4.319395e-03, 4.900000e+01, 5.600000e+01;
  // 4.949800e-03, 5.600000e+01, 6.400000e+01;
  // 5.673713e-03, 6.400000e+01, 7.400000e+01;
  // 6.545409e-03, 7.400000e+01, 8.400000e+01;
  // 7.427902e-03, 8.400000e+01, 9.700000e+01;
  // 8.647304e-03, 9.700000e+01, 1.140000e+02;
  // 1.017103e-02, 1.140000e+02, 1.330000e+02;
  // 1.192320e-02, 1.330000e+02, 1.530000e+02;
  // 1.376966e-02, 1.530000e+02, 1.740000e+02;
  // 1.575296e-02, 1.740000e+02, 1.960000e+02;
  // 1.782699e-02, 1.960000e+02, 2.200000e+02;
  // 2.014816e-02, 2.200000e+02, 2.450000e+02;
  // 2.255938e-02, 2.450000e+02, 2.720000e+02;
  // 2.521716e-02, 2.720000e+02, 3.000000e+02;
  // 2.800307e-02, 3.000000e+02, 3.300000e+02;
  // 3.104837e-02, 3.300000e+02, 3.620000e+02;
  // 3.436577e-02, 3.620000e+02, 3.950000e+02;
  // 3.779190e-02, 3.950000e+02, 4.300000e+02;
  // 4.153333e-02, 4.300000e+02, 4.680000e+02;
  // 4.565891e-02, 4.680000e+02, 5.070000e+02;
  // 4.997736e-02, 5.070000e+02, 5.480000e+02;
  // 5.460218e-02, 5.480000e+02, 5.920000e+02;
  // 5.972602e-02, 5.920000e+02, 6.380000e+02;
  // 6.513817e-02, 6.380000e+02, 6.860000e+02;
  // 7.097763e-02, 6.860000e+02, 7.370000e+02;
  // 7.742613e-02, 7.370000e+02, 8.460000e+02;
  // 9.166592e-02, 8.460000e+02, 1.684000e+03;
  // 9.586801e-04, 1.800000e+01, 2.100000e+01;
  // 1.121244e-03, 2.100000e+01, 2.400000e+01;
  // 1.280398e-03, 2.400000e+01, 2.800000e+01;
  // 1.499323e-03, 2.800000e+01, 3.200000e+01;
  // 1.708473e-03, 3.200000e+01, 3.700000e+01;
  // 1.983315e-03, 3.700000e+01, 4.300000e+01;
  // 2.312776e-03, 4.300000e+01, 4.900000e+01;
  // 2.646720e-03, 4.900000e+01, 5.600000e+01;
  // 3.012052e-03, 5.600000e+01, 6.400000e+01;
  // 3.457834e-03, 6.400000e+01, 7.400000e+01;
  // 4.037938e-03, 7.400000e+01, 8.400000e+01;
  // 4.599814e-03, 8.400000e+01, 9.700000e+01;
  // 5.314140e-03, 9.700000e+01, 1.140000e+02;
  // 6.312425e-03, 1.140000e+02, 1.330000e+02;
  // 7.398199e-03, 1.330000e+02, 1.530000e+02;
  // 8.577277e-03, 1.530000e+02, 1.740000e+02;
  // 9.897978e-03, 1.740000e+02, 1.960000e+02;
  // 1.116435e-02, 1.960000e+02, 2.200000e+02;
  // 1.266586e-02, 2.200000e+02, 2.450000e+02;
  // 1.425817e-02, 2.450000e+02, 2.720000e+02;
  // 1.601715e-02, 2.720000e+02, 3.000000e+02;
  // 1.786442e-02, 3.000000e+02, 3.300000e+02;
  // 1.992947e-02, 3.300000e+02, 3.620000e+02;
  // 2.218864e-02, 3.620000e+02, 3.950000e+02;
  // 2.458264e-02, 3.950000e+02, 4.300000e+02;
  // 2.715063e-02, 4.300000e+02, 4.680000e+02;
  // 3.007401e-02, 4.680000e+02, 5.070000e+02;
  // 3.322121e-02, 5.070000e+02, 5.480000e+02;
  // 3.666103e-02, 5.480000e+02, 5.920000e+02;
  // 4.046824e-02, 5.920000e+02, 6.380000e+02;
  // 4.467488e-02, 6.380000e+02, 6.860000e+02;
  // 4.918692e-02, 6.860000e+02, 7.900000e+02;
  // 5.996411e-02, 7.900000e+02, 1.684000e+03;
  // 5.813668e-04, 1.800000e+01, 2.100000e+01;
  // 6.798892e-04, 2.100000e+01, 2.400000e+01;
  // 7.836975e-04, 2.400000e+01, 2.800000e+01;
  // 9.161333e-04, 2.800000e+01, 3.200000e+01;
  // 1.045291e-03, 3.200000e+01, 3.700000e+01;
  // 1.220038e-03, 3.700000e+01, 4.300000e+01;
  // 1.415476e-03, 4.300000e+01, 4.900000e+01;
  // 1.618170e-03, 4.900000e+01, 5.600000e+01;
  // 1.859247e-03, 5.600000e+01, 6.400000e+01;
  // 2.141280e-03, 6.400000e+01, 7.400000e+01;
  // 2.483203e-03, 7.400000e+01, 8.400000e+01;
  // 2.850001e-03, 8.400000e+01, 9.700000e+01;
  // 3.312549e-03, 9.700000e+01, 1.140000e+02;
  // 3.943816e-03, 1.140000e+02, 1.330000e+02;
  // 4.655653e-03, 1.330000e+02, 1.530000e+02;
  // 5.425006e-03, 1.530000e+02, 1.740000e+02;
  // 6.258805e-03, 1.740000e+02, 1.960000e+02;
  // 7.168202e-03, 1.960000e+02, 2.200000e+02;
  // 8.210800e-03, 2.200000e+02, 2.450000e+02;
  // 9.322299e-03, 2.450000e+02, 2.720000e+02;
  // 1.051796e-02, 2.720000e+02, 3.000000e+02;
  // 1.185194e-02, 3.000000e+02, 3.300000e+02;
  // 1.335701e-02, 3.300000e+02, 3.620000e+02;
  // 1.504748e-02, 3.620000e+02, 3.950000e+02;
  // 1.688223e-02, 3.950000e+02, 4.300000e+02;
  // 1.894307e-02, 4.300000e+02, 4.680000e+02;
  // 2.132506e-02, 4.680000e+02, 5.070000e+02;
  // 2.397109e-02, 5.070000e+02, 5.480000e+02;
  // 2.695367e-02, 5.480000e+02, 5.920000e+02;
  // 3.048928e-02, 5.920000e+02, 6.380000e+02;
  // 3.445020e-02, 6.380000e+02, 6.860000e+02;
  // 3.904010e-02, 6.860000e+02, 1.410000e+03;
  // 3.575703e-04, 1.800000e+01, 2.100000e+01;
  // 4.176482e-04, 2.100000e+01, 2.400000e+01;
  // 4.804413e-04, 2.400000e+01, 2.800000e+01;
  // 5.614078e-04, 2.800000e+01, 3.200000e+01;
  // 6.420933e-04, 3.200000e+01, 3.700000e+01;
  // 7.509800e-04, 3.700000e+01, 4.300000e+01;
  // 8.768687e-04, 4.300000e+01, 4.900000e+01;
  // 1.006723e-03, 4.900000e+01, 5.600000e+01;
  // 1.163303e-03, 5.600000e+01, 6.400000e+01;
  // 1.339278e-03, 6.400000e+01, 7.400000e+01;
  // 1.559134e-03, 7.400000e+01, 8.400000e+01;
  // 1.801676e-03, 8.400000e+01, 9.700000e+01;
  // 2.108824e-03, 9.700000e+01, 1.140000e+02;
  // 2.525561e-03, 1.140000e+02, 1.330000e+02;
  // 3.002230e-03, 1.330000e+02, 1.530000e+02;
  // 3.542622e-03, 1.530000e+02, 1.740000e+02;
  // 4.134644e-03, 1.740000e+02, 1.960000e+02;
  // 4.815718e-03, 1.960000e+02, 2.200000e+02;
  // 5.584594e-03, 2.200000e+02, 2.450000e+02;
  // 6.406899e-03, 2.450000e+02, 2.720000e+02;
  // 7.390371e-03, 2.720000e+02, 3.000000e+02;
  // 8.509150e-03, 3.000000e+02, 3.300000e+02;
  // 9.828682e-03, 3.300000e+02, 3.620000e+02;
  // 1.135391e-02, 3.620000e+02, 3.950000e+02;
  // 1.311169e-02, 3.950000e+02, 4.300000e+02;
  // 1.523028e-02, 4.300000e+02, 4.680000e+02;
  // 1.790272e-02, 4.680000e+02, 5.070000e+02;
  // 2.099758e-02, 5.070000e+02, 5.480000e+02;
  // 2.453239e-02, 5.480000e+02, 1.032000e+03;
  // 2.202036e-04, 1.800000e+01, 2.100000e+01;
  // 2.581491e-04, 2.100000e+01, 2.400000e+01;
  // 2.954210e-04, 2.400000e+01, 2.800000e+01;
  // 3.482940e-04, 2.800000e+01, 3.200000e+01;
  // 3.992407e-04, 3.200000e+01, 3.700000e+01;
  // 4.687384e-04, 3.700000e+01, 4.300000e+01;
  // 5.469676e-04, 4.300000e+01, 4.900000e+01;
  // 6.323534e-04, 4.900000e+01, 5.600000e+01;
  // 7.297971e-04, 5.600000e+01, 6.400000e+01;
  // 8.512513e-04, 6.400000e+01, 7.400000e+01;
  // 1.005413e-03, 7.400000e+01, 8.400000e+01;
  // 1.162369e-03, 8.400000e+01, 9.700000e+01;
  // 1.378690e-03, 9.700000e+01, 1.140000e+02;
  // 1.681696e-03, 1.140000e+02, 1.330000e+02;
  // 2.036906e-03, 1.330000e+02, 1.530000e+02;
  // 2.464163e-03, 1.530000e+02, 1.740000e+02;
  // 2.939006e-03, 1.740000e+02, 1.960000e+02;
  // 3.499306e-03, 1.960000e+02, 2.200000e+02;
  // 4.208564e-03, 2.200000e+02, 2.450000e+02;
  // 5.042398e-03, 2.450000e+02, 2.720000e+02;
  // 6.076901e-03, 2.720000e+02, 3.000000e+02;
  // 7.364254e-03, 3.000000e+02, 3.300000e+02;
  // 8.911763e-03, 3.300000e+02, 3.620000e+02;
  // 1.070827e-02, 3.620000e+02, 3.950000e+02;
  // 1.274154e-02, 3.950000e+02, 4.300000e+02;
  // 1.510773e-02, 4.300000e+02, 7.370000e+02;
  // 1.373663e-04, 1.800000e+01, 2.100000e+01;
  // 1.607552e-04, 2.100000e+01, 2.400000e+01;
  // 1.839888e-04, 2.400000e+01, 2.800000e+01;
  // 2.187672e-04, 2.800000e+01, 3.200000e+01;
  // 2.531392e-04, 3.200000e+01, 3.700000e+01;
  // 2.985629e-04, 3.700000e+01, 4.300000e+01;
  // 3.511622e-04, 4.300000e+01, 4.900000e+01;
  // 4.093714e-04, 4.900000e+01, 5.600000e+01;
  // 4.786728e-04, 5.600000e+01, 6.400000e+01;
  // 5.613578e-04, 6.400000e+01, 7.400000e+01;
  // 6.706556e-04, 7.400000e+01, 8.400000e+01;
  // 7.924434e-04, 8.400000e+01, 9.700000e+01;
  // 9.675314e-04, 9.700000e+01, 1.140000e+02;
  // 1.210841e-03, 1.140000e+02, 1.330000e+02;
  // 1.539389e-03, 1.330000e+02, 1.530000e+02;
  // 1.951574e-03, 1.530000e+02, 1.740000e+02;
  // 2.481215e-03, 1.740000e+02, 1.960000e+02;
  // 3.152265e-03, 1.960000e+02, 2.200000e+02;
  // 3.965797e-03, 2.200000e+02, 2.450000e+02;
  // 4.902888e-03, 2.450000e+02, 2.720000e+02;
  // 6.058142e-03, 2.720000e+02, 3.000000e+02;
  // 7.393679e-03, 3.000000e+02, 4.680000e+02;
  // }
  // --------- fastNLO: Warm-Up run results (end)

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
  B->NScaleDim = 1; // one variable used in scales: pT
  B->Iscale.push_back(0);  // mur=mur(pT), pT = index 0 
  B->Iscale.push_back(0);  // muf=muf(pT), pT = index 0 
  B->ScaleDescript.resize(B->NScaleDim);

  B->ScaleDescript[0].push_back("pT");
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
    fnloBlockBNlojet *refB = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
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
