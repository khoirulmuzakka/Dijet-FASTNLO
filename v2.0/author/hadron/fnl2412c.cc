//
// fastNLO v2 author code for fnl2412c:
//     CMS LHC Dijet Mass Scenario, E_cms = 7 TeV
//     for fastjet anti-kT algo with R=0.7 in E-scheme
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

#include "fj-ak.h"

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
  fj_ak jetclus;
   
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
  double jetsize = 0.7;
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
  double ptjmin = 30.;
  // lowest pT of leading jet for event to be considered (trigger threshold)
  double ptj1min = 60.;
  // highest (pseudo-)rapidity for jets to be considered
  double yjmax  = 3.0;



  
  // Find the two leading jets within phase space cuts
  int njet = 0;
  int ijet[4]     = {0, 0, 0, 0};
  double ptjet[4] = {0.,0.,0.,0.};
  double yjet[4]  = {10.,10.,10.,10.};
  if (nj > 0) {
    // Initialize pointers to the jets, check minimal jet pT and maximal |y,eta|
    for (int i=1; i<=nj; i++) {
      double pti = pj[i].perp();
      double yi  = abs(pj[i].rapidity());
      //DEBUG
      //       cout << "ijet, pti, yi: " << i << ", " << pti << ", " << yi << endl;
      //       cout << "ptlow = " << ptlow << ", raphigh = " << yjmax << endl;
      //DEBUGEND
      if (pti > ptjmin && yi < yjmax) {
	ijet[i]  = 1;
	ptjet[i] = pti; 
	yjet[i]  = yi;
	njet++;
      }
    }
  }

  //DEBUG
  //   cout << "nj, njet = " << nj << ", " << njet << endl;
  //DEBUGEND
  if (njet >= 2) {
    // Get first two jets 
    int ij1 = 0;
    int ij2 = 0;
    for (int i=1; i<=nj; i++) {
      if (ijet[i] == 1) {
	if (ij1 == 0) {
	  ij1 = i;
	} else {
	  if (ij2 == 0) {
	    ij2 = i;
	  }
	}
      }
    }

    // Find the two leading jets
    // Order pt1 >= pt2
    double pt1 = pj[ij1].perp();
    double pt2 = pj[ij2].perp();
    double y1  = abs(pj[ij1].rapidity());
    double y2  = abs(pj[ij2].rapidity());
    if (pt2 > pt1) {
      int itmp     = ij1;
      double pttmp = pt1;
      double ytmp  = y1;
      ij1 = ij2;
      pt1 = pt2;
      y1  = y2;
      ij2 = itmp;
      pt2 = pttmp;
      y2  = ytmp;
    }
    // For 3 jets with pt3 > pt2, exchange 2 and 3
    if (njet > 2 && pj[3].perp() > pt2) {
      ij2 = 3;
      pt2 = pj[ij2].perp();
      y2  = abs(pj[ij2].rapidity());
    }
    // No further check whether pt3 > pt1, i.e. here order not important anymore
    //DEBUG
    //     cout << "ij1, pt1, y1: " << ij1 << ", " << pt1 << ", " << y1 << endl;
    //     cout << "ij2, pt2, y2: " << ij2 << ", " << pt2 << ", " << y2 << endl;
    //DEBUGEND
    
    // Derive dijet variables
    double mjj = sqrt((pj[ij1].T()+pj[ij2].T())*(pj[ij1].T()+pj[ij2].T())
		      -(pj[ij1].X()+pj[ij2].X())*(pj[ij1].X()+pj[ij2].X())
		      -(pj[ij1].Y()+pj[ij2].Y())*(pj[ij1].Y()+pj[ij2].Y())
		      -(pj[ij1].Z()+pj[ij2].Z())*(pj[ij1].Z()+pj[ij2].Z()));
    
    // Determine minimal and maximal (pseudo-)rapidity and maximal pT
    double yjjmin = min(y1,y2);
    double yjjmax = max(y1,y2);
    double ptmax  = max(pt1,pt2);

    // --- Further dijets phase space cuts?
    if ( ptj1min < ptmax && mjjbin[0][0] < mjj ) {
      
      // - set the renormalization and factorization scale to average dijet pT
      double mu = (pt1+pt2)/2.0;

      // --- identify bin number (dim1,dim2) here (mjj,|ymax|)
      int obsbin = -1;
      for (int j = 0; j < A2->GetNObsBin(); j++) {
	if (A2->LoBin[j][0] <= mjj    && mjj    < A2->UpBin[j][0] && 
	    A2->LoBin[j][1] <= yjjmax && yjjmax < A2->UpBin[j][1]) {
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
  } // --- end: dijet events only
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
  table->GetBlockA1()->SetScenName("fnl2412c");  // - fastNLO user: set scenario name
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
  A2->ScDescript.push_back("d2sigma-dijet_dMd|y_max|_(pb_GeV)");
  A2->ScDescript.push_back("CMS_Collaboration");
  A2->ScDescript.push_back("Dijet_Mass");
  A2->ScDescript.push_back("anti-kT_R=0.7");
  A2->ScDescript.push_back("arXiv:1104.1693");

  A2->NScDescript = A2->ScDescript.size();
  A2->Ecms = sqrt(s);
  A2->ILOord = 2;   // --- fastNLO user: power of LO contribution for process
  A2->NDim = 2;     // --- fastNLO user: No of dimensions in which observable is binned
  A2->DimLabel.push_back("pT_[GeV]");  // --- fastNLO user: label of 1st dimension
  A2->IDiffBin.push_back(2);
  A2->DimLabel.push_back("|y_max|");   // --- fastNLO user: label of 2nd dimension
  A2->IDiffBin.push_back(2);

  vector <double> bound;
  bound.resize(2);

  // --- fastNLO user: bin definitions - here in |y_max| and mjj
  const int nrapbins = 5;
  double rapbins[nrapbins+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
  
  const int nmjjbins[nrapbins] = {32, 31, 24, 20, 18};
  
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
  double mjj0[33] = { 156.0, 176.0, 197.0, 220.0, 244.0, 270.0, 296.0, 325.0, 354.0, 386.0, 
		      419.0, 453.0, 489.0, 526.0, 565.0, 606.0, 649.0, 693.0, 740.0, 788.0, 
		      838.0, 890.0, 944.0, 1000.0, 1118.0, 1246.0, 1383.0, 1530.0, 1687.0, 1856.0,
		      2037.0, 2332.0, 2659.0 };
  double mjj1[32] = { 197.0, 220.0, 244.0, 270.0, 296.0, 325.0, 354.0, 386.0, 419.0, 453.0,
		      489.0, 526.0, 565.0, 606.0, 649.0, 693.0, 740.0, 788.0, 838.0, 890.0,
		      944.0, 1000.0, 1118.0, 1246.0, 1383.0, 1530.0, 1687.0, 1856.0, 2037.0, 2332.0,
		      2659.0, 3019.0 };
  double mjj2[25] = { 386.0, 419.0, 453.0, 489.0, 526.0, 565.0, 606.0, 649.0, 693.0, 740.0,
		      788.0, 838.0, 890.0, 944.0, 1000.0, 1118.0, 1246.0, 1383.0, 1530.0, 1687.0,
		      1856.0, 2037.0, 2332.0, 2659.0, 3019.0 };
  double mjj3[21] = { 565.0, 606.0, 649.0, 693.0, 740.0, 788.0, 838.0, 890.0, 944.0, 1000.0,
		      1118.0, 1246.0, 1383.0, 1530.0, 1687.0, 1856.0, 2037.0, 2332.0, 2659.0, 3019.0,
		      3854.0 };
  double mjj4[19] = { 649.0, 693.0, 740.0, 788.0, 838.0, 890.0, 944.0, 1000.0, 1118.0, 1246.0,
		      1383.0, 1530.0, 1687.0, 1856.0, 2037.0, 2332.0, 2659.0, 3019.0, 3854.0 };
  for (int j=0; j<nmjjbins[0]+1; j++) { 
    mjjbins[0][j] = mjj0[j];
  }
  for (int j=0; j<nmjjbins[1]+1; j++) { 
    mjjbins[1][j] = mjj1[j];
  }
  for (int j=0; j<nmjjbins[2]+1; j++) { 
    mjjbins[2][j] = mjj2[j];
  }
  for (int j=0; j<nmjjbins[3]+1; j++) { 
    mjjbins[3][j] = mjj3[j];
  }
  for (int j=0; j<nmjjbins[4]+1; j++) { 
    mjjbins[4][j] = mjj4[j];
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
  // fnl2412c: divide by bin width in pT and |y|

  int nbins = 0;   // --- count total No. bins
  for (int i=0;i<nrapbins;i++){
    for (int j=0;j<nmjjbins[i];j++){
      double binsize = 1.;
      nbins += 1;
      bound[0] = mjjbins[i][j];
      bound[1] = rapbins[i];
      A2->LoBin.push_back(bound);
      bound[0] = mjjbins[i][j+1];
      bound[1] = rapbins[i+1];
      A2->UpBin.push_back(bound);
      //printf(" %d %d  |  %f %f\n",i,j,bound[0],bound[1]);
      binsize = binsize // fnl2412c: Start with preset value = 1 
	* (mjjbins[i][j+1]-mjjbins[i][j]) // ... times dpT
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
  B->NScaleDim = 1; // one variable used in scales: dijet pT average
  B->Iscale.push_back(0);  // mur=mur(pT), pT = index 0 
  B->Iscale.push_back(0);  // muf=muf(pT), pT = index 0 
  B->ScaleDescript.resize(B->NScaleDim);

  B->ScaleDescript[0].push_back("<pT_1,2>");
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
