//
// fastNLO v2.1 author code for fnl5002:
//     ATLAS LHC Dijet Mass Scenario, E_cms = 7 TeV
//     for fastjet anti-kT algo with R=0.6 in E-scheme
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
// Please keep the order of all statements in inittable
// in order to guarantee a working code.
//
// This file contains the following routines:
//   inputfunc		(-> user edits)
//   psinput		(-> user edits)
//   userfunc		(-> user edits)
//   FillEvent		(-> user edits)	
//   inittable		(-> user edits)
//   GetWarmupValues	(-> user edits)
//   DefineBinning	(-> user edits)
//   initfunc		(don't touch)
//   writetable		(don't touch)
//   end_of_event	(don't touch)
//   phys_output	(don't touch)
//   GetEcms		(don't touch)
//   GetNj		(don't touch)
//
// Implementing a new scenario requires to edit:
//  - the jet algorithm ("#include" statement and assignment of "jetclus")
//  - number of jets (inputfunc)
//  - center-of-mass energy (psinput)
//  - compute observable, determine bin no. (userfunc)
//  - declare all variables for table, define bin boundaries (inittable, etc.)
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
  pdf_cteq6_pp pdf_pp;
  pdf_hhc_dummy dummypdf;

  // --- jet algorithm
  fj_ak jetclus;   // fastNLO user: define jet algorithm (consistent with .h file above)
   
  // --- fastNLO definitions (not for user)
  fnloTable *table;
  double nevents;        // No of events calculated so far
  unsigned long nwrite;  // No of events after to write out the table
  string tablefilename;  // The table file to write to
  time_t start_time;
  bool nlo;
  double* xlim;
  double* scale1lo;
  double* scale1hi;
  double* scale2lo;
  double* scale2hi;
   
  void inittable();
  void writetable();
  void GetWarmupValues( fnloBlockBNlojet* B );
  void DefineBinning();
  int FillEvent( double val1 , double val2 , double mu1, double mu2 , const event_hhc& p , const nlo::amplitude_hhc& amp );
  double GetEcms();
  unsigned int GetNj();
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
  //     Reasonable choices for existing accelerators
  //		s =     40000.; // RHIC               200 GeV
  //		s =   3240000.; // TeV Run I         1800 GeV
  //		s =   3841600.; // TeV Run II        1960 GeV
  //		s =    810000.; // LHC Injection Run  900 GeV
  //		s =   5569600.; // LHC Initial Run   2360 GeV
  //		s =  49000000.; // LHC First Run     7000 GeV
  //		s = 100000000.; // LHC Start-up Run 10000 GeV
  //		s = 196000000.; // LHC Design Run   14000 GeV
  s =  49000000.;

  //   You can use your own phase generator. 
  //   Here we use the default.
  ps = 0;
} 

void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
  // --- fastNLO user:
  //     Here is your playground where you compute your observable 
  //     and define the two possible scale variables.
  //     At least one scale variable must be in 'GeV'. All observables
  //     must be in same dimension as your bingrid is defined.

  // --- run the jet algorithm
  double jetsize = 0.6;
  bounded_vector<lorentzvector<double> > pj = jetclus(p,jetsize);
  unsigned int nj = pj.upper(); 

  // --- declare and initialize phase space cut variables
  double pTmin = 20., ymin = 0., ymax = 4.4;
	
  // --- for inclusive jet cross section: loop over all jets
  vector<lorentzvector<double> > pj_ps;
  for(unsigned int i = 1; i <= nj; i++){
    double pt = pj[i].perp(); 
    double rap = fabs(pj[i].rapidity());
    // --- jet in phase space?
    if(pt > pTmin  &&  rap < ymax  &&  rap > ymin) { 
      pj_ps.push_back(pj[i]);
    }
  }
   
  // sort pj_ps array in pt
  vector<lorentzvector<double> > pj_ps_sort;
  lorentzvector<double> temp;
  for(unsigned int i = 0; i < pj_ps.size(); i++) {
    for (unsigned int y=0; y < (pj_ps.size()-i-1); y++){
      if (pj_ps[y].perp() < pj_ps[y+1].perp() ){        
	temp=pj_ps[y];
	pj_ps[y]=pj_ps[y+1];
	pj_ps[y+1]=temp;
      }
    }          
  }

  // identify dijet event
  // identify leading and subleading jet
  lorentzvector<double>* jet1 = NULL;
  lorentzvector<double>* jet2 = NULL;
  if ( pj_ps.size() > 1 ) {
    if ( pj_ps[0].perp() > 30 ) jet1 = &pj_ps[0];
    else return; // leading jet does not fulfills pt max criteria
    if ( pj_ps[1].perp() > 20 ) jet2 = &pj_ps[1];
    else return; // subleading jet does not fulfill pt_2 criteria
  }
  else return; // this is not a dijet event

  // calculate observables
  // this publication is published in TeV
  double ypsstar 	= fabs(jet1->rapidity()-jet2->rapidity())/2.;
  //  double ptmean_TeV	= ( jet1->perp() + jet2->perp() ) /2./1000.;
  //  double ptsum_TeV 	= ( jet1->perp() + jet2->perp() ) /2./1000.;
  double M12_TeV	= ((*jet1)+(*jet2)).mag()/1000.;
  double ptmax_TeV	= jet1->perp()/1000.;
  double pt2_TeV	= jet2->perp()/1000.;
   
  // --- fill fastNLO arrays
  //     Scales must be in GeV or dimensionless (at least one in GeV)
  //     Values must be in same dimension as your binning is defined
  FillEvent( M12_TeV, ypsstar , ptmax_TeV*1000. , pt2_TeV*1000. , p , amp );
} // --- end userfunc

int UserHHC::FillEvent( double val1 , double val2 , double mu1, double mu2 , const event_hhc& p, const nlo::amplitude_hhc& amp ){
  // ---- Fill FastNLO Array ---- //
  // --- fastNLO user: usually nothing to do

  // ---- get x-values ---- //
  double x1 = p[-1].Z()/p[hadron(-1)].Z();
  double x2 = p[0].Z()/p[hadron(0)].Z();   

  // --- fastNLO user: If this is not a single or double differential binning
  //     user has to perform calculation of bin number 'obsbin' by himself.
  int obsbin = table->GetBlockA2()->GetBinNumber( val1 , val2  );
  if (obsbin >= 0) {
    double prefactor = 1./table->GetBlockA2()->BinSize[obsbin]; 
    for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
      ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHCMuVar(obsbin,x1,x2,mu1,mu2,amp,dummypdf,pdf_pp,prefactor);// scales in GeV!
    }
  } // - end: fill fastNLO array
  return 0;
}

// ---- some function that can be used for reference calculations
double Fct_x(double x, double y){ return x; }
double Fct_y(double x, double y){ return y; }
double Fct_xyov2(double x, double y){ return (x+y)/2.; }
double Fct_x2y2ov2(double x, double y){ return sqrt((x*x+y*y)/2.); }
double Fct_x_exp03y(double x, double y){ return x*exp(0.3*y); }

void UserHHC::inittable(){
  // --- fastNLO user: This is the part, where the fastNLO table
  //     and the main initializations are set up. Please refer to the
  //     the documentation which parts have to be changed.
  //
 
  // ---- set up fastNLO and fnloTable---- //
  // --- fastNLO user: nothing to do 
  table		= new fnloTable(tablefilename);
  fnloBlockA1  *A1	= table->GetBlockA1();
  fnloBlockA2  *A2	= table->GetBlockA2();
  fnloBlockBNlojet *B	= new fnloBlockBNlojet(A1,A2);
  table->CreateBlockB(0,B);
  A1->SetContributionHeader();
  A2->Ecms	= GetEcms();
  A2->ILOord	= GetNj();
   
  // ---- set scenario name and description ---- //
  // --- fastNLO user: set scenario name (no whitespaces!)
  A1->SetScenName("fnl5002");

  // --- fastNLO user: up to 20 strings and any number of lines to describe the scenario
  A2->ScDescript.push_back("d2sigma-dijet_dm12dy*_[pb_TeV]");
  A2->ScDescript.push_back("ATLAS_Collaboration");
  A2->ScDescript.push_back("Dijet_Mass");
  A2->ScDescript.push_back("anti-kT_R=0.6");
  A2->ScDescript.push_back("arXiv:1112.6297");

  // --- fastNLO user: Give information about your measurement
  A2->SetIpublunits( 12 );			// --- fastNLO user: set cross section units (negative power of ten), e.g. 'pb' -> 12.
  A2->SetNumDiffBin( 2 );			// --- fastNLO user: No of dimensions in which observable is binned
  bool IsDiffBin = true;			// --- fastNLO user: Are publication units divided by this variable?
  A2->SetDimLabel("m12_[TeV]", 1 , IsDiffBin );// --- fastNLO user: label of 1st dimension
  A2->SetDimLabel("y*",        2 , IsDiffBin );// --- fastNLO user: label of 2nd dimension

  // ---- Define your bingrid in method DefineBinning() ---- //
  // --- fastNLO user: Modifiy function DefineBinning() below according to your bin grid.
  DefineBinning();

  // ---- initialize variables for WarmUp run ---- //
  // --- fastNLO user: Start "Warm-Up" or "Production" run.
  //     See documentation or GetWarmupValues() for more details.
  //     choices for B->SetDoWarmUp((bool))
  //	    -  B->SetDoWarmUp(true)   ->  Do the Warm-Up run
  //	    -  B->SetDoWarmUp(false)  ->  Do a production run
  B->SetDoWarmUp(true);

  // ---- get warm up values or init arrays reasonably ---- //
  // --- fastNLO user: See documentation in GetWarmUpValues for details.
  //     Do not modify this call.
  GetWarmupValues( B );

  // ---- initalize BlockB ---- //
  // --- fastNLO user: nothing to do.
  B->InitLHCConstants(A2,nlo);

  // ---- set number-of x-nodes ---- //
  B->SetNumberOfXNodesPerMagnitude( 8 , xlim );

  // ---- number of scale nodes for mu ---- //
  B->SetNumberOfScaleNodesScale1( 6 );
  B->SetNumberOfScaleNodesScale2( 5 );

  // ---- set names for the two possible scale variables (according to variables used in FillEvent()) ---- //
  B->SetScale1Name( "pT_max_[GeV]" );
  B->SetScale2Name( "pT_2_[GeV]" );

  // ---- Choose function for ScaleNode distances ---- //
  // --- fastNLO user: possibility to choose function which
  //     is used for binning of scale nodes
  //     possible choices
  //        - LogLog Binning ( H(mu) = log(log(mu/0.25)) )
  //        - Linear         ( H(mu) = mu )
  B->InitLogLogScaleNode( A2 , scale1lo , scale1hi , 1 );	// choose function H for scale 1
  B->InitLogLogScaleNode( A2 , scale2lo , scale2hi , 2 );	// choose function H for scale 2

  // ---- Initialize the cross section tables ---- //
  // --- fastNLO user: nothing to do.
  B->ResizeSigmaTildeTables( A2 );

  // ---- Reference tables ---- //
  // --- fastNLO user: You have the possibility to store
  //     three reference calculations in your fastNLO table
  //     which hold the plain nlojet++ results (cteq6m, as=0.1179, etc...).
  //     You can use those three references e.g. for three different scale 
  //     definitions.
  //     Therefore you can specify the functions how to calculate the renormalization
  //     and factorization scale using the scales mu1 and mu2 that you pass to 
  //     FastNLO when calling FillEventHHCMuVar();
  //     Usually you are following this convention:
  //      - 0		-> some 'mixed' or composed scale of scale values 1 and 2
  //      - 1		-> just scale 1 (if this is reasonably defined)
  //      - 2		-> just scale 2 (if this is reasonably defined)
  //
  //     INFO: Don't make any call, if your value does not make any sense (e.g. mu=|y|)
  //     INFO: The calculation of each reference table implies a recalculation of the matrix
  //     elements and slows down your calculation. 
  //     Each reference table needs around +25% computation time
  //
  B->SetFuncMuForReference( Fct_x_exp03y , Fct_x_exp03y , 0 );
  //B->SetFuncMuForReference( Fct_x , Fct_x , 1 );
  //B->SetFuncMuForReference( Fct_y , Fct_y , 2 );
}

void UserHHC::DefineBinning(){
  // ---- define the binning ---- //
  // --- fastNLO user: Define your binning here.
  //     If you use a double differential binning
  //     you should use the method fnloBlockA2::InitBinning().
  //     Therefore you have to define following variables.
  //     1) (int)		number of bins in your first dimension
  //     2) (double*)		array of bin edges in 1st dimension
  //     3) (int*)		array of number of bins for second dimension for each 1st dim bin
  //     4) (vector<double*>)	vector of arrays for bin edges in 2nd dimension for each 1st dim bin

  // ---- initalize the bingrids ---- //
  // --- fastNLO user: bin definitions - here in m12 and y*
  // 1)
  const int ndim2bins = 9;
  // 2)
  double dim2bins[ndim2bins+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 , 3.5 , 4.0, 4.4 };
  // 3)
  const int ndim1bins[ndim2bins] = {20,20,21,20,19,16,12,9,2};
   
  // prepare 4)
  double dim1bins[ndim2bins][22] = {
    { 0.07 , 0.11 , 0.16 , 0.21 , 0.26 , 0.31 , 0.37 , 0.44 , 0.51 , 0.59 , 
      0.67 , 0.76 , 0.85 , 0.95 , 1.06 , 1.18 , 1.31 , 1.45 , 1.60 , 
      1.94 , 2.78 },
    { 0.11 , 0.16 , 0.21 , 0.26 , 0.31 , 0.37 , 0.44 , 0.51 , 0.59 , 
      0.67 , 0.76 , 0.85 , 0.95 , 1.06 , 1.18 , 1.31 , 1.45 , 1.60 , 
      1.76 , 2.12, 2.31 },
    { 0.16 , 0.21 , 0.26 , 0.31 , 0.37 , 0.44 , 0.51 , 0.59 , 
      0.67 , 0.76 , 0.85 , 0.95 , 1.06 , 1.18 , 1.31 , 1.45 , 1.60 , 
      1.76 , 1.94, 2.12, 2.55, 3.61 },
    { 0.26 , 0.31 , 0.37 , 0.44 , 0.51 , 0.59 , 
      0.67 , 0.76 , 0.85 , 0.95 , 1.06 , 1.18 , 1.31 , 1.45 , 1.60 , 
      1.76 , 1.94, 2.12, 2.33, 2.78 , 3.93 },
    { 0.37 , 0.44 , 0.51 , 0.59 , 
      0.67 , 0.76 , 0.85 , 0.95 , 1.06 , 1.18 , 1.31 , 1.45 , 1.60 , 
      1.76 , 1.94, 2.12, 2.33, 2.55 , 3.04 , 4.27 },
    { 0.67 , 0.76 , 0.85 , 0.95 , 1.06 , 1.18 , 1.31 , 1.45 , 1.60 , 
      1.76 , 1.94, 2.12, 2.33, 2.55 , 2.78, 3.31 , 4.64 },
    { 1.18 , 1.31 , 1.45 , 1.60 , 
      1.76 , 1.94, 2.12, 2.33, 2.55 , 2.78, 3.04, 3.61 , 5.04 },
    { 1.76 , 1.94, 2.12, 2.33, 2.55 , 2.78, 3.04, 3.31 , 3.93 , 5.47 },
    { 2.55 , 3.04, 4.27 }
  };

  // 4)
  vector<double*> vdim1bins(ndim2bins);
  for (unsigned int i=0;i<vdim1bins.size();i++) vdim1bins[i]=dim1bins[i]; 
   
  // ---- pass arrays to FnloTable and init bingrid ---- //
  table->GetBlockA2()->InitBinning( ndim2bins , dim2bins , ndim1bins , vdim1bins );
}

void UserHHC::GetWarmupValues( fnloBlockBNlojet* B ){
  // --- fastNLO user: before running a new scenario for the first time,
  //     the first run must be a "Warm-Up Run" (B->IWarmUp=1), which produces
  //     an initialization block as output. This should be copied and pasted
  //     below into GetWarmupValues(). These initialization values must be used for all
  //     production jobs (IWarmUp=0) for a given scenario.
  // --- fastNLO user: Copy the values from the Warm-Up run here
  //     If this IS a Warm-Up run, we initialize the
  //     arrays with reasonable numbers.
  // --- Info: During a warm-up run, you might have to decrease
  //     your number of x-nodes per magnitude, that you do not
  //     run out of memory.

  // ---- Allocate Warmup run arrays ---- //
  const int NObsBin = table->GetBlockA2()->NObsBin;
  xlim = new double[NObsBin];
  scale1lo = new double[NObsBin];
  scale1hi = new double[NObsBin];
  scale2lo = new double[NObsBin];
  scale2hi = new double[NObsBin];

  // Get Warmup Values or initialize arrays with reasonable numbers
  if ( !B->GetDoWarmUp() ) {
    // ---- copy result from warmup run here ---- //
    // --------- fastNLO: Warm-Up run results (end)
  }
  else {
    printf("fastNLO: This is a warm-up run!\n");
    // --- fastNLO user: You can set the number of contributions
    //     after which the WarmUp values are printed
    B->SetWarmUpPrint(10000000);		// default 10000000

    // safe initialziations
    for(int i=0;i<NObsBin;i++){ 
      xlim[i] = 1.1e-07;
      scale1lo[i] = 1.0, scale1hi[i]=9.9e10;
      scale2lo[i] = 1.0, scale2hi[i]=9.9e10;
    }
  }
}
//------ END OF USER DEFINED PARTS, NO USER EDITS BELOW ------

//------ DON'T TOUCH THIS PART! ------
void UserHHC::initfunc(unsigned int)
{
  // --- Initialize event counters
  nevents = 0;
  // Set some defaults
  if (nwrite==0) nwrite = 5000000;
  start_time = std::time(0);
}

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

  nwrite = __save;
  inittable();
}

unsigned int UserHHC::GetNj(){
  unsigned int nj = 0, nu = 0 ,nd = 0;
  inputfunc(nj,nu,nd);
  return nj;
}

double UserHHC::GetEcms(){
  double ecms = 0;
  psinput(NULL,ecms);
  return sqrt(ecms);
}
