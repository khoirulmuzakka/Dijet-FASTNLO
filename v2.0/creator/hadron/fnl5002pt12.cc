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
  B->SetDoWarmUp(false);

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
    // 7191000000 contributions (!= events) in warm-up run
    xlim[ 0 ] = 1.25e-04 , scale1lo[ 0 ] = 30.6292 , scale1hi[ 0 ] = 107.4941 , scale2lo[ 0 ] = 20.0004 , scale2hi[ 0 ] = 65.7875;
    xlim[ 1 ] = 2.47e-04 , scale1lo[ 1 ] = 47.7807 , scale1hi[ 1 ] = 128.6247 , scale2lo[ 1 ] = 24.7222 , scale2hi[ 1 ] = 92.0305;
    xlim[ 2 ] = 5.23e-04 , scale1lo[ 2 ] = 69.6945 , scale1hi[ 2 ] = 151.5707 , scale2lo[ 2 ] = 44.0636 , scale2hi[ 2 ] = 120.7712;
    xlim[ 3 ] = 9.01e-04 , scale1lo[ 3 ] = 91.2532 , scale1hi[ 3 ] = 183.5435 , scale2lo[ 3 ] = 62.6111 , scale2hi[ 3 ] = 149.0746;
    xlim[ 4 ] = 1.38e-03 , scale1lo[ 4 ] = 113.0190 , scale1hi[ 4 ] = 218.9905 , scale2lo[ 4 ] = 81.7582 , scale2hi[ 4 ] = 177.9366;
    xlim[ 5 ] = 1.96e-03 , scale1lo[ 5 ] = 134.7150 , scale1hi[ 5 ] = 261.2707 , scale2lo[ 5 ] = 97.6033 , scale2hi[ 5 ] = 212.8949;
    xlim[ 6 ] = 2.80e-03 , scale1lo[ 6 ] = 160.7031 , scale1hi[ 6 ] = 310.9424 , scale2lo[ 6 ] = 116.5945 , scale2hi[ 6 ] = 252.6097;
    xlim[ 7 ] = 3.95e-03 , scale1lo[ 7 ] = 191.1317 , scale1hi[ 7 ] = 360.2873 , scale2lo[ 7 ] = 138.3406 , scale2hi[ 7 ] = 293.3285;
    xlim[ 8 ] = 5.32e-03 , scale1lo[ 8 ] = 221.4719 , scale1hi[ 8 ] = 416.7222 , scale2lo[ 8 ] = 160.4221 , scale2hi[ 8 ] = 338.9949;
    xlim[ 9 ] = 7.11e-03 , scale1lo[ 9 ] = 255.8731 , scale1hi[ 9 ] = 473.3186 , scale2lo[ 9 ] = 185.7083 , scale2hi[ 9 ] = 384.9539;
    xlim[ 10 ] = 9.17e-03 , scale1lo[ 10 ] = 291.3263 , scale1hi[ 10 ] = 536.7834 , scale2lo[ 10 ] = 210.7437 , scale2hi[ 10 ] = 436.7424;
    xlim[ 11 ] = 1.18e-02 , scale1lo[ 11 ] = 330.2440 , scale1hi[ 11 ] = 600.7053 , scale2lo[ 11 ] = 239.5922 , scale2hi[ 11 ] = 488.1688;
    xlim[ 12 ] = 1.47e-02 , scale1lo[ 12 ] = 369.0773 , scale1hi[ 12 ] = 670.7595 , scale2lo[ 12 ] = 267.7572 , scale2hi[ 12 ] = 545.7452;
    xlim[ 13 ] = 1.84e-02 , scale1lo[ 13 ] = 412.3899 , scale1hi[ 13 ] = 748.3083 , scale2lo[ 13 ] = 298.4413 , scale2hi[ 13 ] = 608.5820;
    xlim[ 14 ] = 2.29e-02 , scale1lo[ 14 ] = 460.0579 , scale1hi[ 14 ] = 833.3715 , scale2lo[ 14 ] = 333.4512 , scale2hi[ 14 ] = 676.4032;
    xlim[ 15 ] = 2.84e-02 , scale1lo[ 15 ] = 512.9704 , scale1hi[ 15 ] = 925.5150 , scale2lo[ 15 ] = 371.2446 , scale2hi[ 15 ] = 753.7032;
    xlim[ 16 ] = 3.50e-02 , scale1lo[ 16 ] = 569.8255 , scale1hi[ 16 ] = 1022.8729 , scale2lo[ 16 ] = 411.9433 , scale2hi[ 16 ] = 832.7206;
    xlim[ 17 ] = 4.29e-02 , scale1lo[ 17 ] = 630.1433 , scale1hi[ 17 ] = 1130.5004 , scale2lo[ 17 ] = 457.2017 , scale2hi[ 17 ] = 919.8861;
    xlim[ 18 ] = 5.22e-02 , scale1lo[ 18 ] = 695.5167 , scale1hi[ 18 ] = 1367.6611 , scale2lo[ 18 ] = 504.4842 , scale2hi[ 18 ] = 1111.8258;
    xlim[ 19 ] = 7.68e-02 , scale1lo[ 19 ] = 842.1419 , scale1hi[ 19 ] = 1959.6845 , scale2lo[ 19 ] = 610.3682 , scale2hi[ 19 ] = 1590.0493;
    xlim[ 20 ] = 3.22e-04 , scale1lo[ 20 ] = 35.0184 , scale1hi[ 20 ] = 119.8452 , scale2lo[ 20 ] = 20.0226 , scale2hi[ 20 ] = 78.7104;
    xlim[ 21 ] = 5.24e-04 , scale1lo[ 21 ] = 50.8418 , scale1hi[ 21 ] = 137.7873 , scale2lo[ 21 ] = 28.3372 , scale2hi[ 21 ] = 102.3712;
    xlim[ 22 ] = 9.02e-04 , scale1lo[ 22 ] = 66.6074 , scale1hi[ 22 ] = 162.5464 , scale2lo[ 22 ] = 39.6320 , scale2hi[ 22 ] = 127.5031;
    xlim[ 23 ] = 1.38e-03 , scale1lo[ 23 ] = 82.7331 , scale1hi[ 23 ] = 193.4670 , scale2lo[ 23 ] = 56.5489 , scale2hi[ 23 ] = 151.8130;
    xlim[ 24 ] = 1.96e-03 , scale1lo[ 24 ] = 98.5833 , scale1hi[ 24 ] = 230.8730 , scale2lo[ 24 ] = 71.2580 , scale2hi[ 24 ] = 180.8726;
    xlim[ 25 ] = 2.80e-03 , scale1lo[ 25 ] = 117.7968 , scale1hi[ 25 ] = 274.4858 , scale2lo[ 25 ] = 85.1112 , scale2hi[ 25 ] = 216.2073;
    xlim[ 26 ] = 3.95e-03 , scale1lo[ 26 ] = 139.9445 , scale1hi[ 26 ] = 318.0191 , scale2lo[ 26 ] = 101.2555 , scale2hi[ 26 ] = 251.8444;
    xlim[ 27 ] = 5.31e-03 , scale1lo[ 27 ] = 161.7571 , scale1hi[ 27 ] = 367.9378 , scale2lo[ 27 ] = 117.3571 , scale2hi[ 27 ] = 290.2984;
    xlim[ 28 ] = 7.11e-03 , scale1lo[ 28 ] = 187.8327 , scale1hi[ 28 ] = 418.0702 , scale2lo[ 28 ] = 135.6106 , scale2hi[ 28 ] = 328.6448;
    xlim[ 29 ] = 9.18e-03 , scale1lo[ 29 ] = 213.1781 , scale1hi[ 29 ] = 474.8127 , scale2lo[ 29 ] = 154.4146 , scale2hi[ 29 ] = 371.5045;
    xlim[ 30 ] = 1.18e-02 , scale1lo[ 30 ] = 241.3849 , scale1hi[ 30 ] = 530.8844 , scale2lo[ 30 ] = 175.2156 , scale2hi[ 30 ] = 416.2504;
    xlim[ 31 ] = 1.47e-02 , scale1lo[ 31 ] = 270.3760 , scale1hi[ 31 ] = 593.4190 , scale2lo[ 31 ] = 195.9036 , scale2hi[ 31 ] = 465.1766;
    xlim[ 32 ] = 1.84e-02 , scale1lo[ 32 ] = 301.5059 , scale1hi[ 32 ] = 660.7679 , scale2lo[ 32 ] = 218.7648 , scale2hi[ 32 ] = 516.2169;
    xlim[ 33 ] = 2.29e-02 , scale1lo[ 33 ] = 337.1121 , scale1hi[ 33 ] = 735.6167 , scale2lo[ 33 ] = 244.0908 , scale2hi[ 33 ] = 577.7795;
    xlim[ 34 ] = 2.84e-02 , scale1lo[ 34 ] = 374.9414 , scale1hi[ 34 ] = 818.2582 , scale2lo[ 34 ] = 272.7432 , scale2hi[ 34 ] = 644.2666;
    xlim[ 35 ] = 3.50e-02 , scale1lo[ 35 ] = 416.6090 , scale1hi[ 35 ] = 904.1803 , scale2lo[ 35 ] = 302.0921 , scale2hi[ 35 ] = 709.6743;
    xlim[ 36 ] = 4.29e-02 , scale1lo[ 36 ] = 460.2385 , scale1hi[ 36 ] = 997.8831 , scale2lo[ 36 ] = 334.7784 , scale2hi[ 36 ] = 783.5086;
    xlim[ 37 ] = 5.22e-02 , scale1lo[ 37 ] = 508.8055 , scale1hi[ 37 ] = 1098.9111 , scale2lo[ 37 ] = 368.5624 , scale2hi[ 37 ] = 861.0552;
    xlim[ 38 ] = 6.32e-02 , scale1lo[ 38 ] = 560.0374 , scale1hi[ 38 ] = 1324.5107 , scale2lo[ 38 ] = 405.1572 , scale2hi[ 38 ] = 1036.5741;
    xlim[ 39 ] = 9.17e-02 , scale1lo[ 39 ] = 674.8661 , scale1hi[ 39 ] = 1440.2839 , scale2lo[ 39 ] = 487.5124 , scale2hi[ 39 ] = 1129.3408;
    xlim[ 40 ] = 7.60e-04 , scale1lo[ 40 ] = 33.3891 , scale1hi[ 40 ] = 115.2252 , scale2lo[ 40 ] = 20.0151 , scale2hi[ 40 ] = 71.0531;
    xlim[ 41 ] = 1.02e-03 , scale1lo[ 41 ] = 43.9278 , scale1hi[ 41 ] = 128.2995 , scale2lo[ 41 ] = 22.3991 , scale2hi[ 41 ] = 88.7183;
    xlim[ 42 ] = 1.38e-03 , scale1lo[ 42 ] = 54.4293 , scale1hi[ 42 ] = 141.8565 , scale2lo[ 42 ] = 30.4315 , scale2hi[ 42 ] = 105.3376;
    xlim[ 43 ] = 1.97e-03 , scale1lo[ 43 ] = 64.7688 , scale1hi[ 43 ] = 168.6778 , scale2lo[ 43 ] = 40.9483 , scale2hi[ 43 ] = 125.5130;
    xlim[ 44 ] = 2.80e-03 , scale1lo[ 44 ] = 77.0442 , scale1hi[ 44 ] = 200.5063 , scale2lo[ 44 ] = 54.4413 , scale2hi[ 44 ] = 149.5184;
    xlim[ 45 ] = 3.96e-03 , scale1lo[ 45 ] = 91.9098 , scale1hi[ 45 ] = 232.6524 , scale2lo[ 45 ] = 65.5890 , scale2hi[ 45 ] = 172.9965;
    xlim[ 46 ] = 5.32e-03 , scale1lo[ 46 ] = 106.4705 , scale1hi[ 46 ] = 268.4398 , scale2lo[ 46 ] = 77.7428 , scale2hi[ 46 ] = 199.9505;
    xlim[ 47 ] = 7.11e-03 , scale1lo[ 47 ] = 123.1344 , scale1hi[ 47 ] = 304.6352 , scale2lo[ 47 ] = 89.5197 , scale2hi[ 47 ] = 227.2180;
    xlim[ 48 ] = 9.17e-03 , scale1lo[ 48 ] = 139.7014 , scale1hi[ 48 ] = 346.2147 , scale2lo[ 48 ] = 101.4882 , scale2hi[ 48 ] = 257.6066;
    xlim[ 49 ] = 1.18e-02 , scale1lo[ 49 ] = 158.5592 , scale1hi[ 49 ] = 387.8554 , scale2lo[ 49 ] = 114.6511 , scale2hi[ 49 ] = 288.5187;
    xlim[ 50 ] = 1.47e-02 , scale1lo[ 50 ] = 177.4103 , scale1hi[ 50 ] = 433.2824 , scale2lo[ 50 ] = 129.4098 , scale2hi[ 50 ] = 321.5270;
    xlim[ 51 ] = 1.84e-02 , scale1lo[ 51 ] = 198.8866 , scale1hi[ 51 ] = 482.2975 , scale2lo[ 51 ] = 143.9425 , scale2hi[ 51 ] = 358.6595;
    xlim[ 52 ] = 2.29e-02 , scale1lo[ 52 ] = 221.2182 , scale1hi[ 52 ] = 538.9812 , scale2lo[ 52 ] = 160.3445 , scale2hi[ 52 ] = 398.6204;
    xlim[ 53 ] = 2.84e-02 , scale1lo[ 53 ] = 246.2698 , scale1hi[ 53 ] = 596.8295 , scale2lo[ 53 ] = 179.0293 , scale2hi[ 53 ] = 443.3181;
    xlim[ 54 ] = 3.50e-02 , scale1lo[ 54 ] = 273.3261 , scale1hi[ 54 ] = 661.8227 , scale2lo[ 54 ] = 197.7605 , scale2hi[ 54 ] = 491.0150;
    xlim[ 55 ] = 4.29e-02 , scale1lo[ 55 ] = 303.3385 , scale1hi[ 55 ] = 727.6952 , scale2lo[ 55 ] = 219.9821 , scale2hi[ 55 ] = 547.0487;
    xlim[ 56 ] = 5.22e-02 , scale1lo[ 56 ] = 334.3821 , scale1hi[ 56 ] = 800.1678 , scale2lo[ 56 ] = 242.6635 , scale2hi[ 56 ] = 595.3741;
    xlim[ 57 ] = 6.32e-02 , scale1lo[ 57 ] = 367.3355 , scale1hi[ 57 ] = 883.8460 , scale2lo[ 57 ] = 267.9655 , scale2hi[ 57 ] = 659.9900;
    xlim[ 58 ] = 7.68e-02 , scale1lo[ 58 ] = 405.6508 , scale1hi[ 58 ] = 962.4906 , scale2lo[ 58 ] = 292.2931 , scale2hi[ 58 ] = 724.5024;
    xlim[ 59 ] = 9.18e-02 , scale1lo[ 59 ] = 443.4268 , scale1hi[ 59 ] = 1156.0837 , scale2lo[ 59 ] = 323.8234 , scale2hi[ 59 ] = 861.6057;
    xlim[ 60 ] = 1.33e-01 , scale1lo[ 60 ] = 533.9291 , scale1hi[ 60 ] = 1638.9824 , scale2lo[ 60 ] = 386.9538 , scale2hi[ 60 ] = 1215.4978;
    xlim[ 61 ] = 1.99e-03 , scale1lo[ 61 ] = 33.9884 , scale1hi[ 61 ] = 109.2662 , scale2lo[ 61 ] = 20.0030 , scale2hi[ 61 ] = 66.8374;
    xlim[ 62 ] = 2.46e-03 , scale1lo[ 62 ] = 40.6420 , scale1hi[ 62 ] = 120.5608 , scale2lo[ 62 ] = 20.1809 , scale2hi[ 62 ] = 79.8859;
    xlim[ 63 ] = 2.95e-03 , scale1lo[ 63 ] = 48.4041 , scale1hi[ 63 ] = 140.8524 , scale2lo[ 63 ] = 26.8088 , scale2hi[ 63 ] = 94.7875;
    xlim[ 64 ] = 3.96e-03 , scale1lo[ 64 ] = 57.6579 , scale1hi[ 64 ] = 154.6625 , scale2lo[ 64 ] = 35.5607 , scale2hi[ 64 ] = 109.7911;
    xlim[ 65 ] = 5.31e-03 , scale1lo[ 65 ] = 66.7962 , scale1hi[ 65 ] = 175.8877 , scale2lo[ 65 ] = 44.0164 , scale2hi[ 65 ] = 126.3224;
    xlim[ 66 ] = 7.12e-03 , scale1lo[ 66 ] = 77.3906 , scale1hi[ 66 ] = 200.2589 , scale2lo[ 66 ] = 49.7228 , scale2hi[ 66 ] = 145.0775;
    xlim[ 67 ] = 9.17e-03 , scale1lo[ 67 ] = 87.8601 , scale1hi[ 67 ] = 225.4865 , scale2lo[ 67 ] = 63.6631 , scale2hi[ 67 ] = 163.6433;
    xlim[ 68 ] = 1.18e-02 , scale1lo[ 68 ] = 99.3630 , scale1hi[ 68 ] = 253.6975 , scale2lo[ 68 ] = 72.0773 , scale2hi[ 68 ] = 182.5772;
    xlim[ 69 ] = 1.47e-02 , scale1lo[ 69 ] = 111.2753 , scale1hi[ 69 ] = 283.3898 , scale2lo[ 69 ] = 80.7677 , scale2hi[ 69 ] = 204.7236;
    xlim[ 70 ] = 1.84e-02 , scale1lo[ 70 ] = 124.0907 , scale1hi[ 70 ] = 315.7611 , scale2lo[ 70 ] = 90.0700 , scale2hi[ 70 ] = 229.8906;
    xlim[ 71 ] = 2.29e-02 , scale1lo[ 71 ] = 138.5918 , scale1hi[ 71 ] = 352.1690 , scale2lo[ 71 ] = 100.7254 , scale2hi[ 71 ] = 253.4699;
    xlim[ 72 ] = 2.84e-02 , scale1lo[ 72 ] = 154.3387 , scale1hi[ 72 ] = 391.0987 , scale2lo[ 72 ] = 111.6317 , scale2hi[ 72 ] = 281.9516;
    xlim[ 73 ] = 3.50e-02 , scale1lo[ 73 ] = 171.3767 , scale1hi[ 73 ] = 433.7890 , scale2lo[ 73 ] = 124.8933 , scale2hi[ 73 ] = 311.6508;
    xlim[ 74 ] = 4.29e-02 , scale1lo[ 74 ] = 190.2490 , scale1hi[ 74 ] = 478.5470 , scale2lo[ 74 ] = 137.9820 , scale2hi[ 74 ] = 346.2940;
    xlim[ 75 ] = 5.23e-02 , scale1lo[ 75 ] = 209.5253 , scale1hi[ 75 ] = 524.0719 , scale2lo[ 75 ] = 151.4654 , scale2hi[ 75 ] = 378.4444;
    xlim[ 76 ] = 6.32e-02 , scale1lo[ 76 ] = 230.8102 , scale1hi[ 76 ] = 579.8252 , scale2lo[ 76 ] = 166.6649 , scale2hi[ 76 ] = 417.2505;
    xlim[ 77 ] = 7.68e-02 , scale1lo[ 77 ] = 253.6047 , scale1hi[ 77 ] = 630.6524 , scale2lo[ 77 ] = 184.1118 , scale2hi[ 77 ] = 456.8418;
    xlim[ 78 ] = 9.17e-02 , scale1lo[ 78 ] = 277.3241 , scale1hi[ 78 ] = 695.8845 , scale2lo[ 78 ] = 200.5424 , scale2hi[ 78 ] = 500.1388;
    xlim[ 79 ] = 1.11e-01 , scale1lo[ 79 ] = 305.1244 , scale1hi[ 79 ] = 826.1241 , scale2lo[ 79 ] = 221.2376 , scale2hi[ 79 ] = 598.9389;
    xlim[ 80 ] = 1.58e-01 , scale1lo[ 80 ] = 365.4398 , scale1hi[ 80 ] = 1168.8924 , scale2lo[ 80 ] = 263.9692 , scale2hi[ 80 ] = 839.8596;
    xlim[ 81 ] = 4.25e-03 , scale1lo[ 81 ] = 30.0766 , scale1hi[ 81 ] = 105.0057 , scale2lo[ 81 ] = 20.0000 , scale2hi[ 81 ] = 58.5259;
    xlim[ 82 ] = 5.11e-03 , scale1lo[ 82 ] = 35.4535 , scale1hi[ 82 ] = 113.8241 , scale2lo[ 82 ] = 20.0016 , scale2hi[ 82 ] = 67.7775;
    xlim[ 83 ] = 6.43e-03 , scale1lo[ 83 ] = 41.2037 , scale1hi[ 83 ] = 122.8324 , scale2lo[ 83 ] = 21.7816 , scale2hi[ 83 ] = 78.7594;
    xlim[ 84 ] = 7.79e-03 , scale1lo[ 84 ] = 47.4054 , scale1hi[ 84 ] = 125.9932 , scale2lo[ 84 ] = 25.8360 , scale2hi[ 84 ] = 89.1032;
    xlim[ 85 ] = 9.19e-03 , scale1lo[ 85 ] = 53.9652 , scale1hi[ 85 ] = 141.4480 , scale2lo[ 85 ] = 33.3106 , scale2hi[ 85 ] = 101.1168;
    xlim[ 86 ] = 1.18e-02 , scale1lo[ 86 ] = 61.3582 , scale1hi[ 86 ] = 158.3457 , scale2lo[ 86 ] = 38.9012 , scale2hi[ 86 ] = 113.2342;
    xlim[ 87 ] = 1.47e-02 , scale1lo[ 87 ] = 68.6745 , scale1hi[ 87 ] = 175.7178 , scale2lo[ 87 ] = 45.4071 , scale2hi[ 87 ] = 126.2610;
    xlim[ 88 ] = 1.84e-02 , scale1lo[ 88 ] = 76.2568 , scale1hi[ 88 ] = 197.5776 , scale2lo[ 88 ] = 54.1979 , scale2hi[ 88 ] = 140.8701;
    xlim[ 89 ] = 2.29e-02 , scale1lo[ 89 ] = 85.4212 , scale1hi[ 89 ] = 219.5032 , scale2lo[ 89 ] = 61.9753 , scale2hi[ 89 ] = 156.7926;
    xlim[ 90 ] = 2.84e-02 , scale1lo[ 90 ] = 95.3465 , scale1hi[ 90 ] = 243.4506 , scale2lo[ 90 ] = 69.1329 , scale2hi[ 90 ] = 174.1224;
    xlim[ 91 ] = 3.50e-02 , scale1lo[ 91 ] = 105.4297 , scale1hi[ 91 ] = 269.4786 , scale2lo[ 91 ] = 76.8534 , scale2hi[ 91 ] = 192.7346;
    xlim[ 92 ] = 4.29e-02 , scale1lo[ 92 ] = 116.6983 , scale1hi[ 92 ] = 298.2567 , scale2lo[ 92 ] = 85.1828 , scale2hi[ 92 ] = 212.6139;
    xlim[ 93 ] = 5.23e-02 , scale1lo[ 93 ] = 128.3493 , scale1hi[ 93 ] = 327.5379 , scale2lo[ 93 ] = 93.9002 , scale2hi[ 93 ] = 233.8658;
    xlim[ 94 ] = 6.32e-02 , scale1lo[ 94 ] = 141.3837 , scale1hi[ 94 ] = 361.5783 , scale2lo[ 94 ] = 103.0720 , scale2hi[ 94 ] = 257.7973;
    xlim[ 95 ] = 7.68e-02 , scale1lo[ 95 ] = 156.7283 , scale1hi[ 95 ] = 395.1566 , scale2lo[ 95 ] = 113.1083 , scale2hi[ 95 ] = 281.7192;
    xlim[ 96 ] = 9.18e-02 , scale1lo[ 96 ] = 171.5501 , scale1hi[ 96 ] = 431.7203 , scale2lo[ 96 ] = 124.3362 , scale2hi[ 96 ] = 309.6231;
    xlim[ 97 ] = 1.11e-01 , scale1lo[ 97 ] = 187.9481 , scale1hi[ 97 ] = 474.4639 , scale2lo[ 97 ] = 136.1222 , scale2hi[ 97 ] = 338.8535;
    xlim[ 98 ] = 1.33e-01 , scale1lo[ 98 ] = 206.3626 , scale1hi[ 98 ] = 562.6249 , scale2lo[ 98 ] = 150.6966 , scale2hi[ 98 ] = 404.6575;
    xlim[ 99 ] = 1.89e-01 , scale1lo[ 99 ] = 244.3430 , scale1hi[ 99 ] = 788.0831 , scale2lo[ 99 ] = 177.0800 , scale2hi[ 99 ] = 567.4124;
    xlim[ 100 ] = 1.17e-02 , scale1lo[ 100 ] = 33.0132 , scale1hi[ 100 ] = 104.7018 , scale2lo[ 100 ] = 20.0023 , scale2hi[ 100 ] = 61.9619;
    xlim[ 101 ] = 1.36e-02 , scale1lo[ 101 ] = 37.6474 , scale1hi[ 101 ] = 110.3996 , scale2lo[ 101 ] = 20.1495 , scale2hi[ 101 ] = 69.2981;
    xlim[ 102 ] = 1.63e-02 , scale1lo[ 102 ] = 42.0423 , scale1hi[ 102 ] = 115.1245 , scale2lo[ 102 ] = 22.4523 , scale2hi[ 102 ] = 77.4556;
    xlim[ 103 ] = 1.99e-02 , scale1lo[ 103 ] = 46.8810 , scale1hi[ 103 ] = 125.7330 , scale2lo[ 103 ] = 26.6438 , scale2hi[ 103 ] = 86.4237;
    xlim[ 104 ] = 2.30e-02 , scale1lo[ 104 ] = 52.2494 , scale1hi[ 104 ] = 134.9770 , scale2lo[ 104 ] = 31.6954 , scale2hi[ 104 ] = 96.1934;
    xlim[ 105 ] = 2.85e-02 , scale1lo[ 105 ] = 58.1374 , scale1hi[ 105 ] = 148.8667 , scale2lo[ 105 ] = 35.6126 , scale2hi[ 105 ] = 106.7954;
    xlim[ 106 ] = 3.51e-02 , scale1lo[ 106 ] = 64.5483 , scale1hi[ 106 ] = 166.0786 , scale2lo[ 106 ] = 45.2121 , scale2hi[ 106 ] = 118.2059;
    xlim[ 107 ] = 4.29e-02 , scale1lo[ 107 ] = 71.0837 , scale1hi[ 107 ] = 182.3080 , scale2lo[ 107 ] = 51.2017 , scale2hi[ 107 ] = 130.4537;
    xlim[ 108 ] = 5.23e-02 , scale1lo[ 108 ] = 79.1065 , scale1hi[ 108 ] = 199.1874 , scale2lo[ 108 ] = 57.0685 , scale2hi[ 108 ] = 143.4683;
    xlim[ 109 ] = 6.33e-02 , scale1lo[ 109 ] = 86.4709 , scale1hi[ 109 ] = 221.2321 , scale2lo[ 109 ] = 62.9322 , scale2hi[ 109 ] = 158.1681;
    xlim[ 110 ] = 7.69e-02 , scale1lo[ 110 ] = 96.0591 , scale1hi[ 110 ] = 242.7160 , scale2lo[ 110 ] = 69.7092 , scale2hi[ 110 ] = 172.8207;
    xlim[ 111 ] = 9.18e-02 , scale1lo[ 111 ] = 104.2005 , scale1hi[ 111 ] = 263.9844 , scale2lo[ 111 ] = 76.2533 , scale2hi[ 111 ] = 189.9631;
    xlim[ 112 ] = 1.11e-01 , scale1lo[ 112 ] = 115.4926 , scale1hi[ 112 ] = 288.6457 , scale2lo[ 112 ] = 83.3108 , scale2hi[ 112 ] = 207.8632;
    xlim[ 113 ] = 1.33e-01 , scale1lo[ 113 ] = 126.1573 , scale1hi[ 113 ] = 316.5331 , scale2lo[ 113 ] = 91.8984 , scale2hi[ 113 ] = 226.6431;
    xlim[ 114 ] = 1.58e-01 , scale1lo[ 114 ] = 138.1136 , scale1hi[ 114 ] = 374.8793 , scale2lo[ 114 ] = 100.3970 , scale2hi[ 114 ] = 269.8449;
    xlim[ 115 ] = 2.24e-01 , scale1lo[ 115 ] = 162.5022 , scale1hi[ 115 ] = 530.4557 , scale2lo[ 115 ] = 119.5640 , scale2hi[ 115 ] = 378.2139;
    xlim[ 116 ] = 3.38e-02 , scale1lo[ 116 ] = 35.6130 , scale1hi[ 116 ] = 106.5675 , scale2lo[ 116 ] = 20.0320 , scale2hi[ 116 ] = 65.0521;
    xlim[ 117 ] = 3.93e-02 , scale1lo[ 117 ] = 39.5026 , scale1hi[ 117 ] = 110.0073 , scale2lo[ 117 ] = 21.0630 , scale2hi[ 117 ] = 71.9962;
    xlim[ 118 ] = 4.74e-02 , scale1lo[ 118 ] = 43.7529 , scale1hi[ 118 ] = 109.3576 , scale2lo[ 118 ] = 24.0936 , scale2hi[ 118 ] = 79.4536;
    xlim[ 119 ] = 5.58e-02 , scale1lo[ 119 ] = 48.2802 , scale1hi[ 119 ] = 120.1992 , scale2lo[ 119 ] = 30.3126 , scale2hi[ 119 ] = 87.4030;
    xlim[ 120 ] = 6.34e-02 , scale1lo[ 120 ] = 52.7437 , scale1hi[ 120 ] = 133.7682 , scale2lo[ 120 ] = 34.1101 , scale2hi[ 120 ] = 96.3283;
    xlim[ 121 ] = 7.69e-02 , scale1lo[ 121 ] = 57.8483 , scale1hi[ 121 ] = 145.3655 , scale2lo[ 121 ] = 41.3341 , scale2hi[ 121 ] = 105.2557;
    xlim[ 122 ] = 9.18e-02 , scale1lo[ 122 ] = 63.9689 , scale1hi[ 122 ] = 159.0486 , scale2lo[ 122 ] = 46.4886 , scale2hi[ 122 ] = 115.6876;
    xlim[ 123 ] = 1.11e-01 , scale1lo[ 123 ] = 69.2713 , scale1hi[ 123 ] = 176.8785 , scale2lo[ 123 ] = 51.3805 , scale2hi[ 123 ] = 126.6045;
    xlim[ 124 ] = 1.33e-01 , scale1lo[ 124 ] = 76.9545 , scale1hi[ 124 ] = 189.4953 , scale2lo[ 124 ] = 55.5122 , scale2hi[ 124 ] = 138.0402;
    xlim[ 125 ] = 1.58e-01 , scale1lo[ 125 ] = 82.7504 , scale1hi[ 125 ] = 209.1856 , scale2lo[ 125 ] = 61.5190 , scale2hi[ 125 ] = 150.9188;
    xlim[ 126 ] = 1.89e-01 , scale1lo[ 126 ] = 91.3963 , scale1hi[ 126 ] = 248.6236 , scale2lo[ 126 ] = 66.6694 , scale2hi[ 126 ] = 179.2475;
    xlim[ 127 ] = 2.66e-01 , scale1lo[ 127 ] = 108.2189 , scale1hi[ 127 ] = 344.4729 , scale2lo[ 127 ] = 80.1106 , scale2hi[ 127 ] = 250.2778;
    xlim[ 128 ] = 8.05e-02 , scale1lo[ 128 ] = 32.2653 , scale1hi[ 128 ] = 91.2924 , scale2lo[ 128 ] = 20.2037 , scale2hi[ 128 ] = 58.5172;
    xlim[ 129 ] = 8.89e-02 , scale1lo[ 129 ] = 35.5462 , scale1hi[ 129 ] = 96.1890 , scale2lo[ 129 ] = 21.3911 , scale2hi[ 129 ] = 63.9276;
    xlim[ 130 ] = 9.99e-02 , scale1lo[ 130 ] = 38.8431 , scale1hi[ 130 ] = 103.9790 , scale2lo[ 130 ] = 22.0052 , scale2hi[ 130 ] = 70.2849;
    xlim[ 131 ] = 1.19e-01 , scale1lo[ 131 ] = 42.6710 , scale1hi[ 131 ] = 106.1034 , scale2lo[ 131 ] = 24.8560 , scale2hi[ 131 ] = 76.8952;
    xlim[ 132 ] = 1.39e-01 , scale1lo[ 132 ] = 46.7662 , scale1hi[ 132 ] = 113.6188 , scale2lo[ 132 ] = 31.7028 , scale2hi[ 132 ] = 83.8383;
    xlim[ 133 ] = 1.62e-01 , scale1lo[ 133 ] = 50.9147 , scale1hi[ 133 ] = 126.2161 , scale2lo[ 133 ] = 36.5104 , scale2hi[ 133 ] = 91.6610;
    xlim[ 134 ] = 1.89e-01 , scale1lo[ 134 ] = 55.7265 , scale1hi[ 134 ] = 136.6608 , scale2lo[ 134 ] = 37.2985 , scale2hi[ 134 ] = 99.8362;
    xlim[ 135 ] = 2.24e-01 , scale1lo[ 135 ] = 60.6561 , scale1hi[ 135 ] = 161.8624 , scale2lo[ 135 ] = 43.8049 , scale2hi[ 135 ] = 118.4705;
    xlim[ 136 ] = 3.15e-01 , scale1lo[ 136 ] = 71.9784 , scale1hi[ 136 ] = 226.2755 , scale2lo[ 136 ] = 54.6285 , scale2hi[ 136 ] = 164.9555;
    xlim[ 137 ] = 2.01e-01 , scale1lo[ 137 ] = 31.4541 , scale1hi[ 137 ] = 87.5030 , scale2lo[ 137 ] = 21.3662 , scale2hi[ 137 ] = 55.6051;
    xlim[ 138 ] = 2.43e-01 , scale1lo[ 138 ] = 37.9726 , scale1hi[ 138 ] = 103.2901 , scale2lo[ 138 ] = 25.8863 , scale2hi[ 138 ] = 78.0877;
    // --------- fastNLO: Warm-Up run results (end)
  }
  else {
    printf("fastNLO: This is a warm-up run!\n");
    // --- fastNLO user: You can set the number of contributions
    //     after which the WarmUp values are printed
    B->SetWarmUpPrint(1000000);		// default 10000000

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
