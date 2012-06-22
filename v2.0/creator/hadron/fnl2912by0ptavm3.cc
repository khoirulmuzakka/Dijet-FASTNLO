//
// fastNLO v2.1 author code for fnl2912b:
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
  double* xlim;
  double* scale1lo;
  double* scale1hi;
  double* scale2lo;
  double* scale2hi;
   
  // --- fastNLO user FYI: steering of these flags is encoded in NLOJet++ run name (option -n name)
  //     if filename matches "deb", "ref", or "wrm" the respective flag is set to true
  bool doDebug, doReference, doWarmUp;
  bool nlo;

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
  //nj = 2U;    
  nj = 3U;

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
  s =  49000000.; // LHC First Run     7000 GeV
  //		s = 100000000.; // LHC Start-up Run 10000 GeV
  //		s = 196000000.; // LHC Design Run   14000 GeV

  //   You can use your own phase generator. 
  //   Here we use the default.
  ps = 0;
} 

void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
  // --- fastNLO: Don't touch this piece of code!
  fnloBlockA2 *A2 =  table->GetBlockA2();

  // --- fastNLO user: Set jet size and run the jet algorithm
  double jetsize = 0.5;
  bounded_vector<lorentzvector<double> > pj = jetclus(p,jetsize);
  unsigned int nj = pj.upper(); 

  // --- give some debug output before selection and sorting
  if ( doDebug ) {
    for (unsigned int i=1; i<=nj; i++) {
      double pti  = pj[i].perp();
      double yi   = pj[i].rapidity();
      double etai = pj[i].prapidity();
      cout << "before cuts: jet # i, pt, y, eta: " << i << ", " << pti << ", " << yi << ", " << etai << endl;
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
  //     Usually, pT and E are in GeV, but this may be changed.
  //     ATTENTION: Scales must always be in GeV!

  // --- declare and initialize phase space cut variables
  const double pTmin = 50., ymin = 0., ymax = 3.0;
	
  // --- loop over all jets and make selection
  vector<lorentzvector<double> > pj_ps;
  for(unsigned int i = 1; i <= nj; i++){
    double pt = pj[i].perp(); 
    double rap = fabs(pj[i].rapidity());
    // --- jet in phase space?
    if(pTmin <= pt && ymin <= rap && rap < ymax) { 
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

  // --- give some debug output after selection and sorting
  unsigned int njet = pj_ps.size();
  if ( doDebug ) {
    cout << "# jets before and after phase space cuts: nj, njet = " << nj << ", " << njet << endl;
    cout << "phase space cuts: yjmin, yjmax, ptjmin: " << ymin << ", " << ymax << ", " << pTmin << endl;
    for (unsigned int i=0; i<njet; i++) {
      double pti  = pj_ps[i].perp();
      double yi   = pj_ps[i].rapidity();
      double etai = pj_ps[i].prapidity();
      cout << "after cuts: jet # i, pt, y, eta: " << i << ", " << pti << ", " << yi << ", " << etai << endl;
    }
  }

  // 3-jet mass requires at least 3 jets
  if (njet > 2) {
    
    // --- declare and initialize additional cut variables
    // minimum pT required for 3rd jet
    const double pT3min = 100.;
    // minimal 3-jet mass for events to be considered
    double m3jmin = A2->LoBin[0][0];
    // maximal |y| of three leading jets
    const double y3jmax = 1.0;
    
    // Derive 3-jet variables
    // 3-jet mass
    lorentzvector<double> pj123 = pj_ps[0] + pj_ps[1] + pj_ps[2]; 
    double m3j = pj123.mag();
    if (m3j < 0.) {cout << "Warning!: Negative Mass" << m3j << endl;}
    
    // Maximal (pseudo-)rapidity
    double y3j = max(abs(pj_ps[0].rapidity()),abs(pj_ps[1].rapidity()));
    y3j = max(y3j,abs(pj_ps[2].rapidity()));
    
    // Average pTs of leading jets
    //    double pT1   = (pj_ps[0].perp()) / 1.0;
    //    double pT12  = (pj_ps[0].perp() + pj_ps[1].perp()) / 2.0;
    double pT123 = (pj_ps[0].perp() + pj_ps[1].perp() + pj_ps[2].perp()) / 3.0;
    
    // pT3
    double pT3 = pj_ps[2].perp();

    // --- Further 3-jet phase space cuts?
    if ( m3jmin <= m3j && y3j < y3jmax && pT3min <= pT3 ) {

      // --- fill fastNLO arrays
      //     Scales must be in GeV or dimensionless (at least one in GeV)
      //     Values must be in same dimension as your binning is defined
      FillEvent( m3j, y3j, pT123, m3j/2., p, amp );
    }
  }
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

// ---- some functions that can be used for reference calculations
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
  A1->SetScenName("fnl2912by0ptavm3");

  // --- fastNLO user: up to 20 strings and any number of lines to describe the scenario
  A2->ScDescript.push_back("d2sigma-3-jet_dM3jd|y3j|_max_[pb_GeV]");
  A2->ScDescript.push_back("CMS_Collaboration");
  A2->ScDescript.push_back("E_cms=7_TeV");
  A2->ScDescript.push_back("3-Jet_Mass");
  A2->ScDescript.push_back("anti-kT_R=0.5");
  A2->ScDescript.push_back("CMS-PAS-QCD-11-003");
  A2->ScDescript.push_back("provided by:");
  A2->ScDescript.push_back("fastNLO_2.1.0");
  A2->ScDescript.push_back("If you use this table, please cite:");
  A2->ScDescript.push_back("  D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch, arXiv:1109.1310");

  // --- fastNLO user: Give information about your measurement
  A2->SetIpublunits( 12 );			// --- fastNLO user: set cross section units (negative power of ten), e.g. 'pb' -> 12.
  A2->SetNumDiffBin( 2 );			// --- fastNLO user: No of dimensions in which observable is binned
  bool IsDiffBin = true;			// --- fastNLO user: Are publication units divided by this variable?
  A2->SetDimLabel("M3j_[GeV]", 1 , IsDiffBin );// --- fastNLO user: label of 1st dimension
  A2->SetDimLabel("|y3j|_max", 2 , IsDiffBin );// --- fastNLO user: label of 2nd dimension

  // ---- Define your bingrid in method DefineBinning() ---- //
  // --- fastNLO user: Modifiy function DefineBinning() below according to your bin grid.
  DefineBinning();

  // ---- initialize variables for WarmUp run ---- //
  // KR: Now written via filename string match to "wrm" into bool variable doWarmUp 
  // --- fastNLO user: Start "Warm-Up" or "Production" run.
  //     See documentation or GetWarmupValues() for more details.
  //     choices for B->SetDoWarmUp((bool))
  //	    -  B->SetDoWarmUp(true)   ->  Do the Warm-Up run
  //	    -  B->SetDoWarmUp(false)  ->  Do a production run
  B->SetDoWarmUp(doWarmUp);
  
  // --- fastNLO user: You can set the number of contributions
  //     after which the WarmUp values are printed
  B->SetWarmUpPrint(1000000);		// default 10000000

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
  B->SetScale1Name( "<pT_1,2,3>_[GeV]" );
  B->SetScale2Name( "M3j_[GeV]/2" );

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
  //     ATTENTION: Please make sure the scales resp. functions to derive scales
  //                actually make sense for your definitions!
  //                Otherwise you might find "NaN"'s in your reference tables
  //                rendering them unusable in the fastNLO merge step.  
  //
  if(doReference){
    //B->SetFuncMuForReference( Fct_x_exp03y , Fct_x_exp03y , 0 );
    //B->SetFuncMuForReference( Fct_x , Fct_x , 1 );
    B->SetFuncMuForReference( Fct_y , Fct_y , 2 );
  }
}

void UserHHC::DefineBinning(){
  // --- fastNLO: Don't touch this piece of code!
  fnloBlockA2 *A2 =  table->GetBlockA2();

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
  // --- fastNLO user: bin definitions - here in M3j and |y3j|_max
  const int ndim2bins = 1;
  const double dim2bins[ndim2bins+1] = {0.0, 1.0};
  
  const int ndim1bins[ndim2bins] = { 37 };
  
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
  const double dim0[38] = {
    313.,   354.,  398.,  445.,  495.,  548.,  604.,  664.,  727.,  794.,
    864.,   938., 1016., 1098., 1184., 1274., 1369., 1469., 1573., 1682.,
    1796., 1916., 2041., 2172., 2309., 2452., 2602., 2758., 2921., 3092.,
    3270., 3456., 3650., 3852., 4063., 4283., 4513., 4753. };
  for (int i=0; i<ndim2bins; i++) {
    for (int j=0; j<ndim1bins[i]+1; j++) { 
      dim1bins[i][j] = dim0[j];
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
  //     define the bin width ("binsize") by which
  //     the cross section is divided to obtain the 
  //     (multi-) differential result.
  //     default: divide by bin width in dim 1 and dim 2
  //              ATTENTION: Don't forget to include a factor of 2 for abs. rapidity |y| !
  // fnl2912: divide by bin width in M3j and |y3j|_max
  // ---- pass arrays to FnloTable and init bingrid ---- //
  const double bwfactor = 2.;
  A2->InitBinningKR( ndim2bins , dim2bins , ndim1bins , dim1bins, bwfactor );
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
  if ( B->GetDoWarmUp() ) {
    cout << endl << "fastNLO: Initializing x limits for warm-up run ..." << endl;
    // safe initialziations
    for(int i=0;i<NObsBin;i++){ 
      xlim[i] = 1.1e-07;
      scale1lo[i] = 1.0, scale1hi[i]=9.9e10;
      scale2lo[i] = 1.0, scale2hi[i]=9.9e10;
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
      // --------- fastNLO: Warm-Up run results (end)
      // Safety check:
      // Count no. of xlim values > 0 (should be equal to NObsBin!)   
      int nxlim = 0;
      for (int i=0; i<NObsBin; i++) {
	if (xlim[i] > 0.) {nxlim++;}
      }
      if (nxlim != NObsBin ) {
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
	int nvar = sscanf(line,"      xlim[ %*u ] = %le , scale1lo[ %*u ] = %le , scale1hi[ %*u ] = %le , scale2lo[ %*u ] = %le , scale2hi[ %*u ] = %le;",
			  &xlim[i],&scale1lo[i],&scale1hi[i],&scale2lo[i],&scale2hi[i]);
	if (nvar != 5) {
	  cerr << "fastNLO: ERROR! x limits line did not match expected format:" << endl;
	  printf("%s",line);
	  exit(1);
	}
	i++;
      }
      if (i != NObsBin ) {
	cerr << "fastNLO: ERROR! Number of x limits read != NObsBin: i = " << i << ", NObsBin = " << NObsBin << endl;
	exit(1);
      }
    }
  }
  // --- print initialized values 
  cout << endl << "fastNLO: Print initialized x and mu limits:" << endl;
  for (int i=0; i<NObsBin; i++) {
    printf("      xlim[ %u ] = %e , scale1lo[ %u ] = %e , scale1hi[ %u ] = %e , scale2lo[ %u ] = %e , scale2hi[ %u ] = %e ;\n",
	   i,xlim[i],i,scale1lo[i],i,scale1hi[i],i,scale2lo[i],i,scale2hi[i]);
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
