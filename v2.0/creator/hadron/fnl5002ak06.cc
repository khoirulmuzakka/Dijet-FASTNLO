//
// fastNLO v2.1 author code for fnl5002ak06:
//     ATLAS LHC Dijets Scenario, E_cms = 7 TeV
//     for fastjet anti-kT algo with R=0.6 in E-scheme
//
// 
// ===================== fastNLO user: ============================
// To create your own scenario, it is recommended to take 
// this code, make a copy and edit the relevant changes.
// Important:
// Edit only those lines which are labeled as "fastNLO user"
// and refer to the documentation ("fastNLO authorcode in 
// NLOJET++") for a detailed explanation of the parameters 
// and variables.
// Please keep the order of all statements in inittable
// in order to guarantee a working code.
//
// This file contains the following routines:
//   inputfunc		(-> user edits)
//   psinput		(-> user edits)
//   userfunc		(-> user edits)
//   inittable		(-> user edits)
//   GetWarmupValues	(-> user edits)
//   DefineBinning	(-> user edits)
//   initfunc		(don't touch)
//   writetable		(don't touch)
//   end_of_event	(don't touch)
//   phys_output	(don't touch)
//   GetEcms		(don't touch)
//   GetNj		(don't touch)
//   FillEvent		(don't touch)	
//
// Implementing a new scenario requires to edit:
//  - the jet algorithm ("#include" statement and assignment of "jetclus")
//  - number of jets (inputfunc)
//  - center-of-mass energy (psinput)
//  - compute observable, determine bin No. (userfunc)
//  - declare all variables for table, define bin boundaries (inittable, etc...)
//  
// ================================================================
// 
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
   //     Possible choices
   //       - nj = 1U
   //       - nj = 2U
   //       - nj = 3U
   nj = 2U;    

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
   double ptmean_TeV	= ( jet1->perp() + jet2->perp() ) /2./1000.;
   double ptsum_TeV 	= ( jet1->perp() + jet2->perp() ) /2./1000.;
   double M12_TeV	= ((*jet1)+(*jet2)).mag()/1000.;
   double ptmax_TeV	= jet1->perp()/1000.;
   
   // --- fill fastNLO arrays
   //     Scales must be in GeV or dimensionless (at least one in GeV)
   //     Values must be in same dimension as your binning is defined
   FillEvent( M12_TeV, ypsstar , ptmax_TeV*1000. , ypsstar , p , amp );

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
   A1->SetScenName("fnl5002-ATLAS-Dijet-06");

   // --- fastNLO user: up to 20 strings and any number of lines to describe the scenario
   A2->ScDescript.push_back("d2sigma/dM12dy* (pb_TeV)");
   A2->ScDescript.push_back("ATLAS Collaboration");
   A2->ScDescript.push_back("Dijet invariant mass");
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
   B->SetScale1Name( "pt_max" );
   B->SetScale2Name( "y*" );


   // ---- Choose function for ScaleNode distances ---- //
   // --- fastNLO user: possibility to choose function which
   //     is used for binning of scale nodes
   //     possible choices
   //        - LogLog Binning ( H(mu) = log(log(mu/0.25)) )
   //        - Linear         ( H(mu) = mu )
   B->InitLogLogScaleNode( A2 , scale1lo , scale1hi , 1 );	// choose function H for scale 1
   B->InitLinearScaleNode( A2 , scale2lo , scale2hi , 2 );	// choose function H for scale 2


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
   // --- fastNLO user: bin definitions - here in y* and m12
   // 1)
   const int nrapbins = 9;
   // 2)
   double rapbins[nrapbins+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 , 3.5 , 4.0, 4.4 };
   // 3)
   const int nptbins[nrapbins] = {20,20,21,20,19,16,12,9,2};
   
   // prepare 4)
   double ptbins[nrapbins][22] = {
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
   vector<double*> vptbins(nrapbins);
   for (unsigned int i=0;i<vptbins.size();i++) vptbins[i]=ptbins[i]; 
   
   // ---- pass arrays to FnloTable and init bingrid ---- //
   table->GetBlockA2()->InitBinning( nrapbins , rapbins , nptbins , vptbins );
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
      // 2090000000 contributions (!= events) in warm-up run
      xlim [0]  = 1.25e-04 , scale1lo [0] =   30.7000 , scale1hi [0] =  107.3000 , scale2lo [0] =    0.0000 , scale2hi [0] =    0.5000;
      xlim [1]  = 2.48e-04 , scale1lo [1] =   47.8000 , scale1hi [1] =  129.0000 , scale2lo [1] =    0.0000 , scale2hi [1] =    0.5000;
      xlim [2]  = 5.23e-04 , scale1lo [2] =   69.5000 , scale1hi [2] =  148.2000 , scale2lo [2] =    0.0000 , scale2hi [2] =    0.5000;
      xlim [3]  = 9.00e-04 , scale1lo [3] =   91.3000 , scale1hi [3] =  183.5000 , scale2lo [3] =    0.0000 , scale2hi [3] =    0.5000;
      xlim [4]  = 1.38e-03 , scale1lo [4] =  112.8000 , scale1hi [4] =  218.6000 , scale2lo [4] =    0.0000 , scale2hi [4] =    0.5000;
      xlim [5]  = 1.96e-03 , scale1lo [5] =  135.0000 , scale1hi [5] =  260.6000 , scale2lo [5] =    0.0000 , scale2hi [5] =    0.5000;
      xlim [6]  = 2.79e-03 , scale1lo [6] =  160.8000 , scale1hi [6] =  310.5000 , scale2lo [6] =    0.0000 , scale2hi [6] =    0.5000;
      xlim [7]  = 3.95e-03 , scale1lo [7] =  191.1000 , scale1hi [7] =  359.7000 , scale2lo [7] =    0.0000 , scale2hi [7] =    0.5000;
      xlim [8]  = 5.31e-03 , scale1lo [8] =  221.9000 , scale1hi [8] =  416.2000 , scale2lo [8] =    0.0000 , scale2hi [8] =    0.5000;
      xlim [9]  = 7.11e-03 , scale1lo [9] =  256.3000 , scale1hi [9] =  472.7000 , scale2lo [9] =    0.0000 , scale2hi [9] =    0.5000;
      xlim [10] = 9.16e-03 , scale1lo [10] =  291.4000 , scale1hi [10] =  536.2000 , scale2lo [10] =    0.0000 , scale2hi [10] =    0.5000;
      xlim [11] = 1.17e-02 , scale1lo [11] =  330.5000 , scale1hi [11] =  599.5000 , scale2lo [11] =    0.0000 , scale2hi [11] =    0.5000;
      xlim [12] = 1.47e-02 , scale1lo [12] =  370.2000 , scale1hi [12] =  669.9000 , scale2lo [12] =    0.0000 , scale2hi [12] =    0.5000;
      xlim [13] = 1.84e-02 , scale1lo [13] =  412.3000 , scale1hi [13] =  747.9000 , scale2lo [13] =    0.0000 , scale2hi [13] =    0.5000;
      xlim [14] = 2.29e-02 , scale1lo [14] =  460.8000 , scale1hi [14] =  831.9000 , scale2lo [14] =    0.0000 , scale2hi [14] =    0.5000;
      xlim [15] = 2.84e-02 , scale1lo [15] =  513.1000 , scale1hi [15] =  924.5000 , scale2lo [15] =    0.0000 , scale2hi [15] =    0.5000;
      xlim [16] = 3.50e-02 , scale1lo [16] =  569.5000 , scale1hi [16] = 1023.8000 , scale2lo [16] =    0.0000 , scale2hi [16] =    0.5000;
      xlim [17] = 4.29e-02 , scale1lo [17] =  630.3000 , scale1hi [17] = 1130.6000 , scale2lo [17] =    0.0000 , scale2hi [17] =    0.5000;
      xlim [18] = 5.22e-02 , scale1lo [18] =  696.0000 , scale1hi [18] = 1368.5000 , scale2lo [18] =    0.0000 , scale2hi [18] =    0.5000;
      xlim [19] = 7.68e-02 , scale1lo [19] =  843.7000 , scale1hi [19] = 1960.7000 , scale2lo [19] =    0.0000 , scale2hi [19] =    0.5000;
      xlim [20] = 3.24e-04 , scale1lo [20] =   35.2000 , scale1hi [20] =  117.9000 , scale2lo [20] =    0.5000 , scale2hi [20] =    1.0000;
      xlim [21] = 5.23e-04 , scale1lo [21] =   51.2000 , scale1hi [21] =  133.4000 , scale2lo [21] =    0.5000 , scale2hi [21] =    1.0000;
      xlim [22] = 9.01e-04 , scale1lo [22] =   67.0000 , scale1hi [22] =  161.6000 , scale2lo [22] =    0.5000 , scale2hi [22] =    1.0000;
      xlim [23] = 1.38e-03 , scale1lo [23] =   82.8000 , scale1hi [23] =  193.4000 , scale2lo [23] =    0.5000 , scale2hi [23] =    1.0000;
      xlim [24] = 1.96e-03 , scale1lo [24] =   98.9000 , scale1hi [24] =  230.1000 , scale2lo [24] =    0.5000 , scale2hi [24] =    1.0000;
      xlim [25] = 2.79e-03 , scale1lo [25] =  118.0000 , scale1hi [25] =  274.7000 , scale2lo [25] =    0.5000 , scale2hi [25] =    1.0000;
      xlim [26] = 3.95e-03 , scale1lo [26] =  139.9000 , scale1hi [26] =  317.5000 , scale2lo [26] =    0.5000 , scale2hi [26] =    1.0000;
      xlim [27] = 5.32e-03 , scale1lo [27] =  162.4000 , scale1hi [27] =  369.0000 , scale2lo [27] =    0.5000 , scale2hi [27] =    1.0000;
      xlim [28] = 7.11e-03 , scale1lo [28] =  187.4000 , scale1hi [28] =  417.2000 , scale2lo [28] =    0.5000 , scale2hi [28] =    1.0000;
      xlim [29] = 9.17e-03 , scale1lo [29] =  213.6000 , scale1hi [29] =  473.4000 , scale2lo [29] =    0.5000 , scale2hi [29] =    1.0000;
      xlim [30] = 1.17e-02 , scale1lo [30] =  242.0000 , scale1hi [30] =  530.0000 , scale2lo [30] =    0.5000 , scale2hi [30] =    1.0000;
      xlim [31] = 1.47e-02 , scale1lo [31] =  269.9000 , scale1hi [31] =  592.4000 , scale2lo [31] =    0.5000 , scale2hi [31] =    1.0000;
      xlim [32] = 1.84e-02 , scale1lo [32] =  301.3000 , scale1hi [32] =  659.6000 , scale2lo [32] =    0.5000 , scale2hi [32] =    1.0000;
      xlim [33] = 2.29e-02 , scale1lo [33] =  337.5000 , scale1hi [33] =  736.5000 , scale2lo [33] =    0.5000 , scale2hi [33] =    1.0000;
      xlim [34] = 2.84e-02 , scale1lo [34] =  375.7000 , scale1hi [34] =  816.2000 , scale2lo [34] =    0.5000 , scale2hi [34] =    1.0000;
      xlim [35] = 3.50e-02 , scale1lo [35] =  417.0000 , scale1hi [35] =  907.0000 , scale2lo [35] =    0.5000 , scale2hi [35] =    1.0000;
      xlim [36] = 4.29e-02 , scale1lo [36] =  460.6000 , scale1hi [36] =  995.2000 , scale2lo [36] =    0.5000 , scale2hi [36] =    1.0000;
      xlim [37] = 5.22e-02 , scale1lo [37] =  508.1000 , scale1hi [37] = 1095.8000 , scale2lo [37] =    0.5000 , scale2hi [37] =    1.0000;
      xlim [38] = 6.32e-02 , scale1lo [38] =  560.9000 , scale1hi [38] = 1325.7000 , scale2lo [38] =    0.5000 , scale2hi [38] =    1.0000;
      xlim [39] = 9.17e-02 , scale1lo [39] =  676.4000 , scale1hi [39] = 1436.2000 , scale2lo [39] =    0.5000 , scale2hi [39] =    1.0000;
      xlim [40] = 7.44e-04 , scale1lo [40] =   33.5000 , scale1hi [40] =  112.5000 , scale2lo [40] =    1.0000 , scale2hi [40] =    1.5000;
      xlim [41] = 1.02e-03 , scale1lo [41] =   44.2000 , scale1hi [41] =  124.7000 , scale2lo [41] =    1.0000 , scale2hi [41] =    1.5000;
      xlim [42] = 1.38e-03 , scale1lo [42] =   54.5000 , scale1hi [42] =  140.8000 , scale2lo [42] =    1.0000 , scale2hi [42] =    1.5000;
      xlim [43] = 1.96e-03 , scale1lo [43] =   65.2000 , scale1hi [43] =  168.5000 , scale2lo [43] =    1.0000 , scale2hi [43] =    1.5000;
      xlim [44] = 2.80e-03 , scale1lo [44] =   77.8000 , scale1hi [44] =  199.5000 , scale2lo [44] =    1.0000 , scale2hi [44] =    1.5000;
      xlim [45] = 3.95e-03 , scale1lo [45] =   92.5000 , scale1hi [45] =  232.3000 , scale2lo [45] =    1.0000 , scale2hi [45] =    1.5000;
      xlim [46] = 5.31e-03 , scale1lo [46] =  106.8000 , scale1hi [46] =  269.1000 , scale2lo [46] =    1.0000 , scale2hi [46] =    1.5000;
      xlim [47] = 7.12e-03 , scale1lo [47] =  123.9000 , scale1hi [47] =  305.6000 , scale2lo [47] =    1.0000 , scale2hi [47] =    1.5000;
      xlim [48] = 9.18e-03 , scale1lo [48] =  140.1000 , scale1hi [48] =  345.5000 , scale2lo [48] =    1.0000 , scale2hi [48] =    1.5000;
      xlim [49] = 1.17e-02 , scale1lo [49] =  158.9000 , scale1hi [49] =  384.5000 , scale2lo [49] =    1.0000 , scale2hi [49] =    1.5000;
      xlim [50] = 1.47e-02 , scale1lo [50] =  177.5000 , scale1hi [50] =  433.1000 , scale2lo [50] =    1.0000 , scale2hi [50] =    1.5000;
      xlim [51] = 1.84e-02 , scale1lo [51] =  198.5000 , scale1hi [51] =  482.7000 , scale2lo [51] =    1.0000 , scale2hi [51] =    1.5000;
      xlim [52] = 2.29e-02 , scale1lo [52] =  221.1000 , scale1hi [52] =  535.9000 , scale2lo [52] =    1.0000 , scale2hi [52] =    1.5000;
      xlim [53] = 2.84e-02 , scale1lo [53] =  247.2000 , scale1hi [53] =  594.0000 , scale2lo [53] =    1.0000 , scale2hi [53] =    1.5000;
      xlim [54] = 3.50e-02 , scale1lo [54] =  275.1000 , scale1hi [54] =  659.7000 , scale2lo [54] =    1.0000 , scale2hi [54] =    1.5000;
      xlim [55] = 4.29e-02 , scale1lo [55] =  302.9000 , scale1hi [55] =  726.2000 , scale2lo [55] =    1.0000 , scale2hi [55] =    1.5000;
      xlim [56] = 5.22e-02 , scale1lo [56] =  334.7000 , scale1hi [56] =  799.5000 , scale2lo [56] =    1.0000 , scale2hi [56] =    1.5000;
      xlim [57] = 6.32e-02 , scale1lo [57] =  368.3000 , scale1hi [57] =  882.5000 , scale2lo [57] =    1.0000 , scale2hi [57] =    1.5000;
      xlim [58] = 7.68e-02 , scale1lo [58] =  406.0000 , scale1hi [58] =  968.4000 , scale2lo [58] =    1.0000 , scale2hi [58] =    1.5000;
      xlim [59] = 9.17e-02 , scale1lo [59] =  445.4000 , scale1hi [59] = 1159.2000 , scale2lo [59] =    1.0000 , scale2hi [59] =    1.5000;
      xlim [60] = 1.32e-01 , scale1lo [60] =  534.6000 , scale1hi [60] = 1636.0000 , scale2lo [60] =    1.0000 , scale2hi [60] =    1.5000;
      xlim [61] = 1.92e-03 , scale1lo [61] =   34.2000 , scale1hi [61] =  109.8000 , scale2lo [61] =    1.5000 , scale2hi [61] =    2.0000;
      xlim [62] = 2.47e-03 , scale1lo [62] =   40.7000 , scale1hi [62] =  117.5000 , scale2lo [62] =    1.5000 , scale2hi [62] =    2.0000;
      xlim [63] = 2.93e-03 , scale1lo [63] =   48.4000 , scale1hi [63] =  133.9000 , scale2lo [63] =    1.5000 , scale2hi [63] =    2.0000;
      xlim [64] = 3.98e-03 , scale1lo [64] =   57.8000 , scale1hi [64] =  151.9000 , scale2lo [64] =    1.5000 , scale2hi [64] =    2.0000;
      xlim [65] = 5.34e-03 , scale1lo [65] =   67.1000 , scale1hi [65] =  174.4000 , scale2lo [65] =    1.5000 , scale2hi [65] =    2.0000;
      xlim [66] = 7.11e-03 , scale1lo [66] =   77.1000 , scale1hi [66] =  200.5000 , scale2lo [66] =    1.5000 , scale2hi [66] =    2.0000;
      xlim [67] = 9.18e-03 , scale1lo [67] =   87.8000 , scale1hi [67] =  225.5000 , scale2lo [67] =    1.5000 , scale2hi [67] =    2.0000;
      xlim [68] = 1.17e-02 , scale1lo [68] =   99.9000 , scale1hi [68] =  252.8000 , scale2lo [68] =    1.5000 , scale2hi [68] =    2.0000;
      xlim [69] = 1.47e-02 , scale1lo [69] =  112.0000 , scale1hi [69] =  281.9000 , scale2lo [69] =    1.5000 , scale2hi [69] =    2.0000;
      xlim [70] = 1.84e-02 , scale1lo [70] =  124.6000 , scale1hi [70] =  313.3000 , scale2lo [70] =    1.5000 , scale2hi [70] =    2.0000;
      xlim [71] = 2.29e-02 , scale1lo [71] =  139.9000 , scale1hi [71] =  350.3000 , scale2lo [71] =    1.5000 , scale2hi [71] =    2.0000;
      xlim [72] = 2.84e-02 , scale1lo [72] =  154.4000 , scale1hi [72] =  390.2000 , scale2lo [72] =    1.5000 , scale2hi [72] =    2.0000;
      xlim [73] = 3.50e-02 , scale1lo [73] =  172.3000 , scale1hi [73] =  431.2000 , scale2lo [73] =    1.5000 , scale2hi [73] =    2.0000;
      xlim [74] = 4.29e-02 , scale1lo [74] =  190.0000 , scale1hi [74] =  475.6000 , scale2lo [74] =    1.5000 , scale2hi [74] =    2.0000;
      xlim [75] = 5.22e-02 , scale1lo [75] =  209.2000 , scale1hi [75] =  526.5000 , scale2lo [75] =    1.5000 , scale2hi [75] =    2.0000;
      xlim [76] = 6.32e-02 , scale1lo [76] =  230.8000 , scale1hi [76] =  575.1000 , scale2lo [76] =    1.5000 , scale2hi [76] =    2.0000;
      xlim [77] = 7.68e-02 , scale1lo [77] =  253.5000 , scale1hi [77] =  627.2000 , scale2lo [77] =    1.5000 , scale2hi [77] =    2.0000;
      xlim [78] = 9.18e-02 , scale1lo [78] =  278.7000 , scale1hi [78] =  691.4000 , scale2lo [78] =    1.5000 , scale2hi [78] =    2.0000;
      xlim [79] = 1.10e-01 , scale1lo [79] =  306.2000 , scale1hi [79] =  824.6000 , scale2lo [79] =    1.5000 , scale2hi [79] =    2.0000;
      xlim [80] = 1.57e-01 , scale1lo [80] =  366.0000 , scale1hi [80] = 1161.7000 , scale2lo [80] =    1.5000 , scale2hi [80] =    2.0000;
      xlim [81] = 4.54e-03 , scale1lo [81] =   30.2000 , scale1hi [81] =   99.0000 , scale2lo [81] =    2.0000 , scale2hi [81] =    2.5000;
      xlim [82] = 5.09e-03 , scale1lo [82] =   35.2000 , scale1hi [82] =  111.5000 , scale2lo [82] =    2.0000 , scale2hi [82] =    2.5000;
      xlim [83] = 6.66e-03 , scale1lo [83] =   41.2000 , scale1hi [83] =  116.5000 , scale2lo [83] =    2.0000 , scale2hi [83] =    2.5000;
      xlim [84] = 7.68e-03 , scale1lo [84] =   48.0000 , scale1hi [84] =  125.0000 , scale2lo [84] =    2.0000 , scale2hi [84] =    2.5000;
      xlim [85] = 9.18e-03 , scale1lo [85] =   54.3000 , scale1hi [85] =  141.2000 , scale2lo [85] =    2.0000 , scale2hi [85] =    2.5000;
      xlim [86] = 1.17e-02 , scale1lo [86] =   61.6000 , scale1hi [86] =  156.4000 , scale2lo [86] =    2.0000 , scale2hi [86] =    2.5000;
      xlim [87] = 1.47e-02 , scale1lo [87] =   68.9000 , scale1hi [87] =  176.0000 , scale2lo [87] =    2.0000 , scale2hi [87] =    2.5000;
      xlim [88] = 1.84e-02 , scale1lo [88] =   77.1000 , scale1hi [88] =  198.2000 , scale2lo [88] =    2.0000 , scale2hi [88] =    2.5000;
      xlim [89] = 2.29e-02 , scale1lo [89] =   85.8000 , scale1hi [89] =  218.4000 , scale2lo [89] =    2.0000 , scale2hi [89] =    2.5000;
      xlim [90] = 2.84e-02 , scale1lo [90] =   95.1000 , scale1hi [90] =  240.2000 , scale2lo [90] =    2.0000 , scale2hi [90] =    2.5000;
      xlim [91] = 3.50e-02 , scale1lo [91] =  105.7000 , scale1hi [91] =  268.9000 , scale2lo [91] =    2.0000 , scale2hi [91] =    2.5000;
      xlim [92] = 4.29e-02 , scale1lo [92] =  116.9000 , scale1hi [92] =  296.9000 , scale2lo [92] =    2.0000 , scale2hi [92] =    2.5000;
      xlim [93] = 5.22e-02 , scale1lo [93] =  129.5000 , scale1hi [93] =  325.2000 , scale2lo [93] =    2.0000 , scale2hi [93] =    2.5000;
      xlim [94] = 6.33e-02 , scale1lo [94] =  142.2000 , scale1hi [94] =  361.1000 , scale2lo [94] =    2.0000 , scale2hi [94] =    2.5000;
      xlim [95] = 7.68e-02 , scale1lo [95] =  156.1000 , scale1hi [95] =  394.5000 , scale2lo [95] =    2.0000 , scale2hi [95] =    2.5000;
      xlim [96] = 9.17e-02 , scale1lo [96] =  171.4000 , scale1hi [96] =  428.3000 , scale2lo [96] =    2.0000 , scale2hi [96] =    2.5000;
      xlim [97] = 1.10e-01 , scale1lo [97] =  188.7000 , scale1hi [97] =  473.7000 , scale2lo [97] =    2.0000 , scale2hi [97] =    2.5000;
      xlim [98] = 1.32e-01 , scale1lo [98] =  205.4000 , scale1hi [98] =  560.6000 , scale2lo [98] =    2.0000 , scale2hi [98] =    2.5000;
      xlim [99] = 1.88e-01 , scale1lo [99] =  247.9000 , scale1hi [99] =  787.4000 , scale2lo [99] =    2.0000 , scale2hi [99] =    2.5000;
      xlim [100] = 1.16e-02 , scale1lo [100] =   33.2000 , scale1hi [100] =   96.8000 , scale2lo [100] =    2.5000 , scale2hi [100] =    3.0000;
      xlim [101] = 1.41e-02 , scale1lo [101] =   37.7000 , scale1hi [101] =  103.1000 , scale2lo [101] =    2.5000 , scale2hi [101] =    3.0000;
      xlim [102] = 1.65e-02 , scale1lo [102] =   41.9000 , scale1hi [102] =  109.9000 , scale2lo [102] =    2.5000 , scale2hi [102] =    3.0000;
      xlim [103] = 2.01e-02 , scale1lo [103] =   47.1000 , scale1hi [103] =  120.0000 , scale2lo [103] =    2.5000 , scale2hi [103] =    3.0000;
      xlim [104] = 2.29e-02 , scale1lo [104] =   52.6000 , scale1hi [104] =  135.5000 , scale2lo [104] =    2.5000 , scale2hi [104] =    3.0000;
      xlim [105] = 2.84e-02 , scale1lo [105] =   58.6000 , scale1hi [105] =  148.7000 , scale2lo [105] =    2.5000 , scale2hi [105] =    3.0000;
      xlim [106] = 3.50e-02 , scale1lo [106] =   65.0000 , scale1hi [106] =  161.9000 , scale2lo [106] =    2.5000 , scale2hi [106] =    3.0000;
      xlim [107] = 4.29e-02 , scale1lo [107] =   72.0000 , scale1hi [107] =  181.7000 , scale2lo [107] =    2.5000 , scale2hi [107] =    3.0000;
      xlim [108] = 5.22e-02 , scale1lo [108] =   79.5000 , scale1hi [108] =  197.2000 , scale2lo [108] =    2.5000 , scale2hi [108] =    3.0000;
      xlim [109] = 6.32e-02 , scale1lo [109] =   86.9000 , scale1hi [109] =  219.8000 , scale2lo [109] =    2.5000 , scale2hi [109] =    3.0000;
      xlim [110] = 7.68e-02 , scale1lo [110] =   95.3000 , scale1hi [110] =  239.8000 , scale2lo [110] =    2.5000 , scale2hi [110] =    3.0000;
      xlim [111] = 9.20e-02 , scale1lo [111] =  105.0000 , scale1hi [111] =  261.4000 , scale2lo [111] =    2.5000 , scale2hi [111] =    3.0000;
      xlim [112] = 1.10e-01 , scale1lo [112] =  114.3000 , scale1hi [112] =  285.4000 , scale2lo [112] =    2.5000 , scale2hi [112] =    3.0000;
      xlim [113] = 1.32e-01 , scale1lo [113] =  125.8000 , scale1hi [113] =  312.8000 , scale2lo [113] =    2.5000 , scale2hi [113] =    3.0000;
      xlim [114] = 1.57e-01 , scale1lo [114] =  137.7000 , scale1hi [114] =  373.9000 , scale2lo [114] =    2.5000 , scale2hi [114] =    3.0000;
      xlim [115] = 2.23e-01 , scale1lo [115] =  163.9000 , scale1hi [115] =  528.1000 , scale2lo [115] =    2.5000 , scale2hi [115] =    3.0000;
      xlim [116] = 3.64e-02 , scale1lo [116] =   35.7000 , scale1hi [116] =  101.5000 , scale2lo [116] =    3.0000 , scale2hi [116] =    3.5000;
      xlim [117] = 4.12e-02 , scale1lo [117] =   39.5000 , scale1hi [117] =  105.6000 , scale2lo [117] =    3.0000 , scale2hi [117] =    3.5000;
      xlim [118] = 4.77e-02 , scale1lo [118] =   43.8000 , scale1hi [118] =  108.5000 , scale2lo [118] =    3.0000 , scale2hi [118] =    3.5000;
      xlim [119] = 5.62e-02 , scale1lo [119] =   48.3000 , scale1hi [119] =  118.5000 , scale2lo [119] =    3.0000 , scale2hi [119] =    3.5000;
      xlim [120] = 6.37e-02 , scale1lo [120] =   53.1000 , scale1hi [120] =  132.2000 , scale2lo [120] =    3.0000 , scale2hi [120] =    3.5000;
      xlim [121] = 7.71e-02 , scale1lo [121] =   58.5000 , scale1hi [121] =  145.8000 , scale2lo [121] =    3.0000 , scale2hi [121] =    3.5000;
      xlim [122] = 9.19e-02 , scale1lo [122] =   63.8000 , scale1hi [122] =  158.6000 , scale2lo [122] =    3.0000 , scale2hi [122] =    3.5000;
      xlim [123] = 1.10e-01 , scale1lo [123] =   70.4000 , scale1hi [123] =  174.0000 , scale2lo [123] =    3.0000 , scale2hi [123] =    3.5000;
      xlim [124] = 1.32e-01 , scale1lo [124] =   77.0000 , scale1hi [124] =  186.5000 , scale2lo [124] =    3.0000 , scale2hi [124] =    3.5000;
      xlim [125] = 1.57e-01 , scale1lo [125] =   83.1000 , scale1hi [125] =  206.1000 , scale2lo [125] =    3.0000 , scale2hi [125] =    3.5000;
      xlim [126] = 1.89e-01 , scale1lo [126] =   91.8000 , scale1hi [126] =  242.0000 , scale2lo [126] =    3.0000 , scale2hi [126] =    3.5000;
      xlim [127] = 2.66e-01 , scale1lo [127] =  109.0000 , scale1hi [127] =  344.1000 , scale2lo [127] =    3.0000 , scale2hi [127] =    3.5000;
      xlim [128] = 8.80e-02 , scale1lo [128] =   32.3000 , scale1hi [128] =   86.5000 , scale2lo [128] =    3.5000 , scale2hi [128] =    4.0000;
      xlim [129] = 9.84e-02 , scale1lo [129] =   35.6000 , scale1hi [129] =   88.3000 , scale2lo [129] =    3.5000 , scale2hi [129] =    4.0000;
      xlim [130] = 9.96e-02 , scale1lo [130] =   38.9000 , scale1hi [130] =   95.2000 , scale2lo [130] =    3.5000 , scale2hi [130] =    4.0000;
      xlim [131] = 1.16e-01 , scale1lo [131] =   42.7000 , scale1hi [131] =  103.8000 , scale2lo [131] =    3.5000 , scale2hi [131] =    4.0000;
      xlim [132] = 1.44e-01 , scale1lo [132] =   46.7000 , scale1hi [132] =  114.0000 , scale2lo [132] =    3.5000 , scale2hi [132] =    4.0000;
      xlim [133] = 1.61e-01 , scale1lo [133] =   51.0000 , scale1hi [133] =  122.3000 , scale2lo [133] =    3.5000 , scale2hi [133] =    4.0000;
      xlim [134] = 1.7e-01 , scale1lo [134] =   55.7000 , scale1hi [134] =  136.6000 , scale2lo [134] =    3.5000 , scale2hi [134] =    4.0000;
      xlim [135] = 2.0e-01 , scale1lo [135] =   60.6000 , scale1hi [135] =  160.0000 , scale2lo [135] =    3.5000 , scale2hi [135] =    4.0000;
      xlim [136] = 3.0e-01 , scale1lo [136] =   71.7000 , scale1hi [136] =  219.4000 , scale2lo [136] =    3.5000 , scale2hi [136] =    4.0000;
      xlim [137] = 1.7e-01 , scale1lo [137] =   31.7000 , scale1hi [137] =   85.6000 , scale2lo [137] =    4.0000 , scale2hi [137] =    4.4000;
      xlim [138] = 2.1e-01 , scale1lo [138] =   37.8000 , scale1hi [138] =   99.0000 , scale2lo [138] =    4.0000 , scale2hi [138] =    4.4000;
      // --------- fastNLO: Warm-Up run results (end)
   }
   else {
      printf("fastNLO: This is a warm-up run!\n");
      // --- fastNLO user: You can set the number of contributions
      //     after which the WarmUp values are printed
      B->SetWarmUpPrint(50000000);		// default 10000000

      // safe initialziations
      for(int i=0;i<NObsBin;i++){ 
	 xlim[i] = 1.1e-07;
	 scale1lo[i] = 1.0, scale1hi[i]=9.9e10;
	 scale2lo[i] = 1.0, scale2hi[i]=9.9e10;
      }
   }
}



// --- fastNLO user: nothing further todo


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



void UserHHC::initfunc(unsigned int)
{
  // --- Initialize event counters 
  nevents = 0;
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
  // --- store table and show calculation time
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

  if(strstr(file,"-born-")!=NULL){
    nlo = false;
    printf("fastNLO: This is a LO run!\n");
  }else{
    if(strstr(file,"-nlo-")!=NULL){
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

