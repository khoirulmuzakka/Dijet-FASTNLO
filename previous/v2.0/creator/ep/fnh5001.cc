//
// fastNLO v2.1 author code for fnh5001:
//     H1, High Q2, HERA-II, Inclusive jets Scenario, kt
//
// 
// ===================== fastNLO user: ============================
// To create your own DIS scenario, it is recommended to take 
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
#include <bits/dis-phasespace.h>
#include <bits/dis-process.h>
#include <bits/dis-jetfunc.h>

//----- used namespaces -----
using namespace nlo;
using namespace std;

// ---- fastNLO ----
#include <fastnlotk/fastNLOCreate.h>

//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
void psinput(phasespace_dis *, double&, double&, double&, double&, double&, double&, double&, double&);
user_base_dis * userfunc();

//----- array of the symbols symbols -----
extern "C"{
  
  struct { 
	const char *name;
	void *address;
  } user_defined_functions[] = 
  {
	//   process index: hhc for e + p --> jets (DIS)
	{"procindex", (void *) "dis"},
  
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
#include <kT_clus.h>
#include "pdf-cteq6.h"
#include "pdf-dis-dummy.h"
#include "fnloTable.h"
#include "fnloBlockBNlojet.h"
//#include "fj-ak.h"   // fastNLO user: .h file for jet algorithm


// ---- fastNLO ----
#include "fastNLOInterfaceNLOJETDIS.cc"


class UserDIS : public basic_user_set<user1d_dis, user1h_dis> 
{
public:
  //   init and user function
   void initfunc(unsigned int);
   void userfunc(const event_dis&, const amplitude_dis&);

   virtual void end_of_event();  
   virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);
  
private:
   //   pdf
   pdf_cteq6dis pdf;
   pdf_dis_dummy dummypdf;

   // algorithms
   kT_clus_long jetclus;
   //fj_ak fjclus;   // fastNLO user: define jet algorithm (consistent with .h file above)
   
   // the jet structore in breit & lab. frame
   bounded_vector<lorentzvector<double> > pjb, pjl; 
   bounded_vector<unsigned int> jet;
   
   //  event in breit frame
   event_dis pbreit;
  
   // boost the jet momenta back to the laboratory frame
   void boost_back_to_lab(const event_dis&);
   
   //fastNLO starts here   
   fnloTable *table;
   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table

   string tablefilename; // The table file to write to
   time_t start_time;
   double* xlim;
   double* scale1lo;
   double* scale1hi;
   double* scale2lo;
   double* scale2hi;

   fastNLOCreate *ftable;
   void InitFastNLO(const std::basic_string<char>& fname);
   
   bool nlo;
   void inittable();
   void writetable();
   unsigned int GetNj();
   double GetEcms();
   void DefineBinning();
   void GetWarmupValues( fnloBlockBNlojet* B );
   void FillEvent( double val1 , double val2 , double mu1, double mu2 , const event_dis& p, const nlo::amplitude_dis& amp , double alem);
   
};
  

user_base_dis * userfunc() {
  return new UserDIS;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
  //  number of jets
   nj = 2U;
  
  //  number of the up and down type flavours
  nu = 2U;
  nd = 3U;
} 

void psinput(phasespace_dis *ps, double& el, double& eh, double& q2min, double& q2max, 
			 double& xmin, double& xmax, double& ymin, double& ymax)
{
  // --- fastNLO user:
  //	 Define your phase space here. Phase space cuts
  //     can still be performed in UserDIS::userfunc.
   
  // energy of the incoming lepton and 
  // hadron in the laboratory frame
  el = 27.6;    // GeV
  eh = 920.0;     // GeV
  
  //  Q^2 cuts
  q2min = 150.0;    // GeV^2
  q2max = 15000.0;  // GeV^2 

  //   xB cuts
  xmin = 0.0;
  xmax = 1.0;
  
  //   y cuts
  ymin = 0.2;
  ymax = 0.7;

  //   You can use your own phase generator. 
  //   Here we use the default.
  ps = 0;
} 


extern"C" double xalfaem_(double *);

double xalpha_em(double mq2) {
   return xalfaem_(&mq2);
}



void UserDIS::userfunc(const event_dis& p, const amplitude_dis& amp)
{

   // --- fastNLO user:
   //     Here is your playground where you compute your observable 
   //     and define the two possible scale variables.
   //     In DIS, the first scale variable must be 'Q'.
   //     All observables must be in same dimension as your bingrid is defined.


   // ---- run the jet finder in bootsed frame ---
   //----- copy the momenta and boost to the breit frame ----
   pbreit = p; lab_to_breit(pbreit);


//    //---- FastJET ---- //
//    const double jetsize = 1.0;
//    bounded_vector<lorentzvector<double> > pjfj = fjclus(pbreit,jetsize);
//    unsigned int njfj = pjfj.upper();



   //----- do the cluster analysis in the breit frame-----
   jetclus.set_up(1,false,2);
   jetclus(pbreit); jetclus.incl(pjb, jet);
   unsigned int nj = pjb.upper();

   //----- jet structure in laboratory frame -----
   pjl = pjb; boost_back_to_lab(p);



   //----- Momentum transfer -----
   double Q2 = -((p[-1] - p[-2]).mag2());
   double Q = sqrt(Q2);


   // ---- el-mag coupling --- 
   double alem = xalpha_em(Q2);


   //----- cuts -----
   const double pTmin = 7.0, pTmax = 50.0;
   const double etamin = -1.0, etamax = 2.5; 

   
   // ---- apply jet cuts and fill jets array----
   unsigned int IndexArray[nj+1];
   double PtArray[nj+1];
   unsigned int jetnum = 0;
   for(unsigned int i = 1; i <= nj; i++) {
      double pt = pjb[i].perp(); 
      double etalab = pjl[i].prapidity();
      if(pt > pTmin && pt < pTmax && etamin < etalab && etalab < etamax) {
	 jetnum++;
 	 PtArray[jetnum] = pt;
	 IndexArray[jetnum] = i;
      }
   }



//    // ---- look for pt max --- //
//    double pt_max = 0.;
//    for ( unsigned int ij = 0 ; ij<jetnum ; ij++ ){
//       if ( PtArray[ij+1] > pt_max ) pt_max = PtArray[ij+1];
//    }   


   /*
   // -------------------------------------------------------------------------------------- //
   // ---- different scale choice ---- //
   // use pt for 1-jet events, and <pT> for 2+-jet events
   // -------------------------------------------------------------------------------------- //

   // ---- dijet cuts -----
   const double pT2min = 5.0, pT2max = 50.0;
   //const double etamin = -1.0, etamax = 2.5; 
   const double M12Min = 16;


   // ---- apply jet cuts and fill jets array----
   double PtArray2[nj+1];
   unsigned int dijetnum = 0;
   unsigned int IndexArray2[nj+1];
   for(unsigned int i = 1; i <= nj; i++) {
      double pt = pjb[i].perp(); 
      double etalab = pjl[i].prapidity();
      if(pt > pT2min && pt < pT2max && etamin < etalab && etalab < etamax) {
	 dijetnum++;
 	 PtArray2[dijetnum] = pt;
	 IndexArray2[dijetnum] = i;
      }
   }

  
   // ---- sort pt-array ---- //
   double t = 0;
   unsigned int y = 0;
   int fred = 0;
   for(unsigned int i = 1; i <= dijetnum; i++) {
      for (y=1; y < (dijetnum+1-i); y++){
	 if (PtArray2[y] < PtArray2[y+1]){        
	    t=PtArray2[y];
	    fred = IndexArray2[y];
	    PtArray2[y]=PtArray2[y+1];
	    IndexArray2[y]=IndexArray2[y+1];
	    PtArray2[y+1]=t;
	    IndexArray2[y+1]=fred;
	 }          
      }
   }


   double pt_mean = 0;
   double M12	= 0;
   if ( dijetnum >= 2 ) {
      // ---- observables ---- //
      pt_mean = ( PtArray2[1] + PtArray2[2] ) / 2;
      lorentzvector<double> psum = pjb[IndexArray2[1]]+ pjb[IndexArray2[2]];
      M12 = psum.mag();
   }
   
   const bool IsDijet = dijetnum >=2 && M12 > M12Min;

   // new scale:
   double muPtLike = IsDijet ?
      pt_mean :
      PtArray[1];
   */

   // -------------------------------------------------------------------------------------- //

   //if ( jetnum < 3 ) return;

   // --- fastNLO v2.2 ----  nlojet-event ----
   vector<fnloEvent> contribs = UsefulNlojetTools::GetFlexibleScaleNlojetContribDIS(p,amp,dummypdf,alem);

   // ---- Inclusive jets: loop over all jets ----
   for ( unsigned int ij = 0 ; ij<jetnum ; ij++ ){
   
      // ---- define observables ----
      double val1 = PtArray[ij+1];
      double val2 = Q2;
      double mu1 = Q;
      double mu2 = PtArray[ij+1];
      //* double mu2 = muPtLike;
      //double mu2 = pt_max;  // pt_max
	   
      // --- fill fastNLO arrays
      //     Scales must be in GeV or dimensionless (at least one in GeV)
      //     Values must be in same dimension as your binning is defined
      FillEvent( val1 , val2 , mu1 , mu2 , p , amp , alem );


      { // fastNLO v2.2
	 // scenario specific quantites
	 fnloScenario scen;
	 scen.SetObservableDimI( Q2 , 0 );
	 scen.SetObservableDimI( PtArray[ij+1] , 1 );
	 scen.SetObsScale1( Q );
	 scen.SetObsScale2( PtArray[ij+1] );
	 
	 ftable->FillAllSubprocesses(contribs,scen); 
      }

   }   
}


void UserDIS::FillEvent( double val1 , double val2 , double mu1, double mu2 , const event_dis& p, const nlo::amplitude_dis& amp , double alem ){
   // ---- Fill FastNLO Array ---- //
   // --- fastNLO user: usually nothing to do

   // ---- get x-value ---- //
   double x = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);
   double Q2 = -((p[-1] - p[-2]).mag2());
   double Q = sqrt(Q2);
   if ( fabs(mu1-Q)>1.e-6 ) {
      printf("UserDIS::FillEvent. Sorry, but scale1 must be 'Q'\n");
      printf("Q = %9.4f, Q2 = %9.4f, mu1 = %9.4f\n",Q,Q2,mu1);
      exit(1);
   }
   int obsbin = table->GetBlockA2()->GetBinNumber( val1 , val2  );

   //---------- fill fastNLO arrays
   if ( obsbin >= 0 ) {
      // --- prefactor: divide by binwidth - and multiply with alpha_elm^2
      double prefactor = 1./table->GetBlockA2()->BinSize[obsbin] * alem*alem; 
      for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
	 ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventDISMuVar(obsbin,x,mu1,mu2,amp,dummypdf,pdf,prefactor);
      }
   } // - end: fill fastNLO array
}


// ---- some function that can be used for reference calculations
double Fct_x(double x, double y){ return x; }
double Fct_y(double x, double y){ return y; }
double Fct_xyov2(double x, double y){ return (x+y)/2.; }
double Fct_x2y2ov2(double x, double y){ return sqrt((x*x+y*y)/2.); }
double Fct_x_exp03y(double x, double y){ return x*exp(0.3*y); }
// --- fastNLO user; Feel free to define any other function here


void UserDIS::inittable(){
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
   A2->Ecms	    = GetEcms();
   A2->ILOord	= GetNj()-1;


   // ---- set scenario name and description ---- //
   // --- fastNLO user: set scenario name (no whitespaces!)
   table->GetBlockA1()->SetScenName("fnh5001");

   // --- fastNLO user: up to 20 strings and any number of lines to describe the scenario
   A2->ScDescript.push_back("Inclusive jets - d2sigma/dpTdQ2 [pb]");
   A2->ScDescript.push_back("H1 Collaboration");
   A2->ScDescript.push_back("Inclusive jet cross sections in DIS");
   A2->ScDescript.push_back("kT-algorithm R=1, High Q2, Hera-II");
   A2->ScDescript.push_back("H1prelim-11-032, Fig. 7");

   
   // --- fastNLO user: Give information about your measurement
   A2->SetIpublunits( 12 );			// --- fastNLO user: set cross section units (negative power of ten), e.g. 'pb' -> 12.
   A2->SetNumDiffBin( 2 );			// --- fastNLO user: No of dimensions in which observable is binned
   bool IsDiffBin = false;			// --- fastNLO user: Are publication units divided by this variable?
   A2->SetDimLabel("p_T", 1 , IsDiffBin );// --- fastNLO user: label of 1st dimension
   A2->SetDimLabel("Q^2",  2 , IsDiffBin );// --- fastNLO user: label of 2nd dimension

   // ---- Define your bingrid in method DefineBinning() ---- //
   // --- fastNLO user: Modifiy function DefineBinning() below according to your bin grid.
   DefineBinning();



   // ---- initalize BlockB ---- //
   // --- fastNLO user: nothing to do.
   // ---- initalize all constants to make fastNLO behave like for DIS scenario ---- //
   B->InitDISConstants(A2,nlo);



   // ---- initialize variables for WarmUp run ---- //
   // --- fastNLO user: Start "Warm-Up" or "Production" run.
   //     See documentation or GetWarmupValues() for more details.
   //     choices for B->IWarmUp
   //	    -  B->SetDoWarmUp(true)   ->  Do the Warm-Up run
   //	    -  B->SetDoWarmUp(false)  ->  Do a production run
   B->SetDoWarmUp(false);


   // ---- get warm up values or init arrays reasonably ---- //
   // --- fastNLO user: See documentation in GetWarmUpValues for details.
   //     Do not modify this call.
   GetWarmupValues( B );
   

   
   // ---- set number-of x-nodes ---- //
   B->SetNumberOfXNodesPerMagnitude(8,xlim);


   // ---- scale description ---- //
   B->SetScale1Name( "Q" );
   B->SetScale2Name( "p_T" );
   //B->SetScale2Name( "pT_likeDijet" );
   
   // ---- number of scale nodes for mu ---- //
   B->SetNumberOfScaleNodesScale1( 5 );
   B->SetNumberOfScaleNodesScale2( 5 );

   
   // ---- Choose function for ScaleNode distances ---- //
   // --- fastNLO user: possibility to choose function which
   //     is used for binning of scale nodes.
   //     possible choices
   //        - LogLog Binning ( H(mu) = log(log(mu/0.25)) )
   //        - Linear         ( H(mu) = mu )
   B->InitLogLogScaleNode( A2 , scale1lo , scale1hi , 1 );	// choose function H for scale 1
   B->InitLogLogScaleNode( A2 , scale2lo , scale2hi , 2 );	// choose function H for scale 2

   
   // ---- final initializations ---- //
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
   //B->SetFuncMuForReference( Fct_x2y2ov2, Fct_x, 0 );
   //B->SetFuncMuForReference( Fct_x , Fct_x , 1 );
   //B->SetFuncMuForReference( Fct_y , Fct_y , 2 );

}

void UserDIS::DefineBinning(){
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
   const int nq2bins = 6;
   // 2)
   double q2bins[nq2bins+1] = {150,200,270,400,700,5000,15000};
   // 3)
   const int netbins[nq2bins] = {4, 4, 4, 4, 4, 4};
   
   // prepare 4)
   double etbins[nq2bins][6] = {
      {7,11,18,30,50},
      {7,11,18,30,50},
      {7,11,18,30,50},
      {7,11,18,30,50},
      {7,11,18,30,50},
      {7,11,18,30,50}
   };
   // 4)
   vector<double*> vetbins(nq2bins);
   for (unsigned int i=0;i<vetbins.size();i++) vetbins[i]=etbins[i]; 
   
   // ---- pass arrays to FnloTable and init bingrid ---- //
   table->GetBlockA2()-> InitBinning( nq2bins , q2bins , netbins , vetbins );
}


void UserDIS::GetWarmupValues( fnloBlockBNlojet* B ){
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
      
      // ---------- standard values: mu1 = Q2, mu2 = p_T
      //        940000000 contributions (!= events) in warm-up run
      xlim [0] = 4.63e-03 , scale1lo [0] =   12.2474 , scale1hi [0] =   14.1421 , scale2lo [0] =    7.0000 , scale2hi [0] =   11.0000;
      xlim [1] = 8.72e-03 , scale1lo [1] =   12.2474 , scale1hi [1] =   14.1421 , scale2lo [1] =   11.0000 , scale2hi [1] =   18.0000;
      xlim [2] = 1.74e-02 , scale1lo [2] =   12.2474 , scale1hi [2] =   14.1421 , scale2lo [2] =   18.0000 , scale2hi [2] =   30.0000;
      xlim [3] = 5.01e-02 , scale1lo [3] =   12.2474 , scale1hi [3] =   14.1421 , scale2lo [3] =   30.0000 , scale2hi [3] =   50.0000;
      xlim [4] = 5.32e-03 , scale1lo [4] =   14.1421 , scale1hi [4] =   16.4317 , scale2lo [4] =    7.0000 , scale2hi [4] =   11.0000;
      xlim [5] = 9.42e-03 , scale1lo [5] =   14.1421 , scale1hi [5] =   16.4317 , scale2lo [5] =   11.0000 , scale2hi [5] =   18.0000;
      xlim [6] = 1.82e-02 , scale1lo [6] =   14.1421 , scale1hi [6] =   16.4317 , scale2lo [6] =   18.0000 , scale2hi [6] =   30.0000;
      xlim [7] = 5.09e-02 , scale1lo [7] =   14.1421 , scale1hi [7] =   16.4317 , scale2lo [7] =   30.0000 , scale2hi [7] =   50.0000;
      xlim [8] = 6.32e-03 , scale1lo [8] =   16.4317 , scale1hi [8] =   20.0000 , scale2lo [8] =    7.0000 , scale2hi [8] =   11.0000;
      xlim [9] = 0.88e-02 , scale1lo [9] =   16.4317 , scale1hi [9] =   20.0000 , scale2lo [9] =   11.0000 , scale2hi [9] =   18.0000;
      xlim [10] = 2.01e-02 , scale1lo [10] =   16.4317 , scale1hi [10] =   20.0000 , scale2lo [10] =   18.0000 , scale2hi [10] =   30.0000;
      xlim [11] = 5.12e-02 , scale1lo [11] =   16.4317 , scale1hi [11] =   20.0000 , scale2lo [11] =   30.0000 , scale2hi [11] =   50.0000;
      xlim [12] = 8.18e-03 , scale1lo [12] =   20.0000 , scale1hi [12] =   26.4575 , scale2lo [12] =    7.0000 , scale2hi [12] =   11.0000;
      xlim [13] = 1.06e-02 , scale1lo [13] =   20.0000 , scale1hi [13] =   26.4575 , scale2lo [13] =   11.0000 , scale2hi [13] =   18.0000;
      xlim [14] = 2.11e-02 , scale1lo [14] =   20.0000 , scale1hi [14] =   26.4575 , scale2lo [14] =   18.0000 , scale2hi [14] =   30.0000;
      xlim [15] = 5.34e-02 , scale1lo [15] =   20.0000 , scale1hi [15] =   26.4575 , scale2lo [15] =   30.0000 , scale2hi [15] =   50.0000;
      xlim [16] = 1.02e-02 , scale1lo [16] =   26.4575 , scale1hi [16] =   70.7107 , scale2lo [16] =    7.0000 , scale2hi [16] =   11.0000;
      xlim [17] = 1.38e-02 , scale1lo [17] =   26.4575 , scale1hi [17] =   70.7107 , scale2lo [17] =   11.0000 , scale2hi [17] =   18.0000;
      xlim [18] = 2.54e-02 , scale1lo [18] =   26.4575 , scale1hi [18] =   70.7107 , scale2lo [18] =   18.0000 , scale2hi [18] =   30.0000;
      xlim [19] = 5.91e-02 , scale1lo [19] =   26.4575 , scale1hi [19] =   70.7107 , scale2lo [19] =   30.0000 , scale2hi [19] =   50.0000;
      xlim [20] = 7.13e-02 , scale1lo [20] =   70.7107 , scale1hi [20] =  122.4745 , scale2lo [20] =    7.0000 , scale2hi [20] =   11.0000;
      xlim [21] = 7.47e-02 , scale1lo [21] =   70.7107 , scale1hi [21] =  122.4745 , scale2lo [21] =   11.0000 , scale2hi [21] =   18.0000;
      xlim [22] = 8.64e-02 , scale1lo [22] =   70.7107 , scale1hi [22] =  122.4745 , scale2lo [22] =   18.0000 , scale2hi [22] =   30.0000;
      xlim [23] = 1.02e-01 , scale1lo [23] =   70.7107 , scale1hi [23] =  122.4745 , scale2lo [23] =   30.0000 , scale2hi [23] =   50.0000;

//       //        120000000 contributions (!= events) in warm-up run
//       //        mu1 = Q2, mu2 = ptLikeDijet
//       xlim [0] = 5.10e-03 , scale1lo [0] =   12.2475 , scale1hi [0] =   14.1421 , scale2lo [0] =    6.0193 , scale2hi [0] =   49.9500;
//       xlim [1] = 9.07e-03 , scale1lo [1] =   12.2475 , scale1hi [1] =   14.1421 , scale2lo [1] =    7.0049 , scale2hi [1] =   49.9736;
//       xlim [2] = 2.04e-02 , scale1lo [2] =   12.2475 , scale1hi [2] =   14.1421 , scale2lo [2] =    7.0076 , scale2hi [2] =   49.8985;
//       xlim [3] = 5.23e-02 , scale1lo [3] =   12.2475 , scale1hi [3] =   14.1421 , scale2lo [3] =    7.0317 , scale2hi [3] =   50.0000;
//       xlim [4] = 5.63e-03 , scale1lo [4] =   14.1421 , scale1hi [4] =   16.4317 , scale2lo [4] =    6.0029 , scale2hi [4] =   49.9253;
//       xlim [5] = 9.82e-03 , scale1lo [5] =   14.1421 , scale1hi [5] =   16.4317 , scale2lo [5] =    7.0013 , scale2hi [5] =   49.7937;
//       xlim [6] = 2.11e-02 , scale1lo [6] =   14.1421 , scale1hi [6] =   16.4317 , scale2lo [6] =    7.0090 , scale2hi [6] =   49.7407;
//       xlim [7] = 5.07e-02 , scale1lo [7] =   14.1421 , scale1hi [7] =   16.4317 , scale2lo [7] =    7.2057 , scale2hi [7] =   50.0000;
//       xlim [8] = 6.69e-03 , scale1lo [8] =   16.4317 , scale1hi [8] =   20.0000 , scale2lo [8] =    6.0196 , scale2hi [8] =   49.8707;
//       xlim [9] = 1.07e-02 , scale1lo [9] =   16.4317 , scale1hi [9] =   20.0000 , scale2lo [9] =    7.0086 , scale2hi [9] =   49.9098;
//       xlim [10] = 2.16e-02 , scale1lo [10] =   16.4317 , scale1hi [10] =   20.0000 , scale2lo [10] =    7.0212 , scale2hi [10] =   49.9484;
//       xlim [11] = 5.40e-02 , scale1lo [11] =   16.4317 , scale1hi [11] =   20.0000 , scale2lo [11] =    7.0652 , scale2hi [11] =   49.9999;
//       xlim [12] = 8.49e-03 , scale1lo [12] =   20.0000 , scale1hi [12] =   26.4575 , scale2lo [12] =    6.0182 , scale2hi [12] =   49.9368;
//       xlim [13] = 1.24e-02 , scale1lo [13] =   20.0000 , scale1hi [13] =   26.4575 , scale2lo [13] =    7.0019 , scale2hi [13] =   49.8723;
//       xlim [14] = 2.44e-02 , scale1lo [14] =   20.0000 , scale1hi [14] =   26.4575 , scale2lo [14] =    7.0174 , scale2hi [14] =   49.9897;
//       xlim [15] = 5.63e-02 , scale1lo [15] =   20.0000 , scale1hi [15] =   26.4575 , scale2lo [15] =    7.0597 , scale2hi [15] =   49.9999;
//       xlim [16] = 1.31e-02 , scale1lo [16] =   26.4575 , scale1hi [16] =   70.7106 , scale2lo [16] =    6.0064 , scale2hi [16] =   49.9353;
//       xlim [17] = 1.68e-02 , scale1lo [17] =   26.4575 , scale1hi [17] =   70.7106 , scale2lo [17] =    7.0001 , scale2hi [17] =   49.9587;
//       xlim [18] = 2.86e-02 , scale1lo [18] =   26.4575 , scale1hi [18] =   70.7107 , scale2lo [18] =    7.0023 , scale2hi [18] =   49.9097;
//       xlim [19] = 6.24e-02 , scale1lo [19] =   26.4575 , scale1hi [19] =   70.7107 , scale2lo [19] =    7.0093 , scale2hi [19] =   50.0000;
//       xlim [20] = 7.44e-02 , scale1lo [20] =   70.7107 , scale1hi [20] =  122.4745 , scale2lo [20] =    6.0142 , scale2hi [20] =   49.9383;
//       xlim [21] = 7.81e-02 , scale1lo [21] =   70.7108 , scale1hi [21] =  122.4745 , scale2lo [21] =    7.0000 , scale2hi [21] =   49.9456;
//       xlim [22] = 9.00e-02 , scale1lo [22] =   70.7107 , scale1hi [22] =  122.4744 , scale2lo [22] =    7.0001 , scale2hi [22] =   49.9741;
//       xlim [23] = 1.24e-01 , scale1lo [23] =   70.7107 , scale1hi [23] =  122.4744 , scale2lo [23] =    7.0080 , scale2hi [23] =   50.0000;

      // --------- fastNLO: Warm-Up run results (end)
   }
   else {
      printf("fastNLO: This is a warm-up run!\n");
      // --- fastNLO user: You can set the number of contributions
      //     after which the WarmUp values are printed
      B->SetWarmUpPrint(100000);		// default 10000000

      // safe initialziations
      for(int i=0;i<NObsBin;i++){ 
	 xlim[i] = 1.1e-07;
	 scale1lo[i] = 1.0, scale1hi[i]=9.9e10;
	 scale2lo[i] = 1.0, scale2hi[i]=9.9e10;
      }
   }
}




// --- fastNLO user: nothing further todo

unsigned int UserDIS::GetNj(){
   unsigned int nj = 0, nu = 0 ,nd = 0;
   inputfunc(nj,nu,nd);
   return nj;
}


double UserDIS::GetEcms(){
   double el = 0, eh = 0, q2min=0, q2max=0,ymin=0,ymax=0,xmin=0,xmax=0;
   psinput(NULL, el, eh, q2min, q2max, xmin, xmax, ymin, ymax);
   double ecms = 4*el*eh;
   return sqrt(ecms);
}

void UserDIS::boost_back_to_lab(const event_dis& p)
{ 
  double x = p[-1].T()/p[hadron(0)].T();
  double bz = (1.0 - x)/(1.0 + x);
  threevector<double> bVec = -((pbreit[-1] + pbreit[hadron(0)]).boostVector());
  lorentzvector<double> p0(pbreit[hadron(0)]);
  
  p0.boost(bVec);
  double phi = p0.phi(), theta = p0.theta();
  unsigned int njet = pjl.upper();
  
  for(unsigned int i = 1; i <= njet; i++) {
    pjl[i].boost(bVec);
    pjl[i].rotateZ(-phi);
    pjl[i].rotateY(-theta); 
    pjl[i].boost(0.0, 0.0, bz);
  }
}

void UserDIS::initfunc(unsigned int)
{
   // ---- Initialize event counters and set some defaults
   nevents = 0;
   if (nwrite==0) nwrite = 5000000;
   start_time = std::time(0);
}

void UserDIS::writetable(){
   table->OpenFileRewrite();
   table->WriteBlockA1();
   table->WriteBlockA2();
   for(int i=0;i< table->GetBlockA1()->GetNcontrib();i++){
      table->WriteBlockBDividebyN(i);
   }
   table->CloseFileWrite();
}

void UserDIS::end_of_event(){
   nevents += 1;
   //-------- store table
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
      printf ("No. events: %.2G writing table....",nevents);
      cout.flush();
      for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
         table->GetBlockB(k)->Nevt = (long long int)nevents;
      }
      writetable();
      printf("fastNLO: Old table written.\n");
      ftable->SetNumberOfEvents(nevents);
      ftable->WriteTable();    
      printf("fastNLO: New table written.\n");
      printf("done.\n");
   }
}

void UserDIS::phys_output(const std::basic_string<char>& __file_name, 
                          unsigned long __save, bool __txt) 
{
   // fastNLO v2.2
   InitFastNLO(__file_name);

   tablefilename.assign(__file_name.c_str());
   tablefilename += ".tab";
   
   //Determine whether we are running LO or NLO
   const char* const file = __file_name.c_str(); 
   if(strstr(file,"-born-")!=NULL){
      nlo = false;
   }else{
      if(strstr(file,"-nlo-")!=NULL){
         nlo = true;
      }else{
         printf("This module can only be run at Born level or at NLO.\n");
         exit(1);
      }
   }
   nwrite = __save;
   inittable();

}

//____________________ fastNLO v2.2 ____________________________
void UserDIS::InitFastNLO(const std::basic_string<char>& __file_name)
{
   // create table and read in steering...
   cout<<"\n ---------------------------------------------------------------\n"<<endl;
   ftable = new fastNLOCreate("fnh5001.str");

   // obtain relevant variables from nlojet
   ftable->SetEcms(UsefulNlojetTools::GetEcms());
   ftable->SetLoOrder(UsefulNlojetTools::GetLoOrder());
   ftable->SetOrderOfAlphasOfCalculation(UsefulNlojetTools::GetOrderOfRun(__file_name));

   // set filename, which is specified through command line
   string tabFilename = __file_name.c_str();
   tabFilename += "_neu.tab";
   ftable->SetFilename(tabFilename);

   // give information to hb.
   //ftable->Print();
}
