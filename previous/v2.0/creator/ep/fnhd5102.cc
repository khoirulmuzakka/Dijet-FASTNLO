//
// fastNLO v2.1 author code for fnhd5102:
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
//#include "h1jet.h"


class UserDIS : public basic_user_set<user1d_dis, user1h_dis> 
{
public:
  //   init and user function
   void initfunc(unsigned int);
   void userfunc(const event_dis&, const amplitude_dis&);
   void lab_to_hcms(event_dis&);  

   virtual void end_of_event();  
   virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);
  
private:
   //   pdf
   pdf_cteq6dis pdf;
   pdf_dis_dummy dummypdf;
   //pdf_dis_flat_dummy flatpdf;

   // algorithms
   kT_clus_long jetclus;
   //h1jet h1_jet;
   
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
  nd = 3U; // just 2 here !
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
  q2min = 4.0;    // GeV^2
  q2max = 80.0;  // GeV^2 

  //   xB cuts
  xmin = 0.0;
  xmax = 1.0;

  //   y cuts
  ymin = 0.1;
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
   // ---- copy the momenta and boost to the breit frame ----
   //  pbreit = p; lab_to_breit(pbreit);
   pbreit = p; lab_to_hcms(pbreit);



   //----- do the cluster analysis in the breit frame-----
   // h1jet pjb = h1_jet(pbreit);
   jetclus.set_up(1,false,2);
   jetclus(pbreit); jetclus.incl(pjb, jet);
   //unsigned int njorg = pjb.upper();

   // ---- jet structure in laboratory frame -----
   pjl = pjb; boost_back_to_lab(p);
   //boost_back_to_hell(p,pbreit);

   unsigned int nj = pbreit.upper(); 
   unsigned int numj = pjb.upper(); 
   //    if ( nj != njorg ) {
   //       cout << "nj = "<<nj<<"\tnjorg = " << njorg <<"\tnumj = " << numj << endl;
   //    }
   nj = (unsigned int)numj;
   
   // ---- Momentum transfer -----
   const double Q2 = -((p[-1] - p[-2]).mag2());
   const double Q = sqrt(Q2);
   
   // ---- el-mag coupling --- 
   const double alem = xalpha_em(Q2);


   int l = p.lower(), u = p.upper();

   // ---- cuts -----
   const double pt1min = 5.5;  
   const double pt2min = 4.;
   const double etamin = -3.0; 
   const double etamax = 0.0;



   // ---- apply jet cuts and fill jets array----
   int hjetidx = -1;	// hard jet index
   int sjetidx = -1;	// soft jet index
   
   
   double tmpptmax = 0.;

   // look for hard jet
   for ( unsigned int j = 1; j <= nj; j++){
      double pt = pjb[j].perp(); 
      double px = pjb[j].X(); 
      double py = pjb[j].Y(); 
      double eta = pjl[j].prapidity();
      
      if ((pt > pt1min) && (eta > etamin) && (eta < etamax) && ((pt > tmpptmax))){
	 tmpptmax = pt;
	 hjetidx = j;
      }
   }
    
   // look for soft jet
   double tmpptmax2 = -10.;
   if (hjetidx > -1){
      for ( unsigned int j = 1; j <= nj; j++){
	 double pt = pjb[j].perp(); 
	 double px = pjb[j].X(); 
	 double py = pjb[j].Y(); 
	 double eta = pjl[j].prapidity();
	 //double etac = pjl[hjetidx].prapidity(); 

	 if ((pt > pt2min) && (eta > etamin) && (eta < etamax)&&(j != hjetidx)&&( pt > tmpptmax2)){
	    tmpptmax2 = pt;
	    sjetidx = j;
	 }
      }
   }
    

   // ---- I Think it is a dijet event ---- //
   // hard jet index hjetidx no longer -1
   // soft jet index hjetidx no longer -1
   if (((hjetidx > 0)&&(sjetidx > 0))){

      double hardpt = pjb[hjetidx].perp();
      double softpt = pjb[sjetidx].perp();
      double hardptlab = pjl[hjetidx].perp();
      double softptlab = pjl[sjetidx].perp();
      
      double meanpt = 0.5*(hardpt + softpt);

      //       double xp    = 0.1;//read_steer::get_xpom();
      //       double logxp = log10(xp);
      //       double hardeta = pjb[hjetidx].prapidity();
      //       double softeta = pjb[sjetidx].prapidity();
      //       double hardetalab = pjl[hjetidx].prapidity();
      //       double softetalab = pjl[sjetidx].prapidity();
      //       double deltaeta = abs(hardeta - softeta);
      //       double eta2lab = pjl[sjetidx].prapidity();
      //       double dphi = abs(pjb[hjetidx].phi() - pjb[sjetidx].phi());
      //       const double pi = 3.14159;
      //       if (dphi > pi) dphi = 2*pi - dphi;
      

       lorentzvector<double> jetsys(pjb[hjetidx]+pjb[sjetidx]);
      
       //        // 4-momentum transfer squared
       //        lorentzvector<double> q(p[-1]-p[-2]);
       //        double y=1.-(p[-2].T()/p[-1].T())*sin(p[-2].theta()/2.)*sin(p[-2].theta()/2.);
       //        double beta = 0.5*Q2/(p[hadron(0)]*q);
       //        double zp = beta*(1+ jetsys.mag2()/Q2);
       double M12 = sqrt(jetsys.mag2());


       // ------------ FastNLO starts here ------------- //
       // ---- define observables ----
       double val1 = hardpt;
       //double val1 = Q2;
       //double val1 = y;
       double val2 = Q2;
       double mu1 = Q;
       double mu2 = meanpt;

       // --- fill fastNLO arrays
       //     Scales must be in GeV or dimensionless (at least one in GeV)
       //     Values must be in same dimension as your binning is defined
       FillEvent( val1 , val2 , mu1 , mu2 , p , amp , alem );

   }

}


void UserDIS::FillEvent( double val1 , double val2 , double mu1, double mu2 , const event_dis& p, const nlo::amplitude_dis& amp , double alem ){
   // ---- Fill FastNLO Array ---- //
   // --- fastNLO user: usually nothing to do

   // ---- get x-value ---- //
   double x = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);
   double Q2 = -((p[-1] - p[-2]).mag2());
   double Q = sqrt(Q2);
   if ( fabs(mu1/Q-1.)>1.e-8 ) {
      printf("UserDIS::FillEvent. Sorry, but scale1 must be 'Q'\n");
      printf("Q = %9.4f, Q2 = %9.4f, mu1 = %9.4f\n",Q,Q2,mu1);
      exit(1);
   }
   int obsbin = table->GetBlockA2()->GetNumDiff() == 1 ? 
      table->GetBlockA2()->GetBinNumber( val1 ):
      table->GetBlockA2()->GetBinNumber( val1 , val2  );

   //cout<<"obsbin =  "<<obsbin <<"\tval1 = " << val1 << "\tval2 = " << val2<< endl;
   
   //---------- fill fastNLO arrays
   if ( obsbin >= 0 ) {
      // --- prefactor: divide by binwidth - and multiply with alpha_elm^2
      double prefactor = 1./table->GetBlockA2()->BinSize[obsbin] * alem*alem; 
      for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
	 ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventDISMuVar(obsbin,x,mu1,mu2,amp,dummypdf,pdf,prefactor);
      }
   } // - end: fill fastNLO array
}



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
   table->GetBlockA1()->SetScenName("fnhd5102");

   // --- fastNLO user: up to 20 strings and any number of lines to describe the scenario
   A2->ScDescript.push_back("H1 Diffractive Dijets - d2sigma/dpt1 [pb]");
   A2->ScDescript.push_back("H1 Collaboration");
   A2->ScDescript.push_back("Diffractive dijet cross sections in DIS");
   A2->ScDescript.push_back("kT-algorithm R=1, Low Q2, Hera-II");
   A2->ScDescript.push_back("test scenario");
   //A2->ScDescript.push_back("Fig. 6a, Tab. 4a");
   A2->ScDescript.push_back("provided by:");
   A2->ScDescript.push_back("fastNLO_2.1.0");
   A2->ScDescript.push_back("If you use this table, please cite:");
   A2->ScDescript.push_back("  D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch, arXiv:1109.1310");

   
   // --- fastNLO user: Give information about your measurement
   A2->SetIpublunits( 12 );			// --- fastNLO user: set cross section units (negative power of ten), e.g. 'pb' -> 12.
   A2->SetNumDiffBin( 1 );			// --- fastNLO user: No of dimensions in which observable is binned
   bool IsDiffBin = false;			// --- fastNLO user: Are publication units divided by this variable?

   A2->SetDimLabel("p_T,1", 1 ,IsDiffBin );// --- fastNLO user: label of 1st dimension
   //A2->SetDimLabel("Q^2", 1 , IsDiffBin );// --- fastNLO user: label of 1st dimension


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
   B->SetNumberOfXNodesPerMagnitude(12,xlim);


   // ---- scale description ---- //
   B->SetScale1Name( "Q" );
   B->SetScale2Name( "<p_T>" );

   
   // ---- number of scale nodes for mu ---- //
   // please watch the warm-up values to see the covered interval
   B->SetNumberOfScaleNodesScale1( 5 );
   B->SetNumberOfScaleNodesScale2( 4 );

   
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

   //  __dQ2__
   //    // 1)
   //    const int nq2bins =
   //    // 2)
   //    double 
   //q2bins[] = {4., 6., 10., 15., 20., 30., 40., 60., 80.}; 
   //  __dpT__
   // 1)
   const int nq2bins = 5;
   // 2)
   double q2bins[nq2bins+1] =     {5.5, 6.5, 7.5, 9., 11., 13.5};



   //    // 4)
   //    vector<double*> vetbins(nq2bins);
   //    for (unsigned int i=0;i<vetbins.size();i++) vetbins[i]=etbins[i]; 
   
   // ---- pass arrays to FnloTable and init bingrid ---- //
   table->GetBlockA2()->InitBinning( nq2bins , q2bins );
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
      //         97000000 contributions (!= events) in warm-up run
      xlim [0] = 1.64e-03 , scale1lo [0] =    2.0000 , scale1hi [0] =    8.9443 , scale2lo [0] =    4.750 , scale2hi [0] =    6.5000;
      xlim [1] = 2.27e-03 , scale1lo [1] =    2.0000 , scale1hi [1] =    8.9443 , scale2lo [1] =    5.250 , scale2hi [1] =    7.5000;
      xlim [2] = 3.02e-03 , scale1lo [2] =    2.0000 , scale1hi [2] =    8.9443 , scale2lo [2] =    5.750 , scale2hi [2] =    9.0000;
      xlim [3] = 4.29e-03 , scale1lo [3] =    2.0000 , scale1hi [3] =    8.9443 , scale2lo [3] =    6.500 , scale2hi [3] =   11.0000;
      xlim [4] = 6.33e-03 , scale1lo [4] =    2.0000 , scale1hi [4] =    8.9443 , scale2lo [4] =    7.500 , scale2hi [4] =   13.5000;

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
      printf("done.\n");
   }
}

void UserDIS::phys_output(const std::basic_string<char>& __file_name, 
                          unsigned long __save, bool __txt) 
{
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


void UserDIS::lab_to_hcms(event_dis& p)
{
   lorentzvector<double> q(p[-1]-p[-2]);
   //double xb = -0.5*q.mag2()/(p[0]*q);

   double tmp_p_in[4]={0.,0.,920.,920.};

   lorentzvector<double> p_in= lorentzvector<double>(tmp_p_in[0],tmp_p_in[1],tmp_p_in[2],tmp_p_in[3]); 
   lorentzvector<double> p_in_hcms;

   threevector<double> bVec = -((q + p_in).boostVector());
   //    threevector<double> bVec = -((q + p[hadron(0)]).boostVector());
   int i, low = p.lower(), up = p.upper();
    

   // p[hadron(0)] ... outgoing parton
   // p[0]   ......... xp*p... incoming

   p_in_hcms = p_in;
   p_in_hcms.boost(bVec);
   //    p_in_hcms.rotateZ(-phi);
   //     p_in_hcms.rotateY(-theta);


   p[hadron(0)].boost(bVec);
   for(i = low; i <= up; i++) 
      p[i].boost(bVec);
    
   //    double phi = p[0].phi(), theta = p[0].theta();
   double phi = p_in_hcms.phi(), theta = p_in_hcms.theta();

    
   p[hadron(0)].rotateZ(-phi);
   p[hadron(0)].rotateY(-theta);  
    
   for(i = low; i <= up; i++) {
      p[i].rotateZ(-phi);
      p[i].rotateY(-theta);
   }
   p_in_hcms.rotateZ(-phi);
   p_in_hcms.rotateY(-theta );
   //    p_in_hcms.rotateZ(-phi);
   //     p_in_hcms.rotateY(-theta);

   //    cout << " p_in_hcms: " << p_in_hcms << endl;
   //     cout << " beam elec: " << p[-1] << endl;
   //     cout << " scat elec: " << p[-2] << endl;

}


