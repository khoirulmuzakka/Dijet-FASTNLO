//
// fastNLO v2.1 author code for fnhd0001:
//   this is the testing module for Diffracive jet production in DIS
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
   void boost_back_to_hell(const event_dis&,event_dis&);

   
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
  q2min = 4.0;    // GeV^2
  q2max = 110.0;  // GeV^2 

  //   xB cuts
  xmin = 0.0;
  xmax = 1.0;

  //   y cuts
  ymin = 0.05;
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

   // ->>>> No BOOST FOR FASTJET implemented !!!

   // ---- jet structure in laboratory frame -----
   pjl = pjb; boost_back_to_hell(p,pbreit);

   unsigned int nj = pbreit.upper(); 
   unsigned int numj = pjb.upper(); 
   //    if ( nj != njorg ) {
   //       cout << "nj = "<<nj<<"\tnjorg = " << njorg <<"\tnumj = " << numj << endl;
   //    }
   
   // ---- Momentum transfer -----
   const double Q2 = -((p[-1] - p[-2]).mag2());
   const double Q = sqrt(Q2);
   
   // ---- el-mag coupling --- 
   const double alem = xalpha_em(Q2);


   int l = p.lower(), u = p.upper();
   int pup = pbreit.upper();


   // ---- cuts -----
   const int ilabcuts  = 1;
   const int ietasort  = 0;
   const int noetacuts = 0;

   const double pt1min = 5.;  
   const double pt2min = 4.;
   const double etamin = -1.; 
   const double etamax = 2.5;

//    int iphas = read_steer::get_phas();
//    if (iphas == 2){
//       pt1min = 5.5;
//       //    etamax = 0.0;
//       etamax = 2.0;
//    }

 


   // ---- apply jet cuts and fill jets array----
   int hjetidx = -1;	// hard jet index
   int sjetidx = -1;	// soft jet index
   
   nj = (unsigned int)numj;
   
   double tmpptmax = 0.;
   if (ietasort)   tmpptmax = 10.;

   /* 
   nlo::weight_dis wt = amp(dummypdf,Q2,Q2,1.0);
   nlo::weight_dis wtf = amp(flatpdf,Q2,Q2,1.0);
   double x = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);
   double yps=1.-(p[-2].T()/p[-1].T())*sin(p[-2].theta()/2.)*sin(p[-2].theta()/2.);
   double s = pow(GetEcms(),2);
   double xsy = x*yps*s;
   double Q21sy = Q2/s/yps;
   printf(" --- Q2 = %9.3f, nj = %d,  y = %6.4f\n",Q2,nj,yps);
   printf("     x  = %7.5f\n",x); 
   printf("     s  = %8.5f, xsy = %9.3f, x=Q2/sy = %7.5f\n",s,xsy,Q21sy);
   printf("     pt1 = %6.3f, pt2 = %6.3f\n",pjb[1].perp(),pjb[2].perp());
   //    printf("     p[1].perp() = %6.3f, p[2].perp() = %6.3f\n",p[1].perp(), p[2].perp());
   //    printf("     p[1].mag2() = %6.3f, p[2].mag2() = %6.3f\n",p[1].mag2(), p[2].mag2());
   //    printf("     p[1].T()    = %6.3f, p[2].T()    = %6.3f\n",p[1].T(), p[2].T());
   //    printf("     p[1].prap() = %6.3f, p[2].prap() = %6.3f\n",p[1].prapidity(), p[2].prapidity());
   //    printf("     pb[1].perp() = %6.3f, pb[2].perp() = %6.3f\n",pbreit[1].perp(), pbreit[2].perp());
   //    printf("     pb[1].mag2() = %6.3f, pb[2].mag2() = %6.3f\n",pbreit[1].mag2(), pbreit[2].mag2());
   //    printf("     pb[1].T()    = %6.3f, pb[2].T()    = %6.3f\n",pbreit[1].T(), pbreit[2].T());
   //    printf("     pb[1].prap() = %6.3f, pb[2].prap() = %6.3f\n",pbreit[1].prapidity(), pbreit[2].prapidity());
   printf("     pjb[1].perp() = %6.3f, pjb[2].perp() = %6.3f\n",pjb[1].perp(), pjb[2].perp());
   //    printf("     pjb[1].mag2() = %6.3f, pjb[2].mag2() = %6.3f\n",pjb[1].mag2(), pjb[2].mag2());
   printf("     pjb[1].T()    = %6.3f, pjb[2].T()    = %6.3f\n",pjb[1].T(), pjb[2].T());
   printf("     pjb[1].prap() = %6.3f, pjb[2].prap() = %6.3f\n",pjb[1].prapidity(), pjb[2].prapidity());

   printf("     pjl[1].perp() = %6.3f, pjl[2].perp() = %6.3f\n",pjl[1].perp(), pjl[2].perp());
   //    printf("     pjl[1].mag2() = %6.3f, pjl[2].mag2() = %6.3f\n",pjl[1].mag2(), pjl[2].mag2());
   printf("     pjl[1].T()    = %6.3f, pjl[2].T()    = %6.3f\n",pjl[1].T(), pjl[2].T());
   printf("     pjl[1].prap() = %6.3f, pjl[2].prap() = %6.3f\n",pjl[1].prapidity(), pjl[2].prapidity());

   printf("     wf = [ %9.3e , %9.3e , %9.3e ],\n     w  = [ %9.3e , %9.3e , %9.3e ]\n",
	  wtf[0],wtf[1],wtf[2],wt[0],wt[1],wt[2]);

   */

   // look for hard jet
   for ( unsigned int j = 1; j <= nj; j++){
      double pt = pjb[j].perp(); 
      double px = pjb[j].X(); 
      double py = pjb[j].Y(); 
      double eta = pjb[j].prapidity();
      
      if (ilabcuts)
	 eta = pjl[j].prapidity();
      
      if (ietasort){
	 if ((pt > pt1min) && (eta > etamin) && (eta < etamax) && ((eta <  tmpptmax))){  // etasort
	    tmpptmax = eta;
	    hjetidx = j;
	 }
      }else{
	 if ((pt > pt1min) && (eta > etamin) && (eta < etamax) && ((pt > tmpptmax))){
	    tmpptmax = pt;
	    hjetidx = j;
	 }
      }
   }
    
   // look for soft jet
   double tmpptmax2 = -10.;
   if (hjetidx > -1){
      for ( unsigned int j = 1; j <= nj; j++){
	 double pt = pjb[j].perp(); 
	 double eta = pjb[j].prapidity();
	 double etac = pjb[hjetidx].prapidity(); 
	 double px = pjb[j].X(); 
	 double py = pjb[j].Y(); 

	 if (ilabcuts){
	    eta = pjl[j].prapidity();
	    etac = pjl[hjetidx].prapidity(); 
	 }

	 if (ietasort){
	    if ((pt > pt2min) && (eta > etamin) && (eta < etamax)&&(j != hjetidx)&&( eta > tmpptmax2)){ // etasort
	       sjetidx = j;
	       tmpptmax2 = eta;
	    }
	 }
	 else{
	    if ((pt > pt2min) && (eta > etamin) && (eta < etamax)&&(j != hjetidx)&&( pt > tmpptmax2)){
	       tmpptmax2 = pt;
	       sjetidx = j;
	    }
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

      double q2pt2 = Q2 + meanpt*meanpt;

      //       double xp    = 0.1;//read_steer::get_xpom();
      //       double logxp = log10(xp);
      double hardeta = pjb[hjetidx].prapidity();
      double softeta = pjb[sjetidx].prapidity();
      double hardetalab = pjl[hjetidx].prapidity();
      double softetalab = pjl[sjetidx].prapidity();
      double deltaeta = abs(hardeta - softeta);
      double eta2lab = pjl[sjetidx].prapidity();
      double dphi = abs(pjb[hjetidx].phi() - pjb[sjetidx].phi());
      const double pi = 3.14159;
      if (dphi > pi) dphi = 2*pi - dphi;
      
      //       if (pr2){
      // 	 cout << " hjetidx: " << hjetidx << " sjetidx: " << sjetidx << endl;
      // 	 cout << " hardetalab: " << hardetalab << " hardptlab: " << hardptlab << endl;
      // 	 cout << " softetalab: " << softetalab << " softptlab: " << softptlab << endl;
      // 	 cout << " hardeta: " << hardeta << " hardpt: " << hardpt << endl;
      // 	 cout << " softeta: " << softeta << " softpt: " << softpt << endl;
      //       } 

       lorentzvector<double> jetsys(pjb[hjetidx]+pjb[sjetidx]);
      
       // 4-momentum transfer squared
       lorentzvector<double> q(p[-1]-p[-2]);
       double y=1.-(p[-2].T()/p[-1].T())*sin(p[-2].theta()/2.)*sin(p[-2].theta()/2.);
       double beta = 0.5*Q2/(p[hadron(0)]*q);
       double zp = beta*(1+ jetsys.mag2()/Q2);
       double M12 = sqrt(jetsys.mag2());

       //        //----- hard scales -----
       //        const double mur = meanpt*meanpt + Q2; // this hsould be correct scale
       //        const double muf = mur;
       //weight_dis wt = amp(&pdf, mur, muf, 389385.730*alem*alem);
       //weight_dis wt3 = amp(&pdfB, mur, muf, 389385.730*alem*alem); // other pdf
       
       //        __fill_hist(1, 0.5, wt);
       //        __fill_hist(2, hardpt, wt); // <- this is what we want!
       //        __fill_hist(3, softpt, wt);
       //        __fill_hist(4, zp, wt);
       //        __fill_hist(5, logxp, wt);
       //        __fill_hist(6, y, wt);
       //        __fill_hist(7, deltaeta, wt);
       //        __fill_hist(8, Q2, wt);
       //        __fill_hist(9, q2pt2, wt);
       //        __fill_hist(19, M12, wt);

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
   table->GetBlockA1()->SetScenName("fnhd1002");

   // --- fastNLO user: up to 20 strings and any number of lines to describe the scenario
   A2->ScDescript.push_back("H1 Diffractive Dijets - d2sigma/dQ2 [pb/GeV2]");
   A2->ScDescript.push_back("H1 Collaboration");
   A2->ScDescript.push_back("Diffractive dijet cross sections in DIS");
   A2->ScDescript.push_back("kT-algorithm R=1, Low Q2, Hera-II");
   A2->ScDescript.push_back("DESY 11-166");
   A2->ScDescript.push_back("Fig. 6a, Tab. 4a");
   //A2->ScDescript.push_back("Fig. 5a");
   A2->ScDescript.push_back("provided by:");
   A2->ScDescript.push_back("fastNLO_2.1.0");
   A2->ScDescript.push_back("If you use this table, please cite:");
   A2->ScDescript.push_back("  D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch, arXiv:1109.1310");

   
   // --- fastNLO user: Give information about your measurement
   A2->SetIpublunits( 12 );			// --- fastNLO user: set cross section units (negative power of ten), e.g. 'pb' -> 12.
   A2->SetNumDiffBin( 1 );			// --- fastNLO user: No of dimensions in which observable is binned
   bool IsDiffBin = true;			// --- fastNLO user: Are publication units divided by this variable?

   A2->SetDimLabel("p*_T,1", 1 ,IsDiffBin );// --- fastNLO user: label of 1st dimension
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
   B->SetNumberOfScaleNodesScale1( 4 );
   B->SetNumberOfScaleNodesScale2( 16 );

   
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
   //    const int nq2bins = 5;
   //    // 2)
   //    double q2bins[nq2bins+1] =       { 4.0 , 8.0 , 16.0 , 32.0 , 60.0 , 110.0 }; // dQ2

   //  __dpT1*__
   // 1)
   const int nq2bins = 4;
   // 2)
   double q2bins[nq2bins+1] =     {5.0, 6.5, 8.5, 12.0 , 150}; // dpT,1
   
   //    // 3)
   //    const int netbins[nq2bins] = {5};
   //    // prepare 4)
   //    double etbins[nq2bins][6] = {
   //       //{5.0, 6.5, 8.5, 12.0 } // dP*T,1
   //       { 4.0 , 8.0 , 16.0 , 32.0 , 60.0 , 110.0 } // dQ2
   //    };


   //   // kind of y
   //    // 3)
   //    const int netbins[nq2bins] = {7};
   //    // prepare 4)
   //    double etbins[nq2bins][18] = {
   //       //{5.0, 6.5, 8.5, 12.0 } // dP*T,1
   //       { 0. , 0.1 , 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 }
   //    };



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

//       //----------------------------
//       //    dy
//       //         21900000 contributions (!= events) in warm-up run
//       xlim [0] = 1.04e-02 , scale1lo [0] =    2.0000 , scale1hi [0] =   10.4881 , scale2lo [0] =    4.5012 , scale2hi [0] =   33.0318;
//       xlim [1] = 5.05e-03 , scale1lo [1] =    2.0000 , scale1hi [1] =   10.4881 , scale2lo [1] =    4.5014 , scale2hi [1] =   66.3293;
//       xlim [2] = 3.47e-03 , scale1lo [2] =    2.0000 , scale1hi [2] =   10.4881 , scale2lo [2] =    4.5061 , scale2hi [2] =   87.2122;
//       xlim [3] = 2.62e-03 , scale1lo [3] =    2.0000 , scale1hi [3] =   10.4880 , scale2lo [3] =    4.5024 , scale2hi [3] =  100.5176;
//       xlim [4] = 2.15e-03 , scale1lo [4] =    2.0000 , scale1hi [4] =   10.4880 , scale2lo [4] =    4.5050 , scale2hi [4] =  112.2086;
//       xlim [5] = 2.11e-03 , scale1lo [5] =    2.0000 , scale1hi [5] =   10.4881 , scale2lo [5] =    4.5174 , scale2hi [5] =  123.2764;
//       xlim [6] = 2.15e-03 , scale1lo [6] =    2.0000 , scale1hi [6] =   10.4880 , scale2lo [6] =    4.5257 , scale2hi [6] =  132.7977;
//       //----------------------------


//       // -------------------------
//       //   dQ2
//       // 9700000 contributions (!= events) in warm-up run
//       xlim [0] = 2.10e-03 , scale1lo [0] =    2.0000 , scale1hi [0] =    2.8284 , scale2lo [0] =    4.50 , scale2hi [0] =   133.;
//       xlim [1] = 2.13e-03 , scale1lo [1] =    2.8284 , scale1hi [1] =    4.0000 , scale2lo [1] =    4.50 , scale2hi [1] =   133.;
//       xlim [2] = 2.28e-03 , scale1lo [2] =    4.0000 , scale1hi [2] =    5.6569 , scale2lo [2] =    4.50 , scale2hi [2] =   133.;
//       xlim [3] = 2.53e-03 , scale1lo [3] =    5.6569 , scale1hi [3] =    7.7460 , scale2lo [3] =    4.50 , scale2hi [3] =   133.;
//       xlim [4] = 2.96e-03 , scale1lo [4] =    7.7460 , scale1hi [4] =   10.4881 , scale2lo [4] =    4.50 , scale2hi [4] =   133.;
//       //----------------------------


//       //----------------------------
//       // dP*T,1
//       //         10000000 contributions (!= events) in warm-up run
//       xlim [0] = 2.07e-03 , scale1lo [0] =    2.0000 , scale1hi [0] =   10.4881 , scale2lo [0] =    4.5020 , scale2hi [0] =    6.5000;
//       xlim [1] = 2.69e-03 , scale1lo [1] =    2.0000 , scale1hi [1] =   10.4881 , scale2lo [1] =    5.2525 , scale2hi [1] =    8.5000;
//       xlim [2] = 4.13e-03 , scale1lo [2] =    2.0000 , scale1hi [2] =   10.4881 , scale2lo [2] =    6.2531 , scale2hi [2] =   12.0000;
//       //----------------------------
      
      
     //----------------------------
     // dP*T,1
     // 37000000 contributions (!= events) in warm-up run
     xlim [0] = 2.07e-03 , scale1lo [0] =    2.0000 , scale1hi [0] =   10.4881 , scale2lo [0] =    5.0000 , scale2hi [0] =    6.5000;
     xlim [1] = 2.69e-03 , scale1lo [1] =    2.0000 , scale1hi [1] =   10.4881 , scale2lo [1] =    6.5000 , scale2hi [1] =    8.5000;
     xlim [2] = 4.13e-03 , scale1lo [2] =    2.0000 , scale1hi [2] =   10.4881 , scale2lo [2] =    8.5000 , scale2hi [2] =   12.0000;
     xlim [3] = 8.20e-03 , scale1lo [3] =    2.0000 , scale1hi [3] =   10.4881 , scale2lo [3] =   12.0000 , scale2hi [3] =  133.2620;
     //       //----------------------------
      
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


void UserDIS::boost_back_to_hell(const event_dis& p1, event_dis& p2) // p1 = initial system, p2 = pbreit
{ 

  double pi = 3.14159265;

  double tmp_p_in[4]={0.,0.,920.,920.000};
  double tmp_e_in[4]={0.,0.,-27.6,27.6};
  double rest_p_lab[4]={0.,0.,0.,1.};
  double test_p_lab[4]={1.,2.,3.,4.};
  double xzlab[4]={1.,0.,1.,sqrt(2.)};

  lorentzvector<double> p_in= lorentzvector<double>(tmp_p_in[0],tmp_p_in[1],tmp_p_in[2],tmp_p_in[3]); 
  lorentzvector<double> e_in= lorentzvector<double>(tmp_e_in[0],tmp_e_in[1],tmp_e_in[2],tmp_e_in[3]); 
  lorentzvector<double> p_rest= lorentzvector<double>(rest_p_lab[0],rest_p_lab[1],rest_p_lab[2],rest_p_lab[3]); 
  lorentzvector<double> p_test= lorentzvector<double>(test_p_lab[0],test_p_lab[1],test_p_lab[2],test_p_lab[3]); 
  lorentzvector<double> xz_lab= lorentzvector<double>(xzlab[0],xzlab[1],xzlab[2],xzlab[3]); 
  lorentzvector<double> p_rest_hcms;
  lorentzvector<double> scat_elec;

  scat_elec = p1[-2];

  //e_in = p1[-1];

  lorentzvector<double> q(p1[-1]-p1[-2]);
  threevector<double> bVec2 = -((q + p_in).boostVector());
  lorentzvector<double> p_in_hcms;
  lorentzvector<double> e_in_hcms;
  lorentzvector<double> xz_hcms;
  lorentzvector<double> p_test_hcms;
  lorentzvector<double> scat_elec_in_hcms;
  lorentzvector<double> p_in_hell;

  p_in_hcms = p_in;
  p_in_hcms.boost(bVec2); 

//   cout << endl;
//   cout << " intermediate p_in_hcms: " << p_in_hcms << endl;

  double fPhi   = p_in_hcms.phi();
  double fTheta = p_in_hcms.theta();
 
  //  cout << " fPhi: " << fPhi*180./pi << "  fTheta: "<< fTheta*180./pi << endl;

  scat_elec_in_hcms = scat_elec;
  scat_elec_in_hcms.boost(bVec2);
  
  scat_elec_in_hcms.rotateZ(-fPhi);
  //  scat_elec_in_hcms.rotateY(pi-fTheta);
  scat_elec_in_hcms.rotateY(-fTheta);

  double fPhi2 = scat_elec_in_hcms.phi();

  p_rest_hcms = p_rest;
  p_rest_hcms.boost(bVec2);
  p_rest_hcms.rotateZ(-fPhi);
//   p_rest_hcms.rotateY(-fTheta + pi);
//   p_rest_hcms.rotateZ(-fPhi2);
  p_rest_hcms.rotateY(-fTheta);

  p_test_hcms = p_test;
  p_test_hcms.boost(bVec2);
  p_test_hcms.rotateZ(-fPhi);
//   p_test_hcms.rotateY(-fTheta + pi);
//   p_test_hcms.rotateZ(-fPhi2);
  p_test_hcms.rotateY(-fTheta);

  xz_hcms = xz_lab;
  xz_hcms.boost(bVec2);
  xz_hcms.rotateZ(-fPhi);
//   xz_hcms.rotateY(-fTheta + pi);
//   xz_hcms.rotateZ(-fPhi2);
  xz_hcms.rotateY(-fTheta );

  e_in_hcms = e_in;
  e_in_hcms.boost(bVec2);
  e_in_hcms.rotateZ(-fPhi);
 //  e_in_hcms.rotateY(-fTheta + pi);
//   e_in_hcms.rotateZ(-fPhi2);
  e_in_hcms.rotateY(-fTheta);


  p_in_hcms.rotateZ(-fPhi);
  p_in_hcms.rotateY(-fTheta );
//   p_in_hcms.rotateZ(-fPhi2);

  // -------------------- now things boosted to HCMS ----------------------------
  // ------------------------ so boost back -------------------------------------

  threevector<double> bVec = -((p_rest_hcms).boostVector());

  lorentzvector<double> xz_relab;
  lorentzvector<double> p_in_relab;
  lorentzvector<double> e_in_relab;
  lorentzvector<double> scat_elec_in_relab;
  lorentzvector<double> p_test_relab;
  lorentzvector<double> p_rest_relab;

  p_in_relab = p_in_hcms;
  p_in_relab.boost(bVec);

  //  cout << " intermediate P_relab: " << p_in_relab << endl;

  double fPhiLab   = p_in_relab.phi();
  double fThetaLab = p_in_relab.theta();

  //  cout << " fPhiLab: " << fPhiLab*180./pi << "  fThetaLab: "<< fThetaLab*180./pi << endl;


  xz_relab = xz_hcms;
  xz_relab.boost(bVec);

  xz_relab.rotateZ(-fPhiLab);
  xz_relab.rotateY(-fThetaLab);

  double fPhiLab2 = xz_relab.phi();

  xz_relab.rotateZ(-fPhiLab2);

  p_test_relab = p_test_hcms;
  p_test_relab.boost(bVec);
  p_test_relab.rotateZ(-fPhiLab);
  p_test_relab.rotateY(-fThetaLab );
    p_test_relab.rotateZ(-fPhiLab2);

  p_in_relab = p_in_hcms;
  p_in_relab.boost(bVec); 
  p_in_relab.rotateZ(-fPhiLab);
  p_in_relab.rotateY(-fThetaLab  );
    p_in_relab.rotateZ(-fPhiLab2);

  e_in_relab = e_in_hcms;
  e_in_relab.boost(bVec); 
  e_in_relab.rotateZ(-fPhiLab);
  e_in_relab.rotateY(-fThetaLab);
    e_in_relab.rotateZ(-fPhiLab2);

  scat_elec_in_relab = scat_elec_in_hcms;
  scat_elec_in_relab.boost(bVec); 
  scat_elec_in_relab.rotateZ(-fPhiLab);
  scat_elec_in_relab.rotateY(-fThetaLab);
    scat_elec_in_relab.rotateZ(-fPhiLab2);

  p_rest_relab.rotateZ(-fPhiLab);
  p_rest_relab.rotateY(-fThetaLab );
  p_rest_relab.rotateZ(-fPhiLab2);

//   cout << endl;
//   cout << " ===== Hell Boost =====" << endl;
//   cout << "          LAB" << endl;
//   cout << " beam prot: " << p_in << endl;
//   cout << " beam elec: " << e_in << endl;
//   cout << " scat elec: " << scat_elec << endl;
//   cout << " rest vect: " << p_rest << endl;
//   cout << " test vect: " << p_test << endl;
//   cout << "          CMS" << endl;
//   cout << " beta vect: " << bVec2 << endl; 
//   cout << " beam prot: " << p_in_hcms << endl;
//   cout << " beam elec: " << e_in_hcms << endl;
//   cout << " scat elec: " << scat_elec_in_hcms << endl;
//   cout << " rest vect: " << p_rest_hcms << endl;
//   cout << " test vect: " << p_test_hcms << endl;
//   cout << " xz   vect: " << xz_hcms << endl;
//   cout << "      LAB reboost" << endl;
//   cout << " beta vect: " << bVec << endl; 
//   cout << " beam prot: " << p_in_relab << endl;
//   cout << " beam elec: " << e_in_relab << endl;
//   cout << " scat elec: " << scat_elec_in_relab << endl;
//   cout << " rest vect: " << p_rest_relab << endl;
//   cout << " test vect: " << p_test_relab << endl;
//   cout << " xz   vect: " << xz_relab << endl;

  // ---- now let's boost back ---



  unsigned int njet = pjl.upper();
  
  for(unsigned int i = 1; i <= njet; i++) {
    pjl[i].boost(bVec);
    pjl[i].rotateZ(-fPhiLab);
    pjl[i].rotateY(-fThetaLab);
    pjl[i].rotateZ(-fPhiLab2);
  }


//   unsigned int njetf = fjets.upper();
//   for(unsigned int i = 1; i <= njetf; i++) {
//     fjetslab[i].boost(bVec);
//     fjetslab[i].rotateZ(-fPhiLab);
//     fjetslab[i].rotateY(-fThetaLab);
//     fjetslab[i].rotateZ(-fPhiLab2);
//   }




//   cout << " hell boost" << endl;
//   cout << p1[-1] << endl;
//   cout << p1[-2] << endl;
//   cout << p_rest << endl;
//   cout << endl;
//   cout << p2[-1] << endl;
//   cout << p2[-2] << endl;
//   cout << p_rest_hcms << endl;
//   cout << endl;

//   //  double x = p[-1].T()/p[hadron(0)].T();
// //   double x = p1[-1].T()/p_in.T();
// //   double bz = (1.0 - x)/(1.0 + x);

//   //   threevector<double> bVec = -((p_in).boostVector());
//   //   threevector<double> bVec = -((pbreit[-1] + p_in_hcms).boostVector());
//   //   threevector<double> bVec = -((pbreit[-1] + pbreit[hadron(0)]).boostVector());
//   //  threevector<double> bVec = -((pbreit[-1] + p_in).boostVector());
//   lorentzvector<double> p0(pbreit[hadron(0)]);

//   p_in_hell = p_rest_hcms;
//   p_in_hell.boost(bVec);

//   p0 = p_in_hcms;

//   //  cout << "***************" << endl;
//   // cout << " p[-1], p[-1].T()" << p[-1] << p[-1].T() << endl;
//   // cout << "bVec=" << bVec << endl;
//   // cout << "pbreit[-1]=" << pbreit[-1] << endl;
//   // cout << "pbreit[hadron(0)]" << pbreit[hadron(0)] << endl;  
//   // cout << "***************" << endl;

//   p0.boost(bVec);
//   double phi = p0.phi(), theta = p0.theta();
//   unsigned int njet = pjl.upper();
//   int  low = plab.lower(), up = plab.upper();
  
// //   cout << p2[0] << endl;

// //   cout << " low: " << low << " up: " << up << endl;
  
//   for(int i = low; i <= up; i++) {
//     p2[i].boost(bVec);
//     //    cout  << "in back boost " <<  p2[i] << endl;
// //     pjl[i].rotateZ(-phi);
// //     //cout  << "in back boost after phi rot " <<  pjl[i] << endl;
// //     ///pjl[i].rotateZ(phi);
// //     pjl[i].rotateY(-theta); 
// //     //cout  << "in back boost after theta rot " <<  pjl[i] << endl;
// //     // pjl[i].rotateY(theta); 
// //     pjl[i].boost(0.0, 0.0, bz);
//   }

//   //  cout << p2[0] << endl;
//   cout << p2[-1] << endl;
//   cout << p2[-2] << endl;
//   cout << p_in_hell << endl;
//   cout << endl;


}
