//
//   fastNLO v2 author code 
//   scenario: fnh501x
//   H1 Trijets (@high Q2)
//   HERA-II
//   for kT algorithm
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
//   inputfunc   (-> user edits)
//   psinput     (-> user edits)
//   initfunc    (don't touch)
//   xalfaem     (don't touch)
//   userfunc    (-> user edits)
//   end_of_event (don't touch)
//   phys_output (don't touch)
//   inittable   (-> user edits)
//   boost_back_to_lab (don't touch)
//
// Implementing a new scenario requires to edit:
//  - the jet algorithm ("#include" statement and assignement of "jetclus")
//  - number of jets (inputfunc)
//  - center-of-mass energy (psinput)
//  - compute observable, determine bin No. (userfunc)
//  - declare all variables for table, define bin boundaries (inittable)
//  
// ================================================================
// 
// last modifications
// 11/09/22 DB update to v2.1
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

#include <kT_clus.h>         // fastNLO user: .h file for jet algorithm

#include "pdf-cteq6.h"
#include "pdf-dis-dummy.h"
#include "fnloTable.h"
#include "fnloBlockBNlojet.h"


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

   // jet algorithm
   kT_clus_long jetclus;   // fastNLO user: define jet algorithm (consistent with .h file above)
   
   // the jet structore in breit & lab. frame
   bounded_vector<lorentzvector<double> > pjb, pjl; 
   bounded_vector<unsigned int> jet;
   
   //  event in breit frame
   event_dis pbreit;
  
  // boost the jet momenta back to the laboratory frame
  void boost_back_to_lab(const event_dis&);
  
   // --- fastNLO definitions (not for user)
   fnloTable *tabledXidQ2;
   fnloTable *tabledXi;
   fnloTable *tabledPtdQ2;
   fnloTable *tabledPt;
   fnloTable *tabledQ2;
   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table
   time_t start_time;
   
   void inittabledXidQ2(string tablefilename);
   void inittabledXi(string tablefilename);
   void inittabledPtdQ2(string tablefilename);
   void inittabledPt(string tablefilename);
   void inittabledQ2(string tablefilename);

   bool nlo;
   bool bWarmup;
   int  fNWarmupPrint;

};
  
user_base_dis * userfunc() {
  return new UserDIS;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
  // --- fastNLO user: select the number of jets for your observable
  nj = 3U;
  
  // --- number of the up and down type flavours
  nu = 2U;
  nd = 3U;
} 

void psinput(phasespace_dis *ps, double& el, double& eh, double& q2min, double& q2max, 
			 double& xmin, double& xmax, double& ymin, double& ymax)
{
  // --- fastNLO user: set the energies of the incoming lepton and hadron in the laboratory frame
  el = 27.6;      // in GeV
  eh = 920.0;     // in GeV
  
  // --- fastNLO user: define the Q^2 boundaries of the phase space
  q2min = 150.0;    // in GeV^2
  q2max = 15000.0;   // in GeV^2 

  // --- fastNLO user: define the x_Bjorken boundaries of the phase space
  xmin = 0.0;
  xmax = 1.0;
  
  // --- fastNLO user: define the y boundaries of the phase space
  ymin = 0.2;
  ymax = 0.7;

  //   You can use your own phase generator. 
  //   Here we use the default.
  ps = 0;
} 


void UserDIS::initfunc(unsigned int)
{
   // ---- Initialize event counters
   nevents = 0;
   // Set some defaults
   if (nwrite==0) nwrite = 10000000; // 5000000

   start_time = std::time(0);
   
}

extern"C" double xalfaem_(double *);

double xalpha_em(double mq2) {
  return xalfaem_(&mq2);
}



void UserDIS::userfunc(const event_dis& p, const amplitude_dis& amp)
{
   
  double alem = 0; // for running alpha EM
  // --- x where PDF is probed ---
  double x = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);
  double Q2 = -((p[-1] - p[-2]).mag2());
  double Q = sqrt(Q2); 
  
  alem = xalpha_em(Q2);

   // --- fastNLO user:
   //     Here is your playground where you compute your observable 
   //     and the bin number ("obsbin") which gets passed to
   //     fastNLO's table filling code.
   //     (all pT, ET and E are in GeV)

   // --- declare and initialize phase space cut variables
   //     here: H1 cuts
   double pt = 0.0;
   double pTmin2 = 5.0, etamin = -1., etamax = 2.5;
   double pTmax = 50.0;
   double Mn_Cut = 16.0;

   double pt_mean = 0.0, pt_sum = 0.0;
   unsigned int jetnum = 0;
   lorentzvector<double> psum2;
   lorentzvector<double> psum3;


   // --- copy the momenta and boost to the breit frame ----
   pbreit = p; lab_to_breit(pbreit);

   // --- do the cluster analysis ----
   jetclus.set_up(1,false,2);
   jetclus(pbreit); jetclus.incl(pjb, jet);
   unsigned int nj = pjb.upper(); 


   // --- jet structure in laboratory frame -----
   pjl = pjb; boost_back_to_lab(p);

   
   unsigned int IndexArray[nj+1];
   double PtArray[nj+1];
   double Theta_pj[nj+1];

   for(unsigned int i = 1; i <= nj; i++) {
     pt = pjb[i].perp();
     double etalab = pjl[i].prapidity();
     if(pt > pTmin2 && pt < pTmax && etamin < etalab && etalab < etamax) {
       jetnum++;
       PtArray[jetnum] = pt;
       Theta_pj[jetnum] = pjb[i].theta();
       IndexArray[jetnum] = i;
    }
   }

  const unsigned int nSim = 3;
  if ( jetnum < nSim) return;
   
  double t = 0;
  unsigned int y = 0;
  int fred = 0;
 
  for(unsigned int i = 1; i <= jetnum; i++) {
    for (y=1; y < (jetnum+1-i); y++){
      if (PtArray[y] < PtArray[y+1]){        
	t=PtArray[y];
	fred = IndexArray[y];
	PtArray[y]=PtArray[y+1];
	IndexArray[y]=IndexArray[y+1];
	PtArray[y+1]=t;
	IndexArray[y+1]=fred;
	t = Theta_pj[y];
	Theta_pj[y] = Theta_pj[y+1];
	Theta_pj[y+1] = t;
      }
    }          
  }
  

  pt_mean = ( PtArray[1] + PtArray[2] + PtArray[3] ) / 3;
  psum2 = pjb[IndexArray[1]] + pjb[IndexArray[2]];
  psum3 = pjb[IndexArray[1]] + pjb[IndexArray[2]] + pjb[IndexArray[3]];

  //----- jet mass cuts -----
  double Mn2 = psum2.mag();
  double Mn3 = psum3.mag();

  if ( Mn2 < Mn_Cut) return;

  double xbj = Q/(2.0*pbreit[hadron(0)].T());
  double xi = xbj*(1+Mn3*Mn3/Q2);


  //  ------------------------
  //    set the scale here !
  double mu = sqrt((pt_mean*pt_mean+Q2)/2); 


  // ----- filling ----- //
  double prefactor = alem*alem;  // 1./A2->BinSize[bin] *  is now in FillEvent
  int bin = -1;

  // ---- dQ2dXi ---- //
  bin = tabledXidQ2->GetBinNumber(xi,Q2);
  tabledXidQ2->FillEventDIS(bin,x,sqrt(Q2),pt_mean,mu,amp,dummypdf,pdf,prefactor);

  // ---- dXi ---- //
  bin = tabledXi->GetBinNumber(xi);
  tabledXi->FillEventDIS(bin,x,sqrt(Q2),pt_mean,mu,amp,dummypdf,pdf,prefactor);

  // ---- dPtdQ2 ---- //
  bin = tabledPtdQ2->GetBinNumber(pt_mean,Q2);
  tabledPtdQ2->FillEventDIS(bin,x,sqrt(Q2),pt_mean,mu,amp,dummypdf,pdf,prefactor);

  // ---- dPt ---- //
  bin = tabledPt->GetBinNumber(pt_mean);
  tabledPt->FillEventDIS(bin,x,sqrt(Q2),pt_mean,mu,amp,dummypdf,pdf,prefactor);

  // ---- dPtdQ2 ---- //
  bin = tabledQ2->GetBinNumber(Q2);
  tabledQ2->FillEventDIS(bin,x,sqrt(Q2),pt_mean,mu,amp,dummypdf,pdf,prefactor);


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
      tabledXidQ2->WriteTable(nevents);
      tabledXi->WriteTable(nevents);
      tabledPtdQ2->WriteTable(nevents);
      tabledPt->WriteTable(nevents);
      tabledQ2->WriteTable(nevents);
      printf("done.\n");
   }
}

void UserDIS::phys_output(const std::basic_string<char>& __file_name, 
                          unsigned long __save, bool __txt) 
{
 
   // --- determine whether we are running LO or NLO
   const char* const file = __file_name.c_str(); 
   if(strstr(file,"born")!=NULL){
      nlo = false;
   }else{
      if(strstr(file,"nlo")!=NULL){
         nlo = true;
      }else{
         printf("This module can only be run at Born level or at NLO.\n");
         exit(1);
      }
   }
   nwrite = __save;


   bWarmup		= false;
   fNWarmupPrint	= 5000000;
   

  string tablefilenamedXidQ2;
  tablefilenamedXidQ2.assign(__file_name.c_str());
  tablefilenamedXidQ2 += "dXidQ2.tab";
  inittabledXidQ2(tablefilenamedXidQ2);

  string tablefilenamedXi;
  tablefilenamedXi.assign(__file_name.c_str());
  tablefilenamedXi += "dXi.tab";
  inittabledXi(tablefilenamedXi);

  string tablefilenamedPtdQ2;
  tablefilenamedPtdQ2.assign(__file_name.c_str());
  tablefilenamedPtdQ2 += "dPtdQ2.tab";
  inittabledPtdQ2(tablefilenamedPtdQ2);

  string tablefilenamedPt;
  tablefilenamedPt.assign(__file_name.c_str());
  tablefilenamedPt += "dPt.tab";
  inittabledPt(tablefilenamedPt);

  string tablefilenamedQ2;
  tablefilenamedQ2.assign(__file_name.c_str());
  tablefilenamedQ2 += "dQ2.tab";
  inittabledQ2(tablefilenamedQ2);




}

void UserDIS::inittabledXidQ2(string tablefilename){

   // --- set up fastNLO
   tabledXidQ2 = new fnloTable(tablefilename);

   // --- fastNLO: fill variable for table header block A1
   tabledXidQ2->GetBlockA1()->SetScenName("H1-Trijets-dXidQ2-HighQ2");  // NO SPACES HERE !!! // - fastNLO user: set scenario name
   tabledXidQ2->GetBlockA1()->SetNcontrib(1);
   tabledXidQ2->GetBlockA1()->SetNmult(0);
   tabledXidQ2->GetBlockA1()->SetNdata(0);
   tabledXidQ2->GetBlockA2()->SetIpublunits(12);  // - fastNLO user: set cross section units
                                            //                 (negative power of ten)
                                            //                 fb:15  pb:12  nb:9  etc.
   // --- fastNLO: fill variables for table header block A2
   fnloBlockA2 *A2 =  tabledXidQ2->GetBlockA2();

   // --- fastNLO user: up to 30 strings to describe the scenario
   A2->ScDescript.push_back("Trijets - d2sigma/dXidQ2 (pb/GeV2)");
   A2->ScDescript.push_back("H1 Collaboration");
   A2->ScDescript.push_back("Trijet cross section in DIS");
   A2->ScDescript.push_back("kT-algorithm, R=1, HighQ2, HERA-II");
   A2->ScDescript.push_back("H1prelim-11-032, Fig. 13");

   A2->Ecms = sqrt(4.*920.*27.6);	// --- fastNLO user: set sqrt(s)
   A2->ILOord = 2;			// --- fastNLO user: power of LO contribution for process // 1 for dijets, 2 for trijets (nj=3)
   A2->NDim = 2;			// --- fastNLO user: No of dimensions in which observable is binned
   A2->DimLabel.push_back("xi");	// --- fastNLO user: label of 1st dimension
   A2->IDiffBin.push_back(2);
   A2->DimLabel.push_back("Q^2");	// --- fastNLO user: label of 2nd dimension
   A2->IDiffBin.push_back(2);


   // --- fastNLO user: bin definitions - here in Q2 and ET
   const int nq2bins = 6;
   double q2bins[nq2bins+1] = {150.,200.,270.,400.,700.,5000.,15000.};

   const int netbins[nq2bins] = {3, 3, 3, 3, 2, 1};
   double etbins[nq2bins][4] = {
     { 0.01 , 0.04 , 0.08 , 0.50 },
     { 0.01 , 0.04 , 0.08 , 0.50 },
     { 0.01 , 0.04 , 0.08 , 0.50 },
     { 0.01 , 0.04 , 0.08 , 0.50 },
     { 0.04 , 0.08 , 0.50 },
     { 0.08 , 0.50 }
   };


   // ---- initalize the bingrids and the normalizations ---- //
   vector<double*> vetbins(nq2bins);
   for (int i=0;i<vetbins.size();i++) vetbins[i]=etbins[i]; 
   A2->InitBinning( nq2bins , q2bins , netbins , vetbins );

   // --- fastNLO table block B
   fnloBlockBNlojet *B = new fnloBlockBNlojet(tabledXidQ2->GetBlockA1(),tabledXidQ2->GetBlockA2());
   tabledXidQ2->CreateBlockB(0,B);

   // ---- initalize all constants to make fastNLO behave like for DIS scenario ---- //
   B->InitDISConstants(A2,nlo);


   // --- arrays for extreme x and (default) scale values (computed in Warm-Up run) ---- //
   double xlim[A2->NObsBin];
   double scale1lo[A2->NObsBin];
   double scale1hi[A2->NObsBin];
   double scale2lo[A2->NObsBin];
   double scale2hi[A2->NObsBin];
   for(int i=0;i<A2->NObsBin;i++){
      xlim[i] = 1.1e-07;
      scale1lo[i] = 1.0, scale1hi[i]=9.9e10;
      scale2lo[i] = 1.0, scale2hi[i]=9.9e10;
   }


   // --- fastNLO user: before running a new scenario for the first time,
   //     the following block (between "start" and "end") should be 
   //     completely removed. The first run must be a "Warm-Up Run" (IWarmUp=1)
   //     which produces an initialization block as output. This should be copied
   //     and pasted below. These initialization values must be used for all
   //     production jobs (IWarmUp=0) for a given scenario (otherwise the result 
   //     tables can not be merged). 
   // --- initialize variables for WarmUp run
   B->IWarmUp = bWarmup;
   if (bWarmup ) cout<<"\nThis is a warmup run.\n"<<endl;      // --- fastNLO user: do the Warm-Up run
   B->IWarmUpPrint = fNWarmupPrint; //3 000 0000
   //
   // --------- fastNLO: Warm-Up run results (start)
   // -749967296 contributions (!= events) in warm-up run
   xlim [0] = 1.000000e-02 , scale1lo [0] =   12.2474 , scale1hi [0] =   14.1421 , scale2lo [0] =    5.0194 , scale2hi [0] =   30.0915;
   xlim [1] = 4.000000e-02 , scale1lo [1] =   12.2474 , scale1hi [1] =   14.1421 , scale2lo [1] =    5.0212 , scale2hi [1] =   38.1906;
   xlim [2] = 8.000001e-02 , scale1lo [2] =   12.2474 , scale1hi [2] =   14.1421 , scale2lo [2] =    5.6748 , scale2hi [2] =   49.8321;
   xlim [3] = 1.000240e-02 , scale1lo [3] =   14.1421 , scale1hi [3] =   16.4317 , scale2lo [3] =    5.0215 , scale2hi [3] =   29.3908;
   xlim [4] = 4.000000e-02 , scale1lo [4] =   14.1421 , scale1hi [4] =   16.4317 , scale2lo [4] =    5.0226 , scale2hi [4] =   37.7886;
   xlim [5] = 8.000000e-02 , scale1lo [5] =   14.1421 , scale1hi [5] =   16.4317 , scale2lo [5] =    5.4767 , scale2hi [5] =   49.7473;
   xlim [6] = 1.071501e-02 , scale1lo [6] =   16.4317 , scale1hi [6] =   20.0000 , scale2lo [6] =    5.0217 , scale2hi [6] =   28.7211;
   xlim [7] = 4.000000e-02 , scale1lo [7] =   16.4317 , scale1hi [7] =   20.0000 , scale2lo [7] =    5.0241 , scale2hi [7] =   37.9151;
   xlim [8] = 8.000000e-02 , scale1lo [8] =   16.4317 , scale1hi [8] =   20.0000 , scale2lo [8] =    5.1759 , scale2hi [8] =   49.9151;
   xlim [9] = 1.261546e-02 , scale1lo [9] =   20.0000 , scale1hi [9] =   26.4575 , scale2lo [9] =    5.0132 , scale2hi [9] =   27.6785;
   xlim [10] = 4.000000e-02 , scale1lo [10] =   20.0000 , scale1hi [10] =   26.4575 , scale2lo [10] =    5.0341 , scale2hi [10] =   38.4812;
   xlim [11] = 8.000000e-02 , scale1lo [11] =   20.0000 , scale1hi [11] =   26.4575 , scale2lo [11] =    5.1155 , scale2hi [11] =   49.8424;
   xlim [12] = 4.000000e-02 , scale1lo [12] =   26.4575 , scale1hi [12] =   70.7070 , scale2lo [12] =    5.0165 , scale2hi [12] =   38.3599;
   xlim [13] = 8.000000e-02 , scale1lo [13] =   26.4575 , scale1hi [13] =   70.7107 , scale2lo [13] =    5.0140 , scale2hi [13] =   49.9122;
   xlim [14] = 8.002628e-02 , scale1lo [14] =   70.7107 , scale1hi [14] =  122.4745 , scale2lo [14] =    5.0077 , scale2hi [14] =   49.9580;

   //    // 17000000 contributions (!= events) in warm-up run
   //    xlim [0] = 1.025726e-02 , scale1lo [0] =   12.2476 , scale1hi [0] =   14.1420 , scale2lo [0] =    5.3104 , scale2hi [0] =   24.2194;
   //    xlim [1] = 4.000406e-02 , scale1lo [1] =   12.2476 , scale1hi [1] =   14.1421 , scale2lo [1] =    5.3441 , scale2hi [1] =   31.6729;
   //    xlim [2] = 8.000420e-02 , scale1lo [2] =   12.2476 , scale1hi [2] =   14.1421 , scale2lo [2] =    6.5696 , scale2hi [2] =   48.4771;
   //    xlim [3] = 1.178356e-02 , scale1lo [3] =   14.1422 , scale1hi [3] =   16.4306 , scale2lo [3] =    5.1763 , scale2hi [3] =   22.9332;
   //    xlim [4] = 4.000667e-02 , scale1lo [4] =   14.1424 , scale1hi [4] =   16.4316 , scale2lo [4] =    5.0278 , scale2hi [4] =   30.1709;
   //    xlim [5] = 8.000019e-02 , scale1lo [5] =   14.1422 , scale1hi [5] =   16.4316 , scale2lo [5] =    6.8610 , scale2hi [5] =   47.1641;
   //    xlim [6] = 1.233040e-02 , scale1lo [6] =   16.4321 , scale1hi [6] =   19.9999 , scale2lo [6] =    5.1331 , scale2hi [6] =   24.6140;
   //    xlim [7] = 4.000115e-02 , scale1lo [7] =   16.4318 , scale1hi [7] =   19.9997 , scale2lo [7] =    5.3304 , scale2hi [7] =   31.8115;
   //    xlim [8] = 8.000148e-02 , scale1lo [8] =   16.4318 , scale1hi [8] =   19.9999 , scale2lo [8] =    6.6321 , scale2hi [8] =   48.8721;
   //    xlim [9] = 1.353028e-02 , scale1lo [9] =   20.0003 , scale1hi [9] =   26.4572 , scale2lo [9] =    5.2592 , scale2hi [9] =   24.3039;
   //    xlim [10] = 4.000513e-02 , scale1lo [10] =   20.0002 , scale1hi [10] =   26.4566 , scale2lo [10] =    5.2742 , scale2hi [10] =   31.0560;
   //    xlim [11] = 8.000024e-02 , scale1lo [11] =   20.0001 , scale1hi [11] =   26.4574 , scale2lo [11] =    5.6309 , scale2hi [11] =   48.6177;
   //    xlim [12] = 4.000021e-02 , scale1lo [12] =   26.4578 , scale1hi [12] =   70.1974 , scale2lo [12] =    5.0329 , scale2hi [12] =   32.1390;
   //    xlim [13] = 8.000233e-02 , scale1lo [13] =   26.4575 , scale1hi [13] =   70.7106 , scale2lo [13] =    5.1892 , scale2hi [13] =   48.8396;
   //    xlim [14] = 8.356003e-02 , scale1lo [14] =   70.7109 , scale1hi [14] =  122.4739 , scale2lo [14] =    5.2018 , scale2hi [14] =   48.3149;
   // xlim[0]=0.01002, mulo[0]=   9.430, muup[0]=  20.535;
   // xlim[1]=0.04000, mulo[1]=   9.535, muup[1]=  26.992;
   // xlim[2]=0.08000, mulo[2]=  10.184, muup[2]=  35.975;
   // xlim[3]=0.01066, mulo[3]=  10.647, muup[3]=  21.103;
   // xlim[4]=0.04000, mulo[4]=  10.771, muup[4]=  27.091;
   // xlim[5]=0.08000, mulo[5]=  11.235, muup[5]=  36.290;
   // xlim[6]=0.01163, mulo[6]=  12.310, muup[6]=  22.084;
   // xlim[7]=0.04000, mulo[7]=  12.315, muup[7]=  27.389;
   // xlim[8]=0.08000, mulo[8]=  12.609, muup[8]=  37.551;
   // xlim[9]=0.01327, mulo[9]=  14.630, muup[9]=  24.855;
   // xlim[10]=0.04000, mulo[10]=  14.737, muup[10]=  29.148;
   // xlim[11]=0.08000, mulo[11]=  14.931, muup[11]=  38.908;
   // xlim[12]=0.04000, mulo[12]=  19.095, muup[12]=  49.803;
   // xlim[13]=0.08000, mulo[13]=  19.268, muup[13]=  60.329;
   // xlim[14]=0.08177, mulo[14]=  50.187, muup[14]=  91.862;

   // --------- fastNLO: Warm-Up run results (end)


   // ---- set number-of x-nodes ---- //
   B->SetNumberOfXNodesPerMagnitude(20,xlim);


   // ---- scale description ---- //
   B->SetScale1Name( "Q" );
   B->SetScale2Name( "<p_T>" );

   
   // ---- number of scale nodes for mu ---- //
   B->SetNumberOfScaleNodesScale1( 15 );
   B->SetNumberOfScaleNodesScale2( 12 );


   // ----- v2.0 tables ---- //
   if ( B->NScaleDep != 3 ){
      // ---- number of scale nodes ---- //
      B->SetNumberOfScaleNodes_v20(13); //! I used 20 previously 13// number of scale nodes for mu
      // ---- set scale variations  ---- //
      B->ScaleFac[0].push_back(1.0);
      // 	B->ScaleFac[0].push_back(0.5); // --- fastNLO user: add any number of
      // 	B->ScaleFac[0].push_back(2.0); // additional scale variations as desired
   }   



   // ---- final initializations ---- //
   B->InitFinalDISValues( A2 , xlim , scale1lo , scale1hi, scale2lo , scale2hi ); 


   // --- fastNLO user: decide whether to include a reference table (for 
   //                   precision studies, not for production jobs)
   //const bool doReference = true;
   const bool doReference = false;

}


void UserDIS::inittabledXi(string tablefilename){

   // --- set up fastNLO
   tabledXi = new fnloTable(tablefilename);

   // --- fastNLO: fill variable for table header block A1
   tabledXi->GetBlockA1()->SetScenName("H1-Trijets-dXi-HighQ2");  // NO SPACES HERE !!! // - fastNLO user: set scenario name
   tabledXi->GetBlockA1()->SetNcontrib(1);
   tabledXi->GetBlockA1()->SetNmult(0);
   tabledXi->GetBlockA1()->SetNdata(0);
   tabledXi->GetBlockA2()->SetIpublunits(12);  // - fastNLO user: set cross section units
                                            //                 (negative power of ten)
                                            //                 fb:15  pb:12  nb:9  etc.
   // --- fastNLO: fill variables for table header block A2
   fnloBlockA2 *A2 =  tabledXi->GetBlockA2();

   // --- fastNLO user: up to 30 strings to describe the scenario
   A2->ScDescript.push_back("Trijets - dsigma/dXi (pb/1)");
   A2->ScDescript.push_back("H1 Collaboration");
   A2->ScDescript.push_back("Trijet cross section in DIS");
   A2->ScDescript.push_back("kT-algorithm, R=1, HighQ2, HERA-II");
   A2->ScDescript.push_back("H1prelim-11-032, Fig. 11c");

   A2->Ecms = sqrt(4.*920.*27.6);	// --- fastNLO user: set sqrt(s)
   A2->ILOord = 2;			// --- fastNLO user: power of LO contribution for process // 1 for dijets, 2 for trijets (nj=3)
   //A2->NDim = 2;			// --- fastNLO user: No of dimensions in which observable is binned
   A2->NDim = 1;			// --- fastNLO user: No of dimensions in which observable is binned
   A2->DimLabel.push_back("xi");	// --- fastNLO user: label of 1st dimension
   A2->IDiffBin.push_back(2);
   // A2->DimLabel.push_back("Q^2");	// --- fastNLO user: label of 2nd dimension
   // A2->IDiffBin.push_back(1);		// I do not want to divide by the Q2-bin width


   // --- fastNLO user: bin definitions - here in Q2 and ET
//    const int nq2bins = 1;
//    double q2bins[nq2bins+1] = { 150.,15000.};

//    const int netbins[nq2bins] = {3};
//    double etbins[nq2bins][4] = {
//      { 0.01 , 0.04 , 0.08 , 0.50 }
//    };

    const int nxibins = 3;
    double xibins[4] = {
        0.01 , 0.04 , 0.08 , 0.50 
    };

   // ---- initalize the bingrids and the normalizations ---- //
   A2->InitBinning( nxibins , xibins );


   // --- fastNLO table block B
   fnloBlockBNlojet *B = new fnloBlockBNlojet(tabledXi->GetBlockA1(),tabledXi->GetBlockA2());
   tabledXi->CreateBlockB(0,B);

   // ---- initalize all constants to make fastNLO behave like for DIS scenario ---- //
   B->InitDISConstants(A2,nlo);



   // --- arrays for extreme x and (default) scale values (computed in Warm-Up run) ---- //
   double xlim[A2->NObsBin];
   double scale1lo[A2->NObsBin];
   double scale1hi[A2->NObsBin];
   double scale2lo[A2->NObsBin];
   double scale2hi[A2->NObsBin];
   for(int i=0;i<A2->NObsBin;i++){
      xlim[i] = 1.1e-07;
      scale1lo[i] = 1.0, scale1hi[i]=9.9e10;
      scale2lo[i] = 1.0, scale2hi[i]=9.9e10;
   }


   // --- initialize variables for WarmUp run
   B->IWarmUp = bWarmup;
   if (bWarmup ) cout<<"\nThis is a warmup run.\n"<<endl;      // --- fastNLO user: do the Warm-Up run
   B->IWarmUpPrint = fNWarmupPrint; //3 000 0000
   // --- fastNLO user: before running a new scenario for the first time,
   //     the following block (between "start" and "end") should be 
   //     completely removed. The first run must be a "Warm-Up Run" (IWarmUp=1)
   //     which produces an initialization block as output. This should be copied
   //     and pasted below. These initialization values must be used for all
   //     production jobs (IWarmUp=0) for a given scenario (otherwise the result 
   //     tables can not be merged). 
   //
   // --------- fastNLO: Warm-Up run results (start)
   xlim [0] = 1.000000e-02 , scale1lo [0] =   12.2474 , scale1hi [0] =   48.6219 , scale2lo [0] =    5.0132 , scale2hi [0] =   30.0915;
   xlim [1] = 4.000000e-02 , scale1lo [1] =   12.2474 , scale1hi [1] =   71.7872 , scale2lo [1] =    5.0165 , scale2hi [1] =   38.4812;
   xlim [2] = 8.000000e-02 , scale1lo [2] =   12.2474 , scale1hi [2] =  122.4745 , scale2lo [2] =    5.0077 , scale2hi [2] =   49.9580;
//    xlim [0] = 1.025726e-02 , scale1lo [0] =   12.2476 , scale1hi [0] =   46.1598 , scale2lo [0] =    5.0872 , scale2hi [0] =   24.6140;
//    xlim [1] = 4.000021e-02 , scale1lo [1] =   12.2476 , scale1hi [1] =   70.1974 , scale2lo [1] =    5.0278 , scale2hi [1] =   32.1390;
//    xlim [2] = 8.000019e-02 , scale1lo [2] =   12.2476 , scale1hi [2] =  122.4739 , scale2lo [2] =    5.1892 , scale2hi [2] =   48.8721;
   // --------- fastNLO: Warm-Up run results (end)

   // ---- set number-of x-nodes ---- //
   B->SetNumberOfXNodesPerMagnitude(20,xlim);

   // ---- scale description ---- //
   B->SetScale1Name( "Q" );
   B->SetScale2Name( "<p_T>" );

   // ---- number of scale nodes for mu ---- //
   B->SetNumberOfScaleNodesScale1( 15 );
   B->SetNumberOfScaleNodesScale2( 12 );

   // ---- final initializations ---- //
   B->InitFinalDISValues( A2 , xlim , scale1lo , scale1hi, scale2lo , scale2hi ); 

   // --- fastNLO user: decide whether to include a reference table (for 
   const bool doReference = false;

}


void UserDIS::inittabledPtdQ2(string tablefilename){

   // --- set up fastNLO
   tabledPtdQ2 = new fnloTable(tablefilename);

   // --- fastNLO: fill variable for table header block A1
   tabledPtdQ2->GetBlockA1()->SetScenName("H1-Trijets-dPtdQ2-HighQ2");  // NO SPACES HERE !!! // - fastNLO user: set scenario name
   tabledPtdQ2->GetBlockA1()->SetNcontrib(1);
   tabledPtdQ2->GetBlockA1()->SetNmult(0);
   tabledPtdQ2->GetBlockA1()->SetNdata(0);
   tabledPtdQ2->GetBlockA2()->SetIpublunits(12);  // - fastNLO user: set cross section units
                                            //                 (negative power of ten)
                                            //                 fb:15  pb:12  nb:9  etc.
   // --- fastNLO: fill variables for table header block A2
   fnloBlockA2 *A2 =  tabledPtdQ2->GetBlockA2();

   // --- fastNLO user: up to 30 strings to describe the scenario
   A2->ScDescript.push_back("Trijets - d2sigma/d<Pt>dQ2 (pb/GeV3)");
   A2->ScDescript.push_back("H1 Collaboration");
   A2->ScDescript.push_back("Trijet cross section in DIS");
   A2->ScDescript.push_back("kT-algorithm, R=1, HighQ2, HERA-II");
   A2->ScDescript.push_back("H1prelim-11-032, Fig. 12");

   A2->Ecms = sqrt(4.*920.*27.6);	// --- fastNLO user: set sqrt(s)
   A2->ILOord = 2;			// --- fastNLO user: power of LO contribution for process // 1 for dijets, 2 for trijets (nj=3)
   A2->NDim = 2;			// --- fastNLO user: No of dimensions in which observable is binned
   A2->DimLabel.push_back("<p_T>");	// --- fastNLO user: label of 1st dimension
   A2->IDiffBin.push_back(2);
   A2->DimLabel.push_back("Q^2");	// --- fastNLO user: label of 2nd dimension
   A2->IDiffBin.push_back(2);


   // --- fastNLO user: bin definitions - here in Q2 and ET
   const int nq2bins = 6;
   double q2bins[nq2bins+1] = {150.,200.,270.,400.,700.,5000.,15000.};
 
  const int netbins[nq2bins] = {3, 3, 3, 3, 3, 3};
   double etbins[nq2bins][4] = {
      {7.,11.,18.,30.},
      {7.,11.,18.,30.},
      {7.,11.,18.,30.},
      {7.,11.,18.,30.},
      {7.,11.,18.,30.},
      {7.,11.,18.,30.}
   };

   // ---- initalize the bingrids and the normalizations ---- //
   vector<double*> vetbins(nq2bins);
   for (int i=0;i<vetbins.size();i++) vetbins[i]=etbins[i]; 
   A2->InitBinning( nq2bins , q2bins , netbins , vetbins );

   // --- fastNLO table block B
   fnloBlockBNlojet *B = new fnloBlockBNlojet(tabledPtdQ2->GetBlockA1(),tabledPtdQ2->GetBlockA2());
   tabledPtdQ2->CreateBlockB(0,B);

   // ---- initalize all constants to make fastNLO behave like for DIS scenario ---- //
   B->InitDISConstants(A2,nlo);


   // --- arrays for extreme x and (default) scale values (computed in Warm-Up run) ---- //
   double xlim[A2->NObsBin];
   double scale1lo[A2->NObsBin];
   double scale1hi[A2->NObsBin];
   double scale2lo[A2->NObsBin];
   double scale2hi[A2->NObsBin];
   for(int i=0;i<A2->NObsBin;i++){
      xlim[i] = 1.1e-07;
      scale1lo[i] = 1.0, scale1hi[i]=9.9e10;
      scale2lo[i] = 1.0, scale2hi[i]=9.9e10;
   }


   // --- fastNLO user: before running a new scenario for the first time,
   //     the following block (between "start" and "end") should be 
   //     completely removed. The first run must be a "Warm-Up Run" (IWarmUp=1)
   //     which produces an initialization block as output. This should be copied
   //     and pasted below. These initialization values must be used for all
   //     production jobs (IWarmUp=0) for a given scenario (otherwise the result 
   //     tables can not be merged). 
   // --- initialize variables for WarmUp run
   B->IWarmUp = bWarmup;
   if (bWarmup ) cout<<"\nThis is a warmup run.\n"<<endl;      // --- fastNLO user: do the Warm-Up run
   B->IWarmUpPrint = fNWarmupPrint; //3 000 0000
   //
   // --------- fastNLO: Warm-Up run results (start)
   // -1014967296 contributions (!= events) in warm-up run
   xlim [0] = 8.959844e-03 , scale1lo [0] =   12.2474 , scale1hi [0] =   14.1421 , scale2lo [0] =    7.0000 , scale2hi [0] =   11.0000;
   xlim [1] = 1.756865e-02 , scale1lo [1] =   12.2474 , scale1hi [1] =   14.1421 , scale2lo [1] =   11.0000 , scale2hi [1] =   18.0000;
   xlim [2] = 4.370402e-02 , scale1lo [2] =   12.2474 , scale1hi [2] =   14.1421 , scale2lo [2] =   18.0000 , scale2hi [2] =   30.0000;
   xlim [3] = 9.720735e-03 , scale1lo [3] =   14.1421 , scale1hi [3] =   16.4317 , scale2lo [3] =    7.0000 , scale2hi [3] =   11.0000;
   xlim [4] = 1.830315e-02 , scale1lo [4] =   14.1421 , scale1hi [4] =   16.4317 , scale2lo [4] =   11.0000 , scale2hi [4] =   18.0000;
   xlim [5] = 4.419733e-02 , scale1lo [5] =   14.1421 , scale1hi [5] =   16.4317 , scale2lo [5] =   18.0000 , scale2hi [5] =   30.0000;
   xlim [6] = 1.071671e-02 , scale1lo [6] =   16.4317 , scale1hi [6] =   20.0000 , scale2lo [6] =    7.0000 , scale2hi [6] =   11.0000;
   xlim [7] = 1.936737e-02 , scale1lo [7] =   16.4317 , scale1hi [7] =   20.0000 , scale2lo [7] =   11.0000 , scale2hi [7] =   18.0000;
   xlim [8] = 4.528831e-02 , scale1lo [8] =   16.4317 , scale1hi [8] =   20.0000 , scale2lo [8] =   18.0000 , scale2hi [8] =   30.0000;
   xlim [9] = 1.261546e-02 , scale1lo [9] =   20.0000 , scale1hi [9] =   26.4575 , scale2lo [9] =    7.0000 , scale2hi [9] =   11.0000;
   xlim [10] = 2.118810e-02 , scale1lo [10] =   20.0000 , scale1hi [10] =   26.4575 , scale2lo [10] =   11.0000 , scale2hi [10] =   18.0000;
   xlim [11] = 4.693139e-02 , scale1lo [11] =   20.0000 , scale1hi [11] =   26.4575 , scale2lo [11] =   18.0000 , scale2hi [11] =   30.0000;
   xlim [12] = 1.705440e-02 , scale1lo [12] =   26.4575 , scale1hi [12] =   70.7107 , scale2lo [12] =    7.0000 , scale2hi [12] =   11.0000;
   xlim [13] = 2.543978e-02 , scale1lo [13] =   26.4575 , scale1hi [13] =   70.7107 , scale2lo [13] =   11.0000 , scale2hi [13] =   18.0000;
   xlim [14] = 5.158560e-02 , scale1lo [14] =   26.4575 , scale1hi [14] =   70.7107 , scale2lo [14] =   18.0000 , scale2hi [14] =   30.0000;
   xlim [15] = 7.858367e-02 , scale1lo [15] =   70.7107 , scale1hi [15] =  122.4745 , scale2lo [15] =    7.0000 , scale2hi [15] =   11.0000;
   xlim [16] = 8.663615e-02 , scale1lo [16] =   70.7107 , scale1hi [16] =  122.4745 , scale2lo [16] =   11.0000 , scale2hi [16] =   18.0000;
   xlim [17] = 1.127987e-01 , scale1lo [17] =   70.7107 , scale1hi [17] =  122.4745 , scale2lo [17] =   18.0000 , scale2hi [17] =   30.0000;

   // --------- fastNLO: Warm-Up run results (end)


   // ---- set number-of x-nodes ---- //
   B->SetNumberOfXNodesPerMagnitude(10,xlim);


   // ---- scale description ---- //
   B->SetScale1Name( "Q" );
   B->SetScale2Name( "<p_T>" );

   
   // ---- number of scale nodes for mu ---- //
   B->SetNumberOfScaleNodesScale1( 15 );
   B->SetNumberOfScaleNodesScale2( 12 );

   // ---- final initializations ---- //
   B->InitFinalDISValues( A2 , xlim , scale1lo , scale1hi, scale2lo , scale2hi ); 

}

void UserDIS::inittabledQ2(string tablefilename){

   // --- set up fastNLO
   tabledQ2 = new fnloTable(tablefilename);

   // --- fastNLO: fill variable for table header block A1
   tabledQ2->GetBlockA1()->SetScenName("H1-Trijets-dQ2-HighQ2");  // NO SPACES HERE !!! // - fastNLO user: set scenario name
   tabledQ2->GetBlockA1()->SetNcontrib(1);
   tabledQ2->GetBlockA1()->SetNmult(0);
   tabledQ2->GetBlockA1()->SetNdata(0);
   tabledQ2->GetBlockA2()->SetIpublunits(12);  // - fastNLO user: set cross section units
                                            //                 (negative power of ten)
                                            //                 fb:15  pb:12  nb:9  etc.
   // --- fastNLO: fill variables for table header block A2
   fnloBlockA2 *A2 =  tabledQ2->GetBlockA2();

   // --- fastNLO user: up to 30 strings to describe the scenario
   A2->ScDescript.push_back("Trijets - d2sigma/dQ2 (pb/GeV2)");
   A2->ScDescript.push_back("H1 Collaboration");
   A2->ScDescript.push_back("Trijet cross section in DIS");
   A2->ScDescript.push_back("kT-algorithm, R=1, HighQ2, HERA-II");
   A2->ScDescript.push_back("H1prelim-11-032, Fig. 11a");

   A2->Ecms = sqrt(4.*920.*27.6);	// --- fastNLO user: set sqrt(s)
   A2->ILOord = 2;			// --- fastNLO user: power of LO contribution for process // 1 for dijets, 2 for trijets (nj=3)
   A2->NDim = 1;			// --- fastNLO user: No of dimensions in which observable is binned
   A2->DimLabel.push_back("Q^2");	// --- fastNLO user: label of 2nd dimension
   A2->IDiffBin.push_back(2);


   // --- fastNLO user: bin definitions - here in Q2 and ET
   const int nq2bins = 6;
   double q2bins[nq2bins+1] = {150.,200.,270.,400.,700.,5000.,15000.};

   // ---- initalize the bingrids and the normalizations ---- //
   A2->InitBinning( nq2bins , q2bins );

   // --- fastNLO table block B
   fnloBlockBNlojet *B = new fnloBlockBNlojet(tabledQ2->GetBlockA1(),tabledQ2->GetBlockA2());
   tabledQ2->CreateBlockB(0,B);

   // ---- initalize all constants to make fastNLO behave like for DIS scenario ---- //
   B->InitDISConstants(A2,nlo);


   // --- arrays for extreme x and (default) scale values (computed in Warm-Up run) ---- //
   double xlim[A2->NObsBin];
   double scale1lo[A2->NObsBin];
   double scale1hi[A2->NObsBin];
   double scale2lo[A2->NObsBin];
   double scale2hi[A2->NObsBin];
   for(int i=0;i<A2->NObsBin;i++){
      xlim[i] = 1.1e-07;
      scale1lo[i] = 1.0, scale1hi[i]=9.9e10;
      scale2lo[i] = 1.0, scale2hi[i]=9.9e10;
   }


   // --- fastNLO user: before running a new scenario for the first time,
   //     the following block (between "start" and "end") should be 
   //     completely removed. The first run must be a "Warm-Up Run" (IWarmUp=1)
   //     which produces an initialization block as output. This should be copied
   //     and pasted below. These initialization values must be used for all
   //     production jobs (IWarmUp=0) for a given scenario (otherwise the result 
   //     tables can not be merged). 
   // --- initialize variables for WarmUp run
   B->IWarmUp = bWarmup;
   if (bWarmup ) cout<<"\nThis is a warmup run.\n"<<endl;      // --- fastNLO user: do the Warm-Up run
   B->IWarmUpPrint = fNWarmupPrint; //3 000 0000
   //
   // --------- fastNLO: Warm-Up run results (start)
   // -474967296 contributions (!= events) in warm-up run
   xlim [0] = 8.959844e-03 , scale1lo [0] =   12.2474 , scale1hi [0] =   14.1421 , scale2lo [0] =    5.0194 , scale2hi [0] =   49.8321;
   xlim [1] = 9.718287e-03 , scale1lo [1] =   14.1421 , scale1hi [1] =   16.4317 , scale2lo [1] =    5.0215 , scale2hi [1] =   49.8526;
   xlim [2] = 1.071501e-02 , scale1lo [2] =   16.4317 , scale1hi [2] =   20.0000 , scale2lo [2] =    5.0217 , scale2hi [2] =   49.9151;
   xlim [3] = 1.261546e-02 , scale1lo [3] =   20.0000 , scale1hi [3] =   26.4575 , scale2lo [3] =    5.0132 , scale2hi [3] =   49.9038;
   xlim [4] = 1.701838e-02 , scale1lo [4] =   26.4575 , scale1hi [4] =   70.7107 , scale2lo [4] =    5.0140 , scale2hi [4] =   49.9694;
   xlim [5] = 7.852251e-02 , scale1lo [5] =   70.7107 , scale1hi [5] =  122.4745 , scale2lo [5] =    5.0077 , scale2hi [5] =   49.9713;

   // --------- fastNLO: Warm-Up run results (end)


   // ---- set number-of x-nodes ---- //
   B->SetNumberOfXNodesPerMagnitude(10,xlim);


   // ---- scale description ---- //
   B->SetScale1Name( "Q" );
   B->SetScale2Name( "<p_T>" );

   
   // ---- number of scale nodes for mu ---- //
   B->SetNumberOfScaleNodesScale1( 15 );
   B->SetNumberOfScaleNodesScale2( 12 );

   // ---- final initializations ---- //
   B->InitFinalDISValues( A2 , xlim , scale1lo , scale1hi, scale2lo , scale2hi ); 

}

void UserDIS::inittabledPt(string tablefilename){

   // --- set up fastNLO
   tabledPt = new fnloTable(tablefilename);

   // --- fastNLO: fill variable for table header block A1
   tabledPt->GetBlockA1()->SetScenName("H1-Trijets-dPt-HighQ2");  // NO SPACES HERE !!! // - fastNLO user: set scenario name
   tabledPt->GetBlockA1()->SetNcontrib(1);
   tabledPt->GetBlockA1()->SetNmult(0);
   tabledPt->GetBlockA1()->SetNdata(0);
   tabledPt->GetBlockA2()->SetIpublunits(12);  // - fastNLO user: set cross section units
                                            //                 (negative power of ten)
                                            //                 fb:15  pb:12  nb:9  etc.
   // --- fastNLO: fill variables for table header block A2
   fnloBlockA2 *A2 =  tabledPt->GetBlockA2();

   // --- fastNLO user: up to 30 strings to describe the scenario
   A2->ScDescript.push_back("Trijets - d2sigma/d<Pt> (pb/GeV)");
   A2->ScDescript.push_back("H1 Collaboration");
   A2->ScDescript.push_back("Trijet cross section in DIS");
   A2->ScDescript.push_back("kT-algorithm, R=1, HighQ2, HERA-II");
   A2->ScDescript.push_back("H1prelim-11-032, Fig. 11b");

   A2->Ecms = sqrt(4.*920.*27.6);	// --- fastNLO user: set sqrt(s)
   A2->ILOord = 2;			// --- fastNLO user: power of LO contribution for process // 1 for dijets, 2 for trijets (nj=3)
   A2->NDim = 1;			// --- fastNLO user: No of dimensions in which observable is binned
   A2->DimLabel.push_back("<p_T>");	// --- fastNLO user: label of 1st dimension
   A2->IDiffBin.push_back(2);


   // --- fastNLO user: bin definitions - here in Q2 and ET
   const int netbins = 3;
   double etbins[4] = {
      7.,11.,18.,30.
   };

   // ---- initalize the bingrids and the normalizations ---- //
   A2->InitBinning( netbins , etbins );

   // --- fastNLO table block B
   fnloBlockBNlojet *B = new fnloBlockBNlojet(tabledPt->GetBlockA1(),tabledPt->GetBlockA2());
   tabledPt->CreateBlockB(0,B);

   // ---- initalize all constants to make fastNLO behave like for DIS scenario ---- //
   B->InitDISConstants(A2,nlo);


   // --- arrays for extreme x and (default) scale values (computed in Warm-Up run) ---- //
   double xlim[A2->NObsBin];
   double scale1lo[A2->NObsBin];
   double scale1hi[A2->NObsBin];
   double scale2lo[A2->NObsBin];
   double scale2hi[A2->NObsBin];
   for(int i=0;i<A2->NObsBin;i++){
      xlim[i] = 1.1e-07;
      scale1lo[i] = 1.0, scale1hi[i]=9.9e10;
      scale2lo[i] = 1.0, scale2hi[i]=9.9e10;
   }


   // --- fastNLO user: before running a new scenario for the first time,
   //     the following block (between "start" and "end") should be 
   //     completely removed. The first run must be a "Warm-Up Run" (IWarmUp=1)
   //     which produces an initialization block as output. This should be copied
   //     and pasted below. These initialization values must be used for all
   //     production jobs (IWarmUp=0) for a given scenario (otherwise the result 
   //     tables can not be merged). 
   // --- initialize variables for WarmUp run
   B->IWarmUp = bWarmup;
   if (bWarmup ) cout<<"\nThis is a warmup run.\n"<<endl;      // --- fastNLO user: do the Warm-Up run
   B->IWarmUpPrint = fNWarmupPrint; //3 000 0000
   //
   // --------- fastNLO: Warm-Up run results (start)
   // -1014967296 contributions (!= events) in warm-up run
   xlim [0] = 8.959844e-03 , scale1lo [0] =   12.2474 , scale1hi [0] =  122.4745 , scale2lo [0] =    7.0000 , scale2hi [0] =   11.0000;
   xlim [1] = 1.756865e-02 , scale1lo [1] =   12.2474 , scale1hi [1] =  122.4745 , scale2lo [1] =   11.0000 , scale2hi [1] =   18.0000;
   xlim [2] = 4.370402e-02 , scale1lo [2] =   12.2474 , scale1hi [2] =  122.4745 , scale2lo [2] =   18.0000 , scale2hi [2] =   30.0000;
   // --------- fastNLO: Warm-Up run results (end)


   // ---- set number-of x-nodes ---- //
   B->SetNumberOfXNodesPerMagnitude(10,xlim);


   // ---- scale description ---- //
   B->SetScale1Name( "Q" );
   B->SetScale2Name( "<p_T>" );

   
   // ---- number of scale nodes for mu ---- //
   B->SetNumberOfScaleNodesScale1( 15 );
   B->SetNumberOfScaleNodesScale2( 12 );

   // ---- final initializations ---- //
   B->InitFinalDISValues( A2 , xlim , scale1lo , scale1hi, scale2lo , scale2hi ); 

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
