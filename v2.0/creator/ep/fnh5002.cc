//
//   fastNLO v2.1 author code 
//   scenario: fnh5002
//   H1 Dijets (@high Q2)
//   for kT algorithm
//   HERA-II
//   H1prelim-11-032
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
//   writetable  (don't touch)
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
// 2011/11/22 DB update to fastNLO v2.1
// 2010/12/20 DB change to DESY-09-032 dijets module
// 2010/10/06 MW final fnh1001 version
// 2010/09/28 MW make user-friendly
// 2009/01/15 TK make code V2.0 compatible
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
   fnloTable *table;
   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table
   string tablefilename;  // The table file to write to
   time_t start_time;
   
   bool nlo;
   void inittable();
   void writetable();
};
  
user_base_dis * userfunc() {
  return new UserDIS;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
  // --- fastNLO user: select the number of jets for your observable
  nj = 2U;
  //nj = 3U;
  
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
   fnloBlockA2 *A2 =  table->GetBlockA2();

   double alem = 0; // for running alpha EM
   // --- x where PDF is probed ---
   double x = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);
   double Q2 = -((p[-1] - p[-2]).mag2());
   alem = xalpha_em(Q2);


   // --- fastNLO user:
   //     Here is your playground where you compute your observable 
   //     and the bin number ("obsbin") which gets passed to
   //     fastNLO's table filling code.
   //     (all pT, ET and E are in GeV)

   // --- declare and initialize phase space cut variables
   //     here: H1 cuts
   double pt = 0.0;
   //double pTmin1 = 7.0;
   double pTmin = 5.0;
   double pTmax = 50.0;
   double etamin = -1.0 , etamax = 2.5;

   double pt_mean = 0.0, pt_sum = 0.0;
   unsigned int jetnum = 0;
   lorentzvector<double> psum;

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
     if(pt >= pTmin && pt < pTmax && etamin <= etalab && etalab < etamax) {
       jetnum++;
       PtArray[jetnum] = pt; // starting from 1
       Theta_pj[jetnum] = pjb[i].theta();
       IndexArray[jetnum] = i;
    }
   }
  const unsigned int nSim = 2;
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
	// and sort the Theta array.
	t = Theta_pj[y];
	Theta_pj[y] = Theta_pj[y+1];
	Theta_pj[y+1] = t;
      }
    }          
  }
  
  pt_mean = ( PtArray[1] + PtArray[2] ) / 2;
  pt_sum = PtArray[1] + PtArray[2];
  psum = pjb[IndexArray[1]]+ pjb[IndexArray[2]];

  double Mn = psum.mag();

  if ( Mn < 16) return;

  int bin = -1;

  //if ( pt_mean<7 || pt_mean>50 ) return;

  // loop over all jets
  for(int j = 0; j < A2->GetNObsBin(); j++) {
    if (pt_mean >= A2->LoBin[j][0]  && pt_mean <  A2->UpBin[j][0] &&
	Q2 >= A2->LoBin[j][1]  && Q2 <  A2->UpBin[j][1]) {
      bin=j;
      break;
    }
  }
  
  //  ------------------------
  //
  //    set the scale here !
  //  
  //double mu = Q;  
  //   double mirkesscale = std::pow ( ( std::sqrt(mirkesscale_jet1) + std::sqrt(mirkesscale_jet2) ) / 2 , 2 ) ; //  mu = (ktb1 + ktb2)/2
  //   double mu = std::sqrt(mirkesscale);

  //   double mirkesscale_jet1 = 2 * PtArray[1]*PtArray[1] / ( 1 + std::cos(Theta_pj[1]) );
  //   double mirkesscale_jet2 = 2 * PtArray[2]*PtArray[2] / ( 1 + std::cos(Theta_pj[2]) );
  //   double mu = ( std::sqrt(mirkesscale_jet1) + std::sqrt(mirkesscale_jet2) ) / 2 ,;
  double mu = sqrt ( ( Q2 +  pt_mean*pt_mean ) / 2  );
  
  //---------- fill fastNLO arrays - don't touch this piece of code!
  if ( bin>=0 ) {
    double prefactor = 1./A2->BinSize[bin] * alem*alem; 
    for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
       if ( table->GetBlockB(k)->NScaleDep == 0 ){
	  if(table->GetBlockB(k)->GetIRef()>0){
	     // --- fastNLO user: don't modify the following calls!
	     ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventDIS(bin,x,mu,amp,pdf,prefactor);
	  }else{
	     ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventDIS(bin,x,mu,amp,dummypdf,prefactor);
	  }
       }
       else if ( table->GetBlockB(k)->NScaleDep == 3 ){
	  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventDISMuVar(bin,x,sqrt(Q2),pt_mean,amp,dummypdf,pdf,prefactor);
       }
    }
  }
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
   inittable();

}

void UserDIS::inittable(){

   // --- set up fastNLO
   table = new fnloTable(tablefilename);

   // --- fastNLO: fill variable for table header block A1
   table->GetBlockA1()->SetScenName("H1-Dijets-HighQ2");  // NO SPACES HERE !!!  - fastNLO user: set scenario name
   table->GetBlockA1()->SetNcontrib(1);
   table->GetBlockA1()->SetNmult(0);
   table->GetBlockA1()->SetNdata(0);
   table->GetBlockA2()->SetIpublunits(12);  // - fastNLO user: set cross section units
                                            //                 (negative power of ten)
                                            //                 fb:15  pb:12  nb:9  etc.
   // --- fastNLO: fill variables for table header block A2
   fnloBlockA2 *A2 =  table->GetBlockA2();

   // --- fastNLO user: up to 30 strings to describe the scenario
   A2->ScDescript.push_back("Dijets - d2sigma/d<pT>dQ2 (pb/GeV3) ");
   A2->ScDescript.push_back("H1 Collaboration");
   A2->ScDescript.push_back("Dijet cross sections in DIS");
   A2->ScDescript.push_back("kT-algorithm R=1, HighQ2, HERA-II");
   A2->ScDescript.push_back("H1prelim-11-032, Fig. 9");

   A2->Ecms = sqrt(4.*920.*27.6);	// --- fastNLO user: set sqrt(s)
   A2->ILOord = 1;    // --- fastNLO user: power of LO contribution for process
   A2->NDim = 2;      // --- fastNLO user: No of dimensions in which observable is binned
   A2->DimLabel.push_back("<p_T>");	// --- fastNLO user: label of 1st dimension
   A2->IDiffBin.push_back(2);		// --- fastNLO user: Is the cross section in the publication divided by this dimension (2: yes, 1:no)
   A2->DimLabel.push_back("Q^2");	// --- fastNLO user: label of 2nd dimension
   A2->IDiffBin.push_back(2);


   // --- fastNLO user: bin definitions - here in Q2 and ET
   const int nq2bins = 6;
   double q2bins[nq2bins+1] = {150.,200.,270.,400.,700.,5000.,15000.};

   const int netbins[nq2bins] = {4, 4, 4, 4, 4, 4};
   double etbins[nq2bins][5] = {
     {7.,11.,18.,30.,50.},
     {7.,11.,18.,30.,50.},
     {7.,11.,18.,30.,50.},
     {7.,11.,18.,30.,50.},
     {7.,11.,18.,30.,50.},
     {7.,11.,18.,30.,50.}
   };


   // ---- initalize the bingrids and the normalizations ---- //
   vector<double*> vetbins(nq2bins);
   for (unsigned int i=0;i<vetbins.size();i++) vetbins[i]=etbins[i]; 
   A2->InitBinning( nq2bins , q2bins , netbins , vetbins );


   // --- fastNLO table block B
   fnloBlockBNlojet *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
   table->CreateBlockB(0,B);
   B->InitDISConstants(A2,nlo);


   // --- initialize variables for WarmUp run
   //B->IWarmUp = 1;     cout<<"\nThis is a warmup run.\n"<<endl;      // --- fastNLO user: do the Warm-Up run
   B->IWarmUp = 0;     //                   or do production run(s)
   
   // - fastNLO user: remember to disable reference-mode in
   //                 Warm-Up run: "doReference = false" (above)
   B->IWarmUpPrint = 1000000; //3 000 0000


   // --- arrays for extreme x and (default) scale values (computed in Warm-Up run)
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
   //
   // --------- fastNLO: Warm-Up run results (start)
   xlim[0]=0.00575, scale1lo [0] =   12.2474 , scale1hi [0] =   14.1421 , scale2lo [0] =    7.0000 , scale2hi [0] =   11.0000;
   xlim[1]=0.00902, scale1lo [1] =   12.2474 , scale1hi [1] =   14.1421 , scale2lo [1] =   11.0000 , scale2hi [1] =   18.0000;
   xlim[2]=0.02051, scale1lo [2] =   12.2474 , scale1hi [2] =   14.1421 , scale2lo [2] =   18.0000 , scale2hi [2] =   30.0000;
   xlim[3]=0.05308, scale1lo [3] =   12.2474 , scale1hi [3] =   14.1421 , scale2lo [3] =   30.0000 , scale2hi [3] =   50.0000;
   xlim[4]=0.00649, scale1lo [4] =   14.1421 , scale1hi [4] =   16.4317 , scale2lo [4] =    7.0000 , scale2hi [4] =   11.0000;
   xlim[5]=0.00973, scale1lo [5] =   14.1421 , scale1hi [5] =   16.4317 , scale2lo [5] =   11.0000 , scale2hi [5] =   18.0000;
   xlim[6]=0.02136, scale1lo [6] =   14.1421 , scale1hi [6] =   16.4317 , scale2lo [6] =   18.0000 , scale2hi [6] =   30.0000;
   xlim[7]=0.05362, scale1lo [7] =   14.1421 , scale1hi [7] =   16.4317 , scale2lo [7] =   30.0000 , scale2hi [7] =   50.0000;
   xlim[8]=0.00743, scale1lo [8] =   16.4317 , scale1hi [8] =   20.0000 , scale2lo [8] =    7.0000 , scale2hi [8] =   11.0000;
   xlim[9]=0.01067, scale1lo [9] =   16.4317 , scale1hi [9] =   20.0000 , scale2lo [9] =   11.0000 , scale2hi [9] =   18.0000;
   xlim[10]=0.02223, scale1lo [10] =   16.4317 , scale1hi [10] =   20.0000 , scale2lo [10] =   18.0000 , scale2hi [10] =   30.0000;
   xlim[11]=0.05467, scale1lo [11] =   16.4317 , scale1hi [11] =   20.0000 , scale2lo [11] =   30.0000 , scale2hi [11] =   50.0000;
   xlim[12]=0.00933, scale1lo [12] =   20.0000 , scale1hi [12] =   26.4575 , scale2lo [12] =    7.0000 , scale2hi [12] =   11.0000;
   xlim[13]=0.01261, scale1lo [13] =   20.0000 , scale1hi [13] =   26.4575 , scale2lo [13] =   11.0000 , scale2hi [13] =   18.0000;
   xlim[14]=0.02419, scale1lo [14] =   20.0000 , scale1hi [14] =   26.4575 , scale2lo [14] =   18.0000 , scale2hi [14] =   30.0000;
   xlim[15]=0.05668, scale1lo [15] =   20.0000 , scale1hi [15] =   26.4575 , scale2lo [15] =   30.0000 , scale2hi [15] =   50.0000;
   xlim[16]=0.01358, scale1lo [16] =   26.4575 , scale1hi [16] =   70.7107 , scale2lo [16] =    7.0000 , scale2hi [16] =   11.0000;
   xlim[17]=0.01682, scale1lo [17] =   26.4575 , scale1hi [17] =   70.7107 , scale2lo [17] =   11.0000 , scale2hi [17] =   18.0000;
   xlim[18]=0.02842, scale1lo [18] =   26.4575 , scale1hi [18] =   70.7107 , scale2lo [18] =   18.0000 , scale2hi [18] =   30.0000;
   xlim[19]=0.06108, scale1lo [19] =   26.4575 , scale1hi [19] =   70.7107 , scale2lo [19] =   30.0000 , scale2hi [19] =   50.0000;
   xlim[20]=0.07445, scale1lo [20] =   70.7107 , scale1hi [20] =  122.4745 , scale2lo [20] =    7.0000 , scale2hi [20] =   11.0000;
   xlim[21]=0.07779, scale1lo [21] =   70.7107 , scale1hi [21] =  122.4745 , scale2lo [21] =   11.0000 , scale2hi [21] =   18.0000;
   xlim[22]=0.08915, scale1lo [22] =   70.7107 , scale1hi [22] =  122.4745 , scale2lo [22] =   18.0000 , scale2hi [22] =   30.0000;
   xlim[23]=0.12192, scale1lo [23] =   70.7107 , scale1hi [23] =  122.4745 , scale2lo [23] =   30.0000 , scale2hi [23] =   50.0000;
   // --------- fastNLO: Warm-Up run results (end)



   // ---- set number-of x-nodes ---- //
   B->SetNumberOfXNodesPerMagnitude(20,xlim);


   // ---- scale description ---- //
   B->SetScale1Name( "Q" );
   B->SetScale2Name( "<p_T>" );

   
   // ---- number of scale nodes for mu ---- //
   B->SetNumberOfScaleNodesScale1( 20 );
   B->SetNumberOfScaleNodesScale2( 15 );


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

   
   // reference table
   if(doReference){
      fnloBlockBNlojet *refB = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
      table->CreateBlockB(1,refB);
      refB->Copy(table->GetBlockB(0));
      refB->InitReferenceTable(A2);
      table->GetBlockA1()->SetNcontrib(2);
   }

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
