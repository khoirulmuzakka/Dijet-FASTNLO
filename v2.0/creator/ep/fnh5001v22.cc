//
// fastNLO v2.1 author code for fnh5001v22:
//     H1, High Q2, HERA-II, Inclusive jets Scenario, kt
//
// 
// ===================== fastNLO user =============================
// To create your own DIS scenario, it is recommended to take 
// this code, make a copy and edit the relevant changes.
//
// This file contains the following routines:
//   inputfunc		(-> user edits)
//   psinput		(-> user edits)
//   userfunc		(-> user edits)
//   initfunc		(don't touch)
//   end_of_event	(don't touch)
//   phys_output	(don't touch)
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
   //pdf_cteq6dis pdf;
   pdf_dis_dummy dummypdf;

   // jet algorithms
   kT_clus_long jetclus;
   //fj_ak fjclus; 
   
   // the jet structore in breit & lab. frame
   bounded_vector<lorentzvector<double> > pjb, pjl; 
   bounded_vector<unsigned int> jet;
   
   //  event in breit frame
   event_dis pbreit;
  
   // boost the jet momenta back to the laboratory frame
   void boost_back_to_lab(const event_dis&);
   
   time_t start_time;
   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table

   // fastNLO
   fastNLOCreate *ftable;
   void InitFastNLO(const std::basic_string<char>& fname);
   
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

// --------- alpha_em running --------- //
extern"C" double xalfaem_(double *);	// 'old' nlojet++ alpha_em-routine
double xalpha_em(double mq2) { return xalfaem_(&mq2); }
// -------- Spiesberger 2012 ---------- //
// extern"C" double aemrun_(double *);
// extern"C" void eprc_init_(int *);
// extern"C" void setpar_(int *);
// extern"C" void readewpar_();
// double xalpha_em(double mq2) {
//    return aemrun_(&mq2);
// }
// ------------------------------------ //


void UserDIS::userfunc(const event_dis& p, const amplitude_dis& amp) {
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


   // --- fastNLO v2.2 ----  nlojet-event ----
   vector<fnloEvent> contribs = UsefulNlojetTools::GetFlexibleScaleNlojetContribDIS(p,amp,dummypdf,alem);

   // ---- Inclusive jets: loop over all jets ----
   for ( unsigned int ij = 0 ; ij<jetnum ; ij++ ){
      // fastNLO v2.2: scenario specific quantites
      fnloScenario scen;
      scen.SetObservableDimI( Q2 , 0 );
      scen.SetObservableDimI( PtArray[ij+1] , 1 );
      scen.SetObsScale1( Q );		// Scale1 MUST be Q !
      scen.SetObsScale2( PtArray[ij+1] );
      
      ftable->FillAllSubprocesses(contribs,scen); 
      
   }   
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
      cout.flush();

      printf("fastNLO: Writing fastNLO table.\n");
      ftable->SetNumberOfEvents(nevents);
      ftable->WriteTable();    
      printf("fastNLO: New table written.\n");
   }
}

void UserDIS::phys_output(const std::basic_string<char>& __file_name, 
                          unsigned long __save, bool __txt) 
{
   // fastNLO v2.2
   nwrite = __save;
   InitFastNLO(__file_name);
}

//____________________ fastNLO v2.2 ____________________________
void UserDIS::InitFastNLO(const std::basic_string<char>& __file_name)
{
   // create table and read in steering...
   cout<<"\n ---------------------------------------------------------------\n"<<endl;
   ftable = new fastNLOCreate("fnh5001v22.str");

   // obtain relevant variables from nlojet
   ftable->SetEcms(UsefulNlojetTools::GetEcms());
   ftable->SetLoOrder(UsefulNlojetTools::GetLoOrder());
   ftable->SetOrderOfAlphasOfCalculation(UsefulNlojetTools::GetOrderOfRun(__file_name));

   // set filename, which is specified through command line
   string tabFilename = __file_name.c_str();
   tabFilename += ".tab";
   ftable->SetFilename(tabFilename);

   // give information to hb.
   //ftable->Print();
}
