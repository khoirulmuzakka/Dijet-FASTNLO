//
// fastNLO v2 creator code for fnl2352:
//     ATLAS LHC Inclusive Jets Scenario, E_cms = 7 TeV
//     for fastjet anti-kT algo with R=0.6 in E-scheme
//
//
// ============== fastNLO user: ===================================
// To create your own scenario, it is recommended to take
// this code, make a copy and edit the relevant changes.
//
// This file contains the following routines:
//   inputfunc    (-> user edits)
//   psinput      (-> user edits)
//   userfunc     (-> user edits)
//   initfunc     (don't touch)
//   end_of_event (don't touch)
//   phys_output  (don't touch)
//
// Implementing a new scenario requires to edit:
//  - the jet algorithm ("#include" statement and assignment of "jetclus")
//  - number of jets (inputfunc)
//  - center-of-mass energy (psinput)
//  - compute observable
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

// ---- fastNLO ----
#include <fastnlotk/fastNLOCreate.h>

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

// --- fastNLO user: include the header file for the jet algorithm
#include "fj-ak.h"

// ---- fastNLO ----
#include "fastNLOInterfaceToNLOJET.cc"

class UserHHC : public basic_user_set<user0d_hhc, user1h_hhc, user2h_hhc>
{
public:
   //   init and user function
   void initfunc(unsigned int);
   void userfunc(const event_hhc&, const amplitude_hhc&);
   virtual void end_of_event();
   virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);

private:
   // --- fastNLO user: define the jet algorithm (consistent with the header file above)
   fj_ak jetclus;

   bounded_vector<lorentzvector<double> > pj;    // the jet structure

   // --- fastNLO definitions (not for user)
   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table

   // --- fastNLO v2.2
   fastNLOCreate *ftable;
   void InitFastNLO(const std::basic_string<char>& fname);
};

user_base_hhc * userfunc() {
   return new UserHHC;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
   // --- fastNLO user: select the number of jets of the LO process for your observable,
   //                   e.g. 2 for incl. jets, 3 for 3-jet mass
   nj = 2U;
   //nj = 3U;

   // --- number of the up and down type flavours (don't touch)
   nu = 2U;
   nd = 3U;
}

void psinput(phasespace_hhc *ps, double& s)
{
   // --- fastNLO user: set the total c.m. energy squared in GeV^2
   s =  49000000.;	// LHC First Run     7000 GeV

   //   You can use your own phase generator.
   //   Here we use the default.
   ps = 0;
}

void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{

   // --- Set the jet size and run the jet algorithm
   double jetsize = 0.6;
   pj = jetclus(p,jetsize);
   unsigned int nj = pj.upper();

   // ---- fastNLO v2.2
   // Analyze inclusive jets in jet loop
   const vector<double>& scalevars = ftable->GetScaleVariations();
   for (unsigned int i = 1; i <= pj.size(); i++) {

      // Get jet quantities
      double pt  = pj[i].perp();
      double rap = abs(pj[i].rapidity());

      // get matrix elements
      double mu = pt;
      //vector<fnloEvent> contribs = UsefulNlojetTools::GetFixedScaleNlojetContribHHC(p,amp,mu);
      vector<vector<fnloEvent> > contribs = UsefulNlojetTools::GetFixedScaleNlojetContribHHC(p,amp,mu,scalevars);

      // scenario specific quantites
      fnloScenario scen;
      scen.SetObservableDimI( rap , 0 );
      scen.SetObservableDimI( pt , 1 );
      scen.SetObsScale1( mu );		// must be consistent with 'mu' from contribs
      ftable->FillAllSubprocesses(contribs,scen); 
   }     
}



void UserHHC::initfunc(unsigned int)
{
   // --- Initialize event counters
   nevents = 0;
   // Set some defaults
   if (nwrite==0) nwrite = 5000000;
}

void UserHHC::end_of_event(){
   nevents += 1;
   // --- store table
   if (( (unsigned long)nevents % nwrite)==0){
      printf ("fastNLO: No. events: %.3G writing table ...\n",nevents);
      ftable->SetNumberOfEvents(nevents);
      ftable->WriteTable();    
   }
}

void UserHHC::phys_output(const std::basic_string<char>& __file_name,
                          unsigned long __save, bool __txt)
{
   nwrite = __save;
   InitFastNLO(__file_name);
}


void UserHHC::InitFastNLO(const std::basic_string<char>& __file_name)
{
   // create table and read in steering...
   cout<<"\n ---------------------------------------------------------------\n"<<endl;
   ftable = new fastNLOCreate("fnl2352v22.str");

   // obtain relevant variables from nlojet
   ftable->SetEcms(UsefulNlojetTools::GetEcms());
   ftable->SetLoOrder(UsefulNlojetTools::GetNj());
   ftable->SetOrderOfAlphasOfCalculation(UsefulNlojetTools::GetOrderOfRun(__file_name));

   // set filename, which is specified through command line
   string tabFilename = __file_name.c_str();
   tabFilename += "_v22.tab";
   ftable->SetFilename(tabFilename);

   // give information to hb.
   //ftable->Print();
}
