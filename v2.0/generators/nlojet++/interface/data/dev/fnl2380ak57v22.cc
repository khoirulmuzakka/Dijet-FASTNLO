//
// fastNLO v2 creator code for fnl2380ak0507y0:
//     CMS LHC Inclusive Jets Scenario, E_cms = 7 TeV
//     for fastjet anti-kT algo with R=0.5 and 0.7 in E-scheme
//
//
// ============== fastNLO user: ===================================
// To create your own scenario, it is recommended to take
// this code, make a copy and edit the relevant changes.
// Important:
// Edit only those lines which are labeled as "fastNLO user"
// and refer to the documentation ("fastNLO creator code in
// NLOJet++") for a detailed explanation of the parameters
// and variables.
// If a code fragment is not explicitely labeled as "fastNLO user",
// it is likely that a modification will interfere with
// the fastNLO routines.
// Please keep the order of all statements in inittable
// in order to guarantee a working code.
//
// This file contains the following routines:
//   inputfunc    (-> user edits)
//   psinput      (-> user edits)
//   userfunc     (-> user edits)
//   inittable    (-> user edits)
//   initfunc     (don't touch)
//   writetable   (don't touch)
//   end_of_event (don't touch)
//   phys_output  (don't touch)
//   GetEcms      (don't touch)
//   GetNj        (don't touch)
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
#include "fastnlotk/fastNLOCreate.h"

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

// --- fastNLO user: include the header file for the jet algorithm
#include "fj-ak.h"
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
   long long int nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table
   double YJMIN;
   double YJMAX;
   double PTJMIN;

   // --- fastNLO user FYI: use of the -n option of NLOJet++
   //     The steering of these flags is encoded in the NLOJet++ run name (option -n name).
   //     If the name matches "deb", "ref", or "wrm", the respective flag is set to true.
   //
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
   s =  49000000.; // LHC First Run     7000 GeV
   //   You can use your own phase generator.
   //   Here we use the default.
   ps = 0;
}

// --- fastNLO user: modify the jet selection in userfunc (default = cutting in |y| min, |y| max and pt min)
//                   (the return value must be true for jets to be UNselected)
struct fNLOSelector {
   fNLOSelector(double ymin, double ymax, double ptmin):
      _ymin (ymin), _ymax (ymax), _ptmin (ptmin){};
   double _ymin, _ymax, _ptmin;
   bool operator() (const lorentzvector<double> &a) {return ! (_ymin <= abs(a.rapidity()) && abs(a.rapidity()) < _ymax && _ptmin <= a.perp());};
};

// --- fastNLO user: modify the jet sorting in userfunc (default = descending in jet pt)
struct fNLOSorter {
   bool operator() (const lorentzvector<double> &a, const lorentzvector<double> &b) {return (a.perp() > b.perp());};
};

void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{

   // --- fastNLO user: in this scenario run jet algo with two different jet sizes R
   static const double R[] = {0.5, 0.7};
   const vector<double> & scalevars = ftable->GetScaleVariations();

   //jet size R loop
   for (unsigned int i=0; i < sizeof(R)/sizeof(R[0]); i++) {
      // --- fastNLO user: set the jet size and run the jet algorithm
      pj = jetclus(p,R[i]);


      // --- select jets in y or eta and ptjmin (failing jets are moved to the end of the jet array pj!)
      static fNLOSelector SelJets(YJMIN,YJMAX,PTJMIN);
      // --- count number of selected jets left at this stage
      size_t njet = std::remove_if(pj.begin(), pj.end(), SelJets) - pj.begin();

      //jet loop
      for (unsigned int j = 1; j <= njet; j++) {
         // Get jet quantities
         double pt  = pj[j].perp();
         // --- preset the renormalization and factorization scale to jet pT
         // --- (can get overwritten below with jet pT bin center -> other scenario)
         double mu = pt;
         //vector<fnloEvent> contribs = UsefulNlojetTools::GetFixedScaleNlojetContribHHC(p,amp,mu);
         vector<vector<fnloEvent> > contribs = UsefulNlojetTools::GetFixedScaleNlojetContribHHC(p,amp,mu,scalevars);
         // scenario specific quantites
         fnloScenario scen;
         scen.SetObservableDimI( R[i] , 0 );
         scen.SetObservableDimI( pt , 1 );
         scen.SetObsScale1( mu );		// must be consistent with 'mu' from contribs
         ftable->FillAllSubprocesses(contribs,scen); 
      } // --- end: jet loop
   } // --- end: jet size loop
} // --- end: fastNLO user playground


//------ DON'T TOUCH THIS PART! ------
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
      printf ("fastNLO: No. events: %.3G writing table ...\n",(double)nevents);
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
   ftable = new fastNLOCreate("fnl2380ak57v22.str");

   // obtain relevant variables from nlojet
   ftable->SetEcms(UsefulNlojetTools::GetEcms());
   ftable->SetLoOrder(UsefulNlojetTools::GetNj());
   ftable->SetOrderOfAlphasOfCalculation(UsefulNlojetTools::GetOrderOfRun(__file_name));

   // set filename, which is specified through command line
   string tabFilename = __file_name.c_str();
   tabFilename += "_v22.tab";
   ftable->SetFilename(tabFilename);

   // Cuts on jets
   // Values are read vom steering file
   // smallest |(pseudo-)rapidity| for jets to be considered
   YJMIN = DOUBLE(YJMIN);
   // largest |(pseudo-)rapidity| for jets to be considered
   YJMAX = DOUBLE(YJMAX);
   // lowest pT for jets to be considered
   PTJMIN = DOUBLE(PTJMIN);


   // give information to hb.
   //ftable->Print();
}

