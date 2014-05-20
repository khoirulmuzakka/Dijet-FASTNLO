//
// fastNLO v2 creator code for fnl5350eta0:
//     CMS LHC Inclusive Jets Scenario, E_cms = 5.02 TeV
//     for fastjet anti-kT algo with R=0.2, 0.3, 0.4 in E-scheme
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
//   inputfunc             (-> user edits)
//   InitFastNLO           (-> user edits)
//   struct fNLOSelector   (-> user edits)
//   struct fNLOSorter     (-> user edits)
//   UserHHC::userfunc     (-> user edits)
//   psinput               (don't touch)
//   userfunc              (don't touch)
//   UserHHC::initfunc     (don't touch)
//   UserHHC::end_of_event (don't touch)
//   UserHHC::phys_output  (don't touch)
//
// Implementing a new scenario requires to edit:
//  - the jet algorithm ("#include" statement and assignment of "jetclus")
//  - number of jets (inputfunc)
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
#include "fastnlotk/fastNLOCreate.h"
#include "fastnlotk/fastNLOEvent.h"

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
#include "fnlo_int_nlojet/fj-ak.h"

// --- fastNLO ---
#include "fnlo_int_nlojet/fnlo_int_hhc_nlojet.h"

// --- fastNLO v2.2
fastNLOCreate *ftable = NULL;
void InitFastNLO(const std::basic_string<char>& fname);

class UserHHC : public basic_user_set<user0d_hhc, user1h_hhc, user2h_hhc>
{
public:
   // --- init and user function
   void initfunc(unsigned int);
   void userfunc(const event_hhc&, const amplitude_hhc&);
   virtual void end_of_event();
   virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);

private:
   // --- fastNLO user: define the jet algorithm (consistent with the header file above)
   fj_ak jetclus;

   // --- define the jet structure
   bounded_vector<lorentzvector<double> > pj;

   // --- fastNLO definitions (not for user)
   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table
};

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
   say::debug["inputfunc"] << "---------- inputfunc called ----------" << endl;
   // --- create fastNLO table and read in steering ... (if not done already)
   if (!ftable) {
      ftable = new fastNLOCreate("fnl5350eta0v22.str",UsefulNlojetTools::GenConsts(),UsefulNlojetTools::ProcConsts_HHC_2Jet());
   }

   // --- fastNLO user: select the number of jets of the LO process for your observable,
   //                   e.g. 2 for incl. jets, 3 for 3-jet mass
   //nj = 1U;
   nj = 2U;
   //nj = 3U;

   // --- fastNLO user: number of the up and down type flavours (usually, you won´t change that)
   nu = 2U;
   nd = 3U;
}

void InitFastNLO(const std::basic_string<char>& __file_name)
{
   say::debug["InitFastNLO"] << "---------- InitFastNLO called ----------" << endl;
   // --- obtain relevant variables from NLOJet++ command line arguments
   ftable->SetOrderOfAlphasOfCalculation(UsefulNlojetTools::GetOrderOfRun(__file_name));

   // --- set fastNLO filename according to NLOJet++ command line arguments
   string tabFilename = __file_name.c_str();
   tabFilename += "_v22.tab";
   ftable->SetFilename(tabFilename);
}

// --- fastNLO user: modify the jet selection in UserHHC::userfunc (default = cutting in |y| min, |y| max and pt min)
//                   (the return value must be true for jets to be UNselected)
// fnl5350eta0: use pseudorapidity eta
struct fNLOSelector {
   fNLOSelector(double ymin, double ymax, double ptmin):
      _ymin (ymin), _ymax (ymax), _ptmin (ptmin){};
   double _ymin, _ymax, _ptmin;
   bool operator() (const lorentzvector<double> &a) {return ! (_ymin <= abs(a.prapidity()) && abs(a.prapidity()) < _ymax && _ptmin <= a.perp());};
};

// --- fastNLO user: modify the jet sorting in UserHHC::userfunc (default = descending in jet pt)
struct fNLOSorter {
   bool operator() (const lorentzvector<double> &a, const lorentzvector<double> &b) {return (a.perp() > b.perp());};
};

void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
   say::debug["UserHHC::userfunc"] << "---------- UserHHC::userfunc called ----------" << endl;

   // --- fastNLO user: in this scenario run jet algo with three different jet sizes R
   const unsigned int ndim1bins = 3;
   const double Rjet[ndim1bins] = { 0.2, 0.3, 0.4 };
   for (unsigned int k=0; k<ndim1bins; k++) {

      // --- fastNLO user: set the jet size and run the jet algorithm
      double jetsize = Rjet[k];
      pj = jetclus(p,jetsize);
      unsigned int nj = pj.upper();

      // --- give some debug output before selection and sorting
      if ( say::debug.GetSpeak() ) {
         for (unsigned int i=1; i<=nj; i++) {
            double pti  = pj[i].perp();
            double yi   = pj[i].rapidity();
            double etai = pj[i].prapidity();
            say::debug["fnl-scenario"] << "Before cuts: jet # i, pt, y, eta: " << i << ", " << pti << ", " << yi << ", " << etai << endl;
         }
      }

      // --- check on maximal no. of jets: 4 (should never be more in NLOJet++)
      if (nj > 4) {
         say::error["fnl-scenario"] << "This scenario is not suited for " << nj <<
            " jets. Aborted!" << endl;
         exit(1);
      }

      // --- fastNLO user:
      //     Here is your playground where you compute your observable
      //     and the bin number ("obsbin"), which gets passed to
      //     fastNLO's table filling code.
      //     Usually, pT and E are in GeV, but this may be changed.
      //     ATTENTION: Scales must always be in GeV!

      // --- declare and initialize phase space cut variables
      // can partially be taken from table binning for this scenario ???
      // smallest |(pseudo-)rapidity| for jets to be considered
      const double yjmin  = 0.0;
      // largest |(pseudo-)rapidity| for jets to be considered
      const double yjmax  = 0.3;
      // lowest pT for jets to be considered
      const double ptjmin = 22.; // A2->LoBin[0][0];

      // --- select jets in y or eta and ptjmin (failing jets are moved to the end of the jet array pj!)
      static fNLOSelector SelJets(yjmin,yjmax,ptjmin);
      // --- count number of selected jets left at this stage
      size_t njet = std::remove_if(pj.begin(), pj.end(), SelJets) - pj.begin();

      // --- sort selected n jets at beginning of jet array pj, by default decreasing in pt
      // fnl5350eta0: Not required for inclusive jets
      //  static fNLOSorter SortJets;
      //  std::sort(pj.begin(), pj.begin() + njet, SortJets);

      // --- give some debug output after selection and sorting
      if ( say::debug.GetSpeak() ) {
         say::debug["fnl-scenario"] << "# jets before and after phase space cuts: nj, njet = " << nj << ", " << njet << endl;
         say::debug["fnl-scenario"] << "phase space cuts: yjmin, yjmax, ptjmin: " << yjmin << ", " << yjmax << ", " << ptjmin << endl;
         for (unsigned int i=1; i<=njet; i++) {
            double pti  = pj[i].perp();
            double yi   = pj[i].rapidity();
            double etai = pj[i].prapidity();
            say::debug["fnl-scenario"] << "After cuts: jet # i, pt, y, eta: " << i << ", " << pti << ", " << yi << ", " << etai << endl;
         }
      }

      // ---- fastNLO v2.2
      // Analyze inclusive jets in jet loop
      const vector<double>& scalevars = ftable->GetScaleVariations();
      for (unsigned int i = 1; i <= njet; i++) {

         // Get jet quantities
         double pt  = pj[i].perp();

         // --- set the renormalization and factorization scale to jet pT
         double mu = pt;

         // get matrix elements
         //vector<fnloEvent> contribs = UsefulNlojetTools::GetFixedScaleNlojetContribHHC(p,amp,mu);
         vector<vector<fnloEvent> > contribs = UsefulNlojetTools::GetFixedScaleNlojetContribHHC(p,amp,mu,scalevars);

         // scenario specific quantites
         fnloScenario scen;
         scen.SetObservableDimI( Rjet[k] , 0 );
         scen.SetObservableDimI( pt , 1 );
         scen.SetObsScale1( mu );   // must be consistent with 'mu' from contribs
         ftable->FillAllSubprocesses(contribs,scen);
      }
   }
}



//------ DON'T TOUCH THIS PART! ------
void psinput(phasespace_hhc *ps, double& s) {
   say::debug["psinput"] << "---------- psinput called ----------" << endl;

   // --- set the center-of-mass energy squared as read from steering file
   s = pow(ftable->GetEcms(),2);

   // --- in principle alternative phase space generators can be used
   // --- we support only the default for now
   ps = 0;
}

user_base_hhc * userfunc() {
   say::debug["userfunc"] << "---------- userfunc called ----------" << endl;
   return new UserHHC;
}

void UserHHC::initfunc(unsigned int)
{
   say::debug["UserHHC::initfunc"] << "---------- UserHHC::initfunc called ----------" << endl;
   // --- Initialize event counters
   nevents = 0;
   // Set some defaults
   if (nwrite==0) nwrite = 1000000;
}

void UserHHC::end_of_event(){
   say::debug["UserHHC::end_of_event"] << "---------- UserHHC::end_of_event called ----------" << endl;
   nevents += 1;
   // --- store table
   if (( (unsigned long)nevents % nwrite)==0){
      ftable->SetNumberOfEvents(nevents);
      ftable->WriteTable();
   }
}

void UserHHC::phys_output(const std::basic_string<char>& __file_name,
                          unsigned long __save, bool __txt)
{
   say::debug["UserHHC::phys_output"] << "---------- UserHHC::phys_output called ----------" << endl;
   nwrite = __save;
   InitFastNLO(__file_name);
}
