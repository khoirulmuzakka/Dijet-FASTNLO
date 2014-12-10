//
// fastNLO v2.2 creator code for inclusive jets
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
// in order to guarantee properly working code.
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
// Implementing a new scenario may imply to change:
//  - the jet algorithm ("#include" statement and assignment of "jetclus")
//  - the number of jets at LO in inputfunc (2-jet or 3-jet observable)
//  - the observable to compute
//  - the scale definition
//
// ================================================================

//------ DON'T TOUCH THIS PART! ------
#include <iostream>
#include <map>
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
#include "fnlo_int_nlojet/fj-jets.h"

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
   // --- fastNLO user: define the jet algorithm (for the choice of included header file above)
   fj_jets  jetclusfj;

   // --- define the jet structure
   bounded_vector<lorentzvector<double> > pj;

   // --- fastNLO definitions (not for user)
   double nevents;        // No. of events calculated so far
   unsigned long nwrite;  // No. of events after which to write out the table
};

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
   //   say::SetGlobalVerbosity(say::INFO);
   say::debug["inputfunc"] << "---------- inputfunc called ----------" << endl;
   // --- create fastNLO table and read in steering ... (if not done already)
   if (!ftable) {
      // --- fastNLO user: adapt the process constants to match the selected number of jets of the LO process, see below.
      //                   Either ProcConsts_HHC_2Jet or ProcConsts_HHC_3Jet
      ftable = new fastNLOCreate("InclusiveJets.str",UsefulNlojetTools::GenConsts(),UsefulNlojetTools::ProcConsts_HHC_2Jet());
   }

   // --- fastNLO user: select the number of jets of the LO process for your observable,
   //                   e.g. 2 for inclusive jets, 3 for 3-jet mass
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
   tabFilename += ".tab";
   ftable->SetFilename(tabFilename);
}

// --- fastNLO user: modify the jet selection in UserHHC::userfunc (default = cutting in |y| min, |y| max and pt min)
//                   (the return value must be true for jets to be UNselected)
struct fNLOSelector {
   fNLOSelector(double ymin, double ymax, double ptmin, bool pseudo=false):
      _ymin (ymin), _ymax (ymax), _ptmin (ptmin), _pseudo (pseudo){};
   double _ymin, _ymax, _ptmin;
   bool _pseudo;
   bool operator() (const lorentzvector<double> &a) {
      if (!_pseudo) return ! (_ymin <= abs(a.rapidity())  && abs(a.rapidity())  < _ymax && _ptmin <= a.perp());
      else          return ! (_ymin <= abs(a.prapidity()) && abs(a.prapidity()) < _ymax && _ptmin <= a.perp());
   };
};

// --- fastNLO user: modify the jet sorting in UserHHC::userfunc (default = descending in jet pt)
struct fNLOSorter {
   bool operator() (const lorentzvector<double> &a, const lorentzvector<double> &b) {return (a.perp() > b.perp());};
};

void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
   say::debug["UserHHC::userfunc"] << "---------- UserHHC::userfunc called ----------" << endl;

   // --- fastNLO user:
   //     Here is your playground where you compute your observable
   //     and the bin number ("obsbin"), which gets passed to
   //     fastNLO's table filling code.
   //     Usually, pT and E are in GeV, but this may be changed.
   //     ATTENTION: Scales must always be in GeV!

   // --- fastNLO user: check and get steering parameters once and store into static vars
   static std::map < std::string, bool > SteeringPars;

   // get general steering parameters needed here from steering file
   SteeringPars["FlexibleScaleTable"] = ftable->TestParameterInSteering("FlexibleScaleTable");
   static bool lFlexibleScaleTable = false; // default
   if ( SteeringPars["FlexibleScaleTable"] ) {
      ftable->GetParameterFromSteering("FlexibleScaleTable",lFlexibleScaleTable);
   }
   SteeringPars["DifferentialDimension"] = ftable->TestParameterInSteering("DifferentialDimension");
   static int read_ndim; // no default
   if ( SteeringPars["DifferentialDimension"] ) {
      ftable->GetParameterFromSteering("DifferentialDimension",read_ndim);
   } else {
      say::error["InclusiveJets"] << "Dimensioning of binning not set, aborted!" << endl;
      exit(1);
   }

   // definition of jet algorithm and overall phase space limits (no defaults)
   //
   // --- fastNLO user: set the jet algorithm and size via steering file
   // fastjet clustering jet algos: 0 = kT, 1 = CA, 2 = anti-kT
   // fastjet cone jet algos: 10 = SISCone, 11 = CDFMidPointCone, 12 = D0RunIICone
   SteeringPars["JetAlgo"] = ftable->TestParameterInSteering("JetAlgo");
   static int jetalgo;
   if ( SteeringPars["JetAlgo"] ) {
      ftable->GetParameterFromSteering("JetAlgo",jetalgo);
      if ( jetalgo < 0 || (2 < jetalgo && jetalgo < 10) || 12 < jetalgo ) {
         say::error["InclusiveJets"] << "Unknown jet algorithm " << jetalgo << ", aborted!" << endl;
         exit(1);
      }
   } else {
      say::error["InclusiveJets"] << "No jet algorithm selected, aborted!" << endl;
      exit(1);
   }
   SteeringPars["Rjet"] = ftable->TestParameterInSteering("Rjet");
   static double jetsize;
   if ( SteeringPars["Rjet"] ) {
      ftable->GetParameterFromSteering("Rjet",jetsize);
   } else {
      say::error["InclusiveJets"] << "Jet size R not defined, aborted!" << endl;
      exit(1);
   }
   SteeringPars["OvThr"] = ftable->TestParameterInSteering("OvThr");
   static double overlapthreshold = 0.5;
   if ( SteeringPars["OvThr"] ) {
      ftable->GetParameterFromSteering("OvThr",overlapthreshold);
   } else if ( jetalgo > 9 ) {
      say::error["InclusiveJets"] << "Overlap threshold not defined for jet algorithm " << jetalgo << ", aborted!" << endl;
      exit(1);
   }
   // --- fastNLO user: declare and initialize overall phase space cuts via steering file
   // overall lowest pT for jets to be considered
   SteeringPars["ptjmin"] = ftable->TestParameterInSteering("ptjmin");
   static double ptjmin;
   if ( SteeringPars["ptjmin"] ) {
      ftable->GetParameterFromSteering("ptjmin",ptjmin);
   } else {
      say::error["InclusiveJets"] << "Minimal jet pT (ptjmin) not defined, aborted!" << endl;
      exit(1);
   }
   // overall highest pT for jets not implemented, since uncritical with respect to CPU time consumption
   // overall smallest |(pseudo-)rapidity| for jets to be considered, use either y or eta but not both
   SteeringPars["yjmin"]   = ftable->TestParameterInSteering("yjmin");
   SteeringPars["etajmin"] = ftable->TestParameterInSteering("etajmin");
   static double yetajmin;
   if ( SteeringPars["yjmin"] && !SteeringPars["etajmin"] ) {
      ftable->GetParameterFromSteering("yjmin",yetajmin);
   } else if ( !SteeringPars["yjmin"] && SteeringPars["etajmin"] ) {
      ftable->GetParameterFromSteering("etajmin",yetajmin);
   } else {
      say::error["InclusiveJets"] << "Minimal jet (pseudo)rapidity (yjmin or etajmin) not uniquely defined, aborted!" << endl;
      exit(1);
   }
   // overall largest |(pseudo-)rapidity| for jets to be considered, use either y or eta but not both
   SteeringPars["yjmax"]   = ftable->TestParameterInSteering("yjmax");
   SteeringPars["etajmax"] = ftable->TestParameterInSteering("etajmax");
   static double yetajmax;
   if ( SteeringPars["yjmax"] && !SteeringPars["etajmax"] ) {
      ftable->GetParameterFromSteering("yjmax",yetajmax);
   } else if ( !SteeringPars["yjmax"] && SteeringPars["etajmax"] ) {
      ftable->GetParameterFromSteering("etajmax",yetajmax);
   } else {
      say::error["InclusiveJets"] << "Maximal jet (pseudo)rapidity (yjmax or etajmax) not uniquely defined, aborted!" << endl;
      exit(1);
   }
   // define logical for decision on cuts in (pseudo-)rapidity, no mixing allowed here
   static bool lpseudo;
   if ( SteeringPars["yjmin"] && SteeringPars["yjmax"] ) {
      lpseudo = false;
   } else if ( SteeringPars["etajmin"] && SteeringPars["etajmax"] ) {
      lpseudo = true;
   } else {
      say::error["InclusiveJets"] << "Phase space cuts mixed in (pseudo-)rapidity, aborted!" << endl;
      say::error["InclusiveJets"] << "Booleans for cut selections are" <<
         " yjmin "    << SteeringPars["yjmin"] <<
         ", yjmax "   << SteeringPars["yjmax"] <<
         ", etajmin " << SteeringPars["etajmin"] <<
         ", etajmax " << SteeringPars["etajmax"] << endl;
      say::error["InclusiveJets"] << "If you really want to mix, the code needs to be adapted." << endl;
      exit(1);
   }

   // apply the jet algorithm to partonic 4-vector array p of NLOJet++
   pj = jetclusfj(p,jetalgo,jetsize,overlapthreshold);
   unsigned int nj = pj.upper();

   // --- check on minimal and maximal no. of jets
   // ATTENTION: In principle, without cuts, there should always be two.
   //            For efficiency reasons though, our interface to the fastjet algorithms
   //            requires a minimal jet pT of 1 GeV. If this is a problem, the ptmin value
   //            in fj-jets.cc needs to be changed.
   // There should never be more than four jets in NLOJet++
   if (nj < 1) {
      say::info["InclusiveJets"] << "This event from NLOJet++ has no jets with pT > 1 GeV. Skipped!" << endl;
      return;
   } else if (nj > 4) {
      say::error["InclusiveJets"] << "This event from NLOJet++ has more than four jets, which should never happen. Aborted!" << endl;
      exit(1);
   }

   // --- give some debug output before selection and sorting
   if ( say::debug.GetSpeak() ) {
      for (unsigned int i=1; i<=nj; i++) {
         double pti  = pj[i].perp();
         double yi   = pj[i].rapidity();
         double etai = pj[i].prapidity();
         say::debug["InclusiveJets"] << "before cuts: jet # i, pt, y, eta: " << i << ", " << pti << ", " << yi << ", " << etai << endl;
      }
   }

   // --- select jets in y or eta and ptjmin (failing jets are moved to the end of the jet array pj!)
   fNLOSelector *SelJets;
   SelJets = new fNLOSelector(yetajmin,yetajmax,ptjmin,lpseudo);

   // --- count number of selected jets left at this stage
   size_t njet = std::remove_if(pj.begin(), pj.end(), *SelJets) - pj.begin();
   if ( njet < 1 ) return; // Nothing to be done

   // --- sort selected n jets at beginning of jet array pj, by default decreasing in pt
   static fNLOSorter SortJets;
   std::sort(pj.begin(), pj.begin() + njet, SortJets);

   // --- give some debug output after selection and sorting
   if ( say::debug.GetSpeak() ) {
      say::debug["InclusiveJets"] << "# jets before and after phase space cuts: nj, njet = " << nj << ", " << njet << endl;
      if ( ! lpseudo ) {
         say::debug["InclusiveJets"] << "phase space cuts: yjmin, yjmax, ptjmin: " << yetajmin << ", " << yetajmax << ", " << ptjmin << endl;
      } else {
         say::debug["InclusiveJets"] << "phase space cuts: etajmin, etajmax, ptjmin: " << yetajmin << ", " << yetajmax << ", " << ptjmin << endl;
      }
      for (unsigned int i=1; i<=njet; i++) {
         double pti  = pj[i].perp();
         double yi   = pj[i].rapidity();
         double etai = pj[i].prapidity();
         say::debug["InclusiveJets"] << "after cuts: jet # i, pt, y, eta: " << i << ", " << pti << ", " << yi << ", " << etai << endl;
      }
   }

   // --- set second choice (only for flexible tables) for the renormalization and factorization scale to max jet pT
   double ptmax  = pj[1].perp();
   double mu2 = ptmax;
   static vector<double> scalevars;
   if ( ! lFlexibleScaleTable ) scalevars = ftable->GetScaleVariations();

   // ---- fastNLO v2.2
   // Analyze inclusive jets in jet loop
   for (unsigned int i = 1; i <= njet; i++) {

      // Get jet quantities
      double pt  = pj[i].perp();
      double yeta;
      if ( ! lpseudo ) {
         yeta = abs(pj[i].rapidity());
      } else {
         yeta = abs(pj[i].prapidity());
      }
      double phi = atan2(pj[i].Y(), pj[i].X());

      // --- set first choice for the renormalization and factorization scale to jet pT
      double mu1 = pt;

      // get matrix elements
      static vector<fnloEvent> contribsflex;
      static vector< vector<fnloEvent> > contribsfix;
      if (lFlexibleScaleTable) {
         contribsflex = UsefulNlojetTools::GetFlexibleScaleNlojetContribHHC(p,amp);
      } else {
         contribsfix  = UsefulNlojetTools::GetFixedScaleNlojetContribHHC(p,amp,mu1,scalevars);
      }

      // scenario specific quantites
      fnloScenario scen;
      if ( read_ndim == 1 ) {        // 1D binning
         scen.SetObservableDimI( pt  , 0 );
      } else if ( read_ndim == 2 ) { // 2D binning
         scen.SetObservableDimI( yeta, 0 );
         scen.SetObservableDimI( pt  , 1 );
      } else if ( read_ndim == 3 ) { // 3D binning
         scen.SetObservableDimI( phi , 0 );
         scen.SetObservableDimI( yeta, 1 );
         scen.SetObservableDimI( pt  , 2 );
      } else {
         say::error["InclusiveJets"] << "More than 3D binning not implemented for inclusive jets, aborted!" << endl;
         say::error["InclusiveJets"] << "DifferentialDimension NDim = " << read_ndim << endl;
         exit(1);
      }
      scen.SetObsScale1( mu1 );   // must be consistent with 'mu' from contribs

      if (lFlexibleScaleTable) {
         scen.SetObsScale2( mu2 );
         ftable->FillAllSubprocesses(contribsflex,scen);
      } else {
         ftable->FillAllSubprocesses(contribsfix,scen);
      }
   }
   delete SelJets;
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
