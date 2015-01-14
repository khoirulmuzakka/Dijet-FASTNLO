//
// fastNLO v2.2 creator code for inclusive dijet events
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
#include "fastnlotk/read_steer.h"

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
   // 'nsave' is initialized by NLOJet++ with 10000 or to the command line number given via --save-after
   // in fastNLO this number is logarithmically increased by factors of 10 up to nwritemax = 10M
   virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);

private:
   // --- fastNLO user: define the jet algorithm (for the choice of included header file above)
   fj_jets jetclusfj;

   // --- define the jet structure
   bounded_vector<lorentzvector<double> > pj;

   // --- fastNLO definitions (not for user)
   double nevents;           // No. of events calculated so far
   unsigned long nwrite;     // Actual no. of events after which to write out the table
   unsigned long nwritemax;  // Maximal no. of events after which to write out the table
};

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
   //   say::SetGlobalVerbosity(say::INFO);
   say::debug["inputfunc"] << "---------- inputfunc called ----------" << endl;
   // --- create fastNLO table and read in steering ... (if not done already)
   if (!ftable) {
      // --- fastNLO user: adapt the process constants to match the selected number of jets of the LO process, see below.
      //                   Either ProcConsts_HHC_2Jet or ProcConsts_HHC_3Jet
      ftable = new fastNLOCreate("InclusiveDijetEvents.str",UsefulNlojetTools::GenConsts(),UsefulNlojetTools::ProcConsts_HHC_2Jet());
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
   static int NDim; // no default
   if ( SteeringPars["DifferentialDimension"] ) {
      ftable->GetParameterFromSteering("DifferentialDimension",NDim);
   } else {
      say::error["InclusiveDijetEvents"] << "Dimensioning of binning not set, aborted!" << endl;
      exit(1);
   }
   // dimension labels
   SteeringPars["DimensionLabels"] = ftable->TestParameterInSteering("DimensionLabels");
   vector<string> DimLabel; // no default
   DimLabel.resize(NDim);
   if ( SteeringPars["DimensionLabels"] ) {
      ftable->GetParameterFromSteering("DimensionLabels",DimLabel);
      //      cout << "DimLabel[0] = " << DimLabel[0] << ", DimLabel[1] = " << DimLabel[1] << endl;
   } else {
      say::error["InclusiveDijetEvents"] << "Dimension labels not set, aborted!" << endl;
      exit(1);
   }
   // dimension labels define the observables
   enum Obs0 { YMAX, YSTAR };
   enum Obs1 { MJJGEV, MJJTEV };
   Obs0 obs0def;
   Obs1 obs1def;
   if ( DimLabel[0] == "|y_max|" ) {
      obs0def = YMAX;
   } else if ( DimLabel[0] == "y_star" ) {
      obs0def = YSTAR;
   } else {
      say::error["InclusiveDijetEvents"] << "Unknown observable, i.e. dimension label, aborted!" << endl;
      say::error["InclusiveDijetEvents"] << "DimLabel[0] = " << DimLabel[0] << endl;
      say::error["InclusiveDijetEvents"] << "Please complement this scenario to include the requested observable." << endl;
      exit(1);
   }
   if ( DimLabel[1] == "Mjj_[GeV]" ) {
      obs1def = MJJGEV;
   } else if ( DimLabel[1] == "Mjj_[TeV]" ) {
      obs1def = MJJTEV;
   } else {
      say::error["InclusiveDijetEvents"] << "Unknown observable, i.e. dimension label, aborted!" << endl;
      say::error["InclusiveDijetEvents"] << "DimLabel[1] = " << DimLabel[1] << endl;
      say::error["InclusiveDijetEvents"] << "Please complement this scenario to include the requested observable." << endl;
      exit(1);
   }
   // scale descriptions
   SteeringPars["ScaleDescriptionScale1"] = ftable->TestParameterInSteering("ScaleDescriptionScale1");
   vector<string> ScaleLabel; // no default
   ScaleLabel.resize(2);
   if ( SteeringPars["ScaleDescriptionScale1"] ) {
      ftable->GetParameterFromSteering("ScaleDescriptionScale1",ScaleLabel[0]);
   } else {
      say::error["InclusiveDijetEvents"] << "No description of scale 1, aborted!" << endl;
      exit(1);
   }
   SteeringPars["ScaleDescriptionScale2"] = ftable->TestParameterInSteering("ScaleDescriptionScale2");
   if ( SteeringPars["ScaleDescriptionScale2"] ) {
      ftable->GetParameterFromSteering("ScaleDescriptionScale2",ScaleLabel[1]);
   } else {
      ScaleLabel[1] = "-";
      say::warn["InclusiveDijetEvents"] << "No description of scale 2, flexible-scale tables not possible!" << endl;
   }
   // scale descriptions define the scales
   enum Scales { PTMAX, PT12AVE, PT123AVE, MJJHALF, PTMAXEXPYSTAR };
   Scales mu1def;
   Scales mu2def = PTMAX;
   if ( ScaleLabel[0] == "pT_max_[GeV]" ) {
      mu1def = PTMAX;
   } else if ( ScaleLabel[0] == "<pT_1,2>_[GeV]" ) {
      mu1def = PT12AVE;
   } else if ( ScaleLabel[0] == "<pT_1,2,3>_[GeV]" ) {
      mu1def = PT123AVE;
   } else if ( ScaleLabel[0] == "Mjj/2_[GeV]" ) {
      mu1def = MJJHALF;
   } else if ( ScaleLabel[0] == "pT_max*exp(0.3*y_star)_[GeV]" ) {
      mu1def = PTMAXEXPYSTAR;
   } else {
      say::error["InclusiveDijetEvents"] << "Unknown scale no. 1, i.e. scale description, aborted!" << endl;
      say::error["InclusiveDijetEvents"] << "ScaleLabel[0] = " << ScaleLabel[0] << endl;
      say::error["InclusiveDijetEvents"] << "Please complement this scenario to include the requested scale." << endl;
      exit(1);
   }
   if ( ScaleLabel[1] == "pT_max_[GeV]" ) {
      mu2def = PTMAX;
   } else if ( ScaleLabel[1] == "<pT_1,2>_[GeV]" ) {
      mu2def = PT12AVE;
   } else if ( ScaleLabel[1] == "<pT_1,2,3>_[GeV]" ) {
      mu2def = PT123AVE;
   } else if ( ScaleLabel[1] == "Mjj/2_[GeV]" ) {
      mu2def = MJJHALF;
   } else if ( ScaleLabel[1] == "pT_max*exp(0.3*y_star)_[GeV]" ) {
      mu2def = PTMAXEXPYSTAR;
   } else {
      say::warn["InclusiveDijetEvents"] << "Unknown scale no. 2, i.e. scale description!" << endl;
      say::warn["InclusiveDijetEvents"] << "ScaleLabel[1] = " << ScaleLabel[1] << endl;
   }

   // definition of jet algorithm and jet phase space limits (no defaults)
   //
   // --- fastNLO user: set the jet algorithm and size via steering file
   // fastjet clustering jet algos: 0 = kT, 1 = CA, 2 = anti-kT
   // fastjet cone jet algos: 10 = SISCone, 11 = CDFMidPointCone, 12 = D0RunIICone
   SteeringPars["JetAlgo"] = ftable->TestParameterInSteering("JetAlgo");
   static int jetalgo;
   if ( SteeringPars["JetAlgo"] ) {
      ftable->GetParameterFromSteering("JetAlgo",jetalgo);
      if ( jetalgo < 0 || (2 < jetalgo && jetalgo < 10) || 12 < jetalgo ) {
         say::error["InclusiveDijetEvents"] << "Unknown jet algorithm " << jetalgo << ", aborted!" << endl;
         exit(1);
      }
   } else {
      say::error["InclusiveDijetEvents"] << "No jet algorithm selected, aborted!" << endl;
      exit(1);
   }
   SteeringPars["Rjet"] = ftable->TestParameterInSteering("Rjet");
   static double jetsize;
   if ( SteeringPars["Rjet"] ) {
      ftable->GetParameterFromSteering("Rjet",jetsize);
   } else {
      say::error["InclusiveDijetEvents"] << "Jet size R not defined, aborted!" << endl;
      exit(1);
   }
   SteeringPars["OvThr"] = ftable->TestParameterInSteering("OvThr");
   static double overlapthreshold = 0.5;
   if ( SteeringPars["OvThr"] ) {
      ftable->GetParameterFromSteering("OvThr",overlapthreshold);
   } else if ( jetalgo > 9 ) {
      say::error["InclusiveDijetEvents"] << "Overlap threshold not defined for jet algorithm " << jetalgo << ", aborted!" << endl;
      exit(1);
   }
   // --- fastNLO user: declare and initialize overall jet phase space cuts via steering file
   // overall lowest pT for jets to be considered
   SteeringPars["ptjmin"] = ftable->TestParameterInSteering("ptjmin");
   static double ptjmin;
   if ( SteeringPars["ptjmin"] ) {
      ftable->GetParameterFromSteering("ptjmin",ptjmin);
   } else {
      say::error["InclusiveDijetEvents"] << "Minimal jet pT (ptjmin) not defined, aborted!" << endl;
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
      say::error["InclusiveDijetEvents"] << "Minimal jet (pseudo)rapidity (yjmin or etajmin) not uniquely defined, aborted!" << endl;
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
      say::error["InclusiveDijetEvents"] << "Maximal jet (pseudo)rapidity (yjmax or etajmax) not uniquely defined, aborted!" << endl;
      exit(1);
   }
   // define logical for decision on cuts in (pseudo-)rapidity, no mixing allowed here
   static bool lpseudo;
   if ( SteeringPars["yjmin"] && SteeringPars["yjmax"] ) {
      lpseudo = false;
   } else if ( SteeringPars["etajmin"] && SteeringPars["etajmax"] ) {
      lpseudo = true;
   } else {
      say::error["InclusiveDijetEvents"] << "Phase space cuts mixed in (pseudo-)rapidity, aborted!" << endl;
      say::error["InclusiveDijetEvents"] << "Booleans for cut selections are" <<
         " yjmin "    << SteeringPars["yjmin"] <<
         ", yjmax "   << SteeringPars["yjmax"] <<
         ", etajmin " << SteeringPars["etajmin"] <<
         ", etajmax " << SteeringPars["etajmax"] << endl;
      say::error["InclusiveDijetEvents"] << "If you really want to mix, the code needs to be adapted." << endl;
      exit(1);
   }
   // minimal number of jets required to be within preselected jet phase space (for dijets this must be two!)
   SteeringPars["Njetmin"] = ftable->TestParameterInSteering("Njetmin");
   static int Njetmin = 2;
   if ( SteeringPars["Njetmin"] ) {
      ftable->GetParameterFromSteering("Njetmin",Njetmin);
   }
   if ( Njetmin < 2 ) {
      say::error["InclusiveDijetEvents"] << "This is a dijet scenario. At least two jets must be required, aborted!" << endl;
      say::error["InclusiveDijetEvents"] << "Please correct the Njetmin requirement. Njetmin = " << Njetmin << endl;
      exit(1);
   }

   // --- fastNLO user: declare and initialize dijet phase space cuts and definitions via steering file
   // overall minimum for observable one, e.g. maximal absolute rapidity |y_max|
   SteeringPars["obs0min"] = ftable->TestParameterInSteering("obs0min");
   static double obs0min = ftable->GetLoBinMin(0); // by default derived from binning in obs0
   if ( SteeringPars["obs0min"] ) {
      ftable->GetParameterFromSteering("obs0min",obs0min);
   }
   // overall maximum for observable one, e.g. maximal absolute rapidity |y_max|
   SteeringPars["obs0max"] = ftable->TestParameterInSteering("obs0max");
   static double obs0max = ftable->GetUpBinMax(0); // by default derived from binning in obs0
   if ( SteeringPars["obs0max"] ) {
      ftable->GetParameterFromSteering("obs0max",obs0max);
   }
   // overall minimum for observable two, e.g. dijet mass mjj
   SteeringPars["obs1min"] = ftable->TestParameterInSteering("obs1min");
   static double obs1min = ftable->GetLoBinMin(1); // by default derived from binning in obs1
   if ( SteeringPars["obs1min"] ) {
      ftable->GetParameterFromSteering("obs1min",obs1min);
   }
   // overall maximum for observable two, e.g. dijet mass mjj
   SteeringPars["obs1max"] = ftable->TestParameterInSteering("obs1max");
   static double obs1max = ftable->GetUpBinMax(1); // by default derived from binning in obs1
   if ( SteeringPars["obs1max"] ) {
      ftable->GetParameterFromSteering("obs1max",obs1max);
   }
   // extra minimal pT requirement for leading jet
   SteeringPars["ptj1min"] = ftable->TestParameterInSteering("ptj1min");
   static double ptj1min = ptjmin; // default is overall jet pT cut
   if ( SteeringPars["ptj1min"] ) {
      ftable->GetParameterFromSteering("ptj1min",ptj1min);
   }
   // extra minimal pT requirement for 2nd leading jet
   SteeringPars["ptj2min"] = ftable->TestParameterInSteering("ptj2min");
   static double ptj2min = ptjmin; // default is overall jet pT cut
   if ( SteeringPars["ptj2min"] ) {
      ftable->GetParameterFromSteering("ptj2min",ptj2min);
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
      say::info["InclusiveDijetEvents"] << "This event from NLOJet++ has no jets with pT > 1 GeV. Skipped!" << endl;
      return;
   } else if (nj > 4) {
      say::error["InclusiveDijetEvents"] << "This event from NLOJet++ has more than four jets, which should never happen. Aborted!" << endl;
      exit(1);
   }

   // --- give some debug output before selection and sorting
   if ( say::debug.GetSpeak() ) {
      for (unsigned int i=1; i<=nj; i++) {
         double pti  = pj[i].perp();
         double yi   = pj[i].rapidity();
         double etai = pj[i].prapidity();
         say::debug["InclusiveDijetEvents"] << "before cuts: jet # i, pt, y, eta: " << i << ", " << pti << ", " << yi << ", " << etai << endl;
      }
   }

   // --- select jets in y (lpseudo = false) or eta (lpseudo = true) and ptjmin
   //     note: failing jets are not deleted but moved to the end of the jet array pj!
   fNLOSelector SelJets (yetajmin,yetajmax,ptjmin,lpseudo);

   // --- count number of selected jets left at this stage
   size_t njet = std::remove_if(pj.begin(), pj.end(), SelJets) - pj.begin();
   if ( (int)njet < Njetmin ) return; // Nothing to be done

   // --- sort selected n jets at beginning of jet array pj, by default decreasing in pt
   static fNLOSorter SortJets;
   std::sort(pj.begin(), pj.begin() + njet, SortJets);

   // --- give some debug output after selection and sorting
   if ( say::debug.GetSpeak() ) {
      say::debug["InclusiveDijetEvents"] << "# jets before and after phase space cuts: nj, njet = " << nj << ", " << njet << endl;
      if ( ! lpseudo ) {
         say::debug["InclusiveDijetEvents"] << "phase space cuts: yjmin, yjmax, ptjmin: " << yetajmin << ", " << yetajmax << ", " << ptjmin << endl;
      } else {
         say::debug["InclusiveDijetEvents"] << "phase space cuts: etajmin, etajmax, ptjmin: " << yetajmin << ", " << yetajmax << ", " << ptjmin << endl;
      }
      for (unsigned int i=1; i<=njet; i++) {
         double pti  = pj[i].perp();
         double yi   = pj[i].rapidity();
         double etai = pj[i].prapidity();
         say::debug["InclusiveDijetEvents"] << "after cuts: jet # i, pt, y, eta: " << i << ", " << pti << ", " << yi << ", " << etai << endl;
      }
   }

   // ---- fastNLO v2.2
   // Analyze inclusive dijet event

   // --- calculate first requested observable
   // rapidities of two leading jets
   double y1 = pj[1].rapidity();
   double y2 = pj[2].rapidity();
   double ystar = abs(y1-y2)/2.;
   double obs0;
   switch(obs0def) {
   case YMAX :
      // maximal rapidity
      obs0 = max(abs(y1),abs(y2));
      break;
   case YSTAR :
      // rapidity separation half
      obs0 = ystar;
      break;
   default :
      say::error["InclusiveDijetEvents"] << "Observable not yet implemented, aborted!" << endl;
      say::error["InclusiveDijetEvents"] << "DimLabel[0] = " << DimLabel[0] << endl;
      say::error["InclusiveDijetEvents"] << "Please complement this scenario to include the requested observable." << endl;
      exit(1);
   }

   // --- calculate second requested observable
   // dijet mass
   lorentzvector<double> pj12 = pj[1] + pj[2];
   double mjj = pj12.mag();
   if (mjj < 0.) {say::warn["InclusiveDijetEvents"] << "Negative mass encountered: " << mjj << endl;}
   double obs1;
   switch(obs1def) {
   case MJJGEV :
      // dijet mass in GeV
      obs1 = mjj;
      break;
   case MJJTEV :
      // dijet mass in TeV
      obs1 = mjj/1000.;
      break;
   default :
      say::error["InclusiveDijetEvents"] << "Observable not yet implemented, aborted!" << endl;
      say::error["InclusiveDijetEvents"] << "DimLabel[1] = " << DimLabel[1] << endl;
      say::error["InclusiveDijetEvents"] << "Please complement this scenario to include the requested observable." << endl;
      exit(1);
   }

   // average pTs of leading jets
   double pT1   = (pj[1].perp()) / 1.0;
   double pT12  = (pj[1].perp() + pj[2].perp()) / 2.0;
   double pT123 = (pj[1].perp() + pj[2].perp() + pj[3].perp()) / 3.0;

   // --- Further dijet phase space cuts?
   if ( obs0min <= obs0 && obs0 < obs0max &&
        obs1min <= obs1 && obs1 < obs1max &&
        ptj1min <= pT1  && ptj2min <= pj[2].perp() ) {

      // --- set the renormalization and factorization scales
      // --- calculate first requested scale
      double mu1;
      switch(mu1def) {
      case PTMAX :
         // maximal jet pT
         mu1 = pT1;
         break;
      case PT12AVE :
         // average pT of leading two jets
         mu1 = pT12;
         break;
      case PT123AVE :
         // average pT of leading three jets
         mu1 = pT123;
         break;
      case MJJHALF :
         // half of dijet mass (must be in GeV!)
         mu1 = mjj/2.;
         break;
      case PTMAXEXPYSTAR :
         // ATLAS definition
         mu1 = pT1 * exp(0.3*ystar);
         break;
      default :
         say::error["InclusiveDijetEvents"] << "Scale not yet implemented, aborted!" << endl;
         say::error["InclusiveDijetEvents"] << "ScaleLabel[0] = " << ScaleLabel[0] << endl;
         say::error["InclusiveDijetEvents"] << "Please complement this scenario to include the requested scale." << endl;
         exit(1);
      }

      // --- calculate second requested scale
      double mu2 = pT1; // default second choice (only for flexible tables)
      if ( lFlexibleScaleTable ) {
         switch(mu2def) {
         case PTMAX :
            // maximal jet pT
            mu2 = pT1;
            break;
         case PT12AVE :
            // average pT of leading two jets
            mu2 = pT12;
            break;
         case PT123AVE :
            // average pT of leading three jets
            mu2 = pT123;
            break;
         case MJJHALF :
            // half of dijet mass (must be in GeV!)
            mu2 = mjj/2.;
            break;
         case PTMAXEXPYSTAR :
            // ATLAS definition
            mu2 = pT1 * exp(0.3*ystar);
            break;
         default :
            say::warn["InclusiveDijetEvents"] << "Scale not yet implemented, aborted!" << endl;
            say::warn["InclusiveDijetEvents"] << "ScaleLabel[1] = " << ScaleLabel[1] << endl;
            say::warn["InclusiveDijetEvents"] << "Please complement this scenario to include the requested scale." << endl;
         }
      }
      static vector<double> scalevars;
      if ( ! lFlexibleScaleTable ) scalevars = ftable->GetScaleVariations();

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
      if ( NDim == 1 ) {        // 1D binning
         // scen.SetObservableDimI( obs0  , 0 );
         say::error["InclusiveDijetEvents"] << "So far only 2D binning implemented here, aborted!" << endl;
         exit(1);
      } else if ( NDim == 2 ) { // 2D binning
         scen.SetObservableDimI( obs0, 0 );
         scen.SetObservableDimI( obs1, 1 );
      } else if ( NDim == 3 ) { // 3D binning
         // scen.SetObservableDimI( obs0, 0 );
         // scen.SetObservableDimI( obs1, 1 );
         // scen.SetObservableDimI( obs2, 2 );
         say::error["InclusiveDijetEvents"] << "So far only 2D binning implemented here, aborted!" << endl;
         exit(1);
      } else {
         say::error["InclusiveDijetEvents"] << "More than 3D binning not implemented for inclusive jets, aborted!" << endl;
         say::error["InclusiveDijetEvents"] << "DifferentialDimension NDim = " << NDim << endl;
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
   nwritemax = 10000000;
}

void UserHHC::end_of_event(){
   say::debug["UserHHC::end_of_event"] << "---------- UserHHC::end_of_event called ----------" << endl;
   nevents += 1;
   // --- store table
   say::debug["UserHHC::end_of_event"] << " nevents = " << nevents << ", nwrite = " << nwrite << endl;
   if (( (unsigned long)nevents % nwrite)==0){
      ftable->SetNumberOfEvents(nevents);
      ftable->WriteTable();
      if ( nwrite < nwritemax ) {
         nwrite *= 10;
      } else {
         nwrite = nwritemax;
      }
   }
}

void UserHHC::phys_output(const std::basic_string<char>& __file_name,
                          unsigned long __save, bool __txt)
{
   say::debug["UserHHC::phys_output"] << "---------- UserHHC::phys_output called ----------" << endl;
   say::debug["UserHHC::phys_output"] << "Before: __save = " << __save << ", nwrite = " << nwrite << endl;
   nwrite = __save;
   InitFastNLO(__file_name);
}
