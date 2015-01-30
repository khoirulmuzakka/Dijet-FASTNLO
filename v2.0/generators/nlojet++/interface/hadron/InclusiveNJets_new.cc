//
// fastNLO v2.2 creator code for inclusive N jets scenarios
//
// ============== fastNLO user: ========================================
// To create your own scenario, it is recommended to take this
// example code, make a copy, and implement the relevant changes.
// Important:
// Edit only those parts which are clearly indicated as modifiable
// "fastNLO user" code. If a code fragment is not explicitely
// labeled as "fastNLO user", it is likely that a modification will
// interfere with the fastNLO routines.
//
// This file contains the following routines:
//   struct fNLOSelector   (-> user edits)
//   struct fNLOSorter     (-> user edits)
//   UserHHC::phys_output  (called once at the start   -> user edits)
//   UserHHC::userfunc     (called once for each event -> user edits)
//   inputfunc             (don't touch)
//   psinput               (don't touch)
//   userfunc              (don't touch)
//   InitfNLO              (don't touch)
//   UserHHC::initfunc     (don't touch)
//   UserHHC::end_of_event (don't touch)
//
// If the provided example routine can not be steered flexibly enough,
// feel free to implement the missing parts.
//
// Implementing a new scenario may imply to change:
//  - the jet algorithm ("#include" statement and assignment of "jetclus..")
//  - the jet observables to compute
//  - the scale definition
//
// =====================================================================

//------ DON'T TOUCH THIS PART! ------
#include <cfloat>
#include <iostream>
#include <map>
#include <bits/hhc-phasespace.h>
#include <bits/hhc-process.h>
#include <bits/hhc-jetfunc.h>

//----- used namespaces -----
using namespace nlo;
using namespace std;

//----- fastNLO -----
#include "fastnlotk/fastNLOEvent.h"
#include "fastnlotk/fastNLOCreate.h"
// --- fastNLO v2.2: include header file of interface to NLOJet++
#include "fastNLOjetpp.h"
#include "fnlo_int_nlojet/fnlo_int_hhc_nlojet.h"

//------ USER DEFINED PART STARTS HERE ------

// --- fastNLO user: include header file for the jet algorithm
#include "fnlo_int_nlojet/fastjet-jets.h"

// --- fastNLO v2.2: define user class to be used with NLOJet++
class UserHHC : public FastNLOUserHHC {
public:
    // --- fastNLO user: evaluate steering file and define physics output (called once before first event)
    virtual void read_steering();
    // --- fastNLO user: analyze parton event (called once for each event)
    virtual void userfunc(const event_hhc&, const amplitude_hhc&);

private:
    // --- fastNLO user: define the jet algorithm (for the choice of included header file above)
    fastjet_jets jetclusfj;

    // --- define the jet structure
    bounded_vector<lorentzvector<double> > pj;

    // --- fastNLO steering
    bool lFlexibleScaleTable; // Fill fixed- or flexible-scale table
    int NDim; // Dimensionality of distributions
    vector<string> DimLabel; // Dimension labels
    // enum to switch between implemented observables (max. of 3 simultaneously)
    enum Obs { PTJETGEV,
               YJET,
               PHIJET };
    Obs obsdef[3];
    double obs[3];
    vector<string> ScaleLabel; // Scale labels
    // enum to switch between implemented scale definitions (max. of 2 simultaneously)
    // (Njet > 1!)
    enum Scales { PTMAX,
                  PTJET };
    Scales mudef[2];
    double mu[2];
    int jetalgo; // Define fastjet jet algorithm
    double jetsize; // Define jet size R
    double overlapthreshold; // Define overlap threshold (for some jet algorithms)
    double ptjmin; // Minimal jet pT (should be >= minimum of 1 GeV specified in interface to fastjet)
    double yetajmin; // Minimal jet (pseudo-)rapidity
    double yetajmax; // Maximal jet (pseudo-)rapidity
    bool lpseudo; // Switch to use either jet rapidity y or jet eta
    int Njetmin; // Minimal number of jets in phase space
    double obsmin[3]; // Minimum in observable in nth dimension (default derived from binning)
    double obsmax[3]; // Maximum in observable in nth dimension (default derived from binning)
};

// --- fastNLO user: modify the jet selection in UserHHC::userfunc (default = cutting in |y| min, |y| max and pt min)
//                   (the return value must be true for jets to be UNselected)
struct fNLOSelector {
    fNLOSelector(double ymin, double ymax, double ptmin, bool pseudo = false)
        : _ymin(ymin)
        , _ymax(ymax)
        , _ptmin(ptmin)
        , _pseudo(pseudo){};
    double _ymin, _ymax, _ptmin;
    bool _pseudo;
    bool operator()(const lorentzvector<double>& a)
    {
        if (!_pseudo)
            return !(_ymin <= abs(a.rapidity()) && abs(a.rapidity()) < _ymax && _ptmin <= a.perp());
        else
            return !(_ymin <= abs(a.prapidity()) && abs(a.prapidity()) < _ymax && _ptmin <= a.perp());
    };
};

FastNLOUserHHC* FastNLOUserHHC::instance = new UserHHC;

// --- fastNLO user: modify the jet sorting in UserHHC::userfunc (default = descending in jet pt)
struct fNLOSorter {
    bool operator()(const lorentzvector<double>& a, const lorentzvector<double>& b) { return (a.perp() > b.perp()); };
};

// --- fastNLO user: check and get steering parameters once and store into static vars
static std::map<std::string, bool> SteeringPars;

// --- fastNLO user: class UserHHC: evaluate steering file and define physics output (called once before first event)
void UserHHC::read_steering()
{
    // --- fastNLO user:
    //     Here is your playground where you can evaluate the steering file settings.
    //     ATTENTION: Some settings are mandatory for the correct functioning!

    // get general steering parameters needed here from steering file
    say::debug["UserHHC::phys_output"] << "Evaluating steering parameters ..." << endl;
    // fixed- or flexible-scale table
    SteeringPars["FlexibleScaleTable"] = ftable->TestParameterInSteering("FlexibleScaleTable");
    lFlexibleScaleTable = false; // default
    if (SteeringPars["FlexibleScaleTable"]) {
        ftable->GetParameterFromSteering("FlexibleScaleTable", lFlexibleScaleTable);
    }
    // dimensionality
    SteeringPars["DifferentialDimension"] = ftable->TestParameterInSteering("DifferentialDimension");
    if (SteeringPars["DifferentialDimension"]) {
        ftable->GetParameterFromSteering("DifferentialDimension", NDim);
    } else {
        say::error["ScenarioCode"] << "Dimensioning of binning not set, aborted!" << endl;
        exit(1);
    }
    if (NDim < 1 || 3 < NDim) {
        say::error["ScenarioCode"] << "Only 1- to 3-dimensional binning implemented, aborted!" << endl;
        say::error["ScenarioCode"] << "Please implement the requested " << NDim << "-dimensional binning." << endl;
        exit(1);
    }
    // dimension labels
    SteeringPars["DimensionLabels"] = ftable->TestParameterInSteering("DimensionLabels");
    DimLabel.resize(NDim);
    if (SteeringPars["DimensionLabels"]) {
        ftable->GetParameterFromSteering("DimensionLabels", DimLabel);
    } else {
        say::error["ScenarioCode"] << "Dimension labels not set, aborted!" << endl;
        exit(1);
    }
    // define the observables according to the dimension labels
    for (int i = 0; i < NDim; i++) {
        if (DimLabel[i] == "|y|") {
            obsdef[i] = YJET;
        } else if (DimLabel[i] == "pT_[GeV]") {
            obsdef[i] = PTJETGEV;
        } else if (DimLabel[i] == "phi") {
            obsdef[i] = PHIJET;
        } else {
            say::error["ScenarioCode"] << "Unknown observable, i.e. dimension label, aborted!" << endl;
            say::error["ScenarioCode"] << "DimLabel[" << i << "] = " << DimLabel[i] << endl;
            say::error["ScenarioCode"] << "Please complement this scenario to include the requested observable." << endl;
            exit(1);
        }
    }
    // scale descriptions
    SteeringPars["ScaleDescriptionScale1"] = ftable->TestParameterInSteering("ScaleDescriptionScale1");
    ScaleLabel.resize(2);
    if (SteeringPars["ScaleDescriptionScale1"]) {
        ftable->GetParameterFromSteering("ScaleDescriptionScale1", ScaleLabel[0]);
    } else {
        say::error["ScenarioCode"] << "No description of scale 1, aborted!" << endl;
        exit(1);
    }
    SteeringPars["ScaleDescriptionScale2"] = ftable->TestParameterInSteering("ScaleDescriptionScale2");
    ScaleLabel[1] = "pT_max_[GeV]"; // default
    if (SteeringPars["ScaleDescriptionScale2"]) {
        ftable->GetParameterFromSteering("ScaleDescriptionScale2", ScaleLabel[1]);
    } else {
        ScaleLabel[1] = "-";
        say::warn["ScenarioCode"] << "No description of scale 2, flexible-scale tables not possible!" << endl;
    }
    // scale descriptions define the scales
    for (unsigned int i = 0; i < 2; i++) {
        if (ScaleLabel[i] == "pT_max_[GeV]") {
            mudef[i] = PTMAX;
        } else if (ScaleLabel[i] == "pT_jet_[GeV]") {
            mudef[i] = PTJET;
        } else {
            say::error["ScenarioCode"] << "Unknown scale, i.e. scale description, aborted!" << endl;
            say::error["ScenarioCode"] << "ScaleLabel[" << i << "] = " << ScaleLabel[i] << endl;
            say::error["ScenarioCode"] << "Please complement this scenario to include the requested scale." << endl;
            exit(1);
        }
    }

    // definition of jet algorithm and jet phase space limits (no defaults)
    //
    // --- fastNLO user: set the jet algorithm and size via steering file
    // fastjet clustering jet algos: 0 = kT, 1 = CA, 2 = anti-kT
    // fastjet cone jet algos: 10 = SISCone, 11 = CDFMidPointCone, 12 = D0RunIICone
    SteeringPars["JetAlgo"] = ftable->TestParameterInSteering("JetAlgo");
    if (SteeringPars["JetAlgo"]) {
        ftable->GetParameterFromSteering("JetAlgo", jetalgo);
        if (jetalgo < 0 || (2 < jetalgo && jetalgo < 10) || 12 < jetalgo) {
            say::error["ScenarioCode"] << "Unknown jet algorithm " << jetalgo << ", aborted!" << endl;
            exit(1);
        }
    } else {
        say::error["ScenarioCode"] << "No jet algorithm selected, aborted!" << endl;
        exit(1);
    }
    SteeringPars["Rjet"] = ftable->TestParameterInSteering("Rjet");
    if (SteeringPars["Rjet"]) {
        ftable->GetParameterFromSteering("Rjet", jetsize);
    } else {
        say::error["ScenarioCode"] << "Jet size R not defined, aborted!" << endl;
        exit(1);
    }
    SteeringPars["OvThr"] = ftable->TestParameterInSteering("OvThr");
    overlapthreshold = 0.5; // default
    if (SteeringPars["OvThr"]) {
        ftable->GetParameterFromSteering("OvThr", overlapthreshold);
    } else if (jetalgo > 9) {
        say::error["ScenarioCode"] << "Overlap threshold not defined for jet algorithm " << jetalgo << ", aborted!" << endl;
        exit(1);
    }
    // --- fastNLO user: declare and initialize overall jet phase space cuts via steering file
    // overall lowest pT for jets to be considered
    SteeringPars["ptjmin"] = ftable->TestParameterInSteering("ptjmin");
    if (SteeringPars["ptjmin"]) {
        ftable->GetParameterFromSteering("ptjmin", ptjmin);
    } else {
        say::error["ScenarioCode"] << "Minimal jet pT (ptjmin) not defined, aborted!" << endl;
        exit(1);
    }
    // overall highest pT for jets not implemented, since uncritical with respect to CPU time consumption
    // overall smallest |(pseudo-)rapidity| for jets to be considered, use either y or eta but not both
    SteeringPars["yjmin"] = ftable->TestParameterInSteering("yjmin");
    SteeringPars["etajmin"] = ftable->TestParameterInSteering("etajmin");
    if (SteeringPars["yjmin"] && !SteeringPars["etajmin"]) {
        ftable->GetParameterFromSteering("yjmin", yetajmin);
    } else if (!SteeringPars["yjmin"] && SteeringPars["etajmin"]) {
        ftable->GetParameterFromSteering("etajmin", yetajmin);
    } else {
        say::error["ScenarioCode"] << "Minimal jet (pseudo)rapidity (yjmin or etajmin) not uniquely defined, aborted!" << endl;
        exit(1);
    }
    // overall largest |(pseudo-)rapidity| for jets to be considered, use either y or eta but not both
    SteeringPars["yjmax"] = ftable->TestParameterInSteering("yjmax");
    SteeringPars["etajmax"] = ftable->TestParameterInSteering("etajmax");
    if (SteeringPars["yjmax"] && !SteeringPars["etajmax"]) {
        ftable->GetParameterFromSteering("yjmax", yetajmax);
    } else if (!SteeringPars["yjmax"] && SteeringPars["etajmax"]) {
        ftable->GetParameterFromSteering("etajmax", yetajmax);
    } else {
        say::error["ScenarioCode"] << "Maximal jet (pseudo)rapidity (yjmax or etajmax) not uniquely defined, aborted!" << endl;
        exit(1);
    }
    // define logical for decision on cuts in (pseudo-)rapidity, no mixing allowed here
    if (SteeringPars["yjmin"] && SteeringPars["yjmax"]) {
        lpseudo = false;
    } else if (SteeringPars["etajmin"] && SteeringPars["etajmax"]) {
        lpseudo = true;
    } else {
        say::error["ScenarioCode"] << "Phase space cuts mixed in (pseudo-)rapidity, aborted!" << endl;
        say::error["ScenarioCode"] << "Booleans for cut selections are"
                                   << " yjmin " << SteeringPars["yjmin"] << ", yjmax " << SteeringPars["yjmax"] << ", etajmin " << SteeringPars["etajmin"] << ", etajmax " << SteeringPars["etajmax"] << endl;
        say::error["ScenarioCode"] << "If you really want to mix, the code needs to be adapted." << endl;
        exit(1);
    }
    // minimal number of jets required to be within preselected jet phase space (for inclusive jets this must be one!)
    SteeringPars["Njetmin"] = ftable->TestParameterInSteering("Njetmin");
    Njetmin = 1;
    if (SteeringPars["Njetmin"]) {
        ftable->GetParameterFromSteering("Njetmin", Njetmin);
    }

    // --- fastNLO user: declare and initialize Njet phase space cuts and definitions via steering file
    // overall minimum & maximum for 1st observable, e.g. maximal absolute rapidity |y_max|
    SteeringPars["obs0min"] = ftable->TestParameterInSteering("obs0min");
    obsmin[0] = ftable->GetLoBinMin(0); // by default derived from binning in obs0
    if (SteeringPars["obs0min"]) {
        ftable->GetParameterFromSteering("obs0min", obsmin[0]);
    }
    SteeringPars["obs0max"] = ftable->TestParameterInSteering("obs0max");
    obsmax[0] = ftable->GetUpBinMax(0); // by default derived from binning in obs0
    if (SteeringPars["obs0max"]) {
        ftable->GetParameterFromSteering("obs0max", obsmax[0]);
    }
    // overall minimum & maximum for 2nd observable, e.g. dijet mass mjj
    obsmin[1] = -DBL_MAX;
    obsmax[1] = +DBL_MAX;
    if (NDim > 1) {
        SteeringPars["obs1min"] = ftable->TestParameterInSteering("obs1min");
        obsmin[1] = ftable->GetLoBinMin(1); // by default derived from binning in obs1
        if (SteeringPars["obs1min"]) {
            ftable->GetParameterFromSteering("obs1min", obsmin[1]);
        }
        SteeringPars["obs1max"] = ftable->TestParameterInSteering("obs1max");
        obsmax[1] = ftable->GetUpBinMax(1); // by default derived from binning in obs1
        if (SteeringPars["obs1max"]) {
            ftable->GetParameterFromSteering("obs1max", obsmax[1]);
        }
    }
    // overall minimum & maximum for 3rd observable
    obsmin[2] = -DBL_MAX;
    obsmax[2] = +DBL_MAX;
    if (NDim > 2) {
        SteeringPars["obs2min"] = ftable->TestParameterInSteering("obs2min");
        obsmin[2] = ftable->GetLoBinMin(2); // by default derived from binning in obs2
        if (SteeringPars["obs2min"]) {
            ftable->GetParameterFromSteering("obs2min", obsmin[2]);
        }
        SteeringPars["obs2max"] = ftable->TestParameterInSteering("obs2max");
        obsmax[2] = ftable->GetUpBinMax(2); // by default derived from binning in obs2
        if (SteeringPars["obs2max"]) {
            ftable->GetParameterFromSteering("obs2max", obsmax[2]);
        }
    }
    jetclusfj.setup(static_cast<fastjet_jets::JetAlgorithm>(jetalgo), jetsize, overlapthreshold);
}

// --- fastNLO v2.2: class UserHHC: analyze parton event (called once for each event)
void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
    say::debug["UserHHC::userfunc"] << "---------- UserHHC::userfunc called ----------" << endl;

    // --- fastNLO user:
    //     Here is your playground where you compute your observables and
    //     scales for each jet or event.
    //     The bin number ("obsbin") gets passed to fastNLO's table filling code.
    //     Usually, pT and E are in GeV, but this may be changed.
    //     ATTENTION: Scales must always be in GeV!

    // apply the jet algorithm to partonic 4-vector array p of NLOJet++
    pj = jetclusfj(p);
    unsigned int nj = pj.upper();

    // --- check on minimal and maximal no. of jets
    // ATTENTION: In principle, without cuts, there should always be two.
    //            For efficiency reasons though, our interface to the fastjet algorithms
    //            requires a minimal jet pT of 1 GeV. If this is a problem, the ptmin value
    //            in fj-jets.cc needs to be changed.
    // There should never be more than four jets in NLOJet++
    if (nj < 1) {
        say::debug["ScenarioCode"] << "This event from NLOJet++ has no jets with pT > 1 GeV. Skipped!" << endl;
        return;
    } else if (nj > 4) {
        say::error["ScenarioCode"] << "This event from NLOJet++ has more than four jets, which should never happen. Aborted!" << endl;
        exit(1);
    }

    // --- give some debug output before selection and sorting
    if (say::info.GetSpeak()) {
        for (unsigned int i = 1; i <= nj; i++) {
            double pti = pj[i].perp();
            double yi = pj[i].rapidity();
            double etai = pj[i].prapidity();
            say::info["ScenarioCode"] << "before cuts: jet # i, pt, y, eta: " << i << ", " << pti << ", " << yi << ", " << etai << endl;
        }
    }

    // --- select jets in y (lpseudo = false) or eta (lpseudo = true) and ptjmin
    //     note: failing jets are not deleted but moved to the end of the jet array pj!
    fNLOSelector SelJets(yetajmin, yetajmax, ptjmin, lpseudo);

    // --- count number of selected jets left at this stage
    size_t njet = std::remove_if(pj.begin(), pj.end(), SelJets) - pj.begin();
    if ((int)njet < Njetmin)
        return; // Nothing to be done

    // --- sort selected n jets at beginning of jet array pj, by default decreasing in pt
    static fNLOSorter SortJets;
    std::sort(pj.begin(), pj.begin() + njet, SortJets);

    // --- give some debug output after selection and sorting
    if (say::info.GetSpeak()) {
        say::info["ScenarioCode"] << "# jets before and after phase space cuts: nj, njet = " << nj << ", " << njet << endl;
        if (!lpseudo) {
            say::info["ScenarioCode"] << "phase space cuts: yjmin, yjmax, ptjmin: " << yetajmin << ", " << yetajmax << ", " << ptjmin << endl;
        } else {
            say::info["ScenarioCode"] << "phase space cuts: etajmin, etajmax, ptjmin: " << yetajmin << ", " << yetajmax << ", " << ptjmin << endl;
        }
        for (unsigned int i = 1; i <= njet; i++) {
            double pti = pj[i].perp();
            double yi = pj[i].rapidity();
            double etai = pj[i].prapidity();
            say::info["ScenarioCode"] << "after cuts: jet # i, pt, y, eta: " << i << ", " << pti << ", " << yi << ", " << etai << endl;
        }
    }

    // ---- fastNLO v2.2
    // Analyze inclusive jets in jet loop
    // set one possible scale choice (Attention: Only correct if jets sorted descending in pT)
    double ptmax = pj[1].perp();
    for (unsigned int k = 1; k <= njet; k++) {

        // derive some jet quantities
        double pt = pj[k].perp();
        double yeta;
        if (!lpseudo) {
            yeta = abs(pj[k].rapidity());
        } else {
            yeta = abs(pj[k].prapidity());
        }
        double phi = atan2(pj[k].Y(), pj[k].X());

        // --- calculate observable of nth dimension
        for (int i = 0; i < NDim; i++) {
            switch (obsdef[i]) {
            case PTJETGEV:
                // jet pT
                obs[i] = pt;
                break;
            case YJET:
                // jet rapidity
                obs[i] = yeta;
                break;
            case PHIJET:
                // jet azimuthal angle
                obs[i] = phi;
                break;
            default:
                say::error["ScenarioCode"] << "Observable not yet implemented, aborted!" << endl;
                say::error["ScenarioCode"] << "DimLabel[" << i << "] = " << DimLabel[i] << endl;
                say::error["ScenarioCode"] << "Please complement this scenario to include the requested observable." << endl;
                exit(1);
            }
        }

        // --- Further Njet phase space cuts?
        if (obsmin[0] <= obs[0] && obs[0] < obsmax[0] && (NDim < 2 || (obsmin[1] <= obs[1] && obs[1] < obsmax[1])) && (NDim < 3 || (obsmin[2] <= obs[2] && obs[2] < obsmax[2]))) {

            // --- set the renormalization and factorization scales
            // --- calculate the requested scales
            for (unsigned int i = 0; i < 2; i++) {
                switch (mudef[i]) {
                case PTMAX:
                    // maximal jet pT
                    mu[i] = ptmax;
                    break;
                case PTJET:
                    // jet pT
                    mu[i] = pt;
                    break;
                default:
                    say::error["ScenarioCode"] << "Scale not yet implemented, aborted!" << endl;
                    say::error["ScenarioCode"] << "ScaleLabel[" << i << "] = " << ScaleLabel[i] << endl;
                    say::error["ScenarioCode"] << "Please complement this scenario to include the requested scale." << endl;
                    exit(1);
                }
            }
            static vector<double> scalevars;
            if (!lFlexibleScaleTable)
                scalevars = ftable->GetScaleVariations();

            // get matrix elements
            static vector<fnloEvent> contribsflex;
            static vector<vector<fnloEvent> > contribsfix;
            if (lFlexibleScaleTable) {
                contribsflex = UsefulNlojetTools::GetFlexibleScaleNlojetContribHHC(p, amp);
            } else {
                contribsfix = UsefulNlojetTools::GetFixedScaleNlojetContribHHC(p, amp, mu[0], scalevars);
            }

            // set and fill scenario specific quantites
            fnloScenario scen;
            if (NDim < 1 || NDim > 3) {
                say::error["ScenarioCode"] << "Less than 1D(?!) or more than 3D binning not implemented here, aborted!" << endl;
                say::error["ScenarioCode"] << "DifferentialDimension NDim = " << NDim << endl;
            }
            for (int i = 0; i < NDim; i++) {
                scen.SetObservableDimI(obs[i], i);
            }

            scen.SetObsScale1(mu[0]); // must be consistent with 'mu' from contribs
            if (lFlexibleScaleTable) {
                scen.SetObsScale2(mu[1]);
                ftable->FillAllSubprocesses(contribsflex, scen);
            } else {
                ftable->FillAllSubprocesses(contribsfix, scen);
            }
        }
    }
}
