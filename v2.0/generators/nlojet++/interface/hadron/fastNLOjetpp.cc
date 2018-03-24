#include "fnlo_int_nlojet/fastNLOjetpp.h"
#include "fastnlotk/fastNLOCreate.h"
#include <iostream>

#include "fnlo_int_nlojet/fnlo_int_hhc_nlojet.h"

//----- declaration of the user defined functions -----
// --- fastNLO v2.2: interface to NLOJet++: read steering file, set LO of
// selected process
void inputfunc(unsigned int &, unsigned int &, unsigned int &);
// --- fastNLO v2.2: interface to NLOJet++: set cms energy and phase space
// generator
void psinput(nlo::phasespace_hhc *, double &);
// --- fastNLO v2.2: interface to NLOJet++: user class
nlo::user_base_hhc *userfunc();

//----- array of the symbols symbols -----
extern "C" {
   struct nlo_symbol {
      const char *name;
      void *address;
   } user_defined_functions[] = {
      //   process index: hhc for hadron-hadron --> jets
      {"procindex", (void *)"hhc"},
      //   input function
      {"inputfunc", (void *)inputfunc},
      //   phase space input function
      {"psinput", (void *)psinput},
      //   user defined functions
      {"userfunc", (void *)userfunc},
      //  end of the list
      {0, 0}};
}

// --- fastNLO v2.2: interface to NLOJet++: read steering file, set LO of
// selected process
void inputfunc(unsigned int &nj, unsigned int &nu, unsigned int &nd) {
   // This is the very first routine called from NLOJet++!
   // It seems that it gets called three times.
   static int ILOord; // no default
   if (!FastNLOUserHHC::instance->ftable) {
      // The following lines should be printed ONCE. Steering parameters have not
      // yet been read.
      std::cout << " # INIT:  [inputfunc] ---------- inputfunc called ----------"
                << std::endl;
      std::cout << " # INIT:  [inputfunc] ---------- initializing ... ----------"
                << std::endl;
      // --- read in steering and create fastNLO table accordingly
      // --- ftable is a global constant
      FastNLOUserHHC::instance->ftable =
         new fastNLOCreate(UsefulNlojetTools::GenConsts(),
                           UsefulNlojetTools::ProcConsts_HHC(),
                           "InclusiveNJets.str");
      if (FastNLOUserHHC::instance->ftable->TestParameterInSteering(
                                                                    "LeadingOrder")) {
         FastNLOUserHHC::instance->ftable->GetParameterFromSteering("LeadingOrder",
                                                                    ILOord);
      } else {
         say::error["ScenarioCode"] << "LO of process not defined, aborted!"
                                    << std::endl;
         exit(1);
      }
      //      ftable->SetLoOrder(ILOord); // Not necessary when set via
      //      LeadingOrder in steering file
      std::cout << " # INIT:  [inputfunc] ---------- LeadingOrder = " << ILOord
                << " ----------" << std::endl;
   }

   // The following line is printed ONLY, when DEBUG print out has been requested
   // by the steering.
   say::debug["inputfunc"] << "---------- inputfunc called ----------" << std::endl;

   // Set the number of jets for the LO process according to the steering,
   // e.g. 2 for hh inclusive jets, 3 for hh 3-jet mass
   // nj = 1U; // only useful for DIS
   switch (ILOord) {
   case 2:
      // hh inclusive jets, hh dijets
      nj = 2U;
      break;
   case 3:
      // hh 3-jets
      nj = 3U;
      break;
   default:
      say::error["ScenarioCode"] << "Unknown LO of process defined, aborted!"
                                 << std::endl;
      exit(1);
   }

   // Set the number of the (massless!) up and down type flavours (usually, you
   // won´t change that)
   nu = 2U;
   nd = 3U;
}

// --- fastNLO v2.2: interface to NLOJet++: set cms energy and phase space
// generator
void psinput(nlo::phasespace_hhc *ps, double &s) {
   //   std::cout << " # INIT:  [psinput] ---------- psinput called ----------" <<
   //   std::endl;
   say::debug["psinput"] << "---------- psinput called ----------" << std::endl;

   // --- set the center-of-mass energy squared as read from steering file
   s = pow(FastNLOUserHHC::instance->ftable->GetEcms(), 2);
   say::debug["psinput"] << "cms energy read from steering: sqrt(s) = "
                         << FastNLOUserHHC::instance->ftable->GetEcms() << std::endl;

   // --- in principle alternative phase space generators can be used
   // --- we support only the default for now
   ps = 0;
}

// --- fastNLO v2.2: interface to NLOJet++: user class
nlo::user_base_hhc *userfunc() {
   say::debug["userfunc"] << "---------- userfunc called ----------" << std::endl;
   return FastNLOUserHHC::instance;
}

FastNLOUserHHC::FastNLOUserHHC() : ftable(NULL) {}

// --- fastNLO user: evaluate steering file and define physics output (called
// once before first event)
void FastNLOUserHHC::phys_output(const std::basic_string<char> &fname,
                                 unsigned long nsave, bool txt) {
   //   std::cout << " # INIT:  [UserHHC::phys_output] ---------- UserHHC::phys_output
   //   called ----------" << std::endl;
   say::debug["UserHHC::phys_output"]
      << "---------- UserHHC::phys_output called ----------" << std::endl;
   say::debug["UserHHC::phys_output"] << "Before: __save = " << nsave
                                      << ", nwrite = " << nwrite << std::endl;
   nwrite = nsave;

   say::debug["InitfNLO"] << "---------- InitfNLO called ----------" << std::endl;
   // --- obtain relevant variables from NLOJet++ command line arguments
   ftable->SetOrderOfAlphasOfCalculation(
                                         UsefulNlojetTools::GetOrderOfRun(fname));

   // --- set fastNLO filename according to NLOJet++ command line arguments
   ftable->SetFilename(fname + ".tab");

   read_steering();
};

// --- fastNLO v2.2: initialize event counter and storage limit (called once)
void FastNLOUserHHC::initfunc(unsigned int) {
   say::debug["UserHHC::initfunc"]
      << "---------- UserHHC::initfunc called ----------" << std::endl;
   nevents = 0;
   nwritemax = 10000000;
};

// --- fastNLO v2.2: count events and store table (called after each event)
void FastNLOUserHHC::end_of_event() {
   say::debug["UserHHC::end_of_event"]
      << "---------- UserHHC::end_of_event called ----------" << std::endl;
   // --- count events
   nevents += 1;
   // --- store table
   say::debug["UserHHC::end_of_event"] << " nevents = " << nevents
                                       << ", nwrite = " << nwrite << std::endl;
   if (((unsigned long)nevents % nwrite) == 0) {
      ftable->SetNumberOfEvents(nevents);
      ftable->WriteTable();
      if (nwrite < nwritemax) {
         nwrite *= 10;
      } else {
         nwrite = nwritemax;
      }
   }
};
