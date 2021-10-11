//!
//! UsefulNlojetTools
//!
//! Collection of useful functions for the fastNLO interface
//! to NLOJET++.
//!
// KR: Drop differentiation in 2jet and 3jet ProcConsts.
//     Default is 2jet, set differences for 3jet via steering!

#include <cstdio>
#include <algorithm> //c++98
#include <utility>   //c++11
#include "fastnlotk/fastNLOEvent.h"
#include "fastnlotk/fastNLOTableConstants.h"
#include "fastnlotk/speaker.h"
#include "fnlo_int_nlojet/pdf-hhc-dummy.h"
#include <bits/hhc-phasespace.h>
#include <bits/hhc-jetfunc.h>

using namespace std;
using namespace nlo;

//----- declaration of the user defined functions -----
// --- fastNLO v2.2: interface to NLOJet++: read steering file, set LO of
// selected process
void inputfunc(unsigned int &, unsigned int &, unsigned int &);
// --- fastNLO v2.2: interface to NLOJet++: set cms energy and phase space
// generator
void psinput(nlo::phasespace_hhc *, double &);
// --- fastNLO v2.2: interface to NLOJet++: user class
nlo::user_base_hhc *userfunc();

namespace UsefulNlojetTools {
   /**
      namespace UsefulNlojetTools

      Collection of useful functions and constant for the interface
      between nlojet++ and fastNLO, if nlojet++ is run in hhc-mode
      (i.e. for pp and ppbar collisions).
   */

   //_______________________________________________________________________
   fastNLO::GeneratorConstants GenConsts() {
      fastNLO::GeneratorConstants GenConsts;
      GenConsts.Name = "NLOJet++_4.1.3";
      GenConsts.UnitsOfCoefficients = 12;
      return GenConsts;
   }

   //_______________________________________________________________________
   fastNLO::ProcessConstants ProcConsts_HHC() {
      fastNLO::ProcessConstants ProcConsts;
      ProcConsts.Name = "";
      ProcConsts.References.push_back(
                                      "Z. Nagy, Phys. Rev. Lett. 88 (2002) 122003,");
      ProcConsts.References.push_back("Z. Nagy, Phys. Rev. D 68 (2003) 094002.");
      // Default: 2-jet observable, i.e. LO power of alpha_S is 2; set to 3 in
      // steering for 3-jet observables!
      ProcConsts.LeadingOrder = 2;
      // Default: 12 ,i.e. pb; if published differently set in steering explicitly.
      //          Note: Internally, NLOJet++ always uses pb.
      ProcConsts.NPDF = 2;
      // KR: Merging q(qbar)g and gq(qbar) subprocesses 5 and 6 into one is
      //     ONLY correct at the statistical limit of large numbers of such events!
      //     Effectively, one generated event is reinterpreted as half of a q(qbar)g
      //     and half
      //     of a gq(qbar) event such that direct comparisons with reference
      //     calculations
      //     are NOT possible any more event-by-event, but only statistically with
      //     large numbers.
      //     On the other hand, the gain in table size is only of the order of -15%.
      //     ==> Therefore:
      //     New default: NSubProcessesLO always 7!
      //     Old default: NSubProcessesLO = 6 at LO for 2-jet observable, 7
      //     otherwise
      // Attention: GetNSubproc(), see below, has to changed accordingly!
      // TODO: Search for more elegant solution here ...
      ProcConsts.NSubProcessesLO = 7;
      ProcConsts.NSubProcessesNLO = 7;
      ProcConsts.NSubProcessesNNLO = 7;
      ProcConsts.IPDFdef1 = 3;
      ProcConsts.IPDFdef2 = 1;
      ProcConsts.IPDFdef3LO = 1;
      ProcConsts.IPDFdef3NLO = 2;
      ProcConsts.IPDFdef3NNLO = 2;
      ProcConsts.NPDFDim = 2;
      // To test full-matrix storage uncomment the following line or set NPDFDim in
      // steering explicitly!
      //      ProcConsts.NPDFDim = 2;
      ProcConsts.AsymmetricProcesses.push_back(make_pair(5, 6));
      ProcConsts.AsymmetricProcesses.push_back(make_pair(6, 5));
      //
      //
      return ProcConsts;
   }

   //_______________________________________________________________________
   unsigned int GetNj() {
      // get order of leading-order
      // inputfunc should be defined in header of nlojet-module beforehand.
      unsigned int nj = 0, nu = 0, nd = 0;
      ::inputfunc(nj, nu, nd);
      return nj;
   }

   double GetEcms() {
      // get center of mass energy
      double ecms = 0;
      ::psinput(NULL, ecms);
      return sqrt(ecms);
   }

   //_______________________________________________________________________
   unsigned int GetLoOrder() { return GetNj(); }

   //_______________________________________________________________________
   unsigned int GetOrderOfRun(const basic_string<char> &__file_name) {
      // --- determine whether we are running LO or NLO
      // and add order of expansion to leading-order .
      static int ord = -1;

      if (ord == -1) {
         const char *const file = __file_name.c_str();
         if (strstr(file, "born") != NULL) {
            // IsNLO = false;
            say::info["GetOrderOfRun"] << "This is a LO run." << endl;
            ord = GetNj();
         } else if (strstr(file, "nlo") != NULL) {
            // IsNLO = true;
            say::info["GetOrderOfRun"] << "This is a NLO run." << endl;
            ord = GetNj() + 1;
         } else if (strlen(file) == 0) {
            // Order is wrongly assumed to be initialized already;
            say::error["GetOrderOfRun"] << "Order of run not properly initialized."
                                        << endl;
            exit(1);
         } else {
            // No alternatives implemented with NLOJet++ version 4
            say::error["GetOrderOfRun"]
               << "This module can only be run at Born level or at NLO!" << endl;
            exit(1);
         }
         say::info["GetOrderOfRun"]
            << "Found order of calculation to be: ord = " << ord << endl;
         return ord;
      } else {
         // Return previously initialized value
         return ord;
      }
   }

   //_______________________________________________________________________
   int NlojetToFastnloIDHHC(int id) {
      if (id == 1)
         return 5;
      else if (id == 2)
         return 6;
      else if (id == 3)
         return 1;
      else if (id == 4)
         return 2;
      else if (id == 5)
         return 3;
      else if (id == 6)
         return 4;
      else
         return id;
   }
   int FastnloIdToNlojetIdHHC(int id) {
      if (id == 1)
         return 3;
      else if (id == 2)
         return 4;
      else if (id == 3)
         return 5;
      else if (id == 4)
         return 6;
      else if (id == 5)
         return 1;
      else if (id == 6)
         return 2;
      else
         return id;
   }

   //_______________________________________________________________________
   vector<fnloEvent> GetFlexibleScaleNlojetContribHHC(const event_hhc &p,
                                                      const amplitude_hhc &amp) {
      static const unsigned int nSubproc = 7; // nlojet must have 7 subprocesses.
      static const double dummyMu2 = 91. * 91.;
      static pdf_hhc_dummy dummypdf;

      // make events
      vector<fnloEvent> ev(nSubproc);

      // fill relevant event quantities
      double x1 = p[-1].Z() / p[hadron(-1)].Z();
      double x2 = p[0].Z() / p[hadron(0)].Z();
      for (unsigned int p = 0; p < nSubproc; p++) {
         ev[p].SetProcessId(p);
         ev[p].SetX1(x1);
         ev[p].SetX2(x2);
      }

      static const double coef = 389385730.;
      // weights
      double weights[7][7]; // weights[amp_i][proc]
      for (int kk = 0; kk < 7; kk++)
         for (unsigned int p = 0; p < nSubproc; p++)
            weights[kk][p] = 0;

      // todo: calculate wt and wtorg from 'weights'
      nlo::weight_hhc wtorg = amp(dummypdf, dummyMu2, dummyMu2, 1.);
      nlo::weight_hhc wt = wtorg;
      // - rearrange subprocesses
      for (int fid = 0; fid < 7; fid++) {
         int nid = FastnloIdToNlojetIdHHC(fid);
         wt[fid] = wtorg[nid];
      }
      wt *= coef;

      // access perturbative coefficients
      nlo::amplitude_hhc::contrib_type itype = amp.contrib();
      nlo::weight_hhc cPDF = dummypdf.pdf(x1, x2, dummyMu2, 2, 3); // 1/x1/x2

      for (int kk = 0; kk < 7; kk++) {
         for (unsigned int fid = 0; fid < nSubproc; fid++) {
            int nid = FastnloIdToNlojetIdHHC(fid);
            weights[kk][fid] = amp._M_fini.amp[kk][nid] * coef * cPDF[nid];
         }
      }

      for (unsigned int p = 0; p < nSubproc; p++) {
         // decompose nlojet-event
         if (itype == nlo::amplitude_hhc::fini) {
            if (amp._M_fini.mode == 0) { // finix1
               ev[p].SetWeight_MuIndependent(weights[0][p]);
               ev[p].SetWeight_log_muf(weights[3][p]);
            } else if (amp._M_fini.mode == 1) { // finix2
               ev[p].SetWeight_MuIndependent(weights[1][p]);
               ev[p].SetWeight_log_muf(weights[4][p]);
            } else if (amp._M_fini.mode == 2) { // fini1
               ev[p].SetWeight_MuIndependent(weights[2][p]);
               ev[p].SetWeight_log_muf(weights[5][p]);
               ev[p].SetWeight_log_mur(weights[6][p]);
            }
         } else { // no fini contribution
            ev[p].AddWeight_MuIndependent(wt[p]);
         }
      }

      return ev;
   }

   // KR: Is it possible to get NSubproc from the Toolkit instead?
   //     This would allow to read the setting as defined in the steering in
   //     contrast to this duplication ...
   //     In addition, this method had set the # of subprocesses for LO 3-jet to 6
   //     which was never tried before ...
   //     ==> Fixed back to 7 now.
   //     New default: NSubProcesses always 7!
   //_______________________________________________________________________
   unsigned int GetNSubproc() {
      static int nSubproc = 7;
      // static int nSubproc = -1;
      // if ( nSubproc == -1 ) {
      //    // Order should have already been initialized --> redemand with empty
      //    string
      //    // If order not yet determined --> stop with error message in
      //    GetOrderOfRun
      //    int nord  = GetOrderOfRun("");
      //    int loord = GetLoOrder();
      //    if ( loord == 2 && nord == loord ) nSubproc = 6 ;
      //    else nSubproc = 7;
      // }
      return nSubproc;
   }

   //_______________________________________________________________________
   vector<fnloEvent> GetFixedScaleNlojetContribHHC(const event_hhc &p,
                                                   const amplitude_hhc &amp,
                                                   double mu) {
      const int nSubproc = GetNSubproc();
      const double Mu2 = mu * mu;
      static pdf_hhc_dummy dummypdf;

      // make events
      vector<fnloEvent> ev(nSubproc);

      // fill relevant event quantities
      double x1 = p[-1].Z() / p[hadron(-1)].Z();
      double x2 = p[0].Z() / p[hadron(0)].Z();
      for (int p = 0; p < nSubproc; p++) {
         ev[p].SetProcessId(p);
         ev[p].SetX1(x1);
         ev[p].SetX2(x2);
      }

      // weights
      // access perturbative coefficients
      static const double coef = 389385730.;

      // todo: calculate wt and wtorg from 'weights'
      nlo::weight_hhc wtorg = amp(dummypdf, Mu2, Mu2, 1.);
      nlo::weight_hhc wt = wtorg;
      // - rearrange subprocesses
      for (int fid = 0; fid < 7; fid++) {
         int nid = FastnloIdToNlojetIdHHC(fid);
         // KR: Check on isnan|isinf now also for 'faster' fill code in toolkit
         //     Commented out here.
         //         if (std::isfinite(wtorg[nid])) {
         wt[fid] = wtorg[nid];
         // } else {
         //    //            wt[fid] = 0.;
         //    wt[fid] = wtorg[nid];
         //    if (std::isnan(wtorg[nid])) {
         //       say::error["GetFixedScaleNlojetContribHHC"]<<"Weight wtorg for nid = " << nid << " is 'nan'!"<<endl;
         //    } else if (std::isinf(wtorg[nid])) {
         //       say::error["GetFixedScaleNlojetContribHHC"]<<"Weight wtorg for nid = " << nid << " is 'inf'!"<<endl;
         //    } else {
         //       say::error["GetFixedScaleNlojetContribHHC"]<<"Weight wtorg for nid = " << nid << " is non-finite!"<<endl;
         //    }
         // }
      }
      wt *= coef;

      // KR, 23.01.2015: This has been moved to the Toolkit part and will disappear
      // here
      //      if(x2>x1){
      // swap subprocesses 6,7
      //         swap(wt[5],wt[6]);
      //      }

      // KR, 23.01.2015: This is part of the NLOJet++ process definitions (for
      // fixed-scale tables) and remains in the interface
      if (nSubproc == 6) {
         wt[5] = (wt[5] + wt[6]) / 2;
         wt[6] = wt[5];
      }

      for (int p = 0; p < nSubproc; p++) {
         ev[p].AddWeight_MuIndependent(wt[p]);
      }
      return ev;
   }

   //_______________________________________________________________________
   vector<vector<fnloEvent> >
   GetFixedScaleNlojetContribHHC(const event_hhc &p, const amplitude_hhc &amp,
                                 double mu, const vector<double> &scalevar) {
      // return contributions for each scale variation
      // return value can directly passed to fastNLOCreate
      vector<vector<fnloEvent> > ctrbs(scalevar.size());
      for (unsigned int is = 0; is < scalevar.size(); is++) {
         double smu = mu * scalevar[is];
         ctrbs[is] = GetFixedScaleNlojetContribHHC(p, amp, smu);
      }
      return ctrbs;
   }
}
