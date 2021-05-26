// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

// MCgrid headers
#include "mcgrid/mcgrid.hh"
#include "mcgrid/mcgrid_binned.hh"

// How many of the 2nd dimension bins should be translated to grids
//#define N_HISTOS 9
#define N_HISTOS 1

namespace Rivet {

   // Dijet mass histograms, R=0.6
   class MCgrid_ATLAS_2012_I1082936 : public Analysis {
   public:

      /// Constructor
      MCgrid_ATLAS_2012_I1082936() : Analysis("MCgrid_ATLAS_2012_I1082936") {}


      /// Book histograms and initialise projections before the run
      void init() {

#if USE_APPL && USE_FNLO
         cerr << "MCgrid error: This analysis is configured to generate both an APPLgrid and a fastNLO table." << endl;
         cerr << "But using multiple binned histograms is not supported yet." << endl;
         exit(-1);
#endif

         // Initialize the projectors:
         const FinalState fs;
         addProjection(fs,"FinalState");

         FastJets fj06(fs,  FastJets::ANTIKT, 0.6);
         fj06.useInvisibles();
         addProjection(fj06, "AntiKT06");

         // Histogram binning in 2nd dimension
         double ystarbins[] = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.4};

         size_t massDsOffset = 4; // 4th series of histograms, i.e. d04-xnn-ymm
         for (size_t i = 0; i < N_HISTOS; i++) {
            _hist_mass.addHistogram(ystarbins[i], ystarbins[i+1], bookHisto1D(massDsOffset, 1, i+1));
         }

         const vector<Histo1DPtr> histos = _hist_mass.getHistograms();
         string subprocFileName("basic");
         // string subprocFileName("MCgrid_ATLAS_2012_I1082936");
#if USE_APPL
         subprocFileName += ".config";
#elif USE_FNLO
         subprocFileName += ".str";
#endif
         MCgrid::subprocessConfig subproc(subprocFileName,
                                          MCgrid::BEAM_PROTON,
                                          MCgrid::BEAM_PROTON);


#if USE_APPL
         MCgrid::applGridConfig config(2, subproc, MCgrid::highPrecAPPLgridArch, 1E-5, 1, 10, 1E7);
#elif USE_FNLO
         MCgrid::fastnloGridArch arch(15, 6, "Lagrange", "Lagrange", "sqrtlog10", "loglog025");
         MCgrid::fastnloConfig config(2, subproc, arch, 7000.0);
#endif

         // Book grids, and push into BinnedGrid instance
         MCgrid::gridPtr grid;
         for (size_t i(0); i < N_HISTOS; i++) {
            grid = MCgrid::bookGrid(histos[i], histoDir(), config);
            _grid_mass.addGrid(ystarbins[i], ystarbins[i+1], grid);
         }
      }

      /// Perform the per-event analysis
      void analyze(const Event& event) {

         // MCgrid event counter
         MCgrid::PDFHandler::HandleEvent(event, histoDir());

         const double weight = event.weight();
         const FastJets &fj  = applyProjection<FastJets>(event, "AntiKT06");
         const Jets& jets    = fj.jetsByPt(Cuts::pT > 20*GeV && Cuts::absrap < 4.4);

         // Identify dijets
         vector<FourMomentum> leadjets;
         foreach (const Jet& jet, jets) {
            const double pT = jet.pT();

            // Make sure we have a leading jet with pT >30 GeV and a second to leading jet with pT>20 GeV
            if (leadjets.size() < 2) {
               if (leadjets.empty() && pT < 30*GeV) continue;
               leadjets.push_back(jet.momentum());
            }
         }
         // Make sure we have the required two leading jets
         if (leadjets.size() < 2) {
            MSG_DEBUG("Could not find two suitable leading jets");
            //            cout << "MSG_DEBUG: Could not find two suitable leading jets" << endl;
         } else {
            //            cout << "pT1 = " << leadjets[0].pT() << ", pT2 = " << leadjets[1].pT() << ", dpT12 = " << leadjets[0].pT() - leadjets[1].pT() << endl;
            const double y1 = leadjets[0].rapidity();
            const double y2 = leadjets[1].rapidity();
            const double ystar = fabs(y1-y2)/2.;
            const double m = (leadjets[0] + leadjets[1]).mass();
            //            cout << "mass = " << m << endl;
            // Fill mass histogram
            _hist_mass.fill(ystar, m/TeV, weight);
            _grid_mass.fill(ystar, m/TeV, event);
         }
      }

      /// Normalise histograms etc., after the run
      void finalize() {
         //         _hist_mass.scale(crossSectionPerEvent()/picobarn, this);
         //         _grid_mass.scale(crossSectionPerEvent()/picobarn);
         _hist_mass.scale(crossSection()/sumOfWeights(), this);
         _grid_mass.scale(crossSection()/sumOfWeights());
         _grid_mass.exportgrids();

         cout << "DEBUG: crossSection() = " << crossSection() << ", sumOfWeights() = " << sumOfWeights() << ", picobarn = " << picobarn << endl;
         cout << "DEBUG: crossSectionPerEvent()        = " << crossSectionPerEvent() << endl;
         cout << "DEBUG: crossSection()/sumOfWeights() = " << crossSection()/sumOfWeights() << endl;

         // Clear MCgrid event counter
         MCgrid::PDFHandler::CheckOutAnalysis(histoDir());
      }

   private:

      /// The di-jet mass spectrum binned in rapidity for akt6 jets
      BinnedHistogram<double>    _hist_mass;
      MCgrid::BinnedGrid<double> _grid_mass;
   };

   // The hook for the plugin system
   DECLARE_RIVET_PLUGIN(MCgrid_ATLAS_2012_I1082936);

}
