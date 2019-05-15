// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {


  class MCgrid_ATLAS_2012_I1082936 : public Analysis {
  public:

    /// Constructor
    MCgrid_ATLAS_2012_I1082936()
      : Analysis("MCgrid_ATLAS_2012_I1082936")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// @todo Initialise and register projections here

      /// @todo Book histograms here, e.g.:
      // _h_XXXX = bookProfile1D(1, 1, 1);
      // _h_YYYY = bookHisto1D(2, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// @todo Do the event by event analysis here

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here

      // scale(_h_YYYY, crossSection()/sumOfWeights()); // norm to cross section
      // normalize(_h_YYYY); // normalize to unity

    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


    /// @name Histograms
    //@{
    Profile1DPtr _h_XXXX;
    Histo1DPtr _h_YYYY;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MCgrid_ATLAS_2012_I1082936);


}
