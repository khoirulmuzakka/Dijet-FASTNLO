#include "fnlo_int_nlojet/fastjet-jets.h"
#include "fastnlotk/speaker.h"
#include "fastjet/ClusterSequence.hh"

#include "fastjet/config.h"
#ifdef FASTJET_ENABLE_PLUGIN_CDFCONES
#include "fastjet/CDFMidPointPlugin.hh"
#endif
#ifdef FASTJET_ENABLE_PLUGIN_D0RUNIICONE
#include "fastjet/D0RunIIConePlugin.hh"
#endif
#ifdef FASTJET_ENABLE_PLUGIN_SISCONE
#include "fastjet/SISConePlugin.hh"
#endif

#include <cstdio>
#include <iostream>

using namespace std;
using namespace nlo;

fastjet_jets::fastjet_jets() : jet_def(NULL), plugin(NULL) {
  _M_pj.resize(1, 0); // this defines an index offset!
}

fastjet_jets::fastjet_jets(const JetAlgorithm jetalgo, const double jetsize,
                           const double overlapthreshold)
    : jet_def(NULL), plugin(NULL) {
  setup(jetalgo, jetsize, overlapthreshold);
}

fastjet_jets::~fastjet_jets() { reset(); }

void fastjet_jets::reset() {
  if (jet_def)
    delete jet_def;
  jet_def = NULL;

  if (plugin)
    delete plugin;
  plugin = NULL;
}

void fastjet_jets::setup(const JetAlgorithm jetalgo, const double jetsize,
                         const double overlapthreshold) {
  reset();
  //----- select fastjet jet algorithm -----
  switch (jetalgo) {
  case fastjet_jets::KT:
  case fastjet_jets::CA:
  case fastjet_jets::antiKT: {
    // fastjet clustering jet algos: 0 = kT, 1 = CA, 2 = anti-kT
    fastjet::Strategy strategy = fastjet::Best;
    fastjet::JetAlgorithm JetAlgo = static_cast<fastjet::JetAlgorithm>(jetalgo);
    jet_def = new fastjet::JetDefinition(JetAlgo, jetsize, strategy);
  } break;

  case fastjet_jets::SISCone:
#ifdef FASTJET_ENABLE_PLUGIN_SISCONE
    plugin = new fastjet::SISConePlugin(jetsize, overlapthreshold);
    jet_def = new fastjet::JetDefinition(plugin);
#endif
    break;

  case fastjet_jets::CDFMidPointCone:
#ifdef FASTJET_ENABLE_PLUGIN_CDFCONES
    plugin = new fastjet::CDFMidPointPlugin(jetsize, overlapthreshold);
    jet_def = new fastjet::JetDefinition(plugin);
#endif
    break;

  case fastjet_jets::D0RunIICone:
// allocate a new plugin for D0RunIICone
#ifdef FASTJET_ENABLE_PLUGIN_D0RUNIICONE
    double min_jet_Et = 6.0;
    plugin =
        new fastjet::D0RunIIConePlugin(jetsize, min_jet_Et, overlapthreshold);
    jet_def = new fastjet::JetDefinition(plugin);
#endif
    break;
  }
  if (jet_def == NULL) {
    cerr << "fj-jets.cc: Selected jet algorithm undefined: jetalgo = "
         << jetalgo << endl;
    cerr << "fj-jets.cc: Aborted!" << endl;
    exit(1);
  }
}

const bounded_vector<lorentzvector<double> > &fastjet_jets::
operator()(const event_hhc &ev) {

  //----- minimal pt to be kept in jet results -----
  double ptmin = 1.;

  //----- fill input objects from NLOJet++ event -----
  if (say::debug.GetSpeak()) {
    say::debug["fj-jets"] << "**************************\n";
    say::debug["fj-jets"] << "Input objects: ip, px, py, pz, E, np\n";
    say::debug["fj-jets"] << "--------------------------\n";
  }
  int np = ev.upper();
  input_objects.clear();
  input_objects.reserve(np);
  for (int ip = 1; ip <= np; ip++) {
    const double px = ev[ip].X();
    const double py = ev[ip].Y();
    const double pz = ev[ip].Z();
    const double E = ev[ip].T();
    input_objects.push_back(fastjet::PseudoJet(px, py, pz, E));
    if (say::debug.GetSpeak()) {
      printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n", ip, px, py, pz, E, np);
    }
  }
  say::debug["fj-jets"] << "**************************\n";

  //----- run the jet clustering with the above jet definition -----
  fastjet::ClusterSequence clust_seq(input_objects, *jet_def);

  //----- get the resulting jets imposing a minimum pt -----
  output_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

  //----- fill output jets from NLOJet++ event -----
  if (say::debug.GetSpeak()) {
    say::debug["fj-jets"] << "**************************\n";
    say::debug["fj-jets"] << "Jet algorithm: " << jet_def->description()
                          << endl;
    say::debug["fj-jets"] << "**************************\n";
    say::debug["fj-jets"]
        << "Output jets (pT > 1 GeV): ij, px, py, pz, E, nj\n";
    say::debug["fj-jets"] << "--------------------------\n";
  }
  unsigned int nj = output_jets.size();
  _M_pj.clear();
  _M_pj.reserve(nj);
  for (unsigned int ij = 0; ij < nj; ij++) {
    const double px = output_jets[ij].px();
    const double py = output_jets[ij].py();
    const double pz = output_jets[ij].pz();
    const double E = output_jets[ij].E();
    _M_pj.push_back(_Lv(px, py, pz, E));
    if (say::debug.GetSpeak()) {
      printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n", ij + 1, px, py, pz, E,
             nj);
    }
  }
  say::debug["fj-jets"] << "**************************\n";

  return _M_pj;
}
