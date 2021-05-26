#ifndef __fastjet_jets_h__
#define __fastjet_jets_h__ 1

#include <bits/hhc-event.h>
#include <bits/hep-bounded_vector.h>
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"

class fastjet_jets {

  //   private types
  typedef nlo::lorentzvector<double> _Lv;

public:
  enum JetAlgorithm {
    // fastjet base algorithms
    KT = fastjet::kt_algorithm,
    CA = fastjet::cambridge_algorithm,
    antiKT = fastjet::antikt_algorithm,
    // jet algorithms included via plugins
    SISCone = 10,
    CDFMidPointCone = 11,
    D0RunIICone = 12
  };

  fastjet_jets();
  fastjet_jets(const JetAlgorithm jetalgo, const double jetsize,
               const double overlapthreshold);
  ~fastjet_jets();
  //   do the clustering and return with the momenta of the jets
  const nlo::bounded_vector<_Lv> &operator()(const nlo::event_hhc &);
  void setup(const JetAlgorithm jetalgo, const double jetsize,
             const double overlapthreshold);

private:
  //   private data members
  nlo::bounded_vector<_Lv> _M_pj;
  void reset();

  fastjet::JetDefinition *jet_def;
  fastjet::JetDefinition::Plugin *plugin;
  std::vector<fastjet::PseudoJet> input_objects;
  std::vector<fastjet::PseudoJet> output_jets;
};

#endif
