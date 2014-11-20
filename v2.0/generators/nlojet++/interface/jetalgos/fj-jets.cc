#include "fnlo_int_nlojet/fj-jets.h"
#include "fastnlotk/speaker.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/CDFMidPointPlugin.hh"
#include "fastjet/D0RunIIConePlugin.hh"
#include "fastjet/SISConePlugin.hh"
#include <cstdio>
#include <iostream>
using namespace std;

const bounded_vector<lorentzvector<double> >&
fj_jets::operator()(const event_hhc& ev, int jetalgo, double jetsize, double overlapthreshold) {

   //----- initialize output array -----
   _M_pj.resize(1,0);

   //----- initialize fastjet pseudojets -----
   vector<fastjet::PseudoJet> input_objects;
   vector<fastjet::PseudoJet> output_jets;

   //----- create null pointers to yet unfilled JetDefinition class and plugin -----
   fastjet::JetDefinition* jet_def = NULL;
   fastjet::JetDefinition::Plugin* plugin = NULL;

   //----- minimal pt to be kept in jet results -----
   double ptmin = 1.;

   //----- fill input objects from NLOJet++ event -----
   if ( say::debug.GetSpeak() ) {
      say::debug["fj-jets"] << "**************************\n";
      say::debug["fj-jets"] << "Input objects: ip, px, py, pz, E, np\n";
      say::debug["fj-jets"] << "--------------------------\n";
   }
   int np = ev.upper();
   for (int ip = 1; ip <= np; ip++) {
      const double px = ev[ip].X();
      const double py = ev[ip].Y();
      const double pz = ev[ip].Z();
      const double E  = ev[ip].T();
      input_objects.push_back(fastjet::PseudoJet(px,py,pz,E));
      if ( say::debug.GetSpeak() ) {
         printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n", ip, px, py, pz, E, np);
      }
   }
   say::debug["fj-jets"] << "**************************\n";

   //----- select fastjet jet algorithm -----
   if ( -1 < jetalgo && jetalgo < 3 ) {
      // fastjet clustering jet algos: 0 = kT, 1 = CA, 2 = anti-kT
      fastjet::Strategy strategy = fastjet::Best;
      fastjet::JetAlgorithm JetAlgo = static_cast<fastjet::JetAlgorithm>(jetalgo);
      jet_def = new fastjet::JetDefinition(JetAlgo, jetsize, strategy);
   } else if ( 9 < jetalgo && jetalgo < 13 ) {
      // fastjet cone jet algos: 10 = SISCone, 11 = CDFMidPointCone, 12 = D0RunIICone
      if ( jetalgo == 10 ) {
         // allocate a new plugin for SISCone
         plugin = new fastjet::SISConePlugin(jetsize, overlapthreshold);
      } else if ( jetalgo == 11 ) {
         // allocate a new plugin for CDFMidPointCone
         plugin = new fastjet::CDFMidPointPlugin(jetsize, overlapthreshold);
      } else if ( jetalgo == 12 ) {
         // allocate a new plugin for D0RunIICone
         double min_jet_Et = 6.0;
         plugin = new fastjet::D0RunIIConePlugin(jetsize, min_jet_Et, overlapthreshold);
      } else {
         cerr << "fj-jets.cc: Selected jet algorithm undefined: jetalgo = " << jetalgo << endl;
         cerr << "fj-jets.cc: Aborted!" << endl;
         exit(1);
      }
      // create a jet-definition based on the plugin
      jet_def = new fastjet::JetDefinition(plugin);
   }

   //----- run the jet clustering with the above jet definition -----
   fastjet::ClusterSequence clust_seq(input_objects, *jet_def);

   //----- get the resulting jets imposing a minimum pt -----
   output_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

   //----- fill output jets from NLOJet++ event -----
   if ( say::debug.GetSpeak() ) {
      say::debug["fj-jets"] << "**************************\n";
      say::debug["fj-jets"] << "Jet algorithm: " << jet_def->description() << endl;
      say::debug["fj-jets"] << "**************************\n";
      say::debug["fj-jets"] << "Output jets (pT > 1 GeV): ij, px, py, pz, E, nj\n";
      say::debug["fj-jets"] << "--------------------------\n";
   }
   unsigned int nj = output_jets.size();
   for (unsigned int ij = 1; ij <= nj; ij++) {
      vector<fastjet::PseudoJet> constituents = clust_seq.constituents(output_jets[ij-1]);
      const double px = output_jets[ij-1].px();
      const double py = output_jets[ij-1].py();
      const double pz = output_jets[ij-1].pz();
      const double E  = output_jets[ij-1].E();
      _M_pj.push_back( _Lv(px, py, pz, E) );
      if ( say::debug.GetSpeak() ) {
         printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n", ij, px, py, pz, E, nj );
      }
   }
   say::debug["fj-jets"] << "**************************\n";

   // todo: it is much more efficient to construct plugin and jet definition only once per run
   if (jet_def)
      delete jet_def;
   if (plugin)
      delete plugin;

   return _M_pj;
}
