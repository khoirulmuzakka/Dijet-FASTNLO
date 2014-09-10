#include "fnlo_int_nlojet/fj-ak.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <cstdio>
#include <iostream>
using namespace std;

const bounded_vector<lorentzvector<double> >&
fj_ak::operator()(const event_hhc& ev, double jetsize) {

   int np = ev.upper();

   //----- initialize the arrays -----
   _M_pj.resize(1,0);

   //----- declare boolean for some initial output
   static bool first = true;

   if ( first ) {
      cout << "**************************\n";
      cout << "Input objects: ip, px, py, pz, E, np\n";
      cout << "--------------------------\n";
   }
   vector<fastjet::PseudoJet> input_objects;
   for (int ip = 1; ip <= np; ip++) {
      const double px = ev[ip].X();
      const double py = ev[ip].Y();
      const double pz = ev[ip].Z();
      const double E  = ev[ip].T();
      input_objects.push_back(fastjet::PseudoJet(px,py,pz,E));
      if ( first ) {
         printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n", ip, px, py, pz, E, np);
      }
   }
   if ( first ) {
      cout << "**************************\n";
   }

   fastjet::Strategy strategy = fastjet::Best;
   fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jetsize, strategy);
   fastjet::ClusterSequence clust_seq(input_objects, jet_def);

   double ptmin = 1.0;
   vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

   if ( first ) {
      cout << "**************************\n";
      cout << "Jet algorithm: " << jet_def.description() << endl;
      cout << "**************************\n";
      cout << "Output jets (pT > 1 GeV): ij, px, py, pz, E, nj\n";
      cout << "--------------------------\n";
   }
   unsigned int nj = inclusive_jets.size();
   for (unsigned int ij = 1; ij <= nj; ij++) {
      vector<fastjet::PseudoJet> constituents = clust_seq.constituents(inclusive_jets[ij-1]);
      const double px = inclusive_jets[ij-1].px();
      const double py = inclusive_jets[ij-1].py();
      const double pz = inclusive_jets[ij-1].pz();
      const double E  = inclusive_jets[ij-1].E();
      _M_pj.push_back( _Lv(px, py, pz, E) );
      if ( first ) {
         printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n", ij, px, py, pz, E, nj );
      }
   }
   if ( first ) {
      cout << "**************************\n";
      //      first = false;
   }

   return _M_pj;
}
