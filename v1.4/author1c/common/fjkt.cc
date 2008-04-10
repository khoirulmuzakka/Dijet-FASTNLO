#include "fjkt.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <cmath>
using namespace std;

const bounded_vector<lorentzvector<double> >&
fjkt::operator()(const event_hhc& ev, double rcone)
{
  int merge = 0, nj = 0, np = 0, nt = ev.upper();
  double dist;

  //----- initialize the arrays -----
  _M_pj.resize(1,0); _M_p.resize(1, nt); 
  _M_ax.resize(1,1);

  vector<fastjet::PseudoJet> input_objects;
  for (int ip = 1; ip <= nt; ip++) {
    const double px = ev[ip].X();
    const double py = ev[ip].Y();
    const double pz = ev[ip].Z();
    const double E  = ev[ip].T();
    cout << "**************************\n";
    cout << "Input objects:\n";
    printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n",
	   ip, px, py, pz, E, nt);
    cout << "**************************\n";
    input_objects.push_back(fastjet::PseudoJet(px,py,pz,E));
  }
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::JetDefinition jet_def(fastjet::kt_algorithm, rcone, strategy);
  fastjet::ClusterSequence clust_seq(input_objects, jet_def);
  
  double ptmin = 0.0;
  vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);
  
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    vector<fastjet::PseudoJet> constituents = clust_seq.constituents(inclusive_jets[i]);

    //    printf("%5u %15.8f %15.8f %15.8f %8u\n",
    //	   i, inclusive_jets[i].rap(),
    //   inclusive_jets[i].phi(),
    //   inclusive_jets[i].perp(),
    //   constituents.size());
    
    const double px = inclusive_jets[i].px();
    const double py = inclusive_jets[i].py();
    const double pz = inclusive_jets[i].pz();
    const double E  = inclusive_jets[i].E();
    cout << "**************************\n";
    cout << "Output jets:\n";
    printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n",
	   i, px, py, pz, E, inclusive_jets.size());
    cout << "**************************\n";
    
    _M_pj.push_back( _Lv() );
    _M_pj[i].X() = px;
    _M_pj[i].Y() = py;
    _M_pj[i].Z() = pz;
    _M_pj[i].T() = E;
  }

  return _M_pj;
}
