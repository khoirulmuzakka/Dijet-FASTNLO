#include "fj-mc-05.h"
#include "fastjet/CDFMidPointPlugin.hh"
// Well, in SISCone the following header is already included at some place
#include "fastjet/ClusterSequence.hh"
#include <cmath>
using namespace std;

const bounded_vector<lorentzvector<double> >&
fj_mc_05::operator()(const event_hhc& ev, double rwhatfor)
{
  int nt = ev.upper();

  //----- initialize the arrays -----
  _M_pj.resize(1,0); _M_p.resize(1, nt); 
  _M_ax.resize(1,1);

  vector<fastjet::PseudoJet> input_objects;
  for (int ip = 1; ip <= nt; ip++) {
    const double px = ev[ip].X();
    const double py = ev[ip].Y();
    const double pz = ev[ip].Z();
    const double E  = ev[ip].T();
//     cout << "**************************\n";
//     cout << "Input objects: ip, px, py, pz, E, np\n";
//     printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n",
//  	   ip, px, py, pz, E, nt);
//     cout << "**************************\n";
    input_objects.push_back(fastjet::PseudoJet(px,py,pz,E));
  }

  // define a generic plugin pointer
  fastjet::JetDefinition::Plugin * plugin = 0;
  
  // allocate a new plugin for CDFMidPoint
  double seed_threshold = 1.0;
  double cone_radius = 0.5;
  double cone_area_fraction = 1.0;
  int    max_pair_size = 2;
  int    max_iterations = 100;
  double overlap_threshold = 0.75;
  plugin = new fastjet::CDFMidPointPlugin (seed_threshold, cone_radius,
                                           cone_area_fraction, max_pair_size,
                                           max_iterations, overlap_threshold);
  
  // create a jet-definition based on the plugin
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
  fastjet::JetDefinition jet_def(plugin);
  fastjet::ClusterSequence clust_seq(input_objects, jet_def);

  double ptmin = 1.0;
  vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);
  
  unsigned int nj = inclusive_jets.size();
  for (unsigned int i = 1; i <= nj; i++) {
    vector<fastjet::PseudoJet> constituents = clust_seq.constituents(inclusive_jets[i-1]);
    const double px = inclusive_jets[i-1].px();
    const double py = inclusive_jets[i-1].py();
    const double pz = inclusive_jets[i-1].pz();
    const double E  = inclusive_jets[i-1].E();
    //    printf("%5u %15.8f %15.8f %15.8f %15.8f\n",
    //    	   i, px, py, pz, E );
    _M_pj.push_back( _Lv(px, py, pz, E) );
  }

//   double Emax1 = 0.;
//   double Emax2 = 0.;
//   unsigned int iord[] = {0,0,0}; 
//   for (unsigned int j = 1; j <= nj; j++) {
//     double E  = _M_pj[j].T();
//     if ( E > Emax1 ) {
//       Emax2 = Emax1;
//       Emax1 = E;
//       iord[2] = iord[1];
//       iord[1] = iord[0];
//       iord[0] = j;
//     } else if ( E > Emax2 ) {
//       Emax2 = E;
//       iord[2] = iord[1];
//       iord[1] = j;
//     } else {
//       iord[2] = j;
//     }
//   }
//   for (unsigned int j = 1; j <= nj; j++) {
//    if ( iord[j-1] > 0 ) {
//       double px = _M_pj[iord[j-1]].X();
//       double py = _M_pj[iord[j-1]].Y();
//       double pz = _M_pj[iord[j-1]].Z();
//       double E  = _M_pj[iord[j-1]].T();
//       double pt = sqrt(px*px + py*py);
//       cout << "**************************\n";
//       cout << "Output jets pt >= 1.0: ijet, px, py, pz, E, pt\n";
//       printf("%5u %15.8f %15.8f %15.8f %15.8f %15.8f\n",
// 	     j, px, py, pz, E, pt);
//       cout << "**************************\n";
//     }
//   }

  delete plugin;
  return _M_pj;
}
