#include "fj-ak.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <cmath>
using namespace std;

const bounded_vector<lorentzvector<double> >&
fj_ak::operator()(const event_hhc& ev, double jetsize)
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

  fastjet::Strategy strategy = fastjet::Best;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jetsize, strategy);
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
// 	     jp, px, py, pz, E, pt);
//       cout << "**************************\n";
//     }
//   }

  return _M_pj;
}
