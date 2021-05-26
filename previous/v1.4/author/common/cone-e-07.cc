#include "cone-e-07.h"
#include <cmath>


const bounded_vector<lorentzvector<double> >&
cone_e_07::operator()(const event_hhc& ev, double rcone)
{
  int merge = 0, nj = 0, np = 0, nt = ev.upper();
  double dist;

  //----- initialize the arrays -----
  _M_pj.resize(1,0); _M_p.resize(1, nt); 
  _M_ax.resize(1,1);

  for(int ip = 1; ip <= nt; ip++) {
    if (ev[ip].perp2() > 1.0e-12) _M_p[++np] = ev[ip];
    const double px = ev[ip].X();
    const double py = ev[ip].Y();
    const double pz = ev[ip].Z();
    const double E  = ev[ip].T();
//     cout << "**************************\n";
//     cout << "Input objects: ip, px, py, pz, E, np\n";
//     printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n",
//  	   ip, px, py, pz, E, nt);
//     cout << "**************************\n";
  }

  nj = np;
  for ( int i = 1; i <= np; i++ ) {
    for ( int j = i + 1; j <= np; j++ ) {     
      if ( merge == 0 && (dist = _M_pair(i,j) ) < 1.0 ) {
	_M_merge(i, j);               //   i<j      Run II E-scheme 
	merge = 1;
	nj--;
	if(j != np) _M_p[j] = _M_p[np];
      }
    }
  }

  // copy merged particles to jet array 
  for (int i = 1; i <=nj; i++) {
    _M_pj.push_back(_M_p[i]);
  }

//   nj++;
//   double Emax1 = 0.;
//   double Emax2 = 0.;
//   unsigned int iord[] = {0,0,0}; 
//   for (unsigned int j = 1; j < nj; j++) {
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
//   unsigned int jp = 0;
//   for (unsigned int j = 1; j < nj; j++) {
//     if ( iord[j-1] > 0 ) {
//       double px = _M_pj[iord[j-1]].X();
//       double py = _M_pj[iord[j-1]].Y();
//       double pz = _M_pj[iord[j-1]].Z();
//       double E  = _M_pj[iord[j-1]].T();
//       double pt = sqrt(px*px + py*py);
//       if ( pt >= 1. ) {
// 	jp++;
//  	cout << "**************************\n";
//  	cout << "Output jets pt >= 1.0: ijet, px, py, pz, E, pt\n";
//  	printf("%5u %15.8f %15.8f %15.8f %15.8f %15.8f\n",
//  	       jp, px, py, pz, E, pt);
//  	cout << "**************************\n";
//       }
//     }
//   }

  return _M_pj;
}

double cone_e_07::_M_pair(int i, int j)
{
  static const double pi = 3.14159265358979323846;
  static const double twopi = 6.28318530717958647692;
  static const double rcone = 0.7;
  static const double rsep = 2.0;
  
   //   - check if distance between both particles is < R_sep*R_cone)
  int k; 
  bool beammerge = false;
  double ei   = _M_p[i].T();
  double ej   = _M_p[j].T();
  double et2i = _M_p[i].perp2();
  double et2j = _M_p[j].perp2();
  double dist  = 0.;
  double distx = 99.9;
  if ( et2i/ei/ei < 1.e-24 && et2j/ej/ej < 1.e-24 ) {
    // merge collinear particles in beam direction if both +z or -z
    if ( _M_p[i].Z() * _M_p[j].Z() > 0. ) { 
      dist = 0.;
      beammerge = true;
    } else {
      dist = 1.e12;
    }
  } else if ( et2i/ei/ei < 1.e-24 || et2j/ej/ej < 1.e-24 ) {
    // do not merge if only one particle in beam direction
    dist = 1.e12;
  } else {
    double dy = _M_p[i].rapidity() - _M_p[j].rapidity();
    double dphi = _M_p[i].phi() - _M_p[j].phi();  
    
    if ( dphi >= pi ) dphi = fmod(pi+dphi, twopi) - pi;
    else if ( dphi < -pi ) dphi = -fmod(pi-dphi, twopi) + pi;
    
    dist = sqrt(dy*dy+dphi*dphi);
  }

  if ( dist < (rcone*rsep) ) {
    double pti = _M_p[i].perp(), ptj = _M_p[j].perp();
    if ( pti < ptj ) {
      k=i;
    } else {
      k=j;
    }
    
    // compute potential jet axis
    _M_ax[1] = _M_p[i] + _M_p[j];
    
    // if two collinear beam particles are merged, set distx to zero
    if ( beammerge ) {
      distx = 0.;
    // check only distance of lower pT particle to jet axis
    //           (higher pT particle always has smaller distance!) 
    } else {
      double dyx = _M_p[k].rapidity() - _M_ax[1].rapidity();
      double dphix = _M_p[k].phi() - _M_ax[1].phi();  
      
      if ( dphix >= pi ) dphix = fmod(pi+dphix, twopi) - pi;
      else if( dphix < -pi ) dphix = -fmod(pi-dphix, twopi) + pi;
      
      distx = sqrt(dyx*dyx+dphix*dphix);
    }
  }
  // test 
  // if (rx>r) cout<<" MW:xxxxxxxxxxxxxx r "<<r<<"  "<<rx<<endl;
  // if (rx/r < 0.5)    cout<<" MW: r/rx "<<r<<"  "<<rx<<"  "<<pti/ptj<<endl;
  
  return (distx/rcone);
}
