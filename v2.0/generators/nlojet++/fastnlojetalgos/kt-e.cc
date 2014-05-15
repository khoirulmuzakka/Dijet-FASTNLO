#include "fastnloja/kt-e.h"
#include <cmath>
using namespace std;

const bounded_vector<lorentzvector<double> >&
kt_e::operator()(const event_hhc& ev, double jetsize)
{
  int imin, jmin, kmin, nt = ev.upper();
  double tmp, pmin, smin, r2 = jetsize*jetsize;

  //----- initialize the arrays -----
  _M_pj.resize(1,0); _M_p.resize(1, nt);

  for(int ip = 1; ip <= nt; ip++) {
    _M_p[ip] = ev[ip];
    const double px = ev[ip].X();
    const double py = ev[ip].Y();
    const double pz = ev[ip].Z();
    const double E  = ev[ip].T();
//     cout << "**************************\n";
//     cout << "Input objects:\n";
//     printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n",
//         ip, px, py, pz, E, nt);
//     cout << "**************************\n";
  }

  unsigned int nj = 1;
  for(int n = nt; n > 1; n--) {
    //----- find the smallest resolution variables -----
    imin = 1; jmin = 2; kmin = 1;
    pmin = _M_pair(1,2); smin = _M_p[1].perp2();
    //    cout << "in pmin " << pmin << "in smin " << smin << endl;
    for(int i = 1; i <= n; i++) {
       if((tmp = _M_p[i].perp2()) < smin) {
        smin = tmp;
        kmin = i;
      }

      for(int j = i + 1; j <= n; j++)
        if((tmp = _M_pair(i, j)) < pmin) {
          pmin = tmp;
          imin = i;
          jmin = j;
        }
    }
    //    cout << "fin pmin " << pmin << " fin smin " << smin << endl;
    //    cout << "fin kmin " << kmin << " fin imin " << imin << " fin jmin " << jmin << endl;
    //    cout << "r2 " << r2 << " smin*r2 " << smin*r2 << " pmin " << pmin << endl;
    //----- do the clustering -----
    if(smin*r2 < pmin) {
      _M_pj.push_back(_M_p[kmin]);
      nj++;
      if(kmin != n) _M_p[kmin] = _M_p[n];
    } else {
      _M_merge(imin, jmin);
      if(jmin != n) _M_p[jmin] = _M_p[n];
    }
  }

  //----- add the last (pseudo)particle to the list of jets -----
  _M_pj.push_back(_M_p[1]);
  nj++;

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
//      jp++;
//      cout << "**************************\n";
//      cout << "Output jets pt >= 1.0: ijet, px, py, pz, E, pt\n";
//      printf("%5u %15.8f %15.8f %15.8f %15.8f %15.8f\n",
//          jp, px, py, pz, E, pt);
//      cout << "**************************\n";
//       }
//     }
//   }

  return _M_pj;
}

double kt_e::_M_pair(int i, int j)
{
  static const double pi = 3.14159265358979323846;
  static const double twopi = 6.28318530717958647692;

  //  double et2 = (_M_p[i].perp2() < _M_p[j].perp2() ? _M_p[i].perp2() : _M_p[j].perp2());
  double ei   = _M_p[i].T();
  double ej   = _M_p[j].T();
  double et2i = _M_p[i].perp2();
  double et2j = _M_p[j].perp2();
  double dist = 0.;
  if ( et2i/ei/ei < 1.e-24 && et2j/ej/ej < 1.e-24 ) {
    // merge collinear particles in beam direction if both +z or -z
    if ( _M_p[i].Z() * _M_p[j].Z() > 0. ) {
      dist = 0.;
    } else {
      dist = 1.e12;
    }
  } else if ( et2i/ei/ei < 1.e-24 || et2j/ej/ej < 1.e-24 ) {
    // do not merge if only one particle in beam direction
    dist = 1.e12;
  } else {
    double et2 = ( et2i < et2j ? et2i : et2j );
    double dy = _M_p[i].rapidity() - _M_p[j].rapidity();
    double dphi = _M_p[i].phi() - _M_p[j].phi();

    if ( dphi >= pi ) dphi = fmod(pi+dphi, twopi) - pi;
    else if( dphi < -pi) dphi = -fmod(pi-dphi, twopi) + pi;

    dist = et2*(dy*dy+dphi*dphi);
    //    cout << "et2: " << et2 << " dy: " << dy << " dphi: " << dphi << endl;
  }

  //  return et2*(dy*dy+dphi*dphi);
  return dist;
}

void kt_e::_M_merge(int i, int j)
{
   _M_p[i] += _M_p[j];
}
