#include "fastnlo_interface_nlojet/cone-et.h"
#include <cmath>

const bounded_vector<lorentzvector<double> >&
cone_et::operator()(const event_hhc& ev, double jetsize)
{
  int merge = 0, nj = 0, np = 0, nt = ev.upper();
  double dist;

  //----- initialize the arrays -----
  _M_pj.resize(1,0); _M_p.resize(1, nt);
  for(int ip = 1; ip <= nt; ip++)
    if(ev[ip].perp2() > 1.0e-12) _M_p[++np] = ev[ip];

  nj = np;
  for(int i = 1; i <= np; i++) {
    for(int j = i + 1; j <= np; j++) {
      if(merge == 0 && (dist = _M_pair(i,j)) < jetsize) {
        _M_merge_Et(i, j);               //   i<j
        merge = 1;
        nj--;
        if(j != np) _M_p[j] = _M_p[np];
        //if (dist > 0.7) cout<<" MW: cone dist= "<<dist<<endl;
      }
    }
  }

  // copy merged particles to jet array
  for (int i = 1; i <=nj; i++) {
    _M_pj.push_back(_M_p[i]);
  }

  return _M_pj;
}

double cone_et::_M_pair(int i, int j)
{
  static const double pi = 3.14159265358979323846;
  static const double twopi = 6.28318530717958647692;
  static const double rsep = 2.0;

  double pti = _M_p[i].perp(), ptj = _M_p[j].perp();
  double pt  = (pti > ptj ? pti : ptj);
  double r0  = (pti+ptj)/pt;
  double r1  = (rsep < r0 ? rsep : r0);

  double deta = _M_p[i].prapidity() - _M_p[j].prapidity();
  double dphi = _M_p[i].phi() - _M_p[j].phi();
  if(dphi >= pi) dphi = fmod(pi+dphi, twopi) - pi;
  else if(dphi < -pi) dphi = -fmod(pi-dphi, twopi) + pi;

  double r   = sqrt(deta*deta+dphi*dphi);

  //if (r < 1 && r/r1 > 0.7) {
  //  cout<<" MW: r/r1= "<<r<<"  "<<r/r1<<endl;
  //  cout<<" MW  ------  p  "<<_M_p[1]<<endl;
  //}

  return (r/r1);
}

void cone_et::_M_merge_Et(int i, int j)
{
  static const double pi = 3.14159265358979323846;
  static const double twopi = 6.28318530717958647692;

  double pt1 = _M_p[i].perp(), pt2 = _M_p[j].perp();
  double eta1 = _M_p[i].prapidity(), eta2 = _M_p[j].prapidity();
  double phi1 = _M_p[i].phi(), phi2 = _M_p[j].phi();
  double dphi = phi1 - phi2;

  if(dphi >= pi) phi1=phi1-twopi;
  if(-dphi >= pi) phi2=phi2-twopi;

  double pts = pt1 + pt2;
  double etas = 0.0;
  double phis = 0.0;
  if(pts>0.0000000000001) {
    etas = (eta1*pt1 + eta2*pt2)/pts;
    phis = (phi1*pt1 + phi2*pt2)/pts;
  }
  else{
    etas = 0.0;
    phis = 0.0;
  }

  double px = pts*std::cos(phis);
  double py = pts*std::sin(phis);
  double pz = pts*std::sinh(etas);
  double ee = pts*std::cosh(etas);
  // MW: fill new result in first vector
  _M_p[i] = _Lv(px , py, pz, ee);
}
