#include "kt-e-10.h"
#include <cmath>


const bounded_vector<lorentzvector<double> >&
kt_e_10::operator()(const event_hhc& ev, double r0)
{
  int imin, jmin, kmin, nt = ev.upper();
  double tmp, pmin, smin, r2 = r0*r0;
  
  //----- initialize the arrays -----
  _M_pj.resize(1,0); _M_p.resize(1, nt);
  for(int ip = 1; ip <= nt; ip++)
    _M_p[ip] = ev[ip];
  
  for(int n = nt; n > 1; n--) {
    //----- find the smallest resolution variables -----
    imin = 1; jmin = 2; kmin = 1;
    pmin = _M_pair(1,2); smin = _M_p[1].perp2();
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
    
    //----- do the clustering -----
    if(smin*r2 < pmin) {
      _M_pj.push_back(_M_p[kmin]);
      if(kmin != n) _M_p[kmin] = _M_p[n];
    } else {
      _M_merge(imin, jmin);
      if(jmin != n) _M_p[jmin] = _M_p[n];
    }
  }
  
  //----- add the last (pseudo)particle to the list of jets ----- 
  _M_pj.push_back(_M_p[1]);

  return _M_pj;
}

double kt_e_10::_M_pair(int i, int j)
{
  static const double pi = 3.14159265358979323846;
  static const double twopi = 6.28318530717958647692;
  
  double et2 = (_M_p[i].perp2() < _M_p[j].perp2() ? _M_p[i].perp2() : _M_p[j].perp2());
  double dy = _M_p[i].rapidity() - _M_p[j].rapidity();
  double dphi = _M_p[i].phi() - _M_p[j].phi();
  
  if(dphi >= pi) dphi = fmod(pi+dphi, twopi) - pi;
  else if(dphi < -pi) dphi = -fmod(pi-dphi, twopi) + pi;

  return et2*(dy*dy+dphi*dphi);
}

void kt_e_10::_M_merge(int i, int j)
{
   _M_p[i] += _M_p[j];

}
