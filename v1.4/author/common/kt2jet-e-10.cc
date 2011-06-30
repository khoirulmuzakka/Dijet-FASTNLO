#include "kt2jet-e-10.h"
#include <cmath>

const bounded_vector<lorentzvector<double> >&
kt2jet_e_10::operator()(const event_hhc& ev, double rcone)
{
  int merge = 0, nj = 0, np = 0, nt = ev.upper();
  double dist;

  //----- initialize the arrays -----
  _M_pj.resize(1,0); _M_p.resize(1, nt); 
  _M_ax.resize(1,1);
  for(int ip = 1; ip <= nt; ip++) 
    if(ev[ip].perp2() > 1.0e-12) _M_p[++np] = ev[ip];

  nj = np;
  for(int i = 1; i <= np; i++) {
    for(int j = i + 1; j <= np; j++) {     
      if(merge == 0 && (dist = _M_pair(i,j)) < 1.0) {
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

  return _M_pj;
}

double kt2jet_e_10::_M_pair(int i, int j)
{
  static const double pi = 3.14159265358979323846;
  static const double twopi = 6.28318530717958647692;
  static const double rcone = 0.7;
  static const double rsep = 1.0;
  int k; 
  double rx;

  
   //   - check if distance between both particles is < R_sep*R_cone)
  double dy = _M_p[i].rapidity() - _M_p[j].rapidity();
  double dphi = _M_p[i].phi() - _M_p[j].phi();  
  if(dphi >= pi) dphi = fmod(pi+dphi, twopi) - pi;
  else if(dphi < -pi) dphi = -fmod(pi-dphi, twopi) + pi;
  double r   = sqrt(dy*dy+dphi*dphi);

  rx=99.9;
  if (r<(rcone*rsep)){
    double pti = _M_p[i].perp(), ptj = _M_p[j].perp();
    if (pti<ptj){
       k=i;
    }else{
       k=j;
    }

    // compute potential jet axis
    _M_ax[1] = _M_p[i] + _M_p[j];

    // check only distance of lower pT particle to jet axis
    //           (higher pT particle always has smaller distance!) 
    double dyx = _M_p[k].rapidity() - _M_ax[1].rapidity();
    double dphix = _M_p[k].phi() - _M_ax[1].phi();  
    if(dphix >= pi) dphix = fmod(pi+dphix, twopi) - pi;
    else if(dphix < -pi) dphix = -fmod(pi-dphix, twopi) + pi;
    rx = sqrt(dyx*dyx+dphix*dphix);

    // test 
    // if (rx>r) cout<<" MW:xxxxxxxxxxxxxx r "<<r<<"  "<<rx<<endl;
    // if (rx/r < 0.5)    cout<<" MW: r/rx "<<r<<"  "<<rx<<"  "<<pti/ptj<<endl;
  }

  return (rx/rcone);
}

