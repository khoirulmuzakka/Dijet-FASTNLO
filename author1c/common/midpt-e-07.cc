#include "midpt-e-07.h"
#include <cmath>

// MW test alphas
//#include "qcdlib.h"

//extern "C" void pxtest_(int *,int *,double vp[], double vj[]);
extern "C" void pxinterface_(int *,int *,double vp[], double vj[]);


const bounded_vector<lorentzvector<double> >&
midpt_e_07::operator()(const event_hhc& ev, double rcone)
{
  int nj = 0, np = 0, nt = ev.upper();
  double part[16],jet[16];

  //----- initialize the arrays -----
  _M_pj.resize(1,0); _M_p.resize(1, nt);
  for(int ip = 1; ip <= nt; ip++) 
    if(ev[ip].perp2() > 1.0e-12) _M_p[++np] = ev[ip];

  // copy particle 4-vectors to 1d array
  for(int ip = 1; ip <= np; ip++) {
    part[(ip*4-4)] =_M_p[ip].X();
    part[(ip*4-3)] =_M_p[ip].Y();
    part[(ip*4-2)] =_M_p[ip].Z();
    part[(ip*4-1)] =_M_p[ip].T();
  }

  // MW: call pxcone
  //pxtest_(&np,&nj,part,jet);
  pxinterface_(&np,&nj,part,jet);
  for(int i = 1; i <= nj; i++) {
    _M_p[i] = _Lv(jet[(i-1)*4],jet[(i-1)*4+1],jet[(i-1)*4+2],jet[(i-1)*4+3]);
    _M_pj.push_back(_M_p[i]);
  }

  return _M_pj;
}

