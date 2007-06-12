#include "cteq5m.h"



extern "C" double ctq5pd_(int *, int *, double *, double *, int *);

weight_dis 
pdf_cteq5m::pdf(double eta, double muf2, unsigned int nu, unsigned int nd)
{
  double mf = sqrt(muf2);
  weight_dis retval;
  int ia, iq, ig = 0, ierror = 0;
  static const int iu[2] = {1,4}, id[3] = {2,3,5};
  
  retval[0] = ctq5pd_(&_M_mode, &ig, &eta, &mf, &ierror);
  
  retval[1] = 0.0;
  for(int u = 0; u < (int) nu && u < 2; u++) {
    ia = -(iq = iu[u]);
    retval[1] += ctq5pd_(&_M_mode, &iq, &eta, &mf, &ierror);
    retval[1] += ctq5pd_(&_M_mode, &ia, &eta, &mf, &ierror);
  }
  
  retval[2] = 0.0;
  for(int d = 0; d < (int) nd && d < 3; d++) {
    ia = -(iq = id[d]);
    retval[2] += ctq5pd_(&_M_mode, &iq, &eta, &mf, &ierror);
    retval[2] += ctq5pd_(&_M_mode, &ia, &eta, &mf, &ierror);
  }

  return retval;
}
