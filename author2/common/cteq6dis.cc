#include "cteq6dis.h"

extern "C" double ctq6pdf_(int *, double *, double *);
extern "C" int setctq6_(int *);

weight_dis pdf_cteq6::
pdf(double eta, double muf2, unsigned int nu, unsigned int nd)
{
  weight_dis retval;
  double mf = sqrt(muf2);
  int ia, iq, ig = 0, im = (int) _M_mode;
  static const int iu[2] = {1,4}, id[3] = {2,3,5};
  static bool first = true;
  static double xmin, xmax;

  if(first) {
    setctq6_(&im);
    cout<<"loading pdf table"<<endl;
    first = false;
    xmin = 1.0;
    xmax = 0.0;
  }

  //----- gluon pdfs -----
  retval[0] = ctq6pdf_(&ig, &eta, &mf);
  
  //---- up type quarks -----
  retval[1] = 0.0;
  for(int u = 0; u < (int) nu && u < 2; u++) {
     ia = -(iq = iu[u]);
     retval[1] += ctq6pdf_(&iq, &eta, &mf);
     retval[1] += ctq6pdf_(&ia, &eta, &mf);
  }


  //----- down type quarks -----
  for(int d = 0; d < (int) nd && d < 3; d++) {
     ia = -(iq = id[d]);
     retval[2] += ctq6pdf_(&iq, &eta, &mf);
     retval[2] += ctq6pdf_(&ia, &eta, &mf);
  }

//    retval[0] = 1.0;
//    retval[1] = 1.0;
//    retval[2] = 1.0;

  return retval;
}
