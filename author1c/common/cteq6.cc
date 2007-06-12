#include "cteq6.h"

extern "C" double ctq6pdf_(int *, double *, double *);
extern "C" int setctq6_(int *);

weight_hhc pdf_cteq6::
pdf(double e1, double e2, double muf2, unsigned int nu, unsigned int nd)
{
  weight_hhc retval;
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
  double A0 = ctq6pdf_(&ig, &e1, &mf);
  double B0 = ctq6pdf_(&ig, &e2, &mf);
  
  //---- up type quarks -----
  double q1, q2, a1, a2;
  double A = 0.0, B = 0.0, Ab = 0.0, Bb = 0.0, D = 0.0, Db = 0.0; 

  for(int u = 0; u < (int) nu && u < 2; u++) {
    ia = -(iq = iu[u]);
    q1 = ctq6pdf_(&iq, &e1, &mf);
    a1 = ctq6pdf_(&ia, &e1, &mf);

    a2 = ctq6pdf_(&iq, &e2, &mf);
    q2 = ctq6pdf_(&ia, &e2, &mf);
    
    A += q1; Ab += a1; B += q2; Bb += a2;
    D += q1*q2 + a1*a2; Db += q1*a2 + a1*q2;
  }

  //----- down type quarks -----
  for(int d = 0; d < (int) nd && d < 3; d++) {
    ia = -(iq = id[d]);
    q1 = ctq6pdf_(&iq, &e1, &mf);
    a1 = ctq6pdf_(&ia, &e1, &mf);

    a2 = ctq6pdf_(&iq, &e2, &mf);
    q2 = ctq6pdf_(&ia, &e2, &mf);
    
    A += q1; Ab += a1; B += q2; Bb += a2;
    D += q1*q2 + a1*a2; Db += q1*a2 + a1*q2;
  }

  // cout<<" MW           e  "<<e1<<"  "<<e2<<endl;

  if(e1<xmin) {
    xmin=e1;
    cout<<" MW   xmin= "<<xmin<<endl;
  }
  if(e2<xmin) {
    xmin=e2;
    cout<<" MW   xmin= "<<xmin<<"   "<<endl;
  }

  //cout<<" MW  in CTEQ6    mu2= "<<muf2<<"  x1,2= "<<e1<<"  "<<e2<<"  "<<xmin<<endl;
  retval[0] = A0*B0;
  retval[1] = (A + Ab)*B0;
  retval[2] = A0*(B + Bb);
  retval[3] = A*B + Ab*Bb - D;
  retval[4] = D;
  retval[5] = Db;
  retval[6] = A*Bb +Ab*B - Db;

  //cout<<" MWcteq6:  "<<retval[0]<<"  "<<retval[1]<<endl;

  return retval;
}
