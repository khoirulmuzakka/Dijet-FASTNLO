#include "pdf-grv-cteq6.h"


extern "C" int grvglo_(float*,float*,float*,float*,float*,float*,float*,float*);
extern "C" int grvgho_(float*,float*,float*,float*,float*,float*,float*,float*);



double pdf_grv_cteq6::photon_in_lepton(double x) 
{
  const double Me2 = 0.00000026112004954086;  //  GeV^2
  double Q2max = 1.0, Q2min = Me2*x*x/(1-x);
  
  return 1.0/137.0*((1+(1-x)*(1-x))/x*log(Q2max/Q2min) 
					+ 2.0*Me2*x*(1.0/Q2max-1.0/Q2min))/6.28318530717958647692; 
}

 
void pdf_grv_cteq6::parton_in_photon(double x, double Q2, double *f)
{
  float __x = x, __q2 = Q2, uh, dh, sh, ch, bh, gh;

  grvgho_(&__x, &__q2, &uh, &dh, &sh, &ch, &bh, &gh);
  
   f[0] = gh; f[1] = dh; f[2] = uh; 
   f[3] = sh; f[4] = ch; f[5] = bh;

}


void pdf_grv_cteq6::parton_in_lepton(double x, double Q2, double *ret_val)
{
  double ymin = x > _M_ymin ? x : _M_ymin;
  
  std::memset(ret_val, 0, 6*sizeof(double));
  if(x >= _M_ymax) return;

  double weight, ph[6];
  double y, jac = 0.5*std::log(_M_ymax/ymin);
  
  for(unsigned int ib = 0; ib < 20; ib++) {
    y = ymin*std::exp((_M_xb[ib]+1)*jac);
    
    weight = jac*_M_wb[ib]*y*(this -> photon_in_lepton(y));
    this -> parton_in_photon(x/y, Q2, ph);
    
    for(unsigned int ip = 0; ip < 6; ip++)
      ret_val[ip] += weight*ph[ip];
  }

}


#define EPS 1.0e-15

void pdf_grv_cteq6::_S_gauleg(unsigned int n, double *x, double *w)
{
  unsigned int m, j, i;
  double z1, z, pp, p3, p2, p1;
  
  m = (n+1)/2;
  for(i = 1; i <= m; i++) {
    z = cos(3.14159265358979323846*(i-0.25)/(n+0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for(j = 1; j <= n; j++) {
		p3 = p2;
		p2 = p1;
		p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp = n*(z*p1-p2)/(z*z-1.0);
      z1 = z;
      z = z1-p1/pp;
    } while(fabs(z-z1) > EPS);
    
    x[i-1] = -z, x[n-i] = z;
    w[i-1] = w[n-i] = 2.0/((1.0-z*z)*pp*pp);
  }
}



