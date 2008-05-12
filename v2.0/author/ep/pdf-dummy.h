#ifndef __pdf_dummy_h__
#define __pdf_dummy_h__ 1


#include <bits/photo-process.h>


//----- used namespaces -----
using namespace nlo;
using namespace std;
using namespace lhpdf;


class pdf_dummy
  : public pdf_and_coupling_photo
{
public:
  //   constructor
   explicit pdf_dummy(unsigned int mem = 0){;}
  
  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) {
     return 1.0;
  }
  
  //  photon in electron (Williams-Weicheker function)
  double photon(double x) 
  {
    const double Me2 = 0.00000026112004954086;  //  GeV^2
    double Q2max = 1.0, Q2min = Me2*x*x/(1-x);

     return 1.0/137.0*((1+(1-x)*(1-x))/x*log(Q2max/Q2min) 
           + 2.0*Me2*x*(1.0/Q2max-1.0/Q2min))/6.28318530717958647692; 
  }
    
  //   the parton distribution function
   void hadron(double x, double Q2, unsigned int, unsigned int, double *f) {
       for(int i=-6; i <= 6; i++) f[i] = 0.;
       f[0] = 1./x;
       f[1] = 1./x;
       f[2] = 1./x;
  }
  
};



#endif
