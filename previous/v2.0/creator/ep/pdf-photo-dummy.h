#ifndef __pdf_photo_dummy_h__
#define __pdf_photo_dummy_h__ 1


#include <bits/photo-process.h>


//----- used namespaces -----
using namespace nlo;
using namespace std;
using namespace lhpdf;


class pdf_photo_dummy
  : public pdf_and_coupling_photo
{
public:
  //   constructor
   explicit pdf_photo_dummy(unsigned int mem = 0){;}
  
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
      printf("pdf_photon_dummy::hadron : PDF in dummy should never be called.\n");
      exit(1);
   }

  weight_photo pdf(double x1, double x2, double mf2, unsigned int nu, unsigned int nd) {
     weight_photo retval;
     //----- gluon pdfs -----
     retval[0] = 1./x2;
     //---- up type quarks -----
     retval[1] = 1./x2;
     //---- down type quarks -----
     retval[2] = 1./x2;
     return retval*(this -> photon(x1));
    
  }  

   
  
};



#endif
