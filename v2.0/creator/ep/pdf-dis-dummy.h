#ifndef __pdf_dis_h__
#define __pdf_dis_h__ 1

#include <algorithm>
#include <bits/dis-process.h>


//----- used namespaces -----
using namespace nlo;
using namespace std;



class pdf_dis_dummy
  : public pdf_and_coupling_dis
{
public:
  //   constructor
  pdf_dis_dummy(unsigned int mem = 0) { ;}
  
  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) {
    return 1.0;
  }
  
  void hadron(double x, double Q2, unsigned int, unsigned int, double *f) 
  {
      printf("pdf_dis_dummy::hadron : PDF in dummy should never be called.\n");
      exit(1);

  }
    
   weight_dis pdf(double x1, double mf2, unsigned int nu, unsigned int nd) 
   {
      weight_dis retval;

      retval[0] = 1./x1;
      retval[1] = 1./x1;
      retval[2] = 1./x1;
    
      return retval;
   }

  

};



#endif
