#ifndef __pdf_hhc_h__
#define __pdf_hhc_h__ 1

#include <algorithm>
#include <bits/hhc-process.h>


//----- used namespaces -----
using namespace nlo;
using namespace std;



class pdf_hhc_dummy
  : public pdf_and_coupling_hhc
{
public:
  //   constructor
  pdf_hhc_dummy(unsigned int mem = 0) { ;}

  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) {
    return 1.0;
  }

  //   parton in lepton distribution function
  void hadronA(double x, double Q2, unsigned int, unsigned int, double *f)
  {
      printf("pdf_hhc_dummy::hadronA : PDF in dummy should never be called.\n");
      exit(1);

  }

   //   the parton in the hadron distribution function
   void hadronB(double x, double Q2, unsigned int, unsigned int, double *f) {
      printf("pdf_hhc_dummy::hadronB : PDF in dummy should never be called.\n");
      exit(1);
   }

   weight_hhc pdf(double x1, double x2, double mf2, unsigned int nu, unsigned int nd)
   {
      weight_hhc retval;

      retval[0] = 1./x1/x2;
      retval[1] = 1./x1/x2;
      retval[2] = 1./x1/x2;
      retval[3] = 1./x1/x2;
      retval[4] = 1./x1/x2;
      retval[5] = 1./x1/x2;
      retval[6] = 1./x1/x2;

      for (int i=0; i<7; i++) {
         if (isnan(retval[i])){
            cout << "fastNLO.pdf_hhc_dummy: WARNING! NaN for 1/x1/x2 in no. " << i << ", x1 = " << x1 << ", x2 = " << x2 << endl;
         }
      }
      return retval;
   }
};
#endif
