#include "fnlo_int_nlojet/pdf-hhc-dummy.h"

using namespace nlo;
using namespace std;

weight_hhc pdf_hhc_dummy::pdf(double x1, double x2, double mf2, unsigned int nu,
                              unsigned int nd) {
  weight_hhc retval;

  retval[0] = 1. / x1 / x2;
  retval[1] = 1. / x1 / x2;
  retval[2] = 1. / x1 / x2;
  retval[3] = 1. / x1 / x2;
  retval[4] = 1. / x1 / x2;
  retval[5] = 1. / x1 / x2;
  retval[6] = 1. / x1 / x2;

  for (int i = 0; i < 7; i++) {
    if (isnan(retval[i])) {
      cout << "fastNLO.pdf_hhc_dummy: WARNING! NaN for 1/x1/x2 in no. " << i
           << ", x1 = " << x1 << ", x2 = " << x2 << endl;
    }
  }
  return retval;
}
