#ifndef __pdf_cteq6_h__
#define __pdf_cteq6_h__ 1


#include <bits/hhc-process.h>
#include <cteq6.h>


//----- used namespaces -----
using namespace nlo;
using namespace std;
using namespace lhpdf;


class pdf_cteq6
  : public pdf_and_coupling_hhc
{
public:
  //   constructor
  explicit pdf_cteq6(unsigned int mem = 0)
    : _M_pdf(mem) {}
  
  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) {
    return _M_pdf(std::sqrt(mr2))/6.28318530717958647692;
  }
  
  //   the parton distribution function
  void hadronA(double x, double Q2, unsigned int, unsigned int, double *f) {
    double __f[13]; _M_pdf(x, sqrt(Q2), __f+6);
    for(int i=-6; i <= 6; i++) f[i] = __f[6+i]/x;
  }
  
  void hadronB(double x, double Q2, unsigned int, unsigned int, double *f) {
    double __f[13]; _M_pdf(x, sqrt(Q2), __f+6);
    for(int i=-6; i <= 6; i++) f[i] = __f[6-i]/x;
  }
  
private:
  lhpdf::cteq6 _M_pdf;
};



#endif
