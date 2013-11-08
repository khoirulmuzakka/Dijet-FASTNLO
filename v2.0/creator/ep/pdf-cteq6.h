#ifndef __pdf_cteq6_h__
#define __pdf_cteq6_h__ 1


#include <bits/photo-process.h>
#include <bits/dis-process.h>
#include <cteq6.h>


//----- used namespaces -----
using namespace nlo;
using namespace std;
using namespace lhpdf;


class pdf_cteq6photo
  : public pdf_and_coupling_photo
{
public:
  //   constructor
  explicit pdf_cteq6photo(unsigned int mem = 0)
    : _M_pdf(mem) {}
  
  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) {
     return _M_pdf(std::sqrt(mr2))/6.28318530717958647692;
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
     double __f[13]; _M_pdf(x, sqrt(Q2), __f+6);
     for(int i=-6; i <= 6; i++) f[i] = __f[6+i]/x;
  }
  
private:
  lhpdf::cteq6 _M_pdf;
};

class pdf_cteq6dis
  : public pdf_and_coupling_dis
{
public:
  //   constructor
  explicit pdf_cteq6dis(unsigned int mem = 0)
    : _M_pdf(mem) {}
  
  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) {
     return _M_pdf(std::sqrt(mr2))/6.28318530717958647692;
  }
  
  //   the parton distribution function
  void hadron(double x, double Q2, unsigned int, unsigned int, double *f) {
    double __f[13]; _M_pdf(x, sqrt(Q2), __f+6);
    for(int i=-6; i <= 6; i++) f[i] = __f[6+i]/x;
  }
  
private:
  lhpdf::cteq6 _M_pdf;
};



#endif
