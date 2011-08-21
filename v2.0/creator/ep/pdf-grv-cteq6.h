#ifndef __pdf_grv_cteq6_h__
#define __pdf_grv_cteq6_h__ 1


#include <algorithm>
#include <bits/hhc-process.h>
#include <cteq6.h>


//----- used namespaces -----
using namespace nlo;
using namespace std;



class pdf_grv_cteq6
  : public pdf_and_coupling_hhc
{
public:
  //   constructor
  pdf_grv_cteq6(unsigned int mem = 0) 
	: _M_pdf(mem), _M_ymin(0.0), _M_ymax(1.0) 
	  { _S_gauleg(20, _M_xb, _M_wb);}
  
  //  set the allowed momentum fraction range of the photon
  void photon_momentum_fraction(double ymin, double ymax) {
	_M_ymin = ymin; _M_ymax = ymax;
  }

  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) {
     return _M_pdf(std::sqrt(mr2))/6.28318530717958647692;
  }
  
  //   parton in lepton distribution function
  void hadronA(double x, double Q2, unsigned int, unsigned int, double *f) 
  {
    double func[6];
    parton_in_lepton(x, Q2, func);
    
    f[0] = func[0]/(137.0*x);
    for(unsigned int i = 1; i < 6; i++) 
       f[i] = f[-i] = func[i]/(137.0*x);
    f[6] = f[-6] = 0.0;
  }
    
  //   the parton in the hadron distribution function
  void hadronB(double x, double Q2, unsigned int, unsigned int, double *f) {
    double __f[13]; _M_pdf(x, sqrt(Q2), __f+6);
    for(int i=-6; i <= 6; i++) f[i] = __f[6+i]/x;
  }
  
private:
  //    proton distribution function (CTEQ6)
  lhpdf::cteq6 _M_pdf;
  
  //    gauss-legendre base point and weights
  //    20 base point must be enough but if you 
  //    think you need more fell free to change 
  double _M_xb[20], _M_wb[20];
  
  //   calculates the base points and weights
  void _S_gauleg(unsigned int, double *, double *);

  //  photon in electron (Williams-Weicheker function)
  double photon_in_lepton(double); 
  double _M_ymin, _M_ymax;

  //   GRV parton distribution in photon (c++ wrapper of the f77 function)
  void parton_in_photon(double, double, double *);

  //  paton_in_lepton = photon_in_lepton \otimes parton_in_photon
  void parton_in_lepton(double, double, double *);
};



#endif
