#ifndef __cteq5m_h__
#define __cteq5m_h__ 1

#include <process.h>

//----- used namespaces -----
using namespace nlo;
using namespace std;



//   Parton distribution function and the strong coupling.
//   In this example we use the mrst99 padf set and the two
//   loop  qcd coupling. This class is inhereted from the 
//   'amplitude_dis::pdf_and_coupling' abstact class and 
//   it must have one compusory member function.
//
//    weight_dis pdf(double x, double mf2, unsigned nu, unsigned nd);
//    returns the parton distributions
//         x = input : momentum fraction
//         mf2 = input : square of the factorization scale
//         nu = input : number of the up type flavours
//         nd = input : number of the down type flavours
//
//    double  alpha_qcd(unsigned nf, double mr2);
//    returns the QCD coupling
//         nf = input : number of the active quark flavours
//         mr2 = input : square of the renormalization scale 
//
class pdf_cteq5m
  : public amplitude_dis::pdf_and_coupling 
{
public:
  //   pdf mode 
  enum mode_type { nlo = 1, lo = 3};

  //   constructor
  explicit pdf_cteq5m(unsigned int __mode = 1, unsigned int __loop = 2)
    : _M_mode(__mode), _M_loop(__loop) {}
  
  //   set mode
  void mode(int __mode) {
    _M_mode = __mode;
  }
  
  void loop(int __loop) {
    _M_loop = __loop;
  }
  
  //   get mode
  int mode() const { return _M_mode;}
  int loop() const { return _M_loop;}
  
  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) {
     //    double al = (_M_mode == nlo ? 0.01878028328484364962 : 0.02021267777267070764);
    double al = 0.01878028328484364962;
    return nlo::alpha_qcd(_M_loop, nf, 91.187, al, std::sqrt(mr2));
  }
      
  //   the parton distribution function
  weight_dis pdf(double, double, unsigned int = 2U, unsigned int = 3U);
  
private:
  int _M_mode, _M_loop;
};


#endif
