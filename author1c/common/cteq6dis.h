#ifndef __cteq6_h__
#define __cteq6_h__ 1


#include <process.h>


//----- used namespaces -----
using namespace nlo;
using namespace std;


//   Parton distribution function and the strong coupling.
//   In this example we use the mrst99 pdf set and the two  -- CTEQ6
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
extern "C" double alfas_(double*,double*,int*);

class pdf_cteq6
  : public amplitude_dis::pdf_and_coupling 
{
public:
  //   pdf mode 
  enum mode_type { nlo = 200, lo = 3};

  //   constructor
  explicit pdf_cteq6(mode_type md = pdf_cteq6::nlo, unsigned int lp = 2)
    : _M_mode(md), _M_loop(lp) {}
  
  //   set mode
  void mode(mode_type md) { _M_mode = md;}
  void loop(int lp) { _M_loop = lp;}
  
  //   get mode
  mode_type mode() const { return _M_mode;}
  int loop() const { return _M_loop;}
  
  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) {
    double al = 0.01878028328484364962;
    //    double al = 0.01971929;     // 5% higher
    // if(_M_mode == lo) al = 0.02021267777267070764;
    return nlo::alpha_qcd(_M_loop, nf, 91.187, al, std::sqrt(mr2));
  }
  
  //   the parton distribution function
  weight_dis pdf(double, double, unsigned int = 2U, unsigned int = 3U);
  
private:
  mode_type _M_mode;
  int _M_loop;
};



#endif
