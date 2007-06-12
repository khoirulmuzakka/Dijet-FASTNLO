#ifndef __kt_et_10_h__
#define __kt_et_10_h__ 1


#include <event.h>

using namespace std;
using namespace nlo;


class kt_et_10
{
  //   private types
  typedef lorentzvector<double> _Lv;
  
  struct _Vec {
    _Vec() {}
    _Vec(const _Lv& p)
      : pt(p.perp()), pt2(p.perp2()),
	eta(p.prapidity()), phi(p.phi()) {}
    
    double pt, pt2, eta, phi;
  };
 
public:
  //   do the clustering and return with the momenta of the jets
  const bounded_vector<_Lv>& operator()(const event_dis&, double = 1.0);
  
private:
  //   private data members
  bounded_vector<_Vec> _M_p;
  bounded_vector<_Lv> _M_pj;

  //  convert : _Vec --> _Lv
  _Lv _Vec_to_Lv(const _Vec& v) {
    return v.pt*_Lv(cos(v.phi), sin(v.phi), sinh(v.eta), cosh(v.eta));
  }
  
  //   private members
  double _M_pair(int, int);
  void  _M_merge(int, int);
};


#endif
