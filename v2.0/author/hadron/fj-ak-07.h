#ifndef __fj_ak_07_h__
#define __fj_ak_07_h__ 1


#include <bits/hhc-event.h>
#include <bits/hep-bounded_vector.h>

using namespace std;
using namespace nlo;

class fj_ak_07
{
  
  //   private types
  typedef lorentzvector<double> _Lv;
  
public:
  //   do the clustering and return with the momenta of the jets
  const bounded_vector<_Lv>& operator()(const event_hhc&, double = 0.7);
  
private:
  //   private data members
  bounded_vector<_Lv> _M_p;
  bounded_vector<_Lv> _M_pj;
  bounded_vector<_Lv> _M_ax;
  
  //   private members
  void _M_merge(int i, int j) { _M_p[i] += _M_p[j];}
};

#endif
