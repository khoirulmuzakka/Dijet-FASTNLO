#ifndef __fj_mc_05_h__
#define __fj_mc_05_h__ 1


#include <event.h>

using namespace std;
using namespace nlo;

class fj_mc_05
{
  
  //   private types
  typedef lorentzvector<double> _Lv;
  
public:
  //   do the clustering and return with the momenta of the jets
  const bounded_vector<_Lv>& operator()(const event_hhc&, double = 0.5);
  
private:
  //   private data members
  bounded_vector<_Lv> _M_p;
  bounded_vector<_Lv> _M_pj;
  bounded_vector<_Lv> _M_ax;
  
  //   private members
  void _M_merge(int i, int j) { _M_p[i] += _M_p[j];}
};

#endif
