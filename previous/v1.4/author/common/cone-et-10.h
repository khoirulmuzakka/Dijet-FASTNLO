#ifndef __cone_et_10_h__
#define __cone_et_10_h__ 1


#include <event.h>

using namespace std;
using namespace nlo;


class cone_et_10
{
  //   private types
  typedef lorentzvector<double> _Lv;
  
public:
  //   do the clustering and return with the momenta of the jets
  const bounded_vector<_Lv>& operator()(const event_hhc&, double = 1.0);
  
private:
  //   private data members
  bounded_vector<_Lv> _M_p;
  bounded_vector<_Lv> _M_pj;
  
  //   private members
  double _M_pair(int, int);
  void _M_merge(int i, int j) { _M_p[i] += _M_p[j];}
  // MW: merge in Et scheme
  void _M_merge_Et(int, int);
};


#endif
