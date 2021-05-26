#ifndef __cone_et_h__
#define __cone_et_h__ 1

#include <bits/hhc-event.h>
#include <bits/hep-bounded_vector.h>

class cone_et {
  //   private types
  typedef nlo::lorentzvector<double> _Lv;

public:
  //   do the clustering and return with the momenta of the jets
  const nlo::bounded_vector<_Lv> &operator()(const nlo::event_hhc &, double);

private:
  //   private data members
  nlo::bounded_vector<_Lv> _M_p;
  nlo::bounded_vector<_Lv> _M_pj;

  //   private members
  double _M_pair(int, int);
  void _M_merge(int i, int j) { _M_p[i] += _M_p[j]; }
  // MW: merge in Et scheme
  void _M_merge_Et(int, int);
};

#endif
