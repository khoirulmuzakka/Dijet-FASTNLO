#ifndef __kt_e_h__
#define __kt_e_h__ 1

#include <bits/hhc-event.h>
#include <bits/hep-bounded_vector.h>

class kt_e {
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
  void _M_merge(int, int);
};

#endif
