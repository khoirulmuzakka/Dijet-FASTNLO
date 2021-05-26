#ifndef __fj_jets_h__
#define __fj_jets_h__ 1

#include <bits/hhc-event.h>
#include <bits/hep-bounded_vector.h>

class fj_jets {

  //   private types
  typedef nlo::lorentzvector<double> _Lv;

public:
  //   do the clustering and return with the momenta of the jets
  const nlo::bounded_vector<_Lv> &operator()(const nlo::event_hhc &, int,
                                             double, double);

private:
  //   private data members
  nlo::bounded_vector<_Lv> _M_pj;
};

#endif
