#ifndef __fj_ak_h__
#define __fj_ak_h__ 1

#include <bits/hhc-event.h>
#include <bits/hep-bounded_vector.h>

using namespace std;
using namespace nlo;

class fj_ak {

   //   private types
   typedef lorentzvector<double> _Lv;

 public:
   //   do the clustering and return with the momenta of the jets
   const bounded_vector<_Lv>& operator()(const event_hhc&, double);

 private:
   //   private data members
   bounded_vector<_Lv> _M_pj;

};

#endif
