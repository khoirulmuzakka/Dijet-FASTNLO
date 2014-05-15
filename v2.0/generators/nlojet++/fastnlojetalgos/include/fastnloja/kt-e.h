#ifndef __kt_e_h__
#define __kt_e_h__ 1


#include <bits/hhc-event.h>
#include <bits/hep-bounded_vector.h>

using namespace std;
using namespace nlo;

class kt_e
{
  typedef lorentzvector<double> _Lv;

public:
   //   do the clustering and return with the momenta of the jets
   const bounded_vector<_Lv>& operator()(const event_hhc&, double);

private:
  //   private data members
   bounded_vector<_Lv> _M_p;
   bounded_vector<_Lv> _M_pj;

   //   private members
   double _M_pair(int, int);
   void  _M_merge(int, int);
};


#endif
