#ifndef __fnloBlockBNlojet__
#define __fnloBlockBNlojet__


#include "fnloBlockB.h"
#include <bits/photo-process.h>
#include <bits/hhc-process.h>


class fnloBlockBNlojet : public fnloBlockB {
 public:
   fnloBlockBNlojet(fnloBlockA1 *blocka1, fnloBlockA2 *blocka2) :fnloBlockB(blocka1,blocka2){_S_gauleg(20, _M_xb, _M_wb);}
   void FillEventPhoto(int ObsBin,double x, double scale1, const nlo::amplitude_photo& amp, nlo::pdf_and_coupling_photo& pdf);
   void FillEventResolved(int ObsBin,double x1, double x2,double scale1, double ymin, double ymax,  double Q2max,const nlo::amplitude_hhc& amp, nlo::pdf_and_coupling_hhc& pdf);

private:
   //    gauss-legendre base point and weights
   //    20 base point must be enough but if you 
   //    think you need more fell free to change 
   double _M_xb[20], _M_wb[20];
   void _S_gauleg(unsigned int n, double *x, double *w);


};


#endif
