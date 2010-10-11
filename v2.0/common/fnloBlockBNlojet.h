#ifndef __fnloBlockBNlojet__
#define __fnloBlockBNlojet__


#include "fnloBlockB.h"
#include <bits/dis-process.h>
#include <bits/photo-process.h>
#include <bits/hhc-process.h>


class fnloBlockBNlojet : public fnloBlockB {
 public:
   fnloBlockBNlojet(fnloBlockA1 *blocka1, fnloBlockA2 *blocka2) :fnloBlockB(blocka1,blocka2){_S_gauleg(20, _M_xb, _M_wb);}
   void FillEventDIS(int ObsBin,double x, double scale1, const nlo::amplitude_dis& amp, nlo::pdf_and_coupling_dis& pdf, double prefactor=1.0);
   void FillEventDIS2Scale(int ObsBin,double x, double scale1, double scale2, const nlo::amplitude_dis& amp, nlo::pdf_and_coupling_dis& pdf, double prefactor=1.0);
   void FillEventPhoto(int ObsBin,double x, double scale1, const nlo::amplitude_photo& amp, nlo::pdf_and_coupling_photo& pdf, double prefactor=1.0);
   void FillEventResolved(int ObsBin,double x1, double x2,double scale1, double ymin, double ymax,  double Q2max,const nlo::amplitude_hhc& amp, nlo::pdf_and_coupling_hhc& pdf, double prefactor=1.0);
   void FillEventHHC(int ObsBin,double x1, double x2,double scale1,const nlo::amplitude_hhc& amp, nlo::pdf_and_coupling_hhc& pdf, double prefactor=1.0);

   //double TransformHx1(double x){return log10(x);}
   //double TransformHx2(double x){return -sqrt(log10(1.0/x));}
   //double TransformHmu(double mu){return log(log(mu/0.25));}
   double PDFwgt(double x){double w=(1.-0.99*x)/sqrt(x); w = w*w*w; return w;}

private:
   // This is for gammaP in ep, for NLOJET++ integration
   //
   //    gauss-legendre base point and weights
   //    20 base point must be enough but if you 
   //    think you need more fell free to change 
   double _M_xb[20], _M_wb[20];
   void _S_gauleg(unsigned int n, double *x, double *w);

   // interpolation kernel
   void Interpol(int nnode, int nmax, double delta, int ikern, int &nmod, vector <double> *kernel);

};


#endif
