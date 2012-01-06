#ifndef __fnloBlockBNlojet__
#define __fnloBlockBNlojet__


#include "fnloBlockB.h"
#include <bits/dis-process.h>
#include <bits/photo-process.h>
#include <bits/hhc-process.h>
#include "fnloBlockA2.h"

class fnloBlockBNlojet : public fnloBlockB {
 public:
   fnloBlockBNlojet(fnloBlockA1 *blocka1, fnloBlockA2 *blocka2) :fnloBlockB(blocka1,blocka2){_S_gauleg(20, _M_xb, _M_wb);counter=0;}
   void FillEventDIS(int ObsBin,double x, double scale1, const nlo::amplitude_dis& amp, nlo::pdf_and_coupling_dis& pdf, double prefactor=1.0);
   void FillEventDISMuVar(int ObsBin, double x, double M1, double M2, const nlo::amplitude_dis& amp, nlo::pdf_and_coupling_dis& dummypdf, nlo::pdf_and_coupling_dis& realpdf, double prefactor=1.0);
   void FillEventPhoto(int ObsBin,double x, double scale1, const nlo::amplitude_photo& amp, nlo::pdf_and_coupling_photo& pdf, double prefactor=1.0);
   void FillEventResolved(int ObsBin,double x1, double x2,double scale1, double ymin, double ymax,  double Q2max,const nlo::amplitude_hhc& amp, nlo::pdf_and_coupling_hhc& pdf, double prefactor=1.0);
   void FillEventHHC(int ObsBin,double x1, double x2,double scale1,const nlo::amplitude_hhc& amp, nlo::pdf_and_coupling_hhc& pdf, double prefactor=1.0);
   void FillEventHHCMuVar(int ObsBin,double x1, double x2, double M1, double M2, const nlo::amplitude_hhc& amp, nlo::pdf_and_coupling_hhc& dummypdf, nlo::pdf_and_coupling_hhc& realpdf, double prefactor=1.0);

   void InitDISConstants( fnloBlockA2* A2 , bool nlo );
   void InitFinalDISValues( fnloBlockA2* A2 , double* xlim , double* scale1lo , double* scale1hi , double* scale2lo = NULL , double* scale2hi = NULL );
   void InitLogLogScaleNode( fnloBlockA2* A2 , double* slo , double* shi , int iScale );
   void InitLinearScaleNode( fnloBlockA2* A2 , double* slo , double* shi , int iScale );
   void ResizeSigmaTildeTables( fnloBlockA2* A2 );
   void InitLHCConstants( fnloBlockA2* A2 , bool nlo );
   void InitReferenceTable( fnloBlockA2* A2 );
   void SetNumberOfXNodesPerMagnitude( int nxPerMagnitude , double* xlim );

   void SetScale1Name( string name );
   void SetScale2Name( string name );

   void SetNumberOfScaleNodesScale1( int nNodes ) { NscalenodeScale1 = nNodes;  };
   void SetNumberOfScaleNodesScale2( int nNodes ) { NscalenodeScale2 = nNodes;  };
   void SetNumberOfScaleNodes_v20( int nNodes ) { Nscalenode.resize(1); Nscalenode[0] = nNodes;  };

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
   enum EInterpolKernel { kCatmulRom = 1 , kLagrangian = 2 };
   void Interpol(int nnode, int nmax, double delta, int ikern, int &nmod, vector <double> *kernel);
   vector < double > Interpol(int nnode, int nmax, double delta, EInterpolKernel eIkern, int &nmod );

   void WarmUp( int ObsBin, double x, double M1, double M2 = 0 , string sx = "xlim", string s1 ="mu", string s2 ="");
   void FillMuVarReferenceTables(int ObsBin, double M1, double M2, const nlo::amplitude_dis& amp, nlo::pdf_and_coupling_dis& realpdf, double prefactor);
   void FillMuVarReferenceTables(int ObsBin, double M1, double M2, const nlo::amplitude_hhc& amp, nlo::pdf_and_coupling_hhc& realpdf, double prefactor);
   unsigned long counter; // here, if you want to have multiple tables in a warm up run
   double* axlo;
   double* a1lo;
   double* a1up;
   double* a2lo;
   double* a2up;
   

   void InitScaleNode( fnloBlockA2* A2, double* slo, double* shi, int iScale  );
   double (fnloBlockBNlojet::*Fct_H_Scale[2])(double); // an array of pointers to member functions
   double (fnloBlockBNlojet::*Fct_H_Scale_Inv[2])(double); // an array of pointers to member functions

   double Function_loglog025( double mu );
   double Function_loglog025_inv( double mu );
   double Function_x( double mu );
   double Function_x_inv( double mu );


public:
   // variables for warm-up run
   int IWarmUp;
   unsigned long IWarmUpCounter;
   unsigned long IWarmUpPrint;
   vector < double > xlo;
   vector < double > scalelo;
   vector < double > scalehi;

   vector < double > scale1hi;
   vector < double > scale1lo;
   vector < double > scale2hi;
   vector < double > scale2lo;

};


#endif
