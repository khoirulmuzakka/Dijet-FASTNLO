// Author: Daniel Britzger
// DESY, 02/04/2012

#ifndef FASTNLODIFFDiffFit
#define FASTNLODIFFDiffFit


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  FastNLODiffUSER                                                     //
//                                                                      //
//  FastNLODiffReader is a standalone code for reading                  //
//  diffractive FastNLO tables of version 2.0 for DIS processes         //
//                                                                      //
//  FastNLO is developed by                                             //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch         //
//    (publication in preparation)                                      //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <string>
#include <cstdio>
#include <vector>
#include "FastNLODiffReader.h"

using namespace std;

class FastNLODiffDiffFit : public FastNLODiffReader {

public:

   FastNLODiffDiffFit(string filename);
   ~FastNLODiffDiffFit(void){;};

protected:

   // inherited functions
   double EvolveAlphas(double Q) const ;
   bool InitPDF();
   vector<double> GetDiffXFX(double xpom, double zpom, double muf) const ;

};



//______________________________________________________________________________


extern "C"{
   void fpdfxq_(int *iset, const double *x, const double *q2, double *pdfs,int *ichk);
   double asfunc_( double* r2, int* nf  , int* ierr);
   //void diffpdf_(double* xpom, double*  zpom, double*  Q2, double *pdfs);

   void evolfg_(int* itype, double(*func)(double*,double*), double (*)[13] , int* iq0, double* epsi) ;
   double func_(double*,double*);
   extern struct{
      double xpout;
   } xpfeed_;
   extern struct{
      int xpselect;
   } xpselection_;

   //common/misc/proton,def,iq0
   extern struct{
      double proton[13];
      double def[12][13];
      int iq0;
   } misc_;

}



//______________________________________________________________________________


FastNLODiffDiffFit::FastNLODiffDiffFit(string filename) : FastNLODiffReader(filename)
{
}


//______________________________________________________________________________


double FastNLODiffDiffFit::EvolveAlphas(double Q) const {
   // --- fastNLO user:
   // Implementation of Alpha_s evolution as function of the
   // factorization scale [and alphas(Mz)].
   //
   double mu2 = Q*Q;
   int ierr = 9876;
   int nf = 9;
   double as = asfunc_( &mu2, &nf  , &ierr);
   if ( ierr > 0 )
      error["EvolveAlphas"]<<"Alphas evolution failed. ierr = "<<ierr<<", Q = "<<Q<<endl;
   return as;
}


//______________________________________________________________________________


bool FastNLODiffDiffFit::InitPDF(){
   // --- fastNLO user:
   //  Initalize PDF parameters if necessary
   //
   // nothing todo!
   // we init the xpom value in the common block
   xpselection_.xpselect = 1;
   xpfeed_.xpout = fxpom;
   //

   //call EVOLFG(1,func,def,iq0,epsi)     ! Evolve dpdf's on the grid
   double epsi = 1.;
   int itype = 1;
   evolfg_(&itype,&func_, misc_.def ,&misc_.iq0,&epsi);
   return true;
}


//______________________________________________________________________________



vector<double> FastNLODiffDiffFit::GetDiffXFX(double xpom, double zpom, double muf) const {
   //
   //  GetDiffXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  xpom, zpom and factorisation scale.
   //

   vector < double > xfx(13);

   double xp = zpom;
   int iqnset = 1;
   int iqnchk = 0;
   double muf2 = muf*muf;
   fpdfxq_(&iqnset, &xp, &muf2, &xfx[0], &iqnchk);

   return xfx;

}


//______________________________________________________________________________


#endif
