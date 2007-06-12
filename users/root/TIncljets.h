#ifndef __TIncljets__
#define __TIncljets__

#include "TObject.h"
#include <vector>
#include <string>
#include <fstream>
#include "nlo-weight.h"
#include "LHAPDFWrap.h"
#include "physconst.h"

using namespace std;

const int kdis   = 1;
const int kpp    = 2;
const int kppbar = 3;

const int kcferror  = 0;
const int kcflo     = 1;
const int kcfnlo    = 2;

const int cmarker = 1234567890; //used to separate section in the table

class TIncljets : public TObject{
 public:
   typedef nlo::weight<3U> weight_dis;
 
   int nQ2;       // No of Q2idity bins 
   double *Q2high;  //[nQ2] array for Q2 boundaries
   int *npt;       //[nQ2] No of pT bins in each Q2 range
   int nptmax;     // maximum number of pt bins
   vector< vector<double> >pthigh;   // array for pT boundaries
   double ptlow;    // lowest pt considered 
   
   int nscalevar;                 // number of scale variations (mu_r, mu_f) for NLO
   vector <double> murscale;        // overall scale factor for renormalization scale
   vector< vector<double> >murval;   // array for renormalization scale values
   vector <double> mufscale;        // overall scale factor for factorization scale
   vector< vector<double> >mufval;    // array for factorization scale values

   int ntot;      // no of xbins 
   vector< vector<double> >xlimit; // array for lower x limits
   vector <vector< vector <  vector <weight_dis >  >  > > disweights; //! array for the weights MAIN ARRAY
   vector <vector< vector < weight_dis >  > > pdf; //! array for the pdfs
   int nsubproc;   // no of subprocesses
       
   vector< double> nevents; // no of events calculated so far for each scale variation
   int nwrite; // no of events after to write out the table

   int ireaction;
   double      s;
   int iproc;
   int ialgo;
   double JetResol1;
   double JetResol2;
   int npow;
   int npowmax;
   int Oalphas;
   int ixscheme;
   int ipdfwgt;

   vector <vector< vector < double  >  > > xsection; // the result

   double q2binning[5];
   double ptbinning[5];

   vector <vector< vector < double  >  > > rebinned; // the result

   vector <vector< vector < double  >  > > reference; 
   

   LHAPDFWrap *lhapdf;
   TIncljets();
   double alphas(double Q, double alphasMZ);
   void ReadTable(string filename);
   void ReadPDF(string filename);
   void SetPDFSet(int set);
   void FillPDFCache(double muffactor);
   void ResetXsection();
   void CalcXsection(double asmz, int scale=0);
   double GetXsection(int Q2, int pt, int order);
   void Rebin();
   double GetRebinned(int Q2, int pt, int order);
   void ResetReference();
   void ReadReference(string filename,int order);
   double GetReference(int Q2, int pt, int order);


   ClassDef(TIncljets,1)         
};
#endif
