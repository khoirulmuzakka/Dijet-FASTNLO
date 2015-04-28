//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO_reader_2.1.0                                                //
//  D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch           //
//                                                                      //
//  The projects web page can be found at:                              //
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//  If you use this code, please cite:                                  //
//    T. Kluge, K. Rabbertz and M. Wobisch, hep-ph/0609285              //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch,        //
//       arXiv:1109.1310                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//
//  fastNLOExtern
//  This class allows to supply the three main abstract methods from an
//  externally linked library.
//
//////////////////////////////////////////////////////////////////////////

#include "fastNLOReader.h"

using namespace std;

// Fortran Functions that get consumed by C++
extern "C" {
   void fastnlo_initpdf_();
   void fastnlo_getxfx_(double *ret, double* x, double* muf, int* f);
   void fastnlo_evolve_as_(double *ret, double* q);
}

class fastNLOExtern : public fastNLOReader {
public:
   fastNLOExtern(std::string tablename) : fastNLOReader(tablename) {}
protected:
   // Call external function to init PDF
   virtual bool InitPDF() {
      fastnlo_initpdf_();
      return true;
   }

   // Call external function to get PDF values
   virtual vector<double> GetXFX(double x, double muf) const {
      vector <double> xfx(13);
      std::cout << muf << std::endl;
      for (int f = -6; f < 7; ++f)
         fastnlo_getxfx_(&(xfx[f + 6]), &x, &muf, &f);
      return xfx;
   }

   // Call external function to get alpha_s value
   virtual double EvolveAlphas(double Q) const {
      double result;
      fastnlo_evolve_as_(&result, &Q);
      return result;
   }
};

std::map<int, fastNLOExtern*> fastNLO_context;

void v2f(vector<double> v, double *result) {
   for (size_t i = 0; i < v.size(); ++i) // memcpy could be dangerous?
      result[i] = v[i];
}

// C++ Functions that get consumed by Fortran
extern "C" {
   void fastnlo_create_(int *ctx, char *tbl_name, int tbl_name_len) {
      fastNLOExtern *reader = new fastNLOExtern(std::string(tbl_name, tbl_name_len));
      *ctx = 0;
      if (fastNLO_context.rbegin() != fastNLO_context.rend())
         *ctx = fastNLO_context.rbegin()->first + 1;
      fastNLO_context[*ctx] = reader;
   }

   void fastnlo_destroy_(int *ctx) {
      fastNLOExtern *reader = fastNLO_context[*ctx];
      delete reader;
      fastNLO_context.erase(*ctx);
   }

   void fastnlo_setscalefactorsmurmuf_(int *ctx, double *muR, double *muF) {
      fastNLO_context[*ctx]->SetScaleFactorsMuRMuF(*muR, *muF);
   }

   void fastnlo_getcrosssection_(int *ctx, double *result) {
      v2f(fastNLO_context[*ctx]->GetCrossSection(), result);
   }

   void fastnlo_getreferencecrosssection_(int *ctx, double *result) {
      v2f(fastNLO_context[*ctx]->GetReferenceCrossSection(), result);
   }

   void fastnlo_getkfactors_(int *ctx, double *result) {
      v2f(fastNLO_context[*ctx]->GetKFactors(), result);
   }

   void fastnlo_getqscales_(int *ctx, int *irelord, double *result) {
      v2f(fastNLO_context[*ctx]->GetQScales(*irelord), result);
   }
}
