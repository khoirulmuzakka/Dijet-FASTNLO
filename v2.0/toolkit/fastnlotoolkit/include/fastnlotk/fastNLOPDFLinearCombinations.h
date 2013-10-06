#ifndef __fastNLOLinearCombinations__
#define __fastNLOLinearCombinations__

#include <string>
#include "speaker.h"
#include "fastNLOCoeffAddBase.h"

using namespace std;

class fastNLOPDFLinearCombinations {

public:
   fastNLOPDFLinearCombinations();
   ~fastNLOPDFLinearCombinations();

   vector<double > CalcPDFLinearCombination(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 = vector<double>(), const vector<double>& pdfx2 = vector<double>() , bool pdf2IsAntiParticle = false) const;

protected:
   vector<double > MakeAntiHadron(const vector<double >& hadron) const;

private: 
   vector<double > CalcPDFLCTwoHadrons(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1, const vector<double>& pdfx2 ) const ;
   vector<double > CalcPDFLCOneHadron(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 ) const;

   vector<double> CalcPDFLinearCombDIS(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1) const;
   vector<double> CalcPDFLinearCombHHC(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const ; // jets in hh
   vector<double> CalcPDFLinearCombttbar(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const ; // ttbar
};
#endif
