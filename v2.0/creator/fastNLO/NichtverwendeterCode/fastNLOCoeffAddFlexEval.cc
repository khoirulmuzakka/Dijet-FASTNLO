#include <cstdlib>
#include <iostream>
#include <cmath>

#include "fastNLOCoeffAddFlexEval.h"


fastNLOCoeffAddFlexEval::fastNLOCoeffAddFlexEval() : fReader(NULL) {
   fastNLOCoeffAddFlex::SetClassName("fastNLOCoeffAddFlexEval");
}

 fastNLOCoeffAddFlexEval::fastNLOCoeffAddFlexEval(int NObsBin, int iLOord) : fastNLOCoeffAddFlex(NObsBin,iLOord) , fReader(NULL) {
    fastNLOCoeffAddFlex::SetClassName("fastNLOCoeffAddFlexEval");
 }

fastNLOCoeffAddFlexEval::fastNLOCoeffAddFlexEval(const fastNLOCoeffBase& base , int iLOord ) : fastNLOCoeffAddFlex(base,iLOord) , fReader(NULL) {
   fastNLOCoeffAddFlex::SetClassName("fastNLOCoeffAddFlexEval");
 }


//________________________________________________________________________________________________________________ //


void fastNLOCoeffAddFlexEval::CalcCrossSection() {
   
}


vector<double> fastNLOCoeffAddFlexEval::GetCrossSection() {
   if ( fCrossSection.empty() ) CalcCrossSection();
   return fCrossSection;
}
