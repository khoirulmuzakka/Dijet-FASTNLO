%module fastnlo
%include "std_string.i"
%include "std_vector.i"
%include <std_pair.i>

namespace std {
   %template(vectori) vector<int>;
   %template(vectors) vector<string>;
   %template(vectord) vector<double>;
   %template(vectord2) vector<vector<double> >;
   %template() pair<double,double>;
   %template(pairvector) vector<pair<double,double> >;
};

%ignore *::operator[];
%ignore *operator>>;
%rename("PrintMessage") speaker::print;

%{
#include "../fastnlotoolkit/include/fastnlotk/fastNLOTable.h"
#include "../fastnlotoolkit/include/fastnlotk/fastNLOReader.h"
#include "../fastnlotoolkit/include/fastnlotk/fastNLOLHAPDF.h"
#include "../fastnlotoolkit/include/fastnlotk/fastNLOAlphas.h"
#include "../fastnlotoolkit/include/fastnlotk/fastNLOCRunDec.h"
#include "../fastnlotoolkit/include/fastnlotk/Alphas.h"
%}
%include "../fastnlotoolkit/include/fastnlotk/speaker.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOPDFLinearCombinations.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOBase.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOTable.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOReader.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOLHAPDF.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOAlphas.h"
%include "../fastnlotoolkit/include/fastnlotk/CRunDec.h"
%include "../fastnlotoolkit/include/fastnlotk/Alphas.h"
