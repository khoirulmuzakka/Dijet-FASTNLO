%module fastnlo

%{
/**
 * This is a workaround for a minor swig bug when building on gcc 4.6.1 and above.
 * Prior to gcc 4.6.1 the STL headers like vector, string, etc. used to 
 * automatically pull in the cstddef header but starting with gcc 4.6.1 they no
 * longer do. This leads to swig generated a file that does not compile so we
 * explicitly include cstddef so the swig generated file will compile.
 */
#include <cstddef>
%}

%include <std_string.i>
%include <std_vector.i>
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
%include "../fastnlotoolkit/include/fastnlotk/fastNLOConstants.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOBase.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOTable.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOReader.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOLHAPDF.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOAlphas.h"
%include "../fastnlotoolkit/include/fastnlotk/CRunDec.h"
%include "../fastnlotoolkit/include/fastnlotk/Alphas.h"
