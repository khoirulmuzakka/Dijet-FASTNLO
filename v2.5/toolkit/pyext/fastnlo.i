// To be tested: Following suggestion from Dmitry Kalinkin, 12.06.2021
%module(directors="1") fastnlo

// generate directors for all classes that have virtual methods
%feature("director");

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
   %template(pairi) pair<int,int>;
   %template(vectorpairi) vector<pair<int,int> >;
   %template(vector2pairi) vector<vector<pair<int,int> > >;
};

%ignore *::operator[];
%ignore *operator>>;
%rename("PrintMessage") speaker::print;


/* Suppress SWIG warning */
#pragma SWIG nowarn=SWIGWARN_PARSE_NESTED_CLASS
/* Redefine nested class in global scope in order for SWIG to generate */
/* a proxy class. Only SWIG parses this definition. */

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
%include "../fastnlotoolkit/include/fastnlotk/fastNLOTable.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOReader.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOLHAPDF.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOAlphas.h"
%include "../fastnlotoolkit/include/fastnlotk/fastNLOCRunDec.h"
%include "../fastnlotoolkit/include/fastnlotk/Alphas.h"
/* Avoid syntax errors because of 'as' by ignoring declarations using '*as' in CRunDec */
/* Let's hope this doesn't break anything important ... */
%ignore mPS2mSI;
%ignore m1S2mSI;
%ignore mRS2mSI;
%ignore mRSp2mSI;
%include "../fastnlotoolkit/include/fastnlotk/CRunDec.h"
