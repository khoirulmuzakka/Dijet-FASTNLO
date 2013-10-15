%module fastnloreader
%include "std_string.i"
%include "std_vector.i"

namespace std {
   %template(vectori) vector<int>;
   %template(vectors) vector<string>;
   %template(vectord) vector<double>;
   %template(vectord2) vector<vector<double> >;
};

%{
#include "../fastnloreader/include/fastnlo/FastNLOReader.h"
#include "../fastnloreader/include/fastnlo/FastNLOUser.h"
#include "../fastnloreader/include/fastnlo/FastNLOLHAPDF.h"
#include "../fastnloreader/include/fastnlo/FastNLOAlphas.h"
#include "../fastnloreader/include/fastnlo/FastNLOCRunDec.h"
#include "../fastnloreader/include/fastnlo/Alphas.h"
%}
%include "../fastnloreader/include/fastnlo/speaker.h"
%include "../fastnloreader/include/fastnlo/FastNLOReader.h"
%include "../fastnloreader/include/fastnlo/FastNLOUser.h"
%include "../fastnloreader/include/fastnlo/FastNLOLHAPDF.h"
%include "../fastnloreader/include/fastnlo/FastNLOAlphas.h"
%include "../fastnloreader/include/fastnlo/CRunDec.h"
%include "../fastnloreader/include/fastnlo/Alphas.h"



