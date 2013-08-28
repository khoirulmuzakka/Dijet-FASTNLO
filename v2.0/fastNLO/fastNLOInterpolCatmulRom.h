// Author: Daniel Britzger
// DESY, 28/06/2013

#ifndef __fastNLOInterpolCatmulRom__
#define __fastNLOInterpolCatmulRom__



#include "speaker.h"
#include <string>
#include <vector>
#include <utility>
#include "fastNLOInterpolBase.h"

namespace fastNLO {
}

using namespace std;
using namespace fastNLO;

class fastNLOInterpolCatmulRom : public fastNLOInterpolBase {

public:

   fastNLOInterpolCatmulRom(double min, double max);
   ~fastNLOInterpolCatmulRom(void);
   
   vector<pair<int,double> > CalcNodeValues(double val);

protected:


private:


};


#endif