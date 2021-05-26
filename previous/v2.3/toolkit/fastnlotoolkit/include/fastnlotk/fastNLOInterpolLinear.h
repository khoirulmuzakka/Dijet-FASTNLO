// Author: Daniel Britzger
// DESY, 29/01/2014

/**
   fastNLOInterpolLinear

   Interpolation routines for linear interpolation.
*/

#ifndef __fastNLOInterpolLinear__
#define __fastNLOInterpolLinear__

#include "speaker.h"
#include <string>
#include <vector>
#include <utility>
#include "fastNLOInterpolBase.h"


class fastNLOInterpolLinear : public fastNLOInterpolBase {
   
public:

   fastNLOInterpolLinear(double min, double max, fastNLOGrid::GridType type);
   ~fastNLOInterpolLinear(void);
   
   //   vector<pair<int,double> > CalcNodeValues(double val);
   void CalcNodeValues(std::vector<std::pair<int,double> >& nodes, double val);

protected:


private:


};


#endif
