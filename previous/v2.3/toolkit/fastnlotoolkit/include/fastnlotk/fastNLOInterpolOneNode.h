// Author: Daniel Britzger
// DESY, 29/01/2014

/**
   fastNLOInterpolOneNode

   No interpolation is performed, but all
   values are stored at one single node.
   The 'x'-value of the node is calculated
   from the H-function.
*/

#ifndef __fastNLOInterpolOneNode__
#define __fastNLOInterpolOneNode__

#include "speaker.h"
#include <string>
#include <vector>
#include <utility>
#include "fastNLOInterpolBase.h"


class fastNLOInterpolOneNode :  public fastNLOInterpolBase {
   
public:

   fastNLOInterpolOneNode(double min, double max, fastNLOGrid::GridType type);
   ~fastNLOInterpolOneNode(void);
   
   //   vector<pair<int,double> > CalcNodeValues(double val);
   void CalcNodeValues(std::vector<std::pair<int,double> >& nodes, double val);

protected:


private:
   std::vector<std::pair<int,double> > fDummyNode;

};


#endif
