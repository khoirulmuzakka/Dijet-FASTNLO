// Author: Daniel Britzger
// DESY, 29/01/2014

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include "fastnlotk/fastNLOInterpolLinear.h"

using namespace std;

//______________________________________________________________________________
fastNLOInterpolLinear::fastNLOInterpolLinear(double min, double max) : fastNLOInterpolBase(min,max,2) {
   debug["fastNLOInterpolLinear"]<<"New fastNLOInterpolLinear instance."<<endl;
}


//______________________________________________________________________________
fastNLOInterpolLinear::~fastNLOInterpolLinear(void) {
}


//______________________________________________________________________________
void fastNLOInterpolLinear::CalcNodeValues(vector<pair<int,double> >& nodes, double x) {
   //! Performs interpolation of value value on grid 'fgrid'.
   //! uses distance measure 'fdm'
   //! returns for for all relevant grid points
   //!  -  the integer number of that node
   //!  -  the value, which this node obtains.
   //!

   // --- relative distance delta - in function fdm H(x)
   // deltascale (Interpol(.,.,.delta,.): relative distance of value to node 'nnode'
   double delta = GetDelta(x);
   // --- get scale interpolation kernel and updated scalenode position: 1 <= nmu < ntot-2
   int nnode = FindLargestPossibleNode(x);
   // --- set nodes
   nodes.resize(2);
   nodes[0] = make_pair(nnode,  1.-delta);
   nodes[1] = make_pair(nnode+1,delta);

   // std::cout<<" *** node: x="<<x<<"\tdelta="<<delta<<"\tnnode="<<nnode
   // 	    <<"\tn0="<<nodes[0].first<<"\tv0="<<nodes[0].second
   // 	    <<"\tn1="<<nodes[1].first<<"\tv0="<<nodes[1].second
   // 	    <<endl;
	    
   if (fLastGridPointWasRemoved ) {
      if ( nodes.back().first==(int)fgrid.size() ) {
	 nodes.resize(1);
	 if ( nodes.back().first==(int)fgrid.size() ) {
	    nodes.resize(0);
	 }
      }
   }
}

