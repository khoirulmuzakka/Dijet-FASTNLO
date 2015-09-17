// Author: Daniel Britzger
// DESY, 28/06/2013

#ifndef __fastNLOInterpolBase__
#define __fastNLOInterpolBase__


#include "speaker.h"
#include <string>
#include <vector>
#include <cmath>
#include <utility>

namespace fastNLOGrid {
   enum GridType {
      kLinear           = 0,            // linear grid
      kLogLog025        = 1,            // loglog grid
      kLog10            = 2,            // log10 grid
      kSqrtLog10        = 3             // sqrt(logarithmic) grid
   };
}


class fastNLOInterpolBase : public PrimalScream {

public:

   fastNLOInterpolBase(double min, double max, int nMinNodes);
   virtual ~fastNLOInterpolBase(void);

   const std::vector<std::pair<int,double> >& GetNodeValues(double val);
   std::vector<std::pair<int,double> >* GetNodeValuesPtr(double val);

   void MakeGrids(fastNLOGrid::GridType type, int nNodes);
   void MakeGridsWithNNodesPerMagnitude(fastNLOGrid::GridType type, int nNodes);
   void RemoveLastNode();

   void PrintGrid();
   std::vector<double> GetGrid() const { return fgrid;}
   const std::vector<double>* GetGridPtr() const { return &fgrid;}
   std::vector<double> GetHGrid() const { return fHgrid;}
   double GetDelta(double);
   bool CheckX(double&);

   static fastNLOGrid::GridType TranslateGridType(std::string in);

protected:

   void SetGrid(std::vector<double> grid);
   void SetHGrid(std::vector<double> grid);
   void MakeGrids(double min, double max, int nNodes);
   std::vector<double> MakeGridFromHGrid(std::vector<double> g);
   std::vector<double> MakeLinearGrid(double min, double max, int nNodes);

   //virtual std::vector<std::pair<int,double> > CalcNodeValues(double val) = 0;
   virtual void CalcNodeValues(std::vector<std::pair<int,double> >& nodes, double val) = 0;

   int FindLargestPossibleNode(double);

   double Function_loglog025( double mu ){
      // function H(mu) = log(log( mu / 0.25 ))
      return log(log(mu/0.25));}
   double Function_loglog025_inv( double mu ){
      // inverse of function H(mu) = log(log( mu / 0.25 ))
      return 0.25*exp(exp(mu));}
   double Function_x( double mu ){
      // function H(mu) = x
      return mu;}
   double Function_x_inv( double mu ){
      // inverse of function H(mu) = x;
      return mu;}
   double Function_log10( double x ){return log10(x);}
   double Function_log10_inv( double x ){return pow(10,x);}
   double Function_sqrtlog10( double x ){return -sqrt(-log10(x));}
   double Function_sqrtlog10_inv( double x ){return pow(10,-pow(x,2));}

   std::vector<double> HGrid_loglog025_inv(std::vector<double> grid);
   std::vector<double> HGrid_log10_inv(std::vector<double> grid);
   std::vector<double> HGrid_sqrtlog10_inv(std::vector<double> grid);
   int GetNMod() const {return fnmod;}
   double GetHx(double);

protected:
   std::vector<std::pair<int,double> > fNodes;

   int fNMinNodes;
   double fvalmin;
   double fvalmax;
   double fLastVal;
   bool fLastGridPointWasRemoved; // odd boolean to agree with original code;
   fastNLOGrid::GridType fdm; // distance measure
   std::vector<double> fgrid;
   std::vector<double> fHgrid;
   int fnmod ; // variable for final nodes. Has to be filled by inherited algorithm

};

#endif
