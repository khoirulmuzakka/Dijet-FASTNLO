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
      kSqrtLog10        = 3,            // sqrt(logarithmic) grid
      kLogLog           = 4,            // loglog grid (not applicable to x)
      k4thrtLog10       = 5,            // log(x)^(1/4) (4th root)
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
   const std::vector<double>& GetGrid() const { return fgrid;}
   const std::vector<double>* GetGridPtr() const { return &fgrid;}
   const std::vector<double>& GetHGrid() const { return fHgrid;}
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

   inline double Function_loglog025( double mu ){
      // function H(mu) = log(log( mu / 0.25 ))
      return log(log(mu/0.25));}
   inline double Function_loglog025_inv( double mu ){
      // inverse of function H(mu) = log(log( mu / 0.25 ))
      return 0.25*exp(exp(mu));}
   inline double Function_loglog( double mu ){
      // function H(mu) = log(log( mu ))
      return log(log(mu));}
   inline double Function_loglog_inv( double mu ){
      // inverse of function H(mu) = log(log( mu ))
      return exp(exp(mu));}
   inline double Function_x( double mu ){
      // function H(mu) = x
      return mu;}
   inline double Function_x_inv( double mu ){
      // inverse of function H(mu) = x;
      return mu;}
   inline double Function_log10( double x ){return log10(x);}
   inline double Function_log10_inv( double x ){return pow(10,x);}
   inline double Function_sqrtlog10( double x ){return -sqrt(-log10(x));}
   inline double Function_sqrtlog10_inv( double x ){return pow(10,-pow(x,2));}
   inline double Function_4thrtlog10( double mu ){return -pow(-log10(mu),0.25);}
   inline double Function_4thrtlog10_inv( double mu ){return pow(10,-pow(mu,4));}

   std::vector<double> HGrid_loglog025_inv(std::vector<double> grid);
   std::vector<double> HGrid_loglog_inv(std::vector<double> grid);
   std::vector<double> HGrid_log10_inv(std::vector<double> grid);
   std::vector<double> HGrid_sqrtlog10_inv(std::vector<double> grid);
   std::vector<double> HGrid_4thrtlog10_inv(std::vector<double> grid);
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
