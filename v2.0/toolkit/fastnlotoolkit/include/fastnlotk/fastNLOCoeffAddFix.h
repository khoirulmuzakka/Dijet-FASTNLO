#ifndef __fastNLOCoeffAddFix__
#define __fastNLOCoeffAddFix__

#include "fastNLOCoeffAddBase.h"
#include "fastNLOConstants.h"


class fastNLOCoeffAddFix : public fastNLOCoeffAddBase {

   friend class fastNLOReader;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddFix() = delete;
   fastNLOCoeffAddFix(int NObsBin);
   fastNLOCoeffAddFix(const fastNLOCoeffBase&);
   virtual ~fastNLOCoeffAddFix(){;}
   virtual fastNLOCoeffAddFix* Clone() const;                                     //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   virtual void Read(std::istream&table);
   void ReadRest(std::istream& table);
   virtual void Write(std::ostream& table);
   virtual void Add(const fastNLOCoeffAddBase& other);
   virtual void Print(int iprint) const;

   // Manipulate coefficient bins
   // Clear all coefficients and event counters
   virtual void Clear();
   // Set number of events to unity and normalize coefficients accordingly
   virtual void NormalizeCoefficients();
   // Multiply all coefficients of all bins by a constant factor
   virtual void MultiplyCoefficientsByConstant(double fact);
   // In the following, iObsIdx is the C++ array index of the concerned bin and
   // not the observable bin no. running from 1 to fNObsBins!
   // Multiply coefficients of one observable bin a factor
   virtual void MultiplyBin(unsigned int iObsIdx, double fact);
   // Erase observable bin from table
   virtual void EraseBin(unsigned int iObsIdx);
   // Catenate observable to table
   virtual void CatBin(const fastNLOCoeffAddFix& other, unsigned int iObsIdx);

   int GetTotalScalevars() const ;
   int GetTotalScalenodes() const ;
   int GetNScaleNode() const { return GetTotalScalenodes(); }
   int GetNScalevar() const { return Nscalevar[0];}
   fastNLO::v1d GetAvailableScaleFactors() const { return ScaleFac[0]; }
   double GetScaleFactor(int iVar) const {
      if ( iVar >= (int)ScaleFac[0].size() )
         this->error["GetScaleFactor"]<<"Scalevariation no. "<<iVar<<" not available. There are only "<<GetNScalevar()<<" available in this table."<< std::endl;
      return ScaleFac[0][iVar];
   }

   double GetSigmaTilde(int iObs, int iSvar, int ix, int is, int iN ) const { return SigmaTilde[iObs][iSvar][ix][is][iN];}
   double GetScaleNode(int iObs, int iSvar, int iNode ) const { return ScaleNode[iObs][0][iSvar][iNode]; }
   std::vector < double > GetScaleNodes(int iObs, int iSvar) const { return ScaleNode[iObs][0][iSvar]; }

   void ResizePdfLC();
   void ResizePdfSplLC();
   void ResizeSigmaTilde();
   bool IsCompatible(const fastNLOCoeffAddFix& other) const;                   //!< Check for compatibility of two contributions for merging/adding
   bool IsCatenable(const fastNLOCoeffAddFix& other) const;        //!< Check for compatibility of two contributions for merging/adding

protected:
   void ReadCoeffAddFix(std::istream& table);

   std::vector < int > Nscalevar;
   //std::vector < int > Nscalenode;
   fastNLO::v2d ScaleFac;
   fastNLO::v4d ScaleNode;
   fastNLO::v5d SigmaTilde; // units are (p)barn * Nevt / BinSize

public:
   fastNLO::v2d AlphasTwoPi_v20;
   fastNLO::v4d PdfLc;
   fastNLO::v4d PdfSplLc1;
   fastNLO::v4d PdfSplLc2;
};

#endif
