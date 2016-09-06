#ifndef __fastNLOCoeffAddFix__
#define __fastNLOCoeffAddFix__

#include "fastNLOCoeffAddBase.h"
#include "fastNLOConstants.h"


class fastNLOCoeffAddFix : public fastNLOCoeffAddBase {

   friend class fastNLOReader;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddFix(int NObsBin);
   fastNLOCoeffAddFix(const fastNLOCoeffBase&);
   virtual ~fastNLOCoeffAddFix(){;}
   virtual fastNLOCoeffBase* Clone() const;                                     //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   virtual void Read(std::istream&table);
   void ReadRest(std::istream& table);
   virtual void Write(std::ostream& table);
   virtual void Add(const fastNLOCoeffAddBase& other);
   virtual void Print(int iprint) const;
   virtual void Clear();                                                        //!< Clear all coefficients and event counters
   virtual void NormalizeCoefficients();                                        //!< Set number of events to 1 and normalize coefficients accordingly.
   virtual void MultiplyCoefficientsByConstant(double coef);                    //!< Multiply all coefficients by constant coef

   // Erase observable bin; iObsIdx is the C++ array index to be removed and
   // not the observable bin no. running from 1 to fNObsBins
   void EraseBin(unsigned int iObsIdx);

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

protected:
   fastNLOCoeffAddFix();
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
