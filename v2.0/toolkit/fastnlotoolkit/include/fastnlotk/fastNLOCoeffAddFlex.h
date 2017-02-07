#ifndef __fastNLOCoeffAddFlex__
#define __fastNLOCoeffAddFlex__

#include "fastNLOCoeffAddBase.h"
#include "fastNLOConstants.h"


class fastNLOCoeffAddFlex : public fastNLOCoeffAddBase {

   friend class fastNLOTable;
   friend class fastNLOReader;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddFlex() = delete;
   fastNLOCoeffAddFlex(int NObsBin, int iLOord);
   explicit fastNLOCoeffAddFlex(const fastNLOCoeffBase& base , int iLOord);
   virtual ~fastNLOCoeffAddFlex(){;}
   virtual fastNLOCoeffAddFlex* Clone() const;                                     //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false) ;
   virtual void Read(std::istream& table);
   void ReadRest(std::istream& table);
   virtual void Write(std::ostream& table);
   virtual void Print(int iprint) const;
   virtual void Add(const fastNLOCoeffAddBase& other);

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
   virtual void CatBin(const fastNLOCoeffAddFlex& other, unsigned int iObsIdx);

   unsigned int GetNScaleNode1(int iObsBin) const { return ScaleNode1[iObsBin].size(); };
   unsigned int GetNScaleNode2(int iObsBin) const { return ScaleNode2[iObsBin].size(); };
   double GetScaleNode1(int iObsBin, int iNode) const { return ScaleNode1[iObsBin][iNode]; };
   double GetScaleNode2(int iObsBin, int iNode) const { return ScaleNode2[iObsBin][iNode]; };
   std::vector < double > GetScaleNodes1(int iObsBin) const { return ScaleNode1[iObsBin]; };
   std::vector < double > GetScaleNodes2(int iObsBin) const { return ScaleNode2[iObsBin]; };
   bool IsCompatible(const fastNLOCoeffAddFlex& other) const;                   //!< check for compatibilty for adding/merging of two tables
   bool IsCatenable(const fastNLOCoeffAddFlex& other) const;        //!< Check for compatibility of two contributions for merging/adding

protected:

   void ReadCoeffAddFlex(std::istream& table);

   int fILOord;   // obtained from Scenario
   int fSTildeDISFormat = 1; // format of sigma-tilde coefficients (0: log(mu2/q2), 1: log(mu2))

   // SigmaTilde [NObsBins] ['n' x-nodes] [n s1-Nodes] [n s2-Nodes] [nsubproc]
   fastNLO::v5d SigmaTildeMuIndep; // units are (p)barn * Nevt / BinSize
   fastNLO::v5d SigmaTildeMuFDep;
   fastNLO::v5d SigmaTildeMuRDep;
   fastNLO::v5d SigmaTildeMuRRDep;
   fastNLO::v5d SigmaTildeMuFFDep;
   fastNLO::v5d SigmaTildeMuRFDep;
   // SigmaRef [NObsBins] [nsubproc]
   fastNLO::v2d SigmaRefMixed;  // units are (p)barn * Nevt / BinSize
   fastNLO::v2d SigmaRef_s1;
   fastNLO::v2d SigmaRef_s2;
   //int NscalenodeScale1;
   //int NscalenodeScale2;
   // ScaleNodeXY [ObsBin] [NscalenodeScaleX]
   fastNLO::v2d ScaleNode1;
   fastNLO::v2d ScaleNode2;

public:
   fastNLO::v3d AlphasTwoPi;
   fastNLO::v5d PdfLcMuVar;

};

#endif
