#ifndef __fastNLOCoeffAddFlex__
#define __fastNLOCoeffAddFlex__

#include "fastNLOCoeffAddBase.h"
#include "fastNLOConstants.h"


class fastNLOCoeffAddFlex : public fastNLOCoeffAddBase {

   friend class fastNLOTable;
   friend class fastNLOReader;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddFlex(int NObsBin, int iLOord);
   fastNLOCoeffAddFlex(const fastNLOCoeffBase& base , int iLOord);
   virtual ~fastNLOCoeffAddFlex(){;}
   virtual fastNLOCoeffBase* Clone() const;                                     //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false) ;
   virtual void Read(std::istream& table);
   void ReadRest(std::istream& table);
   virtual void Write(std::ostream& table);
   virtual void Print() const;
   virtual void Add(const fastNLOCoeffAddBase& other);
   virtual void Clear();                                                        //!< Clear all coefficients and event counters
   virtual void NormalizeCoefficients();                                        //!< Set number of events to 1 and normalize coefficients accordingly.
   virtual void MultiplyCoefficientsByConstant(double coef);                    //!< Multiply all coefficients by constant coef

   unsigned int GetNScaleNode1(int iObsBin) const { return ScaleNode1[iObsBin].size(); };
   unsigned int GetNScaleNode2(int iObsBin) const { return ScaleNode2[iObsBin].size(); };
   double GetScaleNode1(int iObsBin, int iNode) const { return ScaleNode1[iObsBin][iNode]; };
   double GetScaleNode2(int iObsBin, int iNode) const { return ScaleNode2[iObsBin][iNode]; };
   bool IsCompatible(const fastNLOCoeffAddFlex& other) const;                   //!< check for compatibilty for adding/merging of two tables

protected:

   fastNLOCoeffAddFlex();
   void ReadCoeffAddFlex(std::istream& table);

   int fILOord;   // obtained from Scenario

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
