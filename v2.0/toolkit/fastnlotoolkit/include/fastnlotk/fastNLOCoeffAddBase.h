#ifndef __fastNLOCoeffAddBase__
#define __fastNLOCoeffAddBase__


#include "fastNLOCoeffBase.h"
#include "fastNLOConstants.h"


class fastNLOCoeffAddBase : public fastNLOCoeffBase {

   friend class fastNLOTable;
   friend class fastNLOCreate;

public:
   fastNLOCoeffAddBase(int NObsBin);
   fastNLOCoeffAddBase(const fastNLOCoeffBase& base);
   virtual ~fastNLOCoeffAddBase() {;}
   virtual fastNLOCoeffBase* Clone() const;                                     //!< returns 'new' copy of this instance.
   static bool CheckCoeffConstants(const fastNLOCoeffBase* c, bool quiet = false);
   void Read(std::istream& table);
   virtual void Write(std::ostream& table);
   virtual void Add(const fastNLOCoeffAddBase& other);
   virtual void Print(int iprint) const;
   virtual void Clear();                                                        //!< Clear all coefficients and event counters
   virtual void NormalizeCoefficients();                                        //!< Set number of events to 1 and normalize coefficients accordingly.

   // Erase observable bin; iObsIdx is the C++ array index to be removed and
   // not the observable bin no. running from 1 to fNObsBins
   virtual void EraseBin(unsigned int iObsIdx);

   int GetIRef() const {return IRef;}
   double GetNevt() const { return Nevt; }
   double GetNevt(int NObsBin, int NSubproc) const {
      if (Nevt > 0) return Nevt;
      else {std::cout<<"Todo. Preparation for v2.3."<< std::endl; return Nevt;}
   }
   void SetNevt(double nevt) { Nevt = nevt;}                                    //!< Set number of events
   int GetNxmax(int Obsbin) const ;
   int GetXIndex(int Obsbin,int x1bin,int x2bin =0) const ;
   int GetNSubproc() const { return NSubproc;}
   int GetIScaleDep() const { return IScaleDep;}
   int GetNPDF() const {return NPDFPDG.size();}
   int GetPDFPDG(int iPDF) const {return NPDFPDG[iPDF];}
   int GetNPDFDim() const {return NPDFDim;}
   int GetIPDFdef1() const { return IPDFdef1; }
   int GetIPDFdef2() const { return IPDFdef2; }
   int GetIPDFdef3() const { return IPDFdef3; }
   int GetNpow() const {return Npow;}
   int GetNScales() const {return NScales;}
   int GetNScaleDim() const {return NScaleDim;}
   //std::vector<std::string > GetScaleDescript(int iScale=0) const { return ScaleDescript[iScale]; };
   std::string GetScaleDescription(int iScale=0) const { return ScaleDescript[0][iScale]; };         // getter for scale description of scale iScale
   std::vector<std::vector<std::string > > GetScaleDescr() const { return ScaleDescript; }
   int GetNxtot1(int iBin) const { return XNode1[iBin].size(); }
   int GetNxtot2(int iBin) const { return XNode2.size() > 0 ? XNode2[iBin].size() : -1; }

   double GetXNode1(int iObsBin, int iNode) const { return XNode1[iObsBin][iNode]; }
   double GetXNode2(int iObsBin, int iNode) const { return XNode2[iObsBin][iNode]; }
   double GetX1(int iObsBin, int iXnode) const; //! return x value of pdf1 for x-node 1
   double GetX2(int iObsBin, int iXnode) const; //! return x value of pdf1 for x-node 1

   std::vector < double > GetXNodes1(int iObsBin) const { return XNode1[iObsBin]; }
   std::vector < double > GetXNodes2(int iObsBin) const { return XNode2[iObsBin]; }

   bool IsReference() const {return IRef>0;};
   bool IsCompatible(const fastNLOCoeffAddBase& other) const;

   const std::vector<std::vector<std::pair<int,int> > >& GetPDFCoeff() const { return fPDFCoeff;}

protected:
   fastNLOCoeffAddBase();
   void ReadCoeffAddBase(std::istream& table);
   int GetScaledimfromvar(int scalevar) const;

   int IRef;
   int IScaleDep;
   double Nevt;
   int Npow;
   std::vector < int > NPDFPDG;
   int NPDFDim;
   std::vector < int > NFFPDG;
   int NFFDim;
   int NSubproc;
   int IPDFdef1;
   int IPDFdef2;
   int IPDFdef3;
   std::vector<std::vector<std::pair<int,int> > > fPDFCoeff;                                                   //! fPDFCoeff[iSubProc][iPartonPair][pair]
   // Missing: linear PDF combinations for IPDFdef1=0
   std::vector < double > Hxlim1;
   fastNLO::v2d XNode1;
   std::vector < double > Hxlim2;
   fastNLO::v2d XNode2;
   std::vector < int > Nztot;
   std::vector < double > Hzlim;
   fastNLO::v2d ZNode;
   int NScales;
   int NScaleDim;
   std::vector < int > Iscale;                                                                       // not used
   std::vector < std::vector < std::string > > ScaleDescript;

};

#endif
