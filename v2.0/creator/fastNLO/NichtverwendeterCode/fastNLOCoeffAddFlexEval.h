#ifndef __fastNLOCoeffAddFlexEval__
#define __fastNLOCoeffAddFlexEval__

#include "fastNLOCoeffAddFlex.h"
#include "fastNLOCoeffAddEvalBase.h"
#include "fastNLOConstants.h"
#include "fastNLOReader.h"

using namespace std;

class fastNLOCoeffAddFlexEval : public fastNLOCoeffAddFlex , public fastNLOCoeffAddEvalBase {
   friend class fastNLOReader;

public:
   fastNLOCoeffAddFlexEval(int NObsBin, int iLOord);
   fastNLOCoeffAddFlexEval(const fastNLOCoeffBase& base , int iLOord);

   void SetPDFCache( v5d pdfs ) { fPDFLC = pdfs; }
   v5d GetPDFCache() const { return fPDFLC; }

   void SetAlphasCache(v3d AlphasTwoPi ) { fAlphasTwoPi = AlphasTwoPi;}
   v3d GetAlphasCache() const { return fAlphasTwoPi; }

   void CalcCrossSection();
   vector<double> GetCrossSection();
   
   void SetfastNLOReader(fastNLOReader* reader) { fReader = reader;}
   fastNLOReader* GetfastNLOReader() const { return fReader; }

protected:
   fastNLOCoeffAddFlexEval();

public:
   
protected:
   v5d fPDFLC;
   v3d fAlphasTwoPi;
   v1d fCrossSection;
   fastNLOReader* fReader;						// reader to provide pdf and alpha_s interface functions

};

#endif
