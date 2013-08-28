#ifndef __fastNLOCoeffAddEvalBase__
#define __fastNLOCoeffAddEvalBase__

#include "fastNLOCoeffAddBase.h"
#include "fastNLOReader.h"

using namespace std;

class fastNLOCoeffAddEvalBase : public fastNLOCoeffAddBase  {
   friend class fastNLOReader;

public:
   fastNLOCoeffAddEvalBase(int NObsBin, fastNLOReader* reader);
   fastNLOCoeffAddEvalBase(const fastNLOCoeffBase& base , int iLOord, fastNLOReader* reader);

   //    void SetPDFCache( v5d pdfs ) { fPDFLC = pdfs; }
   //    v5d GetPDFCache() const { return fPDFLC; }
   //    void SetAlphasCache(v3d AlphasTwoPi ) { fAlphasTwoPi = AlphasTwoPi;}
   //    v3d GetAlphasCache() const { return fAlphasTwoPi; }

   virtual void CalcCrossSection() = 0;
   virtual v1d GetCrossSection() { return fCrossSection; };
   
   void SetfastNLOReader(fastNLOReader* reader) { fReader = reader;}
   fastNLOReader* GetfastNLOReader() const { return fReader; }

protected:
   fastNLOCoeffAddEvalBase();

public:
   
protected:
   v1d fCrossSection;
   fastNLOReader* fReader;						// reader to provide pdf and alpha_s interface functions

};

#endif
