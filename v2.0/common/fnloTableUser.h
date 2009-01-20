#ifndef __fnloTableUser__
#define __fnloTableUser__


#include "fnloTable.h"


class fnloTableUser : public fnloTable {
public:
   fnloTableUser(string name) : fnloTable(name){alphasmz=0.118;scalevar=0;rescale=1.0;lastusedpdf=-1;}
   int ReadPDF(string filename);
   int ReadPDF2(string filename);
   void SetPDFSet(int set);
   void SetPDFSet2(int set);
   void SetAlphasMZ(double asmz){alphasmz=asmz;}
   void SetRescale(double xmur){rescale=xmur;}
   void SetScalevar(int scale){scalevar = scale;}
   void SetMurScale(int scale);
   void SetMufScale(int scale);
   void FillPDFCache();
   void CalcXsection();
   double GetXsection(int bin1, int order);
   double GetXsectionRef(int bin1, int order);
   double GetSmallestX(int blockb, int ObsBin){return GetBlockB(blockb)->GetSmallestX(ObsBin);}
   double GetSmallestX2(int blockb, int ObsBin){return GetBlockB(blockb)->GetSmallestX2(ObsBin);}


private:
   void GetPdfs(double x, double muf, vector<double> &xfx);
   void GetPdfs2(double x, double muf, vector<double> &xfx);
   string lhapdffilename; // First hadron
   string lhapdffilename2; // Second hadron (for gammaP)
   int lastusedpdf; // when switching between lhapdf and lhapdf2 one needs to reinitialise... unfortunately
   
   int nPDF; // total # of pdfs within series
   int nPDF2;

   int pdfset; // set within pdf series to use 
   int pdfset2;  

   double alphasmz;
   int scalevar1; // used for mur, and for muf if only one scale
   int scalevar2; // used for muf
   int scalevar;  // computed from the scalevar1 and scalevar2

   double rescale;



};


#endif
