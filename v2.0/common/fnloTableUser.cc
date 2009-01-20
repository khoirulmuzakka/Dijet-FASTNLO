#include "fnloTableUser.h"

#include "LHAPDF.h"

int fnloTableUser::ReadPDF(string filename){
   
   lhapdffilename = filename;
   LHAPDF::initPDFSetByName(lhapdffilename);
   lastusedpdf=1;
   nPDF = LHAPDF::numberPDF();

   return 0;
}

int fnloTableUser::ReadPDF2(string filename){
   
   lhapdffilename2 = filename;
   LHAPDF::initPDFSetByName(lhapdffilename2);
   lastusedpdf=2;
   nPDF2 = LHAPDF::numberPDF();

   return 0;
}

void fnloTableUser::SetPDFSet(int set){
   pdfset = set;
   LHAPDF::initPDF(pdfset);
}

void fnloTableUser::SetPDFSet2(int set){
   pdfset2 = set;
   LHAPDF::initPDF(pdfset2);
}

void fnloTableUser::GetPdfs(double x, double muf, vector<double> &xfx){
   if(lastusedpdf!=1){
      LHAPDF::setVerbosity(LHAPDF::SILENT);
      LHAPDF::initPDFSetByName(lhapdffilename);
      LHAPDF::initPDF(pdfset);
      lastusedpdf=1;
   }
   xfx = LHAPDF::xfx(x,muf);
}

void fnloTableUser::GetPdfs2(double x, double muf, vector<double> &xfx){
   if(lastusedpdf!=2){
      LHAPDF::setVerbosity(LHAPDF::SILENT);
      LHAPDF::initPDFSetByName(lhapdffilename2);
      LHAPDF::initPDF(pdfset2);
      lastusedpdf=2;
   }
   xfx = LHAPDF::xfxp(x,muf,0,0);
}

void fnloTableUser::FillPDFCache(){
   if(BlockIndexLO>=0){
      if(GetBlockB(BlockIndexLO)->NPDFDim<2){
         GetBlockB(BlockIndexLO)->FillPDFCache(scalevar, &fnloTableUser::GetPdfs,this);
      }else{
         GetBlockB(BlockIndexLO)->FillPDFCache(scalevar, &fnloTableUser::GetPdfs, &fnloTableUser::GetPdfs2,this);
      }
   }
   if(BlockIndexNLO>=0){
      if(GetBlockB(BlockIndexNLO)->NPDFDim<2){
         GetBlockB(BlockIndexNLO)->FillPDFCache(scalevar, &fnloTableUser::GetPdfs,this);
      }else{
         GetBlockB(BlockIndexNLO)->FillPDFCache(scalevar, &fnloTableUser::GetPdfs,&fnloTableUser::GetPdfs2,this);
      }
   }
}

void fnloTableUser::CalcXsection(){
   if(BlockIndexLO>=0) GetBlockB(BlockIndexLO)->CalcXsection(alphasmz,scalevar,rescale);
   if(BlockIndexNLO>=0) GetBlockB(BlockIndexNLO)->CalcXsection(alphasmz,scalevar,rescale);
   if(BlockIndexLORef>=0) GetBlockB(BlockIndexLORef)->CalcXsection(alphasmz,scalevar,rescale);
   if(BlockIndexNLORef>=0) GetBlockB(BlockIndexNLORef)->CalcXsection(alphasmz,scalevar,rescale);
}

double fnloTableUser::GetXsection(int bin1, int order){
   double xsection = 0;
   int bin = bin1;
   if(order>0 && BlockIndexLO>=0){ // LO
      xsection += GetBlockB(BlockIndexLO)->GetXsection(bin);
   }
   if(order>1 && BlockIndexNLO>=0){ // NLO
      xsection += GetBlockB(BlockIndexNLO)->GetXsection(bin);
   }
   return xsection;
}

double fnloTableUser::GetXsectionRef(int bin1, int order){
   double xsection = 0;
   int bin = bin1;
   if(order>0 && BlockIndexLORef>=0){ // LO
      xsection += GetBlockB(BlockIndexLORef)->GetXsection(bin);
   }
   if(order>1 && BlockIndexNLO>=0){ // NLO
      xsection += GetBlockB(BlockIndexNLORef)->GetXsection(bin);
   }
   return xsection;
}

void fnloTableUser::SetMurScale(int scale){
   scalevar1 = scale;
}

void fnloTableUser::SetMufScale(int scale){
   scalevar2 = scale;
}
