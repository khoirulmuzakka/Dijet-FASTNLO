#include "TIncljets.h"

ClassImp(TIncljets)

TIncljets::TIncljets(){
   q2binning[0] = 150.0;
   q2binning[1] = 200.0;
   q2binning[2] = 300.0;
   q2binning[3] = 600.0;
   q2binning[4] = 5000.0;

   ptbinning[0] = 7.0;
   ptbinning[1] = 11.0;
   ptbinning[2] = 18.0;
   ptbinning[3] = 30.0;
   ptbinning[4] = 50.0;

}

double TIncljets::alphas(double Q, double alphasMZ){
   //Calculating alpha_s at mu_r, given an alpha_s(m_Z) using series approach
   //See Ellis,Sterling,Webber, p26

   // Constants NF,MZ,TWOPI,TWOPISQR are set in physconst.h
   return alphasMZ / (1. + alphasMZ*
                      (  (BETA0/TWOPI) 
                         +(BETA1/TWOPISQR)*alphasMZ  )
                      * log(Q/MZ));
}

void TIncljets::ReadTable(string filename){
   // #define READ(n) table.read(reinterpret_cast<char *>(&n),sizeof(n))
#define READ(n) table >> n

   fstream table(filename.c_str(),ios::in); // open file
   int marker;

   READ(ireaction);
   switch(ireaction){
   case kdis:
      nsubproc = 3;
      break;
   default:
      printf("Ireaction = %d not supported. Stopping.\n",ireaction);      
      exit(1);
      break;
   }
   
   READ(s);
   READ(iproc);
   READ(ialgo);
   READ(JetResol1);
   READ(JetResol2);
   READ(npow);
   READ(Oalphas);
   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after Oalphas in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   nevents.resize(2);
   READ(nevents[0]);
   printf("Reading %.2G LO events.\n",nevents[0]);
   READ(nevents[1]);
   printf("Reading %.2G NLO events.\n",nevents[1]);
   READ(ntot);
   READ(ixscheme);
   READ(ipdfwgt);
   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after ipdfwght in %s. Stopping.\n",filename.c_str());
      exit(1);
   }
   READ(nQ2);
   Q2high = new double[nQ2+1];      //----- array for Q2 boundaries
   for(int i=0;i<nQ2+1;i++){
      READ(Q2high[i]); //Q2[0] ... Q2[NQ2]
   }

   npt     = new  int[nQ2];         // nrap bins in Q2 with each npt[irap] bins in pt
   for(int i=0;i<nQ2;i++){
      READ(npt[i]); //Npt[0] ... Npt[NQ2-1] 
   }

   pthigh.resize(nQ2);
   for(int i=0;i<nQ2;i++){
      pthigh[i].resize(npt[i]+1);
      for(int j=0;j<npt[i]+1;j++){
         READ(pthigh[i][j]); // all pt bins 
      }
   }
    
   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after pthigh in %s. Stopping.\n",filename.c_str());
      exit(1);
   }
   
   xlimit.resize (nQ2);
   for(int i=0;i<nQ2;i++){
      xlimit[i].resize(npt[i]);
      for(int j=0;j<npt[i];j++){
         READ(xlimit[i][j]); // all x limits 
      }
   }

   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after xlimit in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   murval.resize (nQ2);
   for(int i=0;i<nQ2;i++){
      murval[i].resize (npt[i]);
      for(int j=0;j<npt[i];j++){
         READ(murval[i][j]); // all murval 
      }
   }
   
   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after murval in %s. Stopping.\n",filename.c_str());
      exit(1);
   }


   mufval.resize (nQ2);
   for(int i=0;i<nQ2;i++){
      mufval[i].resize (npt[i]);
      for(int j=0;j<npt[i];j++){
         READ(mufval[i][j]); // all mufval 
      }
   }

   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after mufval in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   //scale variations for NLO
   READ(nscalevar);
   murscale.resize(nscalevar);
   mufscale.resize(nscalevar);
   for(int i=0;i<nscalevar;i++){
      READ(murscale[i]);
   }
   for(int i=0;i<nscalevar;i++){
      READ(mufscale[i]);
   }
   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after mufscale in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   if(ireaction == kdis){
      disweights.resize (nQ2);
      pdf.resize (nQ2);
      for(int i=0;i<nQ2;i++){ // Q2
         disweights[i].resize(npt[i]);
         pdf[i].resize(npt[i]);
         for(int j=0;j<npt[i];j++){ // pt
            disweights[i][j].resize(ntot);
            pdf[i][j].resize(ntot);
            for(int k=0;k<ntot;k++){ // x
               disweights[i][j][k].resize(1+nscalevar);
               for(int m=0;m<nsubproc;m++){     //subprocesses 
                  for(int scalevar=0; scalevar<1+nscalevar;scalevar++){ //LO and NLO with scale variations)
                     READ(disweights[i][j][k][scalevar][m]);
                  }
               }
            }
         }
      }
   }else{
      printf("Ireaction = %d not supported. Stopping.\n",ireaction);      
      exit(1);
   }
      
   marker=0;
   READ(marker);// ------------------END of table
   if(marker!=cmarker){
      printf("Misaligned table after weights in %s. Stopping.\n",filename.c_str());
      exit(1);
   }
  
   table.close();
}

void TIncljets::ReadPDF(string filename){
   string path = "/h1/h1gen/lhapdf/LHAPDFv4/PDFsets/";
   path += filename;
   char* buffer = new char[path.length()+1];
   strcpy(buffer,path.c_str());
   lhapdf = new LHAPDFWrap(buffer, 0);
   lhapdf->initPDF(0);
}

void TIncljets::SetPDFSet(int set){
   if(lhapdf)
      lhapdf->initPDF(set);
}

void TIncljets::FillPDFCache(double muffactor){
   vector< double > xfx;
     for(int i=0;i<nQ2;i++){ // Q2
         for(int j=0;j<npt[i];j++){ // pt
            double xlim = xlimit[i][j];
            double hxlim = log10(xlim);
            double muf = mufval[i][j] * muffactor;
             for(int k=0;k<ntot;k++){ // x
                double hx = hxlim *(1. - (k-1.)/(ntot));
                double x = pow(10,hx);
                xfx = lhapdf->xfx(x,muf);
                // 0..5 = tbar, ..., ubar, dbar;
                // 6 = g;
                // 7..12 = d, u, ..., t
                pdf[i][j][k][0] = xfx[6]; //gluon

                pdf[i][j][k][1] = 0.;
                for(int l=0;l<6;l++){
                   pdf[i][j][k][1] += xfx[l] + xfx[l+7]; // sigma
                }
                pdf[i][j][k][2] = 0.;
                for(int l=0;l<13;l++){
                   double x = (l==6 ? 0.0 : xfx[l]);
                   if (!(l&1)) x *= 4.;
                   pdf[i][j][k][2] += x; // delta
               }
                pdf[i][j][k][2] /= 9.;
             }
         }
     }
}

void TIncljets::ResetXsection(){
   xsection.resize(nQ2);
   for(int i=0;i<nQ2;i++){ // Q2
      xsection[i].resize(npt[i]);
      for(int j=0;j<npt[i];j++){ // pt
         xsection[i][j].resize(3); // LO,NLO, Sum
         for(int m=0;m<3;m++) 
            xsection[i][j][m] = 0.;
      }
   }

   rebinned.resize(4);
   for(int i=0;i<4;i++){ // Q2
      rebinned[i].resize(4);
      for(int j=0;j<4;j++){ // pt
         rebinned[i][j].resize(3); // LO,NLO, Sum
         for(int m=0;m<3;m++) 
            rebinned[i][j][m] = 0.;
      }
   }


}

void TIncljets::CalcXsection(double asmz, int scale){
   for(int i=0;i<nQ2;i++){ // Q2
      for(int j=0;j<npt[i];j++){ // pt
         double mur = murval[i][j] * murscale[scale];
         double as = alphas(mur,asmz)/(2*PI); 
         for(int k=0;k<ntot;k++){ // x
            for(int l=0;l<nsubproc;l++){ // subprocesses
               xsection[i][j][0] += disweights[i][j][k][0][l] ; // LO
               //               xsection[i][j][0] += disweights[i][j][k][0][l] * pow(as,1) * pdf[i][j][k][l]; // LO
               xsection[i][j][1] += disweights[i][j][k][scale][l] * pow(as,2) * pdf[i][j][k][l]; // NLO
               xsection[i][j][2] += xsection[i][j][0] + xsection[i][j][1]; // total
            }
         }
      }
   }
}

double TIncljets::GetXsection(int Q2, int pt, int order){
   return xsection[Q2][pt][order];
}

void TIncljets::Rebin(){
   for(int i=0;i<4;i++){ // Q2
      for(int j=0;j<4;j++){ // pt
         for(int k=0;k<nQ2;k++){ // Q2
            double Q2center = (Q2high[k+1]+Q2high[k])/2;
            for(int l=0;l<npt[k];l++){ // pt
               double ptcenter = (pthigh[k][l+1]+pthigh[k][l])/2;
               if(Q2center>q2binning[i] && Q2center<q2binning[i+1]
                  && ptcenter>ptbinning[j] && ptcenter<ptbinning[j+1] ){
                  rebinned[i][j][0] += xsection[k][l][0];
                  rebinned[i][j][1] += xsection[k][l][1];
                  rebinned[i][j][2] += xsection[k][l][2];
               }
            }
         }
         double norm = 1000./(q2binning[i+1]-q2binning[i])/(ptbinning[j+1]-ptbinning[j]);
         rebinned[i][j][0] *= norm;
         rebinned[i][j][1] *= norm;
         rebinned[i][j][2] *= norm;
      }
   }
}

double TIncljets::GetRebinned(int Q2, int pt, int order){
   return rebinned[Q2][pt][order];
}

void TIncljets::ResetReference(){
   reference.resize(4);
   for(int i=0;i<4;i++){ // Q2
      reference[i].resize(4);
      for(int j=0;j<4;j++){ // pt
         reference[i][j].resize(3); // LO,NLO, Sum
         for(int m=0;m<3;m++) 
            reference[i][j][m] = 0.;
      }
   }
}


void TIncljets::ReadReference(string filename,int order){

   fstream table(filename.c_str(),ios::in); // open file
   double binl,binh,temp,cross,uncer;
   

   for (int i=0; i<4; i++) {
      table.ignore(256,'\n');
      for (int j=0; j<4; j++) {
         table >> binl >> temp >> binh >> temp >> temp >> temp
                >> temp >> temp >> temp >> cross >> uncer;
         reference[i][j][order] = cross/(q2binning[i+1]-q2binning[i]);
      }
      table.ignore(256,'\n');
      table.ignore(256,'\n');
   }
}

   
double TIncljets::GetReference(int Q2, int pt, int order){
   return reference[Q2][pt][order];
}

