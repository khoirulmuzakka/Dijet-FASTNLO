// last modification 2005/09/08 - TK

#include <stdio.h>
#include <unistd.h>             
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "gzstream.h"


using namespace std;

typedef vector<double> weight_hhc;
typedef vector<double> weight_dis;

// return values of checkfile 
const int kdis   = 1;
const int kpp    = 2;
const int kppbar = 3;

// return values of checkfile 
const int kcferror  = 0;
const int kcflo     = 1;
const int kcfnlo    = 2;

   const int cmarker = 1234567890; //used to separate section in the table

// binning
   
   int nrap;       // No of rapidity bins 
   double *raphigh;  // array for rapidity boundaries
   int *npt;       // No of pT bins in each y range
   int nptmax;     // maximum number of pt bins
   vector< vector<double> >pthigh;   // array for pT boundaries
   double ptlow;    // lowest pt considered 
   
   int nscalevar;                 // number of scale variations (mu_r, mu_f) for NLO
   vector <double> murscale;        // overall scale factor for renormalization scale
   vector< vector<double> >murval;   // array for renormalization scale values
   vector <double> mufscale;        // overall scale factor for factorization scale
   vector< vector<double> >mufval;    // array for factorization scale values

   int ntot;      // no of xbins 
   vector< vector<double> >xlimit; // array for lower x limits
   vector <vector< vector < vector < vector <weight_hhc> > > > > weights; // array for the weights MAIN ARRAY
   int nsubproc;   // no of subprocesses
   
   
   vector< double> nevents; // no of events calculated so far for each scale variation
   int nwrite; // no of events after to write out the table

   int ireaction;
   double      s;
   int iproc;
   int ialgo;
   double JetResol1;
   double JetResol2;
   int npow;
   int npowmax = 0;
   int Oalphas;
   int ixscheme;
   int ipdfwgt;

bool patchscales = false;


bool cmp(double x1, double x2);
int checkfile(string filename);
int resetbuffer(string filename);
int addfile(string filename);
int writetable(string filename);

int main(int argc,void** argv)
{
  // check parameters
  if(argc<3){
    printf("Usage: nlofast-add [--patch-scales] file1.raw [filex.raw]+ result.txt[.gz]\n");
    return(1);
  }
  // loop over arguments
  int nparts = argc-1;
  vector <string> lofiles;
  vector <string> nlofiles;
  vector <string> outfiles; 

  for(int i=0;i<nparts-1;i++){
     const char* path=(char *)argv[i+1];
     if(strcmp(path,"--patch-scales")==0){
        printf("ATTENTION: murscale[] and mufscale[] will be modified by taking square root.\n");
        patchscales=true;
        continue;
     }
     // File there?
     if (access(path, R_OK) == 0){
        int result = checkfile(path);
        switch(result){
        case kcflo:
           lofiles.push_back(path);
           break;
        case kcfnlo:
           nlofiles.push_back(path);
           break;
        default:
           printf("Invalid format in %s, skipping.\n",path);      
           break;
        }
     }else{
        printf("Cannot access %s, skipping.\n",path);      
     }
  }
  const char* path=(char *)argv[argc-1];
  // check if result file is writable
  if (access(path, F_OK) == 0){
     printf("File for writing the table exists: %s.\n Please remove it first.\n",path);      
     exit(2);
  }else{
     outfiles.push_back(path);
  }

  int losize  =  lofiles.size();
  int nlosize = nlofiles.size();
  int outsize = outfiles.size();

  printf("Found %d LO file(s) and %d NLO file(s).\n",losize, nlosize);
  if(losize>0 || nlosize>0){
     if(losize==0) printf("Using dummy LO data.\n");
     if(nlosize==0) printf("Using dummy NLO data.\n");
     if(outsize>0) printf("File for writing results: %s.\n",outfiles[0].c_str());
     int locount  = 0;
     int nlocount = 0;
     if(nlosize>0){
        resetbuffer(nlofiles[nlocount++]); // Initialise using an NLO file
     }else{
        resetbuffer(lofiles[locount++]); // Initialise using an LO file
     }    
     for(;nlocount<nlosize;nlocount++)
        addfile(nlofiles[nlocount]);
     for(;locount<losize;locount++)
        addfile(lofiles[locount]);

     if(outsize>0){
        printf("Writing table %s in text format. Total number of events: %.2G (LO) %.2G (NLO).\n",outfiles[0].c_str(),nevents[0],nevents[1]);
        writetable(outfiles[0]);
     }
  }else{
      printf("Nothing to add.\n");
  }
  return 0;
}

int resetbuffer(string filename){
#define READ(n) table.read(reinterpret_cast<char *>(&n),sizeof(n))
   fstream table(filename.c_str(),ios::in|ios::binary); // open file
   int marker;

   READ(ireaction);
   switch(ireaction){
   case kpp:
   case kppbar:
      nsubproc = 7;
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
   npowmax = npow;
   if(npow==1){
      READ(nevents[0]);
      printf("Reading LO file %s with %.2G events.\n",filename.c_str(),nevents[0]);
   }else if(npow==2){
      READ(nevents[1]);
      printf("Reading NLO file %s with %.2G events.\n",filename.c_str(),nevents[1]);
   }
   else{
      printf("Npow = %d not supported. Stopping.\n",npow);
      exit(1);
   }

   READ(ntot);
   READ(ixscheme);
   READ(ipdfwgt);
   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after ipdfwght in %s. Stopping.\n",filename.c_str());
      exit(1);
   }
   READ(nrap);
   raphigh = new double[nrap+1];      //----- array for rapidity boundaries
   for(int i=0;i<nrap+1;i++){
      READ(raphigh[i]); //Rap[0] ... Rap[Nrapidity]
   }

   npt     = new  int[nrap];         // nrap bins in rapidity with each npt[irap] bins in pt
   for(int i=0;i<nrap;i++){
      READ(npt[i]); //Npt[0] ... Npt[Nrapidity-1] 
   }

   pthigh.resize(nrap);
   for(int i=0;i<nrap;i++){
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
   
   xlimit.resize (nrap);
   for(int i=0;i<nrap;i++){
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

   murval.resize (nrap);
   for(int i=0;i<nrap;i++){
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


   mufval.resize (nrap);
   for(int i=0;i<nrap;i++){
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
   if(npow==2){
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
   }else{
      // When we're here this means we should use dummy NLO data
      nscalevar=1;
      murscale.resize(nscalevar);
      mufscale.resize(nscalevar);
      murscale[0] = 1.0;
      mufscale[0] = 1.0;
   }

   nevents.resize(1+nscalevar); // Room for LO and  NLO with scale variations
   for(int scalevar=1; scalevar<nscalevar;scalevar++){
      nevents[scalevar+1] = nevents[scalevar]; // Up to now all scale varaitions have same number of events
   }
   
   if(ireaction == kpp || ireaction == kppbar){
      weights.resize (nrap);
      for(int i=0;i<nrap;i++){ // rapidity
         weights[i].resize(npt[i]);
         for(int j=0;j<npt[i];j++){ // pt
            weights[i][j].resize(ntot);
            for(int k=0;k<ntot;k++){ // xmin
               weights[i][j][k].resize(k+1); // half matrix xmin,xmax: (n^2+n)/2
               for(int l=0;l<min(k+1,ntot);l++){ // xmax
                  weights[i][j][k][l].resize(1+nscalevar);
                  for(int scalevar=0; scalevar<1+nscalevar;scalevar++){ //LO & NLO with scale variations
                     weights[i][j][k][l][scalevar].resize(nsubproc); // seven entries in weight
                  }
                  for(int scalevar=0; scalevar<nscalevar;scalevar++){ //either (LO) or (NLO with scale variations)
                     for(int m=0;m<nsubproc;m++){     //subprocesses
                        READ(weights[i][j][k][l][npow-1+scalevar][m]);
                     }
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
   return 0;
}

int addfile(string filename){
    fstream table(filename.c_str(),ios::in|ios::binary); // open file
   int marker;
   int scalevarmax;
  
   int bufferint;
   double bufferdbl;

   READ(bufferint);
   if(bufferint != ireaction){
      printf("Differing ireaction in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   READ(bufferdbl);
   if(!cmp(bufferdbl,s)){
      printf("Differing s in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   READ(bufferint);
   iproc = bufferint;
   if(bufferint != iproc){
      printf("Differing iproc in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   READ(bufferint);
   if(bufferint != ialgo){
      printf("Differing ialgo in %s. Stopping.\n",filename.c_str());
      exit(1);
   }


   READ(bufferdbl);
   if(!cmp(bufferdbl,JetResol1)){
      printf("Differing JetResol1 in %s. Stopping.\n",filename.c_str());
      exit(1);
   }
   READ(bufferdbl);
   if(!cmp(bufferdbl,JetResol2)){
      printf("Differing JetResol2 in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   READ(npow);

   READ(bufferint);
   if(bufferint > Oalphas){
      Oalphas = bufferint; // Write highest found power of a_s to result table
   }

   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after Oalphas in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   READ(bufferdbl);
   if(npow==1){
     nevents[0]+=bufferdbl;
     printf("Reading LO file %s with %.2G events.\n",filename.c_str(),bufferdbl);
   }else{
     if(npow==2){
       for(int scalevar=1; scalevar<nscalevar+1;scalevar++){
	 nevents[scalevar] += bufferdbl;
       }
       printf("Reading NLO file %s with %.2G events.\n",filename.c_str(),bufferdbl);
     }else{
       printf("Npow = %d not supported. Stopping.\n",npow);
       exit(1);
     }
   }
   if (npow>npowmax) npowmax=npow;

   READ(bufferint);
   if(bufferint != ntot){
      printf("Differing ntot in %s. Stopping.\n",filename.c_str());
      exit(1);
   }
   READ(bufferint);
   if(bufferint != ixscheme){
      printf("Differing ixscheme in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   READ(bufferint);
   if(bufferint != ipdfwgt){
      printf("Differing ipdfwgt in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after ipdfwght in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   READ(bufferint);
   if(bufferint != nrap){
      printf("Differing nrap in %s. Stopping.\n",filename.c_str());
      exit(1);
   }


   for(int i=0;i<nrap+1;i++){
      READ(bufferdbl); //Rap[0] ... Rap[Nrapidity]
      if(!cmp(bufferdbl,raphigh[i])){
         printf("Differing raphigh in %s. Stopping.\n",filename.c_str());
         printf("%f %f\n",bufferdbl,raphigh[i]);
         exit(1);
      }
   }

   for(int i=0;i<nrap;i++){
      READ(bufferint); //Npt[0] ... Npt[Nrapidity-1] 
      if(bufferint != npt[i]){
         printf("Differing npt in %s. Stopping.\n",filename.c_str());
         exit(1);
      }
   }

   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i]+1;j++){
         READ(bufferdbl); // all pt bins 
         if(!cmp(bufferdbl,pthigh[i][j])){
            printf("Differing pthigh in %s. Stopping.\n",filename.c_str());
            exit(1);
         }
      }
   }
    
   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after pthigh in %s. Stopping.\n",filename.c_str());
      exit(1);
   }
   
   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         READ(bufferdbl); // all x limits 
         if(!cmp(bufferdbl,xlimit[i][j])){
            printf("Differing xlimit in %s. Stopping.\n",filename.c_str());
            exit(1);
         }
      }
   }

   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after xlimit in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         READ(bufferdbl); // all murval 
         if(!cmp(bufferdbl,murval[i][j])){
            printf("Differing murval in %s. Stopping.\n",filename.c_str());
            exit(1);
         }
      }
   }
   
   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after murval in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         READ(bufferdbl); // all mufval 
         if(!cmp(bufferdbl,mufval[i][j])){
            printf("Differing mufval in %s. Stopping.\n",filename.c_str());
            exit(1);
         }
      }
   }

   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after mufval in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   //scale variations for NLO
   if(npow==2){
      READ(bufferint);
      if(bufferint!=nscalevar){
         printf("Differing nscalevar in %s. Stopping.\n",filename.c_str());
         exit(1);
      }

      for(int i=0;i<nscalevar;i++){
         READ(bufferdbl);
         if(!cmp(bufferdbl,murscale[i])){
            printf("Differing murscale[%d] in %s. Stopping.\n",i,filename.c_str());
            exit(1);
         }
      }
      for(int i=0;i<nscalevar;i++){
         READ(bufferdbl);
         if(!cmp(bufferdbl,mufscale[i])){
            printf("Differing mufscale[%d] in %s. Stopping.\n",i,filename.c_str());
            exit(1);
         }
      }

      READ(marker);// ------------------END of block
      if(marker!=cmarker){
         printf("Misaligned table after mufscale in %s. Stopping.\n",filename.c_str());
         exit(1);
      }
      scalevarmax = nscalevar;
   }
   else{ if(npow==1){
         //LO file
         scalevarmax = 1;
      }else{
         printf("Npow = %d not supported. Stopping.\n",npow);
      }
   }

   if(ireaction == kpp || ireaction == kppbar){
      for(int i=0;i<nrap;i++){ // rapidity
         for(int j=0;j<npt[i];j++){ // pt
            for(int k=0;k<ntot;k++){ // xmin
               for(int l=0;l<min(k+1,ntot);l++){ // xmax
                  for(int scalevar=0; scalevar<scalevarmax;scalevar++){ //LO & NLO scale variations
                     for(int m=0;m<nsubproc;m++){     //subprocesses
                        READ(bufferdbl);
                        weights[i][j][k][l][npow-1+scalevar][m] += bufferdbl;
                     }
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
   return 0;
}

int writetable(string filename){
   //    #define WRITE(n) table.write(reinterpret_cast<char *>(&n),sizeof(n))
    #define WRITE(n) *table << setprecision(40) << n << "\n"
   int marker = 1234567890; //used to separate section in the table

   ostream *table;

   int len = strlen(filename.c_str());
   if(len<3){
      printf ("writetable: strange file name: %s",filename.c_str());
      return 2;
   }
   char extension[4];
   strncpy(extension,(char*)(filename.c_str()+len-3),4);
   if(strcmp(extension,".gz")==0){
      printf("Writing compressed format.\n");
      table = new ogzstream(filename.c_str(),ios::out); // open .gz file
   }else{
      printf("Writing uncompressed format.\n");
      table = new fstream(filename.c_str(),ios::out); // open file
   }

   //really open?
   if(!table->good()){
       printf("Cannot open %s for writing. Stopping.\n",filename.c_str());
       exit(2);
    }

   // ireaction
   WRITE(ireaction);

   // Ecms
   WRITE(s);

   //iproc
   WRITE(iproc);

   //iproc
   WRITE(ialgo);

   //JetResol1
   WRITE(JetResol1);

   //JetResol2
   WRITE(JetResol2);

   //npow
   WRITE(npowmax);

   //Oalphas
   WRITE(Oalphas);
   
   WRITE(marker);// ------------------END of block

   //nevt LO
   WRITE(nevents[0]);

   //nevt NLO
   WRITE(nevents[1]);

   //ntot
   WRITE(ntot);

   //ixscheme
   WRITE(ixscheme);

   //ipdfwgt
   WRITE(ipdfwgt);

   WRITE(marker);// ------------------END of block

   //Nrapidity
   WRITE(nrap);
   for(int i=0;i<nrap+1;i++){
      WRITE(raphigh[i]); //Rap[0] ... Rap[Nrapidity]
   }

   for(int i=0;i<nrap;i++){
      WRITE(npt[i]); //Npt[0] ... Npt[Nrapidity-1] 
   }

   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i]+1;j++){
         WRITE(pthigh[i][j]); // all pt bins 
      }
   }
    
   WRITE(marker);// ------------------END of block
   
    for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         WRITE(xlimit[i][j]); // all x limits 
      }
   }

   WRITE(marker);// ------------------END of block

   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         WRITE(murval[i][j]); // all murval 
      }
   }
   
   WRITE(marker);// ------------------END of block

   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         WRITE(mufval[i][j]); // all mufval 
      }
   }

   WRITE(marker);// ------------------END of block
   WRITE(nscalevar);
   for(int i=0;i<nscalevar;i++){
      if(patchscales){
         murscale[i] = sqrt(murscale[i]);
      }
      WRITE(murscale[i]);
   }
   for(int i=0;i<nscalevar;i++){
      if(patchscales){
         mufscale[i] = sqrt(mufscale[i]);
      }
      WRITE(mufscale[i]);
   }
   WRITE(marker);// ------------------END of block

   for(int scalevar=0;scalevar<1+nscalevar;scalevar++){
      if(nevents[scalevar]==0) nevents[scalevar]=1;
   }
   
   if(ireaction == kpp || ireaction == kppbar){
      for(int i=0;i<nrap;i++){ // rapidity
         for(int j=0;j<npt[i];j++){ // pt
            for(int k=0;k<ntot;k++){ // xmin
               for(int l=0;l<min(k+1,ntot);l++){ // xmax
                  for(int m=0;m<nsubproc;m++){     //subprocesses
                     for(int scalevar=0; scalevar<1+nscalevar;scalevar++){ //LO & NLO with scale variations
                        WRITE(weights[i][j][k][l][scalevar][m]/nevents[scalevar]);
                     }
                  }
               }
            }
         }
      }
   }

   WRITE(marker);// ------------------END of table

   table->flush();
   delete table;
   return 0;
}

int checkfile(string filename){
   fstream table(filename.c_str(),ios::in|ios::binary); // open file
  
   int bufferint,npow,marker;
   double bufferdbl;

   READ(bufferint);
   READ(bufferdbl);
   READ(bufferint);
   READ(bufferint);
   READ(bufferdbl);
   READ(bufferdbl);
   READ(npow);
   READ(bufferint);
   READ(marker);
   if(marker!=cmarker){
      return kcferror;
   }
   if(npow==1){
      return kcflo;
   }
   if(npow==2){
      return kcfnlo;
   }
   return kcferror;
}

bool cmp(double x1, double x2){
   double norm;
   if (x1>0.){
      norm = x1;
   }else{
      norm = 1.; // If x1 is 0, do not try to calculate relative deviation, use absolute
   }

   return(((x1-x2)/norm)<1e-9);
}
