// last modification 2005/08/15 - TK

#include <stdio.h>
#include <unistd.h>             
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "gzstream.h"


using namespace std;

typedef vector<double> weight_hhc;
 
   int marker = 1234567890; //used to separate section in the table

  // binning
   
   int nrap;       // No of rapidity bins 
   double *raphigh;  // array for rapidity boundaries
   int *npt;       // No of pT bins in each y range
   int nptmax;     // maximum number of pt bins
   vector< vector<double> >pthigh;   // array for pT boundaries
   double ptlow;    // lowest pt considered 
   
   double murscale;        // overall scale factor for renormalization scale
   vector< vector<double> >murval; // array for renormalization scale values
   double mufscale;        // overall scale factor for factorization scale
   vector< vector<double> >mufval; // array for factorization scale values

   int ntot;      // no of xbins 
   vector< vector<double> >xlimit; // array for lower x limits
   vector <vector< vector < vector < vector <weight_hhc> > > > > weights; // array for the weights MAIN ARRAY
   
   
   double nevents[6]; // no of events calculated so far for each scale variation
   int nwrite; // no of events after to write out the table
   double xsectsum; // total cross section - internal counter for test output

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



int resetbuffer(const char* filename);
int addfile(const char* filename);
int writetable(const char* filename);

int main(int argc,void** argv)
{
  // check parameters
  if(argc<3){
    printf("Usage: nlofast-add file1.raw [filex.raw]+ result.txt\n");
    return(1);
  }
  // loop over arguments
  int nparts = argc-1;
  bool firstpart = true;
  for(int i=0;i<nparts-1;i++){
     const char* path=(char *)argv[i+1];
     // File there?
    if (access(path, R_OK) == 0){
      if(firstpart){
	resetbuffer(path);
        firstpart = false;
      }else{
	addfile(path);
      }
    }else{
      printf("Cannot access %s, skipping.\n",path);      
    }
  }
  // anything there to write out?
  if(!firstpart){
     const char* path=(char *)argv[argc-1];
     // check if result file is writable
     if (access(path, F_OK) == 0){
        printf("File for writing the table exists: %s.\n Please remove it first.\n",path);      
        return 2;
     }
     printf("Writing table %s in text format. Total number of events: %.0f (LO) %.0f (NLO).\n",path,nevents[0],nevents[1]);      
     writetable(path);
  }
  return 0;
}

int resetbuffer(const char* filename){
#define READ(n) table.read(reinterpret_cast<char *>(&n),sizeof(n))
   const int cmarker = 1234567890; //used to separate section in the table
   fstream table(filename,ios::in|ios::binary); // open file
   int marker;

   READ(ireaction);
   READ(s);
   READ(iproc);
   READ(ialgo);
   READ(JetResol1);
   READ(JetResol2);
   READ(npow);
   READ(Oalphas);
   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after Oalphas in %s. Stopping.\n",filename);
      exit(1);
   }

   npowmax = npow;
   int scalevarmax;
   if(npow==1){
      READ(nevents[0]);
      printf("Reading LO file %s with %.0f events.\n",filename,nevents[0]);
      scalevarmax=1; //LO has no scale variations
   }else if(npow==2){
      READ(nevents[1]);
      for(int scalevar=1; scalevar<5;scalevar++){
         nevents[scalevar+1] = nevents[scalevar];
      }
      printf("Reading NLO file %s with %.0f events.\n",filename,nevents[1]);
      scalevarmax=5; // NLO: central plus 4 variations
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
      printf("Misaligned table after ipdfwght in %s. Stopping.\n",filename);
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
      printf("Misaligned table after pthigh in %s. Stopping.\n",filename);
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
      printf("Misaligned table after xlimit in %s. Stopping.\n",filename);
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
      printf("Misaligned table after murval in %s. Stopping.\n",filename);
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
      printf("Misaligned table after mufval in %s. Stopping.\n",filename);
      exit(1);
   }

  weights.resize (nrap);
   for(int i=0;i<nrap;i++){ // rapidity
      weights[i].resize(npt[i]);
      for(int j=0;j<npt[i];j++){ // pt
         weights[i][j].resize(ntot);
         for(int k=0;k<ntot;k++){ // xmin
            weights[i][j][k].resize(k+1); // half matrix xmin,xmax: (n^2+n)/2
            for(int l=0;l<min(k+1,ntot);l++){ // xmax
               weights[i][j][k][l].resize(6);
               for(int scalevar=0; scalevar<6;scalevar++){ //LO & 5*NLO scale variations
                  weights[i][j][k][l][scalevar].resize(7); // seven entries in weight
               }
               for(int scalevar=0; scalevar<scalevarmax;scalevar++){ //LO & 5*NLO scale variations
                  for(int m=0;m<7;m++){     //subprocesses
                     READ(weights[i][j][k][l][npow-1+scalevar][m]);
                  }
               }
            }
         }
      }
   }

   READ(marker);// ------------------END of table
   if(marker!=cmarker){
      printf("Misaligned table after weights in %s. Stopping.\n",filename);
      exit(1);
   }
  
   table.close();
   return 0;
}

int addfile(const char* filename){
   const int cmarker = 1234567890; //used to separate section in the table
   fstream table(filename,ios::in|ios::binary); // open file
   int marker;
  
   int bufferint;
   double bufferdbl;

   READ(bufferint);
   if(bufferint != ireaction){
      printf("Differing ireaction in %s. Stopping.\n",filename);
      exit(1);
   }

   READ(bufferdbl);
   if(bufferdbl != s){
      printf("Differing s in %s. Stopping.\n",filename);
      exit(1);
   }

   READ(bufferint);
   iproc = bufferint;
   if(bufferint != iproc){
      printf("Differing iproc in %s. Stopping.\n",filename);
      exit(1);
   }

   READ(bufferint);
   if(bufferint != ialgo){
      printf("Differing ialgo in %s. Stopping.\n",filename);
      exit(1);
   }


   READ(bufferdbl);
   if(bufferdbl != JetResol1){
      printf("Differing JetResol1 in %s. Stopping.\n",filename);
      exit(1);
   }
   READ(bufferdbl);
   if(bufferdbl != JetResol2){
      printf("Differing JetResol2 in %s. Stopping.\n",filename);
      exit(1);
   }

   READ(npow);

   READ(bufferint);
   if(bufferint > Oalphas){
      Oalphas = bufferint; // Write highest found power of a_s to result table
   }

   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after Oalphas in %s. Stopping.\n",filename);
      exit(1);
   }
 

   READ(bufferdbl);
   int scalevarmax;
   if(npow==1){
      nevents[0]+=bufferdbl;
      printf("Reading LO file %s with %.0f events.\n",filename,bufferdbl);
      scalevarmax=1; //LO has no scale variations
   }else if(npow==2){
      nevents[1] += bufferdbl;
      for(int scalevar=1; scalevar<5;scalevar++){
         nevents[scalevar+1] += bufferdbl;
      }
      printf("Reading NLO file %s with %.0f events.\n",filename,bufferdbl);
      scalevarmax=5; // NLO: central plus 4 variations
   }
   else{
      printf("Npow = %d not supported. Stopping.\n",npow);
      exit(1);
   }
   if (npow>npowmax) npowmax=npow;

   READ(bufferint);
   if(bufferint != ntot){
      printf("Differing ntot in %s. Stopping.\n",filename);
      exit(1);
   }
   READ(bufferint);
   if(bufferint != ixscheme){
      printf("Differing ixscheme in %s. Stopping.\n",filename);
      exit(1);
   }

   READ(bufferint);
   if(bufferint != ipdfwgt){
      printf("Differing ipdfwgt in %s. Stopping.\n",filename);
      exit(1);
   }

   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after ipdfwght in %s. Stopping.\n",filename);
      exit(1);
   }

   READ(bufferint);
   if(bufferint != nrap){
      printf("Differing nrap in %s. Stopping.\n",filename);
      exit(1);
   }


   for(int i=0;i<nrap+1;i++){
      READ(bufferdbl); //Rap[0] ... Rap[Nrapidity]
      if(bufferdbl != raphigh[i]){
         printf("Differing raphigh in %s. Stopping.\n",filename);
         exit(1);
      }
   }

   for(int i=0;i<nrap;i++){
      READ(bufferint); //Npt[0] ... Npt[Nrapidity-1] 
      if(bufferint != npt[i]){
         printf("Differing npt in %s. Stopping.\n",filename);
         exit(1);
      }
   }

   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i]+1;j++){
         READ(bufferdbl); // all pt bins 
         if(bufferdbl != pthigh[i][j]){
            printf("Differing pthigh in %s. Stopping.\n",filename);
            exit(1);
         }
      }
   }
    
   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after pthigh in %s. Stopping.\n",filename);
      exit(1);
   }
   
   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         READ(bufferdbl); // all x limits 
         if(bufferdbl != xlimit[i][j]){
            printf("Differing xlimit in %s. Stopping.\n",filename);
            exit(1);
         }
      }
   }

   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after xlimit in %s. Stopping.\n",filename);
      exit(1);
   }

   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         READ(bufferdbl); // all murval 
         if(bufferdbl != murval[i][j]){
            printf("Differing murval in %s. Stopping.\n",filename);
            exit(1);
         }
      }
   }
   
   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after murval in %s. Stopping.\n",filename);
      exit(1);
   }

   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         READ(bufferdbl); // all mufval 
         if(bufferdbl != mufval[i][j]){
            printf("Differing mufval in %s. Stopping.\n",filename);
            exit(1);
         }
      }
   }

   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after mufval in %s. Stopping.\n",filename);
      exit(1);
   }

   for(int i=0;i<nrap;i++){ // rapidity
      for(int j=0;j<npt[i];j++){ // pt
         for(int k=0;k<ntot;k++){ // xmin
            for(int l=0;l<min(k+1,ntot);l++){ // xmax
               for(int scalevar=0; scalevar<scalevarmax;scalevar++){ //LO & 5*NLO scale variations
                  for(int m=0;m<7;m++){     //subprocesses
                     READ(bufferdbl);
                     weights[i][j][k][l][npow-1+scalevar][m] += bufferdbl;
                  }
               }
            }
         }
      }
   }

   READ(marker);// ------------------END of table
   if(marker!=cmarker){
      printf("Misaligned table after weights in %s. Stopping.\n",filename);
      exit(1);
   }
 
   table.close();
   return 0;
}

int writetable(const char* filename){
   //    #define WRITE(n) table.write(reinterpret_cast<char *>(&n),sizeof(n))
    #define WRITE(n) *table << setprecision(40) << n << "\n"
   int marker = 1234567890; //used to separate section in the table

   ostream *table;

   int len = strlen(filename);
   if(len<3){
      printf ("writetable: strange file name: %s",filename);
      return 2;
   }
   char extension[4];
   strncpy(extension,(char*)(filename+len-3),4);
   if(strcmp(extension,".gz")==0){
      printf("Writing compressed format.\n");
      table = new ogzstream(filename,ios::out); // open .gz file
   }else{
      printf("Writing uncompressed format.\n");
      table = new fstream(filename,ios::out); // open file
   }

   //really open?
   if(!table->good()){
       printf("Cannot open %s for writing. Stopping.\n",filename);
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

   for(int i=0;i<nrap;i++){ // rapidity
      for(int j=0;j<npt[i];j++){ // pt
         for(int k=0;k<ntot;k++){ // xmin
            for(int l=0;l<min(k+1,ntot);l++){ // xmax
               for(int m=0;m<7;m++){     //subprocesses
                  for(int scalevar=0; scalevar<6;scalevar++){ //LO & 5*NLO scale variations
                     WRITE(weights[i][j][k][l][scalevar][m]/nevents[scalevar]);
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

