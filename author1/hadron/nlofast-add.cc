// last modification 2005/12/07 - TK

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

const int kdis   = 1;
const int kpp    = 2;
const int kppbar = 3;

// return values of checkfile 
const int kcferror  = 0;
const int kcflo     = 1;
const int kcfnlo    = 2;
const int kcfnnlo   = 3;

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

   int nxtot;      // no of xbins 
   vector< vector<double> >xlimit; // array for lower x limits
   vector <vector< vector < vector < vector <weight_hhc> > > > > weights; // array for the weights MAIN ARRAY
   int nsubproc;   // no of subprocesses
   
   
   vector< double> nevents; // no of events calculated so far for each order
   int nwrite; // no of events after to write out the table

   int ireaction;
   double s;
   int iproc;
   int ialgo;
   double JetResol1;
   double JetResol2;
   int nord;
   int nordmax;
   vector< int> npow; // indicates power of alpha_s for each order 1..nord 
   vector< std::basic_string<char> > label; // comments for the orders: "LO","NLO" etc.
   int ixscheme;
   int ipdfwgt;


bool cmp(double x1, double x2);
int checkfile(string filename);
int resetbuffer(string filename);
int addfile(string filename);
int writetable(string filename);

int main(int argc,void** argv)
{
  // check parameters
  if(argc<3){
    printf("Usage: nlofast-add file1.raw [filex.raw]+ result.txt[.gz]\n");
    return(1);
  }
  // loop over arguments
  int nparts = argc-1;
  vector <string> lofiles;
  vector <string> nlofiles;
  vector <string> nnlofiles;
  vector <string> outfiles; 

  for(int i=0;i<nparts-1;i++){
     const char* path=(char *)argv[i+1];
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
         case kcfnnlo:
           nnlofiles.push_back(path);
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
  int nlosize =  nlofiles.size();
  int nnlosize = nnlofiles.size();
  int outsize = outfiles.size();

  printf("Found:\n %d LO file(s)\n %d NLO file(s)\n %d NNLO file(s).\n",losize, nlosize,nnlosize);

  if(losize>0 || nlosize>0 || nnlosize>0){
     if(losize==0){
        printf("Using only higher order tables is not supported. Use at least one LO file.\n");
        return 1;
     }
     if(nlosize==0){
        printf("Processing only LO data.\n");
     }
     if(nnlosize>0 && nlosize==0){
       printf("Using threshold corrections without NLO data is not supported. Use at least one NLO file.\n");
       return 1;
     }
     if(outsize>0) printf("File for writing results: %s.\n",outfiles[0].c_str());
     int locount  = 0;
     int nlocount = 0;
     int nnlocount = 0;

     if(nnlosize>0){
       resetbuffer(nnlofiles[nnlocount++]);  // Initialise using an NNLO file     
     }else{
       if(nlosize>0){
	 resetbuffer(nlofiles[nlocount++]); // Initialise using an NLO file
       }else{
	 resetbuffer(lofiles[locount++]); // Initialise using a LO file
       }    
     }
     for(;nnlocount<nnlosize;nnlocount++)
        addfile(nnlofiles[nnlocount]);
     for(;nlocount<nlosize;nlocount++)
        addfile(nlofiles[nlocount]);
     for(;locount<losize;locount++)
        addfile(lofiles[locount]);

      if(outsize>0){
        printf("Writing table %s in text format. Total number of events: %.2G(LO)",outfiles[0].c_str(),nevents[0]);
 	if(nlosize>0){
           printf(" %.2G(NLO)",nevents[1]);
        }
	if(nnlosize>0){
           printf(" %.2G(NNLO)",nevents[2]);
        }
        printf(".\n");
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
   char bufferchar;

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
   READ(nord);
   if(nord<1 || nord>3){
      printf("nord = %d not supported. Stopping.\n",nord);
      exit(1);
   }
   nordmax = nord; // The highest order is always read in first

   npow.resize(nordmax);
   for(int i=0;i<nordmax;i++) npow[i]=0;

   label.resize(nordmax);
   for(int i=0;i<nordmax;i++) label[i]="not known";

   READ(npow[nord-1]);

   label[nord-1]="";
   READ(bufferchar);
   while(bufferchar!='\n' && (!table.eof()) ){
      label[nord-1] += bufferchar;
      READ(bufferchar);
   }

   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after label in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   nevents.resize(nordmax);

   READ(nevents[nord-1]);
   switch(nord){
   case 1:
      printf("Reading LO file %s with %.2G events.\n",filename.c_str(),nevents[0]);
      break;
   case 2:
      printf("Reading NLO file %s with %.2G events.\n",filename.c_str(),nevents[1]);
      break;
   case 3:
      printf("Reading NNLO file %s with %.2G events.\n",filename.c_str(),nevents[2]);
      break;
   default:
      printf("nord = %d not supported. Stopping.\n",nord);
      exit(1);
   }

   READ(nxtot);
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

   int ncontrib = 0; // Number of orders including scale variations to fill in the table 
   //scale variations for (N)NLO
   if(nord>1){
      READ(nscalevar); // there may be several scale variations for (N)NLO
      ncontrib = nscalevar; // and these are the contributions to fill in the table 
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
      // When we're here this means we have a pure LO file to write
      ncontrib = 1; // one contribution for LO
      nscalevar = 0; // no scale variations
   }
   
   int offset = 0;
   switch(nord){
   case 1: offset =0; break;
   case 2: offset =1; break;
   case 3: offset =1 + nscalevar; break;
   default: printf("nord= %d not supported. Exit.\n",nord); exit(1);
   }
   if(ireaction == kpp || ireaction == kppbar){
      weights.resize (nrap);
      for(int i=0;i<nrap;i++){ // rapidity
         weights[i].resize(npt[i]);
         for(int j=0;j<npt[i];j++){ // pt
            weights[i][j].resize(nxtot);
            for(int k=0;k<nxtot;k++){ // xmin
               weights[i][j][k].resize(k+1); // half matrix xmin,xmax: (n^2+n)/2
               for(int l=0;l<min(k+1,nxtot);l++){ // xmax
		 weights[i][j][k][l].resize((nord-1)*nscalevar+1); // 1 for LO and nscalevar per higher order
                  for(int scalevar=0; scalevar<(nord-1)*nscalevar+1;scalevar++){ 
                     weights[i][j][k][l][scalevar].resize(nsubproc); // seven entries in weight
                  }
                  for(int scalevar=0; scalevar<ncontrib;scalevar++){ //  1 for LO and nscalevar per higher order
                     for(int m=0;m<nsubproc;m++){     //subprocesses
                        READ(weights[i][j][k][l][offset+scalevar][m]);
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
   char bufferchar;
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

   READ(nord);
   if(nord<1 || nord>2){
      printf("nord = %d not supported. Stopping.\n",nord);
      exit(1);
   }

   READ(bufferint);
   if(npow[nord-1]==0){ // not yet filled?
      npow[nord-1] =  bufferint;
   }else if( bufferint !=  npow[nord-1]){
      printf("Differing npow[%d]=%d in %s. Stopping.\n",nord-1,bufferint,filename.c_str());
      exit(1);
   }
   
   label[nord-1]="";
   READ(bufferchar);
   while(bufferchar!='\n' && (!table.eof()) ){
      label[nord-1] += bufferchar;
      READ(bufferchar);
   }

   READ(marker);// ------------------END of block
   if(marker!=cmarker){
      printf("Misaligned table after label in %s. Stopping.\n",filename.c_str());
      exit(1);
   }

   READ(bufferdbl);
   nevents[nord-1] += bufferdbl;

   switch(nord){
   case 1:
      printf("Reading LO file %s with %.2G events.\n",filename.c_str(),bufferdbl);
      break;
   case 2:
      printf("Reading NLO file %s with %.2G events.\n",filename.c_str(),bufferdbl);
      break;
   case 3:
      printf("Reading NNLO file %s with %.2G events.\n",filename.c_str(),bufferdbl);
      break;
   default:
      printf("nord = %d not supported. Stopping.\n",nord);
      exit(1);
   }

   READ(bufferint);
   if(bufferint != nxtot){
      printf("Differing nxtot in %s. Stopping.\n",filename.c_str());
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

   if(nord==1){
      scalevarmax = 1;
   }else if(nord==2){   //scale variations for NLO
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

   int offset = 0;
   switch(nord){
   case 1: offset =0; break;
   case 2: offset =1; break;
   case 3: offset =1 + nscalevar; break;
   default: printf("nord= %d not supported. Exit.\n",nord); exit(1);
   }
   if(ireaction == kpp || ireaction == kppbar){
      for(int i=0;i<nrap;i++){ // rapidity
         for(int j=0;j<npt[i];j++){ // pt
            for(int k=0;k<nxtot;k++){ // xmin
               for(int l=0;l<min(k+1,nxtot);l++){ // xmax
		 for(int scalevar=0; scalevar<scalevarmax;scalevar++){ //1 for LO and nscalevar per higher order
                     for(int m=0;m<nsubproc;m++){     //subprocesses
                        READ(bufferdbl);
                        weights[i][j][k][l][offset+scalevar][m] += bufferdbl;
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

   //nord
   WRITE(nordmax);

   //npow
   for(int i=0;i<nordmax;i++){
      WRITE(npow[i]);
   }

   //label
   for(int i=0;i<nordmax;i++){
      *table << label[i] << endl;
   }

   WRITE(marker);// ------------------END of block

   //nevt
   for(int i=0;i<nordmax;i++){
      WRITE(nevents[i]);
   }

   //nxtot
   WRITE(nxtot);

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
      WRITE(murscale[i]);
   }
   for(int i=0;i<nscalevar;i++){
      WRITE(mufscale[i]);
   }
   WRITE(marker);// ------------------END of block

   for(int i=0;i<nordmax;i++){
      if(nevents[i]==0) nevents[i]=1;
   }
   
   if(ireaction == kpp || ireaction == kppbar){
      for(int i=0;i<nrap;i++){ // rapidity
         for(int j=0;j<npt[i];j++){ // pt
            for(int k=0;k<nxtot;k++){ // xmin
               for(int l=0;l<min(k+1,nxtot);l++){ // xmax
                  for(int m=0;m<nsubproc;m++){     //subprocesses
                     for(int scalevar=0; scalevar<1+(nordmax-1)*nscalevar;scalevar++){  // 1 for LO and nscalevar per higher order
                        //                        WRITE(weights[i][j][k][l][scalevar][m]/nevents[scalevar==0? 0 : 1]);
                        WRITE(weights[i][j][k][l][scalevar][m]/nevents[(int)ceil(scalevar/(float)nscalevar)]);
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
  
   int bufferint,nord,marker;
   double bufferdbl;
   std::basic_string<char> bufferstring;
   char bufferchar;

   READ(bufferint);
   READ(bufferdbl);
   READ(bufferint);
   READ(bufferint);
   READ(bufferdbl);
   READ(bufferdbl);
   READ(nord);
   READ(bufferint);
   
   READ(bufferchar);
   while(bufferchar!='\n'){
      bufferstring += bufferchar;
      READ(bufferchar);
   }
   READ(marker);
   if(marker!=cmarker){
      printf("Misaligned table after npow.\n");
      return kcferror;
   }
   if(nord==1){
      return kcflo;
   }
   if(nord==2){
      return kcfnlo;
   }
   if(nord==3){
      return kcfnnlo;
   }
   printf("Invalid nord=%d.\n",nord);
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
