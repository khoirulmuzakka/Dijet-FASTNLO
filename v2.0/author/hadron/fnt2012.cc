//
//   fastNLO v2 author code 
//   scenario: fnt2012
//   Run IIa D0 M-3jet measurement
//   for Run II midpoint cone algorithm 
//
// 
// ============== fastNLO user: ===================================
// To create your own scenario, it is recommended to take 
// this code, make a copy and edit the relevant changes.
// Important:
// Edit only those lines which are labeled as "fastNLO user"
// and refer to the documentation ("fastNLO authorcode in 
// NLOJET++") for a detailed explanation of the parameters 
// and variables.
// If a code fragment is not explicitely labeled as "fastNLO user",
// it is likely that a modification will interfere with
// the fastNLO routines.
//
// This file contains the following routines:
//   inputfunc   (-> user edits)
//   psinput     (-> user edits)
//   initfunc    (don't touch)
//   userfunc    (-> user edits)
//   writetable  (don't touch)
//   end_of_even (don't touch)
//   phys_output (don't touch)
//   inittable   (-> user edits)
//
// Implementing a new scenario requires to edit:
//  - the jet algorithm ("#include" statement and assignement of "jetclus")
//  - number of jets (inputfunc)
//  - center-of-mass energy (psinput)
//  - compute observable, determine bin No. (userfunc)
//  - declare all variables for table, define bin boundaries (inittable)
//  
// ================================================================
// 
// last modifications
// 2010/09/24 MW make user-friendly
// 2010/09/22 MW implement D0 phase space
// 2009/01/15 TK make code V2.0 compatible
//


//------ DON'T TOUCH THIS PART! ------
#include <bits/hhc-phasespace.h>
#include <bits/hhc-process.h>
#include <bits/hhc-jetfunc.h>

//----- used namespaces -----
using namespace nlo;
using namespace std;

//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
void psinput(phasespace_hhc *, double&);
user_base_hhc * userfunc();

//----- array of the symbols symbols -----
extern "C"{
  struct { 
	const char *name;
	void *address;
  } user_defined_functions[] = 
  {
	//   process index: hhc for hadron-hadron --> jets
	{"procindex", (void *) "hhc"},
	//   input function 
	{"inputfunc", (void *) inputfunc},
	//   phase space input function 
	{"psinput", (void *) psinput},
	//   user defined functions
	{"userfunc",  (void *) userfunc},
	//  end of the list
	{0, 0}
  };
}
//------ END OF THE DO-NOT-TOUCH-PART ------

//------ USER DEFINED PART STARTS HERE ------
#include <algorithm>

#include "fj-sc.h"   // fastNLO user: .h file for jet algorithm

#include "pdf-cteq6.h"
#include "pdf-hhc-dummy.h"
#include "fnloTable.h"
#include "fnloBlockBNlojet.h"

class UserHHC : public basic_user_set<user0d_hhc, user1h_hhc, user2h_hhc>
{
 public:
   //   init and user function
   void initfunc(unsigned int);
   void userfunc(const event_hhc&, const amplitude_hhc&);
   virtual void end_of_event();  
   virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);

 private:
   // --- pdf
   pdf_cteq6 pdf;
   pdf_hhc_dummy dummypdf;

   // --- jet algorithm
   //cone_e_07 jetclus;   // fastNLO user: define jet algorithm (consistent with .h file above)
   fj_sc jetclus;

   bounded_vector<lorentzvector<double> > pj;    // the jet structure 
   
   // --- fastNLO definitions (not for user)
   fnloTable *table;
   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table
   string tablefilename;  // The table file to write to
   time_t start_time;
   
   bool nlo;
   void inittable();
   void writetable();
};

user_base_hhc * userfunc() {
  return new UserHHC;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
  // --- fastNLO user: select the number of jets for your observable
  //nj = 1U;
  //nj = 2U;    
  nj = 3U;

   // --- number of the up and down type flavours (don't touch)
   nu = 2U;
   nd = 3U;
} 


void psinput(phasespace_hhc *ps, double& s)
{
   // --- fastNLO user: set the total c.m. energy square in GeV
   //s = 40000;       // RHIC           200GeV   
   //s = 3240000;     // TeV Run I     1800GeV
   s = 3841600;       // TeV Run II    1960GeV
   //s = 49000000     // LHC           7000GeV
   //s = 100000000    // LHC          10000GeV
   //s = 196000000;   // LHC          14000GeV
  
   //   You can use your own phase generator. 
   //   Here we use the default.
   ps = 0;
} 


void UserHHC::initfunc(unsigned int)
{
   // --- Initialize event counters
   nevents = 0;
   // Set some defaults
   if (nwrite==0) nwrite = 5000000;
   start_time = std::time(0);
}


void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
   fnloBlockA2 *A2 =  table->GetBlockA2();
   double x1 = p[-1].Z()/p[hadron(-1)].Z();
   double x2 = p[0].Z()/p[hadron(0)].Z();

   static const double twopi = 6.28318530717958647692;
   static const double pi = 3.14159265358979323846;

   // --- run the jet algorithm
   double jetsize = 0.7;
   pj = jetclus(p,jetsize);
   unsigned int nj = pj.upper(); 


   // --- fastNLO user:
   //     Here is your playground where you compute your observable 
   //     and the bin number ("obsbin") which gets passed to
   //     fastNLO's table filling code.
   //     (all pT and E are in GeV)

   // --- declare and initialize phase space cut variables
   double pTcut = 150., ymax = 2.4;

   // - require >=3 jets for trijet cross section
   if (nj >= 3) {   
     // - sort jets in descending pT
     int jpnt[6] = {0,1,2,3,4,5};
     for(unsigned int i=1;i<nj;i++){   
       for(unsigned int j=i+1;j<=nj;j++){
         double pt1 = pj[jpnt[i]].perp();
         double pt2 = pj[jpnt[j]].perp();
         if (pt2 > pt1) { // swap pointers
	   //	   printf("jets swapped %d %d %f %f\n",i,j,pt1,pt2);
           int itmp = jpnt[i];
           jpnt[i] = jpnt[j];
           jpnt[j] = itmp;        
         }
       }
     }
     
     // - compute observable and cut variables
     double ptmax = pj[jpnt[1]].perp(); 
     double pt2   = pj[jpnt[2]].perp();
     double ptmin = pj[jpnt[3]].perp(); 
     
     double y1 = pj[jpnt[1]].rapidity();
     double rap1 = fabs(y1);
     double phi1 = pj[jpnt[1]].phi();
     
     double y2 = pj[jpnt[2]].rapidity();
     double rap2 = fabs(y2);
     double phi2 = pj[jpnt[2]].phi();
     
     double y3 = pj[jpnt[3]].rapidity();
     double rap3 = fabs(y3);
     double phi3 = pj[jpnt[3]].phi();
     
     double dphi = fabs(phi1-phi2);
     if (dphi > pi) dphi = twopi - dphi;
     double deltar12 = sqrt((y1-y2)*(y1-y2) + dphi*dphi); 
     dphi = fabs(phi3-phi2);
     if (dphi > pi) dphi = twopi - dphi;
     double deltar23 = sqrt((y3-y2)*(y3-y2) + dphi*dphi); 
     dphi = fabs(phi1-phi3);
     if (dphi > pi) dphi = twopi - dphi;
     double deltar13 = sqrt((y1-y3)*(y1-y3) + dphi*dphi); 
     double deltar = min(deltar12, min(deltar13,deltar23));
     double rapmax = max(rap1,max(rap2,rap3));
     
     double m3jet = 0.001 *  // M3jet in TeV
       sqrt( (pj[jpnt[1]].T()+pj[jpnt[2]].T()+pj[jpnt[3]].T())*(pj[jpnt[1]].T()+pj[jpnt[2]].T()+pj[jpnt[3]].T())
	     -(pj[jpnt[1]].X()+pj[jpnt[2]].X()+pj[jpnt[3]].X())*(pj[jpnt[1]].X()+pj[jpnt[2]].X()+pj[jpnt[3]].X())
	     -(pj[jpnt[1]].Y()+pj[jpnt[2]].Y()+pj[jpnt[3]].Y())*(pj[jpnt[1]].Y()+pj[jpnt[2]].Y()+pj[jpnt[3]].Y())
	    -(pj[jpnt[1]].Z()+pj[jpnt[2]].Z()+pj[jpnt[3]].Z())*(pj[jpnt[1]].Z()+pj[jpnt[2]].Z()+pj[jpnt[3]].Z()));
     
     if(ptmax >= pTcut && ptmin >= 40. && m3jet > 0.4 && rapmax < ymax && deltar>1.4) {
       
       // --- identify bin number for M3jet - between 1 and 10, later check y,pt3 range(s)
       int m3jetbin = -1;
       //for(int j = 0; j < A2->GetNObsBin(); j++) {
       for(int j = 0; j < 10; j++) {
	 if (m3jet >= A2->LoBin[j][0] && m3jet < A2->UpBin[j][0]) {
	   m3jetbin = j;
	   break;
	 }
       }
       if (m3jetbin >= 0) {
	 // - set the renormalization and factorization scale to pTmax
	 //double mu = pt1; 
	 // - set the renormalization and factorization scale to average trijet pT
	 double mu = (ptmax + pt2 + ptmin)/3.; 
	 
	 double prefactor = 1./A2->BinSize[m3jetbin]; // - divide by binwidth
	 
	 // --- contribution to range 3 (pTmin>40, y<2.4)
	 int obsbin = m3jetbin + 20;
	 // --- fill fastNLO arrays - don't touch this piece of code!
	 for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
	   if(table->GetBlockB(k)->GetIRef()>0){
	     // --- fastNLO user: don't modify the following calls!
	     ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,pdf,prefactor);
	   }else{
	     ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,dummypdf,prefactor);
	   }
	 } // - end: fill fastNLO array
	 
	 // --- contribution to range 2 (pTmin>40, y<1.6)
	 if (ymax < 1.6) {
	   int obsbin = m3jetbin + 10;
	   // --- fill fastNLO arrays - don't touch this piece of code!
	   for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
	     if(table->GetBlockB(k)->GetIRef()>0){
	       // --- fastNLO user: don't modify the following calls!
	       ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,pdf,prefactor);
	     }else{
	       ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,dummypdf,prefactor);
	     }
	   } // - end: fill fastNLO array	  
	 } // end - range 2        

	 // --- contribution to range 1 (pTmin>40, y<0.8)
	 if (ymax < 0.8) {
	   int obsbin = m3jetbin;
	   // --- fill fastNLO arrays - don't touch this piece of code!
	   for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
	     if(table->GetBlockB(k)->GetIRef()>0){
	       // --- fastNLO user: don't modify the following calls!
	       ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,pdf,prefactor);
	     }else{
	       ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,dummypdf,prefactor);
	     }
	   } // - end: fill fastNLO array	  
	 } // end - range 1        

	 // --- contribution to range 4 (pTmin>70, y<2.4)
	 if (ptmin > 70.) {
	   int obsbin = m3jetbin;
	   // --- fill fastNLO arrays - don't touch this piece of code!
	   for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
	     if(table->GetBlockB(k)->GetIRef()>0){
	       // --- fastNLO user: don't modify the following calls!
	       ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,pdf,prefactor);
	     }else{
	       ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,dummypdf,prefactor);
	     }
	   } // - end: fill fastNLO array	  
	 } // end - range 4

	 // --- contribution to range 5 (pTmin>100, y<2.4)
	 if (ptmin > 100.) {
	   int obsbin = m3jetbin;
	   // --- fill fastNLO arrays - don't touch this piece of code!
	   for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
	     if(table->GetBlockB(k)->GetIRef()>0){
	       // --- fastNLO user: don't modify the following calls!
	       ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,pdf,prefactor);
	     }else{
	       ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,dummypdf,prefactor);
	     }
	   } // - end: fill fastNLO array	  
	 } // end - range 5

       } // - end: in M3jet histogram range
     } // - end: in analysis 3-jet phase space
     
   } // - end: reqire >=3 jets
   // --- end: fastNLO user playground
   
}

void UserHHC::writetable(){
   table->OpenFileRewrite();
   table->WriteBlockA1();
   table->WriteBlockA2();
   for(int i=0;i< table->GetBlockA1()->GetNcontrib();i++){
      table->WriteBlockBDividebyN(i);
   }
   table->CloseFileWrite();

}

void UserHHC::end_of_event(){
   nevents += 1;
   // --- store table
   if (( (unsigned long)nevents % nwrite)==0){
      time_t hour, min, time = std::time(0) - start_time;
      
      hour = time/3600L;
      time -= hour*3600L;
      min  = time/60L;
      time -= min*60L;
      
      std::cout<<"--->     "
	  <<(hour < 10 ? "0" : "")<<hour
	  <<(min < 10 ? ":0" : ":")<<min
	  <<(time < 10 ? ":0" : ":")<<time<<std::endl;
      printf ("No. events: %.2G writing table....",nevents);
      cout.flush();
      for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
         table->GetBlockB(k)->Nevt = (long long int)nevents;
      }
      writetable();
      printf("done.\n");
   }
}

void UserHHC::phys_output(const std::basic_string<char>& __file_name, 
                          unsigned long __save, bool __txt) 
{
   tablefilename.assign(__file_name.c_str());
   tablefilename += ".tab";
   
   // --- determine whether we are running LO or NLO
   const char* const file = __file_name.c_str(); 
   if(strstr(file,"born")!=NULL){
      nlo = false;
   }else{
      if(strstr(file,"nlo")!=NULL){
         nlo = true;
      }else{
         printf("This module can only be run at Born level or at NLO.\n");
         exit(1);
      }
   }
   nwrite = __save;
   inittable();
}

void UserHHC::inittable(){

   // --- fastNLO user: set the square of the center of mass energy (in GeV)
   //double s = 3240000; // Tevatron RunI
   double s = 3841600; // Tevatron RunII
   //double s = 49000000; // LHC @7TeV 
   //double s = 100000000; // LHC @10TeV
 
   // --- fastNLO user: decide whether to include a reference table (for 
   //                   precision studies, not for production jobs)
   //const bool doReference = true;
   const bool doReference = false;

   // --- set up fastNLO
   table = new fnloTable(tablefilename);

   // --- fastNLO: fill variable for table header block A1
   table->GetBlockA1()->SetScenName("fnt2013denom");  // - fastNLO user: set scenario name
   table->GetBlockA1()->SetNcontrib(1);
   table->GetBlockA1()->SetNmult(0);
   table->GetBlockA1()->SetNdata(0);
   table->GetBlockA2()->SetIpublunits(12);  // - fastNLO user: set cross section units
                                            //                 (negative power of ten)

   // --- fastNLO: fill variables for table header block A2
   fnloBlockA2 *A2 =  table->GetBlockA2();

   // --- fastNLO user: up to 30 strings to describe the scenario
   A2->ScDescript.push_back("dsigma-jet_dpTmax_(pb_GeV)");
   A2->ScDescript.push_back("xxx");
   A2->ScDescript.push_back("D0_Collaboration");
   //A2->ScDescript.push_back("...");
   //A2->ScDescript.push_back("...");
   A2->NScDescript = A2->ScDescript.size();
   A2->Ecms = sqrt(s);
   // A2->ILOord = 2;   // --- fastNLO user: power of LO contribution for process (dijet:2)
   A2->ILOord = 3;   // --- fastNLO user: power of LO contribution for process (trijet:3)
   A2->NDim = 2;     // --- fastNLO user: No of dimensions in which observable is binned
   A2->DimLabel.push_back("pTmax");  // --- fastNLO user: label of 1st dimension
   A2->IDiffBin.push_back(2);
   A2->DimLabel.push_back("pTmin");   // --- fastNLO user: label of 2nd dimension
   A2->IDiffBin.push_back(2);

   vector <double> bound;
   bound.resize(2);

   // --- fastNLO user: bin definitions 
   //  MW:  here rap->rap/pTmin  and  pT->M3-jet(TeV)
   const int nrapbins = 5;
   double rapbins[nrapbins+1] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

   const int nptbins[nrapbins] = {10, 10, 10, 10, 10};
   double ptbins[5][11] = {
     {0.40, 0.45, 0.50, 0.55, 0.61, 0.67, 0.74, 0.81, 0.90, 1.1, 1.5},
     {0.40, 0.45, 0.50, 0.55, 0.61, 0.67, 0.74, 0.81, 0.90, 1.1, 1.5},
     {0.40, 0.45, 0.50, 0.55, 0.61, 0.67, 0.74, 0.81, 0.90, 1.1, 1.5},
     {0.40, 0.45, 0.50, 0.55, 0.61, 0.67, 0.74, 0.81, 0.90, 1.1, 1.5},
     {0.40, 0.45, 0.50, 0.55, 0.61, 0.67, 0.74, 0.81, 0.90, 1.1, 1.5}
   };

   // --- fastNLO user:
   //     define below the bin width ("binsize") by which
   //     the cross section is divided to obtain the 
   //     (multi-) differential result.
   double binsize = 0.;

   int nbins = 0;   // --- count total No. bins
   for(int i=0;i<nrapbins;i++){
      for(int j=0;j<nptbins[i];j++){
 	 nbins += 1;
         bound[0] = ptbins[i][j];
         bound[1] = rapbins[i];
         A2->LoBin.push_back(bound);
         bound[0] = ptbins[i][j+1];
	 bound[1] = rapbins[i+1];
         A2->UpBin.push_back(bound);
	 //printf(" %d %d  |  %f %f\n",i,j,bound[0],bound[1]);

	 binsize = (ptbins[i][j+1]-ptbins[i][j]); // here: only width in M3-jet(TeV)
	 A2->BinSize.push_back(binsize);
      }
   }
   printf(" tot. No. observable bins = %d\n",nbins);

   A2->NObsBin = nbins;
   A2->INormFlag = 0;   // --- fastNLO user: default=0 - set =1 if observable is 
                        //     to be normalized by own integral (in 1st dimension)
                        //     see documentation for details and for other options

   // --- fastNLO table block B
   fnloBlockB *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
   table->CreateBlockB(0,B);
   B->IXsectUnits = 12;    // --- fastNLO user: set to same value as "SetIpublunits"
   B->IDataFlag = 0;
   B->IAddMultFlag = 0;
   B->IContrFlag1 = 1;
   B->IContrFlag3 = 0;
   B->CodeDescript.push_back("NLOJET++ 4.1.3");  // --- fastNLO user: enter NLOJET++ version
   B->NCodeDescr = B->CodeDescript.size();
   B->IRef = 0;
   if (nlo || A2->ILOord > 2) {
     B->NSubproc = 7;
   } else {
     B->NSubproc = 6;
     printf("  this job uses 6 subprocesses \n");
   }
   if(nlo){
      B->CtrbDescript.push_back("NLO");
      B->IContrFlag2 = 2;
      B->IScaleDep = 1;
      B->Npow = A2->ILOord+1;
   }else{
      B->CtrbDescript.push_back("LO");      
      B->IContrFlag2 = 1;
      B->IScaleDep = 0;
      B->Npow = A2->ILOord;
   }
   B->NContrDescr = B->CtrbDescript.size();

   B->NPDF = 2;
   B->NPDFPDG.push_back(2212);   // --- fastNLO user: PDG code for 1st hadron
   B->NPDFPDG.push_back(-2212);  // --- fastNLO user: PDG code for 2nd hadron
   B->NPDFDim = 1;
   B->NFragFunc = 0;
   B->IPDFdef1 = 3;
   B->IPDFdef2 = 1;
   if(B->NSubproc == 7) {
     B->IPDFdef3 = 2;
   } else {
     B->IPDFdef3 = 1;
     printf("  set IPDFdef3=1 consistent with 6 subprocesses \n");
   }
   B->XNode1.resize(A2->NObsBin);

   // --- initialize variables for WarmUp run
   B->IWarmUp = 1;       // --- fastNLO user: do the Warm-Up run
   //B->IWarmUp = 0;     //                   or do production run(s)
   // - fastNLO user: remember to disable reference-mode in
   //                 Warm-Up run: "doReference = false" (above)
   B->IWarmUpPrint = 5000000;
   B->xlo.resize(A2->NObsBin);
   B->scalelo.resize(A2->NObsBin);
   B->scalehi.resize(A2->NObsBin);

   // --- arrays for extreme x and (default) scale values (computed in Warm-Up run)
   double xlim[A2->NObsBin];
   double mulo[A2->NObsBin];
   double muup[A2->NObsBin];
   for(int i=0;i<A2->NObsBin;i++){
     xlim[i] = 1.1e-07, mulo[i] = 3.0, muup[i]=9.9e10; // - safe initializations
   }
   

   // --- fastNLO user: before running a new scenario for the first time,
   //     the following block (between "start" and "end") should be 
   //     completely removed. The first run must be a "Warm-Up Run" (IWarmUp=1)
   //     which produces an initialization block as output. This should be copied
   //     and pasted below. These initialization values must be used for all
   //     production jobs (IWarmUp=0) for a given scenario (otherwise the result 
   //     tables can not be merged). 
   //
   // --------- fastNLO: Warm-Up run results (start)

   // --------- fastNLO: Warm-Up run results (end)


   //printf("* --- xlimits \n");
   for(int i=0;i<A2->NObsBin;i++){
      int nxtot = 12;
      if (i == ((A2->NObsBin)-1)) nxtot = 13; // Darf's etwas mehr sein?
      B->Nxtot1.push_back(nxtot);
 
      // ---- safety factors for ET-scheme / optimized by eta range
      //double hxlim = -sqrt(-log10(xmin*0.81));   // pre Warm-Up
      double hxlim = -sqrt(-log10(xlim[i]));   // use exact value from Warm-Up
      //printf("%d %g %g \n",i,pow(10,-pow(hxlim,2)),xlim[i]);
      B->Hxlim1.push_back(hxlim);
      for(int j=0;j<nxtot;j++){
         double hx = hxlim*( 1.- ((double)j)/(double)nxtot);
         B->XNode1[i].push_back(pow(10,-pow(hx,2))); 
      }
   }

   B->NScales = 2;  // two scales: mur and muf
   B->NScaleDim = 1; // one variable used in scales: pT
   B->Iscale.push_back(0);  // mur=mur(pT), pT = index 0 
   B->Iscale.push_back(0);  // muf=muf(pT), pT = index 0 
   B->ScaleDescript.resize(B->NScaleDim);

   B->ScaleDescript[0].push_back("pT");
   B->NscaleDescript.push_back(B->ScaleDescript[0].size());
   B->Nscalenode.push_back(4); // number of scale nodes for pT

   B->ScaleFac.resize(B->NScaleDim);

   B->ScaleFac[0].push_back(1.0);    // --- fastNLO: central scale (don't change)
   if(nlo){
     B->ScaleFac[0].push_back(0.5);  // --- fastNLO user: add any number of
     B->ScaleFac[0].push_back(2.0);  //             additional scale variations
     B->ScaleFac[0].push_back(0.25); //             as desired.
     //B->ScaleFac[0].push_back(4.0);
     //B->ScaleFac[0].push_back(8.0);
   }
   B->Nscalevar.push_back(B->ScaleFac[0].size());

   //printf("* --- scales (bin,dimension,variation)\n");
   B->ScaleNode.resize(A2->NObsBin);
   for(int i=0;i<A2->NObsBin;i++){
      B->ScaleNode[i].resize(B->NScaleDim);
      for(int j=0;j<B->NScaleDim;j++){
         B->ScaleNode[i][j].resize(B->Nscalevar[j]);
         for(int k=0;k<B->Nscalevar[j];k++){
            B->ScaleNode[i][j][k].resize(B->Nscalenode[j]);
            if(B->Nscalenode[j]==1){
	      B->ScaleNode[i][j][k][0] = B->ScaleFac[0][k] * (muup[i]+mulo[i])/2.; // assume only one scale dimension
	      // B->ScaleNode[i][j][k][0] = B->ScaleFac[0][k]*(A2->UpBin[i][B->Iscale[j]] + A2->LoBin[i][B->Iscale[j]])/2.;// pre Warm-Up
            }else{
               const double mu0scale = .25; // In GeV
               double llscalelo = log(log((B->ScaleFac[0][k]*mulo[i]*0.98)/mu0scale));  // safety factors -> remove when
               double llscalehi = log(log((B->ScaleFac[0][k]*muup[i]*1.02)/mu0scale));  //     interpolation is improved 
	       //double llscalelo = log(log((B->ScaleFac[0][k]*A2->LoBin[i][B->Iscale[j]])/mu0scale)); // pre Warm-Up
               //double llscalehi = log(log((B->ScaleFac[0][k]*A2->UpBin[i][B->Iscale[j]])/mu0scale));// pre Warm-Up
	       //printf("%d %d %d %g %g\n",i,j,k,mu0scale*exp(exp(llscalelo)),mu0scale*exp(exp(llscalehi)));
               for(int l=0;l<B->Nscalenode[j];l++){
		 // later this is the place where the Chebychev nodes will be implemented
		 B->ScaleNode[i][j][k][l] = mu0scale * exp(exp( llscalelo +
								double(l)/double(B->Nscalenode[j]-1)*(llscalehi-llscalelo) ));
               }
            }
         }            
      }
   }

   B->SigmaTilde.resize(A2->NObsBin);
   for(int i=0;i<A2->NObsBin;i++){
      B->SigmaTilde[i].resize(B->Nscalevar[0]);
      for(int k=0;k<B->Nscalevar[0];k++){
         B->SigmaTilde[i][k].resize(B->Nscalenode[0]);
         for(int l=0;l<B->Nscalenode[0];l++){
            int nxmax = B->GetNxmax(i);
            B->SigmaTilde[i][k][l].resize(nxmax);
            for(int m=0;m<nxmax;m++){
               B->SigmaTilde[i][k][l][m].resize(B->NSubproc);
               for(int n=0;n<B->NSubproc;n++){
                  B->SigmaTilde[i][k][l][m][n] = 0.;
               }
            }
         }            
      }
   }
   
   
   // reference table
   if(doReference){
      fnloBlockB *refB = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
      //refB->NSubproc = 7;
      if (nlo || A2->ILOord > 2) {
        refB->NSubproc = 7;
      } else {
        refB->NSubproc = 6;
        printf("  this reference job uses 6 subprocesses \n");
      }
      table->CreateBlockB(1,refB);
      refB->Copy(table->GetBlockB(0));
      refB->IRef = 1;
      refB->Nscalenode[0] = 1;
      refB->Nxtot1.clear();
      refB->Hxlim1.clear();
      for(int i=0;i<A2->NObsBin;i++){
         refB->Nxtot1.push_back(1);
         refB->Hxlim1.push_back(0.);
         refB->XNode1[i].clear(); 
         refB->XNode1[i].push_back(0.); 
      }
      table->GetBlockA1()->SetNcontrib(2);
   }

}
