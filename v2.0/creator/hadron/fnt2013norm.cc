//
//   fastNLO v2 author code 
//   scenario: fnt2013(norm)
//   denominator for Run IIa D0 R3/2 measurement: sigma-dijet
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

#include "cone-e.h"     // fastNLO user: .h file for jet algorithm

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
   cone_e jetclus;   // fastNLO user: define jet algorithm (consistent with .h file above)
   
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
  nj = 2U;    
  //nj = 3U;

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
   double pTcut = 80., ymax = 2.4;

   // - require >=2 jets for dijet cross section
   if (nj >= 2) {   
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
     double pt1 = pj[jpnt[1]].perp(); 
     double y1 = pj[jpnt[1]].rapidity();
     double rap1 = fabs(y1);
     double phi1 = pj[jpnt[1]].phi();
     
     double pt2 = pj[jpnt[2]].perp(); 
     double y2 = pj[jpnt[2]].rapidity();
     double rap2 = fabs(y2);
     double phi2 = pj[jpnt[2]].phi();
     
     double dphi = fabs(phi1-phi2);
     if (dphi > pi) dphi = twopi - dphi;
     double deltar = sqrt((y1-y2)*(y1-y2) + dphi*dphi); 
     
     if(pt1 >= pTcut && rap1 < ymax && rap2 < ymax && pt2 >= 50. && deltar>1.4) {
      
       // - set the renormalization and factorization scale to pTmax
       double mu = pt1; 

       // --- identify bin number for pTmax for pTmin>50
       int obsbin = -1;
       for(int j = 0; j < 12; j++) {
	 if (pt1 >= A2->LoBin[j][0] && pt1 < A2->UpBin[j][0]) {
	   obsbin = j;
	   break;
	 }
       }
       // --- fill fastNLO arrays - don't touch this piece of code!
       if (obsbin >= 0) {
	 double prefactor = 1./A2->BinSize[obsbin]; // - divide by binwidth
	 // --- fastNLO: don't touch the following calls!
	 for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
	   if(table->GetBlockB(k)->GetIRef()>0){
	     ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,pdf,prefactor);
	   }else{
	     ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,dummypdf,prefactor);
	   }
	 } 

	 if (pt2>=70.0 && pt1>=100.0) {
	   obsbin+=11;
	   prefactor = 1./A2->BinSize[obsbin]; // - divide by binwidth
	   // --- fastNLO: don't touch the following calls!
	   for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
	     if(table->GetBlockB(k)->GetIRef()>0){
	       ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,pdf,prefactor);
	     }else{
	       ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,dummypdf,prefactor);
	     }
	   }
	 
	   if (pt2>=90.0 && pt1>=120.0) {
	     obsbin+=10;
	     prefactor = 1./A2->BinSize[obsbin]; // - divide by binwidth
	     // --- fastNLO: don't touch the following calls!
	     for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
	       if(table->GetBlockB(k)->GetIRef()>0){
		 // --- fastNLO user: don't modify the following calls!
		 ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,pdf,prefactor);
	       }else{
		 ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,dummypdf,prefactor);
	       }
	     } 
	   } // pT2>90
	 } // pT2>70
       } // - end: fill fastNLO array  (pT2>50)       
     } // - end: in analysis 2-jet phase space     
   } // - end: reqire >=2 jets
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
   const bool doReference = true;
   //const bool doReference = false;

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
   A2->ScDescript.push_back("dsigma-dijet_dpTmax_(pb_GeV)");
   A2->ScDescript.push_back("xxx");
   A2->ScDescript.push_back("D0_Collaboration");
   //A2->ScDescript.push_back("...");
   //A2->ScDescript.push_back("...");
   A2->NScDescript = A2->ScDescript.size();
   A2->Ecms = sqrt(s);
   A2->ILOord = 2;   // --- fastNLO user: power of LO contribution for process
   A2->NDim = 2;     // --- fastNLO user: No of dimensions in which observable is binned
   A2->DimLabel.push_back("pTmax");  // --- fastNLO user: label of 1st dimension
   A2->IDiffBin.push_back(2);
   A2->DimLabel.push_back("pTmin");   // --- fastNLO user: label of 2nd dimension
   A2->IDiffBin.push_back(2);

   vector <double> bound;
   bound.resize(2);

   // --- fastNLO user: bin definitions 
   //  MW:  here rap->pTmin  and  pT->pTmax
   const int nrapbins = 3;
   double rapbins[nrapbins+1] = {50.0, 70.0, 90.0, 980.0};

   const int nptbins[nrapbins] = {12, 11, 10};
   double ptbins[3][13] = {
     {80., 100., 120., 140., 165., 190., 220., 250., 285., 320., 
      360., 400., 500.},
     {100., 120., 140., 165., 190., 220., 250., 285., 320., 360., 
      400., 500.},
     {120., 140., 165., 190., 220., 250., 285., 320., 360., 400.,
      500.}
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
	 //bound[1] = rapbins[i+1];
         bound[1] = rapbins[3];   // upper bound for pTmin is always 980GeV(=no bound)
         A2->UpBin.push_back(bound);
	 //printf(" %d %d  |  %f %f\n",i,j,bound[0],bound[1]);

	 binsize = (ptbins[i][j+1]-ptbins[i][j]); // here: only width in pTmax
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
   //B->IWarmUp = 1;       // --- fastNLO user: do the Warm-Up run
   B->IWarmUp = 0;     //                   or do production run(s)
   // - fastNLO user: remember to disable reference-mode in
   //                 Warm-Up run: "doReference = false" (above)
   B->IWarmUpPrint = 50000000;
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
   // 900000000 contributions (!= events) in warm-up run 
   xlim[0]=0.00745, mulo[0]=80.0, muup[0]=100.0;
   xlim[1]=0.01042, mulo[1]=100.0, muup[1]=120.0;
   xlim[2]=0.01500, mulo[2]=120.0, muup[2]=140.0;
   xlim[3]=0.02042, mulo[3]=140.0, muup[3]=165.0;
   xlim[4]=0.02830, mulo[4]=165.0, muup[4]=190.0;
   xlim[5]=0.03761, mulo[5]=190.0, muup[5]=220.0;
   xlim[6]=0.05046, mulo[6]=220.0, muup[6]=250.0;
   xlim[7]=0.06514, mulo[7]=250.0, muup[7]=285.0;
   xlim[8]=0.08463, mulo[8]=285.0, muup[8]=320.0;
   xlim[9]=0.10679, mulo[9]=320.0, muup[9]=360.0;
   xlim[10]=0.13498, mulo[10]=360.0, muup[10]=400.0;
   xlim[11]=0.16664, mulo[11]=400.0, muup[11]=500.0;
   xlim[12]=0.01042, mulo[12]=100.0, muup[12]=120.0;
   xlim[13]=0.01503, mulo[13]=120.0, muup[13]=140.0;
   xlim[14]=0.02042, mulo[14]=140.0, muup[14]=165.0;
   xlim[15]=0.02835, mulo[15]=165.0, muup[15]=190.0;
   xlim[16]=0.03761, mulo[16]=190.0, muup[16]=220.0;
   xlim[17]=0.05046, mulo[17]=220.0, muup[17]=250.0;
   xlim[18]=0.06514, mulo[18]=250.0, muup[18]=285.0;
   xlim[19]=0.08463, mulo[19]=285.0, muup[19]=320.0;
   xlim[20]=0.10679, mulo[20]=320.0, muup[20]=360.0;
   xlim[21]=0.13498, mulo[21]=360.0, muup[21]=400.0;
   xlim[22]=0.16664, mulo[22]=400.0, muup[22]=500.0;
   xlim[23]=0.01503, mulo[23]=120.0, muup[23]=140.0;
   xlim[24]=0.02042, mulo[24]=140.0, muup[24]=165.0;
   xlim[25]=0.02835, mulo[25]=165.0, muup[25]=190.0;
   xlim[26]=0.03761, mulo[26]=190.0, muup[26]=220.0;
   xlim[27]=0.05046, mulo[27]=220.0, muup[27]=250.0;
   xlim[28]=0.06514, mulo[28]=250.0, muup[28]=285.0;
   xlim[29]=0.08463, mulo[29]=285.0, muup[29]=320.0;
   xlim[30]=0.10679, mulo[30]=320.0, muup[30]=360.0;
   xlim[31]=0.13498, mulo[31]=360.0, muup[31]=400.0;
   xlim[32]=0.16664, mulo[32]=400.0, muup[32]=500.0;
   // --------- fastNLO: Warm-Up run results (end)


   //printf("* --- xlimits \n");
   for(int i=0;i<A2->NObsBin;i++){
      //int nxtot = 12;
      int nxtot = 15;
      //if (i == ((A2->NObsBin)-1)) nxtot = 13; // Darf's etwas mehr sein?
      B->Nxtot1.push_back(nxtot);
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
   B->Nscalenode.push_back(6); // number of scale nodes for pT

   B->ScaleFac.resize(B->NScaleDim);

   B->ScaleFac[0].push_back(1.0);    // --- fastNLO: central scale (don't change)
   if(nlo){
     //B->ScaleFac[0].push_back(0.5);  // --- fastNLO user: add any number of
     //B->ScaleFac[0].push_back(2.0);  //             additional scale variations
     //B->ScaleFac[0].push_back(0.25); //             as desired.
     //B->ScaleFac[0].push_back(4.0);
     //B->ScaleFac[0].push_back(8.0);
   }
   B->Nscalevar.push_back(B->ScaleFac[0].size());

   const double mu0scale = .25; // In GeV
   B->ScaleNode.resize(A2->NObsBin);
   B->HScaleNode.resize(A2->NObsBin);
   for(int i=0;i<A2->NObsBin;i++){
     B->ScaleNode[i].resize(B->NScaleDim);
     B->HScaleNode[i].resize(B->NScaleDim);
     for(int j=0;j<B->NScaleDim;j++){
       B->ScaleNode[i][j].resize(B->Nscalevar[j]);
       B->HScaleNode[i][j].resize(B->Nscalevar[j]);
       for(int k=0;k<B->Nscalevar[j];k++){
	 B->ScaleNode[i][j][k].resize(B->Nscalenode[j]);
	 B->HScaleNode[i][j][k].resize(B->Nscalenode[j]);
	 if(B->Nscalenode[j]==1){
	   B->ScaleNode[i][j][k][0] = B->ScaleFac[0][k] * (muup[i]+mulo[i])/2.; // assume only one scale dimension
	   B->HScaleNode[i][j][k][0] = log(log((B->ScaleFac[0][k]*(muup[i]+mulo[i])/2.)/mu0scale));
	 }else{
	   double llscalelo = log(log((B->ScaleFac[0][k]*mulo[i])/mu0scale));  
	   double llscalehi = log(log((B->ScaleFac[0][k]*muup[i])/mu0scale));  
	   for(int l=0;l<B->Nscalenode[j];l++){
                 B->HScaleNode[i][j][k][l] = llscalelo +
                   double(l)/double(B->Nscalenode[j]-1)*(llscalehi-llscalelo);
                 B->ScaleNode[i][j][k][l] = mu0scale * exp(exp(B->HScaleNode[i][j][k][l]));
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
   
   
   // --- reference table
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
