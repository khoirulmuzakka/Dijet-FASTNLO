//
// fastNLO v2 author code for fnl2342a:
//     CMS LHC Inclusive Jets Scenario, E_cms = 7 TeV
//     for fastjet anti-kT algo with R=0.5 in E-scheme
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
//  - the jet algorithm ("#include" statement and assignment of "jetclus")
//  - number of jets (inputfunc)
//  - center-of-mass energy (psinput)
//  - compute observable, determine bin No. (userfunc)
//  - declare all variables for table, define bin boundaries (inittable)
//  
// ================================================================
// 
// last modifications
// 2011/01/26 KR Add CMS Inclusive Jets Scenario
// 2011/01/13 KR unify jet sizes into one .cc and .h file
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

#include "fj-ak.h"

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
   fj_ak jetclus;
   
   bounded_vector<lorentzvector<double> > pj;    // the jet structure 
   
   // --- fastNLO definitions (not for user)
   fnloTable *table;
   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table
   string tablefilename;  // The table file to write to
   time_t start_time;
   
   bool doReference;
   bool doWarmUp;
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
   // --- fastNLO user: set the total c.m. energy squared in GeV^2
   //s =     40000.; // RHIC               200 GeV
   //s =   3240000.; // TeV Run I         1800 GeV
   //s =   3841600.; // TeV Run II        1960 GeV
   //s =    810000.; // LHC Injection Run  900 GeV
   //s =   5569600.; // LHC Initial Run   2360 GeV
   s =  49000000.; // LHC First Run     7000 GeV
   //s = 100000000.; // LHC Start-up Run 10000 GeV
   //s = 196000000.; // LHC Design Run   14000 GeV

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

   // --- run the jet algorithm
   double jetsize = 0.5;
   pj = jetclus(p,jetsize);
   unsigned int nj = pj.upper(); 


   // --- fastNLO user:
   //     Here is your playground where you compute your observable 
   //     and the bin number ("obsbin") which gets passed to
   //     fastNLO's table filling code.
   //     (all pT and E are in GeV)

   // --- declare and initialize phase space cut variables
   // lowest pT for jets to be considered
   double ptmin = 10.;
   // highest (pseudo-)rapidity for jets to be considered
   double ymax  = 3.;

   // Analyze inclusive jets in jet loop
   for (unsigned int i = 1; i <= nj; i++) {
     double pt  = pj[i].perp(); 
     double rap = abs(pj[i].rapidity());
     
     // --- jet in phase space?
     if (ptmin < pt && rap < ymax) {

       // - set the renormalization and factorization scale to jet pT
       double mu = pt;

       // --- identify bin number (y,pT)
       int obsbin = -1;
       for (int j = 0; j < A2->GetNObsBin(); j++) {
	 if (A2->LoBin[j][0] <= pt  && pt  < A2->UpBin[j][0] && 
	     A2->LoBin[j][1] <= rap && rap < A2->UpBin[j][1]) {
	   obsbin = j;
	   break;
	 }
       }
	 
	 // --- fill fastNLO arrays - don't touch this piece of code!
	 if (obsbin >= 0) {
 	    double prefactor = 1./A2->BinSize[obsbin]; // - divide by binwidth
	    for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
   	       if(table->GetBlockB(k)->GetIRef()>0){
	          ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,pdf,prefactor);
	       }else{
	 	  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,mu,amp,dummypdf,prefactor);
	       }
	    }
	 } // - end: fill fastNLO array
      }
   }
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
      printf("fastNLO: This is a LO run!\n");
   }else{
      if(strstr(file,"nlo")!=NULL){
         nlo = true;
	 printf("fastNLO: This is a NLO run!\n");
      }else{
         printf("This module can only be run at Born level or at NLO.\n");
         exit(1);
      }
   }

   // --- determine whether this is a warm-up or reference job
   doWarmUp = false;
   if (strstr(file,"wrm")!=NULL) {
     doWarmUp = true;
     printf("fastNLO: This is a warm-up run!\n");
     if ( ! nlo ) {
       printf("fastNLO: WARNING! Warm-up runs are better done at NLO!\n");
     }
   }
   doReference = false;
   if (strstr(file,"ref")!=NULL) {
     doReference = true;
     printf("fastNLO: This is a reference run!\n");
   }
   if ( doWarmUp && doReference ) {
     printf("fastNLO: ERROR! Warm-up and reference runs cannot be done simultaneously:\n");
     printf("         doWarmUp = %i, doReference = %i\n",doWarmUp,doReference);
     exit(2);
   }

   nwrite = __save;
   inittable();
}

void UserHHC::inittable(){

  cout << "Run at NLO: " << nlo << endl; 

   // --- fastNLO user: set the total c.m. energy squared in GeV^2
   //double s =     40000.; // RHIC               200 GeV
   //double s =   3240000.; // TeV Run I         1800 GeV
   //double s =   3841600.; // TeV Run II        1960 GeV
   //double s =    810000.; // LHC Injection Run  900 GeV
   //double s =   5569600.; // LHC Initial Run   2360 GeV
   double s =  49000000.; // LHC First Run     7000 GeV
   //double s = 100000000.; // LHC Start-up Run 10000 GeV
   //double s = 196000000.; // LHC Design Run   14000 GeV

   // --- fastNLO user: decide whether to include a reference table (for 
   //                   precision studies, not for production jobs)
   // KR: Now set via filename string match to "ref"
   //const bool doReference = true;
   //const bool doReference = false;

   // --- set up fastNLO
   table = new fnloTable(tablefilename);

   // --- fastNLO: fill variable for table header block A1
   table->GetBlockA1()->SetScenName("fnl2342a");  // - fastNLO user: set scenario name
   table->GetBlockA1()->SetNcontrib(1);
   table->GetBlockA1()->SetNmult(0);
   table->GetBlockA1()->SetNdata(0);
   table->GetBlockA2()->SetIpublunits(15);  // - fastNLO user: set cross section units
                                            //                 (negative power of ten)

   // --- fastNLO: fill variables for table header block A2
   fnloBlockA2 *A2 =  table->GetBlockA2();

   // --- fastNLO user: up to 20 strings to describe the scenario
   A2->ScDescript.push_back("d2sigma-jet_dpT_d|y|_(fb_GeV)");
   A2->ScDescript.push_back("CMS_Collaboration");
   A2->ScDescript.push_back("CMS-PAP-QCD-10-011");
   //A2->ScDescript.push_back("");

   A2->NScDescript = A2->ScDescript.size();
   A2->Ecms = sqrt(s);
   A2->ILOord = 2;   // --- fastNLO user: power of LO contribution for process
   A2->NDim = 2;     // --- fastNLO user: No of dimensions in which observable is binned
   A2->DimLabel.push_back("pT");  // --- fastNLO user: label of 1st dimension
   A2->IDiffBin.push_back(2);
   A2->DimLabel.push_back("|y|");   // --- fastNLO user: label of 2nd dimension
   A2->IDiffBin.push_back(2);

   vector <double> bound;
   bound.resize(2);

   // --- fastNLO user: bin definitions - here in |y| and pT
   const int nrapbins = 6;
   double rapbins[nrapbins+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};

   const int nptbins[nrapbins] = {49, 49, 49, 46, 39, 32};
   double ptb[50] = { 10., 12., 15., 18., 21., 24., 28., 32., 37., 43., 49., 56., 64., 74., 84., 97.,
		      114., 133., 153., 174., 196., 220., 245., 272., 300., 330., 362., 395., 430., 468.,
		      507., 548., 592., 638., 686., 737., 790., 846., 905., 967.,
		      1032., 1101., 1172., 1248., 1327., 1410., 1497., 1588., 1684., 1784.};
   double ptbins[6][50];
   for (int i=0; i<nrapbins; i++) {
     for (int j=0; j<50; j++) { 
       if (j < nptbins[i]+1 ) {
	 ptbins[i][j] = ptb[j];
       } else {
	 ptbins[i][j] = 0.;
       }
     }
   }

   // --- fastNLO user:
   //     define below the bin width ("binsize") by which
   //     the cross section is divided to obtain the 
   //     (multi-) differential result.
   double binsize = 0.;

   int nbins = 0;   // --- count total No. bins
   for (int i=0;i<nrapbins;i++){
     for (int j=0;j<nptbins[i];j++){
         nbins += 1;
         bound[0] = ptbins[i][j];
         bound[1] = rapbins[i];
         A2->LoBin.push_back(bound);
         bound[0] = ptbins[i][j+1];
         bound[1] = rapbins[i+1];
         A2->UpBin.push_back(bound);
	 //printf(" %d %d  |  %f %f\n",i,j,bound[0],bound[1]);

	 binsize = (ptbins[i][j+1]-ptbins[i][j]) *
	   2.*(rapbins[i+1]-rapbins[i]);  // "*2." to account for |y| vs. y
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
   B->IXsectUnits = 15;    // --- fastNLO user: set to same value as "SetIpublunits"
   B->IDataFlag = 0;
   B->IAddMultFlag = 0;
   B->IContrFlag1 = 1;
   B->IContrFlag3 = 0;
   B->CodeDescript.push_back("NLOJet++ 4.1.3");  // --- fastNLO user: enter NLOJET++ version
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
   B->NPDFPDG.push_back(2212);  // --- fastNLO user: PDG code for 2nd hadron
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
   // KR: Now set via filename string match to "wrm"
   //B->IWarmUp = 1;     // --- fastNLO user: do the Warm-Up run
   B->IWarmUp = 0;     //                   or do production run(s)
   if ( doWarmUp ) {B->IWarmUp = 1;} 
   // KR: This is caught in an error condition now 
   // - fastNLO user: remember to disable reference-mode in
   //                 Warm-Up run: "doReference = false" (above)
   B->IWarmUpPrint = 100000000;
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
   if ( ! doWarmUp ) {
     // 200000000 contributions (!= events) in warm-up run 
     xlim[0]=0.00089, mulo[0]=  10.000, muup[0]=  12.000;
     xlim[1]=0.00106, mulo[1]=  12.000, muup[1]=  15.000;
     xlim[2]=0.00134, mulo[2]=  15.000, muup[2]=  18.000;
     xlim[3]=0.00163, mulo[3]=  18.000, muup[3]=  21.000;
     xlim[4]=0.00189, mulo[4]=  21.000, muup[4]=  24.000;
     xlim[5]=0.00214, mulo[5]=  24.000, muup[5]=  28.000;
     xlim[6]=0.00249, mulo[6]=  28.000, muup[6]=  32.000;
     xlim[7]=0.00290, mulo[7]=  32.000, muup[7]=  37.000;
     xlim[8]=0.00331, mulo[8]=  37.000, muup[8]=  43.000;
     xlim[9]=0.00387, mulo[9]=  43.000, muup[9]=  49.000;
     xlim[10]=0.00434, mulo[10]=  49.000, muup[10]=  56.000;
     xlim[11]=0.00500, mulo[11]=  56.000, muup[11]=  64.000;
     xlim[12]=0.00571, mulo[12]=  64.000, muup[12]=  74.000;
     xlim[13]=0.00659, mulo[13]=  74.000, muup[13]=  84.000;
     xlim[14]=0.00763, mulo[14]=  84.000, muup[14]=  97.000;
     xlim[15]=0.00884, mulo[15]=  97.000, muup[15]= 114.000;
     xlim[16]=0.01019, mulo[16]= 114.000, muup[16]= 133.000;
     xlim[17]=0.01195, mulo[17]= 133.000, muup[17]= 153.000;
     xlim[18]=0.01392, mulo[18]= 153.000, muup[18]= 174.000;
     xlim[19]=0.01579, mulo[19]= 174.000, muup[19]= 196.000;
     xlim[20]=0.01793, mulo[20]= 196.000, muup[20]= 220.000;
     xlim[21]=0.02030, mulo[21]= 220.000, muup[21]= 245.000;
     xlim[22]=0.02263, mulo[22]= 245.000, muup[22]= 272.000;
     xlim[23]=0.02532, mulo[23]= 272.000, muup[23]= 300.000;
     xlim[24]=0.02813, mulo[24]= 300.000, muup[24]= 330.000;
     xlim[25]=0.03132, mulo[25]= 330.000, muup[25]= 362.000;
     xlim[26]=0.03468, mulo[26]= 362.000, muup[26]= 395.000;
     xlim[27]=0.03791, mulo[27]= 395.000, muup[27]= 430.000;
     xlim[28]=0.04176, mulo[28]= 430.000, muup[28]= 468.000;
     xlim[29]=0.04596, mulo[29]= 468.000, muup[29]= 507.000;
     xlim[30]=0.05003, mulo[30]= 507.000, muup[30]= 548.000;
     xlim[31]=0.05505, mulo[31]= 548.000, muup[31]= 592.000;
     xlim[32]=0.05991, mulo[32]= 592.000, muup[32]= 638.000;
     xlim[33]=0.06569, mulo[33]= 638.000, muup[33]= 686.000;
     xlim[34]=0.07160, mulo[34]= 686.000, muup[34]= 736.999;
     xlim[35]=0.07764, mulo[35]= 737.000, muup[35]= 790.000;
     xlim[36]=0.08470, mulo[36]= 790.000, muup[36]= 846.000;
     xlim[37]=0.09206, mulo[37]= 846.000, muup[37]= 905.000;
     xlim[38]=0.10063, mulo[38]= 905.000, muup[38]= 967.000;
     xlim[39]=0.10917, mulo[39]= 967.000, muup[39]=1032.000;
     xlim[40]=0.11908, mulo[40]=1032.000, muup[40]=1101.000;
     xlim[41]=0.12949, mulo[41]=1101.000, muup[41]=1171.999;
     xlim[42]=0.14190, mulo[42]=1172.000, muup[42]=1248.000;
     xlim[43]=0.15418, mulo[43]=1248.000, muup[43]=1326.998;
     xlim[44]=0.16864, mulo[44]=1327.000, muup[44]=1410.000;
     xlim[45]=0.18473, mulo[45]=1410.001, muup[45]=1496.999;
     xlim[46]=0.20119, mulo[46]=1497.001, muup[46]=1588.000;
     xlim[47]=0.22049, mulo[47]=1588.000, muup[47]=1683.999;
     xlim[48]=0.24263, mulo[48]=1684.000, muup[48]=1784.000;
     xlim[49]=0.00054, mulo[49]=  10.000, muup[49]=  12.000;
     xlim[50]=0.00067, mulo[50]=  12.000, muup[50]=  15.000;
     xlim[51]=0.00081, mulo[51]=  15.000, muup[51]=  18.000;
     xlim[52]=0.00097, mulo[52]=  18.000, muup[52]=  21.000;
     xlim[53]=0.00114, mulo[53]=  21.000, muup[53]=  24.000;
     xlim[54]=0.00130, mulo[54]=  24.000, muup[54]=  28.000;
     xlim[55]=0.00152, mulo[55]=  28.000, muup[55]=  32.000;
     xlim[56]=0.00176, mulo[56]=  32.000, muup[56]=  37.000;
     xlim[57]=0.00205, mulo[57]=  37.000, muup[57]=  43.000;
     xlim[58]=0.00235, mulo[58]=  43.000, muup[58]=  49.000;
     xlim[59]=0.00266, mulo[59]=  49.000, muup[59]=  56.000;
     xlim[60]=0.00306, mulo[60]=  56.000, muup[60]=  64.000;
     xlim[61]=0.00356, mulo[61]=  64.000, muup[61]=  74.000;
     xlim[62]=0.00408, mulo[62]=  74.000, muup[62]=  84.000;
     xlim[63]=0.00475, mulo[63]=  84.000, muup[63]=  97.000;
     xlim[64]=0.00543, mulo[64]=  97.000, muup[64]= 114.000;
     xlim[65]=0.00646, mulo[65]= 114.000, muup[65]= 133.000;
     xlim[66]=0.00764, mulo[66]= 133.000, muup[66]= 153.000;
     xlim[67]=0.00876, mulo[67]= 153.000, muup[67]= 174.000;
     xlim[68]=0.00995, mulo[68]= 174.000, muup[68]= 196.000;
     xlim[69]=0.01121, mulo[69]= 196.000, muup[69]= 220.000;
     xlim[70]=0.01272, mulo[70]= 220.000, muup[70]= 245.000;
     xlim[71]=0.01438, mulo[71]= 245.000, muup[71]= 272.000;
     xlim[72]=0.01611, mulo[72]= 272.000, muup[72]= 300.000;
     xlim[73]=0.01803, mulo[73]= 300.000, muup[73]= 330.000;
     xlim[74]=0.02006, mulo[74]= 330.000, muup[74]= 362.000;
     xlim[75]=0.02238, mulo[75]= 362.000, muup[75]= 395.000;
     xlim[76]=0.02469, mulo[76]= 395.000, muup[76]= 430.000;
     xlim[77]=0.02743, mulo[77]= 430.000, muup[77]= 468.000;
     xlim[78]=0.03027, mulo[78]= 468.000, muup[78]= 507.000;
     xlim[79]=0.03331, mulo[79]= 507.000, muup[79]= 548.000;
     xlim[80]=0.03682, mulo[80]= 548.000, muup[80]= 592.000;
     xlim[81]=0.04065, mulo[81]= 592.000, muup[81]= 638.000;
     xlim[82]=0.04501, mulo[82]= 638.000, muup[82]= 686.000;
     xlim[83]=0.04954, mulo[83]= 686.000, muup[83]= 737.000;
     xlim[84]=0.05484, mulo[84]= 737.000, muup[84]= 790.000;
     xlim[85]=0.06031, mulo[85]= 790.000, muup[85]= 846.000;
     xlim[86]=0.06665, mulo[86]= 846.000, muup[86]= 905.000;
     xlim[87]=0.07401, mulo[87]= 905.000, muup[87]= 967.000;
     xlim[88]=0.08209, mulo[88]= 967.000, muup[88]=1032.000;
     xlim[89]=0.09149, mulo[89]=1032.000, muup[89]=1101.000;
     xlim[90]=0.10184, mulo[90]=1101.000, muup[90]=1171.999;
     xlim[91]=0.11356, mulo[91]=1172.000, muup[91]=1248.000;
     xlim[92]=0.12816, mulo[92]=1248.001, muup[92]=1327.000;
     xlim[93]=0.14438, mulo[93]=1327.000, muup[93]=1410.000;
     xlim[94]=0.16305, mulo[94]=1410.000, muup[94]=1496.999;
     xlim[95]=0.18370, mulo[95]=1497.002, muup[95]=1588.000;
     xlim[96]=0.20650, mulo[96]=1588.000, muup[96]=1684.000;
     xlim[97]=0.23256, mulo[97]=1684.001, muup[97]=1784.000;
     xlim[98]=0.00033, mulo[98]=  10.000, muup[98]=  12.000;
     xlim[99]=0.00040, mulo[99]=  12.000, muup[99]=  15.000;
     xlim[100]=0.00048, mulo[100]=  15.000, muup[100]=  18.000;
     xlim[101]=0.00060, mulo[101]=  18.000, muup[101]=  21.000;
     xlim[102]=0.00070, mulo[102]=  21.000, muup[102]=  24.000;
     xlim[103]=0.00080, mulo[103]=  24.000, muup[103]=  28.000;
     xlim[104]=0.00094, mulo[104]=  28.000, muup[104]=  32.000;
     xlim[105]=0.00106, mulo[105]=  32.000, muup[105]=  37.000;
     xlim[106]=0.00125, mulo[106]=  37.000, muup[106]=  43.000;
     xlim[107]=0.00146, mulo[107]=  43.000, muup[107]=  49.000;
     xlim[108]=0.00168, mulo[108]=  49.000, muup[108]=  56.000;
     xlim[109]=0.00193, mulo[109]=  56.000, muup[109]=  64.000;
     xlim[110]=0.00221, mulo[110]=  64.000, muup[110]=  74.000;
     xlim[111]=0.00251, mulo[111]=  74.000, muup[111]=  84.000;
     xlim[112]=0.00294, mulo[112]=  84.000, muup[112]=  97.000;
     xlim[113]=0.00340, mulo[113]=  97.000, muup[113]= 114.000;
     xlim[114]=0.00403, mulo[114]= 114.000, muup[114]= 133.000;
     xlim[115]=0.00475, mulo[115]= 133.000, muup[115]= 153.000;
     xlim[116]=0.00549, mulo[116]= 153.000, muup[116]= 174.000;
     xlim[117]=0.00647, mulo[117]= 174.000, muup[117]= 196.000;
     xlim[118]=0.00724, mulo[118]= 196.000, muup[118]= 220.000;
     xlim[119]=0.00842, mulo[119]= 220.000, muup[119]= 245.000;
     xlim[120]=0.00942, mulo[120]= 245.000, muup[120]= 272.000;
     xlim[121]=0.01058, mulo[121]= 272.000, muup[121]= 300.000;
     xlim[122]=0.01196, mulo[122]= 300.000, muup[122]= 330.000;
     xlim[123]=0.01343, mulo[123]= 330.000, muup[123]= 362.000;
     xlim[124]=0.01516, mulo[124]= 362.000, muup[124]= 395.000;
     xlim[125]=0.01695, mulo[125]= 395.000, muup[125]= 430.000;
     xlim[126]=0.01912, mulo[126]= 430.000, muup[126]= 468.000;
     xlim[127]=0.02152, mulo[127]= 468.000, muup[127]= 507.000;
     xlim[128]=0.02406, mulo[128]= 507.000, muup[128]= 548.000;
     xlim[129]=0.02721, mulo[129]= 548.000, muup[129]= 592.000;
     xlim[130]=0.03050, mulo[130]= 592.000, muup[130]= 638.000;
     xlim[131]=0.03443, mulo[131]= 638.000, muup[131]= 686.000;
     xlim[132]=0.03915, mulo[132]= 686.000, muup[132]= 737.000;
     xlim[133]=0.04476, mulo[133]= 737.000, muup[133]= 790.000;
     xlim[134]=0.05104, mulo[134]= 790.000, muup[134]= 846.000;
     xlim[135]=0.05881, mulo[135]= 846.000, muup[135]= 905.000;
     xlim[136]=0.06719, mulo[136]= 905.000, muup[136]= 967.000;
     xlim[137]=0.07650, mulo[137]= 967.001, muup[137]=1032.000;
     xlim[138]=0.08719, mulo[138]=1032.000, muup[138]=1101.000;
     xlim[139]=0.09929, mulo[139]=1101.000, muup[139]=1171.999;
     xlim[140]=0.11239, mulo[140]=1172.000, muup[140]=1248.000;
     xlim[141]=0.12760, mulo[141]=1248.000, muup[141]=1326.999;
     xlim[142]=0.14438, mulo[142]=1327.000, muup[142]=1410.000;
     xlim[143]=0.16571, mulo[143]=1410.001, muup[143]=1496.998;
     xlim[144]=0.19104, mulo[144]=1497.002, muup[144]=1587.996;
     xlim[145]=0.22049, mulo[145]=1588.002, muup[145]=1683.999;
     xlim[146]=0.26095, mulo[146]=1684.002, muup[146]=1783.996;
     xlim[147]=0.00020, mulo[147]=  10.000, muup[147]=  12.000;
     xlim[148]=0.00024, mulo[148]=  12.000, muup[148]=  15.000;
     xlim[149]=0.00030, mulo[149]=  15.000, muup[149]=  18.000;
     xlim[150]=0.00037, mulo[150]=  18.000, muup[150]=  21.000;
     xlim[151]=0.00043, mulo[151]=  21.000, muup[151]=  24.000;
     xlim[152]=0.00048, mulo[152]=  24.000, muup[152]=  28.000;
     xlim[153]=0.00057, mulo[153]=  28.000, muup[153]=  32.000;
     xlim[154]=0.00066, mulo[154]=  32.000, muup[154]=  37.000;
     xlim[155]=0.00076, mulo[155]=  37.000, muup[155]=  43.000;
     xlim[156]=0.00089, mulo[156]=  43.000, muup[156]=  49.000;
     xlim[157]=0.00103, mulo[157]=  49.000, muup[157]=  56.000;
     xlim[158]=0.00118, mulo[158]=  56.000, muup[158]=  64.000;
     xlim[159]=0.00137, mulo[159]=  64.000, muup[159]=  74.000;
     xlim[160]=0.00159, mulo[160]=  74.000, muup[160]=  84.000;
     xlim[161]=0.00187, mulo[161]=  84.000, muup[161]=  97.000;
     xlim[162]=0.00217, mulo[162]=  97.000, muup[162]= 114.000;
     xlim[163]=0.00256, mulo[163]= 114.000, muup[163]= 133.000;
     xlim[164]=0.00302, mulo[164]= 133.000, muup[164]= 153.000;
     xlim[165]=0.00359, mulo[165]= 153.000, muup[165]= 174.000;
     xlim[166]=0.00421, mulo[166]= 174.000, muup[166]= 196.000;
     xlim[167]=0.00487, mulo[167]= 196.000, muup[167]= 220.000;
     xlim[168]=0.00560, mulo[168]= 220.000, muup[168]= 245.000;
     xlim[169]=0.00650, mulo[169]= 245.000, muup[169]= 272.000;
     xlim[170]=0.00769, mulo[170]= 272.000, muup[170]= 300.000;
     xlim[171]=0.00872, mulo[171]= 300.000, muup[171]= 330.000;
     xlim[172]=0.00985, mulo[172]= 330.000, muup[172]= 362.000;
     xlim[173]=0.01140, mulo[173]= 362.000, muup[173]= 395.000;
     xlim[174]=0.01322, mulo[174]= 395.000, muup[174]= 430.000;
     xlim[175]=0.01530, mulo[175]= 430.000, muup[175]= 468.000;
     xlim[176]=0.01799, mulo[176]= 468.000, muup[176]= 507.000;
     xlim[177]=0.02104, mulo[177]= 507.000, muup[177]= 548.000;
     xlim[178]=0.02460, mulo[178]= 548.000, muup[178]= 592.000;
     xlim[179]=0.02880, mulo[179]= 592.000, muup[179]= 638.000;
     xlim[180]=0.03332, mulo[180]= 638.000, muup[180]= 686.000;
     xlim[181]=0.03852, mulo[181]= 686.000, muup[181]= 736.999;
     xlim[182]=0.04444, mulo[182]= 737.000, muup[182]= 790.000;
     xlim[183]=0.05117, mulo[183]= 790.000, muup[183]= 846.000;
     xlim[184]=0.05910, mulo[184]= 846.000, muup[184]= 905.000;
     xlim[185]=0.06934, mulo[185]= 905.000, muup[185]= 967.000;
     xlim[186]=0.08221, mulo[186]= 967.000, muup[186]=1031.997;
     xlim[187]=0.09794, mulo[187]=1032.000, muup[187]=1100.998;
     xlim[188]=0.12128, mulo[188]=1101.000, muup[188]=1171.995;
     xlim[189]=0.15284, mulo[189]=1172.003, muup[189]=1247.990;
     xlim[190]=0.20629, mulo[190]=1248.002, muup[190]=1326.996;
     xlim[191]=0.28900, mulo[191]=1327.003, muup[191]=1409.982;
     xlim[192]=0.47154, mulo[192]=1410.004, muup[192]=1472.453;
     xlim[193]=0.00013, mulo[193]=  10.000, muup[193]=  12.000;
     xlim[194]=0.00015, mulo[194]=  12.000, muup[194]=  15.000;
     xlim[195]=0.00019, mulo[195]=  15.000, muup[195]=  18.000;
     xlim[196]=0.00023, mulo[196]=  18.000, muup[196]=  21.000;
     xlim[197]=0.00027, mulo[197]=  21.000, muup[197]=  24.000;
     xlim[198]=0.00031, mulo[198]=  24.000, muup[198]=  28.000;
     xlim[199]=0.00036, mulo[199]=  28.000, muup[199]=  32.000;
     xlim[200]=0.00041, mulo[200]=  32.000, muup[200]=  37.000;
     xlim[201]=0.00048, mulo[201]=  37.000, muup[201]=  43.000;
     xlim[202]=0.00056, mulo[202]=  43.000, muup[202]=  49.000;
     xlim[203]=0.00065, mulo[203]=  49.000, muup[203]=  56.000;
     xlim[204]=0.00073, mulo[204]=  56.000, muup[204]=  64.000;
     xlim[205]=0.00090, mulo[205]=  64.000, muup[205]=  74.000;
     xlim[206]=0.00102, mulo[206]=  74.000, muup[206]=  84.000;
     xlim[207]=0.00117, mulo[207]=  84.000, muup[207]=  97.000;
     xlim[208]=0.00141, mulo[208]=  97.000, muup[208]= 114.000;
     xlim[209]=0.00173, mulo[209]= 114.000, muup[209]= 133.000;
     xlim[210]=0.00213, mulo[210]= 133.000, muup[210]= 153.000;
     xlim[211]=0.00254, mulo[211]= 153.000, muup[211]= 174.000;
     xlim[212]=0.00302, mulo[212]= 174.000, muup[212]= 196.000;
     xlim[213]=0.00359, mulo[213]= 196.000, muup[213]= 220.000;
     xlim[214]=0.00427, mulo[214]= 220.000, muup[214]= 245.000;
     xlim[215]=0.00517, mulo[215]= 245.000, muup[215]= 272.000;
     xlim[216]=0.00621, mulo[216]= 272.000, muup[216]= 300.000;
     xlim[217]=0.00740, mulo[217]= 300.000, muup[217]= 330.000;
     xlim[218]=0.00902, mulo[218]= 330.000, muup[218]= 362.000;
     xlim[219]=0.01074, mulo[219]= 362.000, muup[219]= 395.000;
     xlim[220]=0.01282, mulo[220]= 395.000, muup[220]= 429.999;
     xlim[221]=0.01514, mulo[221]= 430.000, muup[221]= 468.000;
     xlim[222]=0.01793, mulo[222]= 468.000, muup[222]= 507.000;
     xlim[223]=0.02119, mulo[223]= 507.000, muup[223]= 548.000;
     xlim[224]=0.02536, mulo[224]= 548.001, muup[224]= 592.000;
     xlim[225]=0.03097, mulo[225]= 592.001, muup[225]= 638.000;
     xlim[226]=0.03800, mulo[226]= 638.000, muup[226]= 685.998;
     xlim[227]=0.04880, mulo[227]= 686.000, muup[227]= 736.998;
     xlim[228]=0.06679, mulo[228]= 737.003, muup[228]= 789.999;
     xlim[229]=0.09428, mulo[229]= 790.001, muup[229]= 845.978;
     xlim[230]=0.15966, mulo[230]= 846.010, muup[230]= 904.906;
     xlim[231]=0.41065, mulo[231]= 905.002, muup[231]= 929.124;
     xlim[232]=0.00008, mulo[232]=  10.000, muup[232]=  12.000;
     xlim[233]=0.00009, mulo[233]=  12.000, muup[233]=  15.000;
     xlim[234]=0.00011, mulo[234]=  15.000, muup[234]=  18.000;
     xlim[235]=0.00014, mulo[235]=  18.000, muup[235]=  21.000;
     xlim[236]=0.00017, mulo[236]=  21.000, muup[236]=  24.000;
     xlim[237]=0.00020, mulo[237]=  24.000, muup[237]=  28.000;
     xlim[238]=0.00022, mulo[238]=  28.000, muup[238]=  32.000;
     xlim[239]=0.00025, mulo[239]=  32.000, muup[239]=  37.000;
     xlim[240]=0.00031, mulo[240]=  37.000, muup[240]=  43.000;
     xlim[241]=0.00036, mulo[241]=  43.000, muup[241]=  49.000;
     xlim[242]=0.00043, mulo[242]=  49.000, muup[242]=  56.000;
     xlim[243]=0.00048, mulo[243]=  56.000, muup[243]=  64.000;
     xlim[244]=0.00056, mulo[244]=  64.000, muup[244]=  74.000;
     xlim[245]=0.00068, mulo[245]=  74.000, muup[245]=  84.000;
     xlim[246]=0.00081, mulo[246]=  84.000, muup[246]=  97.000;
     xlim[247]=0.00097, mulo[247]=  97.000, muup[247]= 114.000;
     xlim[248]=0.00125, mulo[248]= 114.000, muup[248]= 133.000;
     xlim[249]=0.00157, mulo[249]= 133.000, muup[249]= 153.000;
     xlim[250]=0.00200, mulo[250]= 153.000, muup[250]= 174.000;
     xlim[251]=0.00253, mulo[251]= 174.000, muup[251]= 196.000;
     xlim[252]=0.00324, mulo[252]= 196.000, muup[252]= 220.000;
     xlim[253]=0.00405, mulo[253]= 220.000, muup[253]= 245.000;
     xlim[254]=0.00496, mulo[254]= 245.000, muup[254]= 272.000;
     xlim[255]=0.00612, mulo[255]= 272.000, muup[255]= 300.000;
     xlim[256]=0.00748, mulo[256]= 300.000, muup[256]= 330.000;
     xlim[257]=0.00939, mulo[257]= 330.000, muup[257]= 362.000;
     xlim[258]=0.01157, mulo[258]= 362.000, muup[258]= 395.000;
     xlim[259]=0.01494, mulo[259]= 395.000, muup[259]= 429.999;
     xlim[260]=0.02079, mulo[260]= 430.000, muup[260]= 467.998;
     xlim[261]=0.03128, mulo[261]= 468.000, muup[261]= 506.995;
     xlim[262]=0.05266, mulo[262]= 507.001, muup[262]= 547.993;
     xlim[263]=0.15902, mulo[263]= 548.006, muup[263]= 567.949;
   }
   // --------- fastNLO: Warm-Up run results (end)


   //printf("* --- xlimits \n");
   for(int i=0;i<A2->NObsBin;i++){
      int nxtot = 15;
      if (i == ((A2->NObsBin)-1)) nxtot += 1; // Darf's etwas mehr sein?
      B->Nxtot1.push_back(nxtot);
      double hxlim = -sqrt(-log10(xlim[i]));   // use value from Warm-Up run
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
   //B->Nscalenode.push_back(4); // number of scale nodes for pT
   B->Nscalenode.push_back(6); // number of scale nodes for pT

   B->ScaleFac.resize(B->NScaleDim);

   B->ScaleFac[0].push_back(1.0);    // --- fastNLO: central scale (don't change)
   if(nlo){
     B->ScaleFac[0].push_back(0.5);  // --- fastNLO user: add any number of
     B->ScaleFac[0].push_back(2.0);  //             additional scale variations
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
