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
//  - the jet algorithm ("#include" statement and assignement of "jetclus")
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
   //double s =     40000.; // RHIC               200 GeV
   //double s =   3240000.; // TeV Run I         1800 GeV
   //double s =   3841600.; // TeV Run II        1960 GeV
   //double s =    810000.; // LHC Injection Run  900 GeV
   //double s =   5569600.; // LHC Initial Run   2360 GeV
   double s =  49000000.; // LHC First Run     7000 GeV
   //double s = 100000000.; // LHC Start-up Run 10000 GeV
   //double s = 196000000.; // LHC Design Run   14000 GeV

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
   //const bool doReference = true;
   const bool doReference = false;

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
   for (unsigned int i=0; i<nrapbins; i++) {
     for (unsigned int j=0; j<nptbins[i]+1; j++) { 
       pthigh[i][j] = ptb[j];
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
   B->IWarmUp = 1;       // --- fastNLO user: do the Warm-Up run
   //B->IWarmUp = 0;     //                   or do production run(s)
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
   // 500000000 contributions (!= events) in warm-up run 
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
