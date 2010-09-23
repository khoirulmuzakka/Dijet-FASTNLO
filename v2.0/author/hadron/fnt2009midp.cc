//
// fastNLO v2 author code for fnt2009 for midpoint algo (no Rsep)
//  Run IIa D0 incl jets
// 
// last modification
// 2010/09/22 MW implement D0 phase space
// 2009/01/15 TK   make code V2.0 compatible

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
#include "cone-e-07.h"
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
   //   pdf
   pdf_cteq6 pdf;
   pdf_hhc_dummy dummypdf;

   // algorithms
   cone_e_07 jetclus;   // jet algorithm
   
   bounded_vector<lorentzvector<double> > pj;    // the jet structure 
   
   //fastNLO starts here
   
   fnloTable *table;
   
   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table

   string tablefilename; // The table file to write to
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
   //  number of jets 
   //nj = 1U;
   nj = 2U;
   //nj = 3U;

   //  number of the up and down type flavours
   nu = 2U;
   nd = 3U;
} 


void psinput(phasespace_hhc *ps, double& s)
{
   //  total c.m. energy square
   //s = 40000;       // RHIC           200GeV   
   //s = 3240000;       // TeV Run I     1800GeV
   s = 3841600;     // TeV Run II    1960GeV
   //s = 196000000;   // LHC          14000GeV
  
   //   You can use your own phase generator. 
   //   Here we use the default.
   ps = 0;
} 


void UserHHC::initfunc(unsigned int)
{
   // ---- Initialize event counters
   nevents = 0;
   // Set some defaults
   if (nwrite==0) nwrite = 5000000;
    
   start_time = std::time(0);

}


void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
   fnloBlockA2 *A2 =  table->GetBlockA2();

   double pTmin = 50., ymin = 0., ymax = 0.4;  // later 2.4

   double x1 = p[-1].Z()/p[hadron(-1)].Z();
   double x2 = p[0].Z()/p[hadron(0)].Z();

   //----- do the jet analysis -----
   pj = jetclus(p);
   unsigned int nj = pj.upper(); 

   // loop over all jets
   for(unsigned int i = 1; i <= nj; i++){
      double pt = pj[i].perp(); 
      double rap = fabs(pj[i].rapidity());
      // jet valid?
      if(pt > pTmin &&  rap < ymax) {     // && rap > ymin) { 
	// MW: how to search for bin? stupid: linear in all bins? 
         // find corresponding bin number (y,pT)
         int obsbin = -1;
         for(int j = 0; j < A2->GetNObsBin(); j++) {
            if (pt >= A2->LoBin[j][0]  && pt <  A2->UpBin[j][0] && 
		rap >= A2->LoBin[j][1] && rap < A2->UpBin[j][1]) {
               obsbin=j;
               break;
            }
         }
         // --- fill fastNLO arrays
         if ( obsbin>=0 ) {
            for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
               if(table->GetBlockB(k)->GetIRef()>0){
                  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,pt,amp,pdf,0.5); // 0.5 = factor for bin in |y| -> y
               }else{
                  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(obsbin,x1,x2,pt,amp,dummypdf,0.5);
               }
            }
         }        
      }
   }

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
   
   //Determine whether we are running LO or NLO
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

   //double s = 3240000; // Tevatron RunI
   double s = 3841600; // Tevatron RunII
   
   //const bool doReference = true;
   const bool doReference = false;

   //Set up fastNLO
   table = new fnloTable(tablefilename);

   table->GetBlockA1()->SetScenName("fnt2009");
   table->GetBlockA1()->SetNcontrib(1);
   table->GetBlockA1()->SetNmult(0);
   table->GetBlockA1()->SetNdata(0);
   table->GetBlockA2()->SetIpublunits(12);

   fnloBlockA2 *A2 =  table->GetBlockA2();
   A2->ScDescript.push_back("d2sigma-jet_dpT_dy_(pb_GeV)");
   A2->ScDescript.push_back("xxx");
   A2->ScDescript.push_back("D0_Collaboration");
   //A2->ScDescript.push_back("");
   A2->NScDescript = A2->ScDescript.size();
   A2->Ecms = sqrt(s);
   A2->ILOord = 2;
   A2->NDim = 2;
   A2->DimLabel.push_back("pT");
   A2->IDiffBin.push_back(2);
   A2->DimLabel.push_back("y");
   A2->IDiffBin.push_back(2);

   vector <double> bound;
   bound.resize(2);

   const int nptbins = 23;
   double ptbins[nptbins+1] = {
     50., 60., 70., 80., 90., 100., 110., 120., 130., 145., 
     160., 180., 200., 220., 240., 265., 295., 325., 360., 400.,
     445., 490., 540., 665.};
   const int nrapbins = 1;
   double rapbins[nrapbins+1] = {0.,0.4};

   for(int i=0;i<nrapbins;i++){
      for(int j=0;j<nptbins;j++){
         bound[0] = ptbins[j];
         bound[1] = rapbins[i];
         A2->LoBin.push_back(bound);
         bound[0] = ptbins[j+1];
         bound[1] = rapbins[i+1];
         A2->UpBin.push_back(bound);
      }
   }
   A2->NObsBin = nptbins*nrapbins;

   A2->INormFlag = 0;

   // main table
   fnloBlockB *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
   table->CreateBlockB(0,B);
   B->IXsectUnits = 9;
   B->IDataFlag = 0;
   B->IAddMultFlag = 0;
   B->IContrFlag1 = 1;
   B->IContrFlag3 = 0;
   B->CodeDescript.push_back("NLOJET++ 4.1.3");
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
   B->NPDFPDG.push_back(2212);
   B->NPDFPDG.push_back(-2212);
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

   // ===== MW: dirty solution for WarmUp run (but it works)
   // --- initialize variables for WarmUp run
   B->IWarmUp = 1;   // do the Warm-Up run (initialize x,mu limits) 
   //B->IWarmUp = 0;   // production
   B->IWarmUpPrint = 50000000;
   B->xlo.resize(A2->NObsBin);
   B->scalelo.resize(A2->NObsBin);
   B->scalehi.resize(A2->NObsBin);

   // --- arrays for extreme x and (default) scale values
   double xlim[A2->NObsBin];
   double mulo[A2->NObsBin];
   double muup[A2->NObsBin];
   for(int i=0;i<A2->NObsBin;i++){
     xlim[i] = 1.1e-07, mulo[i] = 3.0, muup[i]=9.9e10;
   }
   
   // --------- start: insert here (copy&paste) results from warm-up run

   // --------- end: insert here (copy&paste) results from warm-up run

   //printf("* --- xlimits \n");
   for(int i=0;i<A2->NObsBin;i++){
      int nxtot = 12;
      if (i == ((A2->NObsBin)-1)) nxtot = 13; // Darf's etwas mehr sein?
      B->Nxtot1.push_back(nxtot);
 
      /*   ------------ old code before implementation of Warm-Up Run
      // - Setup the xlimit array - computed from kinematic constraints
      double pt = A2->LoBin[i][0];
      double raphigh = A2->UpBin[i][1];
      double xt = 2*pt/sqrt(s);
      double ymax = log((1.+sqrt(1.-xt*xt))/xt);  // upper kin. y-limit
      if (ymax>raphigh) ymax=raphigh;
      double ymin = A2->LoBin[i][1];
      //   find smallest x by integrating over accessible y-range
      double xmin = 1.0; 
      for (int nr = 0; nr <= 400; nr++) {
	double ytest = ymin + double(nr)*(ymax-ymin)/400.0;
	double xtest = pt*exp(-ytest)/(sqrt(s)-pt*exp(ytest));
	if (xtest<xmin) xmin = xtest;
      }  
      */

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

   B->ScaleDescript[0].push_back("pT of jet");
   B->NscaleDescript.push_back(B->ScaleDescript[0].size());
   B->Nscalenode.push_back(4); // number of scale nodes for pT

   B->ScaleFac.resize(B->NScaleDim);

   B->ScaleFac[0].push_back(1.0);
   if(nlo){
     B->ScaleFac[0].push_back(0.5);
     B->ScaleFac[0].push_back(2.0);
     B->ScaleFac[0].push_back(0.25);
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
