//
// fastNLO author code for fnt1001 for midpoint ago (no Rsep)
//  Run I CDF incl jets:   CDF publication: hep-ph/0102074
//       T. Affolder et al. (CDF Collaboration),
//                    Phys. Rev. D64, 032001 (2001).
//
//  -> the published ansatz was used to derive bins corresponding 
//     to the published CDF bin-centers
// 
// last modification
// 2005/12/07 MW - new xlimit calculation (y-integration)
// 2005/12/.. TK - add "reference" setting
// 2006/01/13 MW - make "reference" option switchable "iref=0/1"
//                 scale variations are now ordered ; and mur=muf
//                 implement bicubic eigenfunctions
// 2006/01/17 MW - divide by binwidth 
//                 store in fb (as in publication) -> variable "unitfactor"
//                 extend tableformat to version 1c
//             --> used for production
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
#include "cone-et-07.h"
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
   cone_et_07 jetclus;   // jet algorithm
   
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
   s = 3240000;       // TeV Run I     1800GeV
   //s = 3841600;     // TeV Run II    1960GeV
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

   double pTmin = 40.5, etamin = 0.1, etamax = 0.7;

   double x1 = p[-1].Z()/p[hadron(-1)].Z();
   double x2 = p[0].Z()/p[hadron(0)].Z();

   //----- do the jet analysis -----
   pj = jetclus(p);
   unsigned int nj = pj.upper(); 

   // loop over all jets
   for(unsigned int i = 1; i <= nj; i++){
      double pt = pj[i].perp(); 
      double prap = fabs(pj[i].prapidity());
      // jet valid?
      if(pt > pTmin &&  prap < etamax && prap > etamin) { 
         // find corresponding bin number (p_T)
         int ptbin = -1;
         for(int j = 0; j < A2->GetNObsBin(); j++) {
            if (pt >= A2->LoBin[j][0]  && pt <  A2->UpBin[j][0]) {
               ptbin=j;
               break;
            }
         }
         //---------- fill fastNLO arrays
         if ( ptbin>=0 ) {
            for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
               if(table->GetBlockB(k)->GetIRef()>0){
                  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(ptbin,x1,x2,pt,amp,pdf,0.5); // 0.5 = factor for bin in |eta| -> eta
               }else{
                  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(ptbin,x1,x2,pt,amp,dummypdf,0.5);
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
   //-------- store table
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

   double s = 3240000; 
   
   const bool doReference = true;

   //Set up fastNLO
   table = new fnloTable(tablefilename);

   table->GetBlockA1()->SetScenName("fnt1001midp");
   table->GetBlockA1()->SetNcontrib(1);
   table->GetBlockA1()->SetNmult(0);
   table->GetBlockA1()->SetNdata(0);
   table->GetBlockA2()->SetIpublunits(9);

   fnloBlockA2 *A2 =  table->GetBlockA2();
   A2->ScDescript.push_back("d2sigma-jet_dET_deta_(nb_GeV)");
   A2->ScDescript.push_back("hep-ph-0102074");
   A2->ScDescript.push_back("CDF_Collaboration");
   A2->ScDescript.push_back("");
   A2->NScDescript = A2->ScDescript.size();
   A2->Ecms = sqrt(s);
   A2->ILOord = 2;
   A2->NDim = 2;
   A2->DimLabel.push_back("E_T");
   A2->IDiffBin.push_back(2);
   A2->DimLabel.push_back("y");
   A2->IDiffBin.push_back(2);

   vector <double> bound;
   bound.resize(2);

   const int nptbins = 33;
   double ptbins[nptbins+1] = {40.5,46.5,52.44,58.26,64.01,69.63,75.19,80.82,86.37,91.8,97.37,102.78,108.38,113.55,119.2,124.31,130.03,135.07,140.86,150.96,162.35,172.42,183.84,193.9,205.54,215.14,237.13,258.36,280.59,301.56,323.9,344.32,383.76,452.71};
   const int nrapbins = 1;
   double rapbins[nrapbins+1] = {0.1,0.7};

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
   B->CodeDescript.push_back("NLOJET++ 4.0.1");
   B->NCodeDescr = B->CodeDescript.size();
   B->IRef = 0;
   B->NSubproc = 7;
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
   B->IPDFdef3 = 1;

   B->XNode1.resize(A2->NObsBin);
   for(int i=0;i<A2->NObsBin;i++){
      int nxtot = 12;
      B->Nxtot1.push_back(nxtot);
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
      // ---- safety factors for ET-scheme / optimized by eta range
      double hxlim = -sqrt(-log10(xmin*0.81));
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

   B->ScaleDescript[0].push_back(" p_T of jet");
   B->NscaleDescript.push_back(B->ScaleDescript[0].size());
   B->Nscalenode.push_back(4); // number of scale nodes for pT

   B->ScaleFac.resize(B->NScaleDim);

   B->ScaleFac[0].push_back(1.0);
   //   B->ScaleFac[0].push_back(0.5);
   //   B->ScaleFac[0].push_back(2.0);
   B->Nscalevar.push_back(B->ScaleFac[0].size());

   B->ScaleNode.resize(A2->NObsBin);
   for(int i=0;i<A2->NObsBin;i++){
      B->ScaleNode[i].resize(B->NScaleDim);
      for(int j=0;j<B->NScaleDim;j++){
         B->ScaleNode[i][j].resize(B->Nscalevar[j]);
         for(int k=0;k<B->Nscalevar[j];k++){
            B->ScaleNode[i][j][k].resize(B->Nscalenode[j]);
            if(B->Nscalenode[j]==1){
               B->ScaleNode[i][j][k][0] = B->ScaleFac[0][k]*(A2->UpBin[i][B->Iscale[j]] + A2->LoBin[i][B->Iscale[j]])/2.;
            }else{
               const double mu0scale = .25; // In GeV
               double llscalelo = log(log((B->ScaleFac[0][k]*A2->LoBin[i][B->Iscale[j]])/mu0scale));
               double llscalehi = log(log((B->ScaleFac[0][k]*A2->UpBin[i][B->Iscale[j]])/mu0scale));
               for(int l=0;l<B->Nscalenode[j];l++){
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
      refB->NSubproc = 7;
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
