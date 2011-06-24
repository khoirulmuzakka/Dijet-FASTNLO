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
   cone_et_07 jetclus;   // fastNLO user: define jet algorithm (consistent with .h file above)
   
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

   double x1 = p[-1].Z()/p[hadron(-1)].Z();
   double x2 = p[0].Z()/p[hadron(0)].Z();

   //----- run the jet algorithm
   pj = jetclus(p);
   unsigned int nj = pj.upper(); 

   // ===== start of the fastNLO user part
   // =====
   double pTmin = 40.5, etamin = 0.1, etamax = 0.7;

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
         // --- fill fastNLO arrays
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

   double s = 3240000; // Tevatron RunI
   
   const bool doReference = true;
   //const bool doReference = false;

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
   //A2->ScDescript.push_back("");
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
   double ptbins[nptbins+1] = {
     40.5, 46.5, 52.44, 58.26, 64.01, 69.63, 75.19, 80.82, 86.37, 91.8,
     97.37, 102.78, 108.38, 113.55, 119.2, 124.31, 130.03, 135.07, 140.86, 150.96,
     162.35, 172.42, 183.84, 193.9, 205.54, 215.14, 237.13, 258.36, 280.59, 301.56,
     323.9, 344.32, 383.76, 452.71};
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
   //B->IWarmUp = 1;   // do the Warm-Up run (initialize x,mu limits) 
   B->IWarmUp = 0;   // production
   B->IWarmUpPrint = 1000000;
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
   // 307000000 contributions (!= events) in warm-up run 
   xlim[0]=0.011738, mulo[0]=40.500001, muup[0]=46.500000;
   xlim[1]=0.013584, mulo[1]=46.500001, muup[1]=52.439998;
   xlim[2]=0.015400, mulo[2]=52.440002, muup[2]=58.260000;
   xlim[3]=0.017213, mulo[3]=58.260004, muup[3]=64.009999;
   xlim[4]=0.019006, mulo[4]=64.010001, muup[4]=69.630000;
   xlim[5]=0.020769, mulo[5]=69.630000, muup[5]=75.189999;
   xlim[6]=0.022643, mulo[6]=75.190002, muup[6]=80.819998;
   xlim[7]=0.024124, mulo[7]=80.820000, muup[7]=86.369999;
   xlim[8]=0.026265, mulo[8]=86.370002, muup[8]=91.799999;
   xlim[9]=0.027668, mulo[9]=91.800001, muup[9]=97.369998;
   xlim[10]=0.030170, mulo[10]=97.370006, muup[10]=102.779999;
   xlim[11]=0.032094, mulo[11]=102.780001, muup[11]=108.379998;
   xlim[12]=0.034015, mulo[12]=108.380006, muup[12]=113.549999;
   xlim[13]=0.035803, mulo[13]=113.550001, muup[13]=119.199998;
   xlim[14]=0.037559, mulo[14]=119.200001, muup[14]=124.309999;
   xlim[15]=0.039124, mulo[15]=124.310001, muup[15]=130.029992;
   xlim[16]=0.041257, mulo[16]=130.030000, muup[16]=135.069999;
   xlim[17]=0.043825, mulo[17]=135.070007, muup[17]=140.859998;
   xlim[18]=0.045953, mulo[18]=140.860003, muup[18]=150.959998;
   xlim[19]=0.048473, mulo[19]=150.960002, muup[19]=162.349999;
   xlim[20]=0.053317, mulo[20]=162.350002, muup[20]=172.419995;
   xlim[21]=0.057052, mulo[21]=172.420003, muup[21]=183.839999;
   xlim[22]=0.061254, mulo[22]=183.840003, muup[22]=193.899998;
   xlim[23]=0.065609, mulo[23]=193.900003, muup[23]=205.539999;
   xlim[24]=0.069768, mulo[24]=205.540007, muup[24]=215.139997;
   xlim[25]=0.074718, mulo[25]=215.140010, muup[25]=237.129994;
   xlim[26]=0.083468, mulo[26]=237.130001, muup[26]=258.359997;
   xlim[27]=0.093144, mulo[27]=258.360003, muup[27]=280.589991;
   xlim[28]=0.104822, mulo[28]=280.590007, muup[28]=301.559996;
   xlim[29]=0.114991, mulo[29]=301.560003, muup[29]=323.899986;
   xlim[30]=0.129211, mulo[30]=323.900038, muup[30]=344.319995;
   xlim[31]=0.141246, mulo[31]=344.320013, muup[31]=383.760000;
   xlim[32]=0.161341, mulo[32]=383.760009, muup[32]=452.709998;
   // --------- end: insert here (copy&paste) results from warm-up run

   //printf("* --- xlimits \n");
   for(int i=0;i<A2->NObsBin;i++){
      int nxtot = 12;
      if (i == ((A2->NObsBin)-1)) nxtot = 13; // Darf's etwas mehr sein?
      B->Nxtot1.push_back(nxtot);
      double hxlim = -sqrt(-log10(xlim[i]));
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
      //  refB->NSubproc = 7;
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
