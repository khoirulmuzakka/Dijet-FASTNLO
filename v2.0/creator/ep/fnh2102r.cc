//------ DON'T TOUCH THIS PART! ------
#include <bits/hhc-phasespace.h>
#include <bits/photo-phasespace.h>
#include <bits/hhc-process.h>
#include <bits/hhc-jetfunc.h>

//----- used namespaces -----
using namespace nlo;
using namespace std;


//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
void psinput(phasespace_photo *, double&, double&);
user_base_hhc * userfunc();

//----- array of the symbols symbols -----
extern "C"{
  
  struct { 
	const char *name;
	void *address;
  } user_defined_functions[] = 
  {
	//   process index: hhc for gamma + p --> jets (photoproduction resolved photon)
	{"procindex", (void *) "photores"},
  
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
#include <kT_clus.h>
#include "pdf-grv-cteq6.h"
#include "pdf-hhc-dummy.h"
#include "fnloTable.h"
#include "fnloBlockBNlojet.h"


//    Defines new sample type to do felavor decomposation in the initial state
const char *mySample_label[8] = {"gg", "qg", "gq", "qr", "qq", "qqb", "qrb","total"};
typedef weight<8U, mySample_label> mySample;

//   Here we specify some base user objects with the new sample type
//     (see photo-jetfunc.h and included header files therein) 
typedef basic_user<jetfunc_hhc, void,        mySample, weight_conversion> myUser0d;

class UserPhoto : public basic_user_set<myUser0d> 
{
public:
  //   init and user function
  void initfunc(unsigned int);
  void userfunc(const event_hhc&, const amplitude_hhc&);

   virtual void end_of_event();  
   virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);
  
private:
   //   pdf
   pdf_grv_cteq6 pdf;
   pdf_hhc_dummy dummypdf;

   // algorithms
   kT_clus_long jetclus;
   
   //  private types
   typedef lorentzvector<double> _Lv;

   // the jet structore
   bounded_vector<_Lv> pj_lab; 
   bounded_vector<unsigned int> jet;
  
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
  return new UserPhoto;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
  //  number of jets
  nj = 2U;
  
  //  number of the up and down type flavours
  nu = 2U;
  nd = 3U;
} 

void psinput(phasespace_photo *ps, double& el, double& eh)
{
  //  energy of incomings
  el = 27.5;
  eh = 920.0;
  
  //   You can use your own phase generator. 
  //   Here we use the default.
  ps = 0;
} 


void UserPhoto::initfunc(unsigned int)
{
   // ---- Initialize event counters
   nevents = 0;
   // Set some defaults
   if (nwrite==0) nwrite = 5000000;
    
   start_time = std::time(0);
  
}


void UserPhoto::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
   fnloBlockA2 *A2 =  table->GetBlockA2();

   //----- ZEUS cuts -----
   double pt = 0.0;
   double pTmin = 17.0, etamin = -1.0, etamax = 2.5;
   double ymin = 0.2; double ymax = 0.85;
   double Q2max = 1.0;

   //----- photon momentum fraction -----
   pdf.photon_momentum_fraction(ymin, ymax);

   //----- Bjorken x -----
   //   double x = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);
   double s = 2.0*(p[hadron(-1)]*p[hadron(0)]);
   double x = 2.0*(p[0]*p[hadron(-1)])/s;
   double x2 = 2.0*(p[-1]*p[hadron(0)])/s;

   //----- do the cluster analysis-----
   jetclus.set_up(1,false,2);
   jetclus(p); jetclus.incl(pj_lab, jet);
   unsigned int nj = pj_lab.upper();

   // loop over all jets
   for(unsigned int i = 1; i <= nj; i++){
      // jet valid?
      if((pt = pj_lab[i].perp()) > pTmin && -pj_lab[i].prapidity() < etamax && -pj_lab[i].prapidity() > etamin) { // pz in wrong direction 
         // find corresponding bin number (E_T)
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
                  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventResolved(ptbin,x,x2,pt,ymin,ymax,Q2max,amp,pdf);
               }else{
                  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventResolved(ptbin,x,x2,pt,ymin,ymax,Q2max,amp,dummypdf);
               }
            }
         }        
      }
   }
}

void UserPhoto::end_of_event(){
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

void UserPhoto::phys_output(const std::basic_string<char>& __file_name, 
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

void UserPhoto::inittable(){

   const bool doReference = true;

   //Set up fastNLO
   table = new fnloTable(tablefilename);

   table->GetBlockA1()->SetScenName("fnh2102r");
   table->GetBlockA1()->SetNcontrib(1);
   table->GetBlockA1()->SetNmult(0);
   table->GetBlockA1()->SetNdata(0);
   table->GetBlockA2()->SetIpublunits(12);

   fnloBlockA2 *A2 =  table->GetBlockA2();
   A2->ScDescript.push_back("Inclusive jet cross sections in photoproduction");
   A2->ScDescript.push_back("ZEUS Collaboration");
   A2->ScDescript.push_back("DESY 02-228");
   A2->ScDescript.push_back("?");
   A2->NScDescript = A2->ScDescript.size();
   A2->Ecms = sqrt(4.*920.*27.5);
   A2->ILOord = 2;
   A2->NDim = 1;
   A2->DimLabel.push_back("E_T");
   A2->IDiffBin.push_back(2);

   vector <double> bound;
   bound.resize(1);
  
   bound[0] = 17.; A2->LoBin.push_back(bound);
   bound[0] = 21.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 25.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 29.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 35.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 41.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 47.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 55.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 71.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 95.; A2->UpBin.push_back(bound);
   A2->NObsBin = 9;

   A2->INormFlag = 0;

   // main table
   fnloBlockB *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
   table->CreateBlockB(0,B);
   B->IXsectUnits = 12;
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
   B->CtrbDescript.push_back("resolved");
   B->NContrDescr = B->CtrbDescript.size();

   B->NPDF = 2;
   B->NPDFPDG.push_back(2212);
   B->NPDFPDG.push_back(20);
   B->NPDFDim = 2;
   B->NFragFunc = 0;
   B->IPDFdef1 = 4;
   B->IPDFdef2 = 1;
   B->IPDFdef3 = 1;

   double xlimits[9] = {0.01,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08};
   B->XNode1.resize(A2->NObsBin);
   for(int i=0;i<A2->NObsBin;i++){
      int nxtot = 14+i;
      B->Nxtot1.push_back(nxtot);
      double xlim = xlimits[i];
      double hxlim = log10(xlim);
      B->Hxlim1.push_back(hxlim);
      for(int j=0;j<nxtot;j++){
         double hx = hxlim*( 1.- ((double)j)/(double)nxtot);
         B->XNode1[i].push_back(pow(10,hx)); 
      }
   }
   
   double xlimits2[9] = {0.02,0.02,0.03,0.04,0.05,0.07,0.10,0.15,0.15};
   B->XNode2.resize(A2->NObsBin);
   for(int i=0;i<A2->NObsBin;i++){
      int nxlow = 15;
      int nxhigh = 20;
      int nxtot = nxlow + nxhigh;
      B->Nxtot2.push_back(nxtot);
      double xlim = xlimits2[i]; 
      double hxlim = log10(xlim);
      double xmax = 0.9;
      double hxmax = log10(xmax);
      B->Hxlim2.push_back(hxlim);
      for(int j=0;j<nxlow;j++){
         double hx = hxlim +  ((double)j)/((double)(nxlow)) * (hxmax-hxlim);
         B->XNode2[i].push_back(pow(10,hx)); 
      }
      for(int j=0;j<nxhigh;j++){
         B->XNode2[i].push_back(xmax +  ((double)j)/((double)(nxhigh))*(1.0-xmax)); 
      }
   }
   

   B->NScales = 2;  // two scales: mur and muf
   B->NScaleDim = 1; // only one variable used in scales: ET
   B->Iscale.push_back(0);  // mur=mur(ET), ET = index 0 
   B->Iscale.push_back(0);  // muf=muf(ET), ET = index 0 
   B->ScaleDescript.resize(B->NScaleDim);
   B->ScaleDescript[0].push_back(" E_T of jet");
   B->NscaleDescript.push_back(B->ScaleDescript[0].size());
   B->Nscalenode.push_back(4); // number of scale nodes

   B->ScaleFac.resize(B->NScaleDim);
   for(int i=0;i<B->NScaleDim;i++){
      B->ScaleFac[i].push_back(1.0);
      B->ScaleFac[i].push_back(0.5);
      B->ScaleFac[i].push_back(2.0);
      B->Nscalevar.push_back(B->ScaleFac[i].size());
   }

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
      refB->Nxtot2.clear();
      refB->Hxlim2.clear();
      for(int i=0;i<A2->NObsBin;i++){
         refB->Nxtot1.push_back(1);
         refB->Hxlim1.push_back(0.);
         refB->XNode1[i].clear(); 
         refB->XNode1[i].push_back(0.); 
         refB->Nxtot2.push_back(1);
         refB->Hxlim2.push_back(0.);
         refB->XNode2[i].clear(); 
         refB->XNode2[i].push_back(0.); 
      }
      table->GetBlockA1()->SetNcontrib(2);
   }

}
 
void UserPhoto::writetable(){
   table->OpenFileRewrite();
   table->WriteBlockA1();
   table->WriteBlockA2();
   for(int i=0;i< table->GetBlockA1()->GetNcontrib();i++){
      table->WriteBlockBDividebyN(i);
   }
   table->CloseFileWrite();
}




