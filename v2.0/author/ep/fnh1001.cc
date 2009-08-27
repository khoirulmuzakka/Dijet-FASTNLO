//------ DON'T TOUCH THIS PART! ------
#include <bits/dis-phasespace.h>
#include <bits/dis-process.h>
#include <bits/dis-jetfunc.h>

//----- used namespaces -----
using namespace nlo;
using namespace std;


//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
void psinput(phasespace_dis *, double&, double&, double&, double&, double&, double&, double&, double&);
user_base_dis * userfunc();

//----- array of the symbols symbols -----
extern "C"{
  
  struct { 
	const char *name;
	void *address;
  } user_defined_functions[] = 
  {
	//   process index: hhc for e + p --> jets (DIS)
	{"procindex", (void *) "dis"},
  
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
#include "pdf-cteq6.h"
#include "pdf-dis-dummy.h"
#include "fnloTable.h"
#include "fnloBlockBNlojet.h"


class UserDIS : public basic_user_set<user1d_dis, user1h_dis> 
{
public:
  //   init and user function
   void initfunc(unsigned int);
   void userfunc(const event_dis&, const amplitude_dis&);

   virtual void end_of_event();  
   virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);
  
private:
   //   pdf
   pdf_cteq6dis pdf;
   pdf_dis_dummy dummypdf;

   // algorithms
   kT_clus_long jetclus;
   
   // the jet structore in breit & lab. frame
   bounded_vector<lorentzvector<double> > pjb, pjl; 
   bounded_vector<unsigned int> jet;
   
   //  event in breit frame
   event_dis pbreit;
  
  // boost the jet momenta back to the laboratory frame
  void boost_back_to_lab(const event_dis&);
  
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
  

user_base_dis * userfunc() {
  return new UserDIS;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
  //  number of jets
  nj = 2U;
  
  //  number of the up and down type flavours
  nu = 2U;
  nd = 3U;
} 

void psinput(phasespace_dis *ps, double& el, double& eh, double& q2min, double& q2max, 
			 double& xmin, double& xmax, double& ymin, double& ymax)
{
  // energy of the incomong lepton and 
  // hadron in the laboratory frame
  el = 27.5;    // GeV
  eh = 820.0;     // GeV
  
  //  Q^2 cuts
  q2min = 150.0;    // GeV^2
  q2max = 5000.0;  // GeV^2 

  //   xB cuts
  xmin = 0.0;
  xmax = 1.0;
  
  //   y cuts
  ymin = 0.2;
  ymax = 0.6;

  //   You can use your own phase generator. 
  //   Here we use the default.
  ps = 0;
} 


void UserDIS::initfunc(unsigned int)
{
   // ---- Initialize event counters
   nevents = 0;
   // Set some defaults
   // KR: Change for testing
   //   if (nwrite==0) nwrite = 5000000;
   if (nwrite==0) nwrite = 5000;

   // KR: Avoid seg faults, initialize table here and not only in
   //     UserDIS::phys_output where it seems to be too late ?
   std::cout << "   before inittable   " << endl;
   inittable();
   std::cout << "   after inittable   " << endl;
    
   start_time = std::time(0);
   
}

extern"C" double xalfaem_(double *);

double xalpha_em(double mq2) {
  return xalfaem_(&mq2);
}



void UserDIS::userfunc(const event_dis& p, const amplitude_dis& amp)
{
   fnloBlockA2 *A2 =  table->GetBlockA2();

   //----- H1 cuts -----
   double pt = 0.0;
   double pTmin = 7.0, etamin = -1.0, etamax = 2.5;

   double alem = 0; // for running alpha EM

   //----- Bjorken x -----
   double x = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);
   double Q2 = -((p[-1] - p[-2]).mag2());
   alem = xalpha_em(Q2);

   //----- copy the momenta and boost to the breit frame ----
   pbreit = p; lab_to_breit(pbreit);

   //----- do the cluster analysis-----
   jetclus.set_up(1,false,2);
   jetclus(pbreit); jetclus.incl(pjb, jet);
   unsigned int nj = pjb.upper(); 
   
   //----- jet structure in laboratory frame -----
   pjl = pjb; boost_back_to_lab(p);

   // loop over all jets
   for(unsigned int i = 1; i <= nj; i++){
      // jet valid?
      if((pt = pjb[i].perp()) > pTmin && pjl[i].prapidity() < etamax && pjl[i].prapidity() > etamin) { 
         // find corresponding bin number (E_T and Q^2)
         int bin = -1;
         for(int j = 0; j < A2->GetNObsBin(); j++) {
            if (pt >= A2->LoBin[j][0]  && pt <  A2->UpBin[j][0] &&
                Q2 >= A2->LoBin[j][1]  && Q2 <  A2->UpBin[j][1]) {
               bin=j;
               break;
            }
         }

         //---------- fill fastNLO arrays
         if ( bin>=0 ) {
            for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
               if(table->GetBlockB(k)->GetIRef()>0){
                  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventDIS(bin,x,pt,amp,pdf,alem * alem);
               }else{
                  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventDIS(bin,x,pt,amp,dummypdf,alem * alem);
               }
            }
         }        
      }
   }

}

void UserDIS::end_of_event(){
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

void UserDIS::phys_output(const std::basic_string<char>& __file_name, 
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

void UserDIS::inittable(){

   const bool doReference = true;
   // KR: Workaround since I couldn't figure out how the filename setting
   //     could work with inittable() in initfunc
   string tablefilename = "fastNLO.tab";

   //Set up fastNLO
   table = new fnloTable(tablefilename);

   table->GetBlockA1()->SetScenName("fnh1001");
   table->GetBlockA1()->SetNcontrib(1);
   table->GetBlockA1()->SetNmult(0);
   table->GetBlockA1()->SetNdata(0);
   table->GetBlockA2()->SetIpublunits(12);

   fnloBlockA2 *A2 =  table->GetBlockA2();
   A2->ScDescript.push_back("Inclusive jet cross sections in DIS");
   A2->ScDescript.push_back("H1 Collaboration");
   A2->ScDescript.push_back("hep-ex/0010054");
   A2->ScDescript.push_back("Eur.Phys.J. C19 (2001) 289-311");
   A2->NScDescript = A2->ScDescript.size();
   A2->Ecms = sqrt(4.*820.*27.5);
   A2->ILOord = 1;
   A2->NDim = 2;
   A2->DimLabel.push_back("E_T");
   A2->IDiffBin.push_back(2);
   A2->DimLabel.push_back("Q^2");
   A2->IDiffBin.push_back(2);

   vector <double> bound;
   bound.resize(2);

   const int netbins = 4;
   const int nq2bins = 4;
   double etbins[netbins+1] = {7.,11.,18.,30.,50.};
   double q2bins[nq2bins+1] = {150.,200.,300.,600.,5000.};

   for(int i=0;i<nq2bins;i++){
      for(int j=0;j<netbins;j++){
         bound[0] = etbins[j];
         bound[1] = q2bins[i];
         A2->LoBin.push_back(bound);
         bound[0] = etbins[j+1];
         bound[1] = q2bins[i+1];
         A2->UpBin.push_back(bound);
      }
   }
   A2->NObsBin = netbins*nq2bins;

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
   B->NSubproc = 3;
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

   B->NPDF = 1;
   B->NPDFPDG.push_back(2212);
   B->NPDFDim = 0;
   B->NFragFunc = 0;
   B->IPDFdef1 = 2;
   B->IPDFdef2 = 1;
   B->IPDFdef3 = 1;

   B->XNode1.resize(A2->NObsBin);
   for(int i=0;i<A2->NObsBin;i++){
      int nxtot = 20;
      B->Nxtot1.push_back(nxtot);
      double s = pow(A2->Ecms,2);
      double xlim = 0.8*(A2->LoBin[i][1]/s/0.6+pow(A2->LoBin[i][0] ,2)/s);
      //      printf("bin#=%d Q2=%f et=%f s=%f xlim=%f\n",i,A2->LoBin[i][1],A2->LoBin[i][0],s,xlim);
      double hxlim = log10(xlim);
      B->Hxlim1.push_back(hxlim);
      for(int j=0;j<nxtot;j++){
         double hx = hxlim*( 1.- ((double)j)/(double)nxtot);
         B->XNode1[i].push_back(pow(10,hx)); 
      }
   }

   B->NScales = 2;  // two scales: mur and muf
   B->NScaleDim = 1; // only one variable used in scales: ET
   B->Iscale.push_back(0);  // mur=mur(ET), ET = index 0 
   B->Iscale.push_back(0);  // muf=muf(ET), ET = index 0 
   B->ScaleDescript.resize(B->NScaleDim);
   B->ScaleDescript[0].push_back(" E_T of jet");
   B->NscaleDescript.push_back(B->ScaleDescript[0].size());
   B->Nscalenode.push_back(10); // number of scale nodes

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
      refB->NSubproc = 3;
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

   if(nlo){
      B->NSubproc = 3;
   }else{
      B->NSubproc = 2;
   }



}

void UserDIS::writetable(){
   table->OpenFileRewrite();
   table->WriteBlockA1();
   table->WriteBlockA2();
   for(int i=0;i< table->GetBlockA1()->GetNcontrib();i++){
      table->WriteBlockBDividebyN(i);
   }
   table->CloseFileWrite();
}


void UserDIS::boost_back_to_lab(const event_dis& p)
{ 
  double x = p[-1].T()/p[hadron(0)].T();
  double bz = (1.0 - x)/(1.0 + x);
  threevector<double> bVec = -((pbreit[-1] + pbreit[hadron(0)]).boostVector());
  lorentzvector<double> p0(pbreit[hadron(0)]);
  
  p0.boost(bVec);
  double phi = p0.phi(), theta = p0.theta();
  unsigned int njet = pjl.upper();
  
  for(unsigned int i = 1; i <= njet; i++) {
    pjl[i].boost(bVec);
    pjl[i].rotateZ(-phi);
    pjl[i].rotateY(-theta); 
    pjl[i].boost(0.0, 0.0, bz);
  }
}
