//
//   fastNLO v2 author code 
//   scenario: fnh1001
//   H1 incl jets (@high Q2)
//   for kT algorithm
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
//   xalfaem     (don't touch)
//   userfunc    (-> user edits)
//   end_of_event (don't touch)
//   phys_output (don't touch)
//   inittable   (-> user edits)
//   writetable  (don't touch)
//   boost_back_to_lab (don't touch)
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
// 2010/10/06 MW final fnh1001 version
// 2010/09/28 MW make user-friendly
// 2009/01/15 TK make code V2.0 compatible
//


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

#include <kT_clus.h>         // fastNLO user: .h file for jet algorithm

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

   // jet algorithm
   kT_clus_long jetclus;   // fastNLO user: define jet algorithm (consistent with .h file above)
   
   // the jet structore in breit & lab. frame
   bounded_vector<lorentzvector<double> > pjb, pjl; 
   bounded_vector<unsigned int> jet;
   
   //  event in breit frame
   event_dis pbreit;
  
  // boost the jet momenta back to the laboratory frame
  void boost_back_to_lab(const event_dis&);
  
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
  
user_base_dis * userfunc() {
  return new UserDIS;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
  // --- fastNLO user: select the number of jets for your observable
  nj = 2U;
  //nj = 3U;
  
  // --- number of the up and down type flavours
  nu = 2U;
  nd = 3U;
} 

void psinput(phasespace_dis *ps, double& el, double& eh, double& q2min, double& q2max, 
			 double& xmin, double& xmax, double& ymin, double& ymax)
{
  // --- fastNLO user: set the energies of the incoming lepton and hadron in the laboratory frame
  el = 27.5;      // in GeV
  eh = 820.0;     // in GeV
  
  // --- fastNLO user: define the Q^2 boundaries of the phase space
  q2min = 150.0;    // in GeV^2
  q2max = 5000.0;   // in GeV^2 

  // --- fastNLO user: define the x_Bjorken boundaries of the phase space
  xmin = 0.0;
  xmax = 1.0;
  
  // --- fastNLO user: define the y boundaries of the phase space
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
   if (nwrite==0) nwrite = 5000000;

   start_time = std::time(0);
   
}

extern"C" double xalfaem_(double *);

double xalpha_em(double mq2) {
  return xalfaem_(&mq2);
}



void UserDIS::userfunc(const event_dis& p, const amplitude_dis& amp)
{
   fnloBlockA2 *A2 =  table->GetBlockA2();

   double alem = 0; // for running alpha EM
   // --- x where PDF is probed ---
   double x = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);
   double Q2 = -((p[-1] - p[-2]).mag2());
   alem = xalpha_em(Q2);


   // --- fastNLO user:
   //     Here is your playground where you compute your observable 
   //     and the bin number ("obsbin") which gets passed to
   //     fastNLO's table filling code.
   //     (all pT, ET and E are in GeV)

   // --- declare and initialize phase space cut variables
   //     here: H1 cuts
   double pt = 0.0;
   double pTmin = 7.0, etamin = -1.0, etamax = 2.5;

   // --- copy the momenta and boost to the breit frame ----
   pbreit = p; lab_to_breit(pbreit);

   // --- do the cluster analysis ----
   jetclus.set_up(1,false,2);
   jetclus(pbreit); jetclus.incl(pjb, jet);
   unsigned int nj = pjb.upper(); 
   
   // --- jet structure in laboratory frame -----
   pjl = pjb; boost_back_to_lab(p);

   // --- loop over all jets
   for(unsigned int i = 1; i <= nj; i++){
      // jet valid?
      if((pt = pjb[i].perp()) > pTmin && pjl[i].prapidity() < etamax && pjl[i].prapidity() > etamin) { 

  	 // --- set the renormalization and factorization scale to jet pT
	 double mu = pt; 

         // --- find corresponding bin number (E_T and Q^2)
         int bin = -1;
         for(int j = 0; j < A2->GetNObsBin(); j++) {
            if (pt >= A2->LoBin[j][0]  && pt <  A2->UpBin[j][0] &&
                Q2 >= A2->LoBin[j][1]  && Q2 <  A2->UpBin[j][1]) {
               bin=j;
               break;
            }
         }

         //---------- fill fastNLO arrays - don't touch this piece of code!
         if ( bin>=0 ) {
	    // --- prefactor: divide by binwidth - and multiply with alpha_elm^2
	    double prefactor = 1./A2->BinSize[bin] * alem*alem; 
            for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
               if(table->GetBlockB(k)->GetIRef()>0){
 		  // --- fastNLO user: don't modify the following calls!
                  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventDIS(bin,x,mu,amp,pdf,prefactor);
               }else{
                  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventDIS(bin,x,mu,amp,dummypdf,prefactor);
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

void UserDIS::inittable(){

   // --- fastNLO user: decide whether to include a reference table (for 
   //                   precision studies, not for production jobs)
   const bool doReference = true;
   //const bool doReference = false;

   // --- set up fastNLO
   table = new fnloTable(tablefilename);

   // --- fastNLO: fill variable for table header block A1
   table->GetBlockA1()->SetScenName("fnh1001");  // - fastNLO user: set scenario name
   table->GetBlockA1()->SetNcontrib(1);
   table->GetBlockA1()->SetNmult(0);
   table->GetBlockA1()->SetNdata(0);
   table->GetBlockA2()->SetIpublunits(12);  // - fastNLO user: set cross section units
                                            //                 (negative power of ten)
                                            //                 fb:15  pb:12  nb:9  etc.
   // --- fastNLO: fill variables for table header block A2
   fnloBlockA2 *A2 =  table->GetBlockA2();

   // --- fastNLO user: up to 30 strings to describe the scenario
   A2->ScDescript.push_back("Inclusive jet cross sections in DIS");
   A2->ScDescript.push_back("H1 Collaboration");
   A2->ScDescript.push_back("hep-ex/0010054");
   A2->ScDescript.push_back("Eur.Phys.J. C19 (2001) 289-311");

   A2->NScDescript = A2->ScDescript.size();
   A2->Ecms = sqrt(4.*820.*27.5);   // --- fastNLO user: set sqrt(s)
   A2->ILOord = 1;    // --- fastNLO user: power of LO contribution for process
   A2->NDim = 2;      // --- fastNLO user: No of dimensions in which observable is binned
   A2->DimLabel.push_back("E_T"); // --- fastNLO user: label of 1st dimension
   A2->IDiffBin.push_back(2);
   A2->DimLabel.push_back("Q^2"); // --- fastNLO user: label of 2nd dimension
   A2->IDiffBin.push_back(2);

   vector <double> bound;
   bound.resize(2);


   // --- fastNLO user: bin definitions - here in Q2 and ET
   const int nq2bins = 4;
   double q2bins[nq2bins+1] = {150.,200.,300.,600.,5000.};

   const int netbins[nq2bins] = {4, 4, 4, 4};
   double etbins[nq2bins][5] = {
     {7.,11.,18.,30.,50.},
     {7.,11.,18.,30.,50.},
     {7.,11.,18.,30.,50.},
     {7.,11.,18.,30.,50.}
   };

   // --- fastNLO user:
   //     define below the bin width ("binsize") by which
   //     the cross section is divided to obtain the 
   //     (multi-) differential result.
   double binsize = 0.;

   int nbins = 0;   // --- count total No. bins
   for(int i=0;i<nq2bins;i++){
      for(int j=0;j<netbins[i];j++){
 	 nbins += 1;
         bound[0] = etbins[i][j];
         bound[1] = q2bins[i];
         A2->LoBin.push_back(bound);
         bound[0] = etbins[i][j+1];
         bound[1] = q2bins[i+1];
         A2->UpBin.push_back(bound);

         binsize = (etbins[i][j+1]-etbins[i][j]) 
           *(q2bins[i+1]-q2bins[i]);  // binsize in ET, Q2
         A2->BinSize.push_back(binsize);
      }
   }
   printf(" tot. No. observable bins = %d\n",nbins);

   A2->NObsBin = nbins;
   A2->INormFlag = 0;    // --- fastNLO user: default=0 - set =1 if observable is 
			 //     to be normalized by own integral (in 1st dimension)
			 //     see documentation for details and for other options

   // --- fastNLO table block B
   fnloBlockBNlojet *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
   table->CreateBlockB(0,B);
   B->IXsectUnits = 12;   // --- fastNLO user: set to same value as "SetIpublunits"
   B->IDataFlag = 0;
   B->IAddMultFlag = 0;
   B->IContrFlag1 = 1;
   //   B->IContrFlag3 = 0;
   B->NScaleDep = 0;
   B->CodeDescript.push_back("NLOJET++ 4.1.3");  // --- fastNLO user: enter NLOJET++ version
   //B->NCodeDescr = B->CodeDescript.size();
   B->IRef = 0;
   if (nlo || A2->ILOord > 1) {
     B->NSubproc = 3;
   } else {
     B->NSubproc = 2;
     printf("  this job uses 2 subprocesses \n");
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
   //B->NContrDescr = B->CtrbDescript.size();

   B->NPDF = 1;
   B->NPDFPDG.push_back(2212);   // --- fastNLO user: PDG code for hadron
   B->NPDFDim = 0;
   B->NFragFunc = 0;
   B->IPDFdef1 = 2;
   B->IPDFdef2 = 1;
   if(B->NSubproc == 3) {
     B->IPDFdef3 = 3;
   } else {
     B->IPDFdef3 = 2;
     printf("  set IPDFdef3=2 consistent with 2 subprocesses \n");
   }
   B->XNode1.resize(A2->NObsBin);

   // --- initialize variables for WarmUp run
   //B->IWarmUp = 1;       // --- fastNLO user: do the Warm-Up run
   B->IWarmUp = 0;     //                   or do production run(s)
   // - fastNLO user: remember to disable reference-mode in
   //                 Warm-Up run: "doReference = false" (above)
   B->IWarmUpPrint = 30000000;
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
   // 180000000 contributions (!= events) in warm-up run 
   xlim[0]=0.00645, mulo[0]=   7.000, muup[0]=  11.000;
   xlim[1]=0.01155, mulo[1]=  11.000, muup[1]=  18.000;
   xlim[2]=0.02522, mulo[2]=  18.000, muup[2]=  30.000;
   xlim[3]=0.06704, mulo[3]=  30.000, muup[3]=  50.000;
   xlim[4]=0.00738, mulo[4]=   7.000, muup[4]=  11.000;
   xlim[5]=0.01280, mulo[5]=  11.000, muup[5]=  18.000;
   xlim[6]=0.02744, mulo[6]=  18.000, muup[6]=  30.000;
   xlim[7]=0.06727, mulo[7]=  30.000, muup[7]=  50.000;
   xlim[8]=0.00928, mulo[8]=   7.000, muup[8]=  11.000;
   xlim[9]=0.01472, mulo[9]=  11.000, muup[9]=  18.000;
   xlim[10]=0.02920, mulo[10]=  18.000, muup[10]=  30.000;
   xlim[11]=0.06772, mulo[11]=  30.000, muup[11]=  50.000;
   xlim[12]=0.01491, mulo[12]=   7.000, muup[12]=  11.000;
   xlim[13]=0.02039, mulo[13]=  11.000, muup[13]=  18.000;
   xlim[14]=0.03482, mulo[14]=  18.000, muup[14]=  30.000;
   xlim[15]=0.07570, mulo[15]=  30.000, muup[15]=  50.000;
   // --------- fastNLO: Warm-Up run results (end)


   for(int i=0;i<A2->NObsBin;i++){
      int nxtot = 20;
      if (i == ((A2->NObsBin)-1)) nxtot = 21; // Darf's etwas mehr sein?
      B->Nxtot1.push_back(nxtot);
      double hxlim = log10(xlim[i]);  // use exact value from Warm-Up run
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
   //B->NscaleDescript.push_back(B->ScaleDescript[0].size());
   B->Nscalenode.push_back(13); // number of scale nodes for mu

   B->ScaleFac.resize(B->NScaleDim);
   for(int i=0;i<B->NScaleDim;i++){
      B->ScaleFac[i].push_back(1.0);
      if(nlo){
	B->ScaleFac[i].push_back(0.5); // --- fastNLO user: add any number of
	B->ScaleFac[i].push_back(2.0); // additional scale variations as desired
	//B->ScaleFac[i].push_back(4.0);
	//B->ScaleFac[i].push_back(8.0);
      }
      B->Nscalevar.push_back(B->ScaleFac[i].size());
   }

   const double mu0scale = 0.25; // --- variable in H(mu) (in GeV)
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
	      // B->ScaleNode[i][j][k][0] = B->ScaleFac[0][k]*(A2->UpBin[i][B->Iscale[j]] + A2->LoBin[i][B->Iscale[j]])/2.;
               B->ScaleNode[i][j][k][0] = B->ScaleFac[0][k]*(muup[i]+mulo[i])/2.;
               B->HScaleNode[i][j][k][0] = log(log((B->ScaleFac[0][k]*(muup[i]+mulo[i])/2.)/mu0scale));
            }else{
               double llscalelo = log(log((B->ScaleFac[0][k]*mulo[i])/mu0scale));
               double llscalehi = log(log((B->ScaleFac[0][k]*muup[i])/mu0scale));
               for(int l=0;l<B->Nscalenode[j];l++){
                 // later this is the place where the Chebychev nodes will be implemented
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
   
   
   // reference table
   if(doReference){
      fnloBlockBNlojet *refB = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
      // refB->NSubproc = 3;
      if (nlo || A2->ILOord > 1) {
        refB->NSubproc = 3;
      } else {
        refB->NSubproc = 2;
        printf("  this reference job uses 2 subprocesses \n");
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
