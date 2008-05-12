//------ DON'T TOUCH THIS PART! ------
#include <bits/photo-phasespace.h>
#include <bits/photo-process.h>
#include <bits/photo-jetfunc.h>
#include <nlo++-module_add.h>

//----- used namespaces -----
using namespace nlo;
using namespace std;


//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
void psinput(phasespace_photo *, double&, double&);
user_base_photo * userfunc();

typedef unsigned long int (*module_add_type)(bool, const list<basic_string<char> >&, const basic_string<char>&);
extern  module_add_type module_add;


//----- array of the symbols symbols -----
extern "C"{
  
  struct { 
	const char *name;
	void *address;
  } user_defined_functions[] = 
  {
	//   process index: hhc for gamma + p --> jets (photoproduction direct photon)
	{"procindex", (void *) "photodir"},
  
	//   input function 
	{"inputfunc", (void *) inputfunc},
  
	//   phase space input function 
	{"psinput", (void *) psinput},
  
	//   user defined functions
	{"userfunc",  (void *) userfunc},
  
	//   module to generate the readable result
	{"main_module_add", (void *) module_add},
  
	//  end of the list
	{0, 0}
  };
}
//------ END OF THE DO-NOT-TOUCH-PART ------



//------ USER DEFINED PART STARTS HERE ------
#include <algorithm>
#include <kT_clus.h>
#include "pdf-cteq6.h"
#include "pdf-dummy.h"
#include "fnloTable.h"
#include "fnloBlockBNlojet.h"


//    Defines new sample type to do felavor decomposation in the initial state
const char *mySample_label[4] = {"g", "u", "d", "total"};
typedef weight<4U, mySample_label> mySample;

namespace nlo {
  template<>
  struct weight_conversion<weight_photo, mySample>
    : public std::unary_function<weight_photo, mySample>
  {
    mySample operator()(const weight_photo& x)
    {
      mySample res;
      
      for(unsigned int i = 0; i < 3; i++)
		res[3] += (res[i] = x[i]);

      return res;
    }
  };
}

//   Here we specify some base user objects with the new sample type
//     (see photo-jetfunc.h and included header files therein) 
typedef basic_user<jetfunc_photo, void,        mySample, weight_conversion> myUser0d;
typedef basic_user<jetfunc_photo, double,      mySample, weight_conversion> myUser1d;
typedef basic_user<jetfunc_photo, histpoint1d, mySample, weight_conversion> myUser1h;
typedef basic_user<jetfunc_photo, histpoint2d, mySample, weight_conversion> myUser2h;


class UserPhoto : public basic_user_set<myUser0d, myUser1h, myUser2h> 
{
public:
  //   init and user function
   void initfunc(unsigned int);
   void userfunc(const event_photo&, const amplitude_photo&);

   virtual void end_of_event();  
   virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);
  
private:
  //   pdf
  pdf_cteq6 pdf;
  pdf_dummy dummypdf;

  // algorithms
  kT_clus_long jetclus;
   
  //  private types
  typedef lorentzvector<double> _Lv;

  // the jet structore
  bounded_vector<_Lv> cj, pj_lab; 
  bounded_vector<unsigned int> jet;
  
   struct pT_sort {
      bool operator()(const _Lv& p1, const _Lv& p2) const {
         return p1.perp2() > p2.perp2();
      }
   };
   
   struct E_sort {
      bool operator()(const _Lv& p1, const _Lv& p2) const {
         return p1.T() > p2.T();
      }
   };

   //fastNLO starts here

   fnloTable *table;
   
   // ===== variables for the b-cubic interpolation =====
   // - the relative distances to the four nearest bins
   vector<double> cm ; 
   // - the weights for the cubic eigenfunctions (1-dim)
   vector<double> cefm; 

   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table

   string tablefilename; // The table file to write to
   bool textoutput; // If true, the table is written in plain ASCII instead of BASE64 encoded doubles (later for XML)
   amplitude_photo::integral_type itype; // Born, NLO etc.  
   time_t start_time;
   
   bool nlo;
   void inittable();
   void writetable();
  
};
  

//----- defines the module to sum up the results of the different runs -----
module_add_type module_add = 
main_module_add<basic_user_result<myUser0d::distbook_type, myUser1h::distbook_type, myUser2h::distbook_type> >; 

user_base_photo * userfunc() {
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
  eh = 820.0;
  
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


void UserPhoto::userfunc(const event_photo& p, const amplitude_photo& amp)
{
   fnloBlockA2 *A2 =  table->GetBlockA2();

   //----- H1 cuts -----
   double pt = 0.0;
   double pTmin = 21.0, etamin = -1.0, etamax = 2.5;

   //----- photon momentum fraction -----
   double y = (p[-1]*p[hadron(0)])/(p[hadron(-1)]*p[hadron(0)]);
   if(y < 0.3 || y > 0.65) return; 

   //----- Bjorken x -----
   double x = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);
   
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
               //               double thisscale = (A2->UpBin[ptbin][0] + A2->LoBin[ptbin][0])/2.;
               if(table->GetBlockB(k)->GetIRef()>0){
                  ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventPhoto(ptbin,x,pt,amp,pdf);
               }else{
                   ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventPhoto(ptbin,x,pt,amp,dummypdf);
               }
            }
         }        
      }
   }

}

void UserPhoto::end_of_event(){
   // let NLOJET++ store its results
   basic_user_set<myUser0d, myUser1h, myUser2h>::end_of_event();
   
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
//    // Suppress output of NLOJET++ files
   basic_user_set<myUser0d, myUser1h, myUser2h>::phys_output("",2000000000,false);

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
   textoutput = __txt;
   nwrite = __save;
   inittable();

}

void UserPhoto::inittable(){

   const bool doReference = true;

   //Set up fastNLO
   table = new fnloTable(tablefilename);

   table->GetBlockA1()->SetScenName("fnh2101d");
   table->GetBlockA1()->SetNcontrib(1);
   table->GetBlockA1()->SetNmult(0);
   table->GetBlockA1()->SetNdata(0);
   table->GetBlockA2()->SetIpublunits(12);

   fnloBlockA2 *A2 =  table->GetBlockA2();
   A2->ScDescript.push_back("Inclusive jet cross sections in photoproduction");
   A2->ScDescript.push_back("H1 Collaboration");
   A2->ScDescript.push_back("DESY 02-225");
   A2->ScDescript.push_back("EPJC ?");
   A2->NScDescript = A2->ScDescript.size();
   A2->Ecms = sqrt(4.*820.*27.5);
   A2->ILOord = 1;
   A2->NDim = 1;
   A2->DimLabel.push_back("E_T");
   A2->IDiffBin.push_back(2);

   vector <double> bound;
   bound.resize(1);
  
   bound[0] = 21.; A2->LoBin.push_back(bound);
   bound[0] = 28.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 35.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 42.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 52.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 62.; A2->LoBin.push_back(bound);  A2->UpBin.push_back(bound);
   bound[0] = 75.; A2->UpBin.push_back(bound);
   A2->NObsBin = 6;

   A2->INormFlag = 0;

   // main table
   fnloBlockB *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
   table->CreateBlockB(0,B);
   B->IXsectUnits = 12;
   B->IDataFlag = 0;
   B->IAddMultFlag = 0;
   B->IContrFlag1 = 1;
   B->IContrFlag2 = 1;
   B->IContrFlag3 = 0;
   if(nlo){
      B->CtrbDescript.push_back("NLO");
   }else{
      B->CtrbDescript.push_back("LO");      
   }
   B->CtrbDescript.push_back("direct");
   B->NContrDescr = B->CtrbDescript.size();
   B->CodeDescript.push_back("NLOJET++ 4.0.1");
   B->NCodeDescr = B->CodeDescript.size();
   B->IRef = 0;
   if(nlo){
      B->IScaleDep = 1;
      B->Npow = 2;
   }else{
      B->IScaleDep = 0;
      B->Npow = 1;
   }

   B->NPDF = 1;
   B->NPDFPDG.push_back(2212);
   B->NPDFDim = 0;
   B->NFragFunc = 0;
   B->NSubproc = 3;
   B->IPDFdef1 = 2;
   B->IPDFdef2 = 1;
   B->IPDFdef3 = 1;

   const int nxtot = 40;
   const double xlim = 5e-4;
   B->XNode1.resize(A2->NObsBin);
   for(int i=0;i<A2->NObsBin;i++){
      B->Nxtot1.push_back(nxtot);
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
   B->ScaleDescript[0].push_back("center of E_T  bin");
   B->NscaleDescript.push_back(B->ScaleDescript[0].size());
   B->Nscalenode.push_back(5); // number of scale nodes

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
               B->ScaleNode[i][j][k][0] = (A2->UpBin[i][B->Iscale[j]] + A2->LoBin[i][B->Iscale[j]])/2.;
            }else{
               for(int l=0;l<B->Nscalenode[j];l++){
                  B->ScaleNode[i][j][k][l] = A2->LoBin[i][B->Iscale[j]] +
                     ((double)l/(double)(B->Nscalenode[j]-1)) *
                     ( A2->UpBin[i][B->Iscale[j]] - A2->LoBin[i][B->Iscale[j]] );
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
            B->SigmaTilde[i][k][l].resize(B->Nxtot1[i]);
            for(int m=0;m<B->Nxtot1[i];m++){
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

void UserPhoto::writetable(){
   table->OpenFileRewrite();
   table->WriteBlockA1();
   table->WriteBlockA2();
   for(int i=0;i< table->GetBlockA1()->GetNcontrib();i++){
      table->WriteBlockBDividebyN(i);
   }
   table->CloseFileWrite();
}
