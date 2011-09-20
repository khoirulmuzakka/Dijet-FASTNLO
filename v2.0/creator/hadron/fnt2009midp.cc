//
//   fastNLO v2 author code 
//   scenario: fnt2009 
//   Run IIa D0 incl jets
//   for Run II midpoint cone algorithm (no Rsep)
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
   double pTmin = 50., ymin = 0., ymax = 2.4;

   // --- for inclusive jet cross section: loop over all jets
   for(unsigned int i = 1; i <= nj; i++){
      double pt = pj[i].perp(); 
      double rap = fabs(pj[i].rapidity());
      // --- jet in phase space?
      if(pt > pTmin  &&  rap < ymax  &&  rap > ymin) { 

         // - set the renormalization and factorization scale to jet pT
         double mu = pt;

         // --- identify bin number (y,pT)
         int obsbin = -1;
         for(int j = 0; j < A2->GetNObsBin(); j++) {
	    if (pt >= A2->LoBin[j][0] &&  pt < A2->UpBin[j][0] && 
	       rap >= A2->LoBin[j][1] && rap < A2->UpBin[j][1]) {
	      obsbin = j;
	      break;
	    }
	 }
	 
	 // --- fill fastNLO arrays - don't touch this piece of code!
	 if (obsbin >= 0) {
 	    double prefactor = 1./A2->BinSize[obsbin]; 
	    for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
   	       if(table->GetBlockB(k)->GetIRef()>0){
  		  // --- fastNLO user: don't modify the following calls!
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
   table->GetBlockA1()->SetScenName("fnt2009");  // - fastNLO user: set scenario name
   table->GetBlockA1()->SetNcontrib(1);
   table->GetBlockA1()->SetNmult(0);
   table->GetBlockA1()->SetNdata(0);
   table->GetBlockA2()->SetIpublunits(12);  // - fastNLO user: set cross section units
                                            //                 (negative power of ten)

   // --- fastNLO: fill variables for table header block A2
   fnloBlockA2 *A2 =  table->GetBlockA2();

   // --- fastNLO user: up to 20 strings to describe the scenario
   A2->ScDescript.push_back("d2sigma-jet_dpT_dy_(pb_GeV)");
   A2->ScDescript.push_back("xxx");
   A2->ScDescript.push_back("D0_Collaboration");
   //A2->ScDescript.push_back("");

   A2->NScDescript = A2->ScDescript.size();
   A2->Ecms = sqrt(s);
   A2->ILOord = 2;   // --- fastNLO user: power of LO contribution for process
   A2->NDim = 2;     // --- fastNLO user: No of dimensions in which observable is binned
   A2->DimLabel.push_back("pT");  // --- fastNLO user: label of 1st dimension
   A2->IDiffBin.push_back(2);
   A2->DimLabel.push_back("y");   // --- fastNLO user: label of 2nd dimension
   A2->IDiffBin.push_back(2);

   vector <double> bound;
   bound.resize(2);

   // --- fastNLO user: bin definitions - here in |y| and pT
   const int nrapbins = 6;
   double rapbins[nrapbins+1] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

   const int nptbins[nrapbins] = {23, 22, 20, 17, 15, 13};
   double ptbins[6][24] = {
     { 50., 60., 70., 80., 90., 100., 110., 120., 130., 145., 
       160., 180., 200., 220., 240., 265., 295., 325., 360., 400.,
       445., 490., 540., 665.},
     { 50., 60., 70., 80., 90., 100., 110., 120., 130., 145., 
       160., 180., 200., 220., 240., 265., 295., 325., 360., 400.,
       445., 495., 635.},
     { 50., 60., 70., 80., 90., 100., 110., 125., 140., 155., 
       170., 190., 210., 230., 250., 270., 300., 335., 375., 415.,
       520.},
     { 50., 60., 70., 80., 90., 100., 110., 125., 140., 155., 
       170., 190., 215., 240., 265., 290., 325., 415.},
     { 50., 60., 70., 80., 90., 100., 110., 125., 140., 160., 
       175., 190., 210., 235., 260., 320.},
     { 50., 60., 70., 80., 90., 100., 110., 120., 130., 145., 
       160., 175., 200., 230.}
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
         bound[1] = rapbins[i+1];
         A2->UpBin.push_back(bound);
	 //printf(" %d %d  |  %f %f\n",i,j,bound[0],bound[1]);

	 binsize = (ptbins[i][j+1]-ptbins[i][j]) *
	   2.*(rapbins[i+1]-rapbins[i]);  // factor "2." to account for |y| vs. y
	 A2->BinSize.push_back(binsize);
      }
   }
   printf(" tot. No. observable bins = %d\n",nbins);

   A2->NObsBin = nbins;
   A2->INormFlag = 0;   // --- fastNLO user: default=0 - set =1 if observable is 
                        //     to be normalized by own integral (in 1st dimension)
                        //     see documentation for details and for other options

   // --- fastNLO table block B
   fnloBlockBNlojet *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
   table->CreateBlockB(0,B);
   B->IXsectUnits = 12;    // --- fastNLO user: set to same value as "SetIpublunits"
   B->IDataFlag = 0;
   B->IAddMultFlag = 0;
   B->IContrFlag1 = 1;
   B->IContrFlag3 = 0;
   B->CodeDescript.push_back("NLOJET++ 4.1.3");  // --- fastNLO user: enter NLOJET++ version
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
   // 1500000000 contributions (!= events) in warm-up run 
   xlim[0]=0.017804, mulo[0]=50.000000, muup[0]=60.000000;
   xlim[1]=0.021519, mulo[1]=60.000000, muup[1]=70.000000;
   xlim[2]=0.025331, mulo[2]=70.000000, muup[2]=80.000000;
   xlim[3]=0.029159, mulo[3]=80.000000, muup[3]=90.000000;
   xlim[4]=0.033093, mulo[4]=90.000000, muup[4]=100.000000;
   xlim[5]=0.037029, mulo[5]=100.000000, muup[5]=110.000000;
   xlim[6]=0.041132, mulo[6]=110.000000, muup[6]=120.000000;
   xlim[7]=0.045214, mulo[7]=120.000000, muup[7]=130.000000;
   xlim[8]=0.049391, mulo[8]=130.000000, muup[8]=145.00;
   xlim[9]=0.055784, mulo[9]=145.000000, muup[9]=160.00;
   xlim[10]=0.062347, mulo[10]=160.000000, muup[10]=180.00;
   xlim[11]=0.071460, mulo[11]=180.000000, muup[11]=200.000000;
   xlim[12]=0.080749, mulo[12]=200.000000, muup[12]=220.00;
   xlim[13]=0.090477, mulo[13]=220.000000, muup[13]=240.00;
   xlim[14]=0.100598, mulo[14]=240.000000, muup[14]=265.000000;
   xlim[15]=0.113865, mulo[15]=265.000000, muup[15]=295.00;
   xlim[16]=0.130207, mulo[16]=295.000000, muup[16]=325.00;
   xlim[17]=0.147893, mulo[17]=325.000000, muup[17]=360.00;
   xlim[18]=0.169812, mulo[18]=360.000000, muup[18]=400.00;
   xlim[19]=0.197073, mulo[19]=400.000000, muup[19]=445.00;
   xlim[20]=0.230718, mulo[20]=445.000000, muup[20]=490.00;
   xlim[21]=0.267692, mulo[21]=490.000000, muup[21]=540.00;
   xlim[22]=0.314192, mulo[22]=540.000000, muup[22]=665.000000;
   xlim[23]=0.012169, mulo[23]=50.000000, muup[23]=60.000000;
   xlim[24]=0.014782, mulo[24]=60.000000, muup[24]=70.000000;
   xlim[25]=0.017446, mulo[25]=70.000000, muup[25]=80.000000;
   xlim[26]=0.020200, mulo[26]=80.000000, muup[26]=90.000000;
   xlim[27]=0.023050, mulo[27]=90.000000, muup[27]=100.000000;
   xlim[28]=0.025898, mulo[28]=100.000000, muup[28]=110.000000;
   xlim[29]=0.028834, mulo[29]=110.000000, muup[29]=120.00;
   xlim[30]=0.031870, mulo[30]=120.000000, muup[30]=130.000000;
   xlim[31]=0.035005, mulo[31]=130.000000, muup[31]=145.000000;
   xlim[32]=0.039869, mulo[32]=145.000000, muup[32]=160.000000;
   xlim[33]=0.044916, mulo[33]=160.000000, muup[33]=180.000000;
   xlim[34]=0.051970, mulo[34]=180.000000, muup[34]=200.000000;
   xlim[35]=0.059422, mulo[35]=200.000000, muup[35]=220.000000;
   xlim[36]=0.067302, mulo[36]=220.000000, muup[36]=240.00;
   xlim[37]=0.075710, mulo[37]=240.000000, muup[37]=265.000000;
   xlim[38]=0.086991, mulo[38]=265.000000, muup[38]=295.00;
   xlim[39]=0.101747, mulo[39]=295.000000, muup[39]=325.00;
   xlim[40]=0.118194, mulo[40]=325.000000, muup[40]=360.00;
   xlim[41]=0.139832, mulo[41]=360.000000, muup[41]=400.00;
   xlim[42]=0.168376, mulo[42]=400.000000, muup[42]=445.000000;
   xlim[43]=0.206314, mulo[43]=445.000000, muup[43]=495.00;
   xlim[44]=0.255313, mulo[44]=495.000000, muup[44]=635.0;
   xlim[45]=0.008467, mulo[45]=50.000000, muup[45]=60.000000;
   xlim[46]=0.010284, mulo[46]=60.000000, muup[46]=70.000000;
   xlim[47]=0.012240, mulo[47]=70.000000, muup[47]=80.000000;
   xlim[48]=0.014243, mulo[48]=80.000000, muup[48]=90.000000;
   xlim[49]=0.016344, mulo[49]=90.000000, muup[49]=100.000000;
   xlim[50]=0.018535, mulo[50]=100.000000, muup[50]=110.000000;
   xlim[51]=0.020789, mulo[51]=110.000000, muup[51]=125.000000;
   xlim[52]=0.024407, mulo[52]=125.000000, muup[52]=140.000000;
   xlim[53]=0.028270, mulo[53]=140.000000, muup[53]=155.000000;
   xlim[54]=0.032371, mulo[54]=155.000000, muup[54]=170.0;
   xlim[55]=0.036745, mulo[55]=170.000000, muup[55]=190.000000;
   xlim[56]=0.043162, mulo[56]=190.000000, muup[56]=210.0;
   xlim[57]=0.050165, mulo[57]=210.000000, muup[57]=230.000000;
   xlim[58]=0.057996, mulo[58]=230.000000, muup[58]=250.0;
   xlim[59]=0.066735, mulo[59]=250.000000, muup[59]=270.0;
   xlim[60]=0.076578, mulo[60]=270.000000, muup[60]=300.0;
   xlim[61]=0.093772, mulo[61]=300.000000, muup[61]=335.000000;
   xlim[62]=0.116925, mulo[62]=335.000000, muup[62]=375.0;
   xlim[63]=0.146453, mulo[63]=375.000000, muup[63]=415.000000;
   xlim[64]=0.179528, mulo[64]=415.000000, muup[64]=520.0;
   xlim[65]=0.005921, mulo[65]=50.000000, muup[65]=60.000000;
   xlim[66]=0.007316, mulo[66]=60.000000, muup[66]=70.000000;
   xlim[67]=0.008787, mulo[67]=70.000000, muup[67]=80.000000;
   xlim[68]=0.010349, mulo[68]=80.000000, muup[68]=90.000000;
   xlim[69]=0.012022, mulo[69]=90.000000, muup[69]=100.000000;
   xlim[70]=0.013810, mulo[70]=100.000000, muup[70]=110.000000;
   xlim[71]=0.015719, mulo[71]=110.000000, muup[71]=125.000000;
   xlim[72]=0.018868, mulo[72]=125.000000, muup[72]=140.000000;
   xlim[73]=0.022352, mulo[73]=140.000000, muup[73]=155.0;
   xlim[74]=0.026303, mulo[74]=155.000000, muup[74]=170.0;
   xlim[75]=0.030744, mulo[75]=170.000000, muup[75]=190.0;
   xlim[76]=0.037704, mulo[76]=190.000000, muup[76]=215.000000;
   xlim[77]=0.048154, mulo[77]=215.000000, muup[77]=240.0;
   xlim[78]=0.060013, mulo[78]=240.000000, muup[78]=265.0;
   xlim[79]=0.073163, mulo[79]=265.000000, muup[79]=290.0;
   xlim[80]=0.087595, mulo[80]=290.000000, muup[80]=325.0;
   xlim[81]=0.111374, mulo[81]=325.000000, muup[81]=415.0;
   xlim[82]=0.004279, mulo[82]=50.000000, muup[82]=60.000000;
   xlim[83]=0.005372, mulo[83]=60.000000, muup[83]=70.000000;
   xlim[84]=0.006594, mulo[84]=70.000000, muup[84]=80.0;
   xlim[85]=0.007938, mulo[85]=80.000000, muup[85]=90.000000;
   xlim[86]=0.009425, mulo[86]=90.000000, muup[86]=100.000000;
   xlim[87]=0.011100, mulo[87]=100.000000, muup[87]=110.000000;
   xlim[88]=0.012996, mulo[88]=110.000000, muup[88]=125.000000;
   xlim[89]=0.016342, mulo[89]=125.000000, muup[89]=140.000000;
   xlim[90]=0.020422, mulo[90]=140.000000, muup[90]=160.000000;
   xlim[91]=0.026665, mulo[91]=160.000000, muup[91]=175.0;
   xlim[92]=0.031906, mulo[92]=175.000000, muup[92]=190.000000;
   xlim[93]=0.037597, mulo[93]=190.000000, muup[93]=210.0;
   xlim[94]=0.046164, mulo[94]=210.000000, muup[94]=235.000000;
   xlim[95]=0.059699, mulo[95]=235.000000, muup[95]=260.0;
   xlim[96]=0.078238, mulo[96]=260.000000, muup[96]=320.0;
   xlim[97]=0.003234, mulo[97]=50.000000, muup[97]=60.0;
   xlim[98]=0.004208, mulo[98]=60.000000, muup[98]=70.000000;
   xlim[99]=0.005357, mulo[99]=70.000000, muup[99]=80.0;
   xlim[100]=0.006751, mulo[100]=80.000000, muup[100]=90.000000;
   xlim[101]=0.008443, mulo[101]=90.000000, muup[101]=100.0;
   xlim[102]=0.010418, mulo[102]=100.000000, muup[102]=110.000000;
   xlim[103]=0.012601, mulo[103]=110.000000, muup[103]=120.000000;
   xlim[104]=0.015003, mulo[104]=120.000000, muup[104]=130.000000;
   xlim[105]=0.017605, mulo[105]=130.000000, muup[105]=145.000000;
   xlim[106]=0.022106, mulo[106]=145.000000, muup[106]=160.0;
   xlim[107]=0.027890, mulo[107]=160.000000, muup[107]=175.0;
   xlim[108]=0.035582, mulo[108]=175.000000, muup[108]=200.000000;
   xlim[109]=0.056449, mulo[109]=200.000000, muup[109]=230.0;
   // --------- fastNLO: Warm-Up run results (end)


   //printf("* --- xlimits \n");
   for(int i=0;i<A2->NObsBin;i++){
      //int nxtot = 12;
     int nxtot = 20;    // --- test
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
   //B->Nscalenode.push_back(4); // number of scale nodes for pT
   B->Nscalenode.push_back(6); // number of scale nodes for pT

   B->ScaleFac.resize(B->NScaleDim);

   B->ScaleFac[0].push_back(1.0);    // --- fastNLO: central scale (don't change)
   if(nlo){
     //B->ScaleFac[0].push_back(2.0);  // --- fastNLO user: add any number of
     //B->ScaleFac[0].push_back(0.5);  //             additional scale variations
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
   
   
   // reference table
   if(doReference){
      fnloBlockBNlojet *refB = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
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
