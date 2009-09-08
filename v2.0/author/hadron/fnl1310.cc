//
// fastNLO author code for fnl1310:
//     CMS LHC test scenario, E_cms = 10 TeV
//     for fastjet kT algo with D=0.6 in E-scheme
// 
// last modification
// 2009/09/05 KR   First V2.0 Version with multiple bins in pT and rap, table header ok
// 2009/09/01 KR   V2.0 Version
//
//# define DEBUG
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
#include "fnloTable.h"
#include "fnloBlockBNlojet.h"
#include "pdf-cteq6.h"
#include "pdf-hhc-dummy.h"
#include "fj-kt-06.h"

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
  fj_kt_06 jetclus;   // jet algorithm
   
  bounded_vector<lorentzvector<double> > pj;    // the jet structure 
   
  //fastNLO starts here
   
  fnloTable *table;
   
  double nevents;        // No of events calculated so far
  unsigned long nwrite;  // No of events after to write out the table

  string tablefilename; // The table file to write to
  time_t start_time;
   
  bool nlo;
  bool ref;
  void inittable();
  void writetable();
};



user_base_hhc * userfunc() {
  return new UserHHC;
}



void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
#ifdef DEBUG
  cout << "Start fastNLO inputfunc ..." << endl; 
#endif /* DEBUG */
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
#ifdef DEBUG
  cout << "Start fastNLO psinput ..." << endl; 
#endif /* DEBUG */
  //  total c.m. energy squared
  //s =     40000.; // RHIC               200 GeV   
  //s =   3240000.; // TeV Run I         1800 GeV
  //s =   3841600.; // TeV Run II        1960 GeV
  //s =    810000.; // LHC Injection Run  900 GeV
  //s =  49000000.; // LHC First Run     7000 GeV
  s = 100000000.; // LHC Start-up Run 10000 GeV
  //s = 196000000.; // LHC Design Run   14000 GeV

  //   You can use your own phase generator. 
  //   Here we use the default.
  ps = 0;
} 



void UserHHC::initfunc(unsigned int)
{
#ifdef DEBUG
  cout << "Start fastNLO initfunc ..." << endl; 
#endif /* DEBUG */
  // ---- Initialize event counters
  nevents = 0;
  // Set some defaults
  // KR: Change for testing
  if (nwrite==0) nwrite = 1000000;
  //  if (nwrite==0) nwrite = 10000;
    
  // KR: Avoid seg faults, initialize table here and not only in
  //     UserHHC::phys_output where it seems to be too late ?
  // Not required any more if phys_output declared virtual in NLOJet++ 
  //  inittable();

  start_time = std::time(0);
}



void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
#ifdef DEBUG
  cout << "Start fastNLO userfunc ..." << endl; 
#endif /* DEBUG */
  fnloBlockA2 *A2 =  table->GetBlockA2();
  
  double x1 = p[-1].Z()/p[hadron(-1)].Z();
  double x2 = p[0].Z()/p[hadron(0)].Z();

  //----- do the jet analysis -----
  //  pj = jetclus(p,0.7);
  pj = jetclus(p);
  unsigned int nj = pj.upper(); 
  
  //----- loop over all jets -----
  for (unsigned int i = 1; i <= nj; i++) {
    int irap = -1;
    int ipt  = -1;
    double rap  = fabs(pj[i].rapidity());
    double pt   = pj[i].perp(); 
    double drap = 1.0;

    //----- first cut on lowest pt and on y
    int nobs = A2->GetNObsBin();
#ifdef DEBUG
    int iobs = A2->RapIndex[1];
    cout << "iobs = " << iobs << endl;
    cout << "ptmin =  " << A2->LoBin[0][0] << ", pt = " << pt <<
      ", ptmax = " << A2->UpBin[iobs-1][0] << endl;
    cout << "ymin =  " << A2->LoBin[0][1] << ", rap = " << rap <<
      ", rapmax = " << A2->UpBin[nobs-1][1] << endl;
#endif /* DEBUG */
    if ( pt >= A2->LoBin[0][0] &&
	 A2->LoBin[0][1] <= rap && rap < A2->UpBin[nobs-1][1] ) {
      
      //----- determine y and pt bin -----
      int nrap = A2->RapIndex.size();
      for ( int i = 0; i < nrap; i++) {        
	int iobs = A2->RapIndex[i];
	if ( A2->LoBin[iobs][1] <= rap && rap < A2->UpBin[iobs][1] ) {
	  irap = i;
	  drap = 2.0* (A2->UpBin[iobs][1] - A2->LoBin[iobs][1]);
#ifdef DEBUG
 	  cout << "Rap bin: i = " << i <<
 	    " , LoBin i,1 = " << A2->LoBin[iobs][1] <<
 	    " , rap = " << rap <<
 	    " , UpBin i,1 = " << A2->UpBin[iobs][1] << endl;
#endif /* DEBUG */
	  break;
	}
      }
      
      if (irap >= 0) { 
	// find corresponding pt bin number
	int npt = nobs;
	int rapindex = A2->RapIndex.size();
	if (irap+1 < rapindex ) {
	  npt = A2->RapIndex[irap+1]; 
	}
	for (int i = A2->RapIndex[irap]; i < npt; i++) {
	  if (A2->LoBin[i][0] <= pt && pt < A2->UpBin[i][0]) {
	    ipt = i;
#ifdef DEBUG
 	    cout << "pt bin: i = " << i <<
 	      " , LoBin i,0 = " << A2->LoBin[i][0] <<
 	      " , pt = " << pt <<
	      " , UpBin i,0 = " << A2->UpBin[i][0] << endl;
#endif /* DEBUG */
	    break;
	  }
	}
      }
    
      //----- fill fastNLO arrays -----
      if ( ipt>=0 ) {
	for (int k=0; k<table->GetBlockA1()->GetNcontrib(); k++){
	  if (table->GetBlockB(k)->GetIRef()>0) {
	    ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(ipt,x1,x2,pt,amp,pdf,0.5); // 0.5 = factor for bin in |eta| -> eta
	  } else {
	    ((fnloBlockBNlojet*)(table->GetBlockB(k)))->FillEventHHC(ipt,x1,x2,pt,amp,dummypdf,0.5);
	  }
	}
      }        
    }
  }
}



void UserHHC::writetable(){
#ifdef DEBUG
  cout << "Start fastNLO writetable ..." << endl; 
#endif /* DEBUG */
  table->OpenFileRewrite();
  table->WriteBlockA1();
  table->WriteBlockA2();
  for(int i=0;i< table->GetBlockA1()->GetNcontrib();i++){
    table->WriteBlockBDividebyN(i);
  }
  table->CloseFileWrite();
}



void UserHHC::end_of_event(){
#ifdef DEBUG
  cout << "Start fastNLO end_of_event ..." << endl; 
#endif /* DEBUG */
  nevents += 1;
  //-------- store table
  if (( (unsigned long)nevents % nwrite)==0){
    time_t hour, min, time = std::time(0) - start_time;
      
    hour = time/3600L;
    time -= hour*3600L;
    min  = time/60L;
    time -= min*60L;
      
    cout     <<"--->     "
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
#ifdef DEBUG
  cout << "Start fastNLO phys_output ..." << endl; 
#endif /* DEBUG */
  tablefilename.assign(__file_name.c_str());
  tablefilename += ".tab";
   
  //Determine whether we are running LO or NLO
  const char* const file = __file_name.c_str(); 
  if (strstr(file,"born") != NULL) {
    nlo = false;
    cout << "fastNLO: Running at Born level ..." << endl;
  } else {
    if (strstr(file,"nlo") != NULL) {
      nlo = true;
      cout << "fastNLO: Running at NLO level ..." << endl;
    } else {
      //      printf("This module can only be run at Born level or at NLO.\n");
      cout <<
	"fastNLO: This module can only run at Born or NLO level, aborted!"
	   << endl;
      exit(1);
    }
  }

  //Determine whether we are running in reference mode or not
  ref = false;
  if (strstr(file,"ref") != NULL) {
    ref = true;
    cout << "fastNLO: Running in reference mode ..." << endl; 
  }
  
  nwrite = __save;
  inittable();
}



void UserHHC::inittable(){
#ifdef DEBUG
  cout << "Start fastNLO inittable ..." << endl; 
#endif /* DEBUG */
  // Replaced by bool ref set in phys_output via filename string match
  //  const bool doReference = false;
  // KR: Workaround since I couldn't figure out how the filename setting
  //     could work with inittable() in initfunc
  //  string tablefilename = "fastNLO.tab";

  //  total c.m. energy squared
  //double s =     40000.; // RHIC               200 GeV   
  //double s =   3240000.; // TeV Run I         1800 GeV
  //double s =   3841600.; // TeV Run II        1960 GeV
  //double s =    810000.; // LHC Injection Run  900 GeV
  //double s =  49000000.; // LHC First Run     7000 GeV
  double s = 100000000.; // LHC Start-up Run 10000 GeV
  //double s = 196000000.; // LHC Design Run   14000 GeV

  // Set up fastNLO
  table = new fnloTable(tablefilename);

  table->GetBlockA1()->SetScenName("fnl1310");
  table->GetBlockA1()->SetNcontrib(1);
  table->GetBlockA1()->SetNmult(0);
  table->GetBlockA1()->SetNdata(0);
  table->GetBlockA2()->SetIpublunits(15);

  // Give basic description
  fnloBlockA2 *A2 = table->GetBlockA2();
  A2->ScDescript.push_back("d2sigma-jet_dpT_dy_(fb_GeV)");
  A2->ScDescript.push_back("CMS-PAS-QCD-08-001");
  A2->ScDescript.push_back("CMS_Collaboration");
  A2->ScDescript.push_back("Inclusive_Jets");
  A2->ScDescript.push_back("fast-kT_R=0.6");
  A2->NScDescript = A2->ScDescript.size();
  A2->Ecms = sqrt(s);
  A2->ILOord = 2;

  // Define phase space
  A2->NDim = 2;
  A2->DimLabel.push_back("pT_jet_(GeV)");
  A2->IDiffBin.push_back(2);
  A2->DimLabel.push_back("y_jet");
  A2->IDiffBin.push_back(2);

  // Fix phase space binning
  vector <double> bound;
  bound.resize(2);
  
  const unsigned int nrap = 6;
  double rapbins[nrap+1] = {0.00,0.55,1.10,1.70,2.50,3.20,5.00};

  const unsigned int npt[nrap] = {37,36,33,27,19,12};
  double ptbins[38] = {  
    53.,   67.,   81.,   97.,  114.,  133.,  153.,  174.,  196.,  220.,
    245.,  272.,  300.,  330.,  362.,  395.,  430.,  468.,  507.,  548.,
    592.,  638.,  686.,  737.,  790.,  846.,  905.,  967., 1032., 1101.,
    1172., 1248., 1327., 1410., 1497., 1588., 1684., 1784.};
  
  int nobs = 0;
  for (unsigned int i=0; i<nrap; i++) {
    // Memorize start of next rap bin in table
    A2->RapIndex.push_back(nobs);
    for (unsigned int j=0; j<npt[i]; j++) {
      nobs++;
      bound[0] = ptbins[j];
      bound[1] = rapbins[i];
      A2->LoBin.push_back(bound);
      bound[0] = ptbins[j+1];
      bound[1] = rapbins[i+1];
      A2->UpBin.push_back(bound);
    }
  }
  A2->NObsBin = nobs;
  A2->INormFlag = 0;

  // Main table
  fnloBlockB *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
  table->CreateBlockB(0,B);
  B->IXsectUnits = 15;
  B->IDataFlag = 0;
  B->IAddMultFlag = 0;
  B->IContrFlag1 = 1;
  B->IContrFlag3 = 0;
  B->CodeDescript.push_back("NLOJet++_4.0.1");
  B->NCodeDescr = B->CodeDescript.size();
  B->IRef = 0;
  //  if (doReference) {B->IRef = 1;}
  if (ref) {B->IRef = 1;}
  if (nlo) {
    B->CtrbDescript.push_back("NLO");
    B->IContrFlag2 = 2;
    B->IScaleDep = 1;
    B->Npow = A2->ILOord+1;
    B->NSubproc = 7;
  } else {
    B->CtrbDescript.push_back("LO");      
    B->IContrFlag2 = 1;
    B->IScaleDep = 0;
    B->Npow = A2->ILOord;
    B->NSubproc = 6;
  }
  B->NContrDescr = B->CtrbDescript.size();

  B->NPDF = 2;
  B->NPDFPDG.push_back(2212);
  B->NPDFPDG.push_back(2212);
  B->NPDFDim = 1;
  B->NFragFunc = 0;
  B->NFFDim = 0;
  B->IPDFdef1 = 3;
  B->IPDFdef2 = 1;
  B->IPDFdef3 = 1;

  B->XNode1.resize(A2->NObsBin);
  for (int i=0; i<A2->NObsBin; i++) {
    int nxtot = 12;
    B->Nxtot1.push_back(nxtot);
    // Setup the xlimit array - computed from kinematic constraints
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

  B->ScaleDescript[0].push_back("pT_jet");
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
  //  if(doReference){
  if(ref){
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
