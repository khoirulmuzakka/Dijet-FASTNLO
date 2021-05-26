//
// fastNLO v2 creator code for fnl5350eta0:
//     CMS LHC Inclusive Jets Scenario, E_cms = 5.02 TeV
//     for fastjet anti-kT algo with R=0.2, 0.3, 0.4 in E-scheme
//
//
// ============== fastNLO user: ===================================
// To create your own scenario, it is recommended to take
// this code, make a copy and edit the relevant changes.
// Important:
// Edit only those lines which are labeled as "fastNLO user"
// and refer to the documentation ("fastNLO creator code in
// NLOJet++") for a detailed explanation of the parameters
// and variables.
// If a code fragment is not explicitely labeled as "fastNLO user",
// it is likely that a modification will interfere with
// the fastNLO routines.
// Please keep the order of all statements in inittable
// in order to guarantee a working code.
//
// This file contains the following routines:
//   inputfunc    (-> user edits)
//   psinput      (-> user edits)
//   userfunc     (-> user edits)
//   inittable    (-> user edits)
//   initfunc     (don't touch)
//   writetable   (don't touch)
//   end_of_event (don't touch)
//   phys_output  (don't touch)
//   GetEcms      (don't touch)
//   GetNj        (don't touch)
//
// Implementing a new scenario requires to edit:
//  - the jet algorithm ("#include" statement and assignment of "jetclus")
//  - number of jets (inputfunc)
//  - center-of-mass energy (psinput)
//  - compute observable, determine bin no. (userfunc)
//  - declare all variables for table, define bin boundaries (inittable, etc.)
//
// ================================================================

//------ DON'T TOUCH THIS PART! ------
#include <iostream>
#include <bits/hhc-phasespace.h>
#include <bits/hhc-process.h>
#include <bits/hhc-jetfunc.h>

//----- used namespaces -----
using namespace nlo;
using namespace std;

//----- declaration of the user defined functions -----
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

// --- fastNLO user: include the header file for the jet algorithm
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

   // --- fastNLO user: define the jet algorithm (consistent with the header file above)
   fj_ak jetclus;

   bounded_vector<lorentzvector<double> > pj;    // the jet structure

   // --- fastNLO definitions (not for user)
   fnloTable *table;
   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table
   string tablefilename;  // The table file to write to
   time_t start_time;

   // --- fastNLO user FYI: use of the -n option of NLOJet++
   //     The steering of these flags is encoded in the NLOJet++ run name (option -n name).
   //     If the name matches "deb", "ref", or "wrm", the respective flag is set to true.
   bool doDebug, doReference, doWarmUp;
   bool nlo;

   void inittable();
   void writetable();
   double GetEcms();
   unsigned int GetNj();
};

user_base_hhc * userfunc() {
   return new UserHHC;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
   // --- fastNLO user: select the number of jets of the LO process for your observable,
   //                   e.g. 2 for incl. jets, 3 for 3-jet mass
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
   //s =     40000.; // RHIC               200 GeV
   //s =   3240000.; // TeV Run I         1800 GeV
   //s =   3841600.; // TeV Run II        1960 GeV
   //s =    810000.; // LHC Injection Run  900 GeV
   //s =   5569600.; // LHC Initial Run   2360 GeV
   //s =   7617600.; // LHC HIpp base Run 2760 GeV
   s =  25200400.; // LHC HIpp base Run 5020 GeV
   //s =  49000000.; // LHC First Run     7000 GeV
   //s =  64000000.; // LHC Second Run    8000 GeV
   //s = 100000000.; // LHC Start-up Run 10000 GeV
   //s = 196000000.; // LHC Design Run   14000 GeV

   //   You can use your own phase generator.
   //   Here we use the default.
   ps = 0;
}

// --- fastNLO user: modify the jet selection in userfunc (default = cutting in |y| min, |y| max and pt min)
//                   (the return value must be true for jets to be UNselected)
// fnl5350: use pseudorapidity eta
struct fNLOSelector {
   fNLOSelector(double ymin, double ymax, double ptmin):
      _ymin (ymin), _ymax (ymax), _ptmin (ptmin){};
   double _ymin, _ymax, _ptmin;
   bool operator() (const lorentzvector<double> &a) {return ! (_ymin <= abs(a.prapidity()) && abs(a.prapidity()) < _ymax && _ptmin <= a.perp());};
};

// --- fastNLO user: modify the jet sorting in userfunc (default = descending in jet pt)
struct fNLOSorter {
   bool operator() (const lorentzvector<double> &a, const lorentzvector<double> &b) {return (a.perp() > b.perp());};
};

void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{

   // --- fastNLO: Don't touch this piece of code!
   fnloBlockA2 *A2 =  table->GetBlockA2();
   double x1 = p[-1].Z()/p[hadron(-1)].Z();
   double x2 = p[0].Z()/p[hadron(0)].Z();

   // --- fastNLO user: in this scenario run jet algo with three different jet sizes R
   const unsigned int ndim2bins = 29;
   const unsigned int ndim1bins = 3;
   const double Rjet[ndim1bins] = { 0.2, 0.3, 0.4 };
   for (unsigned int k=0; k<ndim1bins; k++) {

      // --- fastNLO user: set the jet size and run the jet algorithm
      double jetsize = Rjet[k];
      pj = jetclus(p,jetsize);
      unsigned int nj = pj.upper();

      // --- give some debug output before selection and sorting
      if ( doDebug ) {
         for (unsigned int i=1; i<=nj; i++) {
            double pti  = pj[i].perp();
            double yi   = pj[i].rapidity();
            double etai = pj[i].prapidity();
            cout << "before cuts: jet # i, pt, y, eta: " << i << ", " << pti << ", " << yi << ", " << etai << endl;
         }
      }

      // --- check on maximal no. of jets: 4 (should never be more in NLOJet++)
      if (nj > 4) {
         cout << "fastNLO: ERROR! This scenario is not suited for " << nj <<
            " jets. Aborted!" << endl;
         exit(1);
      }

      // --- fastNLO user:
      //     Here is your playground where you compute your observable
      //     and the bin number ("obsbin"), which gets passed to
      //     fastNLO's table filling code.
      //     Usually, pT and E are in GeV, but this may be changed.
      //     ATTENTION: Scales must always be in GeV!

      // --- declare and initialize phase space cut variables
      // can partially be taken from table binning for this scenario
      // smallest |(pseudo-)rapidity| for jets to be considered
      const double yjmin  = 0.0;
      // largest |(pseudo-)rapidity| for jets to be considered
      const double yjmax  = 0.3;
      // lowest pT for jets to be considered
      const double ptjmin = A2->LoBin[0][0];

      // --- select jets in y or eta and ptjmin (failing jets are moved to the end of the jet array pj!)
      static fNLOSelector SelJets(yjmin,yjmax,ptjmin);
      // --- count number of selected jets left at this stage
      size_t njet = std::remove_if(pj.begin(), pj.end(), SelJets) - pj.begin();

      // --- sort selected n jets at beginning of jet array pj, by default decreasing in pt
      // fnl5350: Not required for inclusive jets
      //  static fNLOSorter SortJets;
      //  std::sort(pj.begin(), pj.begin() + njet, SortJets);

      // --- give some debug output after selection and sorting
      if ( doDebug ) {
         cout << "# jets before and after phase space cuts: nj, njet = " << nj << ", " << njet << endl;
         cout << "phase space cuts: yjmin, yjmax, ptjmin: " << yjmin << ", " << yjmax << ", " << ptjmin << endl;
         for (unsigned int i=1; i<=njet; i++) {
            double pti  = pj[i].perp();
            double yi   = pj[i].rapidity();
            double etai = pj[i].prapidity();
            cout << "after cuts: jet # i, pt, y, eta: " << i << ", " << pti << ", " << yi << ", " << etai << endl;
         }
      }

      // Analyze inclusive jets in jet loop
      for (unsigned int i = 1; i <= njet; i++) {

         // Get jet quantities
         double pt  = pj[i].perp();

         // --- set the renormalization and factorization scale to jet pT
         double mu = pt;

         // --- identify bin number (dim2,dim1) here (pT,R)
         int obsbin = -1;
         for (unsigned int j = 0; j < ndim2bins; j++) {
            if (A2->LoBin[j][0] <= pt  && pt  < A2->UpBin[j][0]) {
               obsbin = j + k*ndim2bins;
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
         } // --- end: fill fastNLO array
      } // --- end: jet loop
   } // --- end: jet size loop
} // --- end: fastNLO user playground

void UserHHC::inittable(){
   // Decide whether to include a reference table
   // (for precision studies, not for production jobs)
   // Set via filename string match to "ref", see also -n option of NLOJet++
   // To set manually use these lines and disable string match
   //const bool doReference = true;
   //const bool doReference = false;

   // --- set up fastNLO
   table = new fnloTable(tablefilename);

   // --- fastNLO: fill variables for table header block A1
   fnloBlockA1 *A1 = table->GetBlockA1();
   A1->SetHeaderDefaults();
   // --- fastNLO user: set the scenario name (no white space)
   A1->SetScenName("fnl5350eta0");

   // --- fastNLO: fill variables for table header block A2
   fnloBlockA2 *A2 = table->GetBlockA2();
   // --- fastNLO user: set the cross section units in barn (negative power of ten)
   A2->SetIpublunits(12);
   // --- fastNLO user: write up to 20 strings to describe the scenario
   A2->ScDescript.push_back("d2sigma-jet_dpTdeta_[pb_GeV]");
   A2->ScDescript.push_back("CMS_Collaboration");
   A2->ScDescript.push_back("Inclusive_Jet_pT");
   A2->ScDescript.push_back("anti-kT_R=0.2,0.3,0.4");
   A2->ScDescript.push_back("CMS-PAS-HIN");
   A2->ScDescript.push_back("provided by:");
   A2->ScDescript.push_back("fastNLO_2.1.0");
   A2->ScDescript.push_back("If you use this table, please cite:");
   A2->ScDescript.push_back("  D. Britzger, K. Rabbertz, F. Stober, M. Wobisch, arXiv:1208.3641");
   A2->NScDescript = A2->ScDescript.size();

   A2->Ecms = GetEcms();
   A2->ILOord = GetNj();

   // --- fastNLO user: no. and names of dimensions, from outer to inner loop, in which the observable is binned
   A2->NDim = 2;
   A2->DimLabel.push_back("pT_[GeV]");
   A2->IDiffBin.push_back(2);
   A2->DimLabel.push_back("R");
   A2->IDiffBin.push_back(2);

   vector <double> bound;
   bound.resize(2);

   // --- fastNLO user: define the binning
   const int ndim1bins = 3;
   const double dim1bins[ndim1bins+1] = { 0.15, 0.25, 0.35, 0.45 };

   const int ndim2bins[ndim1bins] = { 29, 29, 29 };

   cout << endl << "--------------------------------" << endl;
   cout << "Binning in outer-most dimension: " << A2->DimLabel[1] << endl;
   cout << "--------------------------------" << endl;
   for (int i=0; i<ndim1bins+1; i++) {
      cout << "i, dim1bins: " << i << ", " << dim1bins[i] << endl;
   }

   vector< vector<double> >dim2bins;
   dim2bins.resize(ndim1bins);
   for (int i=0; i<ndim1bins; i++) {
      dim2bins[i].resize(ndim2bins[i]+1);
   }
   const double dim0[30] = {
      22.  ,   27.,   33.,   39.,   47.,   55.,   64.,   74.,   84.,   97.,
      114. ,  133.,  153.,  174.,  196.,  220.,  245.,  272.,  300.,  330.,
      362. ,  395.,  430.,  468.,  507.,  548.,  592.,  638.,  790.,  967. };
   for (int i=0; i<ndim1bins; i++) {
      for (int j=0; j<ndim2bins[i]+1; j++) {
         dim2bins[i][j] = dim0[j];
      }
   }

   cout << endl << "--------------------------------" << endl;
   cout << "Binning in inner-most dimension: " << A2->DimLabel[0] << endl;
   cout << "--------------------------------" << endl;
   for (int i=0; i<ndim1bins; i++) {
      for (int j=0; j<ndim2bins[i]+1; j++) {
         cout << "i, j, dim2bins: " << i << ", " << j << ", " << dim2bins[i][j] << endl;
      }
   }
   cout << "================================" << endl;

   // --- fastNLO user:
   //     Define below the bin width ("binsize") by which the cross section is divided
   //     to obtain the (multi-) differential result.
   //     Default: divide by bin width in dim. 1 and dim. 2
   //     ATTENTION: Don't forget to include a factor of 2 e.g. for abs. rapidity |y| !
   //
   //     fnl5350eta0: Divide by bin width in pT (inner dimension)
   //              Do NOT divide by bin width in outer dimension which is jet size R
   //              DO divide by bin width in eta = 2 * 0.3
   //
   int nbins = 0;   // --- count total No. bins
   for (int i=0;i<ndim1bins;i++){
      for (int j=0;j<ndim2bins[i];j++){
         double binsize = 1.; // --- start each bin with preset value = 1.
         nbins += 1;
         bound[0] = dim2bins[i][j];
         bound[1] = dim1bins[i];
         A2->LoBin.push_back(bound);
         bound[0] = dim2bins[i][j+1];
         bound[1] = dim1bins[i+1];
         A2->UpBin.push_back(bound);
         binsize = binsize
            * (dim2bins[i][j+1]-dim2bins[i][j]) // ... times dpT
            * 2. * 0.3; // ... times deta
         A2->BinSize.push_back(binsize);
      }
   }
   printf("fastNLO: Total no. of observable bins = %d\n\n",nbins);
   A2->NObsBin = nbins;

   // --- fastNLO user:
   //     Set this flag to 1 if the observable is to be normalized by its own integral
   //     (in 1st dimension), see documentation for details. Default is 0.
   A2->INormFlag = 0;

   // --- fastNLO table block B
   fnloBlockBNlojet *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
   table->CreateBlockB(0,B);
   B->SetNlojetDefaults();
   B->IXsectUnits = A2->GetIpublunits();
   B->NScaleDep = 0;

   B->IRef = 0;
   if (nlo || A2->ILOord > 2) {
      B->NSubproc = 7;
   } else {
      B->NSubproc = 6;
   }
   printf("fastNLO: This job uses %d subprocesses!\n",B->NSubproc);
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
   // --- fastNLO user: PDG codes for 1st hadron and 2nd hadron
   B->NPDFPDG.push_back(2212);
   B->NPDFPDG.push_back(2212);
   B->NPDFDim = 1;
   B->NFragFunc = 0;
   B->NFFDim = 0;
   B->IPDFdef1 = 3;
   B->IPDFdef2 = 1;
   if(B->NSubproc == 7) {
      B->IPDFdef3 = 2;
   } else {
      B->IPDFdef3 = 1;
   }
   printf("fastNLO: Set IPDFdef3 = %d, consistent with %d subprocesses.\n",B->IPDFdef3,B->NSubproc);
   B->XNode1.resize(A2->NObsBin);

   // --- fastNLO user FYI:
   //     Before calculating final results the binning in x and in the scale has
   //     to be optimized. This is done in so-called warm-up runs. This can
   //     be steered via a name string match to "wrm", see also -n option of NLOJet++.
   //     To set manually use (B->IWarmUp = 1) to do the Warm-Up run or
   //     (B->IWarmUp = 0) to do production run(s). Do not forget to
   //     disable the automatic setting via the string match in this case.
   //
   // --- initialize variables for WarmUp run
   //B->IWarmUp = 1;
   B->IWarmUp = 0;
   if ( doWarmUp ) {B->IWarmUp = 1;}
   // --- fastNLO user FYI:
   //     Remember to disable the reference-mode during Warm-Up runs:
   //     "doReference = false" (above) if setting manually.
   //     Otherwise this is caught automatically in an error condition.
   B->IWarmUpPrint = 1000000;
   //B->IWarmUpPrint = 10000;
   B->xlo.resize(A2->NObsBin);
   B->scalelo.resize(A2->NObsBin);
   B->scalehi.resize(A2->NObsBin);

   // --- arrays for extreme x and (default) scale values (computed in Warm-Up run)
   double xlim[A2->NObsBin];
   double mulo[A2->NObsBin];
   double muup[A2->NObsBin];

   // --- fastNLO user FYI:
   //     Before running a new scenario for the first time, the following block between
   //     "Warm-Up run results (start)" and "Warm-Up run results (end)" must be
   //     completely removed if existing already. The first run must then be
   //     a "Warm-Up Run" (IWarmUp = 1), which produces an initialization block as output file.
   //     This output file is subsequently read by a production job and must be shipped with it.
   //     Alternatively, the initialization block can be copied and pasted into
   //     this scenario, again between "Warm-Up run results (start)" and "Warm-Up run results (end)".
   //     The same initialization values must be used for all production jobs (IWarmUp = 0)
   //     for a given scenario. Otherwise the result tables can not be merged.
   //
   if ( doWarmUp ) {
      cout << endl << "fastNLO: Initializing x limits for warm-up run ..." << endl;
      for(int i=0;i<A2->NObsBin;i++){
         xlim[i] = 1.1e-07, mulo[i] = 3.0, muup[i]=9.9e10; // - safe initializations
      }
   } else {
      cout << endl << "fastNLO: Initializing x limits ..." << endl;
      FILE * infile;
      infile = fopen("fastNLO-warmup.dat","r");
      if ( ! infile ) {
         cout << "fastNLO: WARNING! Could not read x limits from file: fastNLO-warmup.dat" << endl;
         cout << "         Trying to find and use x limits included in scenario creator code ..." << endl;
         // --------- fastNLO: Warm-Up run results (start)
         // fastNLO user: paste warm-up run results here when not using extra file ...
         // --------- fastNLO: Warm-Up run results (end)
         // Safety check:
         // Count no. of xlim values > 0 (should be equal to NObsBin!)
         int nxlim = 0;
         for (int i=0; i<A2->NObsBin; i++) {
            if (xlim[i] > 0.) {nxlim++;}
         }
         if (nxlim != A2->NObsBin ) {
            cerr << "fastNLO: ERROR! Could not find proper x limits in scenario code!" << endl;
            cerr << "         Do a warm-up run first." << endl;
            exit(1);
         }
      } else {
         cout << "fastNLO: Reading x limits from file fastNLO-warmup.dat ..." << endl;
         char line[256];
         // Ignore first documentation line
         if ( ! fgets(line,sizeof(line),infile) ) {
            cerr << "fastNLO: ERROR! Reading empty file: fastNLO-warmup.dat" << endl;
            exit(1);
         }
         printf("fastNLO: Read first, documentation line from file: fastNLO-warmup.dat\n");
         printf("%s",line);
         // Now read all limits
         int i = 0;
         while ( fgets(line,sizeof(line),infile) ) {
            if (i == A2->NObsBin ) {
               cerr << "fastNLO: ERROR! Number of x limits exceed NObsBin: i = " << i << ", NObsBin = " << A2->NObsBin << endl;
               exit(1);
            }
            int nvar = sscanf(line,"      xlim[ %*u ] = %le , mulo[ %*u ] = %le , muup[ %*u ] = %le ;",
                              &xlim[i],&mulo[i],&muup[i]);
            if (nvar != 3) {
               cerr << "fastNLO: ERROR! x limits line did not match expected format:" << endl;
               printf("%s",line);
               exit(1);
            }
            i++;
         }
         if (i != A2->NObsBin ) {
            cerr << "fastNLO: ERROR! Number of x limits read != NObsBin: i = " << i << ", NObsBin = " << A2->NObsBin << endl;
            exit(1);
         }
      }
   }
   // --- print initialized values
   cout << endl << "fastNLO: Print initialized x and mu limits:" << endl;
   for (int i=0; i<A2->NObsBin; i++) {
      printf("      xlim[ %u ] = %e , mulo[ %u ] = %e , muup[ %u ] = %e ;\n",
             i,xlim[i],i,mulo[i],i,muup[i]);
   }

   for(int i=0;i<A2->NObsBin;i++){
      int nxtot = 15;
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
   B->NScaleDim = 1; // one variable used in scales
   B->Iscale.push_back(0);  // mur=mur(pT), pT = index 0
   B->Iscale.push_back(0);  // muf=muf(pT), pT = index 0
   B->ScaleDescript.resize(B->NScaleDim);

   // --- fastNLO user: give the defined process scale a name and units
   B->ScaleDescript[0].push_back("pT_jet_[GeV]");
   // --- fastNLO user: minimal number of scale nodes is 4
   B->Nscalenode.push_back(6); // number of scale nodes

   B->ScaleFac.resize(B->NScaleDim);

   B->ScaleFac[0].push_back(1.0);    // --- fastNLO: central scale (don't change)
   if(nlo){
      // --- fastNLO user: add any number of additional muf scale variations as desired
      B->ScaleFac[0].push_back(0.5);
      B->ScaleFac[0].push_back(2.0);
      //B->ScaleFac[0].push_back(0.25);
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
      fnloBlockBNlojet *refB = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
      if (nlo || A2->ILOord > 2) {
         refB->NSubproc = 7;
      } else {
         refB->NSubproc = 6;
      }
      printf("         This reference job uses %d subprocesses!\n",B->NSubproc);
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
//------ END OF USER DEFINED PARTS, NO USER EDITS BELOW ------

//------ DON'T TOUCH THIS PART! ------
void UserHHC::initfunc(unsigned int)
{
   // --- Initialize event counters
   nevents = 0;
   // Set some defaults
   if (nwrite==0) nwrite = 5000000;
   start_time = std::time(0);
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
      printf ("fastNLO: No. events: %.3G writing table ...\n",nevents);
      cout.flush();
      for (int k=0;k<table->GetBlockA1()->GetNcontrib();k++){
         table->GetBlockB(k)->Nevt = (long long int)nevents;
      }
      writetable();
      printf("fastNLO: Table written.\n");
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
      printf("fastNLO: This is a LO run!\n");
   }else{
      if(strstr(file,"nlo")!=NULL){
         nlo = true;
         printf("fastNLO: This is a NLO run!\n");
      }else{
         printf("fastNLO: ERROR! This module can only be run at Born level or at NLO.\n");
         exit(1);
      }
   }

   // --- determine whether this is a debug, reference, or warm-up run
   doDebug = false;
   if (strstr(file,"deb")!=NULL) {
      doDebug = true;
      printf("fastNLO: This is a debug run. Attention, huge output!\n");
   }
   doReference = false;
   if (strstr(file,"ref")!=NULL) {
      doReference = true;
      printf("fastNLO: This is a reference run!\n");
   }
   doWarmUp = false;
   if (strstr(file,"wrm")!=NULL) {
      doWarmUp = true;
      printf("fastNLO: This is a warm-up run!\n");
      if ( ! nlo ) {
         printf("fastNLO: WARNING! Warm-up runs are better done at NLO!\n");
      }
   }
   if ( doWarmUp && doReference ) {
      printf("fastNLO: ERROR! Warm-up and reference runs cannot be done simultaneously:\n");
      printf("         doWarmUp = %i, doReference = %i\n",doWarmUp,doReference);
      exit(2);
   }

   nwrite = __save;
   inittable();
}

unsigned int UserHHC::GetNj(){
   unsigned int nj = 0, nu = 0 ,nd = 0;
   inputfunc(nj,nu,nd);
   return nj;
}

double UserHHC::GetEcms(){
   double ecms = 0;
   psinput(NULL,ecms);
   return sqrt(ecms);
}
