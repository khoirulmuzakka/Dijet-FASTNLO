//
// fastNLO v2 author code for fnl2932b:
//     CMS LHC 3-jet Mass Scenario, E_cms = 7 TeV
//     for fastjet anti-kT algo with R=0.7 in E-scheme
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

#include "fj-ak.h"   // fastNLO user: .h file for jet algorithm

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
   lorentzvector<double> Boost(lorentzvector<double>&, threevector<double>&);
   virtual void end_of_event();
   virtual void phys_output(const std::basic_string<char>& fname, unsigned long nsave = 10000UL, bool txt = false);
   struct ThreeJetVariables {
      ThreeJetVariables() {}
      double                        mass, xi[3];
      double                        comgij[3], muij[3];
      double                        cosEKa, cosEKb;
      double                        cosTheta, cosThetaB, phi, chi;
      bool                          ok;
   };

private:
   // --- pdf
   pdf_cteq6 pdf;
   pdf_hhc_dummy dummypdf;

   // --- jet algorithm
   fj_ak jetclus;   // fastNLO user: define jet algorithm (consistent with .h file above)

   bounded_vector<lorentzvector<double> > pj;    // the jet structure

   // --- fastNLO definitions (not for user)
   fnloTable *table;
   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table
   string tablefilename;  // The table file to write to
   time_t start_time;

   // --- fastNLO user FYI: steering of these flags is encoded in NLOJet++ run name (option -n name)
   //     if filename matches "deb", "ref", or "wrm" the respective flag is set to true
   bool doDebug, doReference, doWarmUp;
   bool nlo;

   void inittable();
   void writetable();
   double GetEcms();
   unsigned int GetNj();

   //Add y-bin to be used
   unsigned int ybin;
};

user_base_hhc * userfunc() {
   return new UserHHC;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
   // --- fastNLO user: select the number of jets for your observable
   //nj = 1U;
   //nj = 2U;
   nj = 3U;

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
   s =  49000000.; // LHC First Run     7000 GeV
   //s =  64000000.; // LHC Second Run    8000 GeV
   //s = 100000000.; // LHC Start-up Run 10000 GeV
   //s = 196000000.; // LHC Design Run   14000 GeV

   //   You can use your own phase generator.
   //   Here we use the default.
   ps = 0;
}

// --- fastNLO user: modify jet selection in userfunc (default = cutting in |y| min, |y| max and pt min)
//     (return value must be true for jets to be UNselected)
// fnl2932: use rapidity!
struct fNLOSelector {
   fNLOSelector(double ymin, double ymax, double ptmin):
      _ymin (ymin), _ymax (ymax), _ptmin (ptmin){};
   double _ymin, _ymax, _ptmin;
   bool operator() (const lorentzvector<double> &a) {return ! (_ymin <= abs(a.rapidity()) && abs(a.rapidity()) < _ymax && _ptmin <= a.perp());};
};

// --- fastNLO user: modify jet sorting in userfunc (default = descending in jet pt)
struct fNLOSorter {
   bool operator() (const lorentzvector<double> &a, const lorentzvector<double> &b) {return (a.perp() > b.perp());};
};

void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{

   // --- fastNLO: Don't touch this piece of code!
   fnloBlockA2 *A2 =  table->GetBlockA2();
   double x1 = p[-1].Z()/p[hadron(-1)].Z();
   double x2 = p[0].Z()/p[hadron(0)].Z();

   // --- fastNLO user: Set jet size and run the jet algorithm
   double jetsize = 0.5;
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
   //     and the bin number ("obsbin") which gets passed to
   //     fastNLO's table filling code.
   //     Usually, pT and E are in GeV, but this may be changed.
   //     ATTENTION: Scales must always be in GeV!

   // --- declare and initialize phase space cut variables
   // smallest |(pseudo-)rapidity| for jets to be considered
   const double yjmin  = 0.0;
   // largest |(pseudo-)rapidity| for jets to be considered
   const double yjmax  = 2.5;
   // lowest pT for jets to be considered
   const double ptjmin = 50.;

   // --- select jets in y or eta and ptjmin (failing jets are moved to the end of the jet array pj!)
   static fNLOSelector SelJets(yjmin,yjmax,ptjmin);
   // --- count number of selected jets left at this stage
   size_t njet = std::remove_if(pj.begin(), pj.end(), SelJets) - pj.begin();

   // --- sort selected n jets at beginning of jet array pj, by default decreasing in pt
   static fNLOSorter SortJets;
   std::sort(pj.begin(), pj.begin() + njet, SortJets);

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

   // 3-jet mass requires at least 3 jets
   if (njet > 2) {

      // --- declare and initialize additional cut variables
      // maximal |y| of three leading jets
      //    const double y3jmax = 2.5;
      const double yRange[6] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
      // Derive 3-jet variables
      UserHHC::ThreeJetVariables var;
      var.ok = false;
      const double small = 1.e-12;
      // 3-jet mass
      lorentzvector<double> pcms = pj[1] + pj[2] + pj[3];
      var.mass = pcms.mag();
      if (var.mass < 0.) {cout << "Warning!: Negative Mass" << var.mass << endl;}
      threevector<double> cmsboost = pcms.boostVector();
      if(doDebug) std::cout << "boost " << cmsboost << std::endl;

      std::vector<threevector<double> > p3vec;
      std::vector<lorentzvector<double> > p4vec;
      for (int i=0; i<3; i++) {
         lorentzvector<double> ptemp;
         ptemp = Boost(pj[i+1], cmsboost);
         threevector<double> xyz = ptemp.unit();
         p3vec.push_back(xyz);
         p4vec.push_back(ptemp);
      }

      double sumst=0, sint[3];
      for (int i=0; i<3; i++) {
         int j=i+1;
         if (j == 3) j = 0;
         int k=j+1;
         if (k == 3) k = 0;
         double cost = p3vec[i]*p3vec[j];
         sint[k] = std::sqrt(1.-cost*cost);
         sumst += sint[k];
      }
      double x1s = 2.*sint[0]/sumst;
      double x2s = 2.*sint[1]/sumst;
      double x3s = 2.*sint[2]/sumst;
      double xtmp= std::max(x1s,x2s);
      var.xi[0]  = std::max(xtmp,x3s);
      xtmp       = std::min(x1s,x2s);
      var.xi[2]  = std::min(xtmp,x3s);
      var.xi[1]  = 2. - var.xi[0] - var.xi[2];

      int n[3]={0,1,2};
      if (std::abs(var.xi[0]-x1s) < small) {
         if (std::abs(var.xi[1]-x2s) > small) {
            n[1] = 2; n[2] = 1;
         }
      } else if (std::abs(var.xi[0]-x2s) < small) {
         n[0] = 1;
         if (std::abs(var.xi[1]-x1s) < small) {
            n[1] = 0; n[2] = 2;
         } else {
            n[1] = 2; n[2] = 0;
         }
      } else {
         n[0] = 2;
         if (std::abs(var.xi[1]-x1s) < small) {
            n[1] = 0; n[2] = 1;
         } else {
            n[1] = 1; n[2] = 0;
         }
      }

      int ij = 0;
      for (int i=1; i<3; i++) {
         for (int j=0; j<i; j++) {
            var.comgij[ij] = p3vec[n[i]]*p3vec[n[j]];
            var.muij[ij]   = std::sqrt(0.5*(var.xi[i])*(var.xi[j])*(1.-var.comgij[ij]));
            ij++;
         }
      }

      var.cosEKa = (var.xi[0]+2.*var.xi[1]-2.)/(var.xi[0]);
      lorentzvector<double> pcms1 = p4vec[n[1]]+p4vec[n[2]];
      threevector<double> cmboost1 = pcms1.boostVector();
      lorentzvector<double> pv1 = Boost(p4vec[n[0]], cmboost1);
      lorentzvector<double> pv2 = Boost(p4vec[n[1]], cmboost1);
      double Pv1 = pv1.X()*pv1.X()+pv1.Y()*pv1.Y()+pv1.Z()*pv1.Z();
      double Pv2 = pv2.X()*pv2.X()+pv2.Y()*pv2.Y()+pv2.Z()*pv2.Z();
      double pp = Pv1*Pv2;
      double p12= (pv1.X()*pv2.X()+pv1.Y()*pv2.Y()+pv1.Z()*pv2.Z());
      var.cosEKb =-p12/pp;

      var.cosTheta = p3vec[n[0]].Z();
      var.phi      = std::atan2(p3vec[n[0]].Y(),p3vec[n[0]].X());
      threevector<double> xyz = cross(p3vec[n[0]], p3vec[n[1]]);
      var.cosThetaB = xyz.Z()/xyz.mag();
      threevector<double> xyzb = threevector<double>(-p3vec[n[0]].Y(),p3vec[n[0]].X(),0);
      double coschi = xyz*xyzb/((xyz.mag())*(xyzb.mag()));
      double chib = std::acos(coschi);
      if (xyz.Z()>0) chib  =-chib;
      if (chib < 0)  chib += (2.*M_PI);
      var.chi = chib;
      if ( doDebug ) {
         std::cout << "CM Momentum " << pcms << " Mass " << var.mass << "\n";
         for (int i=0; i<3; i++)
            std::cout << "Vector " << i << " i/p " << pj[i+1] << " Transformed "
                      << p4vec[i] << " mag " << p4vec[i].mag() << " Unit vector "
                      << p3vec[i] << "\n";
         std::cout << "sint[0..2] = " << sint[0] << ", " << sint[1] << ", "
                   << sint[2] << " Sum " << sumst << " x " << x1s << ", " << x2s
                   << ", " << x3s << " xi " << var.xi[0] << ", " << var.xi[1]
                   << ", " << var.xi[2] << "\n";
         std::cout << "n[0,1,2] = " << n[0] << ", "<< n[1] << ", "<< n[2] << "\n";
         for (ij=0; ij<3; ij++)
            std::cout << "Combination " << ij << " cos(Omega) " << var.comgij[ij]
                      << " mu " << var.muij[ij] << "\n";
         std::cout << "Cos(E-K) " << var.cosEKa << " Meth2 " << pcms1 << " Vec "
                   << pv1 << " and " << pv2 << " pp/p12/cos " << pp << ", "
                   << p12 << ", " << var.cosEKb << "\n";
         std::cout << "CosTheta " << var.cosTheta << " phi " << var.phi
                   << " Chi " << var.chi << " Normals " << xyz << " and "
                   << xyzb << " cosThetaB " << var.cosThetaB << " coschi "
                   << coschi << "\n";
      }
      var.ok = true;

      // Maximal (pseudo-)rapidity
      double y3j = max(abs(pj[1].rapidity()),abs(pj[2].rapidity()));
      y3j = max(y3j,abs(pj[3].rapidity()));

      // Average pTs of leading jets
      //    double pT1   = (pj[1].perp()) / 1.0;
      //    double pT12  = (pj[1].perp() + pj[2].perp()) / 2.0;
      //    double pT123 = (pj[1].perp() + pj[2].perp() + pj[3].perp()) / 3.0;

      // pT1
      double pT1 = pj[1].perp();

      // --- Further 3-jet phase space cuts?
      if ( y3j >= yRange[ybin] && y3j < yRange[ybin+1]) {

         // --- set the renormalization and factorization scale to m3j / 2
         double mu = var.mass/2.;

         //      UserHHC::ThreeJetVariables jetv = shape3(njet,pj,true);

         // --- identify bin number (dim1,dim2) here (m3j,y3jmax)
         int obsbin = -1;
         for (int j = 0; j < A2->GetNObsBin(); j++) {
            if (A2->LoBin[j][0] <= var.mass && var.mass < A2->UpBin[j][0] &&
                A2->LoBin[j][1] <= pT1 && pT1 < A2->UpBin[j][1]) {
               obsbin = j;
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
      } // --- end: final selection
   } // --- end: 3-jet+ events only
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
   // --- fastNLO user: set scenario name (no white space)
   A1->SetScenName("fnl2932by0m3");

   // --- fastNLO: fill variables for table header block A2
   fnloBlockA2 *A2 = table->GetBlockA2();
   // --- fastNLO user: set cross section units (negative power of ten)
   A2->SetIpublunits(12);
   // --- fastNLO user: up to 20 strings to describe the scenario
   A2->ScDescript.push_back("d2sigma-3-jet_dM3jd|y|_max_[pb_GeV]");
   A2->ScDescript.push_back("CMS_Collaboration");
   A2->ScDescript.push_back("E_cms=7_TeV");
   A2->ScDescript.push_back("3-Jet_Mass");
   A2->ScDescript.push_back("anti-kT_R=0.5");
   A2->ScDescript.push_back("CMS-PAS-QCD-11-003");
   A2->ScDescript.push_back("provided by:");
   A2->ScDescript.push_back("fastNLO_2.1.0");
   A2->ScDescript.push_back("If you use this table, please cite:");
   A2->ScDescript.push_back("  D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch, arXiv:1109.1310");

   A2->NScDescript = A2->ScDescript.size();
   A2->Ecms = GetEcms();
   A2->ILOord = GetNj(); // --- fastNLO user: power of LO contr. for process (2 for incl. jets, 3 for 3-jet mass)
   A2->NDim = 2;     // --- fastNLO user: no. of dimensions in which observable is binned
   A2->DimLabel.push_back("M3j_[GeV]");  // --- fastNLO user: label of 1st dimension
   A2->IDiffBin.push_back(2);
   A2->DimLabel.push_back("|y|_max");   // --- fastNLO user: label of 2nd dimension
   A2->IDiffBin.push_back(2);

   vector <double> bound;
   bound.resize(2);

   ybin = 0;     // ybin to be studied
   // --- fastNLO user: bin definitions - here in M3j and pT
   const int ndim2bins = 5;
   const double dim2bins[ndim2bins+1] = {110., 190., 300., 360., 500., 5000.};

   cout << endl << "------------------------" << endl;
   cout << "Binning in dimension 2: " << A2->DimLabel[1] << endl;
   cout << "------------------------" << endl;
   for (int i=0; i<ndim2bins+1; i++) {
      cout << "i, dim2bins: " << i << ", " << dim2bins[i] << endl;
   }

   vector< vector<double> > dim1bins;
   const int ndim1bins[5*ndim2bins] = {22, 26, 26, 27, 27, 20, 23, 25, 23, 28,
                                       21, 22, 26, 24, 27, 21, 23, 25, 26, 28,
                                       23, 28, 27, 27, 30};
   const int istrtbin[5*ndim2bins] = {7,  8, 12, 14, 16,  7,  9, 11, 15, 15,
                                      8, 10, 11, 15, 17,  8, 10, 13, 15, 17,
                                      8, 10, 14, 16, 16};
   const double dim0[49] = {
      0.,     50.,   65.,   83.,  103.,  125.,  150.,  177.,  207.,  240.,  275.,
      313.,   354.,  398.,  445.,  495.,  548.,  604.,  664.,  727.,  794.,
      864.,   938., 1016., 1098., 1184., 1274., 1369., 1469., 1573., 1682.,
      1796., 1916., 2041., 2172., 2309., 2452., 2602., 2758., 2921., 3092.,
      3270., 3456., 3650., 3852., 4063., 4283., 4513., 4753. }; // Mass
   //  const int ndim1bins[5*ndim2bins] = {20, 20, 20, 20, 20, 19, 19, 20, 20, 20,
   //      19, 20, 20, 20, 20, 19, 20, 20, 20, 20,
   //      20, 20, 20, 20, 20};
   //  const int istrtbin[5*ndim2bins] = {6,  6,  6,  6,  6,  7,  7,  6,  6,  6,
   //     7,  6,  6,  6,  6,  7,  6,   6,  6,  6,
   //     6,  6,  6,  6,  6};
   //  const double dim0[37] = {
   //  0.500, 0.525, 0.550, 0.575, 0.600, 0.625, 0.650, 0.675, 0.700, 0.725,
   //  0.750, 0.775, 0.800, 0.825, 0.850, 0.875, 0.900, 0.910, 0.920, 0.930,
   //  0.940, 0.950, 0.960, 0.970, 0.980, 0.990, 1.000}; // x0(37)
   //  const double dim0[21] = {
   //  0.500, 0.525, 0.550, 0.575, 0.600, 0.625, 0.650, 0.675, 0.700, 0.725,
   //  0.750, 0.775, 0.800, 0.825, 0.850, 0.875, 0.900, 0.925, 0.950, 0.975,
   //  1.000}; // x1(21)
   //  const int ndim1bins[5*ndim2bins] = {18, 19, 20, 20, 20, 18, 18, 19, 19, 20,
   //      17, 19, 19, 20, 20, 18, 19, 19, 20, 20,
   //      18, 19, 20, 20, 20};
   //  const int istrtbin[5*ndim2bins] = {0,  0,  0,  0,  0,  0,  0,   0,  0,  0,
   //     0,  0,  0,  0,  0,  0,  0,   0,  0,  0,
   //     0,  0,  0,  0,  0};
   //  const double dim0[25] = {
   //  0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.250,
   //  0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750,
   //  0.800, 0.850, 0.900, 0.950, 1.000} // x2(25)
   //  const int ndim1bins[5*ndim2bins] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
   //      10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
   //      10, 10, 10, 10, 10};
   //  const int istrtbin[5*ndim2bins] = {0,  0,  0,  0,  0,  0,  0,   0,  0,  0,
   //     0,  0,  0,  0,  0,  0,  0,   0,  0,  0,
   //     0,  0,  0,  0,  0};
   //  const double dim0[11] = {
   //  0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00}; // eka(11)

   dim1bins.resize(ndim2bins);
   for (int i=0; i<ndim2bins; i++) {
      dim1bins[i].resize(ndim1bins[5*ybin+i]+1);
      for (int j=0; j<ndim1bins[5*ybin+i]+1; j++) {
         dim1bins[i][j] = dim0[istrtbin[5*ybin+i]+j];
      }
   }
   cout << endl << "------------------------" << endl;
   cout << "Binning in dimension 1: " << A2->DimLabel[0] << endl;
   cout << "------------------------" << endl;
   for (int i=0; i<ndim2bins; i++) {
      for (int j=0; j<ndim1bins[5*ybin+i]+1; j++) {
         cout << "i, j, dim1bins: " << i << ", " << j << ", " << dim1bins[i][j] << endl;
      }
   }
   cout << "========================" << endl;

   // --- fastNLO user:
   //     define below the bin width ("binsize") by which
   //     the cross section is divided to obtain the
   //     (multi-) differential result.
   //     default: divide by bin width in dim 1 and dim 2
   //              ATTENTION: Don't forget to include a factor of 2 for abs. rapidity |y| !
   // fnl2932: divide by bin width in M3j and |y|_max


   printf("fastNLO: AAAAAA\n");

   int nbins = 0;   // --- count total No. bins
   for (int i=0;i<ndim2bins;i++){
      for (int j=0;j<ndim1bins[5*ybin+i];j++){
         double binsize = 1.; // --- start each bin with preset value = 1.
         nbins += 1;
         bound[0] = dim1bins[i][j];
         bound[1] = dim2bins[i];
         A2->LoBin.push_back(bound);
         bound[0] = dim1bins[i][j+1];
         bound[1] = dim2bins[i+1];
         A2->UpBin.push_back(bound);
         binsize = binsize
            * (dim1bins[i][j+1]-dim1bins[i][j]) // ... e.g. times dM3j
            * 2. * (dim2bins[i+1]-dim2bins[i]); // ... e.g. times d|y|_max
         A2->BinSize.push_back(binsize);
      }
   }
   printf("fastNLO: Total no. of observable bins = %d\n\n",nbins);

   A2->NObsBin = nbins;
   A2->INormFlag = 0;   // --- fastNLO user: default=0 - set =1 if observable is
   //     to be normalized by own integral (in 1st dimension)
   //     see documentation for details and for other options

   // --- fastNLO table block B
   fnloBlockBNlojet *B = new fnloBlockBNlojet(table->GetBlockA1(),table->GetBlockA2());
   table->CreateBlockB(0,B);
   B->SetNlojetDefaults();

   B->IXsectUnits = 12;    // --- fastNLO user: set to same value as "SetIpublunits"
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
   B->NPDFPDG.push_back(2212);   // --- fastNLO user: PDG code for 1st hadron
   B->NPDFPDG.push_back(2212);  // --- fastNLO user: PDG code for 2nd hadron
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

   // --- initialize variables for WarmUp run
   // KR: Now set via filename string match to "wrm"
   //B->IWarmUp = 1;     // --- fastNLO user: do the Warm-Up run
   B->IWarmUp = 0;     //                   or do production run(s)
   if ( doWarmUp ) {B->IWarmUp = 1;}
   // KR: This is caught in an error condition now
   // - fastNLO user: remember to disable reference-mode in
   //                 Warm-Up run: "doReference = false" (above)
   B->IWarmUpPrint = 1000000;
   //B->IWarmUpPrint = 10000;
   B->xlo.resize(A2->NObsBin);
   B->scalelo.resize(A2->NObsBin);
   B->scalehi.resize(A2->NObsBin);

   // --- arrays for extreme x and (default) scale values (computed in Warm-Up run)
   double xlim[A2->NObsBin];
   double mulo[A2->NObsBin];
   double muup[A2->NObsBin];

   // --- fastNLO user: before running a new scenario for the first time,
   //     the following block (between "start" and "end") should be
   //     completely removed. The first run must be a "Warm-Up Run" (IWarmUp=1)
   //     which produces an initialization block as output. This should be copied
   //     and pasted below. These initialization values must be used for all
   //     production jobs (IWarmUp=0) for a given scenario (otherwise the result
   //     tables can not be merged).
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
         cout << "         Trying to find and use x limits included in scenario author code ..." << endl;
         // --------- fastNLO: Warm-Up run results (start)
         // fastNLO user: paste warm-up run results here ...
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
   B->NScaleDim = 1; // one variable used in scales: dijet pT average
   B->Iscale.push_back(0);  // mur=mur(pT), pT = index 0
   B->Iscale.push_back(0);  // muf=muf(pT), pT = index 0
   B->ScaleDescript.resize(B->NScaleDim);

   B->ScaleDescript[0].push_back("M3j/2_[GeV]");
   //B->Nscalenode.push_back(4); // number of scale nodes for pT
   B->Nscalenode.push_back(6); // number of scale nodes for pT

   B->ScaleFac.resize(B->NScaleDim);

   B->ScaleFac[0].push_back(1.0);    // --- fastNLO: central scale (don't change)
   if(nlo){
      B->ScaleFac[0].push_back(0.5);    // --- fastNLO user: add any number of
      B->ScaleFac[0].push_back(2.0);    //             additional scale variations
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

lorentzvector<double> UserHHC::Boost(lorentzvector<double>& a, threevector<double>& b) {
   double bx = b.X(); double by = b.Y(); double bz = b.Z();
   double _M_x = a.X(); double _M_y = a.Y();
   double _M_z = a.Z(); double _M_t = a.T();

   double b2 = bx*bx + by*by + bz*bz;
   double gamma = 1.0/std::sqrt(1.0 - b2);
   double bp = bx*_M_x + by*_M_y + bz*_M_z;
   double gamma2 = (b2 > 0.0 ? (gamma - 1.0)/b2 : 0.0);

   double M_x = _M_x + gamma2*bp*bx - gamma*bx*_M_t;
   double M_y = _M_y + gamma2*bp*by - gamma*by*_M_t;
   double M_z = _M_z + gamma2*bp*bz - gamma*bz*_M_t;
   double M_t = gamma*(_M_t - bp);
   return lorentzvector<double>(M_x, M_y, M_z, M_t);
}
