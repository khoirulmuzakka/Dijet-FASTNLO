//
// fastNLO author code for fnl2442a:
//     CMS LHC dijet mass ratio scenario, E_cms = 7 TeV
//     for fastjet anti-kT algo with R=0.7 in E-scheme
//
//
// last modification
//
//------ DON'T TOUCH THIS PART! ------
#include <phasespace.h>
#include <process.h>
#include <jetfunc.h>
#include <qcdlib.h>

#include <iomanip>              // for ASCII output for table

//----- used namespaces -----
using namespace nlo;
using namespace std;


//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&, double&);
user_hhc * userfunc();

//----- array of the symbols symbols -----
struct { 
   const char *name;
   void *address;
} user_defined_functions[] = 
   {
      //   process index: 3 --> hadron-hadron --> jets
      {"procindex", (void *) "3"},
    
      //   input function 
      {"inputfunc", (void *) inputfunc},
    
      //   user defined functions
      {"userfunc",  (void *) userfunc},
    
      //  end of the list
      {0, 0}
   };
//------ USER DEFINED PART STARTS HERE ------
#include "fj-ak-07.h"
#include "cteq6.h"

class UserHHC : public user_hhc
{
 public:
   //   init and user function
   void initfunc(unsigned int);
   void userfunc(const event_hhc&, const amplitude_hhc&);
   virtual void phys_output(const std::basic_string<char>& __file_name, 
			   unsigned long __save = 10000UL, bool __txt = false);
   void end_of_event();  

 private:

   bool nlo;       // Is the job running at LO or NLO?
   // binning
   
   unsigned int iref;     //  switch for reference mode
   unsigned int refscale; //  select which scale variaton is used in ref-table

   double unitfactor;   // factor to convert from nb to units in publication
   int nrap;       // No of rapidity bins 
   double *raphigh;  // array for rapidity boundaries
   int *npt;       // No of pT bins in each y range
   vector< vector<double> >pthigh;   // array for pT boundaries
   double ptlow;    // lowest pt considered 
   
   int nscalevar;             // number of scale variations (mu_r,mu_f) in NLO
   vector <double> murscale;  // overall scale factor for renormalization scale
   vector <double> mufscale;       // overall scale factor for fact. scale
   int nscalebin;              // No. of  bins in the scale
   vector< vector< vector<double> > >murval; // array for renorm. scale values
   vector< vector< vector<double> > >mufval; // array for fact. scale values
   vector< vector< vector< vector<double> > > >murvaltrans; //...ln(ln(mu/mu0)
   vector< vector< vector< vector<double> > > >mufvaltrans; //...ln(ln(mu/mu0) 
   double mu0scale; // Reference scale for log/log(mu/mu0) transformation

   int nxtot;      // no of xbins 
   vector< vector<double> >xlimit; // array for lower x limits
   vector< vector<double> >hxlim; // array for function h at xlimit
   vector< vector<double> >xsmallest; // array for smallest actual x values
   vector< vector<double> >ptsmallest; // array for smallest actual pT values
   vector< vector<double> >ptlargest; // array for largest actual pT values
   vector< vector<double> >pthi; // array for largest pT prediction
   vector< vector<double> >ptlo; // array for smallest pT prediction
   //    array for the weights MAIN ARRAY
   vector <vector <vector< vector < vector < vector <weight_hhc> > > > > >weights; 
   
   // ===== variables for the bi-cubic x-interpolation =====
   // - the relative distances to the four nearest bins
   vector<double> cmax ;   vector<double> cmin ; 
   // - the weights for the cubic eigenfunctions (1-dim)
   vector<double> cefmin ; vector<double> cefmax ; 
   // - the weights for the bi-cubic eigenfunctions (2-dim)
   vector< vector<double> > bicef;

   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table

   pdf_cteq6 pdf;  //   pdf
   fj_ak_07 jetclus;   // jet algorithm
 
   bounded_vector<lorentzvector<double> > pj;    // the jet structure 
   basic_string<char> tablefilename; // The table file to write to
   bool textoutput; // If true, the table is written in plain ASCII instead of BASE64 encoded doubles (later for XML)
   amplitude_hhc::integral_type itype; // Born, NLO etc.

   time_t start_time;

   void writetable();
};

user_hhc * userfunc() {
   return new UserHHC;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd, double& s)
{
   //  number of jets 
   //nj = 1U;
   nj = 2U;
   //nj = 3U;

   //  total c.m. energy squared
   //s =     40000.; // RHIC               200 GeV
   //s =   3240000.; // TeV Run I         1800 GeV
   //s =   3841600.; // TeV Run II        1960 GeV
   //s =    810000.; // LHC Injection Run  900 GeV
   //s =   5569600.; // LHC Initial Run   2360 GeV
   s =  49000000.; // LHC First Run     7000 GeV
   //s = 100000000.; // LHC Start-up Run 10000 GeV
   //s = 196000000.; // LHC Design Run   14000 GeV

   //  number of the up and down type flavours
   nu = 2U;
   nd = 3U;
} 

void UserHHC::initfunc(unsigned int)
{
   unsigned int nj;
   unsigned int nu;
   unsigned int nd;
   double      s;
   inputfunc(nj,nu,nd,s);

   // ********************************************************
   // ********** switch for reference mode on/off ************
   iref = 0;       //  switch for reference mode
   //iref = 1;       //  switch for reference mode
                  //   0: standard fastNLO table only
                 //    1: include 2nd "reference table" (a_s/PDFs)
   refscale = 2;   // which of the scalevariations is used in ref-table?

   //unitfactor = 1000000.0;  // for fb
   unitfactor = 1000.0;  // for pb   <<< use for CDF dijet mass
   //unitfactor = 1.0;  // for nb    

   //Set up binning  
   // First dimension (histogram numbers xxxxRxx), usually rapidity
   // 
   // # of bins
   nrap = 2;
   double rapb[3] = { 0.0, 1.1, 2.5};
   // In reference mode: Double no. of bins
   nrap = nrap*(iref+1);
   
   // Array for bin boundaries
   raphigh = new double[nrap+1];
   for ( unsigned int i=0; i<nrap/(iref+1)+1; i++) {
     raphigh[i] = rapb[i];
   }
   // In reference mode: Copy high end of bin boundaries
   if ( iref==1 ) {
     for (unsigned int i=0; i<nrap/(iref+1); i++) {
       raphigh[i+nrap/2+1] = raphigh[i+1];
     }
   }

   //DEBUG
   //    for (int i=0; i<nrap+1; i++) {
   //      cout << "i, rapb, raphigh: " << i << ", " << rapb[i] << ", " << raphigh[i] << endl;
   //    }
   //DEBUGEND

   // Second dimension (histogram x axis), usually pT
   // Here:  DijetMass
   // 
   // # of bins npt per irap bin of first dimension
   int nptb[2] = {50, 50};
   npt = new int[nrap];
   for (unsigned int i=0; i<nrap/(iref+1); i++) {
     npt[i] = nptb[i]; 
   }
   // In reference mode: Copy # of bins definition
   if ( iref==1 ) {
     for (unsigned int i=0; i<nrap/(iref+1); i++) {
       npt[i+nrap/2] = npt[i];
     }
   }

   pthigh.resize(nrap);
   for (int i=0; i<nrap; i++) {
     pthigh[i].resize(npt[i]+1);
   }

   // Array for bin boundaries
   double ptb[51] = { 156.0, 176.0, 197.0, 220.0, 244.0, 270.0, 296.0, 325.0, 354.0, 386.0, 
		      419.0, 453.0, 489.0, 526.0, 565.0, 606.0, 649.0, 693.0, 740.0, 788.0, 
		      838.0, 890.0, 944.0, 1000.0, 1058.0, 1118.0, 1181.0, 1246.0, 1313.0, 
		      1383.0, 1455.0, 1530.0, 1607.0, 1687.0, 1770.0, 1856.0, 1945.0,
		      2037.0, 2132.0, 2231.0, 2332.0, 2438.0, 2546.0, 2659.0, 2775.0, 2895.0,
		      3019.0, 3147.0, 3279.0, 3416.0, 3558.0 };
   for (unsigned int i=0; i<nrap/(iref+1); i++) {
     for (int j=0; j<npt[i]+1; j++) { 
       pthigh[i][j] = ptb[j];
     }
   }
   // In reference mode: Copy high end of bin boundaries
   if ( iref==1 ) {
     for ( unsigned int i=0; i<nrap/(iref+1); i++ ) {
       for ( int j=0; j<npt[i]+1; j++) {
         pthigh[i+nrap/2][j] = pthigh[i][j];
       }
     }
   }

   //DEBUG
   //    for (int i=0; i<nrap; i++) {
   //      for (int j=0; j<npt[i]+1; j++) {
   //        cout << "i, j, ptb, pthigh: " << i << ", " << j << ", " << ptb[j] << ", " << pthigh[i][j] << endl;
   //      }
   //    }
   //DEBUGEND

   // lowest pT value for jets to be considered
   ptlow = 30.;

   // Binning in x
   nxtot = 20;

   // no of scalebins, linear interpolation in between for now
   nscalebin = 2;
   mu0scale = 0.25; // 0.25 GeV as reference scale (HERA!)

   // NLO scale variations - for scales in GeV
   nscalevar = 4;
   murscale.resize(nscalevar);
   mufscale.resize(nscalevar);
   murscale[0] = 0.25; mufscale[0] = 0.25;
   murscale[1] = 0.5;  mufscale[1] = 0.5;
   murscale[2] = 1.0;  mufscale[2] = 1.0;
   murscale[3] = 2.0;  mufscale[3] = 2.0;

   xlimit.resize (nrap);
   hxlim.resize (nrap);
   xsmallest.resize (nrap);     // test: find smallest x values
   ptsmallest.resize (nrap);     // test: find smallest pT values in mass bins
   ptlargest.resize (nrap);     // test: find largest pT values in mass bins
   pthi.resize (nrap);     //
   ptlo.resize (nrap);     //
   murval.resize (nrap);
   mufval.resize (nrap);
   murvaltrans.resize (nrap);
   mufvaltrans.resize (nrap);
 
   weights.resize (nrap);
   for( int j = 0; j < nrap; j++) { 
      xlimit[j].resize(npt[j]);
      hxlim[j].resize(npt[j]);
      xsmallest[j].resize(npt[j]);
      ptsmallest[j].resize(npt[j]);
      ptlargest[j].resize(npt[j]);
      pthi[j].resize(npt[j]);
      ptlo[j].resize(npt[j]);
      murval[j].resize (npt[j]);
      mufval[j].resize (npt[j]);
      murvaltrans[j].resize (npt[j]);
      mufvaltrans[j].resize (npt[j]);
      weights[j].resize(npt[j]);
      for( int k = 0; k < npt[j]; k++) {
         // Setup the scale bins
         murval[j][k].resize(nscalebin);
	 mufval[j][k].resize(nscalebin);
         murvaltrans[j][k].resize (nscalevar);
         mufvaltrans[j][k].resize (nscalevar);
         for( int l = 0; l < nscalevar; l++) {
            murvaltrans[j][k][l].resize (nscalebin);
            mufvaltrans[j][k][l].resize (nscalebin);
         }        

         // Setup the weights array
         weights[j][k].resize(nxtot);
         for( int l = 0; l < nxtot; l++) {
            weights[j][k][l].resize(l+1); // half matrix xmin,xmax: (n^2+n)/2
            // scale variation
            for(int i = 0; i < l+1; i++){
               weights[j][k][l][i].resize(nscalevar);
	       for(int m = 0; m < nscalevar; m++){
		 weights[j][k][l][i][m].resize(nscalebin);
	       }
            }
         }
	 
         // - Setup the xlimit array - computed from kinematic constraints
	 //  ----> need good formula for dijets
	 // KR: My x limits
	 double ymax  = raphigh[j+1];
	 double dymax = 2.0*ymax;
	 double ptred = 1./sqrt(2.)/(sqrt(exp(dymax)) + sqrt(exp(-dymax)));
	 
	 // define reasonable upper and lower pT boundaries for each bin
	 // KR: This is for |y| <= 1.0 @ LO
	 //         ptlo[j][k]  = 0.325*pthigh[j][k];
	 //         pthi[j][k]  = 0.502*pthigh[j][k+1];
	 // KR: For |y| <= 1.3 @ NLO (= LO / sqrt(2))
	 ptlo[j][k]  = ptred*pthigh[j][k];
         pthi[j][k]  = 0.502*pthigh[j][k+1];

	 double xtmin = 2.0 * ptlo[j][k] / sqrt(s);
	 double xmin  = exp(-ymax) * xtmin;
	 // MW:
	 //       double xmin = 1.0;
	 // 	 double xminmax = pthigh[j][k]*pthigh[j][k] / s;
	 //          double xmin2 = 1.0/2.71828 *pthigh[j][k] / sqrt(s); // for y<1
	 // 	 if (xminmax > xmin2) {
	 // 	   xmin = xminmax; }
	 // 	 else {
	 // 	   xmin = xmin2; }
	 xlimit[j][k] = xmin * 1.0;

         // ---- safety factors for E-scheme: None
	 xlimit[j][k] = xlimit[j][k]*1.0;
         hxlim[j][k]= -sqrt(-log10(xlimit[j][k]));
         xsmallest[j][k]=0.999;
         ptsmallest[j][k]=999.9;
         ptlargest[j][k]=0.0;

         // - Setup the murval and mufval arrays
	 double dpt = (pthi[j][k]-ptlo[j][k]);
         murval[j][k][0] = ptlo[j][k] +dpt*0.12;
         murval[j][k][1] = pthi[j][k] -dpt*0.4;
	 mufval[j][k][0] = murval[j][k][0];
         mufval[j][k][1] =  murval[j][k][1];

         for( int l = 0; l < nscalevar; l++) {
	   murvaltrans[j][k][l][0]=log(log(murscale[l]*murval[j][k][0]/mu0scale)); 
	   murvaltrans[j][k][l][1]=log(log(murscale[l]*murval[j][k][1]/mu0scale));
	 }
      }
   }


   if (iref==1)      // -> in reference mode: copy all definitions
     for(int j=0;j<nrap/2;j++){
      for( int k = 0; k < npt[j]; k++) {
	 xlimit[j+nrap/2][k] = xlimit[j][k];
         hxlim[j+nrap/2][k]=  hxlim[j][k];
         for(int l=0;l<nscalebin;l++){
	   murval[j+nrap/2][k][l]= murval[j][k][l]; 
	   mufval[j+nrap/2][k][l]= mufval[j][k][l]; 
	   for( int m = 0; m < nscalevar; m++) {
	     murvaltrans[j+nrap/2][k][m][l] = murvaltrans[j][k][m][l];
	     mufvaltrans[j+nrap/2][k][m][l] = mufvaltrans[j][k][m][l];
	   }
	 }
      }
   }

   // print x-limit values at the begin of the job
   printf ("(rapidity, pt) array for this job:\n");
   printf ("#rap #pt xlimit pt_high rap_high \n");
   for( int j = 0; j < nrap; j++) {
      for( int k = 0; k < npt[j]; k++) {
         printf("%3d %3d   %8.6f %8.3f %8.1f \n",j,k,
		xlimit[j][k],pthigh[j][k],raphigh[(j+1)]);
      }
   }

   // ===== variables for the bi-cubic interpolation =====
   // - the relative distances to the four nearest bins
   cmax.resize (4); cmin.resize (4);
   // - the weights for the cubic eigenfunctions (1-dim)
   cefmin.resize (4); cefmax.resize (4);
   // - the weights for the bi-cubic eigenfunctions (2-dim)
   bicef.resize (4);
   for( int i = 0; i < 4; i++) {
     bicef[i].resize (4);
   }

   // ---------------------------------------------------------------
   // ---- Initialize event counters
   nevents = 0;
   // Set some defaults
   if (nwrite==0) nwrite = 5000000;
   if (tablefilename=="") tablefilename = "fastnlotable.raw";

   // Say Hello
   cout << "  " << endl;
   cout << "   *******************************************" << endl;
   cout << "    fastNLO    - initialization" << endl;
   cout << "    Scenario fnl2442a:" << endl;
   cout  <<  "      CMS LHC test scenario, E_cms = 7 TeV,"  <<  endl;
   cout  <<  "      for fastjet anti-kT algo with R=0.7 in E-scheme"  <<  endl; 
   cout << " " << endl;
   cout << "        table file " << tablefilename << endl;
   cout << "        store table after " << nwrite << " events" << endl;
   cout << "        sqrt(s)= " << sqrt(s) << endl;
   cout << "        No. x-bins: " << nxtot << endl;
   cout << "        No. rapidity regions: " << nrap << endl;
   cout << "        No. of pT bins in each rapidity region:" << endl;
   for( int j = 0; j < nrap; j++) {
      cout<<"          rap "<<j<<": "<<npt[j]<<endl;
   }
   cout<<"        No. of scale variations in NLO: "<<nscalevar<<endl;
   for( int j = 0; j < nscalevar; j++) {
     cout<<"          "<<j<<":   (mur/pT) "<<murscale[j]<<
	"      (muf/pT) "<<mufscale[j]<<endl;
   }
   cout<<"  "<<endl;

   if (iref==1) {     // -> in reference mode: make output
     cout<<"  "<<endl;
     cout<<"      *****************************************"<<endl;
     cout<<"      running in  >> reference mode << (iref=1)"<<endl;
     cout<<"      -> filling second table including"<<endl;
     cout<<"         alphas and PDFs in coefficients"<<endl;
     cout<<"         for reference-table use scale No. "<<refscale<<endl;
     cout<<"          -- for precision studies --"<<endl; }
   else {
     cout<<"          -- reference mode disabled --"<<endl;
   }
   cout<<"  "<<endl;
   cout<<"   *******************************************"<<endl;
   cout<<"        "<<endl;

   if (iref==1) {      // -> in reference mode: define PDFs/alphas
     pdf.mode(pdf_cteq6::nlo); pdf.loop(2);
   }

   start_time = ::time(0);
}


void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
   //  typedef lorentzvector<double> _Lv;
   //-------- general stuff
   //  double s = 2.0*(p[hadron(0)]*p[hadron(-1)]);

   //---------------------------------------------------
   //--------- start event processing ------------------
   //---------------------------------------------------

   //----- variables: CM energy and x-values -----------
   double xmin=0.0, xmax=0.0;
   double x1 = p[-1].Z()/p[hadron(-1)].Z();
   double x2 = p[0].Z()/p[hadron(0)].Z();
   if (x1 > x2) {
      xmax = x1;
      xmin = x2; }
   else {
      xmax = x2;
      xmin = x1;
   }
   // - to reweight the Eigenfunctions (universal PDF)
   double reweight = 1/sqrt(xmin)/sqrt(xmax)*(1-0.99*xmin) * (1-0.99*xmax);  
   double hxmin  = -sqrt(-log10(xmin));
   double hxmax  = -sqrt(-log10(xmax));
   double hxone   = 0.0;

   // LO or NLO?
   itype = amp.integral();

   //----- do the jet analysis -----
   pj = jetclus(p);
   int nj = pj.upper();

   // Check on maximal no. of jets: 3
   if (nj > 3) {
     cout << "ERROR: This scenario is not suited for " << nj <<
       " jets. Aborted!" << endl;
     exit(1);
   }

   int njet = 0;
   int ijet[4]     = {0, 0, 0, 0};
   double ptjet[4] = {0.,0.,0.,0.};
   double yjet[4]  = {10.,10.,10.,10.};
   int nrmax = nrap/(iref+1); // - important for iref==1: different No loops
   if (nj > 0) {
// Initialize pointers to the jets, check min pT and max |eta|
     for (int i=1; i<=nj; i++) {
       double pti = pj[i].perp();
       double yi  = abs(pj[i].prapidity());
       //DEBUG
       //       cout << "ijet, pti, yi: " << i << ", " << pti << ", " << yi << endl;
       //       cout << "ptlow = " << ptlow << ", raphigh = " << raphigh[nrmax] << endl;
       //DEBUGEND
       // QUESTION 1: Do we look for the two leading jets in |eta|<2.5 or |eta|<??? 
       // QUESTION 2: Do we have an explicit minimal pT requirement? 
       if (pti > ptlow && yi < raphigh[nrmax]) {
	 ijet[i]  = 1;
	 ptjet[i] = pti; 
	 yjet[i]  = yi;
	 njet++;
       }
     }
   }

   //DEBUG
   //   cout << "nj, njet = " << nj << ", " << njet << endl;
   //DEBUGEND
   if (njet >= 2) {
     // Get first two jets 
     int ij1 = 0;
     int ij2 = 0;
     for (int i=1; i<=nj; i++) {
       if (ijet[i] == 1) {
	 if (ij1 == 0) {
	   ij1 = i;
	 } else {
	   if (ij2 == 0) {
	     ij2 = i;
	   }
	 }
       }
     }

     // Find the two leading jets
     // Order pt1 >= pt2
     double pt1 = pj[ij1].perp();
     double pt2 = pj[ij2].perp();
     double y1  = abs(pj[ij1].prapidity());
     double y2  = abs(pj[ij2].prapidity());
     if (pt2 > pt1) {
       int itmp     = ij1;
       double pttmp = pt1;
       double ytmp  = y1;
       ij1 = ij2;
       pt1 = pt2;
       y1  = y2;
       ij2 = itmp;
       pt2 = pttmp;
       y2  = ytmp;
     }
     // For 3 jets with pt3 > pt2, exchange 2 and 3
     if (njet > 2 && pj[3].perp() > pt2) {
       ij2 = 3;
       pt2 = pj[ij2].perp();
       y2  = abs(pj[ij2].prapidity());
     }
     // No further check whether pt3 > pt1, i.e. here order not important anymore
     //DEBUG
     //     cout << "ij1, pt1, y1: " << ij1 << ", " << pt1 << ", " << y1 << endl;
     //     cout << "ij2, pt2, y2: " << ij2 << ", " << pt2 << ", " << y2 << endl;
     //DEBUGEND

     // Derive dijet variables
     double mjj = sqrt((pj[ij1].T()+pj[ij2].T())*(pj[ij1].T()+pj[ij2].T())
		       -(pj[ij1].X()+pj[ij2].X())*(pj[ij1].X()+pj[ij2].X())
		       -(pj[ij1].Y()+pj[ij2].Y())*(pj[ij1].Y()+pj[ij2].Y())
		       -(pj[ij1].Z()+pj[ij2].Z())*(pj[ij1].Z()+pj[ij2].Z()));
     
     // Determine minimal and maximal pseudorapdity
     double yjjmin = min(y1,y2);
     double yjjmax = max(y1,y2);

     if (mjj > pthigh[0][0]) {
       // --- Later this variable will be the ren./fact. scale
       //     double ptmax = pt1;    // either pTmax
       double ptmax = (pt1+pt2)/2.0; // or the average dijet pT

       // --- Determine eta and pt (= dijet mass) bin 
       // KR: Note, BOTH leading jets are required to be inside
       //     the SAME |rapidity| interval in this analysis!
       // KR: Do normalize to binwidth in eta ...!
       double binwidth = 1.0;
       int rapbin = -1;
       int nloop = nrap;   // - important for iref==1: different No loops
       if (iref == 1) nloop = nrap/2; // -> in reference mode
       for (int j = 0; j < nloop; j++) {        
	 if (raphigh[j] <= yjjmin && yjjmax < raphigh[j+1]) {
	   rapbin=j;
	   binwidth = binwidth * 2.0*(raphigh[(j+1)] - raphigh[j]);
	   break;
	 }
       }
       if (rapbin >=0) {              
	 int ptbin  = -1;
	 for (int j = 0; j < npt[rapbin]; j++) {
	   if (pthigh[rapbin][j] <= mjj && mjj < pthigh[rapbin][j+1]) {
	     ptbin=j;
	     binwidth = binwidth * (pthigh[rapbin][(j+1)]-pthigh[rapbin][j]);
	     break;
	   }
	 }
       
	 //---------- compute weight, fill fastNLO array
	 if ( ptbin>=0 ) {
	   // test if x_min in event is lower than x_limit
	   //      if yes -> make big warning!!! 
	   //      -> need to change x_limit values
	   if (xmin<xlimit[rapbin][ptbin]){
	     printf("Warning: xmin (%f) < xlimit (%f) at pt=%f GeV y=%f \n ",
		    xmin,xlimit[rapbin][ptbin],mjj,yjjmax);
             printf("         in ptbin: %i  rapbin %i \n ", ptbin,rapbin);
	     exit(1);
	   }

	   // ***** identify largest/smallest pT value in each Mjj bin ******
	   //               -> to create the scale bins
	   if (ptmax<ptsmallest[rapbin][ptbin])ptsmallest[rapbin][ptbin]=ptmax;
	   if (ptmax>ptlargest[rapbin][ptbin]) ptlargest[rapbin][ptbin]=ptmax;
	   
	   // *******  identify smallest x value in each pT/y bin  *******
	   //               -> to optimize the x-limit values
	   if (xmin<xsmallest[rapbin][ptbin]*0.98){  // 0.98 reduces output
	     xsmallest[rapbin][ptbin] = xmin;
	     if((itype==amplitude_hhc::lo && nevents> 20000000)|| // 2.5h
	      (itype==amplitude_hhc::nlo && nevents> 15000000)){  // 5hh
	       //if((itype==amplitude_hhc::lo && nevents> 3000000)|| // 2.5h 
	       //(itype==amplitude_hhc::nlo && nevents> 150000)){  // 5hh
	       cout<<" "<<endl;
	       cout<<">>>>>> smaller x found in bin  "<<rapbin<<"  "<<ptbin
		   <<"   -  "<<xmin<<"  "<<xlimit[rapbin][ptbin]<<"  :  "<<
		 (xmin/xlimit[rapbin][ptbin])<<"  >> "<<nevents<<endl;
	       for( int j = 0; j < nrap; j++) {
		 for( int k = 0; k < npt[j]; k++) {
		   cout<<"    bin "<<j<<"  "<<k<<"  "<<xsmallest[j][k]
		       <<"  "<<xlimit[j][k]
		       <<"   "<<pthigh[j][k]<<"  "<<raphigh[(j+1)]<<"  :  "
		       <<(xsmallest[j][k]/xlimit[j][k])<<"  "
		       <<ptsmallest[j][k]<<"  "
		       <<ptlargest[j][k]<<"   "
		       <<ptsmallest[j][k]/ptlo[j][k]<<"  "
		       <<ptlargest[j][k]/pthi[j][k]<<" "
		       <<endl;
		 }
	       }
	     }
	   }
	   // ------------  end: identify smallest x  -------------


	   // **********  determine x_ij position in grid  ************
	   //--- determine fractional contributions for all four x-bins
	   // define the x-bin numbers in the range  [0:nxtot[
	   double hxlimit = hxlim[rapbin][ptbin];
	   int nxmin = int(nxtot *(hxmin-hxlimit)/(hxone-hxlimit));
	   int nxmax = int(nxtot *(hxmax-hxlimit)/(hxone-hxlimit));
	   
	   //-- relative distances in h(xmin), h(xmax): deltamin,deltamax 
	   double delta  = (hxone-hxlimit)/nxtot;
	   double hxi =hxlimit+double(nxmax)/double(nxtot)*(hxone-hxlimit);
	   double hxj =hxlimit+double(nxmin)/double(nxtot)*(hxone-hxlimit);
	   double deltamax = (hxmax-hxi)/delta;
	   double deltamin = (hxmin-hxj)/delta;
	   if(deltamax>1.0 || deltamin>1.0 || deltamax<0.0 || deltamin<0.0){
	     cout<<" -> deltas are off: "<<deltamax<<"  "<<deltamin<<endl;
	   }                             
	   
	   // ===== variables for the bi-cubic interpolation =====
	   // === the relative distances to the four nearest bins
	   cmax[0] = deltamax+1.0;
	   cmax[1] = deltamax;
	   cmax[2] = 1.0-deltamax;
	   cmax[3] = 2.0-deltamax;
	   cmin[0] = deltamin+1.0;
	   cmin[1] = deltamin;
	   cmin[2] = 1.0-deltamin;
	   cmin[3] = 2.0-deltamin;
	   
	   // === the weights for the cubic eigenfunctions (1-dim)
	   //   - linear interpolation in 1st and last =(nxtot-1) bins
	   //   - cubic approximation in the middle 
	   
	   if (nxmax==0 || nxmax==(nxtot-1)) { //linear in 1st and last bin
	     cefmax[0] = 0.0;
	     cefmax[1] = 1.0-cmax[1];
	     cefmax[2] = 1.0-cmax[2];
	     cefmax[3] = 0.0; }
	   else {                              // cubic in the middle
	     cefmax[1]=1.0-2.5*cmax[1]*cmax[1]+1.5*cmax[1]*cmax[1]*cmax[1];
	     cefmax[2]=1.0-2.5*cmax[2]*cmax[2]+1.5*cmax[2]*cmax[2]*cmax[2];
	     cefmax[0]=2.0 - 4.0*cmax[0] + 2.5*cmax[0]*cmax[0]
	       - 0.5*cmax[0]*cmax[0]*cmax[0];
	     cefmax[3]=2.0 - 4.0*cmax[3] + 2.5*cmax[3]*cmax[3]
	       - 0.5*cmax[3]*cmax[3]*cmax[3];
	   }
	   
	   if (nxmin==0 || nxmin==(nxtot-1)) { //linear in 1st and last bin
	     cefmin[0] = 0.0;
	     cefmin[1] = 1.0-cmin[1];
	     cefmin[2] = 1.0-cmin[2];
	     cefmin[3] = 0.0; }
	   else {                              //  cubic in the middle 
	     cefmin[1]=1.0-2.5*cmin[1]*cmin[1]+1.5*cmin[1]*cmin[1]*cmin[1];
	     cefmin[2]=1.0-2.5*cmin[2]*cmin[2]+1.5*cmin[2]*cmin[2]*cmin[2];
	     cefmin[0]=2.0 - 4.0*cmin[0] + 2.5*cmin[0]*cmin[0]
	       - 0.5*cmin[0]*cmin[0]*cmin[0];
	     cefmin[3]= 2.0 - 4.0*cmin[3] + 2.5*cmin[3]*cmin[3]
	       - 0.5*cmin[3]*cmin[3]*cmin[3];
	   }
	   
	   // === the weights for the bi-cubic eigenfunctions (2-dim)
	   for( int i1 = 0; i1 < 4; i1++) {
	     for( int i2 = 0; i2 < 4; i2++) {
	       bicef[i1][i2] = cefmax[i1] * cefmin[i2];
	     }
	   }
	   
	   // loop over scale variations for NLO
	   int scalevarmax;
	   if(itype==amplitude_hhc::lo){
	     scalevarmax=1;
	   }else{
	     scalevarmax=nscalevar;
	   }
	   amp.pdf_and_qcd_coupling(0,1.0);
	   for(int scalevar=0; scalevar<scalevarmax;scalevar++){
	     double mur2 = murscale[scalevar]*murscale[scalevar] *ptmax*ptmax;
	     double muf2 = mufscale[scalevar]*mufscale[scalevar] *ptmax*ptmax;
	     weight_hhc wt = amp(mur2,muf2);
	     wt = wt* reweight*reweight*reweight *389385.730;
	     wt = wt* unitfactor/binwidth;
	     
	     // deal with subprocesses 2 and 3
	     //    - if x1>x2 -> o.k.
	     //    - if x2>x1 -> swap weights for subprocesses 2,3
	     if(x2>x1){
	       double buffer;
	       buffer = wt[1];
	       wt[1] = wt[2];
	       wt[2] = buffer;
	     }

	     //make linear interpolation for scalebins
	     double hmu = log(log(murscale[scalevar]*ptmax/mu0scale));
	     
	     for(int scalebin=0; scalebin<nscalebin;scalebin++){
	       double deltascale;
	       switch(scalebin){
	       case 0:
		 deltascale = 1.0-(hmu-murvaltrans[rapbin][ptbin][scalevar][0])/
		   (murvaltrans[rapbin][ptbin][scalevar][1]-
		    murvaltrans[rapbin][ptbin][scalevar][0]);
		 break;
	       case 1:
		 deltascale = 1.0-(murvaltrans[rapbin][ptbin][scalevar][1]-hmu)/
		   (murvaltrans[rapbin][ptbin][scalevar][1]-
		    murvaltrans[rapbin][ptbin][scalevar][0]);
		 break;
	       default:
		 deltascale = 0.0;
		 printf("Error: not defined scalebin accessed.\n");
		 exit(1);
		 ;
	       }

	       /*
	       cout<<scalebin<<"   dsc "<<deltascale<<" "<<hmu<<" "<<
		 murvaltrans[rapbin][ptbin][scalevar][0]<<
		 murvaltrans[rapbin][ptbin][scalevar][1]<<endl;
	       */

	       // ** ----------------------------------------------------
	       // ** now fill half-table
	       // ** loop over all 16 points that receive contributions
	       for( int i1 = 0; i1 < 4; i1++) {
		 for( int i2 = 0; i2 < 4; i2++) {
		   int imax = nxmax +i1 -1;   // the target index (xmax)
		   int imin = nxmin +i2 -1;  // the target index (xmin)
		   weight_hhc wtmp = wt;  // a working copy of the weights
		   
		   // - check if above diagonal? project back and swap!
		   if (imin>imax) { 
		     int di = imin - imax;
		     imax = imax + di;   // modify indicees
		     imin = imin - di;		
		     double buffer  = wtmp[1]; // swap subprocesses 2,3
		     wtmp[1] = wtmp[2];
		     wtmp[2] = buffer;
		   } 
		   // - ignore if xmin coordinate is < 0 (contrib is zero!)
		   // - ignore if xmax coordinate is >= nxtot
		   //     (contrib. > nxtot are zero / =nxtot PDF is zero)
		   if (imax<nxtot && imin>=0) {
		     weights[rapbin][ptbin][imax][imin][scalevar][scalebin] += 
		       bicef[i1][i2] * deltascale * wtmp;
		   }
		 }
	       }
	     }//-end loop scalebin
	   }//-end loop scalevar
	 
	   // ** ------------------------------------------------------
	   // **   only in "reference mode":
	   // **         -> fill half-table with pdfs and alphas
	   // **         -> only for a single scale (index=refscale)
	   // **         -> don't throw away contributions at x=1
	   // **            but store in last existing bin No. =(nxtot-1)
	   // **                (needed for reference)
	   // **
	   if (iref==1) {  // -> in reference mode
	     
	     double mur2 = murscale[refscale]*murscale[refscale] *ptmax*ptmax;
	     double muf2 = mufscale[refscale]*mufscale[refscale] *ptmax*ptmax;
	     
	     amp.pdf_and_qcd_coupling(pdf, 389385.730);
	     weight_hhc wt = amp(mur2,muf2);
	     wt = wt * unitfactor/binwidth;
	     
	     // deal with subprocesses 2 and 3
	     //    - if x1>x2 -> o.k.
	     //    - if x2>x1 -> swap weights for subprocesses 2,3
	     if(x2>x1){
	       double buffer;
	       buffer = wt[1];
	       wt[1] = wt[2];
	       wt[2] = buffer;
	     }
	     
	     for( int i1 = 0; i1 < 4; i1++) {
	       for( int i2 = 0; i2 < 4; i2++) {
		 int imax = nxmax +i1 -1;   // the target index (xmax)
		 int imin = nxmin +i2 -1;  // the target index (xmin)
		 weight_hhc wtmp = wt;  // a working copy of the weights
		 
		 // - check if above diagonal? project back and swap!
		 if (imin>imax) { 
		   int di = imin - imax;
		   imax = imax + di;   // modify indicees
		   imin = imin - di;		
		   double buffer  = wtmp[1]; // swap subprocesses 2,3
		   wtmp[1] = wtmp[2];
		   wtmp[2] = buffer;
		 }
		 
		 // - ignore if xmin coordinate is < 0 (contrib is zero!)
		 // - for reference: rescue contrib. if coord. is >=nxtot
		 if (imin>=0) {
		   if (imax>=nxtot) imax = nxtot - 1;
		   if (imin>=nxtot) imin = nxtot - 1;
		   weights[rapbin+nrap/2][ptbin][imax][imin][0][0]+=
		     bicef[i1][i2] * wtmp;
		 }
	       }
	     }
	   
	   } //-end reference mode: a_s/PDF table filling 
	 
	 } //-end: IF in pT bins
       }  // - end: IF in rapidity bins
     }   // - ckeck pT range and rapidity range
   }    // - end: nj>=2
}

void UserHHC::writetable(){
#define WRITE(n) table.write(reinterpret_cast<char *>(&n),sizeof(n))
   //#define WRITE(n) table << setprecision(40) << n << endl

   int marker = 1234567890; //used to separate section in the table
   
   //fstream table("table05.raw",ios::out|ios::binary); // open file
   //fstream table("/work/joker-clued0/wobisch/fastNLO/table/run1/run1t99n04.raw",ios::out|ios::binary); // open file
   fstream table(tablefilename.c_str() ,ios::out|ios::binary); // open file
 
   // ireaction
   int ireaction = 2;
   WRITE(ireaction);

   unsigned int nj;
   unsigned int nu;
   unsigned int nd;
   double      s;
   inputfunc(nj,nu,nd,s);

   // Ecms
   s= sqrt(s);
   WRITE(s);

   // five strings with table content
   table << "d2sigma-dijet_dM_dy_(pb_GeV)" << endl;
   table << "CMS-LHC-Test-Scenario" << endl;
   table << "Dijet_Mass" << endl;
   table << "anti-kT_R=0.7" << endl;
   table << "-" << endl;

  //iproc
   int iproc = 2; // dijets
   WRITE(iproc);

   //ialgo
   int ialgo = 5; // anti-kT
   WRITE(ialgo);

   //JetResol1
   double JetResol1 = 0.7; // jet size R
   WRITE(JetResol1);

   //JetResol2
   double JetResol2 = 0.0; // no further parameter
   WRITE(JetResol2);

   // relative order
   int nord = itype+1; // LO: itype=0,  NLO: itype=1
   WRITE(nord);

   // absolute order
   int npow = nord+1; // -> only valid for incl and dijet x-sect
   WRITE(npow);

   switch(nord){
   case 1: 
      table << "LO" << endl;
      break;
   case 2: 
      table << "NLO"  << endl;
      break;
   default:
      table << "not known"  << endl;
   }   
   WRITE(marker);// ------------------END of block

   //nevt
   WRITE(nevents);  // No.of events in table

   //nxtot
   WRITE(nxtot);    // No of x-bins in table

   //ixscheme
   int ixscheme = 2;  //   1 log(x)   2 sqrt(log(1/x)
   WRITE(ixscheme);

   //ipdfwgt
   int ipdfwgt = 1;  //   0 no PDF weighting   1 'standard' weighting
   WRITE(ipdfwgt);

   //iref
   WRITE(iref);     // reference table included?

   WRITE(marker);// ------------------END of block

   //Nrapidity
   WRITE(nrap);
   for(int i=0;i<nrap+1;i++){
      WRITE(raphigh[i]); //Rap[0] ... Rap[Nrapidity]
   }

   for(int i=0;i<nrap;i++){
      WRITE(npt[i]); //Npt[0] ... Npt[Nrapidity-1] 
   }

   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i]+1;j++){
         WRITE(pthigh[i][j]); // all pt bins 
      }
   }
    
   WRITE(marker);// ------------------END of block
   
   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         WRITE(xlimit[i][j]); // all x limits 
      }
   }

   WRITE(marker);// ------------------END of block

   // a brief description how the scale is defined
   table << "pT_dijet_average_(GeV)" << endl;

   WRITE(nscalebin); // No. of Bins in mur,muf

   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         for(int k=0;k<nscalebin;k++){
	   WRITE(murval[i][j][k]); // all murval 
	 }
      }
   }
   
   WRITE(marker);// ------------------END of block

   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         for(int k=0;k<nscalebin;k++){
         WRITE(mufval[i][j][k]); // all mufval 
	 }
      }
   }

   WRITE(marker);// ------------------END of block

   int scalevarmax;
   if(itype==amplitude_hhc::nlo){  // scale variations are only filled above LO
      scalevarmax=nscalevar;
      WRITE(nscalevar);
      for(int i=0;i<nscalevar;i++){
         WRITE(murscale[i]); // all murscale 
      }
      for(int i=0;i<nscalevar;i++){
         WRITE(mufscale[i]); // all mufscale 
      }
      
      WRITE(marker);// ------------------END of block
   }else{
      scalevarmax=1;
   }

   for(int i=0;i<nrap;i++){ // rapidity
     for(int j=0;j<npt[i];j++){ // pt
       for(int k=0;k<nxtot;k++){ // xmax
	 for(int l=0;l<k+1;l++){ // xmin
	   // - new order in tableformat version 1c
	   for(int m=0;m<7;m++){     //subprocesses
	     // loop over scale variations for NLO
	     for(int scalevar=0; scalevar<scalevarmax;scalevar++){
	       for(int scalebin=0; scalebin<nscalebin;scalebin++){
		 WRITE(weights[i][j][k][l][scalevar][scalebin][m]);
	       }
	     }
	   }
	 }
       }
     }
   }

   WRITE(marker);// ------------------END of table

   table.close();
}

void UserHHC::end_of_event(){
   // let NLOJET++ store its results
   user_hhc::end_of_event();
   
   nevents += 1;
   //-------- store table
   if (( (unsigned long)nevents % nwrite)==0){
      time_t hour, min, time = ::time(0) - start_time;
      
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
      writetable();
      printf("done.\n");
   }
}

void UserHHC::phys_output(const std::basic_string<char>& __file_name, 
                          unsigned long __save, bool __txt) 
{
   // Suppress output of NLOJET++ files
   user_hhc::phys_output("",2000000000,false);
   tablefilename = __file_name +".raw";
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
}
