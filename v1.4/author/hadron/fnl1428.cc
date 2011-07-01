//
// fastNLO author code for fnl0005
//         LHC dijet mass distribution   pT>0.1TeV
// 
//   --- pT-hi values are often much too high!
//       >> bin2: 55%   -> bin4: 10%
//        ---> to be changed!!!!
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
#include "kt2jet-e-07.h"
//#include "kt2jet-e-07.h"
//#include "kt-e-07.h"
//#include "cone-et-07.h"
//#include "midpt-e-07.h"
//#include "rsep-et-07.h"
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
   kt2jet_e_07 jetclus;   // jet algorithm
 
   bounded_vector<lorentzvector<double> > pj;    // the jet structure 
   basic_string<char> tablefilename; // The table file to write to
   amplitude_hhc::integral_type itype; // Born, NLO etc.

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

   //  total c.m. energy square
   //s = 40000;       // RHIC           200GeV   
   //s = 3240000;     // TeV Run I     1800GeV
   //s = 3841600;       // TeV Run II    1960GeV
   s = 196000000;   // LHC          14000GeV

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
   iref = 1;       //  switch for reference mode
                  //   0: standard fastNLO table only
                 //    1: include 2nd "reference table" (a_s/PDFs)
   refscale = 1;   // which of the scalevariations is used in ref-table?

   //unitfactor = 1000000.0;  // for fb
   //unitfactor = 1000.0;  // for pb 
   unitfactor = 1.0;  // for nb    <<< use for D0 dijet mass

   //Set up binning  
   nrap   = 5;  //            // No. of bins in rapidity -> here: one
   if (iref==1) nrap=nrap*2;  // -> in reference mode: double No. rap bins

   raphigh = new double[nrap+1];  //----- array for rapidity boundaries
   npt     = new  int[nrap];      // nrap bins in rapid.-each npt[irap] pT bins

   // flexible rap binning
   raphigh[0]=0.0;   // barrel            (0.0-1.0)
   raphigh[1]=0.8;   // extended barrel   (0.8-1.7)
   raphigh[2]=1.6;   // end-cap           (1.5-2.5-3.2) fine granularity
   raphigh[3]=2.5;   //                                 coarse granularity
   raphigh[4]=3.2;   // forward           (3.1-4.9)
   raphigh[5]=4.5;

   if (iref==1)      // -> in reference mode: copy rapidity definitions
     for(int i=0;i<nrap/2;i++){
      raphigh[i+nrap/2+1] = raphigh[i+1];
   }

   //Define binning in pt -> here:  DijetMass
   for (int j=0; j<nrap; j++){
     npt[j]=15;
   }

   if (iref==1)      // -> in reference mode: copy No.pT-bin definitions
     for(int i=0;i<nrap/2;i++){
      npt[i+nrap/2] = npt[i];
   }

   // lowest pT value in sample   -> here: mass
   ptlow = 1000.;          // 

   pthigh.resize(nrap);
   //----- array for pt boundaries
   for(int i=0;i<nrap;i++){
      pthigh[i].resize(npt[i]+1);
   }
   //----- array for pt boundaries
   for (int j=0; j<nrap; j++){
     pthigh[j][0]=1000.0;
     for (int i=1; i<11; i++){    // 0.5 TeV bins -> 6TeV
       pthigh[j][i] = pthigh[j][i-1] + 500.0;
     }  
     for (int i=11; i<16; i++){    // 1TeV bins -> 10TeV
       pthigh[j][i] = pthigh[j][i-1] + 1000.0;
     }  
   }

   if (iref==1)      // -> in reference mode: copy pT-bin definitions
     for(int i=0;i<nrap/2;i++){
       for(int j=0;j<npt[i]+1;j++){
         pthigh[i+nrap/2][j] = pthigh[i][j];
       }
     }

   nxtot = 16;  // 

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
            weights[j][k][l].resize(l+1); // half matrix xmin,xax: (n^2+n)/2
            // scale variation
            for(int i = 0; i < l+1; i++){
               weights[j][k][l][i].resize(nscalevar);
	       for(int m = 0; m < nscalevar; m++){
		 weights[j][k][l][i][m].resize(nscalebin);
	       }
            }
         }

	 // --------- integrate over the accessible phase space 
         // - find xlimit and minimum and maximum pT
	 //cout<<"  --- integrate over phase space (Mjj,y1,y2)"<<endl;
         double xmin = 1.0;
	 double ptmin = 9999999.9;
	 double ptmax = 0.0;
	 int nmax = 47;
         for( int m0 = 0; m0 <=nmax; m0++) {
	   for( int m1 = 0; m1 <=nmax; m1++) {
	     for( int m2 = -nmax; m2 <= nmax; m2++) {
	       double mass = pthigh[j][k]+(pthigh[j][k+1]-pthigh[j][k])*double(m0)/double(nmax);
	       double y1 = raphigh[j]+(raphigh[j+1]-raphigh[j])*double(m1)/double(nmax);
	       double y2 = y1*double(m2)/double(nmax);
	       double pt1= mass/sqrt(2.0*(cosh(fabs(y1-y2))+1.0 ));
	       //double pt1= mass/2.0/cosh(fabs(y1-y2)/2.0);
	       double e1 = (exp(y1)+exp(y2))/(exp(-y1)+exp(-y2));
	       double x1 = mass/sqrt(s)*e1;
	       double x2 = mass/sqrt(s)/e1;
	       if (pt1 > 100.0 && x1 < 1.0 && x2 < 1.0) {
		 //cout<<"   "<<mass<<" "<<y1<<" "<<y2<<" :  "<<pt1<<endl;
		 if (x2 < xmin) xmin = x2;
		 if (x1 < xmin) xmin = x1;
		 if (pt1 < ptmin) ptmin = pt1;
		 if (pt1 > ptmax) ptmax = pt1;
	       }
	     }
	   }
	 }


	 xlimit[j][k] = xmin * 1.0 ;   // use exact result
         //pthi[j][k] = pthigh[j][k+1]/2.0;
         //ptlo[j][k] = pthigh[j][k]/sqrt(2.0*(cosh(2.0*raphigh[j+1])+1.0));
	 //ptlo[j][k] = max(ptlo[j][k], 100.0);
         pthi[j][k] = ptmax; 
	 ptlo[j][k] = ptmin;

         // ---- safety factors for E-scheme: only small safety factor
	 xlimit[j][k] = xlimit[j][k]*0.98;
	 /*
	 if (j==0 && k<8) xlimit[j][k] = xlimit[j][k] *
			    (1.0+0.56*((8.0-double(k))*(8.0-double(k)))/10.0);
	 if (j==1 && k<4) xlimit[j][k] = xlimit[j][k]*
			    (1.0+0.55*((4.0-double(k))*(4.0-double(k)))/10.0);
	 */

	 cout<<" "<<j<<" "<<k<<" - "<<xlimit[j][k]<<"  "<<
	   ptlo[j][k]<<"  "<<pthi[j][k]<<endl;

         hxlim[j][k]= -sqrt(-log10(xlimit[j][k]));
         xsmallest[j][k]=0.999;
         ptsmallest[j][k]=9999999.9;
         ptlargest[j][k]=0.0;

         // - Setup the murval and mufval arrays - for scale bins
	 double dpt = (pthi[j][k]-ptlo[j][k]);
	 double lptlo = ptlo[j][k] + dpt*0.05;   // - starting points
	 double lpthi = pthi[j][k] - dpt*0.3;
	 double llptlo = log(log(lptlo/mu0scale));     // - ln(ln((pt/mu0))
	 double llpthi = log(log(lpthi/mu0scale));

	 for(int m = 0; m < nscalebin; m++){
	   double lpt1 = mu0scale * exp(exp( llptlo +
	      double(m)/double(nscalebin-1)*(llpthi-llptlo) ));
	   murval[j][k][m] = lpt1;  // - assign values for intermed scales
	   mufval[j][k][m] = murval[j][k][m];
	 }

         for( int l = 0; l < nscalevar; l++) {    // -> only really need extreme values
	   for(int m = 0; m < nscalebin; m++){
	     murvaltrans[j][k][l][m]=log(log(murscale[l]*murval[j][k][m]/mu0scale)); 
	     mufvaltrans[j][k][l][m]=log(log(mufscale[l]*mufval[j][k][m]/mu0scale)); 
	   }
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
	 for (int m = 0; m < nscalebin; m++) {
	   cout<<"    mu-bin  "<<m<<"  "<<murval[j][k][m]<<endl;
	 }
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
   cout<<"  "<<endl;
   cout<<"   *******************************************"<<endl;
   cout<<"    fastNLO    - initialization"<<endl;
   cout<<"         *** scenario fnl0005 ***"<<endl;
   cout<<"        (ATLAS dijets - test scenario)"<<endl;
   cout<<" "<<endl;
   cout<<"        table file "<<tablefilename<<endl;
   cout<<"        store table after "<<nwrite<<" events"<<endl;
   cout<<"        sqrt(s)= "<<sqrt(s)<<endl;
   cout<<"        No. x-bins: "<<nxtot<<endl;
   cout<<"        No. rapidity regions: "<<nrap<<endl;
   cout<<"        No. of pT bins in each rapidity region:"<<endl;
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

}


void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
   //  typedef lorentzvector<double> _Lv;
   //-------- general stuff
   //  double s = 2.0*(p[hadron(0)]*p[hadron(-1)]);
   static const double twopi = 6.28318530717958647692;
   static const double pi = 3.14159265358979323846;
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


   // --- logic to identify correct jets --- (from fnt1003)
   // - find two highest pT jets
   //   - 1. compare 1,2 - order
   //   - 2. compare 2 (lower than 1!) and 3 - order
   //   - 3. compare 1,2 order
   // highest: 1,2   (only step 1.)
   // highest: 2,3   (step 1: move 1->2 , step 2: move 1 ->3)
   // highest: 1,3   (step 2: move 3->2)
   // if No.1 is central -> save rap. of No.2 with (1)
   // if No.2 is central -> save rap. of No.1 with (2) -> swap (1,2)
   // -> swapping makes sure that (new, reordered) first jet is central

   // -------------- sort jets in ET
   if (nj>=2) {
     // set pointers to the two leading pT jets
     int ij1=1;
     int ij2=2;
     double pt1=pj[1].perp();
     double pt2=pj[2].perp();
     if (pt2>pt1){
       double pttemp = pt1;
       pt1 = pt2;
       pt2 = pttemp;
       ij1 = 2;
       ij2 = 1;
     }
     if (nj==3 && pj[3].perp()>pt2 ) {
       pt2 = pj[3].perp();
       ij2 = 3;
     }
     if (pt2>pt1){        // can only happen if 3rd orig. jet has highest pT
       double pttemp = pt1;   
       pt1 = pt2;
       pt2 = pttemp;
       int tempij = ij1;
       ij1 = ij2;
       ij2 = tempij;
     }
     // define dijet variables
     // define dijet variables
     double dijmass=sqrt(( pj[ij1].T()+pj[ij2].T())*(pj[ij1].T()+pj[ij2].T())
                         -(pj[ij1].X()+pj[ij2].X())*(pj[ij1].X()+pj[ij2].X())
                         -(pj[ij1].Y()+pj[ij2].Y())*(pj[ij1].Y()+pj[ij2].Y())
                         -(pj[ij1].Z()+pj[ij2].Z())*(pj[ij1].Z()+pj[ij2].Z()));
     //double dijchi = exp( fabs( pj[ij1].rapidity() - pj[ij2].rapidity() ));
     //double yboost = fabs(pj[ij1].rapidity() + pj[ij2].rapidity() )/2.0;
     double rapmax = max(fabs(pj[ij1].rapidity()),fabs(pj[ij2].rapidity()));

     if (pt1>100.0 && dijmass > ptlow ) {
       // --- later this variable will be the ren/fact. scale
       //     double ptmax = pt1;    // either pTmax
       double ptmax = (pt1+pt2)/2.0; // or the average dijet pT

       // --- determine eta and pt(=Mass) bin 
       double binwidth = 1.0;    // 
       int rapbin = -1;
       for( int j = 0; j < nrap; j++) {
	 if (rapmax >= raphigh[j] && rapmax < raphigh[j+1]) {
	   rapbin = j;
	   binwidth = (raphigh[j+1]-raphigh[j])*2.0;
	   break;
	 }
       }
       if (rapbin >=0 ) {              
	 int ptbin  = -1;
	 for( int j = 0; j < npt[rapbin]; j++) {
	   if (dijmass>=pthigh[rapbin][j] && dijmass<pthigh[rapbin][(j+1)]) {
	     ptbin=j;
	     binwidth = binwidth*(pthigh[rapbin][(j+1)]-pthigh[rapbin][j]);
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
		    xmin,xlimit[rapbin][ptbin],dijmass,0.0);
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
	     //if((itype==amplitude_hhc::lo && nevents> 600000000)|| // 7h 
	     //(itype==amplitude_hhc::nlo && nevents> 120000000)){  // 10h
	     if((itype==amplitude_hhc::lo && nevents> 30000000)|| // 2.5h 
		(itype==amplitude_hhc::nlo && nevents> 3000000)){  // 5hh
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
	     
	     // make linear interpolation for scalebins
	     //  now: --> extend to full flexible Bernstein Polynomials
	     double hmu = log(log(murscale[scalevar]*ptmax/mu0scale));

	     // new: define position inside scale range	     <<< for Bernstein
	     double tscale = (hmu-murvaltrans[rapbin][ptbin][scalevar][0]) /
	       (murvaltrans[rapbin][ptbin][scalevar][nscalebin-1] -
		murvaltrans[rapbin][ptbin][scalevar][0]);

	     // determine contributions to different scalebins
	     //  ... = bernst0(z), bernst1(z), ..., bernstn(z)    
	     for(int scalebin=0; scalebin<nscalebin;scalebin++){
	       double deltascale = 0.0;
	       if (nscalebin == 2) {              // n=1
		 switch(scalebin){
		 case 0:
		   deltascale = 1.0-tscale;
		   break;
		 case 1:
		   deltascale = tscale;
		   break;
		 default:
		   deltascale = 0.0;
		   printf("Error: not defined scalebin accessed.\n");
		   exit(1);
		   ;
		 }
	       }
	       else if (nscalebin == 3) {              // n=2
		 switch(scalebin){
		 case 0:
		   deltascale = (1.0-tscale)*(1.0-tscale);
		   break;
		 case 1:
		   deltascale = 2.0*(1.0-tscale)*tscale;
		   break;
		 case 2:
		   deltascale = tscale*tscale;
		   break;
		 default:
		   deltascale = 0.0;
		   printf("Error: not defined scalebin accessed.\n");
		   exit(1);
		   ;
		 }
	       }
	       else if (nscalebin == 4) {              // n=3
		 switch(scalebin){
		 case 0:
		   deltascale = (1.0-tscale)*(1.0-tscale)*(1.0-tscale);
		   break;
		 case 1:
		   deltascale = 3.0*(1.0-tscale)*(1.0-tscale)*tscale;
		   break;
		 case 2:
		   deltascale = 3.0*(1.0-tscale)*tscale*tscale;
		   break;
		 case 3:
		   deltascale = tscale*tscale*tscale;
		   break;
		 default:
		   deltascale = 0.0;
		   printf("Error: not defined scalebin accessed.\n");
		   exit(1);
		   ;
		 }
	       }
	       else if (nscalebin == 5) {              // n=4
		 switch(scalebin){
		 case 0:
		   deltascale = (1.0-tscale)*(1.0-tscale)*(1.0-tscale)*(1.0-tscale);
		   break;
		 case 1:
		   deltascale = 4.0*(1.0-tscale)*(1.0-tscale)*(1.0-tscale)*tscale;
		   break;
		 case 2:
		   deltascale = 6.0*(1.0-tscale)*(1.0-tscale)*tscale*tscale;
		   break;
		 case 3:
		   deltascale = 4.0*(1.0-tscale)*tscale*tscale*tscale;
		   break;
		 case 4:
		   deltascale = tscale*tscale*tscale*tscale;
		   break;
		 default:
		   deltascale = 0.0;
		   printf("Error: not defined scalebin accessed.\n");
		   exit(1);
		   ;
		 }
	       }
	       else if (nscalebin > 5) {              // n>4
		 deltascale = 0.0;
		 printf("Error: nscalebin >4 not implemented.\n");
		 exit(1);
	       }
	     

	       
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
   int ireaction = 3;
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
   table << "d2sigma_dMass_dymax_(nb_GeV)" << endl;
   table << "(unpublished)" << endl;
   table << "ATLAS_test_scenario" << endl;
   table << "-" << endl;
   table << "-" << endl;

  //iproc
   int iproc = 2; // dijets
   WRITE(iproc);

   //ialgo
   int ialgo = 1; // kT/Aachen algo
   WRITE(ialgo);

   //JetResol1
   double JetResol1 = 0.7;  // kT/Aachen: D
   WRITE(JetResol1);

   //JetResol2
   double JetResol2 = 0.0; //   kT/Aachen no 2nd parameter

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
   table << "average_ET_of_dijets" << endl;

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
      printf ("No. events: %.0f writing table\n",nevents);
      writetable();
   }
}

void UserHHC::phys_output(const std::basic_string<char>& __file_name, 
                          unsigned long __save, bool __txt) 
{
   // First let NLOJET do the work
   user_hhc::phys_output(__file_name,__save,__txt);
   tablefilename = __file_name +".raw";
   nwrite = __save;
}
