//
// fastNLO author code for fnt1009
//         Run I   D0 dijet Chi distribution  hep-ex/0012046
//
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
//#include "kt2jet-e-07.h"
//#include "kt2jet-e-07.h"
//#include "kt-e-07.h"
#include "cone-et-07.h"
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
   cone_et_07 jetclus;   // jet algorithm
 
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
   s = 3240000;     // TeV Run I     1800GeV
   //s = 3841600;       // TeV Run II    1960GeV
   //s = 196000000;   // LHC          14000GeV

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
   //unitfactor = 1000.0;  // for pb
   unitfactor = 1.0;  // for nb     <<< use for Chi-distrib

   //Set up binning  
   nrap = 4;  //             // No. of bins in rapidity -> here: Mass
   if (iref==1) nrap=nrap*2;  // -> in reference mode: double No. rap bins

   raphigh = new double[nrap+1];  //----- array for rapidity boundaries
   npt     = new  int[nrap];      // nrap bins in rapid.-each npt[irap] pT bins

   // flexible rap binning
   raphigh[0]=260; 
   raphigh[1]=425; 
   raphigh[2]=475; 
   raphigh[3]=635; 
   raphigh[4]=1400; 

   if (iref==1)      // -> in reference mode: copy rapidity definitions
     for(int i=0;i<nrap/2;i++){
      raphigh[i+nrap/2+1] = raphigh[i+1];
   }

   //Define binning in pt -> here:  Chi bins  1-20: 19bins
   npt[0]=19;
   npt[1]=19;
   npt[2]=13;
   npt[3]=11;

   if (iref==1)      // -> in reference mode: copy No.pT-bin definitions
     for(int i=0;i<nrap/2;i++){
      npt[i+nrap/2] = npt[i];
   }

   // lowest pT value in sample   -> here: ???
   ptlow = 55.;          // absolute lowest pT-cut

   pthigh.resize(nrap);
   //----- array for pt boundaries
   for(int i=0;i<nrap;i++){
      pthigh[i].resize(npt[i]+1);
   }
   //----- array for pt boundaries (chi boundaries)
   for(int i=0;i<nrap;i++){
     pthigh[i][0]=1.0;
     for(int j=0;j<npt[i];j++){
       pthigh[i][j+1]=pthigh[i][j]+1.0;
       cout<<i<<"  "<<j<<"  "<<pthigh[i][j]<<" "<<raphigh[i]<<" "<<
	 raphigh[i+1]<<endl;
     }
   }

   if (iref==1)      // -> in reference mode: copy pT-bin definitions
     for(int i=0;i<nrap/2;i++){
       for(int j=0;j<npt[i]+1;j++){
         pthigh[i+nrap/2][j] = pthigh[i][j];
       }
     }

   nxtot = 12;

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
	 
         // - Setup the xlimit array - computed from kinematic constraints
	 //  ----> need good formula for dijets
	 // - integrate over whole dijet phase space from pT: 10-980 all y
	 // - compute dijet mass & chi 
         // -----------> if in bin -> keep track of pthi/lo

	 double xmin = raphigh[j]*raphigh[j] / s; 
	 xlimit[j][k] = xmin * 1.0;   // exact limit value!!!
	 // define reasonable upper and lower pT boundaries for each bin
	 double xhi = 1.0/sqrt(2.0*(cosh(log(pthigh[j][k]))+1.0));
	 double xlo = 1.0/sqrt(2.0*(cosh(log(pthigh[j][k+1]))+1.0));
         pthi[j][k]  = xhi * raphigh[j+1];
         ptlo[j][k]  = xlo * raphigh[j];

         // ---- safety factors for ET-scheme: bigger safety factors
	 xlimit[j][k] = xlimit[j][k]*0.98;
	 //	 xlimit[j][k] = xlimit[j][k]*0.2;
       
         hxlim[j][k]= -sqrt(-log10(xlimit[j][k]));
         xsmallest[j][k]=0.999;
         ptsmallest[j][k]=999.9;
         ptlargest[j][k]=0.0;

         // - Setup the murval and mufval arrays
	 double dpt = (pthi[j][k]-ptlo[j][k]);
	 cout<<"pT  "<<pthi[j][k]<<"  "<<ptlo[j][k]<<" "<<dpt<<"  "<<
	   (ptlo[j][k]+dpt*0.19)<<" "<<(pthi[j][k]-dpt*0.4)
	     <<endl;

         murval[j][k][0] = ptlo[j][k] +dpt*0.18;   // was .19 
         murval[j][k][1] = pthi[j][k] -dpt*0.4;
         
         if (j==3) { // special definition in last - large - mass bin
	   murval[j][k][0] = ptlo[j][k] +dpt*0.08;     //was .14  .48
	   murval[j][k][1] = pthi[j][k] -dpt*0.65;   // <<- 0.5 was checked!!
	 }
         

         mufval[j][k][0] = murval[j][k][0];
         mufval[j][k][1] = murval[j][k][1];

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
   cout<<"  "<<endl;
   cout<<"   *******************************************"<<endl;
   cout<<"    fastNLO    - initialization"<<endl;
   cout<<"         *** scenario fnt1009 ***"<<endl;
   cout<<"             D0 Run I dijet Chi "<<endl;
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
   //double s = 2.0*(p[hadron(0)]*p[hadron(-1)]);

   //---------------------------------------------------
   //--------- start event processing ------------------
   //---------------------------------------------------

   static const double pi = 3.14159265358979323846;
   static const double twopi = 6.28318530717958647692;

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
     int ij3=3;
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
       ij3 = 2;
     }
     if (pt2>pt1){        // can only happen if 3rd orig. jet has highest pT
       double pttemp = pt1;   
       pt1 = pt2;
       pt2 = pttemp;
       int tempij = ij1;
       ij1 = ij2;
       ij2 = tempij;
     }

     /*    -----  test
     if (pj[ij2].perp() > pj[ij1].perp()) cout<<" 2 > 1 ----------"<<endl;
     if (nj==3){
       if (pj[ij3].perp() > pj[ij1].perp()) cout<<" 3 > 1"<<endl;
       if (pj[ij3].perp() > pj[ij2].perp()) cout<<" 3 > 2"<<endl;
     }
     */

     // define dijet variables
     double delph = abs(pj[ij1].phi() - pj[ij2].phi() );
     if(delph > pi) delph = twopi - delph;
     double dijmass = sqrt(2.0*pj[ij1].perp()*pj[ij2].perp()*(
                      cosh(fabs(pj[ij1].prapidity() - pj[ij2].prapidity() ))
			   -cos(delph) ));
     double dijchi = exp( fabs( pj[ij1].prapidity() - pj[ij2].prapidity() ));
     double yboost = fabs(pj[ij1].prapidity() + pj[ij2].prapidity() )/2.0;

     /*     
     cout<<dijmass<<" "<<dijchi<<" "<<yboost<<"   "<<pj[ij1].rapidity()<<
       " "<<pj[ij2].rapidity()<<"   "<<
       pj[ij1].perp()<<" "<<pj[ij2].perp()<<" "<<xmax<<" "<<xmin<<endl;
     */

     // --- cuts
     if (pt1>ptlow && dijchi>1.0 && yboost<1.5 && dijmass > 260.0) {
       // --- later this variable will be the ren/fact. scale
       //     double ptmax = pt1;    // either pTmax
       double ptmax = (pt1+pt2)/2.0; // or the average dijet pT

       // cout<<ptmax<<"  "<<(pj[ij1].perp()+pj[ij2].perp())*0.5<<endl;

       if (ptmax<2.0) {
	 cout<<"aa "<<pt1<<" "<<pt2<<" "<<dijmass<<" "<<dijchi<<endl;
       }
       
       // --- determine y(=Mass) and pt(=Chi) bin 
       double binwidth = 1.0;    // 
       int rapbin = -1;
       int nloop = nrap;   // - important for iref==1: different No loops
       if (iref==1) nloop = nrap/2; // -> in reference mode
       for( int j = 0; j < nloop; j++) {        
	 if (dijmass >= raphigh[j] && dijmass < raphigh[(j+1)]) {
	   rapbin=j;
	   break;
	 }
       }
       if (rapbin >=0 ) {              
	 int ptbin  = -1;
	 for( int j = 0; j < npt[rapbin]; j++) {
	   if (dijchi >= pthigh[rapbin][j] && dijchi < pthigh[rapbin][(j+1)]) {
	     ptbin=j;
	     binwidth = pthigh[rapbin][(j+1)]-pthigh[rapbin][j];
	     break;
	   }
	 }
	 
	 
	 //---------- compute weight, fill fastNLO array
	 if ( ptbin>=0 ) {
	   // test if x_min in event is lower than x_limit
	   //      if yes -> make big warning!!! 
	   //      -> need to change x_limit values
	   if (xmin<xlimit[rapbin][ptbin]){
	     printf("Warning: xmin (%f) < xlimit (%f) at chi=%f M=%f GeV \n ",
		    xmin,xlimit[rapbin][ptbin],dijchi,dijmass);
	     double s = 2.0*(p[hadron(0)]*p[hadron(-1)]);
	     cout<<rapbin<<" "<<ptbin<<"  "<<
	       (raphigh[rapbin]*raphigh[rapbin]/s)<<endl;
	     cout<<sqrt(( pj[ij1].T()+pj[ij2].T())*(pj[ij1].T()+pj[ij2].T())
                         -(pj[ij1].X()+pj[ij2].X())*(pj[ij1].X()+pj[ij2].X())
                         -(pj[ij1].Y()+pj[ij2].Y())*(pj[ij1].Y()+pj[ij2].Y())
			-(pj[ij1].Z()+pj[ij2].Z())*(pj[ij1].Z()+pj[ij2].Z()))
		 <<" "<<fabs(pj[ij1].prapidity() - pj[ij2].prapidity())<<" "
		 <<(sqrt(2.0*pj[ij1].perp()*pj[ij2].perp()*(cosh(log(dijchi))+1.0)))
		 <<" "
		 <<(sqrt(2.0*pj[ij1].perp()*pj[ij2].perp()*(
							    cosh(fabs(pj[ij1].prapidity() - pj[ij2].prapidity()))
			   -cos(delph))))<<" "
		 <<endl;
	     // ---------- why do the two mass-formulas disagree????

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
	     if((itype==amplitude_hhc::lo && nevents> 700000000)|| // 2.5h 
		(itype==amplitude_hhc::nlo && nevents> 90000000)){  // 5hh
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
   table << "dsigma-dijet_dchi_(nb)" << endl;
   table << "hep-ex/0012046" << endl;
   table << "D0_Collaboration" << endl;
   table << "-" << endl;
   table << "-" << endl;

  //iproc
   int iproc = 2; // dijets
   WRITE(iproc);

   //ialgo
   int ialgo = 2; // midpoint cone algo
   WRITE(ialgo);

   //JetResol1
   double JetResol1 = 0.7;  // midpoint cone: R_cone
   WRITE(JetResol1);

   //JetResol2
   double JetResol2 = 0.5; //midpoint cone: f_overlap (not effective up to NLO)

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
   table << "average_pT_of_dijets" << endl;

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

