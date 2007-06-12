// Fast Computation of Jet Cross Sections in Hadron-Induced Processes
//
//    code for preliminary D0 Run II dijets       2005/11/21 MW
//       as function of mass
//       midpoint cone algo (no Rsep)
//       -> choose mur = dij-mass


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
#include "cone-e-07.h"
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
   
   int nrap;       // No of rapidity bins 
   double *raphigh;  // array for rapidity boundaries
   int *npt;       // No of pT bins in each y range
   vector< vector<double> >pthigh;   // array for pT boundaries
   double ptlow;    // lowest pt considered 
   
   int nscalevar;                 // number of scale variations (mu_r, mu_f) for NLO
   vector <double> murscale;        // overall scale factor for renormalization scale
   vector< vector<double> >murval;   // array for renormalization scale values
   vector <double> mufscale;        // overall scale factor for factorization scale
   vector< vector<double> >mufval;    // array for factorization scale values

   int nxtot;      // no of xbins 
   vector< vector<double> >xlimit; // array for lower x limits
   vector< vector<double> >hxlim; // array for function h at xlimit
   vector< vector<double> >xsmallest; // array for smallest actual x values
   vector< vector<double> >ptsmallest; // array for smallest actual pT values
   vector< vector<double> >ptlargest; // array for largest actual pT values
   vector< vector<double> >pthi; // array for largest pT prediction
   vector< vector<double> >ptlo; // array for smallest pT prediction
   //    array for the weights MAIN ARRAY
   vector <vector< vector < vector < vector <weight_hhc> > > > >weights; 
   
   
   double nevents; // no of events calculated so far
   unsigned long nwrite;  // no of events after to write out the table
   double xsectsum; // total cross section - internal counter for test output

   pdf_cteq6 pdf;  //   pdf
   cone_e_07 jetclus;   // jet algorithm
 
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
   //s = 3240000;       // TeV Run I     1800GeV
   s = 3841600;     // TeV Run II    1960GeV
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

   xsectsum = 0;

   //Set up binning
  
   nrap   = 5;                    // no of bins in rapidity
   // rap bins: 1,2,3,4    0-0.8/1.6-2.4   cc cf ff-opp ff-same
   //           5 integrated over y: 0-2.4 (only the rest, not in 1-4)
   raphigh = new double[nrap+1];  //----- array for rapidity boundaries
   npt     = new  int[nrap];      // nrap bins in rapidity - each npt[irap] pT bins

   //Define binning in rapidity
   for(int i=0;i<nrap+1;i++){
      raphigh[i]=0.0 + i*(2.5/nrap); 
   }

   // -> for new jobs: define "rap" bins in steps of 0.5 (divide by 2*rap)
   raphigh[0]=0.0;
   raphigh[1]=0.5;
   raphigh[2]=1.0;
   raphigh[3]=1.5;
   raphigh[4]=2.0;
   raphigh[5]=2.5;

   //Define binning in mass   - start at mass=100 - total bins: 575
   npt[0]=150; // tot
   npt[1]=150; // cc
   npt[2]=100; // cf
   npt[3]=150; // ff opp
   npt[4]=25; // ff same

   // lowest mass value in sample
   ptlow = 100.0;

   pthigh.resize(nrap);
   //----- array for pt boundaries
   for(int i=0;i<nrap;i++){
      pthigh[i].resize(npt[i]+1);
      pthigh[i][0] = 100.0;
   }
   for (int j=0; j<nrap; j++){
     for (int i=0; i<npt[j]; i++){
       pthigh[j][i+1] = pthigh[j][i]+10.0;
     }  
   }

   nxtot = 25;

   // NLO scale variations - for scales in GeV (factors are not squared!!!)
   nscalevar = 4;
   murscale.resize(nscalevar);
   mufscale.resize(nscalevar);
   
   // mur factors correspond to default scale - variation is done later in usercode
   murscale[0] = 0.25;  mufscale[0] = 0.125;
   murscale[1] = 0.25;  mufscale[1] = 0.25;
   murscale[2] = 0.25;  mufscale[2] = 0.5;
   murscale[3] = 0.25;  mufscale[3] = 1.0;

   xlimit.resize (nrap);
   hxlim.resize (nrap);
   xsmallest.resize (nrap);     // test: find smallest x values
   ptsmallest.resize (nrap);     // test: find smallest pT values in mass bins
   ptlargest.resize (nrap);     // test: find largest pT values in mass bins
   pthi.resize (nrap);     //
   ptlo.resize (nrap);     //
   murval.resize (nrap);
   mufval.resize (nrap);
 
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
      weights[j].resize(npt[j]);
      for( int k = 0; k < npt[j]; k++) {
         // Setup the weights array
         weights[j][k].resize(nxtot);
         for( int l = 0; l < nxtot; l++) {
            weights[j][k][l].resize(l+1); // half matrix xmin,xmax: (n^2+n)/2
            // scale variation
            for(int i = 0; i < l+1; i++){
               weights[j][k][l][i].resize(nscalevar);
            }
         }
         // - Setup the xlimit array - computed from kinematic constraints
	 double y1max = 2.4;
	 double y1min = 0.0;
	 double y2max = 2.4;
	 double y2min = -2.4;
	 if (j==0) {
	   y1min = -2.4;
	   y1max = 2.4;
	   y2min = 0.8;
	   y2max = 1.6;
	   // - second jet is in the experimental "gap" 0.8-1.6
           // ---- new with sqrt o.k. - but need 4% safety
	 }
	 else if (j==1) {
	   y1min = -0.8;
	   y1max = 0.8;
	   y2min = 0.0;
	   y2max = 0.8;
           // ------ new with sqrt o.k. - need 1%v safety
	 }
	 else if (j==2) {
	   y1min = -0.8;
	   y1max = 0.8;
	   y2min = 1.6;
	   y2max = 2.4;
	   // - limits for 15 lowest Mass bins are too low (LO kinematics)
           // - limits for mass bins > 79 are not filled ???
           //         can't be explained by LO kinematics!!!!
           //   should be 0.21 - 0.31 continuously rising
           // ----- now with sqrt almost perfect (need 10% safety)
	 }
	 else if (j==3) {
	   y1min = 1.6;
	   y1max = 2.4;
	   y2min = -2.4;
	   y2max = -1.6;
	   // limits are slightly too low - probably because pT cut not used
           // ---------- now with sqrt o.k.
           // but need solution for 0-16 (no LO phase space)
	 }
	 else if (j==4) {
	   y1min = 1.6;
	   y1max = 2.4;
	   y2min = 1.6;
	   y2max = 2.4;
	   // - something is wrong - current xlimits are always =1
           //  should be from 0.005 to 0.3  continuously rising
           // ------- now with sqrt perfect
	 }
         double pt = pthigh[j][k];     // is here really dijet-mass
	 double xmin = 1.0;
	 double ptmin = 1000.0;
	 double ptmax = 0.0;
	 int nintbin = 200;
         for (int nr1 = 0; nr1 <= nintbin; nr1++) {
           for (int nr2 = 0; nr2 <= nintbin; nr2++) {
             double y1test = y1min +(y1max-y1min)*double(nr1)/double(nintbin);
             double y2test = y2min +(y2max-y2min)*double(nr2)/double(nintbin);
	     double e1 = exp(y1test)+exp(y2test);
	     double e2 = exp(-y1test)+exp(-y2test);
	     double xtest1 = pt/1960.0 * sqrt(e1/e2);
	     double xtest2 = pt/1960.0 * sqrt(e2/e1);
	     double xtest = min(xtest1,xtest2); 
	     double xtestb= max(xtest1,xtest2);
	     double pttest = pt/sqrt( 2.0*(cosh(y1test-y2test)+1));
	     if (pttest>40.0 && xtestb<1 && xtest<xmin) xmin = xtest;
	     
	     if (pttest<ptmin && xtestb<1) ptmin = pttest;
	     if (pttest>ptmax && xtestb<1) ptmax = pttest;
           }
	 }
         xlimit[j][k] = xmin;
	 pthi[j][k]  = max(60.0,ptmax);
	 ptlo[j][k]  = max(60.0,ptmin);

	 // ----- safety margins for NLO kinematics - different ....
	 //pthi[j][k] = pthi[j][k] * 1.6;    // for mu=pTmax
	 //ptlo[j][k] = ptlo[j][k] * 0.8;    // for mu=pTmax
	 pthi[j][k] = pthi[j][k] * 1.28;   // for mu=pT-average
	 ptlo[j][k] = ptlo[j][k] * 0.81;   // for mu=pT-average
	 if (j==2) {
	   xlimit[j][k] = xlimit[j][k] * 0.82; 
	 }
	 else if (j==3) {
	   xlimit[j][k] = xlimit[j][k] * 0.7; 
	   if (k<11) xlimit[j][k] = 0.06 + double(k)*0.001; // empty P.S.
	 }
	 else {
	   xlimit[j][k] = xlimit[j][k] * 0.9; 
	 }

	 hxlim[j][k]  = -sqrt(-log10(xlimit[j][k]));
	 xsmallest[j][k]=0.999;
	 ptsmallest[j][k]=999.9;
	 ptlargest[j][k]=0.0;

	 cout<<j<<" "<<k<<"  "<<pt<<"  "<<xmin<<"    "<<ptmin<<"  "<<ptmax<<
	   "     "<<ptmax/ptmin<<endl;


         // - Setup the murval and mufval arrays
         murval[j][k]= (0.6*pt + 0.4*pthigh[j][k+1]); // Take as mean pt value
         mufval[j][k]= (0.6*pt + 0.4*pthigh[j][k+1]); // at 40% in bin
      }
   }
   cout<<"* "<<endl;
   cout<<"* -------------------------------------------"<<endl;
   cout<<"*  the final limits on x, pt_lo, pt_hi: "<<endl;
   cout<<"* "<<endl;
   // print x-limit values at the begin of the job
   printf ("(rapidity, pt) array for this job:\n");
   printf ("#rap #pt  xlimit       pt_lo   pt_hi  pthi/ptlo  pt_high rap_high \n");
   for( int j = 0; j < nrap; j++) {
      for( int k = 0; k < npt[j]; k++) {
         printf("%3d %3d   %8.6f   %8.1f %8.1f  %5.1f  %8.1f %8.1f \n",j,k,
		xlimit[j][k],ptlo[j][k],pthi[j][k], pthi[j][k]/ptlo[j][k],
		pthigh[j][k],raphigh[(j+1)]);
      }
   }

   // ---------------------------------------------------------------
   // ---- Initialize event counters
   nevents = 0;
   // Set some defaults
   if (nwrite==0) nwrite = 2500000;
   if (tablefilename=="") tablefilename = "table.raw";

   // Say Hello
   cout<<"  "<<endl;
   cout<<"   *******************************************"<<endl;
   cout<<"    fastNLO    - initialization"<<endl;
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
   cout<<"   *******************************************"<<endl;
   cout<<"        "<<endl;


   // Histogram Definitions
   // first test: single bin
   //phys(1, "incl. jets pT>100  y<0.5", 1 , 0.0, 1.0);
   // phys(1, "incl. jets pT>100  y<0.5", 75 , 50.0, 800.0);
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

   //-----  pdf and a_s -----
   //amp.pdf_and_qcd_coupling(pdf, 389385.730);  
   //pdf.mode(pdf_cteq6::nlo); pdf.loop(2);

   // ------------- sort jets in pT
   if (nj>=2) {
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
       double pttemp = pt1;    // not really needed for dijet mass (sort)
       pt1 = pt2;
       pt2 = pttemp;
       int tempij = ij1;
       ij1 = ij2;
       ij2 = tempij;
     }
     // define dijet variables
     double dijmass = sqrt((pj[ij1].T()+pj[ij2].T())*(pj[ij1].T()+pj[ij2].T())
			   -(pj[ij1].X()+pj[ij2].X())*(pj[ij1].X()+pj[ij2].X())
			   -(pj[ij1].Y()+pj[ij2].Y())*(pj[ij1].Y()+pj[ij2].Y())
			   -(pj[ij1].Z()+pj[ij2].Z())*(pj[ij1].Z()+pj[ij2].Z()));
     // --- later this variable will be the ren/fact. scale
     //     double ptmax = pt1;    // either pTmax
     double ptmax = (pt1+pt2)/2.0; // or the average dijet pT

     double rapmax = max(abs(pj[ij1].rapidity()),abs(pj[ij2].rapidity()));
     double rapmin = min(abs(pj[ij1].rapidity()),abs(pj[ij2].rapidity()));
     double sameside = pj[ij1].rapidity() * pj[ij2].rapidity();
     //cout<<dijmass<<" "<<rapmax<<" "<<rapmin<<endl;
     
     int rapbin = -1;
     // check if in rapidity range 1
     if (rapmax<2.4) rapbin=0; // in huge range
     // check if in rapidity ranges 2-5
     if (rapmax<0.8) rapbin=1;                             // cc
     if (rapmin<0.8 && rapmax>1.6 && rapmax<2.4) rapbin=2; // cf
     if (rapmin>1.6 && rapmax<2.4 && sameside<0) rapbin=3; // ff opposite
     if (rapmin>1.6 && rapmax<2.4 && sameside>0) rapbin=4; // ff same
     
     // apply selection cuts
     if (pt2>40 && pt1>60 && dijmass>100 && rapbin>=0) {
       double pt = dijmass; 
       int ptbin  = -1;
       for( int j = 0; j < npt[rapbin]; j++) {
	 if (pt >= pthigh[rapbin][j] && pt < pthigh[rapbin][(j+1)]) ptbin=j;
       }
       
       //---------- compute weight, fill NLOJET histos
       if ( ptbin>=0 ) {
	 // test if x_min in event is lower than x_limit
	 //      if yes -> make big warning!!! 
	 //      -> need to change x_limit values
	 if (xmin<xlimit[rapbin][ptbin]){
	   printf(" Warning: xmin (%f) < xlimit (%f) at pt=%f GeV y=%f \n ",
		  xmin,xlimit[rapbin][ptbin],pt,rapmax);
	   exit(1);
	 }

	 //---------- determine x_ij position in grid
	 //--- determine fractional contributions for all four x-bins
	 // define the x-bin numbers in the range  [0:nxtot[
	 double hxlimit = hxlim[rapbin][ptbin];
	 int nxmin = int(nxtot *(hxmin-hxlimit)/(hxone-hxlimit));
	 int nxmax = int(nxtot *(hxmax-hxlimit)/(hxone-hxlimit));
	 
	 //--- relative distances deltamin,deltamax (= fractional contrib.)
	 double delta  = (hxone-hxlimit)/nxtot;
	 double hxi = hxlimit + double(nxmax)/double(nxtot)*(hxone-hxlimit);
	 double hxj = hxlimit + double(nxmin)/double(nxtot)*(hxone-hxlimit);
	 double deltamax = (hxmax-hxi)/delta;
	 double deltamin = (hxmin-hxj)/delta;
	 if(deltamax>1.0  || deltamin>1 || deltamax<0.0  || deltamin<0.0){
	   cout<<"  ->>> deltas are off: "<<deltamax<<"  "<<deltamin<<endl;
	 }
	 
	 
	 // ***** identify largest/smallest pT value in each Mjj bin ******
	 //               -> to create the scale bins
	 if (ptmax<ptsmallest[rapbin][ptbin])ptsmallest[rapbin][ptbin]=ptmax;
	 if (ptmax>ptlargest[rapbin][ptbin]) ptlargest[rapbin][ptbin]=ptmax;


	 // *****  test:  identify smallest x value in each pT/y bin ******
	 //               -> to optimize the x-limit values
	 if (xmin<xsmallest[rapbin][ptbin]*0.98){  //  0.98 reduces output!
	   xsmallest[rapbin][ptbin] = xmin;
	   if((itype==amplitude_hhc::lo && nevents> 700000000)|| // output
	      (itype==amplitude_hhc::nlo && nevents> 70000000)){  // after2h
	   //if((itype==amplitude_hhc::lo && nevents> 90000000)|| // output
	   //(itype==amplitude_hhc::nlo && nevents> 9000000)){  // after13min
	     cout<<" "<<endl;
	     cout<<">>>>>>> smaller x found in bin  "<<rapbin<<"  "<<ptbin
		 <<"   -  "<<xmin<<"  "<<xlimit[rapbin][ptbin]<<"  :  "<<
	       (xmin/xlimit[rapbin][ptbin])<<"  >> "<<nevents<<endl;
	     for( int j = 0; j < nrap; j++) {
	       for( int k = 0; k < npt[j]; k++) {
		 cout<<"  bin "<<j<<" "<<k<<"  "<<xsmallest[j][k]<<" : "
		     <<(xsmallest[j][k]/xlimit[j][k])<<"  "
		   //<<"  "<<xlimit[j][k]
		   //<<"   "<<pthigh[j][k]<<"  "<<raphigh[(j+1)]
		     <<ptsmallest[j][k]<<"  "
		     <<ptlargest[j][k]<<"   "
		     <<ptsmallest[j][k]/ptlo[j][k]<<"  "
		     <<ptlargest[j][k]/pthi[j][k]<<" "
		     <<endl;
	       }
	     }
	   }
	 }
	 // ****************  end: identify smallest x  **********
	 
	 
	 
	 // loop over scale variations for NLO
	 int scalevarmax;
	 if(itype==amplitude_hhc::lo){
	   scalevarmax=1;
	 }else{
	   scalevarmax=nscalevar;
	 }
	 for(int scalevar=0; scalevar<scalevarmax;scalevar++){
	   
	   double mur2=murscale[scalevar]*murscale[scalevar]
	     *murval[rapbin][ptbin]*murval[rapbin][ptbin];
	   double muf2=mufscale[scalevar]*mufscale[scalevar]
	     *mufval[rapbin][ptbin]*mufval[rapbin][ptbin];
	   weight_hhc wt = amp(mur2,muf2);
	   
	   //physfilld(1, 0.5, wt);
	   //wt = wt/xmin/xmax * 389385.730;   // simple x * PDF
	   //  in reference jobs:  comment the following line
	   wt = wt* reweight*reweight*reweight *389385.730;
	   
	   // deal with subprocesses 2 and 3
	   //    - if x1>x2 -> o.k.
	   //    - if x2>x1 -> swap weights for subprocesses 2,3
	   if(x2>x1){
	     double buffer;
	     buffer = wt[1];
	     wt[1] = wt[2];
	     wt[2] = buffer;
	   }
	   
	   // test: total cross section in analysis phase space
	   if(scalevar==0){
	     xsectsum += wt[0]+wt[1]+wt[2]+wt[3]+wt[4]+wt[5]+wt[6];
	   }
	   
	   // ** ----------------------------------------------------------
	   // **   now fill half-table
	   weights[rapbin][ptbin][nxmax][nxmin][scalevar] += 
	     (1-deltamax)*(1-deltamin)*wt;
	   
	   // -- avoid writing behind  nxmax=(nxtot-1)
	   if (nxmax<(nxtot-1)) {
	     weights[rapbin][ptbin][nxmax+1][nxmin][scalevar] += 
	       deltamax * (1-deltamin)*wt;
	     weights[rapbin][ptbin][nxmax+1][nxmin+1][scalevar] += 
	       deltamax * deltamin*wt;
	   }
	   
	   // -- avoid writing behind  nxmin=nxmax
	   if (nxmin<nxmax) {             // should work with half-table
	     weights[rapbin][ptbin][nxmax][nxmin+1][scalevar] += 
	       (1-deltamax)*deltamin*wt;
	   }
	   
	   // -- compensate for missing contrib. nxmin>nxmax (use symmetry)
	   //     ---> project  (nxmax,nxmin+1) back to (nxmax+1,nxmin)
	   if (nxmin==nxmax && nxmin<(nxtot-1)) { 
	     double buffer;
	     buffer = wt[1];
	     wt[1] = wt[2];
	     wt[2] = buffer; 
	     weights[rapbin][ptbin][nxmax+1][nxmin][scalevar] += 
	       (1-deltamax)*deltamin*wt;
	   }
	 }//-end loop scalevar
       } //- end: IF in pT range
     }  // - end: IF in selection range
   }   // -  end: IF two jets


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
   double JetResol2 = 0.5; // midpoint cone: f_overlap (not effective up to NLO)
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
   WRITE(nevents); // No.of events in table

   //nxtot
   WRITE(nxtot);   // No of x-bins in table

   //ixscheme
   int ixscheme = 2; //   1 log(x)   2 sqrt(log(1/x)
   WRITE(ixscheme);

   //ipdfwgt
   int ipdfwgt = 1;//   0 no PDF weighting   1 'standard' weighting
   WRITE(ipdfwgt);

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


   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         WRITE(murval[i][j]); // all murval 
      }
   }
   
   WRITE(marker);// ------------------END of block

   
   for(int i=0;i<nrap;i++){
      for(int j=0;j<npt[i];j++){
         WRITE(mufval[i][j]); // all mufval 
      }
   }

   WRITE(marker);// ------------------END of block


   int scalevarmax;
   if(itype==amplitude_hhc::nlo){
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
               // loop over scale variations for NLO
               for(int scalevar=0; scalevar<scalevarmax;scalevar++){
                  for(int m=0;m<7;m++){     //subprocesses
                     WRITE(weights[i][j][k][l][scalevar][m]);
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
      cout<<"    total x-section so far:  "<<xsectsum/nevents<<endl;
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

