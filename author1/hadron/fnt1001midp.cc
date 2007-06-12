//
// fastNLO author code for fnt1001 for midpoint ago (no Rsep)
//
// last modification
// 2005/12/07 MW - new xlimit calculation (y-integration)

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
   
   int nrap;       // No of rapidity bins 
   double *raphigh;  // array for rapidity boundaries
   int *npt;       // No of pT bins in each y range
   vector< vector<double> >pthigh;   // array for pT boundaries
   double ptlow;    // lowest pt considered 
   
   int nscalevar;             // number of scale variations (mu_r,mu_f) for NLO
   vector <double> murscale;  // overall scale factor for renormalization scale
   vector< vector<double> >murval; // array for renormalization scale values
   vector <double> mufscale;       // overall scale factor for fact. scale
   vector< vector<double> >mufval; // array for factorization scale values

   int nxtot;      // no of xbins 
   vector< vector<double> >xlimit; // array for lower x limits
   vector< vector<double> >hxlim; // array for function h at xlimit
   vector< vector<double> >xsmallest; // array for smallest actual x values
   //    array for the weights MAIN ARRAY
   vector <vector< vector < vector < vector <weight_hhc> > > > >weights; 
   

   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table
   double xsectsum; // total cross section - internal counter for test output

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
   s = 3240000;       // TeV Run I     1800GeV
   //s = 3841600;     // TeV Run II    1960GeV
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
  

   nrap   = 6;                    // no of bins in rapidity (5 for D0 +1 for CDF)
   raphigh = new double[nrap+1];  //----- array for rapidity boundaries
   npt     = new  int[nrap];      // nrap bins in rapidity - each npt[irap] pT bins

   //Define binning in rapidity
   for(int i=0;i<nrap+1;i++){
     raphigh[i]=0.0 + i*(2.5/nrap); // use for equidistant rap binning
   }
   // flexible rap binning
   raphigh[0]=0.0;        // bins are:  0.0 - 0.5 - 1.0 - 1.5 - 2.0 - 3.0
   raphigh[1]=0.5;        // 
   raphigh[2]=1.0;        // 
   raphigh[3]=1.5;        // 
   raphigh[4]=2.0;        // 
   raphigh[5]=3.0;        // bins are:  0.0 - 0.5 - 1.0 - 1.5 - 2.0 - 3.0
   raphigh[6]=0.7;        // specific for CDF (0.1-0.7)

   //Define binning in pt (90bins)
   npt[0]=24;
   npt[1]=24;
   npt[2]=19;
   npt[3]=15;
   npt[4]=8;
   npt[5]=33;   // CDF

   // lowest pT value in sample
   ptlow = 40.5;          // 60.0 for D0 - 40.5 for CDF

   pthigh.resize(nrap);
   //----- array for pt boundaries
   for(int i=0;i<nrap;i++){
      pthigh[i].resize(npt[i]+1);
   }

   //----- array for pt boundaries
   pthigh[0][0] = 60.0;
   pthigh[0][1] = 70;
   pthigh[0][2] = 80;
   pthigh[0][3] = 90;
   pthigh[0][4] =100;
   pthigh[0][5] =110;
   pthigh[0][6] =120;
   pthigh[0][7] =130;
   pthigh[0][8] =140;
   pthigh[0][9] =150;
   pthigh[0][10]=160;
   pthigh[0][11]=170;
   pthigh[0][12]=180;
   pthigh[0][13]=190;
   pthigh[0][14]=200;
   pthigh[0][15]=210;
   pthigh[0][16]=220;
   pthigh[0][17]=230;
   pthigh[0][18]=250;
   pthigh[0][19]=270;
   pthigh[0][20]=290;
   pthigh[0][21]=320;
   pthigh[0][22]=350;
   pthigh[0][23]=410;
   pthigh[0][24]=560.0;

   pthigh[1][0] = 60.0;
   pthigh[1][1] = 70;
   pthigh[1][2] = 80;
   pthigh[1][3] = 90;
   pthigh[1][4] =100;
   pthigh[1][5] =110;
   pthigh[1][6] =120;
   pthigh[1][7] =130;
   pthigh[1][8] =140;
   pthigh[1][9] =150;
   pthigh[1][10]=160;
   pthigh[1][11]=170;
   pthigh[1][12]=180;
   pthigh[1][13]=190;
   pthigh[1][14]=200;
   pthigh[1][15]=210;
   pthigh[1][16]=220;
   pthigh[1][17]=235;
   pthigh[1][18]=250;
   pthigh[1][19]=270;
   pthigh[1][20]=290;
   pthigh[1][21]=320;
   pthigh[1][22]=350;
   pthigh[1][23]=400;
   pthigh[1][24]=530.0;

   pthigh[2][0] = 60.0;
   pthigh[2][1] = 70;
   pthigh[2][2] = 80;
   pthigh[2][3] = 90;
   pthigh[2][4] =100;
   pthigh[2][5] =110;
   pthigh[2][6] =120;
   pthigh[2][7] =130;
   pthigh[2][8] =140;
   pthigh[2][9] =150;
   pthigh[2][10]=160;
   pthigh[2][11]=170;
   pthigh[2][12]=180;
   pthigh[2][13]=190;
   pthigh[2][14]=200;
   pthigh[2][15]=220;
   pthigh[2][16]=250;
   pthigh[2][17]=290;
   pthigh[2][18]=330;
   pthigh[2][19]=460.0;

   pthigh[3][0] = 60.0;
   pthigh[3][1] = 70.0;
   pthigh[3][2] = 80.0;
   pthigh[3][3] = 90.0;
   pthigh[3][4] =100.0;
   pthigh[3][5] =110;
   pthigh[3][6] =120;
   pthigh[3][7] =130;
   pthigh[3][8] =140;
   pthigh[3][9] =150;
   pthigh[3][10]=160;
   pthigh[3][11]=170;
   pthigh[3][12]=180;
   pthigh[3][13]=200;
   pthigh[3][14]=230;
   pthigh[3][15]=320.0;

   pthigh[4][0] = 60.0;
   pthigh[4][1] = 70.0;
   pthigh[4][2] = 80.0;
   pthigh[4][3] = 90.0;
   pthigh[4][4] =100.0;
   pthigh[4][5] =110.0;
   pthigh[4][6] =130.0;
   pthigh[4][7] =160.0;
   pthigh[4][8] =210.0;

   pthigh[5][0] =  40.5;        // CDF bins from CDF ansatz
   pthigh[5][1] =  46.5;
   pthigh[5][2] =  52.44;
   pthigh[5][3] =  58.26;
   pthigh[5][4] =  64.01;
   pthigh[5][5] =  69.63;
   pthigh[5][6] =  75.19;
   pthigh[5][7] =  80.82;
   pthigh[5][8] =  86.37;
   pthigh[5][9] =  91.8;
   pthigh[5][10] =  97.37;
   pthigh[5][11] =  102.78;
   pthigh[5][12] =  108.38;
   pthigh[5][13] =  113.55;
   pthigh[5][14] =  119.2;
   pthigh[5][15] =  124.31;
   pthigh[5][16] =  130.03;
   pthigh[5][17] =  135.07;
   pthigh[5][18] =  140.86;
   pthigh[5][19] =  150.96;
   pthigh[5][20] =  162.35;
   pthigh[5][21] =  172.42;
   pthigh[5][22] =  183.84;
   pthigh[5][23] =  193.9;
   pthigh[5][24] =  205.54;
   pthigh[5][25] =  215.14;
   pthigh[5][26] =  237.13;
   pthigh[5][27] =  258.36;
   pthigh[5][28] =  280.59;
   pthigh[5][29] =  301.56;
   pthigh[5][30] =  323.9;
   pthigh[5][31] =  344.32;
   pthigh[5][32] =  383.76;
   pthigh[5][33] =  452.71;
   
   nxtot = 25;

   // NLO scale variations - for scales in GeV (factors are not squared!!!)
   nscalevar = 4;
   murscale.resize(nscalevar);
   mufscale.resize(nscalevar);
   // mur factors correspond to default scale - variation is done later in usercode
   murscale[0] = 0.5;  mufscale[0] = 0.5;
   murscale[1] = 0.5;  mufscale[1] = 0.25;
   murscale[2] = 0.5;  mufscale[2] = 1.0;
   murscale[3] = 0.5;  mufscale[3] = 2.0;

   xlimit.resize (nrap);
   hxlim.resize (nrap);
   xsmallest.resize (nrap);     // test: find smallest x values
   murval.resize (nrap);
   mufval.resize (nrap);
 
   weights.resize (nrap);
   for( int j = 0; j < nrap; j++) {
      xlimit[j].resize(npt[j]);
      hxlim[j].resize(npt[j]);
      xsmallest[j].resize(npt[j]);
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
         double pt = pthigh[j][k];
         double xt = 2*pt/sqrt(s);
         double ymax = log((1.+sqrt(1.-xt*xt))/xt);  // upper kin. y-limit
         if (ymax>raphigh[(j+1)]) ymax=raphigh[(j+1)];
	 double ymin = raphigh[(j)];
	 if (j==5) ymin=0.1;   // special case: CDF rap bin in fnt1001

	 //   find smallest x by integrating over accessible y-range
	 double xmin = 1.0; 
	 for (int nr = 0; nr <= 200; nr++) {
	   double ytest = ymin + double(nr)*(ymax-ymin)/200.0;
	   double xtest = pt*exp(-ytest)/(sqrt(s)-pt*exp(ytest));
	   if (xtest<xmin) xmin = xtest;
	 }
	 xlimit[j][k] = xmin;
	 xlimit[j][k] = xlimit[j][k]*0.75; // safety factor for ET-scheme
         
         hxlim[j][k]= -sqrt(-log10(xlimit[j][k]));
         xsmallest[j][k]=0.999;

         // - Setup the murval and mufval arrays - take value at 40% of bin
         murval[j][k]= (0.6*pt + 0.4*pthigh[j][k+1]); 
         mufval[j][k]= (0.6*pt + 0.4*pthigh[j][k+1]); 
      }
   }

   // print x-limit values at the begin of the job
   printf ("(rapidity, pt) array for this job:\n");
   printf ("#rap #pt xlimit pt_high rap_high \n");
   for( int j = 0; j < nrap; j++) {
      for( int k = 0; k < npt[j]; k++) {
         printf("%3d %3d   %8.6f %8.1f %8.1f \n",j,k,
		xlimit[j][k],pthigh[j][k],raphigh[(j+1)]);
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
   // phys(1, "incl. jets pT>100  y<0.5", 1 , 0.0, 1.0);
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

   //-------- analyze inclusive jets in jet loop -------
   for( int i = 1; i <= nj; i++) {

      //------- get jet properties
      double pt = pj[i].perp(); 
      if (pt>ptlow){ // check for low pt first, will fail most of the time
         double rap = abs(pj[i].prapidity());
        
         //------ determine y and pt bin --- for D0 w/ 5 pseudo-rapidity regions
         int rapbin = -1;
         //for( int j = 0; j < nrap; j++) {         // for standard jobs
         for( int j = 0; j < (nrap-1); j++) {       // "-1" to omit CDF bin
            if (rap >= raphigh[j] && rap < raphigh[(j+1)]) rapbin=j;
         }
         // if (rapbin >=0 ) {         // for standard jobs         
         if (rapbin >=0 && pt>60) {    // only in fnt1001: harder pT cut for D0
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
                  printf("fast01.cc: Warning: xmin (%f) < xlimit (%f) at pt=%f GeV y=%f \n ",
                         xmin,xlimit[rapbin][ptbin],pt,rap);
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

               
               // *****  test:  identify smallest x value in each pT/y bin ******
               //               -> to optimize the x-limit values
               if (xmin<xsmallest[rapbin][ptbin]*0.98){  //  0.98 reduces output!
		 xsmallest[rapbin][ptbin] = xmin;
                 if((itype==amplitude_hhc::lo && nevents> 600000000)|| // output
		  (itype==amplitude_hhc::nlo && nevents> 60000000)){  // after3h
		   cout<<" "<<endl;
		   cout<<">>>>>>> smaller x found in bin  "<<rapbin<<"  "<<ptbin
		       <<"   -  "<<xmin<<"  "<<xlimit[rapbin][ptbin]<<"  :  "<<
		     (xmin/xlimit[rapbin][ptbin])<<"  >> "<<nevents<<endl;
                   for( int j = 0; j < nrap; j++) {
                     for( int k = 0; k < npt[j]; k++) {
                       cout<<"    bin "<<j<<"  "<<k<<"  "<<xsmallest[j][k]
			   <<"  "<<xlimit[j][k]
                           <<"   "<<pthigh[j][k]<<"  "<<raphigh[(j+1)]<<"  :  "<<
                         (xsmallest[j][k]/xlimit[j][k])<<endl;
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

                  // -- compensate for missing contrib.: nxmin>nxmax (use symmetry)
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
            } //-end: IF in D0 pT range
         }  // - end: IF in D0 rapidity range


         // ---------------------------------------------------------------------
         // ---------------------------------------------------------------------
         // ------------ code for CDF jets --------------------------------------
         // ---------------------------------------------------------------------
         // ---------------------------------------------------------------------
         if (rap >= 0.1 && rap < 0.7) {
            rapbin=5;
            int ptbin  = -1;
            for( int j = 0; j < npt[rapbin]; j++) {
               if (pt >= pthigh[rapbin][j] && pt < pthigh[rapbin][(j+1)]) ptbin=j;
            }
            //---------- compute weight, fill NLOJET histos
            if ( ptbin>=0 ) {
               if (xmin<xlimit[rapbin][ptbin]){
                  printf("fast01.cc: Warning: (CDF) xmin (%f) < xlimit (%f) at pt=%f GeV y=%f \n ",
                         xmin,xlimit[rapbin][ptbin],pt,rap);
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


               // *****  test:  identify smallest x value in each pT/y bin ******
               //               -> to optimize the x-limit values
               if (xmin<xsmallest[rapbin][ptbin]*0.98){  //  0.98 reduces output!
                 xsmallest[rapbin][ptbin] = xmin;
                 if((itype==amplitude_hhc::lo && nevents> 600000000)|| // output
                  (itype==amplitude_hhc::nlo && nevents> 60000000)){  // after3h
                   cout<<" "<<endl;
                   cout<<">>>>>>> smaller x found in bin  "<<rapbin<<"  "<<ptbin
                       <<"   -  "<<xmin<<"  "<<xlimit[rapbin][ptbin]<<"  :  "<<
                     (xmin/xlimit[rapbin][ptbin])<<"  >> "<<nevents<<endl;
                   for( int j = 0; j < nrap; j++) {
                     for( int k = 0; k < npt[j]; k++) {
                       cout<<"    bin "<<j<<"  "<<k<<"  "<<xsmallest[j][k]
			   <<"  "<<xlimit[j][k]
                           <<"   "<<pthigh[j][k]<<"  "<<raphigh[(j+1)]<<"  :  "<<
                         (xsmallest[j][k]/xlimit[j][k])<<endl;
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

                  double mur2=murscale[scalevar]*murval[rapbin][ptbin]
		    *murval[rapbin][ptbin];
                  double muf2=mufscale[scalevar]*mufval[rapbin][ptbin]
		    *mufval[rapbin][ptbin];
                  weight_hhc wt = amp(mur2,muf2);

                  //  in reference jobs:  comment the following lines
                  //wt = wt/xmin/xmax * 389385.730;   // simple x * PDF
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
                  if (nxmin<nxmax) {             // should work with half-grid
                     weights[rapbin][ptbin][nxmax][nxmin+1][scalevar] += 
		       (1-deltamax)*deltamin*wt;
                  }	    
	    
                  // -- compensate for missing contrib.: nxmin>nxmax (use symmetry)
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
            }// - end: IF in CDF pT range
         }  // - end: IF in CDF rapidity range 
      }     // - end: IF ptlow-cut
   }        // - end jet loop

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
   int iproc = 1; // incl. jets
   WRITE(iproc);

   //ialgo
   int ialgo = 2; // midpoint cone algo
   WRITE(ialgo);

   //JetResol1
   double JetResol1 = 0.7;  // midpoint cone: R_cone
   WRITE(JetResol1);

   //JetResol2
   double JetResol2 = 0.5;  // midpoint cone: f_overlap (not effective up to NLO)
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
   if(itype==amplitude_hhc::nlo){   // scale variations are only filled above LO
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

