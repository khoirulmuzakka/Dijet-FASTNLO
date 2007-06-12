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
#include <time.h>

#include <xmlfastnlo.h>

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

   XMLFastNlo *xml; // Interface to the xml representation
   bool textoutput; // If true, the table is written in plain ASCII instead of BASE64 encoded doubles

   bool nlo;       // Is the job running at LO or NLO?
   // binning
   
   int nrap;       // No of rapidity bins 
   double *raphigh;  // array for rapidity boundaries
   int *npt;       // No of pT bins in each y range
   vector< vector<double> >pthigh;   // array for pT boundaries
   double ptlow;    // lowest pt considered 
   

   int nscalebins; // no of bins in the scale
   vector <vector< vector<double> > >murval; // array for renormalization scale values
   vector <vector < vector<double> > >mufval; // array for factorization scale values

   int nscalevar;             // number of scale variations (mu_r,mu_f) for NLO
   vector <double> murscale;  // overall scale factor for renormalization scale
   vector <double> mufscale;       // overall scale factor for fact. scale

   int nxtot;      // no of xbins 
   vector< vector<double> >xlimit; // array for lower x limits
   vector< vector<double> >hxlim; // array for function h at xlimit
   vector< vector<double> >xsmallest; // array for smallest actual x values



   //    array for the weights MAIN ARRAY
   vector < vector < vector < vector < vector < vector < weight_hhc > > > > > >weights; 
   

   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table

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

   //Set up binning

   nrap   = 5;                    // no of bins in rapidity
   raphigh = new double[nrap+1];  //----- array for rapidity boundaries
   npt     = new  int[nrap];      // nrap bins in rapidity - each npt[irap] pT bins

   // flexible rap binning
   raphigh[0]=0.0;        // bins are:  0.0 - 0.5 - 1.0 - 1.5 - 2.0 - 3.0
   raphigh[1]=0.5;        // 
   raphigh[2]=1.0;        // 
   raphigh[3]=1.5;        // 
   raphigh[4]=2.0;        // 
   raphigh[5]=3.0;        // bins are:  0.0 - 0.5 - 1.0 - 1.5 - 2.0 - 3.0

   //Define binning in pt (90bins)
   npt[0]=24;
   npt[1]=24;
   npt[2]=19;
   npt[3]=15;
   npt[4]=8;

   // lowest pT value in sample
   ptlow = 60.;         

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
   
   nxtot = 25;
   
   nscalebins = 2;

   // NLO scale variations - for scales in GeV (factors are not squared!!!)
   // mur factors correspond to default scale - variation is done later in usercode
   if(nlo){
      nscalevar = 4;
      murscale.resize(nscalevar);
      mufscale.resize(nscalevar);
      murscale[0] = 0.5;  mufscale[0] = 0.5;
      murscale[1] = 0.5;  mufscale[1] = 0.25;
      murscale[2] = 0.5;  mufscale[2] = 1.0;
      murscale[3] = 0.5;  mufscale[3] = 2.0;
   }else{
      nscalevar = 1;
      murscale.resize(nscalevar);
      mufscale.resize(nscalevar);
      murscale[0] = 0.5;  mufscale[0] = 0.5;
      
   }

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
         //Set up the scale bins
         murval[j][k].resize(nscalebins);
         mufval[j][k].resize(nscalebins);

         // Setup the weights array
         weights[j][k].resize(nxtot);
         for( int l = 0; l < nxtot; l++) {
            weights[j][k][l].resize(l+1); // half matrix xmin,xmax: (n^2+n)/2
            // scale variation
            for(int i = 0; i < l+1; i++){
               weights[j][k][l][i].resize(nscalevar);
               for( int m = 0; m < nscalevar; m++) {
                  weights[j][k][l][i][m].resize(nscalebins);
               }               
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
//          murval[j][k][0]= (0.6*pt + 0.4*pthigh[j][k+1]); 
//          mufval[j][k][0]= (0.6*pt + 0.4*pthigh[j][k+1]); 
         
         // The following could be used for nscalebins=2, right?
         murval[j][k][0]= pt; 
         mufval[j][k][0]= pt; 
         murval[j][k][1]= pthigh[j][k+1]; 
         mufval[j][k][1]= pthigh[j][k+1]; 


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
   cout<<"        No. of scale variations : "<<nscalevar<<endl;
   for( int j = 0; j < nscalevar; j++) {
     cout<<"          "<<j<<":   (mur/pT) "<<murscale[j]<<
	"      (muf/pT) "<<mufscale[j]<<endl;
   }
   cout<<"   *******************************************"<<endl;
   cout<<"        "<<endl;

   xml = new XMLFastNlo();
   xml->AddElementwithValue("VersionNo","1");
   //-------------------- Date
   xml->PushElement("Date");
   time_t now = time(0);
   xml->AddElementwithValue("Value",ctime(&now));
   xml->AddElementwithValue("Format","day month date hours:minutes:seconds year");
   xml->PopElement();
   //-------------------- Description 
   xml->PushElement("Description");
   xml->AddElementwithValue("Label","fnt1001midp");
   xml->AddElementwithValue("hepex","?");    

   xml->PushElement("Comment","lines",2); 
   xml->AddElementwithValue("Label","incl. jets");    
   xml->AddElementwithValue("Label","w/o Rsep");    
   xml->PopElement();

   xml->PopElement();
   //--------------------- Definition
   xml->PushElement("Definition");

   xml->PushElement("Reaction","ID",3); 
   xml->AddElementwithValue("Label","ppbar");    
   xml->AddElementwithValue("CenterMassEnergy",sqrt(s));    
   xml->PopElement();

   xml->PushElement("SubReaction","ID",1); 
   xml->AddElementwithValue("nHadrons","2");    
   xml->AddElementwithValue("nContributions","7");
   xml->AddElementwithValue("ContributionType","2");
   xml->AddElementwithValue("Label","?");    
   xml->PopElement();

   xml->PushElement("Process","ID",1); 
   xml->AddElementwithValue("Label","inclusive jet cross section");    
   xml->AddElementwithValue("Observable","d2 sigma/dy dpT  (pb)");    
   xml->PopElement();

   xml->PushElement("JetAlgorithm","ID",2); 
   xml->AddElementwithValue("Label","midpoint cone");    
   xml->AddElementwithValue("nParameters","2");    
   xml->PushElement("Parameter","ID",1); 
   xml->AddElementwithValue("Value","0.7");    
   xml->AddElementwithValue("Label","R_cone");    
   xml->PopElement();
   xml->PushElement("Parameter","ID",2); 
   xml->AddElementwithValue("Value","0.5");    
   xml->AddElementwithValue("Label","f_overlap");    
   xml->PopElement();
   xml->PopElement();

   xml->PopElement();
   //--------------------- Calculation
   xml->PushElement("Calculation");
   xml->AddElementwithValue("nxBins",nxtot);
   xml->AddElementwithValue("xBinStructure",1);

   xml->PushElement("xScheme","ID",2);  
   xml->AddElementwithValue("Label","sqrt(log 1/x)"); 
   xml->PopElement();

   xml->PushElement("pdfWeighting","ID",1); 
   xml->AddElementwithValue("Label","standard weighting"); 
   xml->PopElement();

   xml->PushElement("FactorisationScale","ID",1); 
   xml->AddElementwithValue("Label","pt in GeV"); 
   xml->AddElementwithValue("nVariations",nscalevar); 
   for(int i = 0;i<nscalevar;i++){
      xml->AddElementwithValue("Variation",mufscale[i]); 
   }
   xml->PopElement();

   xml->PushElement("RenormalisationScale","ID",1); 
   xml->AddElementwithValue("Label","pt in GeV"); 
   xml->AddElementwithValue("nVariations","1"); 
   for(int i = 0;i<1;i++){
      xml->AddElementwithValue("Variation",murscale[i]); 
   }
   xml->PopElement();

   int nbins =0;
   for(int i=0;i<nrap;i++){
      nbins += npt[i];
   }
   xml->AddElementwithValue("nBins",nbins); 
   xml->AddElementwithValue("nDimensions",2);

   xml->PushElement("Dimension","ID",1);
   xml->AddElementwithValue("Label","rapidity");    
   xml->AddElementwithValue("nBins",nrap);    
   for(int i=0;i<nrap+1;i++){
      xml->AddElementwithValue("Bound",raphigh[i]);    
   }
   xml->PopElement();

   xml->PushElement("Dimension","ID",2);
   xml->AddElementwithValue("Label","pt in GeV");    
   for(int i=0;i<nrap;i++){
      xml->PushElement("ParentBin","ID",i+1);
      xml->AddElementwithValue("nBins",npt[i]);    
      for(int j=0;j<npt[i]+1;j++){
         xml->AddElementwithValue("Bound",pthigh[i][j]);    
      }
      xml->PopElement();
   }
   xml->PopElement();
   
   xml->PushElement("Bins");
   int bincount = 1;
   for( int j = 0; j < nrap; j++) {
      for( int k = 0; k < npt[j]; k++) {
         xml->PushElement("Bin","ID",bincount++);

         xml->PushElement("Low","Dim",1);
         xml->AddValue(raphigh[j]);
         xml->PopElement();
         xml->PushElement("High","Dim",1);
         xml->AddValue(raphigh[j+1]);
         xml->PopElement();
         xml->PushElement("Low","Dim",2);
         xml->AddValue(pthigh[j][k]);
         xml->PopElement();
         xml->PushElement("High","Dim",2);
         xml->AddValue(pthigh[j][k+1]);
         xml->PopElement();
         xml->AddElementwithValue("xLimit",xlimit[j][k]);
         xml->AddElementwithValue("nScaleBins",nscalebins);
         for( int l = 0; l < nscalebins; l++) {
            xml->PushElement("ScaleBin","ID",l+1);
            xml->AddElementwithValue("muF",mufval[j][k][l]);
            xml->AddElementwithValue("muR",murval[j][k][l]);
            xml->PopElement();
         }

         xml->PopElement();
      }
   }   
   xml->PopElement();
   //--------------------- End


}


void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{
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

   //-------- analyze inclusive jets in jet loop -------
   for( int i = 1; i <= nj; i++) {

      //------- get jet properties
      double pt = pj[i].perp(); 
      if (pt>ptlow){ // check for low pt first, will fail most of the time
         double rap = abs(pj[i].prapidity());
        
         //------ determine y and pt bin ---
         int rapbin = -1;
         for( int j = 0; j < nrap; j++) {
            if (rap >= raphigh[j] && rap < raphigh[(j+1)]) rapbin=j;
         }
         if (rapbin >=0 ) {
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
                  printf("Warning: xmin (%f) < xlimit (%f) at pt=%f GeV y=%f \n ",
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
               

               // loop over scale variations (only >1 for NLO)
               for(int scalevar=0; scalevar<nscalevar;scalevar++){

                  // loop over scale bins
                  for(int scalebin=0; scalebin<nscalebins;scalebin++){
                     
                     double mur2=murscale[scalevar]*murscale[scalevar]
                        *murval[rapbin][ptbin][scalebin]*murval[rapbin][ptbin][scalebin];
                     double muf2=mufscale[scalevar]*mufscale[scalevar]
                        *mufval[rapbin][ptbin][scalebin]*mufval[rapbin][ptbin][scalebin];
                     weight_hhc wt = amp(mur2,muf2);

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
                  weights[rapbin][ptbin][nxmax][nxmin][scalevar][scalebin] += 
		    (1-deltamax)*(1-deltamin)*wt;

                  // -- avoid writing behind  nxmax=(nxtot-1)
                  if (nxmax<(nxtot-1)) {
                     weights[rapbin][ptbin][nxmax+1][nxmin][scalevar][scalebin] += 
		       deltamax * (1-deltamin)*wt;
                     weights[rapbin][ptbin][nxmax+1][nxmin+1][scalevar][scalebin] += 
		       deltamax * deltamin*wt;
                  }

                  // -- avoid writing behind  nxmin=nxmax
                  if (nxmin<nxmax) {             // should work with half-table
                     weights[rapbin][ptbin][nxmax][nxmin+1][scalevar][scalebin] += 
		       (1-deltamax)*deltamin*wt;
                  }

                  // -- compensate for missing contrib.: nxmin>nxmax (use symmetry)
                  //     ---> project  (nxmax,nxmin+1) back to (nxmax+1,nxmin)
                  if (nxmin==nxmax && nxmin<(nxtot-1)) { 
                     double buffer;
                     buffer = wt[1];
                     wt[1] = wt[2];
                     wt[2] = buffer; 
                     weights[rapbin][ptbin][nxmax+1][nxmin][scalevar][scalebin] += 
		       (1-deltamax)*deltamin*wt;
                  }
                  }//-end loop scalebin
               }//-end loop scalevar
            } //-end: IF in pT range
         }  // - end: IF in rapidity range
      }     // - end: IF ptlow-cut
   }        // - end jet loop

}

void UserHHC::writetable(){

   xml->ReplaceElement("Contribution","ID",1);
   xml->AddElementwithValue("Order",itype+1);
   xml->AddElementwithValue("AlphasPower",itype+2);
   char buffer[64];
   switch(itype){
   case amplitude_hhc::lo: 
      strcpy(buffer,"LO");
      break;
   case amplitude_hhc::nlo: 
      strcpy(buffer,"NLO");
       break;
    default:
      strcpy(buffer,"not known");
   }   
   xml->AddElementwithValue("Label",buffer);
   xml->AddElementwithValue("nEvents",nevents);
   xml->PopElement();

   xml->ReplaceElement("Table");
    
   if(textoutput){
      xml->PushElement("Format","ID",1);
      xml->AddValue("Numbers as text speparated by spaces");
      xml->PopElement();     
      for(int i=0;i<nrap;i++){ // rapidity
         for(int j=0;j<npt[i];j++){ // pt
            for(int k=0;k<nxtot;k++){ // xmax
               for(int l=0;l<k+1;l++){ // xmin
                  for(int scalevar=0; scalevar<nscalevar;scalevar++){ // scalevar
                     for(int m=0;m<7;m++){     //subprocesses
                        for(int n=0;n<nscalebins;n++){     //scalebin
                           sprintf(buffer,"%.15g ",weights[i][j][k][l][scalevar][n][m]);
                           xml->AddValue(buffer);
                        }
                      }
                  }
               }
            }
         }
      }
   }else{
      xml->PushElement("Format","ID",2);
      xml->AddValue("Numbers as BASE64 encoded IEEE 64-bit floating point");
      xml->PopElement();
      unsigned int base64length;
      XMLByte *base64buffer; 
      for(int i=0;i<nrap;i++){ // rapidity
         for(int j=0;j<npt[i];j++){ // pt
            for(int k=0;k<nxtot;k++){ // xmax
               for(int l=0;l<k+1;l++){ // xmin
                  for(int scalevar=0; scalevar<nscalevar;scalevar++){  // scalevar
                     for(int m=0;m<7;m++){     //subprocesses
                        for(int n=0;n<nscalebins;n++){     //scalebin
                           double number = weights[i][j][k][l][scalevar][n][m];
                           base64buffer = Base64::encode((XMLByte *)&number,sizeof(number),&base64length);
                           base64buffer[base64length-1] = 0;
                           xml->AddValue((char*)base64buffer);
                           XMLString::release(&base64buffer);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   xml->PopElement();
   xml->Write(tablefilename.c_str());

}

void UserHHC::end_of_event(){
   // let NLOJET++ store its results
   user_hhc::end_of_event();
   
   nevents += 1;
   //-------- store table
   if (( (unsigned long)nevents % nwrite)==0){
      const unsigned long startMillis = XMLPlatformUtils::getCurrentMillis();
      printf ("No. events: %.2G writing table... ",nevents);
      writetable();
      const unsigned long endMillis = XMLPlatformUtils::getCurrentMillis();
      printf(" done (%ld ms).\n",endMillis-startMillis);
   }
 }

void UserHHC::phys_output(const std::basic_string<char>& __file_name, 
                          unsigned long __save, bool __txt) 
{
   // Suppress output of NLOJET++ files
   user_hhc::phys_output("",2000000000,false);

   // Our outputfile
   tablefilename = __file_name +".xml";

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

