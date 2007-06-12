// Fast Computation of Jet Cross Sections in Hadron-Induced Processes

// last modification 2006/02/06 - 6pm(GE) TK

//------ DON'T TOUCH THIS PART! ------
#include <phasespace.h>
#include <process.h>
#include <jetfunc.h>
#include <qcdlib.h>
#include "cteq6dis.h"

#include <iomanip>              // for ASCII output for table

//----- used namespaces -----
using namespace nlo;
using namespace std;


//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
user_dis * userfunc();
genevent_dis::input * pscuts();

//----- array of the symbols symbols -----
struct { 
   const char *name;
   void *address;
} user_defined_functions[] = 
   {
      //   process index: 1 --> dis (only jets in the final state)
      {"procindex", (void *) "1"},
    
      //   input function 
      {"inputfunc", (void *) inputfunc},
    
      //   user defined functions
      {"userfunc",  (void *) userfunc},

      //   function for defining the phase space cuts
      {"pscuts",    (void *) pscuts},
    
      //  end of the list
      {0, 0}
   };
class physcuts : public genevent_dis::input
{
public:
  //  constructor
  physcuts() { 
    setcuts();
    checkcuts(); 
  }
  
  // generate the bjorken variables 'xB and y'
  double operator()(phasespace_dis *, double&, double&);
  
  // energy of the incomong lepton and hadron in the laboratory frame
  double energy_lepton() const { return El;}
  double energy_hadron() const { return Eh;}
  double getymax() const { return ymax;}

  
private:
  //   data members
  double El, Eh;
  double Q2min, Q2max;
  double xmin, xmax;
  double ymin, ymax;
  
  //   set the phase space cuts
  void setcuts();
  void checkcuts();
};

physcuts *thecuts;


void physcuts::checkcuts()   
{
  double s = 4.0*El*Eh;
  
  if(Q2min <= 0.0)  throw "Q2min must be greater than zero";
  if(Q2min > Q2max) throw "Q2min must be less than Q2max";
  if(xmin > xmax)   throw "xmin must be less than xmax";
  if(ymin > ymax)   throw "ymin must be less than ymax";
  
  if(xmax == xmin && ymax == ymin && Q2max == Q2min)
    throw "at most two variables can be constrained";
  
  if(s*xmax*ymax < Q2min || s*xmin*ymin > Q2max)
    throw "no phase space aviable";
  
  if(s*xmin*ymin > Q2min) Q2min = s*xmin*ymin;
  if(s*xmax*ymax < Q2max) Q2max = s*xmax*ymax;
}

double physcuts::operator()(phasespace_dis *ps, double& x, double& y) 
{
  double Q2, yjac, qjac, s = 4.0*El*Eh;
  
  if(Q2min == Q2max) {
    Q2 = Q2min;
    qjac = 1.0;
  } else {
    qjac  = std::log(Q2max/Q2min);
    Q2 = Q2min*std::exp(ps -> random()*qjac);
    qjac *= Q2;
  } 
  
  double yn = ymin;
  double yx = ymax;
  
  if(ymin*xmax*s < Q2) yn = Q2/(s*xmax);
  if(ymax*xmin*s > Q2) yx = Q2/(s*xmin);
  
  if(yn == yx) {
    y = yn;
    yjac = 1.0;
  } else if(yn < yx) {
    yjac = std::log(yx/yn);
    y = yn*std::exp(ps -> random()*yjac);
  } else throw "no phase space aviable";
  
  x = Q2/(s*y);
  
  return qjac*yjac/s;
}

//------ USER ROUTINE ------
#include "kt-et-10-dis.h"

class UserDIS : public user_dis
{
 public:
   //   init and user function
   void initfunc(unsigned int);
   void userfunc(const event_dis&, const amplitude_dis&);
   virtual void phys_output(const std::basic_string<char>& __file_name, 
			   unsigned long __save = 10000UL, bool __txt = false);
   void end_of_event();  

 private:
   // algorithms
   kt_et_10_dis jetclus;
  
   // the jet structore in breit & lab. frame
   bounded_vector<lorentzvector<double> > pjb, pjl; 
  
   //  event in breit frame
   event_dis pbreit;

   bool nlo;       // Is the job running at LO or NLO?

   unsigned int iref;     //  switch for reference mode
   unsigned int refscale; //  select which scale variaton is used in ref-table

   double unitfactor;   // factor to convert from nb to units in publication
   int nrap;       // No of rapidity bins 
   double *raphigh;  // array for rapidity boundaries
   int *npt;       // No of pT bins in each y range
   vector< vector<double> >pthigh;   // array for pT boundaries
   double ptlow;    // lowest pt considered 
   
   int nscalevar;             // number of scale variations (mu_r,mu_f) in NLO
   int nscalebin;             // number of points on the scale axis, where to calculate
   vector <double> murscale;  // overall scale factor for renormalization scale
   vector< vector< vector<double> > >murval; // array for renormalization scale values
   vector< vector< vector< vector<double> > > >murvaltrans; // array for renormalization scale values after log/log(mu/mu0) transformation
   vector <double> mufscale;       // overall scale factor for fact. scale
   vector< vector< vector<double> > >mufval; // array for factorization scale values 
   vector< vector< vector< vector<double> > > >mufvaltrans; // array for factorization scale values after log/log(mu/mu0) transformation
   double mu0scale; // Reference scale for og/log(mu/mu0) transformation

   int nxtot;      // no of xbins 
   vector< vector<double> >xlimit; // array for lower x limits
   vector< vector<double> >hxlim; // array for lower x limits, after transformation
   vector< vector<double> >xsmallest; // array for smallest actual x values
   //    array for the weights MAIN ARRAY
   vector <vector <vector< vector < vector < weight_dis > > > > >weights; 
   
  // ===== variables for the b-cubic interpolation =====
  // - the relative distances to the four nearest bins
  vector<double> cm ; 
  // - the weights for the cubic eigenfunctions (1-dim)
  vector<double> cefm; 

   double nevents;        // No of events calculated so far
   unsigned long nwrite;  // No of events after to write out the table

   pdf_cteq6 pdf;  //   pdf
 
   basic_string<char> tablefilename; // The table file to write to
   bool textoutput; // If true, the table is written in plain ASCII instead of BASE64 encoded doubles (later for XML)
   amplitude_dis::integral_type itype; // Born, NLO etc.
   
   time_t start_time;
   
  // boost the jet momenta back to the laboratory frame
   void boost_back_to_lab(const event_dis&);
   void writetable();
};

user_dis * userfunc() {
   return new UserDIS;
}

genevent_dis::input * pscuts() {
   thecuts = new physcuts;
   return thecuts;
}

void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
  //  number of jets
  nj = 3U;

  //  number of the up and down type flavours
  nu = 2U;
  nd = 3U;
} 

void physcuts::setcuts()
{
  // energy of the incomong lepton and 
  // hadron in the laboratory frame
  El = 27.5;    // GeV
  Eh = 820.0;     // GeV
  
  //  Q^2 cuts
  Q2min = 5.0;   // GeV^2
  Q2max = 85.0;  // GeV^2 

  //   xB cuts
  xmin = 1.e-4;
  xmax = 4.e-3;
  
  //   y cuts
  ymin = 0.1;
  ymax = 0.7;
}

void UserDIS::initfunc(unsigned int)
{

   unsigned int nj;
   unsigned int nu;
   unsigned int nd;
   inputfunc(nj,nu,nd);
   double El = thecuts->energy_lepton();
   double Eh = thecuts->energy_hadron();     
   double s = 4.0*El*Eh;

   // ********************************************************
   // ********** switch for reference mode on/off ************
   iref = 0;       //  switch for reference mode
   //   iref = 1;       //  switch for reference mode
                  //   0: standard fastNLO table only
                 //    1: include 2nd "reference table" (a_s/PDFs)
   refscale = 1;   // which of the scalevariations is used in ref-table?

   //unitfactor = 1000000.0;  // for fb
   unitfactor = 1000.0;  // for pb << used in H1 publication
   //   unitfactor = 1.0;  // for nb     

   //Set up binning  
   nrap   = 11;                 // no of bins in deltaeta2 aka rapidity
   if (iref==1) nrap *= 2;  // -> in reference mode: double No. rap bins

   raphigh = new double[nrap+1];  //----- array for rapidity boundaries
   npt     = new  int[nrap];      // nrap bins in rapid.-each npt[irap] pT bins

   //Define binning in deta2
   raphigh[0]=   0.0;         
   raphigh[1]=   0.6;         
   raphigh[2]=   1.4;         
   raphigh[3]=   3.0;         
   //For deta1<1
   raphigh[4]=   0.0;         
   raphigh[5]=   0.6;         
   raphigh[6]=   1.4;         
   raphigh[7]=   3.0;         
   //For deta1>1
   raphigh[8]=   0.0;         
   raphigh[9]=   0.6;         
   raphigh[10]=   1.4;         
   raphigh[11]=   2.5;         

   if (iref==1)      // -> in reference mode: copy rapidity definitions
     for(int i=0;i<nrap/2;i++){
      raphigh[i+nrap/2+1] = raphigh[i+1];
   }

   pthigh.resize(nrap);
   for(int i=0;i<nrap;i++){  
      //Define binning in pt (mean pt of dijets)
      npt[i]=1;
      pthigh[i].resize(npt[i]+1);
      //----- array for pt boundaries
      pthigh[i][0] = 6.0;
      pthigh[i][1] = 100.;
   }
   // lowest pT value in sample
   ptlow = 6.0;
        
   // Bins in x
   nxtot = 30;

   // NLO scale variations - for scales in GeV
   nscalevar = 3;
   // no of scalebins = degree of Bernstein polynomials + 1
   nscalebin = 4;
   
   mu0scale = .25; // .25 Gev as reference scale

   murscale.resize(nscalevar);
   mufscale.resize(nscalevar);
   murscale[0] = 0.5;  mufscale[0] = 0.5;  
   murscale[1] = 1.0;  mufscale[1] = 1.0;
   murscale[2] = 2.0;  mufscale[2] = 2.0;

   xlimit.resize (nrap);
   hxlim.resize (nrap);
   xsmallest.resize (nrap);     // test: find smallest x values
   murval.resize (nrap);
   mufval.resize (nrap);
   murvaltrans.resize (nrap);
   mufvaltrans.resize (nrap);

   weights.resize (nrap);
   for( int j = 0; j < nrap; j++) {
      xlimit[j].resize(npt[j]);
      hxlim[j].resize(npt[j]);
      xsmallest[j].resize(npt[j]);
      murval[j].resize (npt[j]);
      mufval[j].resize (npt[j]);
      murvaltrans[j].resize (npt[j]);
      mufvaltrans[j].resize (npt[j]);
      weights[j].resize(npt[j]);
      for( int k = 0; k < npt[j]; k++) {
         murval[j][k].resize (nscalebin);
         mufval[j][k].resize (nscalebin);
         murvaltrans[j][k].resize (nscalevar);
         mufvaltrans[j][k].resize (nscalevar);
         for( int l = 0; l < nscalevar; l++) {
            murvaltrans[j][k][l].resize (nscalebin);
            mufvaltrans[j][k][l].resize (nscalebin);
         }        
         // Setup the weights array
         weights[j][k].resize(nxtot);
         for( int l = 0; l < nxtot; l++) {
            // scale variation
            weights[j][k][l].resize(nscalevar);
            for(int scalevar=0; scalevar<nscalevar;scalevar++){
               weights[j][k][l][scalevar].resize(nscalebin);
            }
         }

         // - Setup the xlimit array - computed from kinematic constraints
         xlimit[j][k]= 1.e-2; 
         hxlim[j][k] = -sqrt(-log10(xlimit[j][k]));
         xsmallest[j][k]=0.999;

         // - Setup the murval and mufval arrays - for scale bins
	 double dpt = (pthigh[j][k+1]-pthigh[j][k]);
	 double lptlo = pthigh[j][k] + dpt*0.00;   // - starting points
	 double lpthi = pthigh[j][k+1] - dpt*0.9;
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
	 hxlim[j+nrap/2][k] = hxlim[j][k];
         for(int l=0;l<nscalebin;l++){
            murval[j+nrap/2][k][l]= murval[j][k][l]; 
            mufval[j+nrap/2][k][l]= mufval[j][k][l];
            for( int m = 0; m < nscalevar; m++) {
               murvaltrans[j+nrap/2][k][m][l] =  murvaltrans[j][k][m][l];
               mufvaltrans[j+nrap/2][k][m][l] =  mufvaltrans[j][k][m][l];
            }
         }
      }
   }

   // print x-limit values at the begin of the job
   printf ("(Q2, pt) array for this job:\n");
   printf ("#Q2 #pt xlimit pt_high rap_high \n");
   for( int j = 0; j < nrap; j++) {
      for( int k = 0; k < npt[j]; k++) {
         printf("%3d %3d   %8.6f %8.1f %8.1f \n",j,k,
		xlimit[j][k],pthigh[j][k],raphigh[(j+1)]);
      }
   }

   // ===== variables for the bi-cubic interpolation =====
   // - the relative distances to the four nearest bins
   cm.resize (4);
   // - the weights for the cubic eigenfunctions (1-dim)
   cefm.resize (4);

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
   cout<<"         *** scenario fnh1001 ***"<<endl;
   cout<<" "<<endl;
   cout<<"        table file "<<tablefilename<<endl;
   cout<<"        store table after "<<nwrite<<" events"<<endl;
   cout<<"        sqrt(s)= "<<sqrt(s)<<endl;
   cout<<"        No. x-bins: "<<nxtot<<endl;
   cout<<"        No. Q2 regions: "<<nrap<<endl;
   cout<<"        No. of pT bins in each Q2 region:"<<endl;
   for( int j = 0; j < nrap; j++) {
      cout<<"          Q2 "<<j<<": "<<npt[j]<<endl;
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

extern"C" double xalfaem_(double *);

double xalpha_em(double mq2) {
  return xalfaem_(&mq2);
}


void UserDIS::userfunc(const event_dis& p, const amplitude_dis& amp)
{
  //----- cuts in the laboratory & breit frame -----
   const double degtorad = 3.141592654/180.;
   //cuts on electron in lab frame
   const double eletamin = -log(tan(175.*degtorad/2.));
   const double eletamax = -log(tan(156.*degtorad/2.));
   if(p[-2].T()<10.) return; // Cut on Electron energy in lab
   if(p[-2].prapidity()<eletamin || p[-2].prapidity()>eletamax) return; // Cut on Electron theta in lab


   double alem = 0;
   double x=0.;  // pdf-x
   double hx   = 0.0; //  x after transformation
   double hxone   = 0.0; // x=1 after transformation

   const int rapbin  =  0; 
   int rapbin1 = -1; // In this scenario, up to 2 plots can be filled per event
   int rapbin2 = -1;

   double Q2 =0.;
   double binwidth1 = 0.0;
   double binwidth2 = 0.0;
   // - to reweight the Eigenfunctions (universal PDF)
   double reweight = 1.0;  

   //---------------------------------------------------
   //--------- start event processing ------------------
   //---------------------------------------------------


   // LO or NLO?
   itype = amp.integral();

   //----- copy the momenta and boost to the breit frame ----
   pbreit = p; lab_to_breit(pbreit);
   //----- jet structure in breit frame -----
   pjb = jetclus(pbreit);
   int nj = pjb.upper(); 

   if(nj>4){
      printf("This usermodule is not prepared for more than 4 partons.\n");
      exit(1);
   }
   if(nj<3) return; // At least 3 jets needed
   //----- jet structure in laboratory frame -----
   pjl = pjb; boost_back_to_lab(p);

   const double etamin = p[-2].prapidity(); // eta_jets > eta_electron
   const double etamax = -log(tan(7.*degtorad/2.)); //  pseudo-rapidity cut on jets, Lab frame

   int jetindex[4];
   int jetcount = 0;
   for( int i = 1; i <= nj; i++) { // identify jets fulfilling eta and E_t cut 
      double eta = pjl[i].prapidity();
      double pt = pjl[i].perp();
      if(eta > etamin && eta < etamax && pt>ptlow ){
         jetindex[jetcount++] = i;
      }
   }

   if(jetcount<3) return; // at least 3 jets needed


   for(int i = 0; i < jetcount; i++) { // order in eta
      for(int j = 0; j < i; j++) {
         if(pjl[jetindex[i]].prapidity()>pjl[jetindex[j]].prapidity()){ // swap
            int temp = jetindex[i];
            jetindex[i] = jetindex[j]; 
            jetindex[j] = temp; 
         }
      }
   }

   int fwdjet = -1; // Index of forward jet
   for(int i = 0; i < jetcount; i++) { // look for "forward" jet
      int j = jetindex[i];
      if(pjl[j].T()>(0.035*p[hadron(0)].T()) && pjl[j].prapidity()>-log(tan(20.*degtorad/2.))){
         fwdjet = i;
         break;
      }
   }
   if (fwdjet<0) return; // No forward jet found

   int jet0 = -1; // Index of forward jet
   int jet1 = -1; // Index of other jet
   int jet2 = -1; // Index of other jet
   int jet3 = -1; // Index of other jet, not used in analysis

   // how many jets are left between the forward jet and the electron?
   switch(jetcount-fwdjet-1){
   case 0:
      return;
      break;
   case 1:
      return;
      break;
   case 2:
      jet0 = jetindex[fwdjet];
      jet1 = jetindex[fwdjet+1];
      jet2 = jetindex[fwdjet+2];
      break;
   case 3:
      jet0 = jetindex[fwdjet];
      jet1 = jetindex[fwdjet+1];
      jet2 = jetindex[fwdjet+2];
      jet3 = jetindex[fwdjet+3];
      // take 2 jets with highest pt
      if(pjl[jet3].perp()<pjl[jet2].perp() && pjl[jet3].perp()<pjl[jet1].perp() ){
         ;} // Already ok
      else{
         if(pjl[jet2].perp()<pjl[jet1].perp() && pjl[jet2].perp()<pjl[jet3].perp() ){ // Swap 2 and 3
            int temp = jet2;
            jet2 = jet3; 
            jet3 = temp; 
         }else{
            if(pjl[jet1].perp()<pjl[jet2].perp() && pjl[jet1].perp()<pjl[jet3].perp() ){ // Put 1 to the end
               int temp = jet1;
               jet1 = jet2; 
               jet2 = jet3;
               jet3 = temp;
            }else{
               printf("Inconsistency when finding dijets, jetcount=%d.\n",jetcount);
               exit(1);
            }
         }
      }
      // After this procedure the dijets are jet1 and jet2, jet3 will be not used 
      break;
   default:
      printf("Inconsistency when finding dijets, jetcount=%d.\n",jetcount);
      exit(1);
      break;
   }
   double meanpt = sqrt(1./3. *(pow(pjl[jet0].perp(),2)+pow(pjl[jet1].perp(),2)+pow(pjl[jet2].perp(),2))); // the scale
   double deta1 = pjl[jet1].prapidity()-pjl[jet2].prapidity(); // delta rapidity between dijets
   double deta2 = pjl[jet0].prapidity()-pjl[jet1].prapidity(); // delta rapidity between forward and backward jets
 
  for( int j = 0; j < 3; j++) {
      if (deta2 >= raphigh[j] && deta2 < raphigh[(j+1)]){
         rapbin1=j;
         binwidth1 = raphigh[(j+1)] - raphigh[j];
         break;
      }
  }
  if(deta1<1.0){
     for( int j = 4; j < 7; j++) {
        if (deta2 >= raphigh[j] && deta2 < raphigh[(j+1)]){
           rapbin2=j;
           binwidth2 = raphigh[(j+1)] - raphigh[j];
           break;
        }
     }
  }else{
     for( int j = 8; j < 11; j++) {
        if (deta2 >= raphigh[j] && deta2 < raphigh[(j+1)]){
           rapbin2=j;
           binwidth2 = raphigh[(j+1)] - raphigh[j];
           break;
        }
     }
  }
      
  if(rapbin1==-1 && rapbin2==-1) return; // Anything to fill?

  Q2 = 2*(p[-2]*p[-1]);
  alem = xalpha_em(Q2);
  x = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);
  hx   = -sqrt(-log10(x)); // x binning
  hxone   = 0.0;
  reweight = 1.0/sqrt(x) * (1.0-0.99*x);  
  const int ptbin  = 0; // this scenatio has no pt bins

  //---------- compute weight, fill fastNLO array

  // test if x_min in event is lower than x_limit
  //      if yes -> make big warning!!! 
  //      -> need to change x_limit values
  if (x<xlimit[rapbin][ptbin]){
     printf("Error: x (%f) < xlimit (%f) at pt=%f GeV Q2=%f, Q2#=%d , pt#=%d \n ",
            x,xlimit[rapbin][ptbin],meanpt,Q2,rapbin,ptbin);
     exit(1);
  }

  // *******  identify smallest x value in each pT/y bin  *******
  //               -> to optimize the x-limit values
//   if (x<xsmallest[rapbin][ptbin]*0.98){  // 0.98 reduces output
//      xsmallest[rapbin][ptbin] = x;
//      if((itype==amplitude_dis::lo && nevents> 1000000000)|| // 7h 
//         (itype==amplitude_dis::nlo && nevents> 120000000)){  // 10h
//         //             if((itype==amplitude_dis::lo && nevents> 600000)|| // 
//         //                (itype==amplitude_dis::nlo && nevents> 60000)){  // 
//         cout<<" "<<endl;
//         cout<<">>>>>> smaller x found in bin  "<<rapbin<<"  "<<ptbin
//             <<"   -  "<<x<<"  "<<xlimit[rapbin][ptbin]<<"  :  "<<
//            (x/xlimit[rapbin][ptbin])<<"  >> "<<nevents<<endl;
//         for( int j = 0; j < nrap; j++) {
//            for( int k = 0; k < npt[j]; k++) {
//               cout<<"    bin "<<j<<"  "<<k<<"  "<<xsmallest[j][k]
//                   <<"  "<<xlimit[j][k]
//                   <<"   "<<pthigh[j][k]<<"  "<<raphigh[(j+1)]<<"  :  "
//                   <<(xsmallest[j][k]/xlimit[j][k])<<endl;
//            }
//         }
//      }
//   }
  // ------------  end: identify smallest x  -------------


  // **********  determine x_ij position in grid  ************
  //--- determine fractional contributions
  double hxlimit = hxlim[rapbin][ptbin];

  // define the x-bin number in the range  [0:ntot[
  int nx = int(nxtot *(hx-hxlimit)/(hxone-hxlimit));

  //-- relative distances in h(x): deltam
  double delta  = (hxone-hxlimit)/nxtot;
  double hxi =hxlimit+double(nx)/double(nxtot)*(hxone-hxlimit);
  double deltam = (hx-hxi)/delta;
  if(deltam>1.0 || deltam<0.0 ){
     cout<<" -> delta is off: "<<deltam<<endl;
  }                             


  // ===== variables for the bi-cubic interpolation =====
  // === the relative distances to the four nearest bins
  cm[0] = deltam+1.0;
  cm[1] = deltam;
  cm[2] = 1.0-deltam;
  cm[3] = 2.0-deltam;

  // === the weights for the cubic eigenfunctions (1-dim)
  //   - linear interpolation in 1st and last =(nxtot-1) bins
  //   - cubic approximation in the middle 

  if (nx==0 || nx==(nxtot-1)) { //linear in 1st and last bin
     cefm[0] = 0.0;
     cefm[1] = 1.0-cm[1];
     cefm[2] = 1.0-cm[2];
     cefm[3] = 0.0; }
  else {                              // cubic in the middle
     cefm[1]=1.0-2.5*cm[1]*cm[1]+1.5*cm[1]*cm[1]*cm[1];
     cefm[2]=1.0-2.5*cm[2]*cm[2]+1.5*cm[2]*cm[2]*cm[2];
     cefm[0]=2.0 - 4.0*cm[0] + 2.5*cm[0]*cm[0]
        - 0.5*cm[0]*cm[0]*cm[0];
     cefm[3]=2.0 - 4.0*cm[3] + 2.5*cm[3]*cm[3]
        - 0.5*cm[3]*cm[3]*cm[3];
  }
	       
  // loop over scale variations for NLO
  int scalevarmax;
  if(itype==amplitude_dis::lo){
     scalevarmax=1;
  }else{
     scalevarmax=nscalevar;
  }
  amp.pdf_and_qcd_coupling(0,1.0);
  for(int scalevar=0; scalevar<scalevarmax;scalevar++){
     double mur2 = murscale[scalevar]*murscale[scalevar]
        * meanpt*meanpt;
     double muf2 = mufscale[scalevar]*mufscale[scalevar]
        * meanpt*meanpt;
     weight_dis wt = amp(mur2,muf2);
     //--------------------------------
     // transform to MWs pdf linear combinations
     double pdfsigma = 1./3.*(-wt[1] + 4.*wt[2]);
     double pdfdelta = 3.*(wt[1]-wt[2]);
     wt[1] = pdfsigma;
     wt[2] = pdfdelta;
                     
     wt *= alem*alem*389385.730;
     wt *= reweight*reweight*reweight;
     wt *= unitfactor;
     //make Bernstein interpolation for scalebins
     double hmu   = log(log(murscale[scalevar]*meanpt/mu0scale));
     double tscale = (hmu-murvaltrans[rapbin][ptbin][scalevar][0]) /
        (murvaltrans[rapbin][ptbin][scalevar][nscalebin-1] -
         murvaltrans[rapbin][ptbin][scalevar][0]); 
     for(int scalebin=0; scalebin<nscalebin;scalebin++){
        double deltascale=0.;
        switch(nscalebin){ // nscalebin = degree of bernstein polynomial+1
        case 1:
           deltascale = 1.0;
           break;
        case 2:
           switch(scalebin){                  
           case 0:
              deltascale = (1.0-tscale);
              break;
           case 1:
              deltascale = tscale;
              break;
           default: printf("Error: scalebin = %d not existant.\n",scalebin);exit(1);
           }
           break;
        case 3:
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
           default: printf("Error: scalebin = %d not existant.\n",scalebin);exit(1);
           }
           break;
        case 4:
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
           default: printf("Error: scalebin = %d not existant.\n",scalebin);exit(1);
           }                  
           break;
        default:
           printf("Error: nscalebin = %d not supported.\n",nscalebin);
           exit(1);
        }
        // ** ----------------------------------------------------
        // ** now fill table
        // ** loop over all 4 points that receive contributions
        for( int i1 = 0; i1 < 4; i1++) {
           int im = nx +i1 -1;   // the target index
               
           // - ignore if x coordinate is < 0 (contrib is zero!)
           // - ignore if x coordinate is >= nxtot
           //     (contrib. > nxtot are zero / =nxtot PDF is zero)
           if (im<nxtot && im>=0) {
               if(rapbin1>=0)
                  weights[rapbin1][ptbin][im][scalevar][scalebin] +=  cefm[i1]  * deltascale / binwidth1 * wt;
               if(rapbin2>=0)
                  weights[rapbin2][ptbin][im][scalevar][scalebin] +=  cefm[i1]  * deltascale / binwidth2 * wt;
           }
        }
     }

                     
  }//-end loop scalevar

  // ** ------------------------------------------------------
  // **   only in "reference mode":
  // **         -> fill table with pdfs and alphas
  // **         -> only for a single scale (index=refscale)
  // **         -> don't throw away contributions at x=1
  // **            store in first bin No. = 0
  // **                (needed for reference)
  // **
  if (iref==1) {  // -> in reference mode
     double mur2 = murscale[refscale]*murscale[refscale]
        * meanpt*meanpt; // "real" scale is the jet pt
     double muf2 = mufscale[refscale]*mufscale[refscale]
        * meanpt*meanpt;
     amp.pdf_and_qcd_coupling(pdf, 389385.730);
     weight_dis wt = amp(mur2,muf2);

     wt *= alem*alem*unitfactor;
      if(rapbin1>=0)
         weights[rapbin1+nrap/2][ptbin][0][0][0] +=   wt/ binwidth1 ;
      if(rapbin2>=0)
         weights[rapbin2+nrap/2][ptbin][0][0][0] +=   wt/ binwidth2;
  } //-end reference mode: a_s/PDF table filling 
         
}

void UserDIS::boost_back_to_lab(const event_dis& p)
{ 
  double x = p[-1].T()/p[hadron(0)].T();
  double bz = (1.0 - x)/(1.0 + x);
  threevector<double> bVec = -((pbreit[-1] + pbreit[hadron(0)]).boostVector());
  lorentzvector<double> p0(pbreit[hadron(0)]);
  
   p0.boost(bVec);
   double phi = p0.phi(), theta = p0.theta();
   unsigned int njet = pjl.upper();
    for(unsigned int i = 1; i <= njet; i++) {
       pjl[i].boost(bVec);
       pjl[i].rotateZ(-phi);
       pjl[i].rotateY(-theta); 
       pjl[i].boost(0.0, 0.0, bz);
    }
}

void UserDIS::writetable(){
#define WRITE(n) table.write(reinterpret_cast<char *>(&n),sizeof(n))
   //#define WRITE(n) table << setprecision(40) << n << endl

   int marker = 1234567890; //used to separate section in the table
   
   fstream table(tablefilename.c_str() ,ios::out|ios::binary); // open file
 
   // ireaction
   int ireaction = 1;
   WRITE(ireaction);

   double El = thecuts->energy_lepton();
   double Eh = thecuts->energy_hadron();     
   double s = 4.0*El*Eh;

   // Ecms
   s= sqrt(s);
   WRITE(s);

   // five strings with table content
   table << "dsigma-deta2_(pb)" << endl;
   table << "hep-ex/" << endl;
   table << "H1_Collaboration" << endl;
   table << "-" << endl;
   table << "-" << endl;

  //iproc
   int iproc = 1; // incl. jets
   WRITE(iproc);

   //ialgo
   int ialgo = 1; // incl. kt
   WRITE(ialgo);

   //JetResol1
   double JetResol1 = 1.0; // D parameter
   WRITE(JetResol1);

   //JetResol2
   double JetResol2 = 0.0; // not used
   WRITE(JetResol2);

   // relative order
   int nord = itype+1; // LO: itype=0,  NLO: itype=1
   WRITE(nord);

   // absolute order
   int npow = nord+1; // 3-jets in DIS -> LO=alpha_s^2,  NLO=alpha_s^3 
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
   table << "mean_ET_of_jets" << endl;

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
   if(itype==amplitude_dis::nlo){   // scale variations are only filled above LO
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

   for(int i=0;i<nrap;i++){ // Q2
      for(int j=0;j<npt[i];j++){ // pt
         for(int k=0;k<nxtot;k++){ // x
	      for(int m=0;m<3;m++){     //subprocesses
		// loop over scale variations for NLO
                 for(int scalevar=0; scalevar<scalevarmax;scalevar++){
                    for(int scalebin=0; scalebin<nscalebin;scalebin++){
                       WRITE(weights[i][j][k][scalevar][scalebin][m]);
                    }
                 }
              }
         }
      }
   }
   
   WRITE(marker);// ------------------END of table

   table.close();
}

void UserDIS::end_of_event(){
   // let NLOJET++ store its results
   user_dis::end_of_event();
   
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

void UserDIS::phys_output(const std::basic_string<char>& __file_name, 
                          unsigned long __save, bool __txt) 
{
   // Suppress output of NLOJET++ files
   user_dis::phys_output("",2000000000,false);
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

