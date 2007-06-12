// Fast Computation of Jet Cross Sections in Hadron-Induced Processes

// last modification 2005/11/18 - 6pm(GE) TK

//------ DON'T TOUCH THIS PART! ------
#include <phasespace.h>
#include <process.h>
#include <jetfunc.h>
#include <qcdlib.h>
#include <iomanip>

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
#include "kt-et-10.h"

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

   //   pdf
   //   pdf_cteq5m pdf;
   
   // algorithms
   kt_et_10 jetclus;
  
   // the jet structore in breit & lab. frame
   bounded_vector<lorentzvector<double> > pjb, pjl; 
  
   //  event in breit frame
   event_dis pbreit;

   bool nlo;       // Is the job running at LO or NLO?
   // binning
   
   int nQ2;       // No of Q2 bins 
   double *Q2high;  // array for Q2 boundaries
   int *npt;       // No of pT bins in each Q2 range
   vector< vector<double> >pthigh;   // array for pT boundaries
   double ptlow;    // lowest pt considered 
   
   int nscalevar;                 // number of scale variations (mu_r, mu_f) for NLO
   vector <double> murscale;        // overall scale factor for renormalization scale
   vector< vector<double> >murval;   // array for renormalization scale values
   vector <double> mufscale;        // overall scale factor for factorization scale
   vector< vector<double> >mufval;    // array for factorization scale values

   int ntot;      // no of xbins 
   vector< vector<double> >xlimit; // array for lower x limits
   vector< vector<double> >xsmallest; // array for smallest actual x values
   vector <vector< vector < vector <weight_dis> > >  >weights; // array for the weights MAIN ARRAY
   
   double nevents; // no of events calculated so far
   unsigned long nwrite;  // no of events after to write out the table
   double xsectsum; // total cross section - internal counter for test output

   basic_string<char> tablefilename; // The table file to write to

   amplitude_dis::integral_type itype; // Born, NLO etc.

  
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
  nj = 2U;

  //  number of the up and down type flavours
  nu = 2U;
  nd = 3U;
} 

void physcuts::setcuts()
{
  // energy of the incomong lepton and 
  // hadron in the laboratory frame
  El = 27.5;    // GeV
  Eh = 920.0;     // GeV
  
  //  Q^2 cuts
  Q2min = 100.0;    // GeV^2
  Q2max = 15000.0;  // GeV^2 

  //   xB cuts
  xmin = 0.0;
  xmax = 1.0;
  
  //   y cuts
  ymin = 0.2;
  ymax = 0.7;
}


void UserDIS::initfunc(unsigned int njet)
{

   unsigned int nj;
   unsigned int nu;
   unsigned int nd;
   inputfunc(nj,nu,nd);
   double El = thecuts->energy_lepton();
   double Eh = thecuts->energy_hadron();     
   double s = 4.0*El*Eh;
   double ymax = thecuts->getymax();

   xsectsum = 0;

   nQ2    = 18;                    // no of bins in Q2
   Q2high = new double[nQ2+1];      //----- array for Q2 boundaries
   npt    = new  int[nQ2];         // nQ2 bins in Q2 with each npt[iQ2] bins in pt

   //Define binning in Q2
   Q2high[0]=   100.0;         
   Q2high[1]=   120.0;         
   Q2high[2]=   150.0;         
   Q2high[3]=   170.0;         
   Q2high[4]=   200.0;         
   Q2high[5]=   240.0;         
   Q2high[6]=   300.0;         
   Q2high[7]=   400.0;         
   Q2high[8]=   600.0;         
   Q2high[9]=   800.0;         
   Q2high[10]= 1000.0;         
   Q2high[11]= 1200.0;
   Q2high[12]= 1600.0;
   Q2high[13]= 2000.0;         
   Q2high[14]= 2600.0;         
   Q2high[15]= 3200.0;         
   Q2high[16]= 4000.0;         
   Q2high[17]= 5000.0;         
   Q2high[18]=15000.0;         

 
   pthigh.resize(nQ2);
   for(int i=0;i<nQ2;i++){  
      //Define binning in pt
      npt[i]=21;
      pthigh[i].resize(npt[i]+1);
      //----- array for pt boundaries
      pthigh[i][0] = 5.0;
      pthigh[i][1] = 6.0;
      pthigh[i][2] = 7.0;
      pthigh[i][3] = 8.0;
      pthigh[i][4] = 9.0;
      pthigh[i][5] = 10.0;
      pthigh[i][6] = 11.0;
      pthigh[i][7] = 12.0;
      pthigh[i][8] = 13.0;
      pthigh[i][9] = 14.0;
      pthigh[i][10] = 16.0;
      pthigh[i][11] = 18.0;
      pthigh[i][12] = 20.0;
      pthigh[i][13] = 22.0;
      pthigh[i][14] = 24.0;
      pthigh[i][15] = 27.0;
      pthigh[i][16] = 30.0;
      pthigh[i][17] = 33.0;
      pthigh[i][18] = 36.0;
      pthigh[i][19] = 40.0;
      pthigh[i][20] = 44.0;
      pthigh[i][21] = 50.0;
   }
   // lowest pT value in sample
   ptlow = 5.0;     
    
   // Bins in x
   ntot = 75;

   // NLO scale variations
   nscalevar = 3;
   murscale.resize(nscalevar);
   mufscale.resize(nscalevar);
   
   // Array used for scale variations, mur factors correspond to default scale - variation is done later in usercode
   murscale[0] = 1.0; mufscale[0] = 1.0;
   murscale[1] = 1.0; mufscale[1] = 0.5;
   murscale[2] = 1.0; mufscale[2] = 2.0;

   xlimit.resize (nQ2);
   xsmallest.resize (nQ2);     // test: find smallest x values
   murval.resize (nQ2);
   mufval.resize (nQ2);
 
   weights.resize (nQ2);
   for( int j = 0; j < nQ2; j++) {
      xlimit[j].resize(npt[j]);
      xsmallest[j].resize(npt[j]);   // test: find smallest x values
      murval[j].resize (npt[j]);
      mufval[j].resize (npt[j]);
      weights[j].resize(npt[j]);
      double Q2 = Q2high[j];
      for( int k = 0; k < npt[j]; k++) {
         // Setup the weights array
         weights[j][k].resize(ntot);
         for( int l = 0; l < ntot; l++) {
           // scale variation
            weights[j][k][l].resize(nscalevar);
         }
         // - Setup the xlimit array
         double pt = pthigh[j][k];
         xlimit[j][k]= Q2/s/ymax+pow(pthigh[j][k+1],2)/s;
	 xsmallest[j][k]=0.999;
	 
         // - Setup the murval and mufval arrays
         murval[j][k]= (0.6*pt + 0.4*pthigh[j][k+1]); // Take as mean pt value 40% right of left bin boundary 
         mufval[j][k]= sqrt((0.7*Q2 + 0.3*Q2high[j+1])); 
      }
   }
   // Initialise event counters
   nevents = 0;
   // Set some defaults
   if (nwrite==0) nwrite = 2500000;
   if (tablefilename=="") tablefilename = "table.raw";

   // Say Hello
   cout<<"  "<<endl;
   cout<<"   *******************************************"<<endl;
   cout<<"    fastNLO    - initialiation"<<endl;
   cout<<"        table file "<<tablefilename<<endl;
   cout<<"        store table after "<<nwrite<<" events"<<endl;
   cout<<"        sqrt(s)= "<<sqrt(s)<<endl;
   cout<<"        No. x-bins: "<<ntot<<endl;
   cout<<"        No. Q2 bins: "<<nQ2<<endl;
   cout<<"        No. of pT bins in each Q2 bin region:"<<endl;
   for( int j = 0; j < nQ2; j++) {
      cout<<"         Q2 "<<j<<": "<<npt[j]<<endl;
   }
   cout<<"        No. of scale variations in NLO: "<<nscalevar<<endl;
   for( int j = 0; j < nscalevar; j++) {
     cout<<"          "<<j<<":   (mur/pT) "<<murscale[j]<<
	"      (muf/Q2) "<<mufscale[j]<<endl;
   }
   cout<<"   *******************************************"<<endl;
   cout<<"        "<<endl;

}

extern"C" double xalfaem_(double *);

double xalpha_em(double mq2) {
  return xalfaem_(&mq2);
}

void UserDIS::userfunc(const event_dis& p, const amplitude_dis& amp)
{
   typedef lorentzvector<double> _Lv;
   double alem = 0;
   bool firstjet = true;
   double x=0.;
   double Q2 =0;
   int Q2bin = -1;

   // LO or NLO?
   itype = amp.integral();

  
//   //----- copy the momenta and boost to the breit frame ----
   pbreit = p; lab_to_breit(pbreit);
  
  //----- do the cluster analysis-----

  //----- jet structure in breit frame -----
  pjb = jetclus(pbreit);
  int nj = pjb.upper(); 
  
  //----- jet structure in laboratory frame -----
  pjl = pjb; boost_back_to_lab(p);

  
  //----- cuts in the laboratory & breit frame -----
  const double emin = -1.0;
  const double emax = 2.5;        //  pseudo-rapidity cut

  //-------- analyze inclusive jets in jet loop -------
  for(int i = 1; i <= nj; i++) {
     double pt = pjb[i].perp();
     if(pt > ptlow){
        double eta = pjl[i].prapidity();
        if(eta > emin && eta < emax){
           if(firstjet){
              Q2 = -((p[-1] - p[-2]).mag2());
              alem = xalpha_em(Q2);
              //----- set the pdf function and the overall factors -----
              //              pdf.mode(pdf_cteq5m::nlo); pdf.loop(2);
              //              amp.pdf_and_qcd_coupling(pdf, alem*alem*389385.730);
              x = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);

              //              alem = xalpha_em(Q2);
              //------ determine Q2  bin 
              for( int j = 0; j < nQ2; j++) {
                 if (Q2 >= Q2high[j] && Q2 < Q2high[(j+1)]) Q2bin=j;
              }
              firstjet = false;
           }

           if (Q2bin >=0 ) {
              //------ determine pt bin 
              int ptbin  = -1;
              for( int j = 0; j < npt[Q2bin]; j++) {
                 if (pt >= pthigh[Q2bin][j] && pt < pthigh[Q2bin][(j+1)]) ptbin=j;
              }
              //---------- compute weight, fill NLOJET histos
              if ( ptbin>=0 ) {

                 // loop over scale variations for NLO
                 int scalevarmax;
                 if(itype==amplitude_dis::lo){
                    scalevarmax=1;
                 }else{
                    scalevarmax=nscalevar;
                 }
                 for(int scalevar=0; scalevar<scalevarmax;scalevar++){
                 
                    double mur2=murscale[scalevar]*murscale[scalevar]*murval[Q2bin][ptbin]*murval[Q2bin][ptbin];
                    double muf2=mufscale[scalevar]*mufscale[scalevar]*mufval[Q2bin][ptbin]*mufval[Q2bin][ptbin];
                    weight_dis wt = amp(mur2,muf2);


                    wt *= 1./x*alem*alem*389385.730;


                    // test if x_min in event is lower than x_limit
                    //      if yes -> make big warning!!! 
                    //      -> need to change x_limit calculation
                    if (x<xlimit[Q2bin][ptbin]){
                       printf("Warning: x (%f) < xlimit (%f) at pt=%f GeV Q2=%f GeV2 \n ",
                              x,xlimit[Q2bin][ptbin],pt,Q2);
                       exit(1);
                    }
              
                    // *****  test:  identify smallest x value in each pT/Q2 bin ********************
                    //               -> to optimize the x-limit values
//                     if (x<xsmallest[Q2bin][ptbin]*0.98){
//                        xsmallest[Q2bin][ptbin] = x;
//                        if (nevents>1000000){
//                           cout<<" "<<endl;
//                           cout<<">>>>>>>>>>> smaller x found in bin  "<<Q2bin<<"  "<<ptbin
//                               <<"   -  "<<x<<"  "<<xlimit[Q2bin][ptbin]<<"  :  "<<
//                              (x/xlimit[Q2bin][ptbin])<<"  >> "<<nevents<<endl;
//                           for( int j = 0; j < nQ2; j++) {
//                              for( int k = 0; k < npt[j]; k++) {
//                                 printf ("%f\n",xsmallest[j][k]);
//                              }
//                           }
//                 }
//                    }
                    // ****************  end: identify smallest x  *********************************

                    //---------- determine x_i position in grid
                    //--- determine fractional contributions for x-bin
                    double hx   = log10(x);
                    double hxlimit = log10(xlimit[Q2bin][ptbin]);
                    double hxone   = 0.0;
             
                    // define the x-bin number in the range  [0:ntot[
                    int nx = int(ntot *(hx-hxlimit)/(hxone-hxlimit));

                    //--- relative distance delta (= fractional contributions)
                    double delta  = (hxone-hxlimit)/ntot;
                    double hxi = hxlimit + double(nx)/double(ntot) * (hxone-hxlimit);
                    double deltam = (hx-hxi)/delta;
                    if(deltam>1.0  || deltam<0.0){
                       cout<<"  ->>> delta is off  : "<<deltam<<endl;
                    }

                    // transform to MWs pdf linear combinations
                    double pdfsigma = 1./3.*(-wt[1] + 4.*wt[2]);
                    double pdfdelta = 3.*(wt[1]-wt[2]);
                    wt[1] = pdfsigma;
                    wt[2] = pdfdelta;
 
                    // ** ----------------------------------------------------------
                    // **   now fill table
                    weights[Q2bin][ptbin][nx][scalevar]     += (1-deltam)*wt;
                    if(nx<ntot-1)
                       weights[Q2bin][ptbin][nx+1][scalevar]   += (deltam)*wt;


                 }
              }
           }
        }
     }
  }
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
   //  #define WRITE(n) table << setprecision(40) << n << endl

   int marker = 1234567890; //used to separate section in the table
   
   fstream table(tablefilename.c_str() ,ios::out); // open file
   // fstream table(tablefilename.c_str() ,ios::out|ios::binary); // open file
 
   // ireaction
   int ireaction = 1;
   WRITE(ireaction);

   double El = thecuts->energy_lepton();
   double Eh = thecuts->energy_hadron();     
   double s = 4.0*El*Eh;

   // Ecms
   s= sqrt(s);
   WRITE(s);

   //iproc
   int iproc = 1; // incl. jets
   WRITE(iproc);

   //ialgo
   int ialgo = 1;
   WRITE(ialgo);

   //JetResol1
   double JetResol1 = 1.;
   WRITE(JetResol1);

   //JetResol2
   double JetResol2 = 1.;
   WRITE(JetResol2);

   // relative order
   int nord = itype+1; // LO: itype=0, nord=1,  NLO: itype=1,nord=2
   WRITE(nord);

   // absolute order
   int npow = nord; // incl. jets and dijets in DIS -> LO=alpha_s^1,  NO=alpha_s^2 
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
   WRITE(nevents);

   //ntot
   WRITE(ntot);

   //ixscheme
   int ixscheme = 1;
   WRITE(ixscheme);

   //ipdfwgt
   int ipdfwgt = 0;
   WRITE(ipdfwgt);

   WRITE(marker);// ------------------END of block

   //NQ2
   WRITE(nQ2);
   for(int i=0;i<nQ2+1;i++){
      WRITE(Q2high[i]); //Q2[0] ... Q2[NQ2]
   }

   for(int i=0;i<nQ2;i++){
      WRITE(npt[i]); //Npt[0] ... Npt[NQ2-1] 
   }

   for(int i=0;i<nQ2;i++){
      for(int j=0;j<npt[i]+1;j++){
         WRITE(pthigh[i][j]); // all pt bins 
      }
   }
    
   WRITE(marker);// ------------------END of block
   
    for(int i=0;i<nQ2;i++){
      for(int j=0;j<npt[i];j++){
         WRITE(xlimit[i][j]); // all x limits 
      }
   }

   WRITE(marker);// ------------------END of block



   for(int i=0;i<nQ2;i++){
      for(int j=0;j<npt[i];j++){
         WRITE(murval[i][j]); // all murval 
      }
   }
   
   WRITE(marker);// ------------------END of block

   
   for(int i=0;i<nQ2;i++){
      for(int j=0;j<npt[i];j++){
         WRITE(mufval[i][j]); // all mufval 
      }
   }

   WRITE(marker);// ------------------END of block

   int scalevarmax;
   if(itype==amplitude_dis::nlo){
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

   for(int i=0;i<nQ2;i++){ // Q2
      for(int j=0;j<npt[i];j++){ // pt
         for(int k=0;k<ntot;k++){ // x
            // loop over scale variations for NLO
            for(int scalevar=0; scalevar<scalevarmax;scalevar++){
               for(int m=0;m<3;m++){     //subprocesses
                  WRITE(weights[i][j][k][scalevar][m]);
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
      printf ("No. events: %.2G writing table\n",nevents);
     writetable();
      printf("done.\n");
   }
}

void UserDIS::phys_output(const std::basic_string<char>& __file_name, 
                            unsigned long __save, bool __txt) 
{
   // First let NLOJET do the work
   user_dis::phys_output(__file_name,__save,__txt);
   tablefilename = __file_name +".raw";
   nwrite = __save;
}
