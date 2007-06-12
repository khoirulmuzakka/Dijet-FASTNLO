// Nlojet++ standalone Cross Check of Fast Computation

// 2005/09/14 - KR, adapted from Markus' file

//------ DON'T TOUCH THIS PART! ------
#include <phasespace.h>
#include <process.h>
#include <jetfunc.h>
#include <qcdlib.h>

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
#include "rsep-et-07.h"
#include "cteq6.h"

class UserHHC : public user_hhc {
public:
  // init and user function
  void initfunc(unsigned int);
  void userfunc(const event_hhc&, const amplitude_hhc&);
  
private:
  // define binning
  int nrap;                       // no. of rapidity bins 
  double *raphigh;                // array for rapidity boundaries
  int *npt;                       // no. of pT bins in each y range
  vector< vector<double> >pthigh; // array for pT boundaries
  double ptlow;                   // lowest pt considered 
  vector< vector<double> >ptmid;  // array for pT mid points
  
  // pdf and jet algorithm
  pdf_cteq6 pdf;
  rsep_et_07 jetclus;
  
  // the jet structure 
  bounded_vector<lorentzvector<double> > pj; 
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
  s = 3240000;    // Run I
  //  s = 3841600;
  
  //  number of the up and down type flavours
  nu = 2U;
  nd = 3U;
} 



void UserHHC::initfunc(unsigned int) {
  // initialize binning
  nrap    = 6;                  // no. of bins in rapidity (+1 for CDF)
  raphigh = new double[nrap+1]; // array for rapidity boundaries
  npt     = new  int[nrap];     // nrap bins in rapidity with each npt[irap] bins in pt

  // initialize binning in rapidity
  raphigh[0] = 0.0;        // bins are:  0.0 - 0.5 - 1.0 - 1.5 - 2.0 - 3.0
  raphigh[1] = 0.5;        // 
  raphigh[2] = 1.0;        // 
  raphigh[3] = 1.5;        // 
  raphigh[4] = 2.0;        // 
  raphigh[5] = 3.0;        // bins are:  0.0 - 0.5 - 1.0 - 1.5 - 2.0 - 3.0
  raphigh[6] = 0.7;        // specific for CDF (0.1-0.7)

  // initialize binning in pt (90 + 33 bins in total)
  npt[0]=24;
  npt[1]=24;
  npt[2]=19;
  npt[3]=15;
  npt[4]=8;
  npt[5]=33;               // CDF

  // lowest pT value in sample
  ptlow = 40.5;            // 60.0 for D0 - 40.5 for CDF
  
  // arrays for pt boundaries and mid points:
  ptmid.resize(nrap); 
  pthigh.resize(nrap);
  for(int i=0;i<nrap;i++){
    ptmid[i].resize(npt[i]+1);
    pthigh[i].resize(npt[i]+1);
  }
   
  // pt boundaries
  // D0 rapidity bin 1
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
  // D0 rapidity bin 2
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
  // D0 rapidity bin 3
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
  // D0 rapidity bin 4
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
  // D0 rapidity bin 5
  pthigh[4][0] = 60.0;
  pthigh[4][1] = 70.0;
  pthigh[4][2] = 80.0;
  pthigh[4][3] = 90.0;
  pthigh[4][4] =100.0;
  pthigh[4][5] =110.0;
  pthigh[4][6] =130.0;
  pthigh[4][7] =160.0;
  pthigh[4][8] =210.0;
  // bins from CDF ansatz, rapidity bin 0.1 - 0.7
  pthigh[5][0] =  40.5;
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
  
  // prepare histogram naming
  char histname[256]; // histogram name for phys
  int ires;           // temporary string manip variables
  char buffer[128];

  // Histogram Definitions (phys)
  for (int jrap=0;jrap<nrap;jrap++) {
    //    cout << "y bin no. = " << jrap << endl;
    point_type ptbins[npt[jrap]];
    for (int ipt=0;ipt<npt[jrap];ipt++) {
      //      cout << "pt bin no. = " << ipt << endl;
      ptmid[jrap][ipt] = (pthigh[jrap][ipt]+pthigh[jrap][ipt+1])/2;
      ptbins[ipt].xmin = pthigh[jrap][ipt];
      ptbins[ipt].xmid = ptmid[jrap][ipt];
      ptbins[ipt].xmax = pthigh[jrap][ipt+1];
    }
    strcpy(histname,"Scenario fnt1001: ");
    ires = sprintf(buffer,"Rapidity bin no. %i : %.1f < |y| < %.1f",
		   jrap,raphigh[jrap],raphigh[jrap+1]);
    if (jrap < 5) {
      strcat(histname,buffer);
    } else {
      strcat(histname,"CDF rapidity bin   : 0.1 < |y| < 0.7");
    }
    cout << "Booking histogram: " << histname << endl;
    phys(jrap+1,histname,npt[jrap],ptbins,0.25);
  }
}



void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp) {
  typedef lorentzvector<double> _Lv;
  
  // define logic array for rapidity acceptance  
  bool lrap[nrap];
  
  // do the jet analysis
  pj = jetclus(p);
  unsigned int nj = pj.upper(); 
  
  // cuts (D0)
  double s = 2.0*(p[hadron(0)]*p[hadron(-1)]);
  double ymin = 0.0, ymax = 3.0, ptcut = 60.0; 
  double pt_min = sqrt(s), pt_max = 0.0, rap, pt, ht = 0.0;
  
  unsigned int jetnum = 0;
  
  for (int irap=0; irap < nrap; irap++) {
    lrap[irap] = false;
    //    cout << "1. irap, lrap: " << irap << ", " << lrap[irap] << endl;
  }

  for (unsigned int i = 1; i <= nj; i++) {
    pt = pj[i].perp(); rap = abs(pj[i].rapidity());
    if (rap < ymax && rap >= ymin && pt > ptcut) {
      if (pt < pt_min) pt_min = pt;
      if (pt > pt_max) pt_max = pt;
      jetnum++; ht += pt;
    }
  }
  
  // Type of the contribution & pdf
  amplitude_hhc::contrib_type itype = amp.contrib();
  amp.pdf_and_qcd_coupling(pdf, 389385.730);  
  pdf.mode(pdf_cteq6::nlo); pdf.loop(2);
  
  // Scenario fnt1001 jets standalone histograms
  for (unsigned int i = 1; i <= nj; i++) {
    pt = pj[i].perp(); rap = abs(pj[i].prapidity());
    // D0 bins
    for (int jrap = 0; jrap < nrap-1; jrap++) {
      lrap[jrap] = false;
      if ( raphigh[jrap] <= rap && rap < raphigh[jrap+1] &&
	   pt > 60.0 ) {
	lrap[jrap] = true;
      }
      //      cout << "D0 jrap, lrap, pt: " << jrap << ", " << lrap[jrap] << ", " << pt << endl;
    }
    // CDF bin
    lrap[nrap-1] = false;
    if ( 0.1 <= rap && rap < raphigh[nrap] &&
	 pt > 40.5 ) {
      lrap[nrap-1] = true;
    }
    //    cout << "CDF nrap-1, lrap, pt: " << nrap-1  << ", " << lrap[nrap-1] << ", " << pt << endl;
    
    for (int jrap = 0; jrap < nrap; jrap++) {
      //      cout << "rap, jrap, lrap: " << rap << ", " << jrap << ", " << lrap[jrap] << endl; 
      if (lrap[jrap]) {
	double mu2 = pt*pt/4.0;
	weight_hhc wt = amp(mu2,mu2);
	// Divide by delta y and multiply by 1000. to convert to pb
	// A factor of two is there due to |y|!
	if ( jrap < nrap-1 ) {
	  // D0 bins
	  wt = 1000. * wt / (raphigh[jrap+1]-raphigh[jrap]) / 2.; 
	} else {
	  // CDF bin
	  wt = 1000. * wt / 1.2;
	} 
	physfilld(jrap+1, pt, wt);
      }
    }
  }
}
