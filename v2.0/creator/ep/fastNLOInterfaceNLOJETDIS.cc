//
// UsefulNlojetTools
// 
// Collection of useful functions for the fastNLO interface
// to NLOJET++.
//

#include <cstdio>
#include <algorithm> //c++98
#include <utility>   //c++11
#include <fastnlotk/fastNLOCreate.h>

namespace UsefulNlojetTools {

   //_______________________________________________________________________
   unsigned int GetNj(){
      // get order of leading-order
      // inputfunc should be defined in header of nlojet-module beforehand.
      unsigned int nj = 0, nu = 0 ,nd = 0;
      inputfunc(nj,nu,nd);
      return nj;
   }  

   double GetEcms(){
      // get center of mass energy
      //psinput(NULL,ecms); // hhc
      double devnull, el, eh;
      psinput(NULL, el, eh, devnull, devnull, devnull, devnull, devnull, devnull);
      double ecms = 4*el*eh;
      return sqrt(ecms);
   }
   //_______________________________________________________________________
   unsigned int GetLoOrder(){
      // get order of leading-order
      // inputfunc should be defined in header of nlojet-module beforehand.
      unsigned int nj = 0, nu = 0 ,nd = 0;
      inputfunc(nj,nu,nd);
      return nj - 1;
   }  
   
   //_______________________________________________________________________
   unsigned int GetOrderOfRun(const std::basic_string<char>& __file_name){ 
      // --- determine whether we are running LO or NLO
      // and add order of expansion to leading-order .
      const char* const file = __file_name.c_str(); 
      static int ord = -1;
      
      if ( ord == -1 ) {
	 if(strstr(file,"born")!=NULL){
	    //IsNLO = false;
	    printf("\nfastNLO: This is a LO run!\n");
	    ord = GetLoOrder();
	 }else{
	    if(strstr(file,"nlo")!=NULL){
	       //IsNLO = true;
	       printf("\nfastNLO: This is a NLO run!\n");
	       ord = GetLoOrder();
	       ord += 1;
	    }else{
	       // it is 'nlojet'
	       printf("fastNLO: ERROR! This module can only be run at Born level or at NLO.\n");
	       exit(1);
	    }
	 }
	 cout<<"\nfastNLO: Found order of calculation to be: ord="<<ord<<endl;
	 return ord;
      }
      else return ord;
   }

   //_______________________________________________________________________
   unsigned int GetNSubproc() {
      static int nSubproc = -1;
      if ( nSubproc == -1 ) {
	 int nord = GetOrderOfRun("bla");
	 int loord = GetLoOrder() ;
	 if ( nord == loord ) nSubproc = 2 ;
	 else nSubproc = 3;
      }
      return nSubproc;
   }

   //_______________________________________________________________________
   vector<fnloEvent> GetFlexibleScaleNlojetContribDIS(const event_dis& p , const amplitude_dis& amp , nlo::pdf_and_coupling_dis& dummypdf, double alem){
      const int nSubproc = GetNSubproc() ; // nlojet must have 3 subprocesses.
      static const double dummyMu2 = 100.*100.;
      
      // make events
      vector<fnloEvent> ev(nSubproc);

      // fill relevant event quantities
      double x = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);
      //double Q2 = -((p[-1] - p[-2]).mag2());
      for ( int p = 0 ; p<nSubproc ; p++ )  {
	 ev[p].SetProcessId( p );
	 ev[p].SetX1( x );
      }

      
      double coef = 389385730.*alem*alem;

      // ---- calcualte weights for fini contributions ---- //
      nlo::weight_dis wt = amp(dummypdf,dummyMu2*dummyMu2,dummyMu2*dummyMu2,coef); // M1*M1=Q2 is a 'dummy'-scales
      // - order is 0:Delta, 1:Gluon, 2:Sigma
      double pdfsigma = 1.0/3.0*(-wt[1] + 4.*wt[2]);
      double pdfdelta = 3.0*(wt[1]-wt[2]);
      wt[1] = wt[0];
      wt[0] = pdfdelta;
      wt[2] = pdfsigma;
      //wt *= 389385730.;
      //if(IXsectUnits!=12)  wt *= pow(10.,(IXsectUnits-12)) ;
      
      nlo::weight_dis cPDF = dummypdf.pdf(x,dummyMu2*dummyMu2,2,3); // M1*M1 is a dummy scale
      nlo::amplitude_dis::contrib_type itype = amp.contrib();

      // weights
      double weights[5][nSubproc]; // weights[amp_i][proc]
      if ( nSubproc > 2 ) {
	 for ( int kk = 0 ; kk<5 ; kk ++ )	 for ( int p = 0 ; p<nSubproc ; p ++ ) weights[kk][p] = 0; 
	 if(itype == nlo::amplitude_dis::fini) {
	    for ( int kk = 0 ; kk<5 ; kk ++ ){
	       double wSigma = 1.0/3.0*(-amp._M_fini.amp[kk][1]*coef*cPDF[1] + 4.*amp._M_fini.amp[kk][2]*coef*cPDF[2]);
	       double wDelta = 3.0*(amp._M_fini.amp[kk][1]*coef*cPDF[1] - amp._M_fini.amp[kk][2]*coef*cPDF[2]);
	       double wGluon = amp._M_fini.amp[kk][0]*coef*cPDF[0];
	       weights[kk][0] = wDelta;
	       weights[kk][1] = wGluon;
	       weights[kk][2] = wSigma;
	    }
	 }
      }
      
      for(int p=0 ; p<nSubproc ; p++){
	 // decompose nlojet-event
	 if(itype == nlo::amplitude_dis::fini) {
	    if (amp._M_fini.mode==0) { //finix1
	       ev[p].SetWeight_MuIndependent( weights[0][p] );
 	       ev[p].SetWeight_log_muf( weights[2][p] );		  
	    }
	    else if (amp._M_fini.mode==1) { //finix2
 	       ev[p].SetWeight_MuIndependent( weights[1][p] );
 	       ev[p].SetWeight_log_muf( weights[3][p] );		  
 	       ev[p].SetWeight_log_mur( weights[4][p] );		  
	    }
	 }
	 else { // no fini contribution
	    ev[p].AddWeight_MuIndependent( wt[p] );
	 }
      }

      return ev;
   }

}
