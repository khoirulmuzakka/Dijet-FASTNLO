//
// UsefulNlojetTools
// 
// Collection of useful functions for the fastNLO interface
// to NLOJET++.
//

#include <cstdio>
#include <algorithm> //c++98
#include <utility>   //c++11
#include "fastnlotk/fastNLOCreate.h"

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
      double ecms = 0;
      psinput(NULL,ecms);
      return sqrt(ecms);
   }
   
   //_______________________________________________________________________
   unsigned int GetOrderOfRun(const std::basic_string<char>& __file_name){ 
      // --- determine whether we are running LO or NLO
      // and add order of expansion to leading-order .
      const char* const file = __file_name.c_str(); 
      
      if(strstr(file,"born")!=NULL){
	 //IsNLO = false;
	 printf("fastNLO: This is a LO run!\n");
	 return GetNj();
      }else{
	 if(strstr(file,"nlo")!=NULL){
	    //IsNLO = true;
	    printf("fastNLO: This is a NLO run!\n");
	    return GetNj() + 1;
	 }else{
	    // it is 'nlojet'
	    printf("fastNLO: ERROR! This module can only be run at Born level or at NLO.\n");
	    exit(1);
	 }
      }
   }

   //_______________________________________________________________________
   int NlojetToFastnloIDHHC(int id) {
      if ( id == 1 ) return 5;
      else if ( id == 2 ) return 6;
      else if ( id == 3 ) return 1;
      else if ( id == 4 ) return 2;
      else if ( id == 5 ) return 3;
      else if ( id == 6 ) return 4;
      else return id;
   }
   int FastnloIdToNlojetIdHHC(int id) {
      if ( id == 1 ) return 3;
      else if ( id == 2 ) return 4;
      else if ( id == 3 ) return 5;
      else if ( id == 4 ) return 6;
      else if ( id == 5 ) return 1;
      else if ( id == 6 ) return 2;
      else return id;
   }


   //_______________________________________________________________________
   vector<fnloEvent> GetFlexibleScaleNlojetContribHHC(const event_hhc& p , const amplitude_hhc& amp , nlo::pdf_and_coupling_hhc& dummypdf){
      static const int nSubproc = 7 ; // nlojet must have 7 subprocesses.
      static const double dummyMu2 = 100.*100.;
      
      // make events
      vector<fnloEvent> ev(nSubproc);

      // fill relevant event quantities
      double x1 = p[-1].Z()/p[hadron(-1)].Z();
      double x2 = p[0].Z()/p[hadron(0)].Z();
      for ( int p = 0 ; p<nSubproc ; p++ )  {
	 ev[p].SetProcessId( p );
	 ev[p].SetX1( x1 );
	 ev[p].SetX2( x2 );
      }

      // weights
      double weights[7][7]; // weights[amp_i][proc]
      for ( int kk = 0 ; kk<7 ; kk ++ )	 for ( int p = 0 ; p<nSubproc ; p ++ ) weights[kk][p] = 0; 
      
      // access perturbative coefficients
      nlo::amplitude_hhc::contrib_type itype = amp.contrib();
      nlo::weight_hhc cPDF = dummypdf.pdf(x1,x2,dummyMu2,2,3); // 1/x1/x2
      static const double coef = 389385730.;

      for ( int kk = 0 ; kk<7 ; kk ++ ) {
	 for ( int fid = 0 ; fid<nSubproc ; fid ++ ) {
	    int nid = FastnloIdToNlojetIdHHC(fid);
	    weights[kk][fid] = amp._M_fini.amp[kk][nid]*coef*cPDF[nid];
	 }
      }

      // todo: calculate wt and wtorg from 'weights'
      nlo::weight_hhc wtorg = amp(dummypdf,dummyMu2,dummyMu2, 1.);
      nlo::weight_hhc wt = wtorg;
      // - rearrange subprocesses
      for ( int fid = 0 ; fid<7 ; fid ++ ) {
	 int nid = FastnloIdToNlojetIdHHC(fid);
	 wt[fid] = wtorg[nid];
      }
      wt *= coef;
   
      if(x2>x1){
	 // swap subprocesses 6,7
	 swap(wt[5],wt[6]);
	 for ( int kk = 0 ; kk<7 ; kk ++ )
	    swap(weights[kk][5],weights[kk][6]);
      }	       

      for(int p=0 ; p<nSubproc ; p++){
	 // decompose nlojet-event
	 if(itype == nlo::amplitude_hhc::fini) {
	    if (amp._M_fini.mode==0) { //finix1
	       ev[p].SetWeight_MuIndependent( weights[0][p] );
	       ev[p].SetWeight_log_muf( weights[3][p] );		  
	    }
	    else if (amp._M_fini.mode==1) { //finix2
	       ev[p].SetWeight_MuIndependent( weights[1][p] );
	       ev[p].SetWeight_log_muf( weights[4][p] );		  
	    }
	    else if(amp._M_fini.mode==2){ //fini1
	       ev[p].SetWeight_MuIndependent( weights[2][p] );
	       ev[p].SetWeight_log_muf( weights[5][p] );		  
	       ev[p].SetWeight_log_mur( weights[6][p] );		  
	    }
	 }
	 else { // no fini contribution
	    ev[p].AddWeight_MuIndependent( wt[p] );
	 }
      }

      return ev;
   }

}
