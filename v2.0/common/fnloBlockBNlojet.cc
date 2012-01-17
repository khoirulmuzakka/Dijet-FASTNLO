#include <cstdlib>
#include <iostream>

#include "fnloBlockBNlojet.h"

//
// note (Oct 8, 2010) - MW:
// the two routines FillEventDIS, and FillEventHHC are running
// the other routines (DIS2scale and photoproduction) still need some work
//

fnloBlockBNlojet::fnloBlockBNlojet(fnloBlockA1 *blocka1, fnloBlockA2 *blocka2) :fnloBlockB(blocka1,blocka2) {
   _S_gauleg(20, _M_xb, _M_wb);
   counter=0;
   for ( int i = 0 ; i < 3 ; i++ ) Fct_MuR_Ref[i] = NULL;
   for ( int i = 0 ; i < 3 ; i++ ) Fct_MuF_Ref[i] = NULL;
}


void fnloBlockBNlojet::FillEventDIS(int ObsBin, double x, double scale1, const nlo::amplitude_dis& amp, nlo::pdf_and_coupling_dis& pdf, double prefactor){
   fnloBlockA2 *A2 =  BlockA2;

   // ---
   // --- Warm-Up Run to identify extreme x,mu values
   // ---
   // KR: Add file output for later automatic read in
   if (IWarmUp == 1) {
      WarmUp( ObsBin , x , scale1 , 0, "xlim" , "mu" );
      return;
   }

   // --- select interpolation kernel for x and for mu 
   //             1:CatmulRom   2:Lagrangian
   const int ikernx = 1;     
   const int ikernmu = 1;    

   if(this->IRef>0){
      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){
         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         nlo::weight_dis wt = amp(pdf,mu2,mu2,prefactor);

	 // --- rearrange elements of wt for reference table
         //     store gluon in [1] - add quark contributions in [0]
	 //     this is _not_ consistent with final format, but the best solution
         //     at least the gluon (wt[1]) is at the same position
         //     and at order(alphas) also wt[0] has the same content
	 double tmp = wt[1]+wt[2];
	 wt[1] = wt[0];
	 wt[0] = tmp;
	 wt[2] = 0;

         wt *= 389385730.;
         if(IXsectUnits!=12){
            wt *= pow(10.,(IXsectUnits-12)) ;
         }
         
         for(int proc=0;proc<NSubproc;proc++){
            SigmaTilde[ObsBin][scalevar][0][0][proc] += wt[proc];
         }
     }
   }else{

      if (x<XNode1[ObsBin][0]){
         printf("fnloBlockBNlojet::FillEventDIS: find: x (%f) smaller than lowest x-node (%f) for bin #%d .\n",
                x,XNode1[ObsBin][0],ObsBin);
      }

      // **********  determine x_ij position in grid  ************
      // --- determine fractional contributions
      double hx = log10(x);
      double hxone = 0.0;

      // --- define the x-node number in the range: 0 <= nxnode < ntot
      double hxlimit = Hxlim1[ObsBin];
      int nxnode = int(Nxtot1[ObsBin] *(hx-hxlimit)/(hxone-hxlimit));
      if (nxnode < 0) nxnode = 0;  // move into available range (allow to extrapolate for x<xlimit)
      int nx = 0;   // - variable for final node (after modification in kernel)

      // --- relative distance in h(x): deltam
      double delta  = (hxone-hxlimit)/Nxtot1[ObsBin];
      double hxi = hxlimit+double(nxnode)/double(Nxtot1[ObsBin])*(hxone-hxlimit);
      double deltam = (hx-hxi)/delta;

      // --- get x-interpolation kernel and updated x-node position: 1 <= nx < ntot-1
      vector<double> cefm(4) ; 
      Interpol(nxnode, (Nxtot1[ObsBin]-1), deltam, ikernx, nx, &cefm);

      // --- PDF reweighting - compute weights, modify cefm[.] 
      //     but only those within grid, there are no nodes at x=1
      double pdfwgtm = PDFwgt(x);
      for( int i1 = 0; i1 < 4; i1++) {
        if ((nx-1+i1) >= 0 && (nx-1+i1) < Nxtot1[ObsBin] ) {
          cefm[i1] *= pdfwgtm/PDFwgt(XNode1[ObsBin][nx-1+i1]);
	}
      }

      // --- loop over all scale variations
      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){

	 // --- compute renormalization=factorization scale squared 
         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;

	 // --- get coefficients - convert to our 3 subprocesses         
         nlo::weight_dis wt = amp(pdf,mu2,mu2,prefactor);
         double pdfsigma = 1.0/3.0*(-wt[1] + 4.*wt[2]);
         double pdfdelta = 3.0*(wt[1]-wt[2]);
	 // - order is 0:Delta, 1:Gluon, 2:Sigma
	 wt[1] = wt[0];
	 wt[0] = pdfdelta;
	 wt[2] = pdfsigma;

         wt *= 389385730.;
         if(IXsectUnits!=12){
            wt *= pow(10.,(IXsectUnits-12)) ;
         }

         // --- define the scale-bin number in the range:  0 <= scalenode < nscalebin-2
         int scalenode = Nscalenode[0]-2;  // --- initialize with largest possible value
	 int nscale = 0;                   // --- variable for final scale node
	 double mu0scale = 0.25;           // --- parameter in transformation function H(mu)
         vector<double> cefscale(4);
	 cefscale[0] = 0.0;
	 cefscale[1] = 0.5;
	 cefscale[2] = 0.5;
	 cefscale[3] = 0.0;

         if(Nscalenode[0]>1){
	    // --- find scale position in range:  0 <= scalenode < nscalenode-1
            for(int i=1;i<Nscalenode[0];i++){
	      //printf("test node %d %d %f %f \n",ObsBin,i,scale1,ScaleNode[ObsBin][0][scalevar][i]);
               if (ScaleFac[0][scalevar]*scale1<ScaleNode[ObsBin][0][scalevar][i]){
                  scalenode=i-1;
                  break;
               }
            }

	    // --- relative distance delta - in function H(mu)
	    double deltascale = (log(log(ScaleFac[0][scalevar]*scale1/mu0scale)) - HScaleNode[ObsBin][0][scalevar][scalenode])/
               (HScaleNode[ObsBin][0][scalevar][scalenode+1]-HScaleNode[ObsBin][0][scalevar][scalenode]);

	    // --- get scale interpolation kernel and updated scalenode position: 1 <= nmu < ntot-2
	    Interpol(scalenode, (Nscalenode[0]-2), deltascale, ikernmu, nscale, &cefscale);
         }


         // --- loop over all 4 scale nodes that receive contributions
         for(int i3 = 0; i3 < 4; i3++){
            int is = nscale + (i3-1);     // --- the target scale index

	    // --- check: all elements should end up within grid
            if (is < 0){
	       printf(" is<0     %d %d %f \n",nscale,i3,cefscale[i3]);
	       exit(1);
            } else if (is > Nscalenode[0]-1){
	       printf(" is>max   %d %d %f \n",nscale,i3,cefscale[i3]);
	       exit(1);
            }
            // --- loop over all 4 x-nodes that receive contributions
            for( int i1 = 0; i1 < 4; i1++) {           
               int im = nx + (i1-1);     // --- the target x index
	       if (im < Nxtot1[ObsBin]) {
		  for(int proc=0;proc<NSubproc;proc++){
                  //                  printf("%d %d %d %d %d %f\n",ObsBin,scalevar,is,im,proc,cefscale[i3]);
		    SigmaTilde[ObsBin][scalevar][is][im][proc] += cefm[i1] * cefscale[i3] * wt[proc];
		  }
               }
            }
         }
      }
   }
}



// ___________________________________________________________________________________________________ //



void fnloBlockBNlojet::FillEventDISMuVar(int ObsBin, double x, double M1, double M2, const nlo::amplitude_dis& amp, nlo::pdf_and_coupling_dis& dummypdf, nlo::pdf_and_coupling_dis& realpdf, double prefactor){
   // ------------------------------------------------------------------------------------------------------- //
   // 
   // Method for use in nlojet++-module for filling FastNLO sigmatilde-tables for v2.1 tables in DIS
   //
   //    ObsBin:		bin number of this event in the measurement binning
   //    x			x value of this event
   //    M1, M2			values M1 and M2 of your scales (scale1 and scale2). Usually 'Q' and 'p_T'
   //    amp			nlojet++ amplitude object
   //    dummypdf		dummypdf as nlojet++-pdf-object
   //    realpdf		like dummypdf but will evaluate real pdf table for reference tables
   //    prefactor		sth like:  w_PS * 1/x * alpha_em^2 h-bar^2 c^2
   //
   // ------------------------------------------------------------------------------------------------------- //

   fnloBlockA2 *A2 =  BlockA2;
   
   // ----- Set Constants ---- //
   // --- select interpolation kernel for x and for mu 
   //             kCatmulRom   kLagrangian
   const EInterpolKernel ekernx = kCatmulRom;     

   // ------------- interpolation kernel ------------ //
   const EInterpolKernel ekernmu = kCatmulRom;    
   const double mu0scale = 0.25;           // --- parameter in transformation function H(mu)
   
   //NEventsBin[ObsBin]++;

   
   // ---------------------------------------------------------------------------------------- //
   // ------------------------ Warm-Up Run to identify extreme x,mu values ------------------- //
   // ---------------------------------------------------------------------------------------- //

   if (IWarmUp == 1) {
      WarmUp( ObsBin , x , M1 , M2 , "xlim" , "scale1" , "scale2");
      return;
   }

   
   // ---------------------------------------------------------------------------------------- //
   // --------------------------------- reference tables ------------------------------------- //
   // ---------------------------------------------------------------------------------------- //
   FillMuVarReferenceTables( ObsBin, M1, M2, amp, realpdf, prefactor);



   // ---- warning, if your event is outside the 'x'-range of your x-grid ---- //
   if (x<XNode1[ObsBin][0]){
      printf("fnloBlockBNlojet::FillEventDISMuVar: find: x (%f) smaller than lowest x-node (%f) for bin #%d .\n",
	     x,XNode1[ObsBin][0],ObsBin);
   }


   // ---------------------------------------------------------------------------------------- //
   //  ------------------- determine x_ij position in grid ----------------------------------- //
   // ---------------------------------------------------------------------------------------- //
   // --- determine fractional contributions
   double hx = (this->*Fct_H_XNode)(x);
   double hxone = 0.0;
   
   // --- define the x-node number in the range: 0 <= nxnode < ntot
   double hxlimit = Hxlim1[ObsBin];
   int nxnode = int(Nxtot1[ObsBin] *(hx-hxlimit)/(hxone-hxlimit));
   if (nxnode < 0) nxnode = 0;  // move into available range (allow to extrapolate for x<xlimit)
   
   // --- relative distance in h(x): deltam
   double delta  = (hxone-hxlimit)/Nxtot1[ObsBin];
   double hxi = hxlimit+double(nxnode)/double(Nxtot1[ObsBin])*(hxone-hxlimit);
   double deltam = (hx-hxi)/delta;
   
   // --- get x-interpolation kernel and updated x-node position: 1 <= nx < ntot-1
   int nx = 0;   // - variable for final node (after modification in kernel)
   vector < double > cefm = Interpol(nxnode, (Nxtot1[ObsBin]-1), deltam, ekernx, nx);

   // --- PDF reweighting - compute weights, modify cefm[.] 
   //     but only those within grid, there are no nodes at x=1
   double pdfwgtm = PDFwgt(x);
   for( int i1 = 0; i1 < 4; i1++) {
      if ((nx-1+i1) >= 0 && (nx-1+i1) < Nxtot1[ObsBin] ) {
	 cefm[i1] *= pdfwgtm/PDFwgt(XNode1[ObsBin][nx-1+i1]);
      }
   }
   


   if( NscalenodeScale1<=2 || NscalenodeScale2<=2){ 
      printf("fnloBlockBNlojet::FillEventDISMuVar(). Error. Sorry, but you need some more scale nodes!\n");exit(1);
   }
   // ---------------------------------------------------------------------------------------- //
   // ------------------------------------ scale1 nodes -------------------------------------- //
   // ---------------------------------------------------------------------------------------- //
   // --- find scale position in range:  0 <= scalenode1 < nscalenode1-1
   int scalenode1 = NscalenodeScale1-2;  // --- initialize with largest possible value
   for( int iNode=1 ; iNode<NscalenodeScale1 ; iNode++ ){
      if ( M1 < ScaleNode1[ObsBin][iNode]){
	 scalenode1=iNode-1;
	 break;
      }
   }
      
   // --- relative distance delta - in function H(mu)
   // deltascale (Interpol(.,.,.delta,.): relative distance of value to node 'nnode'
   double deltascale1 = 
      ( (this->*Fct_H_Scale[0])(M1) - HScaleNode1[ObsBin][scalenode1] ) /
      (HScaleNode1[ObsBin][scalenode1+1] - HScaleNode1[ObsBin][scalenode1]);
      
   // --- get scale interpolation kernel and updated scalenode position: 1 <= nmu < ntot-2
   int nscale1 = 0;                      // --- variable for final scale2-node
   vector<double> cefscale1 = Interpol(scalenode1, NscalenodeScale1-2, deltascale1, ekernmu, nscale1);
  
      
   // ---------------------------------------------------------------------------------------- //
   // ------------------------------------ scale2 nodes -------------------------------------- //
   // ---------------------------------------------------------------------------------------- //
   // --- find scale position in range:  0 <= scalenode2 < nscalenode2-1
   int scalenode2 = NscalenodeScale2-2;  // --- initialize with largest possible value
   for( int iNode=1 ; iNode<NscalenodeScale2 ; iNode++ ){
      if ( M2 < ScaleNode2[ObsBin][iNode]){
	 scalenode2=iNode-1;
	 break;
      }
   }
	
   // --- relative distance delta - in function H(mu)
   // deltascale (Interpol(.,.,.delta,.): relative distance of value to node 'nnode'
   double deltascale2 = 
      ( (this->*Fct_H_Scale[1])(M2) - HScaleNode2[ObsBin][scalenode2] ) /
      (HScaleNode2[ObsBin][scalenode2+1] - HScaleNode2[ObsBin][scalenode2]);
      
   // --- get scale interpolation kernel and updated scalenode position: 1 <= nmu < ntot-2
   int nscale2 = 0;                      // --- variable for final scale2-node
   vector<double> cefscale2 = Interpol(scalenode2, NscalenodeScale2-2, deltascale2, ekernmu, nscale2);
   

   
   // ---------------------------------------------------------------------------------------- //
   // ----------------- calculate the matrix element and get weights ------------------------- //
   // ---------------------------------------------------------------------------------------- //
   // --- get coefficients - convert to our 3 subprocesses ---- //
   nlo::weight_dis wt = amp(dummypdf,M1*M1,M1*M1,prefactor); // M1*M1=Q2 is a 'dummy'-scales
   double pdfsigma = 1.0/3.0*(-wt[1] + 4.*wt[2]);
   double pdfdelta = 3.0*(wt[1]-wt[2]);
   // - order is 0:Delta, 1:Gluon, 2:Sigma
   wt[1] = wt[0];
   wt[0] = pdfdelta;
   wt[2] = pdfsigma;
   
   wt *= 389385730.;
   if(IXsectUnits!=12)  wt *= pow(10.,(IXsectUnits-12)) ;
   
   // ---- calcualte weights for fini contributions ---- //
   nlo::amplitude_dis::contrib_type itype = amp.contrib();
   double coef = 389385730.*prefactor;
   if(IXsectUnits!=12)  coef *= pow(10.,(IXsectUnits-12)) ;
   
   nlo::weight_dis cPDF = dummypdf.pdf(x,M1*M1,2,3);	// M1*M1 is a dummy scale
   
   //double weights[5][3];		// would be nice to make this a nlo::weight_dis[5] object
   nlo::weight_dis weights[5];
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
   

   
   // ---------------------------------------------------------------------------------------- //
   // -------------------------------- update sigmatilde tables ------------------------------ //
   // ---------------------------------------------------------------------------------------- //
   
   // --- loop over all 4x 4x4  nodes that receive contributions
   for(int iNode1 = 0; iNode1 < 4; iNode1++ ){
      int idxQ = nscale1 + (iNode1-1);     // --- the target scale1 index
      for(int iNode2 = 0; iNode2 < 4; iNode2++ ){
	 int idxPt = nscale2 + (iNode2-1);     // --- the target scale2 index
	 
	 // --- check: all elements should end up within grid ---- //
	 if (idxQ < 0){
	    printf(" idxQ<0     %d %d %f \n",nscale1,iNode1,cefscale1[iNode1]);
	    exit(1);
	 } else if (idxQ > NscalenodeScale1-1){
	    printf(" idxQ>max   %d %d %f \n",nscale1,iNode1,cefscale1[iNode1]);
	    exit(1);
	 }
	 if (idxPt < 0){
	    printf(" idxPt<0     %d %d %f \n",nscale2,iNode2,cefscale2[iNode2]);
	    exit(1);
	 } else if (idxPt > NscalenodeScale2-1){
	    printf(" idxPt>max   %d %d %f \n",nscale2,iNode2,cefscale2[iNode2]);
	    exit(1);
	 }
	 // ---- check done ---- //
			    
	 // todo: check if the 4x4 grid of cefscale1 und cefscale2 is in sum 1.
	  

	 // ---- loop over all 4 _x_-nodes that receive contributions ---- //
	 for( int iX = 0; iX < 4; iX++) {           
	    int idxX = nx + (iX-1);     // --- the target x index
	    if (idxX < Nxtot1[ObsBin]) {
	       
	       // ---- if fini amplitude? then use 'weights' ---- //
	       if(itype == nlo::amplitude_dis::fini) { 
		  if (amp._M_fini.mode==0) {//finix
		     for(int proc=0;proc<NSubproc;proc++){
			SigmaTildeMuIndep[ObsBin][idxX][idxQ][idxPt][proc] += cefm[iX] * cefscale1[iNode1] * cefscale2[iNode2] * weights[0][proc]; //amp._M_fini.amp[0][proc]*coef;//  wt[proc];
			SigmaTildeMuFDep[ObsBin][idxX][idxQ][idxPt][proc] += cefm[iX] * cefscale1[iNode1] * cefscale2[iNode2] * weights[2][proc];//  wt[proc];
		     }
		  }		
		  if(amp._M_fini.mode==1){ //fini1
		     for(int proc=0;proc<NSubproc;proc++){
			SigmaTildeMuIndep[ObsBin][idxX][idxQ][idxPt][proc] += cefm[iX] * cefscale1[iNode1] * cefscale2[iNode2] * weights[1][proc]; //amp._M_fini.amp[1][proc]*coef;//  wt[proc];
			SigmaTildeMuFDep[ObsBin][idxX][idxQ][idxPt][proc] += cefm[iX] * cefscale1[iNode1] * cefscale2[iNode2] * weights[3][proc]; //amp._M_fini.amp[3][proc]*coef;//  wt[proc];
			SigmaTildeMuRDep[ObsBin][idxX][idxQ][idxPt][proc] += cefm[iX] * cefscale1[iNode1] * cefscale2[iNode2] * weights[4][proc]; //amp._M_fini.amp[4][proc]*coef;//  wt[proc];
		     }
		  }
	       }
	       // ---- no fini contribution -> use 'wt' ---- //
	       else { 
		  for(int proc=0;proc<NSubproc;proc++){
		     SigmaTildeMuIndep[ObsBin][idxX][idxQ][idxPt][proc] += cefm[iX] * cefscale1[iNode1] * cefscale2[iNode2] * wt[proc];
		  }
	       } // if fini
	    } // if idxX
	 } // of 4-xnodes
      } // for 4 scale2 nodes
   } // for 4 scale2 nodes  
}


//________________________________________________________________________________________________________________ //


void fnloBlockBNlojet::FillEventHHCMuVar(int ObsBin, double x1, double x2, double M1, double M2, const nlo::amplitude_hhc& amp, nlo::pdf_and_coupling_hhc& dummypdf, nlo::pdf_and_coupling_hhc& realpdf, double prefactor){
   // ------------------------------------------------------------------------------------------------------- //
   // 
   // Method for use in nlojet++-module for filling FastNLO sigmatilde-tables for v2.1 tables in pp and ppbar
   //
   //    ObsBin:		bin number of this event in the measurement binning
   //    x1, x2			x1 and x2 value of this event
   //    M1, M2			values M1 and M2 of your scales (scale1 and scale2). Usually 'Q' and 'p_T'
   //    amp			nlojet++ amplitude object
   //    dummypdf		dummypdf as nlojet++-pdf-object
   //    realpdf		like dummypdf but will evaluate real pdf table for reference tables
   //    prefactor		sth like:  w_PS * 1/x * alpha_em^2 h-bar^2 c^2
   //
   // ------------------------------------------------------------------------------------------------------- //

   // --- select interpolation kernel for x and for mu 
   //             1:CatmulRom   2:Lagrangian
   const EInterpolKernel eIkernx  = kCatmulRom;    // kLagrangian ?!?!
   const EInterpolKernel eIkernmu = kCatmulRom;    // kLagrangian ?!?!
   
   if(this->IRef>0){
     printf("fnloBlockBNlojet::FillEventHHCMuVar should not be called for refernce tables (Ref>0)\n");
     return;
   }

   // ---------------------------------------------------------------------------------------- //
   // ------------------------ Warm-Up Run to identify extreme x,mu values ------------------- //
   // ---------------------------------------------------------------------------------------- //

   if (IWarmUp == 1) {
      WarmUp( ObsBin , min(x1,x2) , M1 , M2 , "xlim" , "scale1" , "scale2");
      return;
   }

   // ---------------------------------------------------------------------------------------- //
   // --------------------------------- reference tables ------------------------------------- //
   // ---------------------------------------------------------------------------------------- //
   FillMuVarReferenceTables( ObsBin, M1, M2, amp, realpdf, prefactor);
   
   
   // ---- Half matrix x variables ---- //
   double xmin=0.0, xmax=0.0;
   if (x1 > x2) {
     xmax = x1;
     xmin = x2; }
   else {
     xmax = x2;
     xmin = x1;
   }

   // ---------------------------- check validity of call ------------------------------------ //
   if(this->NPDFDim != 1){
     printf("fnloBlockBNlojet::FillEventHHC: Error, only NPDFDim=1 (half matrix) implemented so far.\n");
     exit(1);   }
   if (xmin<XNode1[ObsBin][0]){
     printf("fnloBlockBNlojet::FillEventHHC: find: xmin (%f) smaller than lowest x-node (%f) for bin #%d .\n",
	    xmin,XNode1[ObsBin][0],ObsBin);   }
   if( NscalenodeScale1<=3 || NscalenodeScale2<=3){
      printf("fnloBlockBNlojet::FillEventHHCMuVar(). Error. Sorry, but you need some more scale nodes!\n");exit(1);
   }
   // ---------------------------------------------------------------------------------------- //

      
   // ---------------------------------------------------------------------------------------- //
   // **********  determine x_ij position in grid  ************
   // ---------------------------------------------------------------------------------------- //
   // --- determine fractional contributions
   double hxmin  = (this->*Fct_H_XNode)(xmin);
   double hxmax  = (this->*Fct_H_XNode)(xmax);
   double hxone   = 0.0;

   // --- define the x-node numbers in the range: 0 <= nxnode < Nxtot1
   double hxlimit = Hxlim1[ObsBin];
   int nxmin = int(Nxtot1[ObsBin] *(hxmin-hxlimit)/(hxone-hxlimit));
   int nxmax = int(Nxtot1[ObsBin] *(hxmax-hxlimit)/(hxone-hxlimit));
   if (nxmin < 0) nxmin = 0;  // move into available range
   if (nxmax < 0) nxmax = 0;

   // --- relative distances in h(xmin), h(xmax): deltamin,deltamax 
   double delta  = (hxone-hxlimit)/Nxtot1[ObsBin];
   double hxj = hxlimit+double(nxmin)/double(Nxtot1[ObsBin])*(hxone-hxlimit);
   double hxi = hxlimit+double(nxmax)/double(Nxtot1[ObsBin])*(hxone-hxlimit);
   double deltamin = (hxmin-hxj)/delta;
   double deltamax = (hxmax-hxi)/delta;

   // --- x-interpolation kernels and updated x-node positions: 1 <= nx < ntot-1
   int nxminf = 0;  // - variables for final nodes (after modification in kernel)
   int nxmaxf = 0;
   vector<double> cefmin = Interpol(nxmin, (Nxtot1[ObsBin]-1), deltamin, eIkernx, nxminf);
   vector<double> cefmax = Interpol(nxmax, (Nxtot1[ObsBin]-1), deltamax, eIkernx, nxmaxf);

   // --- PDF reweighting - compute weights, modify cefmax[.], cefmin[.] 
   //     but only those within grid, there are no nodes at x=1
   double pdfwgtmax = PDFwgt(xmax);
   for( int i1 = 0; i1 < 4; i1++) {
     if ((nxmaxf-1+i1) >= 0 && (nxmaxf-1+i1) < Nxtot1[ObsBin] ) {
       cefmax[i1] *= pdfwgtmax/PDFwgt(XNode1[ObsBin][nxmaxf-1+i1]);
     }
   }
   double pdfwgtmin = PDFwgt(xmin);
   for( int i2 = 0; i2 < 4; i2++) {
     if ((nxminf-1+i2) >= 0 && (nxminf-1+i2) < Nxtot1[ObsBin] ) {
       cefmin[i2] *= pdfwgtmin/PDFwgt(XNode1[ObsBin][nxminf-1+i2]);
     } 
   }

   // === the weights for the bi-cubic eigenfunctions (2-dim)
   double bicef[4][4];      
   for( int i1 = 0; i1 < 4; i1++) {
     for( int i2 = 0; i2 < 4; i2++) {
       bicef[i1][i2] = cefmax[i1] * cefmin[i2];
     }
   }




   // ---- calculate the matrix element ---- //
   // --- get coefficients - convert to our linear combinations ---- //
   double mu2 = M1*M1; // just a dummy scale

   // ---------------------------------------------------
   // remember from nlojet++ proc-hhc/weight.cc
   //   const char *weight_label_hhc[7] = {"gg", "qg", "gq", "qr", "qq", "qqb", "qrb"};
   // ---------------------------------------------------
   // remember from nlojet++ proc-hhc/process.cc
   //    retval[0] = A0*B0;
   //    retval[1] = (A + Ab)*B0;
   //    retval[2] = A0*(B + Bb);
   //    retval[3] = A*B + Ab*Bb - D;
   //    retval[4] = D;
   //    retval[5] = Db;
   //    retval[6] = A*Bb +Ab*B - Db;
   // ---------------------------------------------------
   //
   //

   //nlo::weight_hhc wt = amp(pdf,mu2,mu2,prefactor);
   nlo::weight_hhc wtorg = amp(dummypdf,mu2,mu2,prefactor);
   // - rearrange subprocesses
   nlo::weight_hhc wt = wtorg;
   wt[0] = wtorg[0];
   wt[1] = wtorg[3];
   wt[2] = wtorg[4];
   wt[3] = wtorg[5];
   wt[4] = wtorg[6];
   wt[5] = wtorg[1];
   wt[6] = wtorg[2];
   // -- case NSubproc=6: see below
   
   wt *= 389385730.;
   if(IXsectUnits!=12){
     wt *= pow(10.,(IXsectUnits-12)) ;
   }
   
   // deal with subprocesses 2 and 3 -> moved to 6 and 7
   //    - if x1>x2 -> o.k.
   //    - if x2>x1 -> swap weights for subprocesses 2,3 -> now 6,7
   if(x2>x1){
      double buffer;
      //       buffer = wt[1];
      //       wt[1] = wt[2];
      //       wt[2] = buffer;
      buffer  = wt[5]; // swap subprocesses 6,7
      wt[5] = wt[6];
      wt[6] = buffer;
   }
   //    // --- combine subprocesses 5,6 here after possible swapping
   //    if (NSubproc == 6) {   
   //      wt[5] = (wt[5]+wt[6])/2.;
   //      wt[6] = wt[5];    // set equal in case of swapping below diagonal (at bottom) 
   //    }
   
   

   // ---- calcualte weights for fini contributions ---- //
   nlo::amplitude_hhc::contrib_type itype = amp.contrib();
   double coef = 389385730.*prefactor;
   if(IXsectUnits!=12)  coef *= pow(10.,(IXsectUnits-12)) ;
   nlo::weight_hhc cPDF = dummypdf.pdf(x1,x2,M1*M1,2,3); // 1/x1/x2
   
   

   double weights[7][7]; // weights[amp_i][proc]
   //nlo::weight_hhc weights[7]; nicer, but I do net get the syntax correct
   if(itype == nlo::amplitude_hhc::fini) { 
     for ( int kk = 0 ; kk<7 ; kk ++ ){ // contrib amp_i
	weights[kk][0] = amp._M_fini.amp[kk][0]*coef*cPDF[0];//wtorg[0];
	weights[kk][1] = amp._M_fini.amp[kk][3]*coef*cPDF[3];
	weights[kk][2] = amp._M_fini.amp[kk][4]*coef*cPDF[4];
	weights[kk][3] = amp._M_fini.amp[kk][5]*coef*cPDF[5];
	weights[kk][4] = amp._M_fini.amp[kk][6]*coef*cPDF[6];
	weights[kk][5] = amp._M_fini.amp[kk][1]*coef*cPDF[1];
	weights[kk][6] = amp._M_fini.amp[kk][2]*coef*cPDF[2];
	//        for ( int nproc = 0 ; nproc<7 ; nproc ++ ){ // nsubproc
	// 	 weights[kk][nproc] = amp._M_fini.amp[kk][nproc]*coef*cPDF[nproc];//wtorg[0];
	//        }
     }     
   }

   // neu
   if(x2>x1){
     for ( int kk = 0 ; kk<7 ; kk ++ ){ // contrib amp_i
	double buffer;
	buffer = weights[kk][5];
	weights[kk][5] = weights[kk][6];
	weights[kk][6] = buffer;
     }
   }


   /*
   // ---------------------- cross check ------------------------ //
   //
   //    we calculate here some control cross section with a weired scale (2*pt)^2
   //        - 'wtctr'
   //    further we perform the xs-calculation with single amplitude contributions
   //        - 'wself'
   //    we look, if wtctr is the same as wself
   //
   // ----------------------------------------------------------- //

   // ----- calculate weights for cross checks ---- //
   mu2 = 2*M1*2*M1;
   nlo::weight_hhc wtctrlorg = amp(dummypdf, mu2, mu2 ,prefactor); // mu is a 'dummy'-sca

   //    // - rearrange subprocesses
   nlo::weight_hhc wtctrl = wtctrlorg;
   //    wtctrl[0] = wtctrlorg[0];
   //    wtctrl[1] = wtctrlorg[3];
   //    wtctrl[2] = wtctrlorg[4];
   //    wtctrl[3] = wtctrlorg[5];
   //    wtctrl[4] = wtctrlorg[6];
   //    wtctrl[5] = wtctrlorg[1];
   //    wtctrl[6] = wtctrlorg[2];
   //    // no further shuffling...

   wtctrl *= 389385730.;
   if(IXsectUnits!=12)  wtctrl *= pow(10.,(IXsectUnits-12)) ;

   amp(dummypdf, M1*M1, M1*M1 ,prefactor); // call it with another scale for proper cross check!

   // ----- calculate weights by myself for cross checks ---- //
   double wself[7] = {0};
 
   if(itype == nlo::amplitude_hhc::fini) { 
	
     if (amp._M_fini.mode==0) { //finix1
       for(int proc=0;proc<NSubproc;proc++){
 	    wself[proc] += weights[0][proc];
 	    wself[proc] += weights[3][proc] * std::log(mu2);
 	    printf("fini.mode=0: proc=[%d], MuIndep=%9.5f, MuFDep=%9.5f, ln(mu2)=%9.5f, mu2=%9.5f, w*ln=%9.5f \n",
 		   proc,weights[0][proc] , weights[3][proc] , std::log(mu2) , mu2, weights[3][proc] * std::log(mu2)); 
 	  }
 	}		
     if (amp._M_fini.mode==1) { //finix2
       for(int proc=0;proc<NSubproc;proc++){
 	    wself[proc] += weights[1][proc];
 	    wself[proc] += weights[4][proc] * std::log(mu2);
 	    printf("fini.mode=1: proc=[%d], MuIndep=%9.5f, MuFDep=%9.5f, ln(mu2)=%9.5f, mu2=%9.5f, w*ln=%9.5f \n",
 		   proc,weights[1][proc] , weights[4][proc] , std::log(mu2) , mu2, weights[4][proc] * std::log(mu2)); 
 	  }
 	}		
     if(amp._M_fini.mode==2){ //fini1
 	  for(int proc=0;proc<NSubproc;proc++){
 	    wself[proc] += weights[2][proc];
 	    wself[proc] += weights[5][proc] * std::log(mu2);
 	    wself[proc] += weights[6][proc] * std::log(mu2);
 	    printf("fini.mode=2: proc=[%d], MuIndep=%9.5f, MuFDep=%9.5f, MuRDep=%9.5f, ln(mu2)=%9.5f, mu2=%9.5f, w_f*ln=%9.5f, w_r*ln=%9.5f \n",
 		   proc,weights[2][proc],weights[5][proc],weights[6][proc], std::log(mu2) , mu2, weights[5][proc] * std::log(mu2), weights[6][proc] * std::log(mu2)); 
 	  }
     }
   }
   else { // no fini contribution
     printf("no fini contribution. wt==wself\n"); 
     for(int proc=0;proc<NSubproc;proc++){
       wself[proc] += wt[proc]; // mu should not matter here!
     }
   }
    
   for(int proc=0;proc<NSubproc;proc++){
     printf("proc=[%d]: wtFNLO=%9.6f, wself=%9.6f\n",proc,wtctrl[proc],wself[proc]);
   }

   return;
   // ------------------------ cross check ende ----------------------------------- //
   */


   // ---------------------------------------------------------------------------------------- //
   // ------------------------------ calculation of the scale nodes -------------------------- //
   // ---------------------------------------------------------------------------------------- //
   int scalenode1 = NscalenodeScale1-2;  // --- initialize with largest possible value
   
   // --- find scale position in range:  0 <= scalenode1 < nscalenode1-1
   for( int iNode=1 ; iNode<NscalenodeScale1 ; iNode++ ){
      if ( M1 < ScaleNode1[ObsBin][iNode]){
	 scalenode1=iNode-1;
	 break;
      }
   }
   // --- relative distance delta - in function H(mu)
   // deltascale (Interpol(.,.,.delta,.): relative distance of value to node 'nnode'
   double deltascale1 = 
      ( (this->*Fct_H_Scale[0])(M1) - HScaleNode1[ObsBin][scalenode1] ) /
      (HScaleNode1[ObsBin][scalenode1+1] - HScaleNode1[ObsBin][scalenode1]);
   
   // --- get scale interpolation kernel and updated scalenode position: 1 <= nmu < ntot-2
   int nscale1 = 0;                      // --- variable for final scale1-node
   vector<double> cefscale1 =  Interpol(scalenode1, NscalenodeScale1-2, deltascale1, eIkernmu, nscale1);

      
      
   // ---------------------------------------------------------------------------------------- //
   // ---------------------------- calculate second scale  interpolation --------------------- //
   // ---------------------------------------------------------------------------------------- //
   int scalenode2 = NscalenodeScale2-2;  // --- initialize with largest possible value

   // --- find scale position in range:  0 <= scalenode2 < nscalenode2-1
   for( int iNode=1 ; iNode<NscalenodeScale2 ; iNode++ ){
      //printf("test node %d %d %f %f \n",ObsBin,i,scale2,ScaleNode[ObsBin][i]);
      if ( M2 < ScaleNode2[ObsBin][iNode]){
	 scalenode2=iNode-1;
	 break;
      }
   }
   // --- relative distance delta - in function H(mu)
   // deltascale (Interpol(.,.,.delta,.): relative distance of value to node 'nnode'
   double deltascale2 = 
      ( (this->*Fct_H_Scale[1])(M2) - HScaleNode2[ObsBin][scalenode2] ) /
      (HScaleNode2[ObsBin][scalenode2+1] - HScaleNode2[ObsBin][scalenode2]);
   
   // --- get scale interpolation kernel and updated scalenode position: 1 <= nmu < ntot-2
   int nscale2 = 0;                      // --- variable for final scale2-node
   vector<double> cefscale2 =  Interpol(scalenode2, NscalenodeScale2-2, deltascale2, eIkernmu, nscale2);
   


   // ---------------------------------------------------------------------------------------- //
   // -------------- loop over all 4x4x4x4  nodes that receive contributions ----------------- //
   // ---------------------------------------------------------------------------------------- //
   for(int iNode1 = 0; iNode1 < 4; iNode1++ ){
      int idx1 = nscale1 + (iNode1-1);     // --- the target scale1 index
      for(int iNode2 = 0; iNode2 < 4; iNode2++ ){
	 int idx2 = nscale2 + (iNode2-1);     // --- the target scale2 index
       
	 // --- check: all elements should end up within grid
	 // todo
              
	 // --- loop over all 16 xmin,xmax points that receive contributions
	 for( int i1 = 0; i1 < 4; i1++) {           
	    for (int i2 = 0; i2 < 4; i2++){
	       int xmaxbin = nxmaxf + (i1 -1);   // - the target x-max index
	       int xminbin = nxminf + (i2 -1);   // - the target x-min index
	       if (xmaxbin < Nxtot1[ObsBin] && xminbin < Nxtot1[ObsBin]) {

		  nlo::weight_hhc wtmp = wt;  // - a working copy of the weights
		  // --- check if above diagonal? project back and swap qg<->gq
		  if (xminbin>xmaxbin) { 
		     int di = xminbin - xmaxbin;
		     xmaxbin = xmaxbin + di;   // modify indicees
		     xminbin = xminbin - di;		

		     double buffer;
// 		     buffer = wtmp[1]; // swap subprocesses 2,3 // if you use nlojet-type PDF-LCs
// 		     wtmp[1] = wtmp[2];
// 		     wtmp[2] = buffer;
		     buffer  = wtmp[5]; // swap subprocesses 6,7
		     wtmp[5] = wtmp[6];
		     wtmp[6] = buffer;
		     for ( int kk = 0 ; kk<7 ; kk ++ ){ // contrib amp_i
			double weightsbuffer; 
// 			weightsbuffer	= weights[kk][1];
// 			weights[kk][1]	= weights[kk][2];
// 			weights[kk][2]	= weightsbuffer;
			weightsbuffer	= weights[kk][5];
			weights[kk][5]	= weights[kk][6];
			weights[kk][6]	= weightsbuffer;
		     }
		     // 	       // new ordering (1,2)->(5,6)	 
		     // 	       buffer  = wtmp[5]; // swap subprocesses 6,7
		     // 	       wtmp[5] = wtmp[6];
		     // 	       wtmp[6] = buffer;
		     // 	       // shuffle also the 'weights'
		     // 	       for ( int kk = 0 ; kk<7 ; kk ++ ){ // contrib amp_i
		     // 		 double weightsbuffer; 
		     // 		 //weightsbuffer	= weights[kk][1];
		     // 		 //weights[kk][1]	= weights[kk][2];
		     // 		 //weights[kk][2]	= weightsbuffer;
		     // 		 weightsbuffer	= weights[kk][5];
		     // 		 weights[kk][5]	= weights[kk][6];
		     // 		 weights[kk][6]	= weightsbuffer;
		     //              }
		  }
	     
		  int im = GetXIndex(ObsBin,xminbin,xmaxbin);
		  // printf("fastNLO: index %d in xmaxbin #%d xminbin #%d\n",im,xmaxbin,xminbin);

		  double fac = bicef[i1][i2] * cefscale1[iNode1] * cefscale2[iNode2];
		  if(itype == nlo::amplitude_hhc::fini) { 
		     if (amp._M_fini.mode==0) { //finix1
			for(int proc=0;proc<NSubproc;proc++){
			   SigmaTildeMuIndep[ObsBin][im][idx1][idx2][proc] +=  fac * weights[0][proc];
			   SigmaTildeMuFDep [ObsBin][im][idx1][idx2][proc] +=  fac * weights[3][proc];
			}
		     }		
		     else if (amp._M_fini.mode==1) { //finix2
			for(int proc=0;proc<NSubproc;proc++){
			   SigmaTildeMuIndep[ObsBin][im][idx1][idx2][proc] +=  fac * weights[1][proc];
			   SigmaTildeMuFDep [ObsBin][im][idx1][idx2][proc] +=  fac * weights[4][proc];
			}
		     }		
		     else if(amp._M_fini.mode==2){ //fini1
			for(int proc=0;proc<NSubproc;proc++){
			   SigmaTildeMuIndep[ObsBin][im][idx1][idx2][proc] +=  fac * weights[2][proc];
			   SigmaTildeMuFDep [ObsBin][im][idx1][idx2][proc] +=  fac * weights[5][proc];
			   SigmaTildeMuRDep [ObsBin][im][idx1][idx2][proc] +=  fac * weights[6][proc];
			}
		     }
		  }
		  else { // no fini contribution
		     for(int proc=0;proc<NSubproc;proc++){
			SigmaTildeMuIndep[ObsBin][im][idx1][idx2][proc] +=  fac * wtmp[proc];
		     }
		  }
	       }
	    }
	 }
      }
   }

}


// ___________________________________________________________________________________________________ //


void fnloBlockBNlojet::WarmUp( int ObsBin, double x, double M1, double M2, string sx, string s1, string s2 ){

   // -------------------------------------------------------------------------- //
   //  
   //  WarmUp(). Simple method needed for 'warm-up' runs. 
   //  Can be called within filling routines for IWarmUp==1 values
   //  
   //  Looks and writes out for:
   //     - smallest x-value in each bin
   //     - smallest and largest scale value in each bin (for up to two scales)
   //
   //  Write-Out is done
   //     - on screen
   //     - into file called: "fastNLO-warmup.dat"
   //
   //  Input:
   //     - x		x-value: repeat for many events
   //     - M1, M2.	scale values (M2 is optional)
   //     - sx, s1, s2	names that are used in the write-out of this run for each variable (x, M1, M2)
   //
   // -------------------------------------------------------------------------- //

   //static unsigned long counter = 0;
   //    static double* axlo = NULL;
   //    static double* a1lo = NULL;
   //    static double* a1up = NULL;
   //    static double* a2lo = NULL;
   //    static double* a2up = NULL;

   // init arrays
   if ( counter == 0 ){
      axlo = new double[BlockA2->GetNObsBin()];
      a1lo = new double[BlockA2->GetNObsBin()];
      a1up = new double[BlockA2->GetNObsBin()];
      a2lo = new double[BlockA2->GetNObsBin()];
      a2up = new double[BlockA2->GetNObsBin()];
      for ( int i = 0 ; i<BlockA2->GetNObsBin() ; i++ ){
	 axlo[i] = 0; 
	 a1lo[i] = 0;
	 a1up[i] = 0;
	 a2lo[i] = 0;
	 a2up[i] = 0;
     }
   }

   // init starting values
   if ( axlo[ObsBin] == 0.) axlo[ObsBin] = x;
   if ( a1lo[ObsBin] == 0.) a1lo[ObsBin] = M1;
   if ( a2lo[ObsBin] == 0.) a2lo[ObsBin] = M2;
   
   // look for smallest/largest values
   if ( axlo[ObsBin] > x)  axlo[ObsBin] = x;
   if ( a1lo[ObsBin] > M1) a1lo[ObsBin] = M1;
   if ( a1up[ObsBin] < M1) a1up[ObsBin] = M1;
   if ( a2lo[ObsBin] > M2) a2lo[ObsBin] = M2;
   if ( a2up[ObsBin] < M2) a2up[ObsBin] = M2;
      
   // ---- count ---- //
   counter++;
      
   // ---- write out --- //
   if ( (counter % IWarmUpPrint) == 0) {
      // ---- printout on screen ---- //
      printf(" // %16lu contributions (!= events) in warm-up run \n",counter);
      for (int i=0;i<BlockA2->GetNObsBin();i++){
	 printf("	%s [%d] = %8.2e", sx.data(), i, axlo[i] ); // xmin
	 if ( a1lo[0] != 0 ) {	    printf(" , %slo [%d] = %9.4f , %shi [%d] = %9.4f", 	   s1.data(),  i, a1lo[i], s1.data(),  i, a1up[i] );	 }
	 if ( a2lo[0] != 0 ) {	    printf(" , %slo [%d] = %9.4f , %shi [%d] = %9.4f", 	   s2.data(),  i, a2lo[i], s2.data(),  i, a2up[i] );	 }
	 printf(";\n");
      }
      
      // ---- write out to file ---- //
      FILE * ofile;
      ofile = fopen("fastNLO-warmup.dat","w");
      fprintf(ofile,"      // %lu contributions (!= events) in warm-up run \n",counter);
      for (unsigned int i=0;i<BlockA2->GetNObsBin();i++){
	 fprintf(ofile,"	%s [%d] = %8.2e", sx.data(), i, axlo[i] );
	 if ( a1lo[0] != 0 ) {	    fprintf(ofile," , %slo [%d] = %9.4f , %shi [%d] = %9.4f", 	   s1.data(),  i, a1lo[i], s1.data(),  i, a1up[i] );	 }
	 if ( a2lo[0] != 0 ) {	    fprintf(ofile," , %slo [%d] = %9.4f , %shi [%d] = %9.4f", 	   s2.data(),  i, a2lo[i], s2.data(),  i, a2up[i] );	 }
	 fprintf(ofile,";\n");
      }
      fclose(ofile);
   }

   return;

}


// ___________________________________________________________________________________________________ //


void fnloBlockBNlojet::FillMuVarReferenceTables(int ObsBin, double M1, double M2, const nlo::amplitude_dis& amp, nlo::pdf_and_coupling_dis& realpdf, double prefactor){
   // -------------------------------------------------------------------------- //
   //  
   //  FillMuVarReferenceTables(). Method for filling the three DIS
   //  reference tables 
   //  Fills:
   //     - SigmaRefMixed
   //     - SigmaRef_s1
   //     - SigmaRef_s2
   //
   // -------------------------------------------------------------------------- //

   // simplify calls:
   vector<vector<double> >* SigmaRef[3] = { &SigmaRefMixed, &SigmaRef_s1, &SigmaRef_s2 };

   // ---- SigmaRefMixed ---- //
   for ( int iR = 0 ; iR < 3 ; iR++ ) {
      if ( Fct_MuR_Ref[iR] && Fct_MuF_Ref[iR] ){
	 double mur2 = pow((Fct_MuR_Ref[iR])(M1,M2),2);
	 double muf2 = pow((Fct_MuF_Ref[iR])(M1,M2),2);
	 if ( mur2 < 1. ){
	    printf("fnloBlockBNlojet::FillMuVarReferenceTables. Sorry, but your composite scale for reference cross section %d is only %7.4f GeV small. This seems to be unphysical and leads to 'nan'.\n",iR,sqrt(mur2));
	    exit(1);
	 }
		nlo::weight_dis wt = amp(realpdf,mur2,muf2,prefactor); // the REAL pdf
	   double tmp = wt[1]+wt[2];
	   wt[1] = wt[0];
	   wt[0] = tmp;
	   wt[2] = 0;
	   wt *= 389385730.;
	   if(IXsectUnits!=12)  wt *= pow(10.,(IXsectUnits-12)) ;
	 for(int proc=0;proc<NSubproc;proc++){
	    (*SigmaRef[iR])[ObsBin][proc] += wt[proc];
	 }
      }      
   }
}


// ___________________________________________________________________________________________________ //


void fnloBlockBNlojet::FillMuVarReferenceTables(int ObsBin, double M1, double M2, const nlo::amplitude_hhc& amp, nlo::pdf_and_coupling_hhc& realpdf, double prefactor){
   // -------------------------------------------------------------------------- //
   //  
   //  FillMuVarReferenceTables(). Method for filling the three
   //  reference tables that come with the MuVar tables
   //
   //  Fills:
   //     - SigmaRefMixed
   //     - SigmaRef_s1
   //     - SigmaRef_s2
   //
   // -------------------------------------------------------------------------- //

   // simplify calls:
   vector<vector<double> >* SigmaRef[3] = { &SigmaRefMixed, &SigmaRef_s1, &SigmaRef_s2 };

   // ---- SigmaRefMixed ---- //
   for ( int iR = 0 ; iR < 3 ; iR++ ) {
      if ( Fct_MuR_Ref[iR] && Fct_MuF_Ref[iR] ){
	 double mur2 = pow((Fct_MuR_Ref[iR])(M1,M2),2);
	 double muf2 = pow((Fct_MuF_Ref[iR])(M1,M2),2);
	 if ( mur2 < 1. ){
	    printf("fnloBlockBNlojet::FillMuVarReferenceTables. Sorry, but your composite scale for reference cross section %d is only %7.4f GeV small. This seems to be unphysical and leads to 'nan'.\n",iR,sqrt(mur2));
	    exit(1);
	 }
	 nlo::weight_hhc wt = amp(realpdf,mur2,muf2,prefactor); // the REAL pdf
	 wt *= 389385730.;
	 if(IXsectUnits!=12)  wt *= pow(10.,(IXsectUnits-12)) ;
	 for(int proc=0;proc<NSubproc;proc++){
	    (*SigmaRef[iR])[ObsBin][proc] += wt[proc];
	 }
      }      
   }
}


// ___________________________________________________________________________________________________ //


void fnloBlockBNlojet::SetFuncMuForReference( double (*func_mur)(double,double) , double (*func_muf)(double,double) , int iRefTable ){
   printf(" *  fnloBlockBNlojet::SetFuncMuForReference(). Test call for renormalization scale.\n");
   printf(" *    Scale1 = 1 ,      Scale2 = 1        ->  mu = func(1,1)             = %9.4f\n",(*func_mur)(1,1));
   printf(" *    Scale1 = 91.1876, Scale2 = 91.1876  ->  mu = func(91.1876,91.1876) = %9.4f\n",(*func_mur)(91.1876,91.1876));
   printf(" *    Scale1 = 1,       Scale2 = 91.1876  ->  mu = func(1,91.1876)       = %9.4f\n",(*func_mur)(1,91.1876));
   printf(" *    Scale1 = 91.1876, Scale2 = 1        ->  mu = func(91.1876,1)       = %9.4f\n",(*func_mur)(91.1876,1));
   printf(" *  fnloBlockBNlojet::SetFuncMuForReference(). Test call for factorization scale.\n");
   printf(" *    Scale1 = 1 ,      Scale2 = 1        ->  mu = func(1,1)             = %9.4f\n",(*func_muf)(1,1));
   printf(" *    Scale1 = 91.1876, Scale2 = 91.1876  ->  mu = func(91.1876,91.1876) = %9.4f\n",(*func_muf)(91.1876,91.1876));
   printf(" *    Scale1 = 1,       Scale2 = 91.1876  ->  mu = func(1,91.1876)       = %9.4f\n",(*func_muf)(1,91.1876));
   printf(" *    Scale1 = 91.1876, Scale2 = 1        ->  mu = func(91.1876,1)       = %9.4f\n",(*func_muf)(91.1876,1));
   
   if ( iRefTable >= 3 || iRefTable < 0 ){
      printf(" *  fnloBlockBNlojet::SetFuncMuForReference(). Error. iRefTable must be 0, 1 or 2, but is %d.\n",iRefTable);
      exit(1);
   }

   Fct_MuR_Ref[iRefTable]	= func_mur;
   Fct_MuF_Ref[iRefTable]	= func_muf;
}


// ___________________________________________________________________________________________________ //


void fnloBlockBNlojet::FillEventPhoto(int ObsBin, double x, double scale1, const nlo::amplitude_photo& amp, nlo::pdf_and_coupling_photo& pdf, double prefactor){

   fnloBlockA2 *A2 =  BlockA2;
   
   if(this->IRef>0){

      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){
         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         nlo::weight_photo wt = amp(pdf,mu2,mu2,prefactor);

         double binsize = 1.0;
         for(int dim=0; dim<A2->NDim; dim++){
            binsize *= (A2->UpBin[ObsBin][dim] - A2->LoBin[ObsBin][dim]);
         }
         
         wt *= 389385730./137.0/binsize;
         if(IXsectUnits!=12){
            wt *= pow(10.,(IXsectUnits-12)) ;
         }
         
         for(int proc=0;proc<NSubproc;proc++){
            SigmaTilde[ObsBin][scalevar][0][0][proc] += wt[proc];
         }
      }
   }else{

      if (x<XNode1[ObsBin][0]){
         printf("fnloBlockBNlojet::FillEventPhoto: Error: x (%f) smaller than lowest point (%f) at bin #%d .\n",
                x,XNode1[ObsBin][0],ObsBin);
         exit(1);
      }

      // **********  determine x_ij position in grid  ************
      //--- determine fractional contributions
      double hxlimit = Hxlim1[ObsBin];
      double hx = log10(x);
      double hxone = 0.0;

      // define the x-bin number in the range  [0:ntot[
      int nx = int(Nxtot1[ObsBin] *(hx-hxlimit)/(hxone-hxlimit));

      //-- relative distances in h(x): deltam
      double delta  = (hxone-hxlimit)/Nxtot1[ObsBin];
      double hxi = hxlimit+double(nx)/double(Nxtot1[ObsBin])*(hxone-hxlimit);
      double deltam = (hx-hxi)/delta;

      // ===== variables for the bi-cubic interpolation =====
      // === the relative distances to the four nearest bins
      vector<double> cm(4) ; 
      cm[0] = deltam+1.0;
      cm[1] = deltam;
      cm[2] = 1.0-deltam;
      cm[3] = 2.0-deltam;

      vector<double> cefm(4) ; 
      if (nx==0 || nx==(Nxtot1[ObsBin]-1)) { //linear in 1st and last bin
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

      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){

         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         nlo::weight_photo wt = amp(pdf,mu2,mu2,prefactor);

         double pdfsigma = 1.0/3.0*(-wt[1] + 4.*wt[2]);
         double pdfdelta = 3.*(wt[1]-wt[2]);
         //wt[1] = pdfdelta;
         //wt[2] = pdfsigma;
	 // MW: order should be 0:Delta, 1:Gluon, 2:Sigma
	 wt[1] = wt[0];
	 wt[0] = pdfdelta;
	 wt[2] = pdfsigma;

         double binsize = 1.0;
         for(int dim=0; dim<A2->NDim; dim++){
            binsize *= (A2->UpBin[ObsBin][dim] - A2->LoBin[ObsBin][dim]);
         }

         wt *= 389385730./137.0/binsize;
         if(IXsectUnits!=12){
            wt *= pow(10.,(IXsectUnits-12)) ;
         }

         // define the scale-bin number in the range  [0:nscalebin[
         int scalenode = 0;
         double cefscale[4] = {0.0,0.5,0.5,0.0};
         if(Nscalenode[0]>1){
            for(int i=1;i<Nscalenode[0];i++){
               if (ScaleFac[0][scalevar]*scale1<ScaleNode[ObsBin][0][scalevar][i]){
                  scalenode=i;
                  break;
               }
            }
            double deltascale = (ScaleFac[0][scalevar]*scale1 - ScaleNode[ObsBin][0][scalevar][scalenode-1])/
               (ScaleNode[ObsBin][0][scalevar][scalenode]-ScaleNode[ObsBin][0][scalevar][scalenode-1]);
            vector<double> cscale(4) ; 
            cscale[0] = deltascale+1.0;
            cscale[1] = deltascale;
            cscale[2] = 1.0-deltascale;
            cscale[3] = 2.0-deltascale;

            if(scalenode<2 || scalenode>ScaleNode[ObsBin][0][scalevar][scalenode-1]-3){
               cefscale[0] = 0.0;
               cefscale[1] = 1.0-cscale[1];
               cefscale[2] = 1.0-cscale[2];
               cefscale[3] = 0.0;
            }else{
               cefscale[1]=1.0-2.5*cscale[1]*cscale[1]+1.5*cscale[1]*cscale[1]*cscale[1];
               cefscale[2]=1.0-2.5*cscale[2]*cscale[2]+1.5*cscale[2]*cscale[2]*cscale[2];
               cefscale[0]=2.0 - 4.0*cscale[0] + 2.5*cscale[0]*cscale[0]
                  - 0.5*cscale[0]*cscale[0]*cscale[0];
               cefscale[3]=2.0 - 4.0*cscale[3] + 2.5*cscale[3]*cscale[3]
                  - 0.5*cscale[3]*cscale[3]*cscale[3];
            }
         }

         // ** loop over all 4 scale nodes that receive contributions
         for(int i3 = 0; i3 < 4; i3++){
            int is = scalenode +i3 -2;  // the target scale index
            if(is<0){
               is = 0;
            }
            if (is> Nscalenode[0]-1){
               is = Nscalenode[0]-1;
            }
            // ** loop over all 4 x points that receive contributions
            for( int i1 = 0; i1 < 4; i1++) {           
               int x1bin = (nx +i1 -1);
               if(x1bin<0) x1bin = 0;
               if(x1bin>Nxtot1[ObsBin]-1) x1bin =Nxtot1[ObsBin]-1;
               int im = x1bin;   // the target x index
               for(int proc=0;proc<NSubproc;proc++){
                  //                     printf("%d %d %d %d %d %f\n",ObsBin,scalevar,is,im,proc,cefscale[i3]);
                  SigmaTilde[ObsBin][scalevar][is][im][proc] +=  cefm[i1]  * cefscale[i3] * wt[proc];
               }
            }
         }
      }
   }
}

void fnloBlockBNlojet::FillEventResolved(int ObsBin, double x, double x2, double scale1, double ymin, double ymax, double Q2max, const nlo::amplitude_hhc& amp, nlo::pdf_and_coupling_hhc& pdf, double prefactor){
   fnloBlockA2 *A2 =  BlockA2;
   
   if(this->IRef>0){

      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){
         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         nlo::weight_hhc wt = amp(pdf,mu2,mu2,prefactor);

         double binsize = 1.0;
         for(int dim=0; dim<A2->NDim; dim++){
            binsize *= (A2->UpBin[ObsBin][dim] - A2->LoBin[ObsBin][dim]);
         }
         
         wt *= 389385730./binsize;
         if(IXsectUnits!=12){
            wt *= pow(10.,(IXsectUnits-12)) ;
         }
         
         for(int proc=0;proc<NSubproc;proc++){
            SigmaTilde[ObsBin][scalevar][0][0][proc] += wt[proc];
         }
      }
   }else{

      if (x<XNode1[ObsBin][0]){
         printf("fnloBlockBNlojet::FillEventResolved: Error: x (%f) smaller than lowest point (%f) at bin #%d .\n",
                x,XNode1[ObsBin][0],ObsBin);
         exit(1);
      }
      //--- determine fractional contribution x1
      double hxlimit = Hxlim1[ObsBin];
      double hx = log10(x);
      double hxone = 0.0;

      // define the x-bin number in the range  [0:ntot[
      int nx = int(Nxtot1[ObsBin] *(hx-hxlimit)/(hxone-hxlimit));

      //-- relative distances in h(x): deltam
      double delta  = (hxone-hxlimit)/Nxtot1[ObsBin];
      double hxi = hxlimit+double(nx)/double(Nxtot1[ObsBin])*(hxone-hxlimit);
      double deltam = (hx-hxi)/delta;

      // ===== variables for the bi-cubic interpolation =====
      // === the relative distances to the four nearest bins
      vector<double> cm(4) ; 
      cm[0] = deltam+1.0;
      cm[1] = deltam;
      cm[2] = 1.0-deltam;
      cm[3] = 2.0-deltam;

      vector<double> cefm(4) ; 
      if (nx==0 || nx>=(Nxtot1[ObsBin]-3)) { //linear in 1st and last bin
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


      double thisymin = x2 > ymin ? x2 : ymin;
      if(x2 >= ymax) return;
      double y, jac = 0.5*std::log(ymax/thisymin);
      // loop over y-range
      for(unsigned int ib = 0; ib < 20; ib++) {
         y = thisymin*std::exp((_M_xb[ib]+1)*jac);
         double x2prime = x2/y;
         double photoweight = jac*_M_wb[ib]*y/x2;
         const double Me2 = 0.00000026112004954086;  //  GeV^2
         double Q2min = Me2*y*y/(1-y);
         photoweight *=  1.0/137.0*((1+(1-y)*(1-y))/y*log(Q2max/Q2min) 
                                    + 2.0*Me2*y*(1.0/Q2max-1.0/Q2min))/6.28318530717958647692; 

         if (x2prime<XNode2[ObsBin][0]){
            printf("fnloBlockBNlojet::FillEventResolved: Error: x2 (%f) smaller than lowest point (%f) at bin #%d .\n",
                   x2prime,XNode2[ObsBin][0],ObsBin);
            exit(1);
         }

         // define the x-bin number in the range  [0:ntot[
         int nx2;
         //-- relative distances in h(x): deltam
         double deltam;

         if(x2prime<0.9){
            //--- determine fractional contribution x2
            double hxlimit = Hxlim2[ObsBin];
            double hx = log10(x2prime);
            double hxone = log10(0.9);

            // define the x-bin number in the range  [0:ntot[
            nx2 = int( (Nxtot2[ObsBin]-20) *(hx-hxlimit)/(hxone-hxlimit));

            //-- relative distances in h(x): deltam
            double delta  = (hxone-hxlimit)/ (Nxtot2[ObsBin]-20);
            double hxi = hxlimit+double(nx2)/double(Nxtot2[ObsBin]-20)*(hxone-hxlimit);
            deltam = (hx-hxi)/delta;
         }else{
            nx2 = int( (Nxtot2[ObsBin]-20) + 20 *(x2prime-0.9)/(1.0-0.9));
            deltam = x2prime - XNode2[ObsBin][nx2];
         }

         //         printf("x2= %g nx2=%d deltam=%g\n",x2prime,nx2,deltam);
         // ===== variables for the bi-cubic interpolation =====
         // === the relative distances to the four nearest bins
         vector<double> cm2(4) ; 
         cm2[0] = deltam+1.0;
         cm2[1] = deltam;
         cm2[2] = 1.0-deltam;
         cm2[3] = 2.0-deltam;

         vector<double> cefm2(4) ; 
         if (nx2==0 || nx2>=(Nxtot2[ObsBin]-3)) { //linear in 1st and last bin
            cefm2[0] = 0.0;
            cefm2[1] = 1.0-cm2[1];
            cefm2[2] = 1.0-cm2[2];
            cefm2[3] = 0.0; }
         else {                              // cubic in the middle
            cefm2[1]=1.0-2.5*cm2[1]*cm2[1]+1.5*cm2[1]*cm2[1]*cm2[1];
            cefm2[2]=1.0-2.5*cm2[2]*cm2[2]+1.5*cm2[2]*cm2[2]*cm2[2];
            cefm2[0]=2.0 - 4.0*cm2[0] + 2.5*cm2[0]*cm2[0]
               - 0.5*cm2[0]*cm2[0]*cm2[0];
            cefm2[3]=2.0 - 4.0*cm2[3] + 2.5*cm2[3]*cm2[3]
               - 0.5*cm2[3]*cm2[3]*cm2[3];
         }

         for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){

            double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
            nlo::weight_hhc wt = amp(pdf,mu2,mu2,prefactor);

            double binsize = 1.0;
            for(int dim=0; dim<A2->NDim; dim++){
               binsize *= (A2->UpBin[ObsBin][dim] - A2->LoBin[ObsBin][dim]);
            }

            wt *= 389385730./binsize;
            if(IXsectUnits!=12){
               wt *= pow(10.,(IXsectUnits-12)) ;
            }

            // define the scale-bin number in the range  [0:nscalebin[
            int scalenode = 0;
            double cefscale[4] = {0.0,0.5,0.5,0.0};
            if(Nscalenode[0]>1){
               for(int i=1;i<Nscalenode[0];i++){
                  if (ScaleFac[0][scalevar]*scale1<ScaleNode[ObsBin][0][scalevar][i]){
                     scalenode=i;
                     break;
                  }
               }
               double deltascale = (ScaleFac[0][scalevar]*scale1 - ScaleNode[ObsBin][0][scalevar][scalenode-1])/
                  (ScaleNode[ObsBin][0][scalevar][scalenode]-ScaleNode[ObsBin][0][scalevar][scalenode-1]);
               vector<double> cscale(4) ; 
               cscale[0] = deltascale+1.0;
               cscale[1] = deltascale;
               cscale[2] = 1.0-deltascale;
               cscale[3] = 2.0-deltascale;

               if(scalenode<2 || scalenode>ScaleNode[ObsBin][0][scalevar][scalenode-1]-3){
                  cefscale[0] = 0.0;
                  cefscale[1] = 1.0-cscale[1];
                  cefscale[2] = 1.0-cscale[2];
                  cefscale[3] = 0.0;
               }else{
                  cefscale[1]=1.0-2.5*cscale[1]*cscale[1]+1.5*cscale[1]*cscale[1]*cscale[1];
                  cefscale[2]=1.0-2.5*cscale[2]*cscale[2]+1.5*cscale[2]*cscale[2]*cscale[2];
                  cefscale[0]=2.0 - 4.0*cscale[0] + 2.5*cscale[0]*cscale[0]
                     - 0.5*cscale[0]*cscale[0]*cscale[0];
                  cefscale[3]=2.0 - 4.0*cscale[3] + 2.5*cscale[3]*cscale[3]
                     - 0.5*cscale[3]*cscale[3]*cscale[3];
               }
            }

            // ** loop over all 4 scale nodes that receive contributions
            for(int i3 = 0; i3 < 4; i3++){
               int is = scalenode +i3 -2;  // the target scale index
               if(is<0){
                  is = 0;
               }
               if (is> Nscalenode[0]-1){
                  is = Nscalenode[0]-1;
               }
               // ** loop over all 4 x1 points that receive contributions
               for( int i1 = 0; i1 < 4; i1++) {           
                  for (int i2 = 0; i2 < 4; i2++){
                     int x1bin = (nx +i1 -1);
                     if(x1bin<0) x1bin = 0;
                     if(x1bin>Nxtot1[ObsBin]-1) x1bin =Nxtot1[ObsBin]-1;
                     int x2bin = (nx2 +i2 -1);
                     if(x2bin<0) x2bin = 0;
                     if(x2bin>Nxtot2[ObsBin]-1) x2bin =Nxtot2[ObsBin]-1;
                     int im = GetXIndex(ObsBin,x1bin,x2bin);
                     //                        if(i3==0) printf("fastNLO: filled at index %d in x2bin #%d x1bin #%d at x=%f\n",im,(nx2 +i2-1),(nx +i1 -1), XNode2[ObsBin][nx2+i2-1]);
                     for(int proc=0;proc<NSubproc;proc++){
                        //                           printf("%d %d %d %d %d %g %g\n",ObsBin,scalevar,is,im,proc,cefscale[i3],photoweight);
                        SigmaTilde[ObsBin][scalevar][is][im][proc] +=  cefm[i1]  * cefm2[i2] * photoweight * cefscale[i3] * wt[proc];
                     }
                  }
               }
            }
         }
      }
   }
}

#define EPS 1.0e-15

void fnloBlockBNlojet::_S_gauleg(unsigned int n, double *x, double *w)
{
   unsigned int m, j, i;
  double z1, z, pp, p3, p2, p1;
  
  m = (n+1)/2;
  for(i = 1; i <= m; i++) {
    z = cos(3.14159265358979323846*(i-0.25)/(n+0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for(j = 1; j <= n; j++) {
		p3 = p2;
		p2 = p1;
		p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp = n*(z*p1-p2)/(z*z-1.0);
      z1 = z;
      z = z1-p1/pp;
    } while(fabs(z-z1) > EPS);
    
    x[i-1] = -z, x[n-i] = z;
    w[i-1] = w[n-i] = 2.0/((1.0-z*z)*pp*pp);
  }
}



void fnloBlockBNlojet::FillEventHHC(int ObsBin, double x1, double x2, double scale1, const nlo::amplitude_hhc& amp, nlo::pdf_and_coupling_hhc& pdf, double prefactor){
   fnloBlockA2 *A2 =  BlockA2;

   // ---
   // --- Warm-Up Run to identify the extreme x,mu values
   // ---
   // KR: Add file output for later automatic read in
   if (IWarmUp == 1) {
     if (xlo[ObsBin] == 0.) xlo[ObsBin] = min(x1,x2);
     if (scalelo[ObsBin] == 0.) scalelo[ObsBin] = scale1;
     if (xlo[ObsBin] > min(x1,x2)) xlo[ObsBin] = min(x1,x2);
     if (scalelo[ObsBin] > scale1) scalelo[ObsBin] = scale1;
     if (scalehi[ObsBin] < scale1) scalehi[ObsBin] = scale1;
     IWarmUpCounter++;
     // KR: Create file already at start time, avoids problems when testing grid functionality  
     if ( IWarmUpCounter == 1 ) {
       FILE * ofile;
       ofile = fopen("fastNLO-warmup.dat","w");
       fprintf(ofile,"      // %lu contributions (!= events) in warm-up run \n",IWarmUpCounter);
       fclose(ofile);
     }
     if ( (IWarmUpCounter % IWarmUpPrint) == 0) {
       FILE * ofile;
       ofile = fopen("fastNLO-warmup.dat","w");
       printf("      // %lu contributions (!= events) in warm-up run \n",IWarmUpCounter);
       fprintf(ofile,"      // %lu contributions (!= events) in warm-up run \n",IWarmUpCounter);
       for (unsigned int i=0;i<BlockA2->GetNObsBin();i++){
         printf("      xlim[ %u ] = %e , mulo[ %u ] = %e , muup[ %u ] = %e ;\n",
		i,xlo[i],i,scalelo[i],i,scalehi[i]);
	 fprintf(ofile,"      xlim[ %u ] = %e , mulo[ %u ] = %e , muup[ %u ] = %e ;\n",
		 i,xlo[i],i,scalelo[i],i,scalehi[i]);
       }
       fclose(ofile);
     }
     return;
   }

   // --- select interpolation kernel for x and for mu 
   //             1:CatmulRom   2:Lagrangian
   const int ikernx = 2;     
   const int ikernmu = 2;    

   if(this->IRef>0){

      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){
         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         //nlo::weight_hhc wt = amp(pdf,mu2,mu2,prefactor);
         nlo::weight_hhc wtorg = amp(pdf,mu2,mu2,prefactor);
         nlo::weight_hhc wt;
	 wt[0] = wtorg[0];
	 wt[1] = wtorg[3];
	 wt[2] = wtorg[4];
	 wt[3] = wtorg[5];
	 wt[4] = wtorg[6];
	 wt[5] = wtorg[1];
	 wt[6] = wtorg[2];

	 if (NSubproc == 6) {
	   wt[5] += wt[6];    // -- sum is correct in ref-mode
	 }

         wt *= 389385730.;
         if(IXsectUnits!=12){
            wt *= pow(10.,(IXsectUnits-12)) ;
         }
         
         for(int proc=0;proc<NSubproc;proc++){
            SigmaTilde[ObsBin][scalevar][0][0][proc] += wt[proc];
         }
      }
   }else{
      if(this->NPDFDim != 1){
         printf("fnloBlockBNlojet::FillEventHHC: Error, only NPDFDim=1 (half matrix) implemented so far.\n");
         exit(1);
      }

      // Half matrix x variables
      double xmin=0.0, xmax=0.0;
      if (x1 > x2) {
         xmax = x1;
         xmin = x2; }
      else {
         xmax = x2;
         xmin = x1;
      }
      if (xmin<XNode1[ObsBin][0]){
         printf("fnloBlockBNlojet::FillEventHHC: find: xmin (%f) smaller than lowest x-node (%f) for bin #%d .\n",
                xmin,XNode1[ObsBin][0],ObsBin);
      }

      // **********  determine x_ij position in grid  ************
      // --- determine fractional contributions
      double hxmin  = -sqrt(-log10(xmin));
      double hxmax  = -sqrt(-log10(xmax));
      double hxone   = 0.0;

      // --- define the x-node numbers in the range: 0 <= nxnode < Nxtot1
      double hxlimit = Hxlim1[ObsBin];
      int nxmin = int(Nxtot1[ObsBin] *(hxmin-hxlimit)/(hxone-hxlimit));
      int nxmax = int(Nxtot1[ObsBin] *(hxmax-hxlimit)/(hxone-hxlimit));
      if (nxmin < 0) nxmin = 0;  // move into available range
      if (nxmax < 0) nxmax = 0;
      int nxminf = 0;  // - variables for final nodes (after modification in kernel)
      int nxmaxf = 0;

      // --- relative distances in h(xmin), h(xmax): deltamin,deltamax 
      double delta  = (hxone-hxlimit)/Nxtot1[ObsBin];
      double hxj = hxlimit+double(nxmin)/double(Nxtot1[ObsBin])*(hxone-hxlimit);
      double hxi = hxlimit+double(nxmax)/double(Nxtot1[ObsBin])*(hxone-hxlimit);
      double deltamin = (hxmin-hxj)/delta;
      double deltamax = (hxmax-hxi)/delta;

      // --- x-interpolation kernels and updated x-node positions: 1 <= nx < ntot-1
      vector<double> cefmin(4) ;
      vector<double> cefmax(4) ;
      Interpol(nxmin, (Nxtot1[ObsBin]-1), deltamin, ikernx, nxminf, &cefmin);
      Interpol(nxmax, (Nxtot1[ObsBin]-1), deltamax, ikernx, nxmaxf, &cefmax);

      // --- PDF reweighting - compute weights, modify cefmax[.], cefmin[.] 
      //     but only those within grid, there are no nodes at x=1
      double pdfwgtmax = PDFwgt(xmax);
      for( int i1 = 0; i1 < 4; i1++) {
	if ((nxmaxf-1+i1) >= 0 && (nxmaxf-1+i1) < Nxtot1[ObsBin] ) {
	  cefmax[i1] *= pdfwgtmax/PDFwgt(XNode1[ObsBin][nxmaxf-1+i1]);
	}
      }
      double pdfwgtmin = PDFwgt(xmin);
      for( int i2 = 0; i2 < 4; i2++) {
	if ((nxminf-1+i2) >= 0 && (nxminf-1+i2) < Nxtot1[ObsBin] ) {
	  cefmin[i2] *= pdfwgtmin/PDFwgt(XNode1[ObsBin][nxminf-1+i2]);
	} 
      }

      // === the weights for the bi-cubic eigenfunctions (2-dim)
      double bicef[4][4];      
      for( int i1 = 0; i1 < 4; i1++) {
         for( int i2 = 0; i2 < 4; i2++) {
            bicef[i1][i2] = cefmax[i1] * cefmin[i2];
         }
      }

      // --- loop over all scale variations
      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){

	 // --- compute renormalization=factorization scale squared 
         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         //nlo::weight_hhc wt = amp(pdf,mu2,mu2,prefactor);
	 nlo::weight_hhc wtorg = amp(pdf,mu2,mu2,prefactor);
	 // - rearrange subprocesses
	 nlo::weight_hhc wt;
	 wt[0] = wtorg[0];
	 wt[1] = wtorg[3];
	 wt[2] = wtorg[4];
	 wt[3] = wtorg[5];
	 wt[4] = wtorg[6];
	 wt[5] = wtorg[1];
	 wt[6] = wtorg[2];
	 // -- case NSubproc=6: see below

         wt *= 389385730.;
         if(IXsectUnits!=12){
            wt *= pow(10.,(IXsectUnits-12)) ;
         }

         // deal with subprocesses 2 and 3 -> moved to 6 and 7
         //    - if x1>x2 -> o.k.
         //    - if x2>x1 -> swap weights for subprocesses 2,3 -> now 6,7
         if(x2>x1){
            double buffer;
            //buffer = wt[1];
            //wt[1] = wt[2];
            //wt[2] = buffer;
            buffer = wt[5];
            wt[5] = wt[6];
            wt[6] = buffer;
         }
	 // --- combine subprocesses 5,6 here after possible swapping
	 if (NSubproc == 6) {   
	   wt[5] = (wt[5]+wt[6])/2.;
	   wt[6] = wt[5];    // set equal in case of swapping below diagonal (at bottom) 
	 }

         // --- define the scale-bin number in the range:  0 <= scalenode < nscalebin-2
         int scalenode = Nscalenode[0]-2;  // --- initialize with largest possible value
         int nscale = 0;                   // --- variable for final scale node
         double mu0scale = 0.25;           // --- parameter in transformation function H(mu)
         vector<double> cefscale(4);
         cefscale[0] = 0.0;
         cefscale[1] = 0.5;
         cefscale[2] = 0.5;
         cefscale[3] = 0.0;

         if(Nscalenode[0]>1){
            // --- find scale position in range:  0 <= scalenode < nscalenode-1
            for(int i=1;i<Nscalenode[0];i++){
               if (ScaleFac[0][scalevar]*scale1<ScaleNode[ObsBin][0][scalevar][i]){
                  scalenode=i-1;
                  break;
               }
            }

	    // --- relative distance delta - in function H(mu)
            double deltascale = (log(log(ScaleFac[0][scalevar]*scale1/mu0scale)) 
				 - HScaleNode[ObsBin][0][scalevar][scalenode])/
	      (HScaleNode[ObsBin][0][scalevar][scalenode+1]-HScaleNode[ObsBin][0][scalevar][scalenode]);

            // --- get scale interpolation kernel and updated scalenode position: 1 <= nmu < ntot-2
            Interpol(scalenode, (Nscalenode[0]-2), deltascale, ikernmu, nscale, &cefscale);
         }

         // --- loop over all 4 scale nodes that receive contributions
         for(int i3 = 0; i3 < 4; i3++){
            int is = nscale + (i3-1);  // --- the target scale index

            // --- check: all elements should end up within grid
            if (is < 0){
	      printf(" is<0     %d %d %f \n",nscale,i3,cefscale[i3]);
	      exit(1);
            } else if (is > Nscalenode[0]-1){
	      printf(" is>max   %d %d %f \n",nscale,i3,cefscale[i3]);
	      exit(1);
            }
            // --- loop over all 16 xmin,xmax points that receive contributions
            for( int i1 = 0; i1 < 4; i1++) {           
               for (int i2 = 0; i2 < 4; i2++){
		  int xmaxbin = nxmaxf + (i1 -1);   // - the target x-max index
                  int xminbin = nxminf + (i2 -1);   // - the target x-min index
		  if (xmaxbin < Nxtot1[ObsBin] && xminbin < Nxtot1[ObsBin]) {
		    nlo::weight_hhc wtmp = wt;  // - a working copy of the weights
		    // --- check if above diagonal? project back and swap qg<->gq
		    if (xminbin>xmaxbin) { 
		      int di = xminbin - xmaxbin;
		      xmaxbin = xmaxbin + di;   // modify indicees
		      xminbin = xminbin - di;		
		      double buffer;
		      //buffer = wtmp[1]; // swap subprocesses 2,3
		      //wtmp[1] = wtmp[2];
		      //wtmp[2] = buffer;
		      // new ordering (1,2)->(5,6)	 
		      buffer  = wtmp[5]; // swap subprocesses 6,7
		      wtmp[5] = wtmp[6];
		      wtmp[6] = buffer;
		    }
		    int im = GetXIndex(ObsBin,xminbin,xmaxbin);
		    // printf("fastNLO: index %d in xmaxbin #%d xminbin #%d\n",im,xmaxbin,xminbin);
		    for(int proc=0;proc<NSubproc;proc++){
		      // printf("%d %d %d %d %d %g %g\n",ObsBin,scalevar,is,im,proc,bicef[i1][i2],cefscale[i3]);
		      SigmaTilde[ObsBin][scalevar][is][im][proc] +=  bicef[i1][i2] * cefscale[i3] * wtmp[proc];
		    }
		  }
               }
            }
         }
      }
   }
}


//________________________________________________________________________________________________________________ //


void fnloBlockBNlojet::Interpol(int nnode, int nmax, double delta, int ikern, int &nmod, vector < double > *kernel){
  // nnode: number of the next node to the left of the current value
  // nmax:  number of the last node which could lie to the left of a potential value
  // delta: relative distance of value to node 'nnode'
  // ikern: select interpolation kernel   1:Catmul Rom  2: Lagrange
  // nmod:  modified number of next node to the left (to be used for storage - relevant only at boundaries)
  // kernel: array(4) containing the interpolation kernel

  // --- distances to all nodes
  vector <double> dist(4);
  dist[0] = 1. + delta;
  dist[1] = 0. + delta;
  dist[2] = 1. - delta;
  dist[3] = 2. - delta;

  vector<double> &kern = *kernel;   // --- to make acces via "[]" more easy

  if (ikern == 1) {         // --- Catmul Rom interpolation kernel
    if (nnode == 0 ) { // --- left boundary
      kern[0] = 1.0 - 7.0/6.0*delta - 1.0/6.0*delta*delta + 1.0/3.0*delta*delta*delta;
      kern[1] = 4.0/3.0*delta + 1.0/3.0*delta*delta - 2.0/3.0*delta*delta*delta;
      kern[2] = -1.0/6.0*delta - 1.0/6.0*delta*delta + 1.0/3.0*delta*delta*delta;
      kern[3] = 0.0;
      nmod = nnode + 1;
    } else if (nnode == nmax) { // --- right boundary
      kern[0] = 0.0;
      kern[1] = -1.0/6.0*dist[2] - 1.0/6.0*dist[2]*dist[2] + 1.0/3.0*dist[2]*dist[2]*dist[2];
      kern[2] = 4.0/3.0*dist[2] + 1.0/3.0*dist[2]*dist[2] - 2.0/3.0*dist[2]*dist[2]*dist[2];
      kern[3] = 1.0 -7.0/6.0*dist[2] -1.0/6.0*dist[2]*dist[2] +1.0/3.0*dist[2]*dist[2]*dist[2];
      nmod = nnode - 1;
    } else { // --- central region
      kern[0] = 2.0 - 4.0*dist[0] + 2.5*dist[0]*dist[0] - 0.5*dist[0]*dist[0]*dist[0];
      kern[1] = 1.0 - 2.5*dist[1]*dist[1] + 1.5*dist[1]*dist[1]*dist[1];
      kern[2] = 1.0 - 2.5*dist[2]*dist[2] + 1.5*dist[2]*dist[2]*dist[2];
      kern[3] = 2.0 - 4.0*dist[3] + 2.5*dist[3]*dist[3] - 0.5*dist[3]*dist[3]*dist[3];
      nmod = nnode;
    }

  } else if (ikern == 2) {  // --- Lagrange interpolation kernel
    if (nnode == 0 ) { // --- left boundary
      kern[0] = 1.0 - 11./6.*delta + delta*delta - 1./6.*delta*delta*delta;
      kern[1] = 3.0*delta - 2.5*delta*delta + 0.5*delta*delta*delta;
      kern[2] = -1.5*delta + 2.0*delta*delta - 0.5*delta*delta*delta;
      kern[3] = 1./3.*delta - 0.5*delta*delta + 1./6.*delta*delta*delta;
      nmod = nnode + 1;
    } else if (nnode == nmax) { // --- right boundary
      kern[0] = 1./3.*dist[2] - 0.5*dist[2]*dist[2] + 1./6.*dist[2]*dist[2]*dist[2];
      kern[1] = -1.5*dist[2] + 2.0*dist[2]*dist[2] - 0.5*dist[2]*dist[2]*dist[2];
      kern[2] = 3.0*dist[2] - 2.5*dist[2]*dist[2] + 0.5*dist[2]*dist[2]*dist[2];
      kern[3] = 1.0 - 11./6.*dist[2] + dist[2]*dist[2] - 1./6.*dist[2]*dist[2]*dist[2];
      nmod = nnode - 1;
    } else { // --- central region - 
      kern[0] = 1.0 - 11./6.*dist[0] + 1.*dist[0]*dist[0] - 1./6.*dist[0]*dist[0]*dist[0];
      kern[1] = 1.0 - 0.5*dist[1] - 1.0*dist[1]*dist[1] + 0.5*dist[1]*dist[1]*dist[1];
      kern[2] = 1.0 - 0.5*dist[2] - 1.0*dist[2]*dist[2] + 0.5*dist[2]*dist[2]*dist[2];
      kern[3] = 1.0 - 11./6.*dist[3] + 1.*dist[3]*dist[3] - 1./6.*dist[3]*dist[3]*dist[3];
      nmod = nnode;
    }

  } else {               // --- unknown interpolation kernel
    printf("fnloBlockBNlojet::Interpol: Error  interpolation kernel (%d) is not defined.\n",ikern);
    exit(1);
  }
}


//________________________________________________________________________________________________________________ //



vector < double > fnloBlockBNlojet::Interpol(int nnode, int nmax, double delta, EInterpolKernel eIkern, int &nmod ){
   // simplified call for ::Interpol 
   // 
   // input:
   //    nnode		number of the next node to the left of the current value
   //    nmax		number of the last node which could lie to the left of a potential value
   //    delta		relative distance of value to node 'nnode'
   //    eIkern		interpolation kernel
   // 'returns' :
   //     nmod		modified number of next node to the left (to be used for storage - relevant only at boundaries)
   //     kern  array(4) containing the interpolation kernel (pleas pass a suitable array)

   // return vector
   vector < double > kern(4);

   if ( eIkern == kCatmulRom) {         // --- Catmul Rom interpolation kernel
      Interpol(nnode, nmax, delta, 1, nmod , &kern );
   } 

   else if ( eIkern == kLagrangian ) {  // --- Lagrange interpolation kernel
      Interpol(nnode, nmax, delta, 2, nmod , &kern );
   } 
   
   else {               // --- unknown interpolation kernel
      printf("fnloBlockBNlojet::Interpol: Error  interpolation kernel (%d) is not defined.\n",eIkern);
      exit(1);
   }
   
   return kern;
}


//________________________________________________________________________________________________________________ //



void fnloBlockBNlojet::InitDISConstants( fnloBlockA2* A2 , bool nlo ){
   // -------------------------------------------------------------------------- //
   //  
   //  InitDISConstants(). Method for initalizing all necessary fastNLO values
   //  for a reasonable v2.1 table.
   //
   //  This method is only for v2.1 tables
   //  
   //
   // -------------------------------------------------------------------------- //
   

   IXsectUnits	= A2->Ipublunits;
   IDataFlag = 0;
   IAddMultFlag = 0;
   IContrFlag1 = 1;
   //   IContrFlag3 = 0;
   CodeDescript.push_back("NLOJET++ 4.1.3");  // --- fastNLO user: enter NLOJET++ version
   IRef = 0;
   
   if (nlo || A2->ILOord > 1) {
      NSubproc = 3;
      printf("  this job uses 3 subprocesses \n");
   } else {
      NSubproc = 2;
      printf("  this job uses 2 subprocesses \n");
   }
   
   if(nlo){
      CtrbDescript.push_back("NLO");
      IContrFlag2 = 2;
      IScaleDep = 1;
      Npow = A2->ILOord+1;
   }else{
      CtrbDescript.push_back("LO");      
      IContrFlag2 = 1;
      IScaleDep = 0;
      Npow = A2->ILOord;
   }
  
   NPDF = 1;
   NPDFPDG.push_back(2212);	// --- fastNLO user: PDG code for hadron
   NPDFDim = 0;			// DIS
   NFragFunc = 0;
   IPDFdef1 = 2;
   IPDFdef2 = 1;

   if(NSubproc == 3) {
      IPDFdef3 = 3;
   } else {
      IPDFdef3 = 2;
      printf("  set IPDFdef3=2 consistent with 2 subprocesses \n");
   }

   IWarmUp = 0;			// no warm-up run -> production run.
   IWarmUpPrint = 10000000 ;

   NScaleDep = 3;
   
   // ---- resize the vectors ---- //
   XNode1.resize(A2->NObsBin);

   scale1lo.resize(A2->NObsBin);
   scale1hi.resize(A2->NObsBin);
   scale2lo.resize(A2->NObsBin);
   scale2hi.resize(A2->NObsBin);

   // ---- those numbers are partly not perfectly defined ---- //
   NScales = 1;		// 
   Iscale.resize(1);
   Iscale[0] = 0;		// mur=mur(ET), ET = index 0 

   NScaleDim = 1;		// NEVER SET NScaleDim TO ANY OTHER VALUE THAN 1 !!!
   ScaleDescript.resize(1);
   ScaleDescript[0].resize(2);
   ScaleDescript[0][0] = ("Q");
   ScaleDescript[0][1] = ("pt");
   
   if ( NScaleDep != 3 ){
      ScaleFac.resize(1);	// 1 = NScaleDim
   }

}


//________________________________________________________________________________________________________________ //



void fnloBlockBNlojet::InitLHCConstants( fnloBlockA2* A2 , bool nlo ){
   // -------------------------------------------------------------------------- //
   //  
   //  InitDISConstants(). Method for initalizing all necessary fastNLO values
   //  for a reasonable v2.1 table.
   //
   //  This method is only for v2.1 tables
   //  
   //
   // -------------------------------------------------------------------------- //
   
   CodeDescript.push_back("NLOJET++ 4.1.3");  // --- fastNLO user: enter NLOJET++ version

   IXsectUnits	= A2->Ipublunits;
   IDataFlag	= 0;
   IAddMultFlag	= 0;
   IContrFlag1	= 1;
   IContrFlag3	= 0;
   IRef		= 0;
   // -> v2.0 
   // NSubproc	= (nlo || A2->ILOord > 2) ? 7 : 6;
   // printf("  this job uses %d subprocesses \n",NSubproc);
   NSubproc	= 7;
   printf("  this job uses ALWAYS (also for LO) %d subprocesses \n",NSubproc);
   
   if(nlo){
      CtrbDescript.push_back("NLO");
      IContrFlag2 = 2;
      IScaleDep = 1;
      Npow = A2->ILOord+1;
   }else{
      CtrbDescript.push_back("LO");      
      IContrFlag2 = 1;
      IScaleDep = 0;
      Npow = A2->ILOord;
   }
  
   NPDF		= 2;
   NPDFPDG.push_back(2212);	// --- fastNLO user: PDG code for 1st hadron
   NPDFPDG.push_back(2212);	// --- fastNLO user: PDG code for 2nd hadron
   NPDFDim	= 1;			// pp
   NFragFunc	= 0;
   IPDFdef1	= 3;
   IPDFdef2	= 1;
   IPDFdef3	= NSubproc == 7 ? 2 : 1;
   printf("         Set IPDFdef3 = %d, consistent with %d subprocesses.\n",IPDFdef3,NSubproc);
   

   //IWarmUp	= 0;			// no warm-up run -> production run.
   IWarmUpPrint	= 10000000 ;

   NScaleDep	= 3;
   
   // ---- resize the vectors ---- //
   XNode1.resize(A2->NObsBin);

   scale1lo.resize(A2->NObsBin);
   scale1hi.resize(A2->NObsBin);
   scale2lo.resize(A2->NObsBin);
   scale2hi.resize(A2->NObsBin);

   // ---- those numbers are partly not ambigously defined in v2.1 ---- //
   NScales = 1;		// 
   Iscale.resize(1);
   Iscale[0] = 0;		// mur=mur(ET), ET = index 0 

   NScaleDim = 1;		// NEVER SET NScaleDim TO ANY OTHER VALUE THAN 1 !!!
   ScaleDescript.resize(1);
   ScaleDescript[0].resize(2);

}


//________________________________________________________________________________________________________________ //



void fnloBlockBNlojet::InitFinalDISValues( fnloBlockA2* A2 , double* xlim , double* scale1lo , double* scale1hi , double* scale2lo , double* scale2hi ){
   // -------------------------------------------------------------------------- //
   //  
   //  InitFinalDISValues(). Method for initalizing all necessary fastNLO values
   //  for a reasonable v2.1 table.
   //
   //   this method could be used for v2.1 and for v2.0 tables
   //   for v2.0 tables just pass scale1lo and scale1hi arrays
   //  
   //
   // -------------------------------------------------------------------------- //
   
   if ( NScaleDim != 1 ) cout << "Error! NScaleDim is not supposed to be <= one." << endl;

   if ( NScaleDep != 3 ){
      const double mu0scale = 0.25; // --- variable in H(mu) (in GeV)
      Nscalevar.push_back(ScaleFac[0].size());

      ResizeTable( &ScaleNode  , A2->NObsBin , NScaleDim , Nscalevar[0] , Nscalenode[0] );
      ResizeTable( &HScaleNode , A2->NObsBin , NScaleDim , Nscalevar[0] , Nscalenode[0] );
     
      // init ScaleNode and HScaleNode
      for(int i=0;i<A2->NObsBin;i++){
	 int j = 0; // this was once the NScaleDim loop...
	 for(int k=0;k<Nscalevar[j];k++){ 
	    if(Nscalenode[j]==1){
	       ScaleNode  [i][j][k][0]  = ScaleFac[0][k]*(scale1hi[i]+scale1lo[i])/2.;
	       HScaleNode [i][j][k][0] = log(log((ScaleFac[0][k]*(scale1hi[i]+scale1lo[i])/2.)/mu0scale));
	    }else{
	       double llscale1lo = log(log((ScaleFac[0][k]*scale1lo[i])/mu0scale));
	       double llscale1hi = log(log((ScaleFac[0][k]*scale1hi[i])/mu0scale));
	       for(int l=0;l<Nscalenode[j];l++){
		  // 1) later this is the place where the Chebychev nodes will be implemented
		  // 2) here also: llscale1lo is supposed to be scale 1
		  // 3) ScaleNode'X' is for 2ScaleInterpolationTables (NScaleDep==2)
		  // 4) The one without number is for crosschecks and that FNLO is still working (NScaleDep == 1)
		  // 		HScaleNode [i][j][k][l] = llscalelo +  double(l)/double(Nscalenode[j]-1)*(llscalehi-llscalelo);
		  // 		ScaleNode [i][j][k][l] = mu0scale * exp(exp(HScaleNode [i][j][k][l]));
		  HScaleNode [i][j][k][l] = llscale1lo +  double(l)/double(Nscalenode[j]-1)*(llscale1hi-llscale1lo);
		  ScaleNode [i][j][k][l] = mu0scale * exp(exp(HScaleNode [i][j][k][l]));
	       }
	    }
	 }
      }

      // ---- init sigma tilde ---- //
      int XmaxFromI[1] = {0};
      ResizeTable( &SigmaTilde , A2->NObsBin , Nscalevar[0] , Nscalenode[0] , XmaxFromI , NSubproc );

   } // NScaleDep != 3 (v2.0 tables)
   else if ( NScaleDep == 3 ){
      if ( scale2hi==NULL ) printf("Error.\n");
      InitLogLogScaleNode( A2 , scale1lo , scale1hi , 1 );
      InitLogLogScaleNode( A2 , scale2lo , scale2hi , 2 );
      ResizeSigmaTildeTables( A2 );
   }   

}


//________________________________________________________________________________________________________________ //


void fnloBlockBNlojet::InitScaleNode( fnloBlockA2* A2, double* slo , double* shi, int iScale  ){
   //
   // InitScaleNode
   // call this function after initializing Fct_H_Scale[iScale-1]
   // 
   // This function initializes and resizes
   //    HScaleNode
   //    ScaleNode
   //    

   // ---- get 'universal' pointers ---- //
   vector < vector < double > >* node	= (iScale == 1) ? &ScaleNode1 : &ScaleNode2;
   vector < vector < double > >* Hnode	= (iScale == 1) ? &HScaleNode1 : &HScaleNode2;
   int Nscalenodes			= (iScale == 1) ? NscalenodeScale1 : NscalenodeScale2 ;

   // ---- init scale nodes ---- //
   ResizeTable( node   , A2->NObsBin , Nscalenodes );
   ResizeTable( Hnode  , A2->NObsBin , Nscalenodes );
   
   // ---- make binning ---- //
   for(int i=0;i<A2->NObsBin;i++){
      double llslo = (this->*Fct_H_Scale[iScale-1])(slo[i]);
      double llshi = (this->*Fct_H_Scale[iScale-1])(shi[i]);
      for(int l=0;l<Nscalenodes;l++){ 
	 (*Hnode)[i][l]   = llslo +  double(l)/double(Nscalenodes-1)*(llshi-llslo);
	 (*node) [i][l]   = (this->*Fct_H_Scale_Inv[iScale-1])((*Hnode)[i][l]);
	 //printf("bin %3d: Hnode[%3d][%2d] = %7.4f, ScaleNode[%3d][%2d] = %7.4f\n",i,i,l,(*Hnode)[i][l],i,l,(*node) [i][l]);
      }
   }
}


//________________________________________________________________________________________________________________ //


void fnloBlockBNlojet::InitLogLogScaleNode( fnloBlockA2* A2 , double* slo , double* shi , int iScale ){
   // -------------------------------------------------------------------------- //
   //  
   //  InitLogLogScaleNode(). Initialize scale nodes in a loglog distance 
   //  
   //  Input values:
   //       slo, shi:	arrays from warm-up run (extreme values)
   //       iScale:		Scale 1 or 2
   //
   // -------------------------------------------------------------------------- //

   // ---- check validity of call ---- //
   if ( iScale != 1 && iScale != 2 ){
      printf("fnloBlockBNlojet::InitLogLogScaleNode. Error. iScale must be 1 or 2.\n");
      exit(1);
   }
   if ( NScaleDim != 1 ) cout << "Error! NScaleDim is not supposed to be <= one." << endl;
   if ( NScaleDep != 3 ){
      printf("fnloBlockBNlojet::InitLogLogScaleNode. Error. This method works sofar only for v2.1 tables.\n");
   } 

   // ---- init pointers to correspoding functions ---- //
   Fct_H_Scale[iScale-1]	= &fnloBlockBNlojet::Function_loglog025;
   Fct_H_Scale_Inv[iScale-1]	= &fnloBlockBNlojet::Function_loglog025_inv;

   // ---- init the scale nodes ---- //
   InitScaleNode( A2, slo, shi,  iScale ) ;   

}


//________________________________________________________________________________________________________________ //


void fnloBlockBNlojet::InitLinearScaleNode( fnloBlockA2* A2 , double* slo , double* shi , int iScale ){
   // -------------------------------------------------------------------------- //
   //  
   //  InitLogLogScaleNode(). Initialize scale nodes in a linear distance 
   //  
   //  Input values:
   //       slo, shi:	arrays from warm-up run (extreme values)
   //       iScale:		ScaleNode 1 or 2
   //
   // -------------------------------------------------------------------------- //

   // ---- check validity of call ---- //
   if ( iScale != 1 && iScale != 2 ){
      printf("fnloBlockBNlojet::InitLogLogScaleNode. Error. iScale must be 1 or 2.\n");
      exit(1);
   }
   if ( NScaleDim != 1 ) cout << "Error! NScaleDim is not supposed to be <= one." << endl;
   if ( NScaleDep != 3 ){
      printf("fnloBlockBNlojet::InitLogLogScaleNode. Error. This method works sofar only for v2.1 tables.\n");
   } 

   // ---- init pointers to correspoding vectors ---- //
   Fct_H_Scale[iScale-1]	= &fnloBlockBNlojet::Function_x;
   Fct_H_Scale_Inv[iScale-1]	= &fnloBlockBNlojet::Function_x_inv;

   // ---- init the scale nodes ---- //
   InitScaleNode( A2, slo, shi, iScale ) ;   
}


//________________________________________________________________________________________________________________ //


double fnloBlockBNlojet::Function_loglog025( double mu ){
   // function H(mu) = log(log( mu / 0.25 ))
   return log(log(mu/0.25));}
double fnloBlockBNlojet::Function_loglog025_inv( double mu ){
   // inverse of function H(mu) = log(log( mu / 0.25 ))
   return 0.25*exp(exp(mu));}

double fnloBlockBNlojet::Function_x( double mu ){
   // function H(mu) = x
   return mu;}
double fnloBlockBNlojet::Function_x_inv( double mu ){
   // inverse of function H(mu) = x;
   return mu;}

double fnloBlockBNlojet::Function_log10( double x ){
   return log10(x);}
double fnloBlockBNlojet::Function_log10_inv( double x ){
   return pow(10,x);}

double fnloBlockBNlojet::Function_sqrtlog10( double x ){
   return -sqrt(-log10(x));}
double fnloBlockBNlojet::Function_sqrtlog10_inv( double x ){
   return pow(10,-pow(x,2));}


//________________________________________________________________________________________________________________ //


void fnloBlockBNlojet::ResizeSigmaTildeTables( fnloBlockA2* A2 ){
   int XmaxFromI[1] = {0};
   // ---- init sigma tilde for the subprocess dependent table  ---- //
   ResizeTable( &SigmaTildeMuIndep , A2->NObsBin , XmaxFromI , NscalenodeScale1 , NscalenodeScale2 , NSubproc );
   ResizeTable( &SigmaTildeMuFDep  , A2->NObsBin , XmaxFromI , NscalenodeScale1 , NscalenodeScale2 , NSubproc );
   ResizeTable( &SigmaTildeMuRDep  , A2->NObsBin , XmaxFromI , NscalenodeScale1 , NscalenodeScale2 , NSubproc );
   
   ResizeTable( &SigmaRefMixed     , A2->NObsBin , NSubproc );
   ResizeTable( &SigmaRef_s1       , A2->NObsBin , NSubproc );
   ResizeTable( &SigmaRef_s2       , A2->NObsBin , NSubproc );
}


//________________________________________________________________________________________________________________ //


void fnloBlockBNlojet::InitReferenceTable( fnloBlockA2* A2 ){
   // -------------------------------------------------------------------------- //
   //  
   //  InitReferenceTable(). Method for initalizing all necessary fastNLO values
   //  for a reasonable v2.0 and v2.1 reference table.
   //
   // -------------------------------------------------------------------------- //

   // do I need sth. like this here?
   //          // refB->NSubproc = 3;
   //       if (nlo || A2->ILOord > 1) {
   //         refB->NSubproc = 3;
   //       } else {
   //         refB->NSubproc = 2;
   //         printf("  this reference job uses 2 subprocesses \n");
   //       }

   IRef = 1; // IS Reference
   Nscalenode[0] = 1;
   Nxtot1.clear();
   Hxlim1.clear();
   NScaleDep = 1; //! the 2-scale-interplation mode has automatic reference cross sections.
   for(int i=0;i<A2->NObsBin;i++){
      Nxtot1.push_back(1);
      Hxlim1.push_back(0.);
      XNode1[i].clear(); 
      XNode1[i].push_back(0.); 
   }
   
   int XmaxFromI[1] = {0};
   ResizeTable( &ScaleNode  , A2->NObsBin , NScaleDim , Nscalevar[0] , Nscalenode[0] );
   ResizeTable( &HScaleNode , A2->NObsBin , NScaleDim , Nscalevar[0] , Nscalenode[0] );
   ResizeTable( &SigmaTilde , A2->NObsBin , Nscalevar[0] , Nscalenode[0] , XmaxFromI , NSubproc );
      
}


//________________________________________________________________________________________________________________ //



void fnloBlockBNlojet::SetScale1Name( string name ){
   if ( ScaleDescript.empty() ) printf("fnloBlockBNlojet::SetScale1Name. Error.\n");
   if ( ScaleDescript[0].empty() ) printf("fnloBlockBNlojet::SetScale1Name. Error.\n");
   if ( ScaleDescript[0].size() != 2 ) printf("fnloBlockBNlojet::SetScale1Name. Error.\n");
   ScaleDescript[0][0]	 = name;  
}

void fnloBlockBNlojet::SetScale2Name( string name ){
   if ( ScaleDescript.empty() ) printf("fnloBlockBNlojet::SetScale2Name. Error.\n");
   if ( ScaleDescript[0].empty() ) printf("fnloBlockBNlojet::SetScale2Name. Error.\n");
   if ( ScaleDescript[0].size() != 2 ) printf("fnloBlockBNlojet::SetScale2Name. Error.\n");
   ScaleDescript[0][1]	 = name;  
}


//________________________________________________________________________________________________________________ //



void fnloBlockBNlojet::SetNumberOfXNodesPerMagnitude( int nxPerMagnitude , double* xlim ){
   // -------------------------------------------------------------------------- //
   //  
   //  Set number of x-noder per order of magnitude in x.
   //
   //  This method is tested only for v2.1 tables
   //  XNode1 must be 'resized' to the right number of ObsBins before usage.
   //  This function makes use of Fct_H_XNode, which is implemented consistently 
   //  only for v2.1 tables yet.
   //
   //  Input
   //     - nxPerMagnitude	number of x-nodes for each order
   //	  - xlim	array (from warm-up run) holding lowest x-limit in each bin
   //  
   //
   // -------------------------------------------------------------------------- //
   printf("   fnloBlockBNlojet::SetNumberOfXNodesPerMagnitude(). Info. Number of x-nodes in each ObsBin (%d per order of magnitude, or at least %d nodes ):\n  *  ",nxPerMagnitude,nxPerMagnitude);

   // set functions how the x-nodes are binned
   if ( NPDF == 1) {
      Fct_H_XNode		= &fnloBlockBNlojet::Function_log10;
      Fct_H_XNode_Inv		= &fnloBlockBNlojet::Function_log10_inv;
   }
   else if ( NPDF == 2 ){
      Fct_H_XNode		= &fnloBlockBNlojet::Function_sqrtlog10;
      Fct_H_XNode_Inv		= &fnloBlockBNlojet::Function_sqrtlog10_inv;
   }

   for(int i=0;i<XNode1.size();i++){
      if ( (xlim[i] < 1.e-6 || xlim[i]>1) && IWarmUp == 0 ) {
	 printf("fnloBlockBNlojet::SetNumberOfXNodesPerMagnitude. Warning. You might have not initialized your xlim-array properly. Please do this before calling SetNumberOfXNodesPerMagnitude. xlim[%d] = %8.3e, IWarmUp = %d.\n",i,xlim[i],IWarmUp);
      }
      int nxtot	= (int)(fabs(log10(xlim[i]))*nxPerMagnitude);
      if ( nxtot < nxPerMagnitude ) nxtot = nxPerMagnitude; // at least nxPerMagnitude points

      // printing...
      printf("%d: %d, ",i,nxtot);
      if ( i==XNode1.size()-1 )printf("\n");
      
      // init x-nodes
      Nxtot1.push_back(nxtot);
      double hxlim = (this->*Fct_H_XNode)(xlim[i]);
      Hxlim1.push_back(hxlim);
      for(int j=0;j<nxtot;j++){
         double hx = hxlim*( 1.- ((double)j)/(double)nxtot);
	 XNode1[i].push_back((this->*Fct_H_XNode_Inv)(hx));
      }
   }
}
