#include <cstdlib>
#include <iostream>

#include "fnloBlockBNlojet.h"

//
// note (Oct 8, 2010) - MW:
// the two routines FillEventDIS, and FillEventHHC are running
// the other routines (DIS2scale and photoproduction) still need some work
//


void fnloBlockBNlojet::FillEventDIS(int ObsBin, double x, double scale1, const nlo::amplitude_dis& amp, nlo::pdf_and_coupling_dis& pdf, double prefactor){
   fnloBlockA2 *A2 =  BlockA2;

   // ---
   // --- Warm-Up Run to identify extreme x,mu values
   // ---
   // KR: Add file output for later automatic read in
   if (IWarmUp == 1) {
     if (xlo[ObsBin] == 0.) xlo[ObsBin] = x;
     if (scalelo[ObsBin] == 0.) scalelo[ObsBin] = scale1;
     if (xlo[ObsBin] > x) xlo[ObsBin] = x;
     if (scalelo[ObsBin] > scale1) scalelo[ObsBin] = scale1;
     if (scalehi[ObsBin] < scale1) scalehi[ObsBin] = scale1;
     IWarmUpCounter++;
     // KR: Create file already at start time, avoids problems when testing grid functionality  
     if ( IWarmUpCounter == 1 ) {
       FILE * ofile;
       ofile = fopen("fastNLO-warmup.dat","w");
       fprintf(ofile," // %lu contributions (!= events) in warm-up run \n",IWarmUpCounter);
       fclose(ofile);
     }
     if ( (IWarmUpCounter % IWarmUpPrint) == 0) {
       FILE * ofile;
       ofile = fopen("fastNLO-warmup.dat","w");
       printf(" // %lu contributions (!= events) in warm-up run \n",IWarmUpCounter);
       fprintf(ofile," // %lu contributions (!= events) in warm-up run \n",IWarmUpCounter);
       for (int i=0;i<BlockA2->GetNObsBin();i++){
         printf(" xlim[%d]=%7.5f, mulo[%d]=%8.3f, muup[%d]=%8.3f;\n",i,xlo[i],i,scalelo[i],i,scalehi[i]);
	 //	 fprintf(ofile,"%7.5f, %8.3f, %8.3f;\n",xlo[i],scalelo[i],scalehi[i]);
	 fprintf(ofile,"%e, %e, %e;\n",xlo[i],scalelo[i],scalehi[i]);
       }
       fclose(ofile);
     }
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
void fnloBlockBNlojet::FillEventDIS2Scale(int ObsBin, double x, double scale1, double scale2,  const nlo::amplitude_dis& amp, nlo::pdf_and_coupling_dis& pdf, double prefactor){
   fnloBlockA2 *A2 =  BlockA2;

   if(this->IRef>0){

      for(int scalevar=0; scalevar<GetTotalScalevars(); scalevar++){
         int scalei1 = scalevar / Nscalevar[1];
         int scalei2 = scalevar % Nscalevar[1];

         double mur2 = ScaleFac[0][scalei1]*ScaleFac[0][scalei1]*scale1*scale1;
         double muf2 = ScaleFac[1][scalei2]*ScaleFac[1][scalei2]*scale2*scale2;
         //         printf("%d %d %d %f %f %f %f\n",scalevar,scalei1,scalei2,ScaleFac[0][scalei1],ScaleFac[1][scalei2],mur2,muf2);
         nlo::weight_dis wt = amp(pdf,mur2,muf2,prefactor);

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
         printf("fnloBlockB::FillEventDIS: Error: x (%f) smaller than lowest point (%f) at bin #%d .\n",
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

      for(int scalevar=0; scalevar<GetTotalScalevars(); scalevar++){

         int scalei1 = scalevar / Nscalevar[1];
         int scalei2 = scalevar % Nscalevar[1];

         double mur2 = ScaleFac[0][scalei1]*ScaleFac[0][scalei1]*scale1*scale1;
         double muf2 = ScaleFac[1][scalei2]*ScaleFac[1][scalei2]*scale2*scale2;

         nlo::weight_dis wt = amp(pdf,mur2,muf2,prefactor);

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

         wt *= 389385730./binsize;
         if(IXsectUnits!=12){
            wt *= pow(10.,(IXsectUnits-12)) ;
         }

         // define the scale-bin for dimension 1
         int scalenode1 = 0;
         double cefscale1[4] = {0.0,0.5,0.5,0.0};
         if(Nscalenode[0]>1){
            for(int i=1;i<Nscalenode[0];i++){
               if (ScaleFac[0][scalei1]*scale1<ScaleNode[ObsBin][0][scalei1][i]){
                  scalenode1=i;
                  break;
               }
            }
            double deltascale = (ScaleFac[0][scalei1]*scale1 - ScaleNode[ObsBin][0][scalei1][scalenode1-1])/
               (ScaleNode[ObsBin][0][scalei1][scalenode1]-ScaleNode[ObsBin][0][scalei1][scalenode1-1]);
            vector<double> cscale(4) ; 
            cscale[0] = deltascale+1.0;
            cscale[1] = deltascale;
            cscale[2] = 1.0-deltascale;
            cscale[3] = 2.0-deltascale;

            if(scalenode1<2 || scalenode1>Nscalenode[0]-3){
               cefscale1[0] = 0.0;
               cefscale1[1] = 1.0-cscale[1];
               cefscale1[2] = 1.0-cscale[2];
               cefscale1[3] = 0.0;
            }else{
               cefscale1[1]=1.0-2.5*cscale[1]*cscale[1]+1.5*cscale[1]*cscale[1]*cscale[1];
               cefscale1[2]=1.0-2.5*cscale[2]*cscale[2]+1.5*cscale[2]*cscale[2]*cscale[2];
               cefscale1[0]=2.0 - 4.0*cscale[0] + 2.5*cscale[0]*cscale[0]
                  - 0.5*cscale[0]*cscale[0]*cscale[0];
               cefscale1[3]=2.0 - 4.0*cscale[3] + 2.5*cscale[3]*cscale[3]
                  - 0.5*cscale[3]*cscale[3]*cscale[3];
            }
         }
 
        // define the scale-bin for dimension 2
         int scalenode2 = 0;
         double cefscale2[4] = {0.0,0.5,0.5,0.0};
         if(Nscalenode[1]>1){
            for(int i=1;i<Nscalenode[1];i++){
               if (ScaleFac[1][scalei2]*scale2<ScaleNode[ObsBin][1][scalei2][i]){
                  scalenode2=i;
                  break;
               }
            }
            double deltascale = (ScaleFac[1][scalei2]*scale2 - ScaleNode[ObsBin][1][scalei2][scalenode2-1])/
               (ScaleNode[ObsBin][1][scalei2][scalenode2]-ScaleNode[ObsBin][1][scalei2][scalenode2-1]);
            vector<double> cscale(4) ; 
            cscale[0] = deltascale+1.0;
            cscale[1] = deltascale;
            cscale[2] = 1.0-deltascale;
            cscale[3] = 2.0-deltascale;

            if(scalenode2<2 || scalenode2>Nscalenode[1]-3){
               cefscale2[0] = 0.0;
               cefscale2[1] = 1.0-cscale[1];
               cefscale2[2] = 1.0-cscale[2];
               cefscale2[3] = 0.0;
            }else{
               cefscale2[1]=1.0-2.5*cscale[1]*cscale[1]+1.5*cscale[1]*cscale[1]*cscale[1];
               cefscale2[2]=1.0-2.5*cscale[2]*cscale[2]+1.5*cscale[2]*cscale[2]*cscale[2];
               cefscale2[0]=2.0 - 4.0*cscale[0] + 2.5*cscale[0]*cscale[0]
                  - 0.5*cscale[0]*cscale[0]*cscale[0];
               cefscale2[3]=2.0 - 4.0*cscale[3] + 2.5*cscale[3]*cscale[3]
                  - 0.5*cscale[3]*cscale[3]*cscale[3];
            }
         }

         // ** loop over all scale nodes that receive contributions
         for(int i3 = 0; i3 < 4; i3++){
            int is1 = scalenode1 +i3 -2;  // the target scale index
            if(is1<0){
               is1 = 0;
            }
            if (is1> Nscalenode[0]-1){
               is1 = Nscalenode[0]-1;
            }
            for(int i4 = 0; i4 < 4; i4++){
               int is2 = scalenode2 +i4 -2;  // the target scale index
               if(is2<0){
                  is2 = 0;
               }
               if (is2> Nscalenode[1]-1){
                  is2 = Nscalenode[1]-1;
               }
               int is = is1 * Nscalenode[1] + is2;
               // ** loop over all 4 x points that receive contributions
               for( int i1 = 0; i1 < 4; i1++) {           
                  int x1bin = (nx +i1 -1);
                  if(x1bin<0) x1bin = 0;
                  if(x1bin>Nxtot1[ObsBin]-1) x1bin =Nxtot1[ObsBin]-1;
                  int im = x1bin;   // the target x index
                  for(int proc=0;proc<NSubproc;proc++){
                     //                     printf("%d %d %d %d %d %f %f %f \n",ObsBin,scalevar,is,im,proc,cefm[i1],cefscale1[i3],cefscale2[i4]);
                     SigmaTilde[ObsBin][scalevar][is][im][proc] +=  cefm[i1]  * cefscale1[i3] * cefscale2[i4] *  wt[proc];
                  }
               }
            }
         }
      }
   }
}


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

   // --- Warm-Up Run to identify the extreme x,mu values
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
       fprintf(ofile," // %lu contributions (!= events) in warm-up run \n",IWarmUpCounter);
       fclose(ofile);
     }
     if ( (IWarmUpCounter % IWarmUpPrint) == 0) {
       FILE * ofile;
       ofile = fopen("fastNLO-warmup.dat","w");
       printf(" // %lu contributions (!= events) in warm-up run \n",IWarmUpCounter);
       fprintf(ofile," // %lu contributions (!= events) in warm-up run \n",IWarmUpCounter);
       for (int i=0;i<BlockA2->GetNObsBin();i++){
         printf(" xlim[%d]=%7.5f, mulo[%d]=%8.3f, muup[%d]=%8.3f;\n",i,xlo[i],i,scalelo[i],i,scalehi[i]);
	 //	 fprintf(ofile,"%7.5f, %8.3f, %8.3f;\n",xlo[i],scalelo[i],scalehi[i]);
	 fprintf(ofile,"%e, %e, %e;\n",xlo[i],scalelo[i],scalehi[i]);
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
