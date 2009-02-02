#include "fnloBlockBNlojet.h"

void fnloBlockBNlojet::FillEventDIS(int ObsBin, double x, double scale1, const nlo::amplitude_dis& amp, nlo::pdf_and_coupling_dis& pdf, double prefactor){
   fnloBlockA2 *A2 =  BlockA2;

   if(this->IRef>0){
      amp.pdf_and_qcd_coupling(pdf, prefactor);

      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){
         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         nlo::weight_dis wt = amp(mu2,mu2);

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

      amp.pdf_and_qcd_coupling(pdf, prefactor);

      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){

         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         
         nlo::weight_dis wt = amp(mu2,mu2);

         double pdfsigma = 1.0/3.0*(-wt[1] + 4.*wt[2]);
         double pdfdelta = 3.*(wt[1]-wt[2]);
         wt[1] = pdfdelta;
         wt[2] = pdfsigma;

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

            if(scalenode<2 || scalenode>Nscalenode[0]-3){
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
                  //                  printf("%d %d %d %d %d %f\n",ObsBin,scalevar,is,im,proc,cefscale[i3]);
                  SigmaTilde[ObsBin][scalevar][is][im][proc] +=  cefm[i1]  * cefscale[i3] * wt[proc];
               }
            }
         }
      }
   }
}
void fnloBlockBNlojet::FillEventDIS2Scale(int ObsBin, double x, double scale1, double scale2,  const nlo::amplitude_dis& amp, nlo::pdf_and_coupling_dis& pdf, double prefactor){
   fnloBlockA2 *A2 =  BlockA2;

   if(this->IRef>0){
      amp.pdf_and_qcd_coupling(pdf, prefactor);

      for(int scalevar=0; scalevar<GetTotalScalevars(); scalevar++){
         int scalei1 = scalevar / Nscalevar[1];
         int scalei2 = scalevar % Nscalevar[1];

         double mur2 = ScaleFac[0][scalei1]*ScaleFac[0][scalei1]*scale1*scale1;
         double muf2 = ScaleFac[1][scalei2]*ScaleFac[1][scalei2]*scale2*scale2;
         //         printf("%d %d %d %f %f %f %f\n",scalevar,scalei1,scalei2,ScaleFac[0][scalei1],ScaleFac[1][scalei2],mur2,muf2);
         nlo::weight_dis wt = amp(mur2,muf2);

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

      amp.pdf_and_qcd_coupling(pdf, prefactor);

      for(int scalevar=0; scalevar<GetTotalScalevars(); scalevar++){

         int scalei1 = scalevar / Nscalevar[1];
         int scalei2 = scalevar % Nscalevar[1];

         double mur2 = ScaleFac[0][scalei1]*ScaleFac[0][scalei1]*scale1*scale1;
         double muf2 = ScaleFac[1][scalei2]*ScaleFac[1][scalei2]*scale2*scale2;

         nlo::weight_dis wt = amp(mur2,muf2);

         double pdfsigma = 1.0/3.0*(-wt[1] + 4.*wt[2]);
         double pdfdelta = 3.*(wt[1]-wt[2]);
         wt[1] = pdfdelta;
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
      amp.pdf_and_qcd_coupling(pdf, prefactor);

      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){
         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         nlo::weight_photo wt = amp(mu2,mu2);

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
         printf("fnloBlockB::FillEventPhoto: Error: x (%f) smaller than lowest point (%f) at bin #%d .\n",
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

      amp.pdf_and_qcd_coupling(pdf, prefactor);
      
      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){

         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         nlo::weight_photo wt = amp(mu2,mu2);

         double pdfsigma = 1.0/3.0*(-wt[1] + 4.*wt[2]);
         double pdfdelta = 3.*(wt[1]-wt[2]);
         wt[1] = pdfdelta;
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
      amp.pdf_and_qcd_coupling(pdf, prefactor);

      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){
         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         nlo::weight_hhc wt = amp(mu2,mu2);

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
         printf("fnloBlockB::FillEventResolved: Error: x (%f) smaller than lowest point (%f) at bin #%d .\n",
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
            printf("fnloBlockB::FillEventResolved: Error: x2 (%f) smaller than lowest point (%f) at bin #%d .\n",
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

         amp.pdf_and_qcd_coupling(pdf, prefactor);
      
         for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){

            double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
            nlo::weight_hhc wt = amp(mu2,mu2);

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
                     int im = x1bin + x2bin * Nxtot1[ObsBin];   // the target x index
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
   
   if(this->IRef>0){
      amp.pdf_and_qcd_coupling(pdf, prefactor);

      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){
         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         nlo::weight_hhc wt = amp(mu2,mu2);

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

      if (x1<XNode1[ObsBin][0]){
         printf("fnloBlockB::FillEventHHC: Error: x1 (%f) smaller than lowest point (%f) at bin #%d .\n",
                x1,XNode1[ObsBin][0],ObsBin);
         exit(1);
      }
      //--- determine fractional contribution x1
      double hxlimit = Hxlim1[ObsBin];
      double hx = log10(x1);
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


         if (x2<XNode2[ObsBin][0]){
            printf("fnloBlockB::FillEventHHC: Error: x2 (%f) smaller than lowest point (%f) at bin #%d .\n",
                   x2,XNode2[ObsBin][0],ObsBin);
            exit(1);
         }

         // define the x-bin number in the range  [0:ntot[
         int nx2;
         //-- relative distances in h(x): deltam

         //--- determine fractional contribution x2
         hxlimit = Hxlim2[ObsBin];
         hx = log10(x2);
         hxone = 0.0;

         // define the x-bin number in the range  [0:ntot[
         nx2 = int(Nxtot2[ObsBin] *(hx-hxlimit)/(hxone-hxlimit));

         //-- relative distances in h(x): deltam
         delta  = (hxone-hxlimit)/ Nxtot2[ObsBin];
         hxi = hxlimit+double(nx2)/double(Nxtot2[ObsBin])*(hxone-hxlimit);
         deltam = (hx-hxi)/delta;

         //         printf("x2= %g nx2=%d deltam=%g\n",x2,nx2,deltam);
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

         amp.pdf_and_qcd_coupling(pdf, prefactor);
      
         for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){

            double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
            nlo::weight_hhc wt = amp(mu2,mu2);

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
                     int im = x1bin + x2bin * Nxtot1[ObsBin];   // the target x index
                     //                        if(i3==0) printf("fastNLO: filled at index %d in x2bin #%d x1bin #%d at x=%f\n",im,(nx2 +i2-1),(nx +i1 -1), XNode2[ObsBin][nx2+i2-1]);
                     for(int proc=0;proc<NSubproc;proc++){
                        //                           printf("%d %d %d %d %d %g %g\n",ObsBin,scalevar,is,im,proc,cefscale[i3],photoweight);
                        SigmaTilde[ObsBin][scalevar][is][im][proc] +=  cefm[i1]  * cefm2[i2] * cefscale[i3] * wt[proc];
                     }
                  }
               }
         }
      }
   }
}
