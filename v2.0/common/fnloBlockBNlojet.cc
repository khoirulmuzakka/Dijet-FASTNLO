#include "fnloBlockBNlojet.h"

void fnloBlockBNlojet::FillEventPhoto(int ObsBin, double x, double scale1, const nlo::amplitude_photo& amp, nlo::pdf_and_coupling_photo& pdf){

   fnloBlockA2 *A2 =  BlockA2;
   
   if(this->IRef>0){
      amp.pdf_and_qcd_coupling(pdf, 1.0);

      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){
         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         nlo::weight_photo wt = amp(mu2,mu2);

         double binsize = 1.0;
         for(int dim=0; dim<A2->NDim; dim++){
            binsize *= (A2->UpBin[ObsBin][dim] - A2->LoBin[ObsBin][dim]);
         }
         
         wt *= 389385730./137.0/binsize;
         if(IXsectUnits!=12){
            wt *= pow(10.,(12-IXsectUnits)) ;
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
      if(deltam>1.0 || deltam<0.0 ){
         cout<<" -> delta is off: "<<deltam<<endl;
      }                             

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

      amp.pdf_and_qcd_coupling(pdf, 1.0);
      
      for(int scalevar=0; scalevar<Nscalevar[0]; scalevar++){

         double mu2 = ScaleFac[0][scalevar]*ScaleFac[0][scalevar]*scale1*scale1;
         nlo::weight_photo wt = amp(mu2,mu2);

         double pdfsigma = 1.0/3.0*(-wt[1] + 4.*wt[2]);
         double pdfdelta = 3.*(wt[1]-wt[2]);
         wt[1] = pdfsigma;
         wt[2] = pdfdelta;

         double binsize = 1.0;
         for(int dim=0; dim<A2->NDim; dim++){
            binsize *= (A2->UpBin[ObsBin][dim] - A2->LoBin[ObsBin][dim]);
         }

         wt *= 389385730./137.0/binsize;
         if(IXsectUnits!=12){
            wt *= pow(10.,(12-IXsectUnits)) ;
         }

         // define the scale-bin number in the range  [0:nscalebin[
         int scalenode = 0;
         double cefscale[2] = {0.5,0.5};
         if(Nscalenode[0]>1){
            scalenode = int((Nscalenode[0]-1) *(scale1-ScaleNode[ObsBin][0][scalevar][0])/
                            (ScaleNode[ObsBin][0][scalevar][Nscalenode[0]-1]-ScaleNode[ObsBin][0][scalevar][0]));
            double deltascale = (scale1 - ScaleNode[ObsBin][0][scalevar][scalenode])/
               (ScaleNode[ObsBin][0][scalevar][scalenode+1]-ScaleNode[ObsBin][0][scalevar][scalenode]);
            //               printf("%f %f %f %f\n",deltascale,scale1,ScaleNode[ObsBin][0][scalevar][scalenode+1],ScaleNode[ObsBin][0][scalevar][scalenode]);
            cefscale[0] = 1.-deltascale;
            cefscale[1] = deltascale;
         }

         //         printf("pt=%f nnodes=%d lowest=%f highest=%f bin=%d\n",scale1,Nscalenode[0],ScaleNode[ObsBin][0][scalevar][0],ScaleNode[ObsBin][0][scalevar][Nscalenode[0]-1],scalenode);

         // ** loop over all 2 scale nodes that receive contributions
         for(int i3 = 0; i3 < 2; i3++){
            int is = scalenode +i3 -1;  // the target scale index
            if(is<0){
               is = 0;
            }
            if (is> Nscalenode[0]-1){
               is = Nscalenode[0]-1;
            }
            // ** loop over all 4 x points that receive contributions
            for( int i1 = 0; i1 < 4; i1++) {           
               int im = nx +i1 -1;   // the target x index
               if (im>=0 && im<Nxtot1[0]) {
                  for(int proc=0;proc<NSubproc;proc++){
                     //                     printf("%d %d %d %d %d %f\n",ObsBin,scalevar,is,im,proc,cefscale[i3]);
                     SigmaTilde[ObsBin][scalevar][is][im][proc] +=  cefm[i1]  * cefscale[i3] * wt[proc];
                  }
               }
            }
         }
      }
   }
}
