// Author: Daniel Britzger
// DESY, 02/04/2012

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  FastNLODiffReader                                                   //
//  D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch           //
//                                                                      //
//  The projects web page can be found at:                              // 
//    http://projects.hepforge.org/fastnlo                              //
//                                                                      //
//  If you use this code, please cite:                                  //
//    T. Kluge, K. Rabbertz and M. Wobisch, hep-ph/0609285              //
//    D. Britzger, T. Kluge, K. Rabbertz, F. Stober, M. Wobisch,        //
//       arXiv:1109.1310                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "FastNLODiffReader.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cfloat>

using namespace std;

//______________________________________________________________________________


FastNLODiffReader::FastNLODiffReader(string filename) : FastNLOReader(filename)
{
   cout << "yeah! it's diffraction!"<<endl;
   fzmin = 0; fzmax = 1.;
   fxpom = 0.01;
}


//______________________________________________________________________________



void FastNLODiffReader::SetXPomSlicing(int nSlice, double* xpom, double* dxpom)
{
   fxPoms.resize(nSlice);
   fdxPoms.resize(nSlice);
   for ( int i = 0 ; i<nSlice ; i++ ){
      fxPoms[i]  = xpom[i];
      fdxPoms[i] = dxpom[i];
   }
}


//______________________________________________________________________________


void FastNLODiffReader::SetXPomLogSlicing(int nStep, double xpommin, double xpommax )
{

   if ( xpommin < 1.e-4 ) {
      printf("FastNLODiffReader::SetXPomLogSlicing. xpommin should not be too small in order to have sufficent nodes.\n");
      if ( xpommin == 0 ) {
	 printf("FastNLODiffReader::SetXPomLogSlicing. xpommin should espc. not be '0'!\n");
	 exit(1);
      }
   }

   // new array
   double *binning = new double[nStep+1];
   double *dxpom   = new double[nStep+1];
   double *xpom    = new double[nStep+1];
   
   double delta_x_log = (log10(xpommax)-log10(xpommin))/nStep;
   
   binning[0]=xpommin;
   // put equidistant binwith on a logarithmic scale
   for(int i=1;i<=nStep;++i){
      binning[i] = pow(10.,(log10(binning[i-1])+ delta_x_log));
      dxpom[i-1] = binning[i] - binning[i-1];
      xpom[i-1]  = pow(10.,((log10(binning[i-1])+log10(binning[i]))/2.));
      //cout << "binning[i] = "<<binning[i]<<"\tdxpom = "<<dxpom[i-1] << "\txpom = " << xpom[i-1] << endl;
   }

   
   SetXPomSlicing( nStep, xpom, dxpom);
}

//______________________________________________________________________________


void FastNLODiffReader::SetXPomExpSlicing(int nStep, double xpommin, double xpommax )
{

   // new array
   double *binning = new double[nStep+1];
   double *dxpom   = new double[nStep+1];
   double *xpom    = new double[nStep+1];
   
   double delta_x_log = (exp(xpommax)-exp(xpommin))/nStep;
   
   binning[0]=xpommin;
   // put equidistant binwith on a logarithmic scale
   for(int i=1;i<=nStep;++i){
      binning[i] = log(exp(binning[i-1])+ delta_x_log);
      dxpom[i-1] = binning[i] - binning[i-1];
      xpom[i-1]  = log((exp(binning[i-1])+exp(binning[i]))/2.);
      cout << "binning[i] = "<<binning[i]<<"\tdxpom = "<<dxpom[i-1] << "\txpom = " << xpom[i-1] << endl;
   }
   
   SetXPomSlicing( nStep, xpom, dxpom);
}

//______________________________________________________________________________


void FastNLODiffReader::SetXPomLinSlicing(int nStep, double xpommin, double xpommax )
{
   // new array
   double *binning = new double[nStep+1];
   double *dxpom   = new double[nStep+1];
   double *xpom    = new double[nStep+1];
   double delta_x_log = (xpommax-xpommin)/nStep;
   binning[0]=xpommin;
   // put equidistant binwith on a linear scale
   for(int i=1;i<=nStep;++i){
      binning[i] = binning[i-1]+ delta_x_log;
      dxpom[i-1] = binning[i] - binning[i-1];
      xpom[i-1]  = (binning[i-1]+binning[i])/2.;
   }
   SetXPomSlicing( nStep, xpom, dxpom);
}

//______________________________________________________________________________
void FastNLODiffReader::FillPDFCache( bool ReCalcCrossSection ){
   printf("FastNLODiffReader::FillPDFCache(). ERROR. PDF Cache cannot be filled in diffractive version, since xpom integration has still to be performed\n");
   printf("  Please access directly FastNLODiffReader::GetDiffCrossSection()\n");
   exit(1);
}

//______________________________________________________________________________
void FastNLODiffReader::CalcCrossSection(){
   printf("FastNLODiffReader::CalcCrossSection(). ERROR. This method is not valid for diffractive tables.\n");
   printf("  Please access directly FastNLODiffReader::GetDiffCrossSection()\n");
   exit(1);
}

//______________________________________________________________________________

vector<double> FastNLODiffReader::GetReferenceCrossSection( ){
   printf("FastNLODiffReader::GetReferenceCrossSection(). No reference cross sections in diffractive version\n");
   return vector<double>();
}

//______________________________________________________________________________

vector < double > FastNLODiffReader::GetCrossSection( ){
   return GetDiffCrossSection();
}

vector < double > FastNLODiffReader::GetDiffCrossSection( ){
   // Get fast calculated NLO cross section

   vector < double > xs(NObsBin);
   
   double interv = 0;
   // do the xpom integration
   //fprintf(stderr," * Integration xpom in %d slices. [",fxPoms.size());
   printf(" * Integrating xpom in %d slices. [",fxPoms.size());
   fflush(stdout);
   FillAlphasCache();
   for ( unsigned int ixp = 0 ; ixp<fxPoms.size() ; ixp++ ){
      fxpom = fxPoms[ixp];
      // always recalculate cross section
      FastNLOReader::FillPDFCache(false);
      FastNLOReader::CalcCrossSection();

      double slicesize = fdxPoms[ixp];
      interv+=slicesize;
      
      printf(".");
      fflush(stdout);
      
      for ( unsigned int i = 0 ; i<NObsBin ; i++ ){
	 xs[i] += XSection[i] * slicesize ;
      }
   }
   printf("]\n");
   cout << " * Integrated interval in xpom: " << interv << endl;

   // set this cross section also to FastNLO mother class
   XSection = xs;

   // k-factors
   for ( int i = 0 ; i<NObsBin ; i++ ){
      FastNLOReader::kFactor[i] = FastNLOReader::XSection[i] / FastNLOReader::XSection_LO[i];
   }


   return xs;
}


//______________________________________________________________________________
