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

#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include "fastnlo/FastNLODiffReader.h"

using namespace std;

//______________________________________________________________________________


FastNLODiffReader::FastNLODiffReader(string filename) : FastNLOReader(filename)
{
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
      debug["SetXPomlicing"]<<"[i]="<<i<<"\tfxPomx[i]="<<fxPoms[i]<<"\tfdxPoms[i]="<<fdxPoms[i]<<endl;
   }
}


//______________________________________________________________________________


void FastNLODiffReader::SetXPomLogSlicing(int nStep, double xpommin, double xpommax )
{

   if ( xpommin < 1.e-4 ) {
      warn["SetXPomLogSlicing"]<<"xpommin should not be too small in order to have sufficent nodes."<<endl;
      if ( xpommin == 0 ) {
	 warn["SetXPomLogSlicing"]<<"xpommin should espc. not be '0'!"<<endl;
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
      debug["SetXPomExpSlicing"]<< "binning[i] = "<<binning[i]<<"\tdxpom = "<<dxpom[i-1] << "\txpom = " << xpom[i-1] << endl;
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
   error["FillPDFCache"]<<"PDF Cache cannot be filled in diffractive version, since xpom integration has still to be performed"<<endl;
   error>>"  Please access directly FastNLODiffReader::GetDiffCrossSection()"<<endl;
   exit(1);
}

//______________________________________________________________________________
void FastNLODiffReader::CalcCrossSection(){
   error["CalcCrossSection"]<<"This method is not valid for diffractive tables."<<endl;
   error>>"  Please access directly FastNLODiffReader::GetDiffCrossSection()\n";
   exit(1);
}

//______________________________________________________________________________

vector<double> FastNLODiffReader::GetReferenceCrossSection( ){
   error["GetReferenceCrossSection"]<<"No reference cross sections in diffractive version"<<endl;
   return vector<double>();
}

//______________________________________________________________________________

void FastNLODiffReader::PrintCrossSectionsWithReference(){
   error["PrintCrossSectionsWithReference"]<<"No reference cross sections in diffractive version"<<endl;
   return;
}

//______________________________________________________________________________

vector < double > FastNLODiffReader::GetCrossSection( ){
   return GetDiffCrossSection();
}

vector < double > FastNLODiffReader::GetDiffCrossSection( ){
   // Get fast calculated NLO cross section

   vector < double > xs(NObsBin);
   vector < double > xsLO(NObsBin);
   if ( fxPoms.empty() ){
      error["GetDiffCrossSection"]<<"No xpom slicing given."<<endl;
      return xs;
   }

   double interv = 0;
   // do the xpom integration
   info["GetDiffCrossSection"]<<"Integrating xpom in "<<fxPoms.size()<<" slices. [";
   fflush(stdout);
   for ( unsigned int ixp = 0 ; ixp<fxPoms.size() ; ixp++ ){
      fxpom = fxPoms[ixp];
      // always recalculate cross section
      FastNLOReader::CalcCrossSection(); // this calls GetXFX() very very often!

      for ( int i = 0 ; i<NObsBin ; i++ ){
	 if (i==0)debug["GetDiffCrossSection"]<<"i="<<i<<"\tixp="<<ixp<<"\tfxpom="<<fxpom<<"\tXSection[i]="<<XSection[i]<<"\tfdxPoms[ixp]="<<fdxPoms[ixp]<<endl;;
	 xs[i] += XSection[i] * fdxPoms[ixp] ;
	 xsLO[i] += XSection_LO[i] * fdxPoms[ixp] ;
      }
      info>>"."; fflush(stdout);
      interv+=fdxPoms[ixp];
   }
   info>>"]"<<endl;
   info["GetDiffCrossSection"]<< "Integrated interval in xpom: " << interv << endl;

   // set this cross section also to FastNLO mother class
   XSection = xs;
   XSection_LO = xsLO;

   // k-factors
   for ( int i = 0 ; i<NObsBin ; i++ ){
      FastNLOReader::kFactor[i] = FastNLOReader::XSection[i] / FastNLOReader::XSection_LO[i];
   }

   return xs;
}


//______________________________________________________________________________

vector<double> FastNLODiffReader::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  x-proton and factorisation scale.
   //

   // get pdf
   double zpom = xp/fxpom;
   vector < double > a(13);
   if ( zpom > fzmin && zpom < fzmax ){
      // find x-node index
      int nx = -1;
      int nb = -1;
      for ( unsigned int ib = 0 ; nb == -1 && ib<BBlocksSMCalc[0][0]->XNode1.size() ; ib++ ) {
	 if ( BBlocksSMCalc[0][0]->GetNxmax(ib) != BBlocksSMCalc[0][1]->GetNxmax(ib) )
	    error["FastNLODiffReader::GetXFX"]<<"LO and NLO tables must have same number of x-bins."<<endl;
	 for ( unsigned int ix = 0 ; nx == -1 && ix<BBlocksSMCalc[0][0]->XNode1[ib].size(); ix++ ) {
	    if ( BBlocksSMCalc[0][0]->XNode1[ib][ix] != BBlocksSMCalc[0][1]->XNode1[ib][ix] )
	       error["FastNLODiffReader::GetXFX"]<<"LO and NLO tables must have idnetical x-bins."<<endl;
	    if ( xp == BBlocksSMCalc[0][0]->XNode1[ib][ix] ) {
	       nx = ix;
	       nb = ib;
	    }
	 }
      }

      // check if this is the 'last' or 'first' xnode
      bool IsLastX  = nx == (int)BBlocksSMCalc[0][0]->XNode1[nb].size()-1 ;
      bool IsFirstX = nx == 0 ;

      if ( nx == -1 || nb == -1 ) {
	 //printf("Warning. Could not find x-node index for xp = %12.8e.\n",xp);
	 for ( int ib = 0 ; nb == -1 && ib<(int)BBlocksSMCalc[0][0]->XNode1.size() ; ib++ ) {
	    for ( int ix = 0 ; nx == -1 && ix<(int)BBlocksSMCalc[0][0]->XNode1[ib].size(); ix++ ) {
	       if ( xp == BBlocksSMCalc[0][0]->XNode1[ib][ix] ) {
		  nx = ix;
		  nb = ib;
	       }
	       if (  fabs ( 1. - xp / BBlocksSMCalc[0][0]->XNode1[ib][ix] ) < 1.e-6 ){
		  warn["FastNLODiffReader::GetXFX"]<<"Could not find x-node index for xp = "<<xp<<endl;
		  warn>>"   but a quite close one: xp = "<<xp<<", xnode = "<<BBlocksSMCalc[0][0]->XNode1[ib][ix]<<endl;
	       }
	    }
	 }
	 //exit(1);
	 IsLastX = true;
	 IsFirstX = true;
      }

      a = GetDiffXFX( fxpom, zpom, muf);

      // calc reweight at integration edges
      if ( !IsLastX && !IsFirstX ){
	 const double x2 = BBlocksSMCalc[0][0]->XNode1[nb][nx+1];// next node
	 const double x1 = BBlocksSMCalc[0][0]->XNode1[nb][nx-1];// prev. node
	 const double zpom2 = x2/fxpom;
	 const double zpom1 = x1/fxpom;
	 double xSpan = 1.;
	 // wenn jetzt der naechste bin nicht mehr in fzmax ist, dann wird gewichtet
	 if ( zpom2 > fzmax && zpom < fzmax ){
	    double xmax = fzmax*fxpom;
	    double ldelx = log10(xmax) - log10(xp);
	    double ldelx0 = log10(x2) - log10(xp);
	    xSpan *= ldelx/ldelx0 + 0.5 ;
	 }
	 if ( zpom1 < fzmin && zpom > fzmin ){
	    double xmin = fzmin*fxpom;
	    double ldelx = log10(xp) - log10(xmin);
	    double ldelx0 = log10(xp) - log10(x1);
	    xSpan *= ldelx/ldelx0 + 0.5 ;
	 }
	 if ( xSpan != 1. ){
	    for ( unsigned  int i = 0 ; i<a.size() ; i++ ) a[i]*=xSpan;
	 }
      }
   }

   return a;
}


//______________________________________________________________________________
