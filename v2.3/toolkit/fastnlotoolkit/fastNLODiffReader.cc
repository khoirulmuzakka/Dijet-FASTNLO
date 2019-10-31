// Author: Daniel Britzger
// DESY, 02/04/2012
// Update for toolkit: DB: 12/05/14

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLODiffReader                                                   //
//                                                                      //
//  fastNLODiffReader for the evaluation of fastNLO tables in           //
//  diffractive DIS.                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include "fastnlotk/fastNLODiffReader.h"

using namespace std;

//______________________________________________________________________________


fastNLODiffReader::fastNLODiffReader(string filename) : fastNLOReader(filename) {
   fzmin = 0;
   fzmax = 1.;
   fxpom = 0.01;
   fProtonE = 920.;
}


//______________________________________________________________________________



void fastNLODiffReader::SetXPomSlicing(int nSlice, double* xpom, double* dxpom) {
   fxPoms.resize(nSlice);
   fdxPoms.resize(nSlice);
   for (int i = 0 ; i<nSlice ; i++) {
      fxPoms[i]  = xpom[i];
      fdxPoms[i] = dxpom[i];
      logger.debug["SetXPomlicing"]<<"[i]="<<i<<"\tfxPomx[i]="<<fxPoms[i]<<"\tfdxPoms[i]="<<fdxPoms[i]<<endl;
   }
   fPDFCached = 0; // reset cache.
}


//______________________________________________________________________________


void fastNLODiffReader::SetXPomLogSlicing(int nStep, double xpommin, double xpommax) {

   if (xpommin < 1.e-4) {
      logger.warn["SetXPomLogSlicing"]<<"xpommin should not be too small in order to have sufficent nodes."<<endl;
      if (xpommin == 0) {
         logger.warn["SetXPomLogSlicing"]<<"xpommin should espc. not be '0'!"<<endl;
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
   for (int i=1; i<=nStep; ++i) {
      binning[i] = pow(10.,(log10(binning[i-1])+ delta_x_log));
      dxpom[i-1] = binning[i] - binning[i-1];
      xpom[i-1]  = pow(10.,((log10(binning[i-1])+log10(binning[i]))/2.));
   }
   SetXPomSlicing(nStep, xpom, dxpom);
}

//______________________________________________________________________________


void fastNLODiffReader::SetXPomExpSlicing(int nStep, double xpommin, double xpommax) {

   // new array
   double *binning = new double[nStep+1];
   double *dxpom   = new double[nStep+1];
   double *xpom    = new double[nStep+1];

   double delta_x_log = (exp(xpommax)-exp(xpommin))/nStep;

   binning[0]=xpommin;
   // put equidistant binwith on a logarithmic scale
   for (int i=1; i<=nStep; ++i) {
      binning[i] = log(exp(binning[i-1])+ delta_x_log);
      dxpom[i-1] = binning[i] - binning[i-1];
      xpom[i-1]  = log((exp(binning[i-1])+exp(binning[i]))/2.);
      logger.debug["SetXPomExpSlicing"]<< "binning[i] = "<<binning[i]<<"\tdxpom = "<<dxpom[i-1] << "\txpom = " << xpom[i-1] << endl;
   }

   SetXPomSlicing(nStep, xpom, dxpom);
}

//______________________________________________________________________________


void fastNLODiffReader::SetXPomLinSlicing(int nStep, double xpommin, double xpommax) {
   // new array
   double *binning = new double[nStep+1];
   double *dxpom   = new double[nStep+1];
   double *xpom    = new double[nStep+1];
   double delta_x_log = (xpommax-xpommin)/nStep;
   binning[0]=xpommin;
   // put equidistant binwith on a linear scale
   for (int i=1; i<=nStep; ++i) {
      binning[i] = binning[i-1]+ delta_x_log;
      dxpom[i-1] = binning[i] - binning[i-1];
      xpom[i-1]  = (binning[i-1]+binning[i])/2.;
   }
   SetXPomSlicing(nStep, xpom, dxpom);
}

//______________________________________________________________________________
void fastNLODiffReader::FillPDFCache(bool ReCalcCrossSection) {
   logger.error["FillPDFCache"]<<"PDF Cache cannot be filled in diffractive version, since xpom integration has still to be performed"<<endl;
   logger.error>>"  Please access directly fastNLODiffReader::GetDiffCrossSection()"<<endl;
   exit(1);
}

//______________________________________________________________________________
void fastNLODiffReader::CalcCrossSection() {
   logger.error["CalcCrossSection"]<<"This method is not valid for diffractive tables."<<endl;
   logger.error>>"  Please access directly fastNLODiffReader::GetDiffCrossSection()\n";
   exit(1);
}

//______________________________________________________________________________

vector<double> fastNLODiffReader::GetReferenceCrossSection() {
   logger.error["GetReferenceCrossSection"]<<"No reference cross sections in diffractive version"<<endl;
   return vector<double>();
}

//______________________________________________________________________________

void fastNLODiffReader::PrintCrossSectionsWithReference() {
   logger.error["PrintCrossSectionsWithReference"]<<"No reference cross sections in diffractive version"<<endl;
   return;
}

//______________________________________________________________________________

vector < double > fastNLODiffReader::GetCrossSection() {
   return GetDiffCrossSection();
}

vector < double > fastNLODiffReader::GetDiffCrossSection() {
   // Get fast calculated NLO cross section

   // check cached values:
   const double asNew = CalcReferenceAlphas();
   bool IsAsCached = fAlphasCached == CalcReferenceAlphas();
   bool IsPDFCached = fPDFCached == CalcNewPDFChecksum();
   if ( IsAsCached && IsPDFCached ) {
      logger.debug["GetDiffCrossSection"]<<"No need for re-calculation of cross section. All values cached. asNew = "<<asNew<<endl;
      return XSection;
   }

   vector < double > xs(NObsBin);
   vector < double > xsLO(NObsBin);
   XSection_LO.resize(NObsBin);
   if (fxPoms.empty()) {
      logger.error["GetDiffCrossSection"]<<"No xpom slicing given."<<endl;
      return xs;
   }

   double interv = 0;
   // do the xpom integration
   logger.info["GetDiffCrossSection"]<<"Integrating xpom in "<<fxPoms.size()<<" slices. [";
   fflush(stdout);
   
   fXSection_vs_xIPzIP.clear();
   fXSection_vs_xIPzIP.resize(NObsBin);
   if ( fPrintxIPzIP ) {
      logger.info>>"]"<<endl;
      printf("%8s%14s%14s%14s\n","ObsBin","xIP","zIP","cs");
   }
   for (unsigned int ixp = 0 ; ixp<fxPoms.size() ; ixp++) {
      fxpom = fxPoms[ixp];
      // always recalculate cross section
      fastNLOReader::CalcCrossSection(); // this calls GetXFX() very very often!

      for (unsigned int i = 0 ; i<NObsBin ; i++) {
         if (i==0) logger.debug["GetDiffCrossSection"]<<"i="<<i<<"\tixp="<<ixp<<"\tfxpom="<<fxpom<<"\tXSection[i]="<<XSection[i]<<"\tfdxPoms[ixp]="<<fdxPoms[ixp]<<endl;
         xs[i] += XSection[i] * fdxPoms[ixp] ;
         xsLO[i] += XSection_LO[i] * fdxPoms[ixp] ;
	 
	 //vector<map<double,double> > fXSection_vsX1;
	 //std::vector < std::map< std::pair<double, double>, double > > fXSection_vs_xIPzIP;
	 if ( fPrintxIPzIP ) {
	    for ( auto xc : fXSection_vsX1[i] ) {
	       printf("%8d%14.6f%14.6f%14.6f\n",i,fxpom,xc.first,xc.second*fdxPoms[ixp]);
	       fXSection_vs_xIPzIP[i][make_pair(fxpom,xc.first)] = xc.second*fdxPoms[ixp];
	       //cout<<i<<"\t"<<fxpom<<"\t"<<xc.first<<"\t"<<xc.second*fdxPoms[ixp]<<endl;;
	    }
	 }
	 // for ( auto xc : fXSection_vsX1[i] ) {
	 //    fXSection_vs_xIPzIP[i][make_pair(xc.first,fxpom)] = xc.second*fdxPoms[ixp];
	 // }
	 // for ( auto d : fXSection_vsX1[0] ) cout<<" | "<<d.first<<", "<<d.second;
	 // cout<<endl;
      }

      // radek
      for(unsigned i = 0; i < NObsBin; ++i) {
         //cout << "The basic cross section in bin "<<i <<" : " << XSection[i] << endl;                                                                                            
         double sumMy=0;
         for(auto &x :  fXSection_vsQ2[i]) {
            sumMy += x.second;
            x.second *= fdxPoms[ixp];
         }
      }
      xsQ2[fxpom] = fXSection_vsQ2;

      if ( !fPrintxIPzIP ) logger.info>>".";
      fflush(stdout);
      interv+=fdxPoms[ixp];
   }
   logger.info>>"]"<<endl;
   logger.info["GetDiffCrossSection"]<< "Integrated interval in xpom: " << interv << endl;
   
   // set this cross section also to FastNLO mother class
   XSection = xs;
   XSection_LO = xsLO;

   // k-factors
   fastNLOReader::kFactor.resize(NObsBin);
   for (unsigned int i = 0 ; i<NObsBin ; i++) {
      fastNLOReader::kFactor[i] = fastNLOReader::XSection[i] / fastNLOReader::XSection_LO[i];
   }

   return xs;
}


//______________________________________________________________________________

vector<double> fastNLODiffReader::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pdf-interface. It should return a vector of 13
   //  parton flavors from tbar to t at a certain
   //  x-proton and factorisation scale.
   //

   // get pdf
   double zpom = xp/(fxpom * fProtonE/920. );
   //double zpom = xp/fxpom;
   vector < double > a(13);
   if (zpom > fzmin && zpom < fzmax) {
      // find x-node index
      int nx = -1;
      int nb = -1;
      for (int ib = 0 ; nb == -1 && ib<B_LO()->GetNObsBin() ; ib++) {
	 if (B_NLO() && (B_LO()->GetNxmax(ib) != B_NLO()->GetNxmax(ib)))
             logger.error["fastNLODiffReader::GetXFX"]<<"LO and NLO tables must have same number of x-bins."<<endl;
          for (int ix = 0 ; nx == -1 && ix<B_LO()->GetNxtot1(ib); ix++) {
             if (B_NLO() && (B_LO()->GetXNode1(ib,ix) != B_NLO()->GetXNode1(ib,ix)))
                logger.error["fastNLODiffReader::GetXFX"]<<"LO and NLO tables must have idnetical x-bins."<<endl;
             if ( xp == B_LO()->GetXNode1(ib,ix) ) {
                nx = ix;
                nb = ib;
             }
          }
      }

      // check if this is the 'last' or 'first' xnode
      bool IsLastX  = nx == (int)B_LO()->GetNxtot1(nx)-1 ;
      bool IsFirstX = nx == 0 ;

      if (nx == -1 || nb == -1) {
         //printf("Warning. Could not find x-node index for xp = %12.8e.\n",xp);
         for (int ib = 0 ; nb == -1 && ib<(int)B_LO()->GetNObsBin() ; ib++) {
            for (int ix = 0 ; nx == -1 && ix<(int)B_LO()->GetNxtot1(ib); ix++) {
               if (xp == B_LO()->GetXNode1(ib,ix)) {
                  nx = ix;
                  nb = ib;
               }
               if (fabs(1. - xp / B_LO()->GetXNode1(ib,ix)) < 1.e-6) {
                  logger.warn["fastNLODiffReader::GetXFX"]<<"Could not find x-node index for xp = "<<xp<<endl;
                  logger.warn>>"   but a quite close one: xp = "<<xp<<", xnode = "<<B_LO()->GetXNode1(ib,ix)<<endl;
               }
            }
         }
         //exit(1);
         IsLastX = true;
         IsFirstX = true;
      }

      a = GetDiffXFX(fxpom, zpom, muf);

      // calc reweight at integration edges
      if (!IsLastX && !IsFirstX) {
         const double x2 = B_LO()->GetXNode1(nb,nx+1) ;// next node
         const double x1 = B_LO()->GetXNode1(nb,nx-1) ;// prev. node
         const double zpom2 = x2/fxpom;
         const double zpom1 = x1/fxpom;
         double xSpan = 1.;
         // wenn jetzt der naechste bin nicht mehr in fzmax ist, dann wird gewichtet
         if (zpom2 > fzmax && zpom < fzmax) {
            double xmax = fzmax*fxpom;
            double ldelx = log10(xmax) - log10(xp);
            double ldelx0 = log10(x2) - log10(xp);
            xSpan *= ldelx/ldelx0 + 0.5 ;
         }
         if (zpom1 < fzmin && zpom > fzmin) {
            double xmin = fzmin*fxpom;
            double ldelx = log10(xp) - log10(xmin);
            double ldelx0 = log10(xp) - log10(x1);
            xSpan *= ldelx/ldelx0 + 0.5 ;
         }
         if (xSpan != 1.) {
            for (unsigned  int i = 0 ; i<a.size() ; i++) a[i]*=xSpan;
         }
      }
   }

   return a;
}


//______________________________________________________________________________
