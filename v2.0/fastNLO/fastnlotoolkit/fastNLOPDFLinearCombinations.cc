#include <cmath>
#include <cstdlib>
#include "fastnlotk/fastNLOConstants.h"
#include "fastnlotk/fastNLOPDFLinearCombinations.h"

using namespace std;
using namespace fastNLO;


fastNLOPDFLinearCombinations::fastNLOPDFLinearCombinations() {
   //: PrimalScream("fastNLOPDFLinearCombinations")  {
}


fastNLOPDFLinearCombinations::~fastNLOPDFLinearCombinations(){
}


//______________________________________________________________________________


vector<double > fastNLOPDFLinearCombinations::CalcPDFLinearCombination(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1, const vector<double>& pdfx2 , bool pdf2IsAntiParticle ) const {
   //
   // return vector of PDF linear combinations which can be used
   // to calculate the cross section.
   //
   // bool pdf2IsAntiParticle specifies, if pdfx2 still has to be 'inverted' or not
   //

   switch ( c->GetNPDF() ) {
   case 0:	// no PDF involved in process; e.g. e+e-
      return vector<double >();
      break;
   case 1:	// one PDF invovled in process: e.g. DIS
      return CalcPDFLCOneHadron(c,pdfx1);
      break;
   case 2:	// two PDFs involved in process: e.g. pp, ppbar
      if ( !pdf2IsAntiParticle ) {
	 vector<double> Antipdf2 = MakeAntiHadron(pdfx2);
	 return CalcPDFLCTwoHadrons(c,pdfx1,Antipdf2);
      }
      else return CalcPDFLCTwoHadrons(c,pdfx1,pdfx2);
      break;
   default:
      //error["CalcPDFLinearCombination"]<<"Unknown number of PDFs involved in process. NPDF="<<c->GetNPDF()<<endl;
      say::error<<"[CalcPDFLinearCombination] Unknown number of PDFs involved in process. NPDF="<<c->GetNPDF()<<endl;
      exit(1);
      return vector<double >();
   }
}


//______________________________________________________________________________


vector<double > fastNLOPDFLinearCombinations::CalcPDFLCOneHadron(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 ) const {
   //
   // Calculate PDF linear combinations for DIS processes
   //

   // ---- check for ep-DIS ---- //
   bool IsONEPDF = ( c->GetNPDF() == 1 );
   bool IsDIS    = ( c->GetIPDFdef1() == 2 );
   bool IsNCDIS  = ( c->GetIPDFdef2() == 1 );
   bool IsProton = ( c->GetPDFPDG(0) == 2212 );
   if ( IsDIS && IsONEPDF && IsNCDIS && IsProton ) return CalcPDFLinearCombDIS(c,pdfx1);
   // ---- unknown process ---- //
   else {
      //error["CalcPDFLCDIS"]<<"Could not identify process. Printing and exiting"<<endl;
      say::error<<"Error. Could not identify process. Printing and exiting"<<endl;
      c->Print();
      exit(1);
      return vector<double >();
   }
}



//______________________________________________________________________________



vector<double > fastNLOPDFLinearCombinations::CalcPDFLCTwoHadrons(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1, const vector<double>& pdfx2 ) const {
   //
   // Calculate PDF linear combinations for processes with two hadrons
   //
   // -----------------------------------------------
   // implement other processes e.g. like this:
   //    bool IsDrellYan = (...)   // check for the identifier(s) (most likely IPDFDef2
   //    if ( IsDrellYan ) CalcPDFLinearCombDrellYan(c,pdfx1,pdfx2);
   // and implement a new function called CalcPDFLinearCombDrellYan(...)
   // -----------------------------------------------

   // ---- check for jet-productions  ---- //
   bool IsTwoPDF = ( c->GetNPDF() == 2 );
   bool IsTwoIdenticHadrons = (c->GetIPDFdef1() == 3  &&  c->GetPDFPDG(0) == fabs(c->GetPDFPDG(1)) );
   bool IsHHJets = ( c->GetIPDFdef2() == 1 );

   if ( IsTwoPDF && IsTwoIdenticHadrons && IsHHJets ) {
      return CalcPDFLinearCombHHC(c,pdfx1,pdfx2);
   }
   // else if (...)  //space for other processes
   // ---- (yet) unknown process ---- //
   else {
      say::error<<"[CalcPDFLinearCombination] Could not identify process. Printing and exiting..."<<endl;
      say::error<<"IsTwoPDF = "<<IsTwoPDF<<endl;
      say::error<<"IsTwoIdenticHadrons = "<<IsTwoIdenticHadrons<<endl;
      say::error<<"IsHHJets = "<<IsHHJets<<endl;
      c->Print();
      exit(1);
      return vector<double >();
   }
}


//______________________________________________________________________________


vector<double > fastNLOPDFLinearCombinations::MakeAntiHadron(const vector<double>& xfx ) const {
   vector < double > xfxbar(13);
   for (unsigned int p = 0 ; p<13 ; p++) {
      xfxbar[p] = xfx[12-p];
   }
   return xfxbar;
}

//______________________________________________________________________________


vector<double> fastNLOPDFLinearCombinations::CalcPDFLinearCombDIS(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1) const {
   //
   //  Internal method.
   //  CalcPDFLinearComb is used to calculate the
   //  linear combinations of different parton flavors
   //  according to the used subprocesses of your
   //  calculation/table.
   //

   int NSubproc = c->GetNSubproc();

   vector < double > pdflc(3);
   pdflc[1] = pdfx1[6]; //gluon
   for (int l=0; l<13; l++) {
      double temp = (l==6 ? 0.0 : pdfx1[l]);
      if (!(l&1)) temp *= 4.;
      pdflc[0] += temp; // delta
   }
   pdflc[0] /= 9.;
   if (NSubproc>2) { // only from NLO
      for (int l=0; l<6; l++) {
         pdflc[2] += pdfx1[5-l] + pdfx1[l+7]; // sigma
      }
   }
   return pdflc;
}

//______________________________________________________________________________


vector<double> fastNLOPDFLinearCombinations::CalcPDFLinearCombHHC(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const {
   //
   //  Internal method.
   //  CalcPDFLinearComb is used to calculate the
   //  linear combinations of different parton flavors
   //  according to the used subprocesses of your
   //  calculation/table.
   //

   int NSubproc = c->GetNSubproc();

   double SumQ1  = 0;
   double SumQB1 = 0;
   double SumQ2  = 0;
   double SumQB2 = 0;
   vector <double> Q1(6);
   vector <double> QB1(6);
   vector <double> Q2(6);
   vector <double> QB2(6);
   for (int k = 0 ; k<6 ; k++) {
      Q1[k]  = pdfx1[k+7];  //! read 1st PDF at x1
      QB1[k] = pdfx1[5-k];
      SumQ1  += Q1[k];
      SumQB1 += QB1[k];
      Q2[k]  = pdfx2[k+7];//  ! read 2nd PDF at x2
      QB2[k] = pdfx2[5-k];
      SumQ2  += Q2[k];
      SumQB2 += QB2[k];
   }
   double G1     = pdfx1[6];
   double G2     = pdfx2[6];

   //   - compute S,A
   double S = 0;
   double A = 0;
   for (int k = 0 ; k<6 ; k++) {
      S += (Q1[k]*Q2[k]) + (QB1[k]*QB2[k]);
      A += (Q1[k]*QB2[k]) + (QB1[k]*Q2[k]);
   }

   //c   - compute seven combinations
   vector <double> H(7);
   H[0]  = G1*G2;
   H[1] = SumQ1*SumQ2 + SumQB1*SumQB2 - S;
   H[2] = S;
   H[3] = A;
   H[4] = SumQ1*SumQB2 + SumQB1*SumQ2 - A;
   H[5] = (SumQ1+SumQB1)*G2;
   H[6] = G1*(SumQ2+SumQB2);

   if (NSubproc == 6) {
      H[5] += H[6];
      H.resize(6);
   }
   return H;

}


//______________________________________________________________________________


vector<double> fastNLOPDFLinearCombinations::CalcPDFLinearCombttbar(const fastNLOCoeffAddBase* c, const vector<double>& pdfx1 , const vector<double>& pdfx2) const {
   //
   // Calculate pdf-lcs for ttbar cross sections in pp/ppbar
   // Used by NNLO generator of Marco Guzzi
   //
   // pdf[0] = gg
   // pdf[1] = qq
   //
   // DB 24.08.13
   //


   vector <double> pdflc(2);
   pdflc[0] += pdfx1[6]*pdfx2[6]; // gg
   // qq
   for (int k = 0 ; k<6 ; k++) {
      pdflc[1] += pdfx1[k]*pdfx2[12-k];
      pdflc[1] += pdfx1[12-k]*pdfx2[k];
   }
   return pdflc;

}
