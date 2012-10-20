// Author: Daniel Britzger
// DESY, 11/08/2012

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO_reader_2.1.0                                                //
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

#ifndef FASTNLOCRUNDEC
#define FASTNLOCRUNDEC

#include "FastNLOReader.h"
#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include "CRunDec.h"

using namespace std;


class FastNLOCRunDec : public FastNLOReader {
private:
public:
   FastNLOCRunDec(string name);
   FastNLOCRunDec(string name, string LHAPDFfile, int PDFset = 0);

   // ----- Printout ---- //
   void PrintPDFInformation() const ;		// Print LHAPDF settings
   void PrintRunDecValues();			// Print values, which are passed to CRunDec for alpha_s evolution
   
   // ---- getters and setters LHAPDF variables ---- //
   void SetLHAPDFfilename( string filename );
   void SetLHAPDFset( int set );
   int GetIPDFSet() const {return fiPDFSet;};
   int GetNPDFSets() const {return fnPDFs;};	

   // ---- getters and setters CRunDec variables ---- //
   void   SetAlphasMz( double AlphasMz , bool ReCalcCrossSection = false);
   double GetAlphasMz() const { return fAlphasMz; };
   void   SetMz( double Mz , bool ReCalcCrossSection = false );
   double GetMz() const { return fMz; };
   void   SetNf( double nf , bool ReCalcCrossSection = false );
   int    GetNf() const { return fNf; };
   void   SetNloop( double nloop , bool ReCalcCrossSection = false );
   int    GetNloop() const {return fNloop;};

protected:
   // inherited functions
   double EvolveAlphas(double Q) const ;
   void InitPDF();
   vector<double> GetXFX(double xp, double muf) const ;

   // ---- LHAPDF vars ---- //
   string fLHAPDFfilename;
   int fnPDFs;
   int fiPDFSet;
   
   // ---- CRunDec ---- //
   // declare variables as static, that all instances use same alpha_s evolution
   static double fAlphasMz;
   static double fMz;
   static int fNf;
   static int fNloop;
   void InitReasonableRunDecValues();
   static CRunDec fcrundec;

};


double FastNLOCRunDec::fAlphasMz=0.1180000654;
double FastNLOCRunDec::fMz=0.1180000654;
int FastNLOCRunDec::fNf = 5;
int FastNLOCRunDec::fNloop=2;
CRunDec FastNLOCRunDec::fcrundec=CRunDec(FastNLOCRunDec::fNf);


//______________________________________________________________________________


FastNLOCRunDec::FastNLOCRunDec(string name) : FastNLOReader(name) {
   warn["FastNLOCRunDec"]<<"Please initialize a PDF set using SetLHAPDFfilename( PDFFile )!"<<std::endl;
   warn["FastNLOCRunDec"]<<"Also do not forget to fill the PDF cache afterwards via FillPDFCache()!"<<std::endl;
}


//______________________________________________________________________________


FastNLOCRunDec::FastNLOCRunDec(string name, string LHAPDFfile, int PDFset) : FastNLOReader(name){
   SetLHAPDFfilename(LHAPDFfile);
   SetLHAPDFset(PDFset);
   // do cross sections calculation, since everything is ready
   InitReasonableRunDecValues();
   FillAlphasCache();
   FillPDFCache();
   CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCRunDec::SetAlphasMz( double AlphasMz , bool ReCalcCrossSection ){
   debug["SetAlphasMz"]<<"Setting alpha_s(Mz)="<<AlphasMz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the alpha_s value at M_Z
   //
   fAlphasMz    = AlphasMz;             // new alpha_s value
   FillAlphasCache();
   if ( ReCalcCrossSection ) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCRunDec::SetMz( double Mz , bool ReCalcCrossSection ){
   debug["SetMz"]<<"Setting MZ-"<<Mz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the Z-Boson mass
   //
   fMz    = Mz;             // new alpha_s value
   FillAlphasCache();
   if ( ReCalcCrossSection ) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCRunDec::SetNloop( double nloop, bool ReCalcCrossSection ){
   debug["SetNloop"]<<"Setting n-loop="<<nloop<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set n loop calculation
   //
   fNloop    = nloop;             // new alpha_s value
   FillAlphasCache();
   if ( ReCalcCrossSection ) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCRunDec::SetNf( double Nf , bool ReCalcCrossSection ){
   debug["SetNf"]<<"Setting number of flavors to "<<Nf<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the number of flavors
   //
   fNf    = Nf;             // new alpha_s value
   //fcrundec.SetConstants(fNf);
   FillAlphasCache();
   if ( ReCalcCrossSection ) CalcCrossSection();
}

//______________________________________________________________________________

void FastNLOCRunDec::InitReasonableRunDecValues(){
   fAlphasMz = 0.11800000313;
   fMz = 91.1876;
   SetNf(5);
   fNloop=2;
   if ( info.GetSpeak() ) {
      info["InitReasonableRunDecValues"]<<"Printing initialized CRunDecValues."<<endl;
      PrintRunDecValues();
   }
}

   
//______________________________________________________________________________

void FastNLOCRunDec::PrintRunDecValues(){
   static const string csep41("#########################################");
   cout<<csep41<<csep41<<endl;
   cout<<"CRunDec Values: Alphas(Mz)="<<fAlphasMz
	<<"\tMZ="<<fMz
	<<"\tn-flavors="<<fNf
	<<"\tn-loop="<<fNloop<<endl;
   cout<<csep41<<csep41<<endl;
}


//______________________________________________________________________________

double FastNLOCRunDec::EvolveAlphas(double Q) const {
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   //
   // Implementation of Alpha_s evolution as function of Mu_r.
   //
   // alpha_s evolution as given by the CRunDec program by
   // Barbara Schmidt, Matthias Steinhauser;
   // see also description in
   //   K.~G.~Chetyrkin, J.~H.~Kuhn and M.~Steinhauser,
   //   ``RunDec: A Mathematica package for running and decoupling of the strong
   //   coupling and quark masses,'' Comput.\ Phys.\ Commun.\  {\bf 133} (2000) 43
   //   [arXiv:hep-ph/0004189].

//    // - do some initial print out
//    const string csep41("#########################################");
//    const string cseps = csep41 + csep41;
//    static bool first = true;
//    if ( first ) {
//       first = false;
//       cout << endl << " " << cseps << endl;
//       printf(" # ALPHAS-CRUNDEC: First call:\n");
//       cout << " " << cseps << endl;
//       //    printf(" # ALPHAS-CRUNDEC: PI              = %#18.15f\n",twopi/2.);
//       printf(" # ALPHAS-CRUNDEC: M_Z/GeV         = %#9.6f\n",Alphas::GetMz());
//       printf(" # ALPHAS-CRUNDEC: a_s(M_Z)        = %#9.6f\n",Alphas::GetAlphasMz());
//       printf(" # APLHAS-CRUNDEC: a_s loop        = %2i\n",Alphas::GetNLoop());
//       //    printf(" # APLHAS-CRUNDEC: flavor-matching = %s\n",(bFlavorMatching?"true":"false"));
//       printf(" # APLHAS-CRUNDEC: nf (M_Z)        = %2d\n",Alphas::GetNf());
//       cout << " " << cseps << endl;
//    }

   double as_crundec = fcrundec.AlphasExact(fAlphasMz, fMz,Q, fNf, fNloop);
   return as_crundec; 
}


//______________________________________________________________________________


void FastNLOCRunDec::InitPDF(){
   //
   //  Initalize some necessary LHAPDF parameters
   //
   // security, if multiple instance with different pdfs are instantiated.
   // we always reinizialized the set PDF-set.

   //LHAPDF::setVerbosity(LHAPDF::SILENT);
   LHAPDF::setVerbosity(LHAPDF::LOWKEY);
   if ( fLHAPDFfilename == ""){
     error["InitPDF"]<<"Empty LHAPDF filename! Please define a PDF set here!\n";
     exit(1);
   } else {
     // Do not use the ByName feature, destroys ease of use on the grid without LHAPDF
     //LHAPDF::initPDFSetByName(fLHAPDFfilename);
     //cout << "PDF set name " << fLHAPDFfilename << endl;
     LHAPDF::initPDFSet(fLHAPDFfilename);
     fnPDFs = LHAPDF::numberPDF();
     if ( fnPDFs < fiPDFSet ){
       error["InitPDF"]<<"There are only "<<fnPDFs<<" pdf sets within this LHAPDF file. You were looking for set number "<<fiPDFSet<<std::endl;
     }
     LHAPDF::initPDF(fiPDFSet);
   }
}


//______________________________________________________________________________



vector<double> FastNLOCRunDec::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   return LHAPDF::xfx(xp,muf);
}


//______________________________________________________________________________


void FastNLOCRunDec::SetLHAPDFfilename( string filename ) {
   fLHAPDFfilename = filename;
   // reset pdfset
   fiPDFSet = 0;
   InitPDF();
}


//______________________________________________________________________________


void FastNLOCRunDec::SetLHAPDFset( int set ) {
   fiPDFSet = set;
   InitPDF();
}


//______________________________________________________________________________


void FastNLOCRunDec::PrintPDFInformation() const {
   //
   // print out the information about the currently used LHAPDF file.
   // unfortunately there is no getter for lhapdf-filename or
   // used pdf-member-id available.
   // One must take care, that one is always using the desired pdf.
   //
   // e.g. If one has two FastNLOReader instances and one initalizes the
   // second instance with another pdf. Then also the first one is using this
   // pdf when evaluating CalcCrossSection (after a PDFCacheRefilling).
   //
   printf(" ##################################################################################\n");
   printf(" #  FastNLOCRunDec::PrintCurrentLHAPDFInformation.\n");
   printf(" #      Your currently initalized pdf is called:\n");
   LHAPDF::getDescription();
   printf(" #      Information about current PDFSet in current LHAPDF-file cannot be displayed.\n");
   printf(" #      Please use FastNLOReader::SetLHAPDFset(int) to choose a pdf-set.\n");
   printf(" ##################################################################################\n");
}



#endif
