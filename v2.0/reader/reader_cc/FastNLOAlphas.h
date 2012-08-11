// Author: Daniel Britzger
// DESY, 20/04/2012

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

#ifndef FASTNLOALPHAS
#define FASTNLOALPHAS

#include "FastNLOReader.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
#include "speaker.h"
#include "Alphas.h"

using namespace std;


class FastNLOAlphas : public FastNLOReader {

private:
   FastNLOAlphas(string name);
public:
   FastNLOAlphas(string name, string LHAPDFfile, int PDFset = 0);

   void SetLHAPDFfilename( string filename );
   void SetLHAPDFset( int set );
   int GetIPDFSet() const {return fiPDFSet;};
   int GetNPDFSets() const {return fnPDFs;};
   void PrintPDFInformation() const ;

   // ---- Alphas vars ---- //
   void SetAlphasMz( double AlphasMz , bool ReCalcCrossSection = false );
   double GetAlphasMz() const { return fAlphasMz; };
   void SetGRVtoPDG2011_2loop();


protected:
   // inherited functions
   //double EvolveAlphas(double Q, double alphasMz ) const ;
   double EvolveAlphas(double Q) const ;
   void InitPDF();
   vector<double> GetXFX(double xp, double muf) const ;

   // ---- LHAPDF vars ---- //
   string fLHAPDFfilename;
   int fnPDFs;
   int fiPDFSet;

   // ---- Alphas vars ---- //
   double fAlphasMz;

};



//______________________________________________________________________________


FastNLOAlphas::FastNLOAlphas(string name) : FastNLOReader(name) {
   // --- fastNLO user: if you have interface FastNLOAlphas::EvolveAlphas(double,double)
   //     it is convenient to automatically interface it here.
   warn["FastNLOAlphas"]<<"Please set LHAPDFfilename and init the LHAPDF::PDFset."<<std::endl;
}



//______________________________________________________________________________


FastNLOAlphas::FastNLOAlphas(string name, string LHAPDFfile, int PDFset) : FastNLOReader(name){
   // we need a PDF
   SetLHAPDFfilename(LHAPDFfile);
   SetLHAPDFset(PDFset);
   // set some reasonable starting values;
   SetGRVtoPDG2011_2loop();
   fAlphasMz = 0.118000000321;
   // do cross sections calculation, since everything is yet ready
   FillAlphasCache();
   FillPDFCache();
   CalcCrossSection();
}


//______________________________________________________________________________


void FastNLOAlphas::SetAlphasMz( double AlphasMz , bool ReCalcCrossSection ){
   debug["SetAlphasMz"]<<"Setting alpha_s(Mz)="<<AlphasMz<<" and RecalculateCrossSection="<<(ReCalcCrossSection?"Yes":"No")<<endl;
   //
   //  Set the alpha_s value at M_Z
   //
   fAlphasMz    = AlphasMz;             // new alpha_s value
   FillAlphasCache();
   if ( ReCalcCrossSection ) CalcCrossSection();
}


//______________________________________________________________________________


double FastNLOAlphas::EvolveAlphas(double Q) const {
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'alphasMz' is not used here!
   //
   return Alphas::CalcAlphasMu ( Q , fAlphasMz );
}


//______________________________________________________________________________


void FastNLOAlphas::SetGRVtoPDG2011_2loop(){
   info["SetGrVtoPDF2011"]<<"Resetting GRV Alphas::Alphas evolution."<<endl;
   Alphas::SetMz(91.1876); // PDG 2011
   Alphas::SetNf(5);
   Alphas::SetNLoop(2);
   Alphas::SetFlavorMatchingOn(true);
   if ( info.GetSpeak() ) {
      info<<"Calling Alphas::PrintInfo()."<<endl;
      info<<"Alpha_s(Mz) value is taken from FastNLOAlphas, instead of Alphas::Alphas."<<endl;
      Alphas::PrintInfo();
   }
}


//______________________________________________________________________________


void FastNLOAlphas::InitPDF(){
   //
   //  Initalize some necessary LHAPDF parameters
   //
   // security, if multiple instance with different pdfs are instantiated.
   // we always reinizialized the set PDF-set.

   if ( fLHAPDFfilename == ""){
      error["InitPDF"]<<"You must specify a LHAPDF filename first.\n"; exit(1);
   }

   LHAPDF::setVerbosity(LHAPDF::SILENT);
   //LHAPDF::setVerbosity(LHAPDF::LOWKEY);
   //cout << " * LHAPDF version: " << LHAPDF::getVersion() <<endl;
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


//______________________________________________________________________________



vector<double> FastNLOAlphas::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   return LHAPDF::xfx(xp,muf);
}


//______________________________________________________________________________


void FastNLOAlphas::SetLHAPDFfilename( string filename ) {
   fLHAPDFfilename = filename;
   // reset pdfset
   fiPDFSet = 0;
   InitPDF();
}


//______________________________________________________________________________


void FastNLOAlphas::SetLHAPDFset( int set ) {
   fiPDFSet = set;
   InitPDF();
}


//______________________________________________________________________________


void FastNLOAlphas::PrintPDFInformation() const {
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
   printf(" #  FastNLOAlphas::PrintCurrentLHAPDFInformation.\n");
   printf(" #      Your currently initalized pdf is called:\n");
   LHAPDF::getDescription();
   printf(" #      Information about current PDFSet in current LHAPDF-file cannot be displayed.\n");
   printf(" #      Please use FastNLOReader::SetLHAPDFset(int) to choose a pdf-set.\n");
   printf(" ##################################################################################\n");
}



#endif
