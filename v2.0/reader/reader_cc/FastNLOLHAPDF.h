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

#ifndef FASTNLOLHAPDF
#define FASTNLOLHAPDF

#include "FastNLOReader.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
#include "speaker.h"

using namespace std;


class FastNLOLHAPDF : public FastNLOReader {

private:
   FastNLOLHAPDF(string name);
public:
   FastNLOLHAPDF(string name, string LHAPDFfile, int PDFset = 0);

   void SetLHAPDFfilename( string filename );
   void SetLHAPDFset( int set );
   int GetIPDFSet() const {return fiPDFSet;};
   int GetNPDFSets() const {return fnPDFs;};
   void PrintPDFInformation() const ;

protected:
   // inherited functions
   double EvolveAlphas(double Q) const ;
   void InitPDF();
   vector<double> GetXFX(double xp, double muf) const ;

   // ---- LHAPDF vars ---- //
   string fLHAPDFfilename;
   int fnPDFs;
   int fiPDFSet;
  
};



//______________________________________________________________________________


FastNLOLHAPDF::FastNLOLHAPDF(string name) : FastNLOReader(name) {
   // --- fastNLO user: if you have interface FastNLOLHAPDF::EvolveAlphas(double,double)
   //     it is convenient to automatically interface it here.
   warn["FastNLOLHAPDF"]<<"Please set LHAPDFfilename and init the LHAPDF::PDFset."<<std::endl;
}



//______________________________________________________________________________


FastNLOLHAPDF::FastNLOLHAPDF(string name, string LHAPDFfile, int PDFset) : FastNLOReader(name){
   SetLHAPDFfilename(LHAPDFfile);
   SetLHAPDFset(PDFset);
   // do cross sections calculation, since everything is yet ready
   FillAlphasCache();
   FillPDFCache();
   CalcCrossSection();
}



//______________________________________________________________________________


double FastNLOLHAPDF::EvolveAlphas(double Q) const {
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'alphasMz' is not used here!
   //
   return LHAPDF::alphasPDF(Q);
}


//______________________________________________________________________________


void FastNLOLHAPDF::InitPDF(){
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



vector<double> FastNLOLHAPDF::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   return LHAPDF::xfx(xp,muf);
}


//______________________________________________________________________________


void FastNLOLHAPDF::SetLHAPDFfilename( string filename ) {
   fLHAPDFfilename = filename;
   // reset pdfset
   fiPDFSet = 0;
   InitPDF();
}


//______________________________________________________________________________


void FastNLOLHAPDF::SetLHAPDFset( int set ) {
   fiPDFSet = set;
   InitPDF();
}


//______________________________________________________________________________


void FastNLOLHAPDF::PrintPDFInformation() const {
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
   printf(" #  FastNLOLHAPDF::PrintCurrentLHAPDFInformation.\n");
   printf(" #      Your currently initalized pdf is called:\n");
   LHAPDF::getDescription();
   printf(" #      Information about current PDFSet in current LHAPDF-file cannot be displayed.\n");
   printf(" #      Please use FastNLOReader::SetLHAPDFset(int) to choose a pdf-set.\n");
   printf(" ##################################################################################\n");
}



#endif
