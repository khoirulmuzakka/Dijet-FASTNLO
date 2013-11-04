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

#ifndef FASTNLOLHAPDF6
#define FASTNLOLHAPDF6

#include "FastNLOReader.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
#include "speaker.h"


using namespace std;


class FastNLOLHAPDF6 : public FastNLOReader {

private:
public:
   FastNLOLHAPDF6(string name);
   FastNLOLHAPDF6(string name, string LHAPDFfile, int PDFSet = 0);

   void SetLHAPDFFilename(string filename);
   void SetLHAPDFMember(int set);
   int GetIPDFMember() const {
      return fiPDFMember;
   };
   int GetNPDFMembers() const {
      return fnPDFs;
   };
   int GetNPDFMaxMember() const {
      return fnPDFs-1;
   };
   void PrintPDFInformation() const ;

protected:
   // inherited functions
   double EvolveAlphas(double Q) const ;
   bool InitPDF();
   vector<double> GetXFX(double xp, double muf) const ;

   // ---- LHAPDF vars ---- //
   string fLHAPDFFilename;
   LHAPDF::PDFSet* PDFSet;
   LHAPDF::PDF* PDF;
   int fnPDFs;
   int fiPDFMember;

   double fchksum;


};



//______________________________________________________________________________


FastNLOLHAPDF6::FastNLOLHAPDF6(string name) : FastNLOReader(name) , fnPDFs(0) , fiPDFMember(0) , fchksum(0.) {
   info["FastNLOLHAPDF"]<<"Please initialize a PDF file using SetLHAPDFFilename( PDFFile ) and a PDF set using SetLHAPDFMember(int PDFMember)"<<std::endl;
}


//______________________________________________________________________________


FastNLOLHAPDF6::FastNLOLHAPDF6(string name, string LHAPDFFile, int PDFMember) : FastNLOReader(name) , fchksum(0.) {
   SetLHAPDFFilename(LHAPDFFile);
   SetLHAPDFMember(PDFMember);
}



//______________________________________________________________________________


double FastNLOLHAPDF6::EvolveAlphas(double Q) const {
   //debug<<"EvolveAlphas with Q="<<Q<<endl;
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'alphasMz' is not used here!
   //
   return PDF->alphasQ(Q);
}



//______________________________________________________________________________



vector<double> FastNLOLHAPDF6::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   vector <double> xfx;
   for (int id=-6; id<7; id++) {
      xfx.push_back(PDF->xfxQ(id, xp, muf));
   }
   return xfx;
}


//______________________________________________________________________________


void FastNLOLHAPDF6::SetLHAPDFFilename(string filename) {
   PDFSet = new LHAPDF::PDFSet(filename);
   fnPDFs = PDFSet->size();
}


//______________________________________________________________________________


void FastNLOLHAPDF6::SetLHAPDFMember(int set) {
   PDF = PDFSet->mkPDF(set);
}


//______________________________________________________________________________


void FastNLOLHAPDF6::PrintPDFInformation() const {
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
   cout << PDFSet->description();
}

bool FastNLOLHAPDF6::InitPDF() {
return true;
}

#endif
