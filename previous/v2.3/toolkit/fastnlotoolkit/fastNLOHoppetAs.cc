// Author: Daniel Britzger
// DESY, 20/04/2012

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  fastNLO_toolkit                                                     //
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

//////////////////////////////////////////////////////////////////////////
//
//  fastNLOAlphas
//  This class inherits the PDF interface from
//  fastNLOLHAPDF, while the alpha_s evolution
//  is superseeded by the Alphas.h class.
//
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <LHAPDF/LHAPDF.h>
#include "fastnlotk/fastNLOLHAPDF.h"
#include "fastnlotk/fastNLOHoppetAs.h"
#include "hoppet_v1.h"

using namespace std;



//______________________________________________________________________________
//
fastNLOHoppetAs::fastNLOHoppetAs(std::string name, std::string LHAPDFFile, int PDFMem) : fastNLOHoppet(name,LHAPDFFile,PDFMem) {
   // Without PDF info use PDG values as default
   //   SetPDGValues();
   //   PrintParmValues();
   // Set initial values via LHAPDF6 info system
   SetLHAPDFValues(LHAPDFFile, PDFMem);
   // Print out values for checking
   //   PrintParmValues();
   fastNLOHoppet::InitPDF();
};



// Evolution
std::vector<double> fastNLOHoppetAs::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the pre-defined pdf-interface.
   //
   return fastNLOLHAPDF::GetXFX(xp, muf);
}
