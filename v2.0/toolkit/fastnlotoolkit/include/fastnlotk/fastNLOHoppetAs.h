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

#ifndef FASTNLOHOPPETAS
#define FASTNLOHOPPETAS

#include "fastNLOHoppet.h"



class fastNLOHoppetAs : public fastNLOHoppet {

 public:
   // Only allow constructor with information on LHAPDF set and member
   // as needed for HOPPET initialisation.
   fastNLOHoppetAs(std::string name) = delete;
   fastNLOHoppetAs(std::string name, std::string LHAPDFFile, int PDFMem = 0);

 protected:
   virtual std::vector<double> GetXFX(double xp, double muf) const;

};

#endif
