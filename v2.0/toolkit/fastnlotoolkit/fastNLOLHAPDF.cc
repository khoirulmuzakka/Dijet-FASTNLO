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

#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "fastnlotk/fastNLOReader.h"
#include "fastnlotk/fastNLOLHAPDF.h"

using namespace std;



//______________________________________________________________________________


fastNLOLHAPDF::fastNLOLHAPDF(string name) : fastNLOReader(name), fnPDFs(0) , fiPDFMember(0) , fchksum(0.) {
   logger.info["fastNLOLHAPDF"]<<"Please initialize a PDF file using SetLHAPDFFilename( PDFFile ) and a PDF set using SetLHAPDFMember(int PDFMember)"<<std::endl;

   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   PDFSet = NULL;
   PDF = NULL;
   #endif
}
//______________________________________________________________________________


fastNLOLHAPDF::~fastNLOLHAPDF() {
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   delete PDFSet;
   delete PDF;
   #endif
}

//______________________________________________________________________________


fastNLOLHAPDF::fastNLOLHAPDF(string name, string LHAPDFFile, int PDFMember) : fastNLOReader(name), fnPDFs(0) , fiPDFMember(0) , fchksum(0.) {
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   PDFSet = NULL;
   PDF = NULL;
   #endif
   SetLHAPDFFilename(LHAPDFFile);
   SetLHAPDFMember(PDFMember);
   // Call additional initialization. Not necessary for LHAPDF.
   InitEvolveAlphas();
   // Everything set. Do cross sections calculation.
   CalcCrossSection();
}

   // Getters
   int fastNLOLHAPDF::GetIPDFMember() const {
      return fiPDFMember;
   };
   int fastNLOLHAPDF::GetNPDFMembers() const {
      return fnPDFs;
   };
   int fastNLOLHAPDF::GetNPDFMaxMember() const {
      return fnPDFs-1;
   };


//______________________________________________________________________________


double fastNLOLHAPDF::EvolveAlphas(double Q) const {
   //debug<<"EvolveAlphas with Q="<<Q<<endl;
   //
   // Implementation of Alpha_s evolution as function of Mu_r only.
   //
   // the alpha_s evolution is done within LHAPDF.
   //
   // WARNING: You cannot change alpha_s(Mz), but is is
   // defined with the pdf. 'alphasMz' is not used here!
   //
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   return PDF->alphasQ(Q);
   #else
   return LHAPDF::alphasPDF(Q);
   #endif
}


//______________________________________________________________________________


bool fastNLOLHAPDF::InitPDF() {
   //
   //  Initalize some necessary LHAPDF parameters
   //  return true, if successful initialization
   //  return false, if PDF initialization failed
   //
   // LHAPDF interface:
   // security, if multiple instance with different pdfs are instantiated.
   // we always reinitialized the set PDF-set.

   if (fLHAPDFFilename == "") {
      logger.warn["InitPDF"]<<"Empty LHAPDF filename! Please define a PDF set here!\n";
      return false;
   }

   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   // Suppress LHAPDF6 output at each member change
   LHAPDF::setVerbosity(0);
   return true;

   #else
   //LHAPDF::setVerbosity(LHAPDF::SILENT);
   LHAPDF::setVerbosity(LHAPDF::LOWKEY);
   // Do not use the ByName feature, destroys ease of use on the grid without LHAPDF
   //LHAPDF::initPDFSetByName(fLHAPDFFilename);
   //cout << "PDF set name " << fLHAPDFFilename << endl;
   if (fchksum == 0 || fchksum != CalcChecksum(1.)) {
      // need to reset LHAPDF.
      logger.debug["InitPDF"]<<"Need to reset lhapdf. fchksum="<<fchksum<<"\tCalcChecksum(1.)="<<CalcChecksum(1.)<<endl;
      LHAPDF::initPDFSet(fLHAPDFFilename);
      fnPDFs = LHAPDF::numberPDF()+1; // LHAPDF counts 0-44 and returns, 44 which must be 45
      if (fnPDFs < fiPDFMember+1) {
         logger.error["InitPDF"]<<"There are only "<<fnPDFs<<" pdf sets within this LHAPDF file. You were looking for set number "<<fiPDFMember<<std::endl;
         return false;
      }
      LHAPDF::initPDF(fiPDFMember);
   }
   fchksum = CalcChecksum(1.);
   return true;
   #endif
}


//______________________________________________________________________________



vector<double> fastNLOLHAPDF::GetXFX(double xp, double muf) const {
   //
   //  GetXFX is used to get the parton array from the
   //  pre-defined pdf-interface.
   //
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   // vector<double> xfx(13);
   // PDF->xfxQ(xp, muf, xfx);
   vector <double> xfx(13);
   PDF->xfxQ(xp,muf,xfx);
   return xfx;
   #else
   return LHAPDF::xfx(xp,muf);
   #endif
}


//______________________________________________________________________________


void fastNLOLHAPDF::SetLHAPDFFilename(string filename) {
   if (filename != fLHAPDFFilename) fchksum = 0;
   fLHAPDFFilename = filename;
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   if ( PDFSet ) delete PDFSet;
   PDFSet = new LHAPDF::PDFSet(filename);
   fnPDFs = PDFSet->size();
   SetLHAPDFMember(0);
   #else
   // Reset pdfset member to zero
   fiPDFMember = 0;
   // KR: Reactivated this. Why was it switched off?
   // --> Mass settings etc. can be read from LHAPDF after setting the filename, i.e. the set.
   InitPDF();
   #endif
}


//______________________________________________________________________________


void fastNLOLHAPDF::SetLHAPDFMember(int set) {
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   if ( PDF ) delete PDF;
   PDF = PDFSet->mkPDF(set);
   #else
   fiPDFMember = set;
   if (fchksum == CalcChecksum(1.)) {  // nothin has changed? we set only the pdfmember
      logger.debug["SetLHAPDFMember"]<<"Changing only pdfmember!"<<endl;
      LHAPDF::initPDF(fiPDFMember);
      fchksum = CalcChecksum(1.);
   } else  {
      logger.debug["SetLHAPDFMember"]<<"Demanding full re-initalization of PDF."<<endl;
      fchksum = 0;
   }
   //InitPDF();
   #endif
}


//______________________________________________________________________________


void fastNLOLHAPDF::PrintPDFInformation() const {
   //
   // print out the information about the currently used LHAPDF file.
   // unfortunately there is no getter for lhapdf-filename or
   // used pdf-member-id available.
   // One must take care, that one is always using the desired pdf.
   //
   // e.g. If one has two fastNLOReader instances and one initalizes the
   // second instance with another pdf. Then also the first one is using this
   // pdf when evaluating CalcCrossSection (after a PDFCacheRefilling).
   //

   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   cout << PDFSet->description();
   #else
   printf(" ##################################################################################\n");
   printf(" #  fastNLOLHAPDF::PrintCurrentLHAPDFInformation.\n");
   printf(" #      Your currently initalized pdf is called:\n");
   LHAPDF::getDescription();
   printf(" #      Information about current PDFMember in current LHAPDF-file cannot be displayed.\n");
   printf(" #      Please use fastNLOReader::SetLHAPDFMember(int) to choose a pdf-set.\n");
   printf(" ##################################################################################\n");
   #endif
}

void fastNLOLHAPDF::SetMz(double Mz) {
   logger.warn["SetMz"]<<"WARNING! The Z mass cannot be changed in alpha_s evolution of LHAPDF!"<<endl;
}

void fastNLOLHAPDF::SetQMass(int pdgid, double mq) {
   logger.warn["SetQMass"]<<"WARNING! The quark masses cannot be changed in alpha_s evolution of LHAPDF!"<<endl;
}

void fastNLOLHAPDF::SetNFlavor(int nflavor) {
   logger.warn["SetNFlavor"]<<"WARNING! The no. of active flavors cannot be changed in alpha_s evolution of LHAPDF!"<<endl;
}

void fastNLOLHAPDF::SetNLoop(int nloop) {
   logger.warn["SetNLoop"]<<"WARNING! The no. of loops cannot be changed in alpha_s evolution of LHAPDF!"<<endl;
}

void fastNLOLHAPDF::SetAlphasMz(double AlphasMz, bool ReCalcCrossSection) {
   logger.warn["SetAlphasMz"]<<"WARNING! alpha_s(M_Z) cannot be changed in alpha_s evolution of LHAPDF!"<<endl;
}

void fastNLOLHAPDF::InitEvolveAlphas() {
   // For LHAPDF do nothing
}

double fastNLOLHAPDF::GetQMass(int pdgid) const {
   if (pdgid < 1 || pdgid > 6 ) {
      logger.error["GetQMass"]<<"PDG code out of quark range 1-6! Aborted\n";
      exit(1);
   }
   return LHAPDF::getQMass(pdgid);
}

int fastNLOLHAPDF::GetNLoop() const {
   return (LHAPDF::getOrderAlphaS() + 1);
}

int fastNLOLHAPDF::GetNFlavor() const {
   return (LHAPDF::getNf());
}

double fastNLOLHAPDF::GetAlphasMz() const {
   #if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   return PDF->alphasQ(91.1876);
   #endif
   return LHAPDF::alphasPDF(91.1876);
}


#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
//______________________________________________________________________________
vector<LHAPDF::PDFUncertainty>  fastNLOLHAPDF::GetPDFUncertaintyLHAPDF(double cl, bool alternative) {
   //! Calculate PDF uncertainty on cross sections.
   //! Formulae for PDF uncertainties are taken from LHAPDF6.
   //! Function returns struct with all values.
   //!   PDFUncertainty.central [new central value]
   //!   PDFUncertainty.errplus
   //!   PDFUncertainty.errminus
   //!   PDFUncertainty.errsymm
   //!   PDFUncertainty.scale   [scalefactor]
   //! More documentation. See docu of LHAPDF::PDFSet::uncertainty()
   //!
   //! 'cl' is used to rescale uncertainties to a particular confidence level
   //!
   //! If the PDF set is given in the form of replicas, then optional argument
   //! 'alternative' equal to true (default: false) will construct a confidence
   //! interval from the probability distribution of replicas, with the central
   //! value given by the median.

   vector<LHAPDF::PDFUncertainty> PDFUnc;
   const unsigned int nMem = GetNPDFMembers();
   const unsigned int nObsBins = GetNObsBin();

   int iMem0 = fiPDFMember;
   vector<vector<double> > CSs(nObsBins); // [bin][iMem]
   for ( unsigned int iObs = 0 ; iObs<nObsBins ; iObs++ ) CSs[iObs].resize(nMem);
   for ( unsigned int iMem = 0 ; iMem<nMem ; iMem++ ) {
      SetLHAPDFMember(iMem);
      CalcCrossSection();
      vector<double> xs = GetCrossSection();
      for ( unsigned int iObs = 0 ; iObs<nObsBins ; iObs++ ) CSs[iObs][iMem] = xs[iObs];
   }
   for ( unsigned int iObs = 0 ; iObs<nObsBins ; iObs++ ) {
      PDFUnc.push_back(PDFSet->uncertainty(CSs[iObs],cl,alternative));
   }
   SetLHAPDFMember(iMem0);
   CalcCrossSection();
   return PDFUnc;
}


//______________________________________________________________________________
vector<double>  fastNLOLHAPDF::CalcPDFUncertaintyMinus(const vector<LHAPDF::PDFUncertainty>& PDFUnc) const{
   //!get vector<double> for PDF-minus uncertainty
   //! uncertainties are 'positive'!
   vector<double> ret(GetNObsBin());
   for ( unsigned int iObs = 0 ; iObs<GetNObsBin() ; iObs++ )
      ret[iObs] = PDFUnc[iObs].errminus;
   return ret;
}


//______________________________________________________________________________
vector<double>  fastNLOLHAPDF::CalcPDFUncertaintyPlus(const vector<LHAPDF::PDFUncertainty>& PDFUnc) const {
   //!get vector<double> for PDF-plus uncertainty
   vector<double> ret(GetNObsBin());
   for ( unsigned int iObs = 0 ; iObs<GetNObsBin() ; iObs++ )
      ret[iObs] = PDFUnc[iObs].errplus;
   return ret;
}


//______________________________________________________________________________
vector<double>  fastNLOLHAPDF::CalcPDFUncertaintyRelMinus(const vector<LHAPDF::PDFUncertainty>& PDFUnc) const{
   //!get vector<double> for PDF-minus uncertainty
   //! uncertainties are 'negative'!
   //! Relative uncertainties are calcuated with 'median' cross section in case of MC uncertainties
   vector<double> ret(GetNObsBin());
   for ( unsigned int iObs = 0 ; iObs<GetNObsBin() ; iObs++ )
      ret[iObs] = PDFUnc[iObs].errminus*-1 / PDFUnc[iObs].central;
   return ret;
}


//______________________________________________________________________________
vector<double>  fastNLOLHAPDF::CalcPDFUncertaintyRelPlus(const vector<LHAPDF::PDFUncertainty>& PDFUnc) const {
   //!get vector<double> for PDF-plus uncertainty
   //! Relative uncertainties are calcuated with 'median' cross section in case of MC uncertainties
   vector<double> ret(GetNObsBin());
   for ( unsigned int iObs = 0 ; iObs<GetNObsBin() ; iObs++ )
      ret[iObs] = PDFUnc[iObs].errplus / PDFUnc[iObs].central;
   return ret;
}


//______________________________________________________________________________
vector<double>  fastNLOLHAPDF::CalcPDFUncertaintyCentral(const vector<LHAPDF::PDFUncertainty>& PDFUnc) const {
   //!get vector<double> for new central value
   vector<double> ret(GetNObsBin());
   for ( unsigned int iObs = 0 ; iObs<GetNObsBin() ; iObs++ )
      ret[iObs] = PDFUnc[iObs].central;
   return ret;
}


//______________________________________________________________________________
vector<double>  fastNLOLHAPDF::CalcPDFUncertaintySymm(const vector<LHAPDF::PDFUncertainty>& PDFUnc) const {
   //!get vector<double> for symmetrized PDF uncertainty
   vector<double> ret(GetNObsBin());
   for ( unsigned int iObs = 0 ; iObs<GetNObsBin() ; iObs++ )
      ret[iObs] = PDFUnc[iObs].errsymm;
   return ret;
}
#endif

//______________________________________________________________________________
fastNLOReader::XsUncertainty fastNLOLHAPDF::GetPDFUncertainty(const EPDFUncertaintyStyle ePDFUnc) {
   fastNLOReader::XsUncertainty XsUnc;
   unsigned int NObsBin = GetNObsBin();
   unsigned int nMem = GetNPDFMaxMember();
   vector < double > xs0;
   vector < double > dxseigu;
   vector < double > dxseigl;

   // Check input
   logger.debug["GetPDFUncertainty"]<<"ePDFUnc = "<<ePDFUnc<<endl;
   if ( ePDFUnc == kPDFNone ) {
      logger.info["GetPDFUncertainty"]<<"No PDF uncertainty, only averaged cross section result evaluated (correct for NNPDF, wrong otherwise!)."<<endl;
   } else if ( ePDFUnc == kHessianSymmetric ) {
      logger.info["GetPDFUncertainty"]<<"Calculating symmetric Hessian PDF uncertainties."<<endl;
   } else if ( ePDFUnc == kHessianAsymmetric ) {
      logger.info["GetPDFUncertainty"]<<"Calculating asymmetric Hessian PDF uncertainties."<<endl;
   } else if ( ePDFUnc == kHessianAsymmetricMax ) {
      logger.info["GetPDFUncertainty"]<<"Calculating asymmetric Hessian PDF uncertainties considering maximal pairwise deviations per eigenvector."<<endl;
   } else if ( ePDFUnc == kHessianCTEQCL68 ) {
      logger.info["GetPDFUncertainty"]<<"Calculating pairwise asymmetric Hessian PDF uncertainties rescaled to CL68 (for CTEQ PDFs)."<<endl;
   } else if ( ePDFUnc == kMCSampling ) {
      logger.info["GetPDFUncertainty"]<<"Calculating statistical sampling PDF uncertainties."<<endl;
   } else if ( ePDFUnc == kHeraPDF10 ) {
      logger.error["GetPDFUncertainty"]<<"ERROR! HERAPDF1.0 uncertainty style not yet implemented, exiting."<<endl;
      logger.error["GetPDFUncertainty"]<<"ePDFUnc = "<<ePDFUnc<<endl;
      exit(1);
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
   } else if ( ePDFUnc == kLHAPDF6 ) {
      logger.info["GetPDFUncertainty"]<<"Calculating LHAPDF6 PDF uncertainties."<<endl;
#endif
   } else {
      logger.error["GetPDFUncertainty"]<<"ERROR! Selected PDF uncertainty style not yet implemented, exiting."<<endl;
      logger.error["GetPDFUncertainty"]<<"ePDFUnc = "<<ePDFUnc<<endl;
      exit(1);
   }
   unsigned int nEig = nMem/2;
   logger.info["GetPDFUncertainty"]<<"Info: Number of highest PDF set member nMem = " << nMem << endl;
   logger.info["GetPDFUncertainty"]<<"Info: Guessed number of eigen vectors nMem/2 = " << nEig << endl;
   if ( nMem < 2 ) {
      logger.error["GetPDFUncertainty"]<<"ERROR! This PDF set has only one or two members: nMem = " << nMem << endl;
      logger.error["GetPDFUncertainty"]<<"PDF uncertainty calculation impossible, aborted!" << endl;
      exit(1);
   } else if ( nMem%2 == 1 && (ePDFUnc == kHessianAsymmetric || ePDFUnc == kHessianAsymmetricMax || ePDFUnc == kHessianCTEQCL68) ) {
      logger.error["GetPDFUncertainty"]<<"ERROR! Odd number of PDF members found: nMem = " << nMem << endl;
      logger.error["GetPDFUncertainty"]<<"This cannot work with selected asymmetric Hessian uncertainties, aborted!" << endl;
      exit(1);
   }

   if ( ePDFUnc != kLHAPDF6 ) {

      // Evaluate central/zeroth PDF set member
      SetLHAPDFMember(0);
      CalcCrossSection();
      xs0 = GetCrossSection();
      // Initialize cross section xs, uncertainties dxs, and temporary vector dxseig to zero
      for ( unsigned int iobs = 0; iobs < NObsBin; iobs++ ) {
         XsUnc.xs.push_back(0);
         XsUnc.dxsu.push_back(0);
         XsUnc.dxsl.push_back(0);
         dxseigu.push_back(0);
         dxseigl.push_back(0);
      }

      // Pairwise loop over other PDF set members, sum up weights and shifted weights squared ((weights-xs0)^2)
      // Shifted weights allow more precise variance calculation
      for ( unsigned int ieig = 1; ieig <= nEig; ieig++ ) {
         SetLHAPDFMember(2*ieig-1);
         CalcCrossSection();
         for ( unsigned int iobs = 0; iobs < NObsBin; iobs++ ) {
            XsUnc.xs[iobs] = XsUnc.xs[iobs] + XSection[iobs];
            dxseigu[iobs] = max(0.,XSection[iobs]-xs0[iobs]);
            dxseigl[iobs] = min(0.,XSection[iobs]-xs0[iobs]);
         }
         SetLHAPDFMember(2*ieig);
         CalcCrossSection();
         for ( unsigned int iobs = 0; iobs < NObsBin; iobs++ ) {
            XsUnc.xs[iobs] = XsUnc.xs[iobs] + XSection[iobs];
            // Take only maximal one in case of one-sided deviations for one eigenvector
            if ( ePDFUnc == kHessianAsymmetricMax || ePDFUnc == kHessianCTEQCL68 ) {
               dxseigu[iobs] = pow(max(dxseigu[iobs],XSection[iobs]-xs0[iobs]),2.);
               dxseigl[iobs] = pow(min(dxseigl[iobs],XSection[iobs]-xs0[iobs]),2.);
            }
            // Add up both even in case of one-sided deviations for one eigenvector
            else {
               dxseigu[iobs] = pow(dxseigu[iobs],2.) + pow(max(0.,XSection[iobs]-xs0[iobs]),2.);
               dxseigl[iobs] = pow(dxseigl[iobs],2.) + pow(min(0.,XSection[iobs]-xs0[iobs]),2.);
            }
            // Fill into dxs
            XsUnc.dxsu[iobs] = XsUnc.dxsu[iobs] + dxseigu[iobs];
            XsUnc.dxsl[iobs] = XsUnc.dxsl[iobs] + dxseigl[iobs];
         }
      }

      // Evaluate last PDF set member in case of odd numbers (not possible with asymmetric Hessian!)
      if ( nMem%2 == 1 ) {
         SetLHAPDFMember(nMem);
         CalcCrossSection();
         for ( unsigned int iobs = 0; iobs < NObsBin; iobs++ ) {
            XsUnc.xs[iobs] = XsUnc.xs[iobs] + XSection[iobs];
            XsUnc.dxsu[iobs] = XsUnc.dxsu[iobs] + pow(max(0.,XSection[iobs]-xs0[iobs]),2.);
            XsUnc.dxsl[iobs] = XsUnc.dxsl[iobs] + pow(min(0.,XSection[iobs]-xs0[iobs]),2.);
         }
      }

      // Derive chosen relative uncertainties
      for ( unsigned int iobs = 0; iobs < NObsBin; iobs++ ) {
         // No PDF uncertainty, only averaged cross section result evaluated (Correct for NNPDF, wrong otherwise!).
         if ( ePDFUnc == kPDFNone ) {
            XsUnc.xs[iobs]   = XsUnc.xs[iobs]/nMem;
            XsUnc.dxsu[iobs] = 0.;
            XsUnc.dxsl[iobs] = 0.;
         }
         // Uncertainty is +- sqrt (sum of all summed squared deviations)
         else if ( ePDFUnc == kHessianSymmetric ) {
            XsUnc.xs[iobs]   = xs0[iobs];
            XsUnc.dxsu[iobs] = +sqrt(XsUnc.dxsu[iobs] + XsUnc.dxsl[iobs]);
            XsUnc.dxsl[iobs] = -XsUnc.dxsu[iobs];
         }
         // Uncertainty is sqrt (separately summed upper and lower squared deviations)
         else if ( ePDFUnc == kHessianAsymmetric || ePDFUnc == kHessianAsymmetricMax || ePDFUnc == kHessianCTEQCL68 ) {
            XsUnc.xs[iobs]   = xs0[iobs];
            XsUnc.dxsu[iobs] = +sqrt(XsUnc.dxsu[iobs]);
            XsUnc.dxsl[iobs] = -sqrt(XsUnc.dxsl[iobs]);
            if ( ePDFUnc == kHessianCTEQCL68 ) {
               XsUnc.dxsu[iobs] = XsUnc.dxsu[iobs]/TOCL90;
               XsUnc.dxsl[iobs] = XsUnc.dxsl[iobs]/TOCL90;
            }
         }
         // Central value xs0 is replaced by sampling average; uncertainty is sqrt of sampling variance
         else if ( ePDFUnc == kMCSampling ) {
            XsUnc.xs[iobs]   = XsUnc.xs[iobs]/nMem;
            XsUnc.dxsu[iobs] = +sqrt((nMem/(nMem-1)) * ( (XsUnc.dxsu[iobs]+XsUnc.dxsl[iobs])/nMem -
                                                         pow(XsUnc.xs[iobs],2.) - pow(xs0[iobs],2.) + 2.*XsUnc.xs[iobs]*xs0[iobs] ));
            XsUnc.dxsl[iobs] = -XsUnc.dxsu[iobs];
         }
         // HERAPDF not yet implemented
         else {
            logger.error["GetPDFUncertainty"]<<"ERROR! Selected PDF uncertainty style not yet implemented, exiting."<<endl;
            logger.error["GetPDFUncertainty"]<<"ePDFUnc = "<<ePDFUnc<<endl;
            exit(1);
         }
         // Give back +- relative uncertainties
         if ( abs(XsUnc.xs[iobs]) > DBL_MIN ) {
            XsUnc.dxsu[iobs] = XsUnc.dxsu[iobs] / XsUnc.xs[iobs];
            XsUnc.dxsl[iobs] = XsUnc.dxsl[iobs] / XsUnc.xs[iobs];
         } else {
            XsUnc.dxsu[iobs] = 0.;
            XsUnc.dxsl[iobs] = 0.;
         }
         logger.debug["GetPDFUncertainty"]<<"iobs = " << iobs << "dxsl = " << XsUnc.dxsl[iobs] << ", dxsu = " << XsUnc.dxsu[iobs] <<endl;
      }
   } else {
#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6
      // Default is no alternative for NNPDF replica uncertainties
      //      vector<LHAPDF::PDFUncertainty> L6PDFUnc = GetPDFUncertaintyLHAPDF(100*erf(1/sqrt(2)),true);
      vector<LHAPDF::PDFUncertainty> L6PDFUnc = GetPDFUncertaintyLHAPDF();
      XsUnc.xs   = CalcPDFUncertaintyCentral(L6PDFUnc);
      XsUnc.dxsu = CalcPDFUncertaintyRelPlus(L6PDFUnc);
      XsUnc.dxsl = CalcPDFUncertaintyRelMinus(L6PDFUnc);
#endif
   }

   logger.info["GetPDFUncertainty"]<<"Setting PDF member back to default of zero."<<endl;
   SetLHAPDFMember(0);

   return XsUnc;
}
