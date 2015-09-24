#ifndef __fnlogeneratorconstants__
#define __fnlogeneratorconstants__

#include <string>
#include <vector>

namespace fastNLO {

   struct GeneratorConstants {
      //! GeneratorConstants
      //!
      //! Collection of generator specific constants.
      //! These are:
      //!  - name and version of generator
      //!  - references for generator
      //!  - (additional information about generator may be included in References)
      std::string Name;                      //!< Name and version of generator
      std::vector<std::string> References;   //!< References for generator. Include additional information here (e.g. 'run-mode' or process).
      //! Get 'CodeDescription' usable for fastNLO table
      std::vector<std::string > GetCodeDescription() {
         std::vector<std::string > CodeDescr(References.size()+1);
         CodeDescr[0] = Name;
         for ( unsigned int i = 0 ; i<References.size() ; i++ )
            CodeDescr[i+1] = References[i];
         return CodeDescr;
      }
   };


   struct ProcessConstants {
      //! ProcessConstants
      //!
      //! Collection of process specific constants.
      //! Please see fastNLO table format definition for a detailed explanation.
      //!
      int LeadingOrder;          //! Order in alpha_s of leading order process
      int UnitsOfCoefficients;   //! X section units of coefficients passed to fastNLO (neg. power of 10: pb->12, fb->15)
      int NPDF;                  //! No. of PDFs involved
      int NSubProcessesLO;       //! No. of LO   subprocesses
      int NSubProcessesNLO;      //! No. of NLO  subprocesses
      int NSubProcessesNNLO;     //! No. of NNLO subprocesses
      int IPDFdef1;              //! Flag 1 to define PDF linear combinations of partonic subprocesses (e.g. hh --> jets: 3)
      int IPDFdef2;              //! Flag 2 to define PDF linear combinations (dep. on IPDFdef1; for 3 e.g. 1 for jet specific LCs, 121 for generic 11x11 matrix)
      int IPDFdef3LO;            //! Flag 3 to define PDF LCs at   LO (dep. on IPDFdef1, IPDFdef2; for 3, 1 e.g. 6 subprocesses, ignored for IPDFdef2==121)
      int IPDFdef3NLO;           //! Flag 3 to define PDF LCs at  NLO (dep. on IPDFdef1, IPDFdef2; for 3, 1 e.g. 7 subprocesses, ignored for IPDFdef2==121)
      int IPDFdef3NNLO;          //! Flag 3 to define PDF LCs at NNLO (dep. on IPDFdef1, IPDFdef2; for 3, 1 e.g. 7 subprocesses, ignored for IPDFdef2==121)
      int NPDFDim;               //! Define internal storage mode for PDF LCs (dep. on NPDF; e.g. for 1: 0 for linear, for 2: 1 for half- or 2 for full-matrix)
      std::vector<std::vector<std::pair<int,int> > > PDFCoeffLO;   //! PDF Linear combinations for   LO calculation (used only if IPDFdef2==0)
      std::vector<std::vector<std::pair<int,int> > > PDFCoeffNLO;  //! PDF Linear combinations for  NLO calculation (used only if IPDFdef2==0)
      std::vector<std::vector<std::pair<int,int> > > PDFCoeffNNLO; //! PDF Linear combinations for NNLO calculation (used only if IPDFdef2==0)
      std::vector<std::pair<int,int> > AsymmetricProcesses;        //! Specify processes that need to be exchanged in half-matrix notation, when xmin>xmax (only if NPDFDim==1)
      std::string Name;                      //!< More precise description for specific contribution (e.g. LO, pp -> 2 jets; also can add 'run-mode' and further details)
      std::vector<std::string> References;   //!< References for process (also other plain text lines can be included here)
      std::vector<std::string > GetProcessDescription() {
         //! Get 'ContrDescription' usable for fastNLO table
         unsigned int iadd = 0;
         if ( Name != "" ) iadd=1;
         std::vector<std::string > ProcDescr(References.size()+iadd);
         if ( iadd != 0 ) ProcDescr[0] = Name;
         for ( unsigned int i = 0 ; i<References.size() ; i++ )
            ProcDescr[i+iadd] = References[i];
         return ProcDescr;
      }
   };

};

#endif
