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
      //!  - name of generator
      //!  - references for generator
      //!  - units of coefficients (typically pb -> 12)
      //!  - (additional information about generator may 
      //!    be included in References)
      int UnitsOfCoefficients; //!< Units of coeffients as passed to fastNLO (negative power of 10: pb->12, fb->15
      std::string Name; //!< Name of generator
      std::vector<std::string> References; //!< References for generator. Include additional information here (e.g. 'run-mode' or process)
      std::vector<string > GetCodeDescription() {
	 //! Get 'CodeDescription' usable for fastNLO table
	 std::vector<string > CodeDescr(References.size()+1);
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
      //! Please see fastNLO table format definition for 
      //! a detailed explanation.
      //! 
      //! AsymmetricProcesses may only be required if NPDFDim=1 is chosen.
      int NPDF; //!< Number of PDFs involved
      int NSubProcessesLO; //!< Number of subprocesses of the considered process in LO run.
      int NSubProcessesNLO; //!< Number of subprocesses of the considered process in NLO run.
      int NSubProcessesNNLO; //!< Number of subprocesses of the considered process in NNLO run.
      int IPDFdef1; //!< Define PDF linear combinations corresponding to partonic subprocesses (hadron-hadron: 3
      int IPDFdef2; //!< Flag to define PDF linear combinations (dependent on IPDFdef1. Use 1 for jet-production in pp/ppbar)
      int IPDFdef3LO; //!< Unique identifier (dependent on IPDFdef1, IPDFdef2) for specifying the PDF linear combinations.
      int IPDFdef3NLO; //!< Unique identifier (dependent on IPDFdef1, IPDFdef2) for specifying the PDF linear combinations.
      int IPDFdef3NNLO; //!< Unique identifier (dependent on IPDFdef1, IPDFdef2) for specifying the PDF linear combinations. 
      int NPDFDim; //!< Internal way to store PDF linear combinations. Use 1 (half-matrix storage) or 2 (full-matrix storage) for hadron-hadron collisions.
      std::vector<std::pair<int,int> > AsymmetricProcesses; //!< (if NPDFDim=1) Specify processes that need to be exchanged in half-matrix notation, when xmin>xmax
   };
      
};

#endif
