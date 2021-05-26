#ifndef __fnlotableconstants__
#define __fnlotableconstants__

// NEVER EVER include a project's internal config.h in installable header files!
// Use for conditional compilation only in .cc source code files.
// Otherwise conflicts with other linked projects are to be expected.
#include <iostream>
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
      // KR: Use C++11 possibility of non-static member initialization
      std::string Name{"Undefined"};         //!< Name and version of generator
      std::vector<std::string> References;   //!< References for generator. Potentially include additional information.
      int UnitsOfCoefficients {12};          //!< Prefactor of x section in barns for coefficients passed to fastNLO (neg. power of 10: pb->12, fb->15)
      //! Transform these constants into 'CodeDescription' usable with fastNLO table
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
      int LeadingOrder{-1};      //!< Order in alpha_s of leading order process
      int NPDF{-1};              //!< No. of PDFs involved
      int NSubProcessesLO{-1};   //!< No. of LO   subprocesses
      int NSubProcessesNLO{-1};  //!< No. of NLO  subprocesses
      int NSubProcessesNNLO{-1}; //!< No. of NNLO subprocesses
      int IPDFdef1{-1};          //!< Flag 1 to define PDF linear combinations of partonic subprocesses (e.g. hh --> jets: 3)
      int IPDFdef2{-1};          //!< Flag 2 to define PDF linear combinations (dep. on IPDFdef1; for 3 e.g. 1 for jet specific LCs, 121 for generic 11x11 matrix)
      int IPDFdef3LO{-1};        //!< Flag 3 to define PDF LCs at   LO (dep. on IPDFdef1, IPDFdef2; for 3, 1 e.g. 6 subprocesses, ignored for IPDFdef2==121)
      int IPDFdef3NLO{-1};       //!< Flag 3 to define PDF LCs at  NLO (dep. on IPDFdef1, IPDFdef2; for 3, 1 e.g. 7 subprocesses, ignored for IPDFdef2==121)
      int IPDFdef3NNLO{-1};      //!< Flag 3 to define PDF LCs at NNLO (dep. on IPDFdef1, IPDFdef2; for 3, 1 e.g. 7 subprocesses, ignored for IPDFdef2==121)
      int NPDFDim{-1};           //!< Define internal storage mode for PDF LCs (dep. on NPDF; e.g. for 1: 0 for linear, for 2: 1 for half- or 2 for full-matrix)
      std::vector<std::vector<std::pair<int,int> > > PDFCoeffLO;   //!< PDF Linear combinations for   LO calculation (used only if IPDFdef2==0)
      std::vector<std::vector<std::pair<int,int> > > PDFCoeffNLO;  //!< PDF Linear combinations for  NLO calculation (used only if IPDFdef2==0)
      std::vector<std::vector<std::pair<int,int> > > PDFCoeffNNLO; //!< PDF Linear combinations for NNLO calculation (used only if IPDFdef2==0)
      std::vector<std::vector<int> > PDFLiCoInLO;   //!< PDF Linear combinations for   LO calculation (used only if IPDFdef2==0) [definition as in steering] (used if PDFCoeffLO is empty)
      std::vector<std::vector<int> > PDFLiCoInNLO;  //!< PDF Linear combinations for  NLO calculation (used only if IPDFdef2==0) [definition as in steering]
      std::vector<std::vector<int> > PDFLiCoInNNLO; //!< PDF Linear combinations for NNLO calculation (used only if IPDFdef2==0) [definition as in steering]
      std::vector<std::pair<int,int> > AsymmetricProcesses; //!< Specify processes that need to be exchanged in half-matrix notation, when xmin>xmax (only if NPDFDim==1)
      std::string Name{"Undefined"};       //!<< More precise description for specific contribution (e.g. LO, pp -> 2 jets; also can add 'run-mode' and further details)
      std::vector<std::string> References; //!<< References for process (also other plain text lines can be included here)
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


   struct ScenarioConstants {
      //! fastNLO Scenario constants
      //! Steering parameters for a fastNLO scenario.
      //! Contains mostly the binning and description of the specific scenario

      std::string ScenarioName{"Undefined"};        //!< Name of the scenario
      std::vector<std::string> ScenarioDescription; //!< Description of the scenario
      int PublicationUnits{12}; //!< Unit of data cross sections (negative power of 10, e.g. 12->pb, 15->fb)
      int DifferentialDimension{0};             //!< Dimensionality of binning (1: single-differential, 2: double-differential; also decides if SingleDifferentialBinning or DoubleDifferentialBinning is used)
      std::vector<std::string> DimensionLabels; //!< Labels (symbol and unit) for the measurement dimensions (from outer to inner "loop"), e.g. "|y|" and "p_T [GeV]". This may also help to define the observables to be calculated in an automatized way!
      std::vector<int> DimensionIsDifferential; //!< Specify for each dimension whether:  0: the cross section is NOT differential,  i.e. there are two bin borders (but NO division (normalization) by bin width);  1 : the cross section is point-wise differential, i.e. only one point is given; 2 : the cross section is bin-wise differential,   i.e. there are two bin borders and division by bin width
      bool CalculateBinSize{true}; //!< Calculate bin width from lower and upper bin boundaries
      double BinSizeFactor{1};     //!< Possibility to provide additional normalization factor, e.g. of 2 for bins in |y|
      std::vector<double> BinSize; //!< If 'CalculateBinSize' is 'false' provide table with bin widths 'by hand' for normalization. If the calculation should not be divided by bin width, then use 'DimensionIsDifferential' equal '0', and set 'CalculateBinSize' 'true' for each dimension.
      std::string ScaleDescriptionScale1{"Undefined"}; //!< "<pT_1,2>_[GeV]" # This defines the scale to be used (Note: The 1st scale should always be in units of [GeV]!)
      std::string ScaleDescriptionScale2{"Undefined"}; //!< "pT_max_[GeV]"   # Specify 2nd scale name and unit (ONLY for flexible-scale tables)
      std::vector<double> SingleDifferentialBinning; //!< Observable binning Use either 'SingleDifferentialBinning' or 'DoubleDifferentialBinning' or 'TripleDifferentialBinning' in accordance with 'DifferentialDimension' above
      std::vector<std::vector<double> > DoubleDifferentialBinning; //!< Observable binning
      std::vector<std::vector<double> > TripleDifferentialBinning; //!< Observable binning
      double CenterOfMassEnergy{7000.}; //!< Center-of-mass energy in GeV. LHC Next Run II: 13000
      int PDF1{2212};                   //!< PDF of 1st hadron (following PDG convention: proton 2212).
      int PDF2{2212};                   //!< PDF of 2nd hadron (following PDG convention: proton 2212).
      std::string OutputFilename{"fastnlo"}; //!< Filename of fastNLO output table
      int OutputPrecision{8};                //!< Number of decimal digits to store in output table (def.=8).
      bool OutputCompression{true};          //!< If zlib available, gzip output table.
      // KR: Cache experimental, switch off by default
      int CacheType{0};    //!< Cache type: 1 or 2, 0 for deactivation
      int CacheMax{0};     //!< maximum size of cache
      int CacheCompare{0}; //!< number of elements to be compared with new entry
      // Flex-scale tables; should be future default ...
      bool FlexibleScaleTable{false};  //!< Create table fully flexible in mu_f (larger size, and requires scale independent weights during creation), true, or table with fixed number of mu_f scale factors, def.=false.
      int NFlexScales{2};              //!< No. of flexible scales to fill simultaneously; 1 or 2
      double FlexConstScale2{91.1876}; //!< Constant value set for 2nd "flexible" scale
      // Fixed-scale tables
      std::vector<double> ScaleVariationFactors; //!< Factorization scale variations (only needed for fixed-scale tables), List of scale factors must include factor '1', Scale factors will be ordered according to fastNLO convention: (1, min, ... , max). Defaults: {0.5, 1, 2}
      bool ReadBinningFromSteering{false};    //!< Specify if binning is read from fScenConst or from warmup
      bool IgnoreWarmupBinningCheck{false};   //!< Don't check warmup binning to avoid too many floating precision issues
      bool ApplyPDFReweighting{true};         //!<  Apply reweighting of pdfs for an optimized interpolation, def.=true.
      bool CheckScaleLimitsAgainstBins{true}; //!< For warmup-run! Set limits for scale nodes to bin borders, if possible
      // KR: Attention: Meaning of InclusiveJets keyword changed between nnlo-bridge code 0.0.40 and 0.0.46 to NNLOJET; see str files for info
      bool InclusiveJets{false}; //!< Flag to store one entry per jet, not just one per event; not stored in table; only used with NNLOJET so far
      // KR: Feature experimental, switch off by default
      double ReduceXmin{0.}; //!< Reduce xmin by n nodes (no change in number of x nodes)
      /**# -------------------------------------------------------------------- #
         #   Choose fastNLO interpolation kernels and distance measures
         # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         #   Currently implemented interpolation kernels
         #     Catmull
         #     Lagrange
         #     OneNode
         #     Linear
         #
         #   Currently implemented distance measures
         #     linear
         #     loglog025        eq. to (log(log(4*x)))
         #     log10
         #     sqrtlog10        eq. to sqrt(log_10(x))
         # -------------------------------------------------------------------- #
      */
      std::string X_Kernel{"Lagrange"};             //!< "Lagrange"
      std::string X_DistanceMeasure{"sqrtlog10"};   //!< "sqrtlog10"
      int X_NNodes{15};                             //!< 15
      std::string X_NNodeCounting{"NodesPerBin"};   //!< "NodesPerBin" ("NodesMax", "NodesPerMagnitude")

      std::string Mu1_Kernel{"Lagrange"};           //!< "Lagrange"
      std::string Mu1_DistanceMeasure{"loglog025"}; //!< "loglog025"
      int Mu1_NNodes{6};                            //!< 6

      std::string Mu2_Kernel{"Lagrange"};           //!< "Lagrange"; Scale2 not used for fixed-scale tables
      std::string Mu2_DistanceMeasure{"loglog025"}; //!< "loglog025"
      int Mu2_NNodes{6};                            //!< 6

      // KR: Defaults can be set using SetScenConstsDefaults().
      //     If the struct is directly created within another project,
      //     uninitialised struct members are possible!
   };


   struct WarmupConstants {
      //! Variables from warmup-run
      //! Initalize WarmupConstants with ScenarioConstants
      //! for consistency.
      //! Furthermore needed for full initialization:
      //!   + OrderInAlphasOfWarmupRunWas
      //!   + Binning
      //!   + Values
      //!   + headerValues
      int OrderInAlphasOfWarmupRunWas;
      bool CheckScaleLimitsAgainstBins;
      std::string ScaleDescriptionScale1;
      std::string ScaleDescriptionScale2;
      int DifferentialDimension;
      std::vector<std::string> DimensionLabels;
      std::vector<int> DimensionIsDifferential;
      std::vector<std::vector<double> > Values;
      std::vector<std::string> headerValues;
      std::vector<std::vector<double> > Binning;
   public:
      WarmupConstants(const ScenarioConstants& scenario) {
         Init();
         ScaleDescriptionScale1 = scenario.ScaleDescriptionScale1;
         ScaleDescriptionScale2 = scenario.ScaleDescriptionScale2;
         CheckScaleLimitsAgainstBins = scenario.CheckScaleLimitsAgainstBins;
         DifferentialDimension = scenario.DifferentialDimension;
         DimensionLabels = scenario.DimensionLabels;
         DimensionIsDifferential = scenario.DimensionIsDifferential;
         std::cout<<"Warning [WarmupConstants]. Binning has not be taken over from ScenarioConstants (not implemented.)"<<std::endl;
         //exit(4);
         // Binning = ;
         // Values = ;
         // headerValues = ;
      }
      void Init(){
         OrderInAlphasOfWarmupRunWas=-1;
         //    ObsBin      x_min      x_max    173.3GeV_min    173.3GeV_max           y_min           y_max
         headerValues.clear();
         // headerValues.push_back("ObsBin");
         // headerValues.push_back("x_min");
         // headerValues.push_back("x_max");
      };

      WarmupConstants() {
         Init();
      }
   };
};

#endif
