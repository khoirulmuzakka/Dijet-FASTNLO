// Author: Daniel Britzger
// DESY, 29/07/2013
// ___________________________________________________________________________________________________
//! The fastNLOCreate class
/*!
  This class can generate/fill one single contribution for a fastNLO table.
  It supports fixed scale and flexible-scale tables.
  fastNLOCreate inherits from fastNLOTable, but is only able to hold one
  coefficient table as member.
  fastNLOCreate no enables to fill this coefficient table, i.e.  to add further
  contributions, e.g. from a MC generator.

  fastNLOCreate works only with a steering file. Example steering files
  are provided together with this class.
  The steering specifies, which kind of process are stored and also
  the binning is read in.

  fastNLOCreate also handles the warmup runs. If no warmup table is found
  it automatically runs in warmup mode. If warmup values are available, it
  runs in 'production' mode.

  In order ot obtain a full fastNLO table, i.e. a table with LO and NLO contributions,
  several contributions have to be merged, using the fnlo-merge (or fnlo-tk-merge)
  tools.

  For example applications, please contact the authors.

*/
// ___________________________________________________________________________________________________

#include <algorithm>
#include <cfloat>
#include <string>
#include <unistd.h>
#include "fastnlotk/fastNLOCreate.h"
#include "fastnlotk/fastNLOTools.h"
#include "fastnlotk/fastNLOCoeffAddFlex.h"
#include "fastnlotk/fastNLOCoeffAddFix.h"
#include "fastnlotk/fastNLOInterpolCatmullRom.h"
#include "fastnlotk/fastNLOInterpolLagrange.h"
#include "fastnlotk/fastNLOInterpolLinear.h"
#include "fastnlotk/fastNLOInterpolOneNode.h"
#include "fastnlotk/read_steer.h"
#include <unistd.h>


using namespace std;

// ___________________________________________________________________________________________________
fastNLOCreate::fastNLOCreate() {
   logger.SetClassName("fastNLOCreate");
}


// ___________________________________________________________________________________________________
fastNLOCreate::fastNLOCreate(const string& steerfile, const fastNLO::GeneratorConstants& GenConsts, const fastNLO::ProcessConstants& ProcConsts) {
   //! Constructor of fastNLOCreate
   //!
   //! Pass steering parameters through
   //! GeneratorConstants, ProcessConstants and the steering file
   //! (see GeneratorConstants.h file for details)
   //!
   //!

   logger.SetClassName("fastNLOCreate");
   logger.debug["fastNLOCreate"]<<"Create table from GenConsts, ProcConsts, and steering file"<<endl;

   //! Set constants from arguments
   fGenConsts  = GenConsts;
   fProcConsts = ProcConsts;

   //! Steering file settings take precedence over settings in code
   //! The WarmupFilename is set from steering
   ReadSteering(steerfile);
   ReadGenAndProcConstsFromSteering();

   bool check = CheckProcConsts();
   if (!check) {
      logger.warn["fastNLOCreate"]<<"Process constants not properly initialized! Please check your steering."<<endl;
   }

   logger.debug["fastNLOCreate"]<<"Instantiate from GenConsts, ProcConsts, and steering file"<<endl;
   Instantiate();
}


// ___________________________________________________________________________________________________
fastNLOCreate::fastNLOCreate(const string& warmupfile, const fastNLO::GeneratorConstants& GenConsts,
                             const fastNLO::ProcessConstants& ProcConsts, const fastNLO::ScenarioConstants& ScenConsts) {
   //! Constructor of fastNLOCreate
   //!
   //! Pass all needed steering paramters through
   //! GeneratorConstants, ProcessConstants and ScenarioConstants
   //! (see GeneratorConstants.h file for details)
   //!
   //! warmupfile: filename to be written out or read in for production run, if already existent
   //!
   //! No steering file is read in
   //!

   logger.SetClassName("fastNLOCreate");
   logger.debug["fastNLOCreate"]<<"Create table from GenConsts, ProcConsts, ScenConsts, and warmup file"<<endl;
   logger.debug["fastNLOCreate"]<<"The warmup filename set from function call is: " << warmupfile << endl;

   //! Set constants from arguments
   fGenConsts  = GenConsts;
   fScenConsts = ScenConsts;
   fProcConsts = ProcConsts;
   if ( !CheckProcConsts() ) {
      logger.warn["fastNLOCreate"]<<"Process constants not properly initialized! Please check the ProcessConstants."<<endl;
      exit(1);
   }

   //! The WarmupFilename is set from argument
   fSteerfile = warmupfile;
   fWarmupFilename = warmupfile;

   // --- check and transform parton combinations
   if ( fProcConsts.IPDFdef2 == 0 ) {
      if ( fProcConsts.IPDFdef3LO > 0 && fProcConsts.PDFCoeffLO.empty() )
         fProcConsts.PDFCoeffLO   = ReadPartonCombinations(0,fProcConsts.PDFLiCoInLO);
      if ( fProcConsts.IPDFdef3NLO > 0 && fProcConsts.PDFCoeffNLO.empty()  )
         fProcConsts.PDFCoeffNLO  = ReadPartonCombinations(1,fProcConsts.PDFLiCoInNLO);
      if ( fProcConsts.IPDFdef3NNLO > 0 && fProcConsts.PDFCoeffNNLO.empty() )
         fProcConsts.PDFCoeffNNLO = ReadPartonCombinations(2,fProcConsts.PDFLiCoInNNLO);
   }

   logger.debug["fastNLOCreate"]<<"Instantiate table from GenConsts, ProcConsts, ScenConsts, and warmup file"<<endl;
   Instantiate();
}


// ___________________________________________________________________________________________________
fastNLOCreate::fastNLOCreate(const fastNLO::GeneratorConstants& GenConsts, const fastNLO::ProcessConstants& ProcConsts,
                             const fastNLO::ScenarioConstants& ScenConsts, const fastNLO::WarmupConstants& WarmupConsts ) {
   //! Constructor of fastNLOCreate
   //!
   //! Pass all needed steering paramters through
   //! GeneratorConstants, ProcessConstants, ScenarioConstants and WarmupConstants
   //! (see GeneratorConstants.h file for details)
   //!
   //! No steering or warmup file is read in
   //!

   logger.SetClassName("fastNLOCreate");
   logger.debug["fastNLOCreate"]<<"Create table from GenConsts, ProcConsts, ScenConsts, and WarmupConsts"<<endl;

   //! Set constants from arguments
   fGenConsts  = GenConsts;
   fScenConsts = ScenConsts;
   fProcConsts = ProcConsts;
   fWarmupConsts = WarmupConsts;
   if ( !CheckProcConsts() ) {
      logger.warn["fastNLOCreate"]<<"Process constants not properly initialized! Please check the ProcessConstants."<<endl;
      exit(1);
   }

   //! No WarmupFile required, a pseudo-WarmupFilename is defined here
   fSteerfile = "NoSteeringFileMode";//warmupfile;
   fWarmupFilename = fSteerfile;
   logger.debug["fastNLOCreate"]<<"Warmup set from code; the pseudo(!)-warmup filename is: " << fWarmupFilename << endl;

   // --- check and transform parton combinations
   if ( fProcConsts.IPDFdef2 == 0 ) {
      if ( fProcConsts.IPDFdef3LO > 0 && fProcConsts.PDFCoeffLO.empty() )
         fProcConsts.PDFCoeffLO   = ReadPartonCombinations(0,fProcConsts.PDFLiCoInLO);
      if ( fProcConsts.IPDFdef3NLO > 0 && fProcConsts.PDFCoeffNLO.empty()  )
         fProcConsts.PDFCoeffNLO  = ReadPartonCombinations(1,fProcConsts.PDFLiCoInNLO);
      if ( fProcConsts.IPDFdef3NNLO > 0 && fProcConsts.PDFCoeffNNLO.empty() )
         fProcConsts.PDFCoeffNNLO = ReadPartonCombinations(2,fProcConsts.PDFLiCoInNNLO);
   }

   logger.debug["fastNLOCreate"]<<"Instantiate table from GenConsts, ProcConsts, ScenConsts, and WarmupConsts"<<endl;
   Instantiate();
}


// ___________________________________________________________________________________________________
fastNLOCreate::fastNLOCreate(const string& steerfile, string steeringNameSpace, bool shouldReadSteeringFile) {
   //!
   //! Constructor for fastNLOCreate
   //!
   //! steerfile:                Steering file containing all required flags
   //! steeringNameSpace:        Alternative read_steer namespace to be used (optional)
   //! shouldReadSteeringFile:   Force to NOT read the steering file (e.g. the file has already been read in)  (optional).
   //!

   logger.SetClassName("fastNLOCreate");
   logger.debug["fastNLOCreate"]<<"Create table from steering file"<<endl;

   //! Set constants from defaults
   SetGenConstsDefaults();
   SetProcConstsDefaults();

   //! Steering file settings take precedence over defaults
   //! The WarmupFilename is set from steering
   ReadSteering(steerfile, steeringNameSpace, shouldReadSteeringFile);
   ReadGenAndProcConstsFromSteering();

   bool check = CheckProcConsts();
   if (!check) {
      logger.warn["fastNLOCreate"]<<"Process constants not properly initialized! Please check your steering."<<endl;
   }

   logger.debug["fastNLOCreate"]<<"Instantiate table from steering file"<<endl;
   Instantiate();
}


// ___________________________________________________________________________________________________
void fastNLOCreate::SetGenConstsDefaults() {
   //! Set default values for generator constants
   fGenConsts.Name = "Theory";
   fGenConsts.UnitsOfCoefficients = 12;   //!< X section units of coefficients passed to fastNLO (neg. power of 10: pb->12, fb->15)
}


// ___________________________________________________________________________________________________
void fastNLOCreate::SetProcConstsDefaults() {
   //! Set default values for process constants
   fProcConsts.LeadingOrder        = -1;   //!< Order in alpha_s of leading order process
   fProcConsts.NPDF                = -1;   //!< No. of PDFs involved
   fProcConsts.NSubProcessesLO     = -1;   //!< No. of LO   subprocesses
   fProcConsts.NSubProcessesNLO    = -1;   //!< No. of NLO  subprocesses
   fProcConsts.NSubProcessesNNLO   = -1;   //!< No. of NNLO subprocesses
   fProcConsts.IPDFdef1            = -1;   //!< Flag 1 to define PDF linear combinations of partonic subprocesses (e.g. hh --> jets: 3)
   fProcConsts.IPDFdef2            = -1;   //!< Flag 2 to define PDF linear combinations (dep. on IPDFdef1; for 3 e.g. 1 for jet specific LCs, 121 for generic 11x11 matrix)
   fProcConsts.IPDFdef3LO          = -1;   //!< Flag 3 to define PDF LCs at   LO (dep. on IPDFdef1, IPDFdef2; for 3, 1 e.g. 6 subprocesses, ignored for IPDFdef2==121)
   fProcConsts.IPDFdef3NLO         = -1;   //!< Flag 3 to define PDF LCs at  NLO (dep. on IPDFdef1, IPDFdef2; for 3, 1 e.g. 7 subprocesses, ignored for IPDFdef2==121)
   fProcConsts.IPDFdef3NNLO        = -1;   //!< Flag 3 to define PDF LCs at NNLO (dep. on IPDFdef1, IPDFdef2; for 3, 1 e.g. 7 subprocesses, ignored for IPDFdef2==121)
   fProcConsts.NPDFDim             = -1;   //!< Define internal storage mode for PDF LCs (dep. on NPDF; e.g. for 1: 0 for linear, for 2: 1 for half- or 2 for full-matrix)
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::CheckProcConsts() {
   //! Check that reasonable values different from the defaults have been set
   logger.debug["CheckProcConsts"]<<"Checking ProcConsts: "<<endl;

   bool checkok = true;
   if (fProcConsts.LeadingOrder < 0) {
      logger.warn["CheckProcConsts"]<<"Order in alpha_s of leading order process not properly set: "<<fProcConsts.LeadingOrder<<endl;
      checkok = false;
   }
   if (fGenConsts.UnitsOfCoefficients < 0) {
      logger.warn["CheckProcConsts"]<<"Power of X section units of coefficients not properly set: "<<fGenConsts.UnitsOfCoefficients<<endl;
      checkok = false;
   }
   if (fProcConsts.NPDF < 1) {
      logger.warn["CheckProcConsts"]<<"No. of PDFs not properly set: "<<fProcConsts.NPDF<<endl;
      checkok = false;
   }
   if (fProcConsts.NSubProcessesLO < 1) {
      logger.warn["CheckProcConsts"]<<"No. of LO subprocesses not properly set: "<<fProcConsts.NSubProcessesLO<<endl;
      checkok = false;
   }
   if (fProcConsts.NSubProcessesNLO < 1) {
      logger.warn["CheckProcConsts"]<<"No. of NLO subprocesses not properly set: "<<fProcConsts.NSubProcessesNLO<<endl;
      checkok = false;
   }
   if (fProcConsts.NSubProcessesNNLO < 1) {
      logger.warn["CheckProcConsts"]<<"No. of NNLO subprocesses not properly set: "<<fProcConsts.NSubProcessesNNLO<<endl;
      checkok = false;
   }
   if (fProcConsts.IPDFdef1 < 0) {
      logger.warn["CheckProcConsts"]<<"Flag 1 to define PDF linear combination not properly set: "<<fProcConsts.IPDFdef1<<endl;
      checkok = false;
   }
   if (fProcConsts.IPDFdef2 < 0) {
      logger.warn["CheckProcConsts"]<<"Flag 2 to define PDF linear combination not properly set: "<<fProcConsts.IPDFdef2<<endl;
      checkok = false;
   }
   if (fProcConsts.IPDFdef3LO < 0) {
      logger.warn["CheckProcConsts"]<<"Flag 3 LO to define PDF linear combination not properly set: "<<fProcConsts.IPDFdef3LO<<endl;
      checkok = false;
   }
   if (fProcConsts.IPDFdef3NLO < 0) {
      logger.warn["CheckProcConsts"]<<"Flag 3 NLO to define PDF linear combination not properly set: "<<fProcConsts.IPDFdef3NLO<<endl;
      checkok = false;
   }
   if (fProcConsts.IPDFdef3NNLO < 0) {
      logger.warn["CheckProcConsts"]<<"Flag 3 NNLO to define PDF linear combination not properly set: "<<fProcConsts.IPDFdef3NNLO<<endl;
      checkok = false;
   }
   if (fProcConsts.NPDFDim < 0) {
      logger.warn["CheckProcConsts"]<<"Internal storage mode for PDF LCs not properly set: "<<fProcConsts.NPDFDim<<endl;
      checkok = false;
   }
   //    std::vector<std::vector<std::pair<int,int> > > PDFCoeffLO; //! PDF Linear combinations for LO calculation (used only if IPDFdef2==0)
   //    std::vector<std::vector<std::pair<int,int> > > PDFCoeffNLO; //! PDF Linear combinations for NLO calculation (used only if IPDFdef2==0)
   //    std::vector<std::vector<std::pair<int,int> > > PDFCoeffNNLO; //! PDF Linear combinations for NNLO calculation (used only if IPDFdef2==0)
   //    std::vector<std::pair<int,int> > AsymmetricProcesses; //!< (if NPDFDim=1) Specify processes that need to be exchanged in half-matrix notation, when xmin>xmax
    //         std::vector<std::string > ProcDescr(References.size()+iadd);

   return checkok;
}


// ___________________________________________________________________________________________________
void fastNLOCreate::ReadGenAndProcConstsFromSteering() {
   //! ReadGenAndProcConstsFromSteering()
   //! If generator and process constants have not been set in
   //! the constructor, then obtain these values from the steering
   //! file.

   logger.debug["ReadGenAndProcConstsFromSteering"]<<endl;

   // Generator constants
   if (EXIST_NS(CodeDescription,fSteerfile)) {
      vector<string > CodeDescr = STRING_ARR_NS(CodeDescription,fSteerfile);
      fGenConsts.Name = CodeDescr[0];
      if (CodeDescr.size() > 1) {
         fGenConsts.References.resize(CodeDescr.size()-1);
         for (unsigned int i = 0 ; i< fGenConsts.References.size() ; i++)
            fGenConsts.References [i] = CodeDescr[i+1];
      }
   }
   if (EXIST_NS(UnitsOfCoefficients,fSteerfile)) fGenConsts.UnitsOfCoefficients = INT_NS(UnitsOfCoefficients,fSteerfile);

   // Process constants
   if (EXIST_NS(LeadingOrder,fSteerfile))        fProcConsts.LeadingOrder = INT_NS(LeadingOrder,fSteerfile);
   if (EXIST_NS(NPDF,fSteerfile))                fProcConsts.NPDF = INT_NS(NPDF,fSteerfile);
   if (EXIST_NS(NSubProcessesLO,fSteerfile))     fProcConsts.NSubProcessesLO = INT_NS(NSubProcessesLO,fSteerfile);
   if (EXIST_NS(NSubProcessesNLO,fSteerfile))    fProcConsts.NSubProcessesNLO = INT_NS(NSubProcessesNLO,fSteerfile);
   if (EXIST_NS(NSubProcessesNNLO,fSteerfile))   fProcConsts.NSubProcessesNNLO = INT_NS(NSubProcessesNNLO,fSteerfile);
   if (EXIST_NS(IPDFdef1,fSteerfile))            fProcConsts.IPDFdef1 = INT_NS(IPDFdef1,fSteerfile);
   if (EXIST_NS(IPDFdef2,fSteerfile))            fProcConsts.IPDFdef2 = INT_NS(IPDFdef2,fSteerfile);
   if (EXIST_NS(IPDFdef3LO,fSteerfile))          fProcConsts.IPDFdef3LO = INT_NS(IPDFdef3LO,fSteerfile);
   if (EXIST_NS(IPDFdef3NLO,fSteerfile))         fProcConsts.IPDFdef3NLO = INT_NS(IPDFdef3NLO,fSteerfile);
   if (EXIST_NS(IPDFdef3NNLO,fSteerfile))        fProcConsts.IPDFdef3NNLO = INT_NS(IPDFdef3NNLO,fSteerfile);

   if ( fProcConsts.IPDFdef2 == 0 ) {
      if ( fProcConsts.IPDFdef3LO > 0 && fProcConsts.IPDFdef3LO != (int)fProcConsts.PDFCoeffLO.size()  )
         fProcConsts.PDFCoeffLO   = ReadPartonCombinations(0,INT_TAB_NS(PartonCombinationsLO,fSteerfile));
      if ( fProcConsts.IPDFdef3NLO > 0 && fProcConsts.IPDFdef3NLO != (int)fProcConsts.PDFCoeffNLO.size() )
         fProcConsts.PDFCoeffNLO  = ReadPartonCombinations(1,INT_TAB_NS(PartonCombinationsNLO,fSteerfile));
      if ( fProcConsts.IPDFdef3NNLO > 0 && fProcConsts.IPDFdef3NNLO != (int)fProcConsts.PDFCoeffNNLO.size() )
         fProcConsts.PDFCoeffNNLO = ReadPartonCombinations(2,INT_TAB_NS(PartonCombinationsNNLO,fSteerfile));
      // check again
      // --- check and transform parton combinations
      if ( fProcConsts.IPDFdef3LO > 0 && fProcConsts.PDFCoeffLO.empty() )
         fProcConsts.PDFCoeffLO   = ReadPartonCombinations(0,fProcConsts.PDFLiCoInLO);
      if ( fProcConsts.IPDFdef3NLO > 0 && fProcConsts.PDFCoeffNLO.empty()  )
         fProcConsts.PDFCoeffNLO  = ReadPartonCombinations(1,fProcConsts.PDFLiCoInNLO);
      if ( fProcConsts.IPDFdef3NNLO > 0 && fProcConsts.PDFCoeffNNLO.empty() )
         fProcConsts.PDFCoeffNNLO = ReadPartonCombinations(2,fProcConsts.PDFLiCoInNNLO);
   }

   if (EXIST_NS(NPDFDim,fSteerfile))             fProcConsts.NPDFDim = INT_NS(NPDFDim,fSteerfile);

   // read asymmetric processes if half-matrix notation is requested
   if ( fProcConsts.NPDFDim == 1 ) {
      if ( fProcConsts.NPDF == 2 && fProcConsts.IPDFdef1 == 3 && ( fProcConsts.IPDFdef2 == 121 || fProcConsts.IPDFdef2 == 169 ) ) {
         const int np = fProcConsts.IPDFdef2==121 ? 11:13;
         int p1 = 0;
         int p2 = 0;
         for ( int p = 0 ; p<fProcConsts.IPDFdef2 ; p++ ) {
            int pid = p1*(np)+p2;
            int asympid = p2*(np)+p1;
            if ( pid != asympid )  // actually not needed necessarily
               fProcConsts.AsymmetricProcesses.push_back(make_pair(pid,asympid));
            p2++;
            if ( p2 == np ) {
               p2=0;
               p1++;
            }
         }
      }
      else if ( fProcConsts.NPDF == 2 ) {
         if (EXIST_NS(AsymmetricProcesses,fSteerfile)) {
            fProcConsts.AsymmetricProcesses.clear();
            vector<vector<int> > asym = INT_TAB_NS(AsymmetricProcesses,fSteerfile);
            for (unsigned int i = 0 ; i<asym.size() ; i++) {
               if ( asym[i].size() != 2 ) {
                  logger.error["ReadGenAndProcConstsFromSteering"]<<"Asymmetric process "<<asym[i][0]<<", must have exactly one counter process."<<endl;
                  exit(1);
               }
               fProcConsts.AsymmetricProcesses.push_back(make_pair(asym[i][0],asym[i][1]));
            }
         }
      }
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::Instantiate() {
   //! Instantiate all internal members
   //! and prepare for filling
   logger.debug["Instantiate"]<<"Instantiate all internal members and prepare for filling " << endl;

   // init member variables
   fReader = NULL;

   fCacheMax = 20; // currently not working!
   fWarmupXMargin = 4; // was 4
   fWarmupNDigitMu1 = 1; //1 by purpose
   fWarmupNDigitMu2 = 2; //2 by purpose

   // Try to get warm-up values.
   // Otherwise a warm-up run will be initialized.
   GetWarmupValues();

   // --- ScenConsts sanity test
   // ...

   ILOord = fProcConsts.LeadingOrder;
   fIOrd = ILOord; // initialize with LO

   // -------------------------
   // header
   SetScenName(fScenConsts.ScenarioName);
   SetItabversion(fastNLO::tabversion);

   // ---- scenario specific flags
   Ipublunits   = fScenConsts.PublicationUnits;
   ScDescript   = fScenConsts.ScenarioDescription;

   Ecms         = fScenConsts.CenterOfMassEnergy;   // is often superseeded by generator-specific code.
   INormFlag    = 0;

   string filename = fScenConsts.OutputFilename;
   if ( filename.find(".gz") != string::npos ) fScenConsts.OutputCompression = true;
   else if ( fScenConsts.OutputCompression ) filename += ".gz";
   SetFilename(filename);

   fIsFlexibleScale  = fScenConsts.FlexibleScaleTable;
   fApplyPDFReweight = fScenConsts.ApplyPDFReweighting;
   SetOutputPrecision(fScenConsts.OutputPrecision);

   //if ( !fIsFlexibleScale )  ReadScaleFactors() ; // is called only when setting order of calculation


   // ---- init bin grid
   if ( fIsWarmup )  ReadBinningFromScenarioConsts();
   else if ( fScenConsts.ReadBinningFromSteering ) {
      ReadBinningFromScenarioConsts();
      CheckWarmupConsistency();
   }
   else {
      UseBinGridFromWarmup();
   }

   // now create one coefficient tasble.
   InitCoeffTable();
   // no info output for the following calls
   bool vol = logger.info.GetSpeak();
   logger.info.DoSpeak(false);
   SetOrderOfAlphasOfCalculation(fIOrd);
   logger.info.DoSpeak(vol);//reset verbosity level

   // Init interpolation kernels
   if (!fIsWarmup) {
      InitInterpolationKernels();
      InitGrids();
   }
}


// ___________________________________________________________________________________________________
fastNLOCreate::~fastNLOCreate() {
   // todo. cleanup arrays of kernels.
}


// ___________________________________________________________________________________________________
void fastNLOCreate::ReadSteering(string steerfile, string steeringNameSpace, bool shouldReadSteeringFile) {
   //! read in steering file
   //! The filename of the steering file
   //! is used as the 'namespace' of keys in read_steer
   //! if there is no steeringNameSpace given explicitly

   //! Remove extension from steerfile to define default steering namespace
   string steerbase = steerfile.substr(0, steerfile.find_last_of("."));

   logger.debug["ReadSteering"]<<"Steerfile = "<<steerfile<<endl;
   if (steeringNameSpace.empty()) {
      steeringNameSpace = steerbase;
   }
   fSteerfile =  steeringNameSpace;
   if ( shouldReadSteeringFile ) {
      READ_NS(steerfile,fSteerfile);
   }

   //! Set verbosity from steering or to default WARNING
   if (EXIST_NS(GlobalVerbosity,fSteerfile) )
      SetGlobalVerbosity(STRING_NS(GlobalVerbosity,fSteerfile));
   else SetGlobalVerbosity("WARNING");

   //! Set WarmupFilename from steering
   if (EXIST_NS(WarmUpFilename,fSteerfile)) {
      SetWarmupTableFilename(STRING_NS(WarmUpFilename,fSteerfile));
      logger.debug["fastNLOCreate"]<<"The warmup filename set in steering file is: " << STRING_NS(WarmUpFilename,fSteerfile) << endl;
   } else {
      //      string fWarmUpFile = STRING_NS(ScenarioName,fSteerfile) + "_" + steeringNameSpace + "_warmup.txt";
      // Default for NLOJet++ is like this; to be unified across generators and create methods
      // TODO: KR: Consider keeping all warmup filenames to .wrm
      //      string fWarmUpFile = steeringNameSpace + "_" + STRING_NS(ScenarioName,fSteerfile) + "_warmup.txt";
      string fWarmUpFile = steeringNameSpace + "_" + STRING_NS(ScenarioName,fSteerfile) + ".wrm";
      SetWarmupTableFilename(fWarmUpFile);
      logger.debug["fastNLOCreate"]<<"The warmup filename derived from steering is: " << fWarmUpFile << endl;
   }


   if (logger.info.GetSpeak())
      PRINTALL();

   // ---- set/read scenario constants (fScenConsts)
   fScenConsts.PublicationUnits = INT_NS(PublicationUnits,fSteerfile);
   fScenConsts.ScenarioName = STRING_NS(ScenarioName,fSteerfile);
   fScenConsts.ScenarioDescription = STRING_ARR_NS(ScenarioDescription,fSteerfile);
   fScenConsts.CenterOfMassEnergy = DOUBLE_NS(CenterOfMassEnergy,fSteerfile);
   fScenConsts.OutputFilename = STRING_NS(OutputFilename,fSteerfile);
   fScenConsts.OutputCompression = BOOL_NS(OutputCompression,fSteerfile);
   fScenConsts.FlexibleScaleTable = BOOL_NS(FlexibleScaleTable,fSteerfile);
   fScenConsts.ApplyPDFReweighting =  BOOL_NS(ApplyPDFReweighting,fSteerfile);
   fScenConsts.OutputPrecision = INT_NS(OutputPrecision,fSteerfile);
   fScenConsts.ReadBinningFromSteering = BOOL_NS(ReadBinningFromSteering,fSteerfile);

   // --- (pre-)set flexible scale flag
   fIsFlexibleScale  = fScenConsts.FlexibleScaleTable;

   // ---- set/read further scenario constants (fScenConsts)
   if ( fScenConsts.ReadBinningFromSteering ) {
      // read values from 'main steering'
      fScenConsts.DifferentialDimension = INT_NS(DifferentialDimension,fSteerfile);
      NDim = fScenConsts.DifferentialDimension;
      fScenConsts.DimensionIsDifferential = INT_ARR_NS(DimensionIsDifferential,fSteerfile);;
      fScenConsts.DimensionLabels = STRING_ARR_NS(DimensionLabels,fSteerfile);
      fScenConsts.ScaleDescriptionScale1 = STRING_NS(ScaleDescriptionScale1,fSteerfile);
      if (!fIsFlexibleScale)
         fScenConsts.ScaleDescriptionScale2 = STRING_NS(ScaleDescriptionScale2,fSteerfile);

      if ( NDim==1 ) fScenConsts.SingleDifferentialBinning = DOUBLE_ARR_NS(SingleDifferentialBinning,fSteerfile);
      else if ( NDim==2 ) fScenConsts.DoubleDifferentialBinning = DOUBLE_TAB_NS(DoubleDifferentialBinning,fSteerfile);
      else if ( NDim==3 ) fScenConsts.TripleDifferentialBinning = DOUBLE_TAB_NS(TripleDifferentialBinning,fSteerfile);

      fScenConsts.CalculateBinSize = BOOL_NS(CalculateBinSize,fSteerfile);
      fScenConsts.BinSizeFactor = DOUBLE_NS(BinSizeFactor,fSteerfile);
      if ( !fScenConsts.CalculateBinSize )
         fScenConsts.BinSize = DOUBLE_ARR_NS(BinSize,fSteerfile);

      //ReadBinSize();
   }
   else {
      // warmup file not yet read in
      // read in needed flags in GetWarmupValues()

      // fScenConsts.DifferentialDimension = INT_NS(Warmup.DifferentialDimension,fSteerfile);
      // fScenConsts.DimensionIsDifferential = INT_ARR_NS(Warmup.DimensionIsDifferential,fSteerfile);
      // fScenConsts.DimensionLabels = STRING_ARR_NS(Warmup.DimensionLabels,fSteerfile);
      // fScenConsts.ScaleDescriptionScale1 = STRING_NS(Warmup.ScaleDescriptionScale1,fSteerfile);
      // if (!fIsFlexibleScale)
      //         fScenConsts.ScaleDescriptionScale2 = STRING_NS(Warmup.ScaleDescriptionScale2,fSteerfile);

      fScenConsts.CalculateBinSize = false;

      // NDim = fScenConsts.DifferentialDimension;
      // if ( NDim==1 ) fScenConsts.SingleDifferentialBinning = ;
      // else if ( NDim==2 ) fScenConsts.DoubleDifferentialBinning = ;
      // else if ( NDim==3 ) fScenConsts.TripleDifferentialBinning = ;
   }

   fScenConsts.CheckScaleLimitsAgainstBins = BOOL_NS(CheckScaleLimitsAgainstBins,fSteerfile);
   fScenConsts.X_Kernel = STRING_NS(X_Kernel,fSteerfile);
   fScenConsts.X_DistanceMeasure = STRING_NS(X_DistanceMeasure,fSteerfile);
   fScenConsts.X_NNodes = INT_NS(X_NNodes,fSteerfile);
   fScenConsts.X_NNodeCounting="NodesPerBin";
   if ( EXIST_NS(X_NodeCounting,fSteerfile) )
        fScenConsts.X_NNodeCounting = STRING_NS(X_NodeCounting,fSteerfile);
   if ( EXIST_NS(X_NoOfNodesPerMagnitude,fSteerfile) )
        logger.warn[""]<<"Key 'X_NoOfNodesPerMagnitude' found. This is no longer supported and ignored. Please use steering key 'X_NNodeCounting' instead."<<endl;
        //fScenConsts.X_NoOfNodesPerMagnitude = BOOL_NS(X_NoOfNodesPerMagnitude,fSteerfile);

   fScenConsts.Mu1_Kernel = STRING_NS(Mu1_Kernel,fSteerfile);
   fScenConsts.Mu1_DistanceMeasure = STRING_NS(Mu1_DistanceMeasure,fSteerfile);
   fScenConsts.Mu1_NNodes = INT_NS(Mu1_NNodes,fSteerfile);
   if ( fIsFlexibleScale ) {
      fScenConsts.Mu2_Kernel = STRING_NS(Mu2_Kernel,fSteerfile);
      fScenConsts.Mu2_DistanceMeasure = STRING_NS(Mu2_DistanceMeasure,fSteerfile);
      fScenConsts.Mu2_NNodes = INT_NS(Mu2_NNodes,fSteerfile);
   }


   if ( fProcConsts.NPDF > 0 ) fScenConsts.PDF1 = INT_NS(PDF1,fSteerfile);
   if ( fProcConsts.NPDF > 1 ) fScenConsts.PDF2 = INT_NS(PDF2,fSteerfile);

   fScenConsts.ScaleVariationFactors = DOUBLE_ARR_NS(ScaleVariationFactors,fSteerfile);

}


// ___________________________________________________________________________________________________
vector<vector<pair<int,int> > > fastNLOCreate::ReadPartonCombinations(int ord, const vector<vector<int> >& PartonCombinations) {
   //! Read PDF linear combinations from steering file
   //! and convert to internal format


   //vector<vector<int> > PartonCombinations;
   if ( ord==0 ) {
      //PartonCombinations = INT_TAB_NS(PartonCombinationsLO,fSteerfile);
      if ( (int)PartonCombinations.size() != fProcConsts.NSubProcessesLO ) {
         logger.error["ReadPartonCombinations"]<<"Number of parton combinations for LO processes must be identical to number of subprocesses. NSubProcessesLO="
                                        <<fProcConsts.NSubProcessesLO<<", # parton combinations="<<PartonCombinations.size()<<". Exiting." <<endl;
         exit(1);
      }
   }
   else if ( ord==1 ) {
      //PartonCombinations = INT_TAB_NS(PartonCombinationsNLO,fSteerfile);
      if ( (int)PartonCombinations.size() != fProcConsts.NSubProcessesNLO ) {
         logger.error["ReadPartonCombinations"]<<"Number of parton combinations for NLO processes must be identical to number of subprocesses.  NSubProcessesNLO="
                                        <<fProcConsts.NSubProcessesNLO<<", # parton combinations="<<PartonCombinations.size()<<". Exiting." <<endl;
         exit(1);
      }
   }
   else if ( ord==2 ) {
      //PartonCombinations = INT_TAB_NS(PartonCombinationsNNLO,fSteerfile);
      if ( (int)PartonCombinations.size() != fProcConsts.NSubProcessesNNLO ) {
         logger.error["ReadPartonCombinations"]<<"Number of parton combinations for NNLO processes must be identical to number of subprocesses.  NSubProcessesNNLO="
                                        <<fProcConsts.NSubProcessesNNLO<<", # parton combinations="<<PartonCombinations.size()<<". Exiting." <<endl;
         exit(1);
      }
   }

   string sord[3] = {"LO","NLO","NNLO"};
   vector<vector<pair<int,int> > > PDFCoeff(PartonCombinations.size());
   // check if all partons are used
   vector<bool> b1(13);
   vector<bool> b2(13);
   for ( int i = 0 ; i<13 ; i++ ) {
      b1[i] = false;
      b2[i] = false;
   }
   for ( unsigned int k=0 ; k<PartonCombinations.size() ; k++ ) {
      if ( PartonCombinations[k].empty() ) {
         logger.error["ReadPartonCombinations"]<<"Row "<<k<<" PartonCombinations"<<sord[ord]<<" does not contain any information. Exiting."<<endl;
         exit(1);
      }
      int iSubProc = PartonCombinations[k][0];
      if ( iSubProc >= (int)PDFCoeff.size() ) {
         logger.error["ReadPartonCombinations"]<<"Subprocess "<<iSubProc<<" in row "<<k+1<<" of PartonCombinations"<<sord[ord]<<" is larger than the total number of subprocesses. Exiting."<<endl;
         exit(1);
      }
      if ( PartonCombinations[k].size()%2 != 1 || PartonCombinations[k].size()<=1) {
         logger.error["ReadPartonCombinations"]<<"Row "<<k<<" of PartonCombinations"<<sord[ord]<<" does not fit format: 'iProc [pair0] [pair1] ... [pairN]. Exiting"<<endl;
         exit(1);
      }
      if ( !PDFCoeff[iSubProc].empty() ) {
         logger.error["ReadPartonCombinations"]<<"Subprocess "<<iSubProc<<" appears twice in the PartonCombinations"<<sord[ord]<<". Exiting."<<endl;
         exit(1);
      }

      for ( unsigned int i=1 ; i<PartonCombinations[k].size() ; i+=2 ) {
         logger.debug["ReadPartonCombinations"]<<"Adding to subprocess "<<iSubProc<<" parton pair (" << PartonCombinations[k][i]<<","<<PartonCombinations[k][i+1]<<")."<<endl;
         int iPart1 = PartonCombinations[k][i];
         int iPart2 = PartonCombinations[k][i+1];
         if ( abs(iPart1) > 6 || abs(iPart2) > 6 ) {
            logger.error["ReadPartonCombinations"]<<"Parton flavor is larger than 6. There is nothing beyond the top-quark. Exiting."<<endl;
            exit(1);
         }
         // KR: Switch off warning. Using the same flavour of a parton in multiple combinations is no problem.
         // if (  b1[iPart1+6]  ) logger.warn["ReadPartonCombinations"]<<"Parton "<<iPart1<<" of hadron 1 is used multiple times in PartonCombinations"<<sord[ord]<<"."<<endl;
         // if (  b2[iPart2+6]  ) logger.warn["ReadPartonCombinations"]<<"Parton "<<iPart2<<" of hadron 2 is used multiple times in PartonCombinations"<<sord[ord]<<"."<<endl;

         b1[iPart1+6] = true;
         b2[iPart2+6] = true;
         PDFCoeff[iSubProc].push_back(std::make_pair(iPart1,iPart2));
      }
   }

   // check if all subprocesses are filled
   for ( unsigned int k=1 ; k<PDFCoeff.size() ; k++ ) {
      if ( PDFCoeff[k].empty() ) {
         logger.error["ReadPartonCombinations"]<<"PartonCombinations"<<sord[ord]<<" does not conatain any information about PDF for subprocess "<<k<<". Exiting."<<endl;
         exit(1);
      }
   }
   // check if all partons are used
   // gluons and quarks up to b, bbar must be present --> error
   for ( int p = 1 ; p<12 ; p++ ) {
      if ( !b1[p] ) {
         logger.error["ReadPartonCombinations"]<<"Parton "<<p-6<<" of hadron 1 is not used in PartonCombinations"<<sord[ord]<<". Exiting."<<endl;
         exit(1);
      }
      if ( !b2[p] ) {
         logger.error["ReadPartonCombinations"]<<"Parton "<<p-6<<" of hadron 2 is not used in PartonCombinations"<<sord[ord]<<". Exiting."<<endl;
         exit(1);
      }
   }
   // t, tbar might be absent --> info
   int p = 0;
   if ( !b1[p] ) logger.info["ReadPartonCombinations"]<<"Parton "<<p-6<<" of hadron 1 is not used in PartonCombinations"<<sord[ord]<<"."<<endl;
   if ( !b2[p] ) logger.info["ReadPartonCombinations"]<<"Parton "<<p-6<<" of hadron 2 is not used in PartonCombinations"<<sord[ord]<<"."<<endl;
   p = 12;
   if ( !b1[p] ) logger.info["ReadPartonCombinations"]<<"Parton "<<p-6<<" of hadron 1 is not used in PartonCombinations"<<sord[ord]<<"."<<endl;;
   if ( !b2[p] ) logger.info["ReadPartonCombinations"]<<"Parton "<<p-6<<" of hadron 2 is not used in PartonCombinations"<<sord[ord]<<"."<<endl;;

   return PDFCoeff;
}


// ___________________________________________________________________________________________________
void fastNLOCreate::SetGlobalVerbosity(string sverb) {
   if (sverb=="DEBUG" || sverb=="Debug" || sverb=="debug")
      speaker::SetGlobalVerbosity(say::DEBUG);
   else if (sverb=="MANUAL" || sverb=="Manual" || sverb=="manual")
      speaker::SetGlobalVerbosity(say::MANUAL);
   else if (sverb=="INFO" || sverb=="Info" || sverb=="info")
      speaker::SetGlobalVerbosity(say::INFO);
   else if (sverb=="WARNING" || sverb=="Warning" || sverb=="warning")
      speaker::SetGlobalVerbosity(say::WARNING);
   else if (sverb=="ERROR" || sverb=="Error" || sverb=="error")
      speaker::SetGlobalVerbosity(say::ERROR);
   else if (sverb=="SILENT" || sverb=="Silent" || sverb=="silent")
      speaker::SetGlobalVerbosity(say::SILENT);
   else
      speaker::SetGlobalVerbosity(say::INFO);
}


// ___________________________________________________________________________________________________
void fastNLOCreate::ReadScaleFactors() {
   //! read scale factors from steering
   //! and init member fScaleFac

   if (fIsFlexibleScale) {logger.warn["ReadScaleFactors"]<<"This function is only reasonable for fixed-scale tables!"<<endl;}
   vector<double> svar = fScenConsts.ScaleVariationFactors;// DOUBLE_ARR_NS(ScaleVariationFactors,fSteerfile);
   fScaleFac.resize(svar.size());
   if (svar.empty()) {
      // 'ScaleVariationFactors' not found -> using default
      logger.warn["ReadScaleFactors"]<<"No list of scale-factors found in steering file. Using only scale-factor of '1'."<<endl;
      fScaleFac.push_back(1.0);
   } else if (GetTheCoeffTable() && GetTheCoeffTable()->IsLO()) {
      // scale-factors not needed -> using default of 1.0
      logger.info["ReadScaleFactors"]<<"This is a leading-order run. There is no MuF scale dependence in LO. Using only scale factor of 1."<<endl;
      fScaleFac.resize(1);
      fScaleFac[0] = 1.0;
   } else if (fIsWarmup) {
      // scale-factors not needed -> using default of 1.0
      logger.info["ReadScaleFactors"]<<"This is a warmup run. Using only scale factor of 1."<<endl;
      fScaleFac.resize(1);
      fScaleFac[0] = 1.0;
   } else  {
      // sort list according to fastNLO conventions
      //  0:  1.0
      //  1-n:  nmin...nmax (without 1.0)
      vector<double> vtemp;
      bool foundunity = false;
      for (unsigned int k = 0 ; k<svar.size(); k++) {
         if (fabs(svar[k]-1.0)<1.e-8  && foundunity) {
            logger.info["ReadScaleFactors"]<<"Found scale factor 1.0 two times in list ScaleVariationFactors. Ignoring second appearance."<<endl;
            fScaleFac.resize(fScaleFac.size()-1);
         } else if (fabs(svar[k]-1.0)<1.e-8) foundunity = true;
         else vtemp.push_back(svar[k]);
      }
      if (!foundunity) {
         logger.error["ReadScaleFactors"]<<"Could not find scale factor of 1.0 in list ScaleVariationFactors. Exiting."<<endl;
         exit(1);
      }
      sort(vtemp.begin(),vtemp.end());

      fScaleFac[0] = 1.0;
      int s = 0;
      for (unsigned int k = 0 ; k<vtemp.size(); k++) {
         fScaleFac[s+1] = vtemp[k];
         if (fScaleFac[s+1] == fScaleFac[s]) {
            logger.info["ReadScaleFactors"]<<"Found scale factor '"<<fScaleFac[k+1]<<"' two times in list ScaleVariationFactors. Ignoring second appearance."<<endl;
            fScaleFac.resize(fScaleFac.size()-1);
            s--;
         }
         s++;
      }
   }
   logger.info["ReadScaleFactors"] << "Using the following scale factors:" << endl;
   for (unsigned int k = 0 ; k<fScaleFac.size(); k++)
      logger.info["ReadScaleFactors"] << "ScaleVar " << k << ": " << fScaleFac[k] << endl;
}


// ___________________________________________________________________________________________________
void fastNLOCreate::InitCoeffTable() {
   logger.debug["InitCoeffTable"]<<endl;
   //! create a coeff table
   //CreateCoeffTable(0);
   CreateCoeffTable();

   //! set 'usual' variables for perturbative calculations
   InitVariablesInCoefficientTable();

   //! read in process specific variables
   ReadCoefficientSpecificVariables();
}


// ___________________________________________________________________________________________________
int fastNLOCreate::CreateCoeffTable() {
   logger.debug["CreateCoeffTable"]<<endl;
   if (!fCoeff.empty()) {
      logger.error["CreateCoeffAddFix"]<<"Vector of coefficients must be empty, since only one coefficient table is allowed."<<endl;
      exit(1);
   }
   if (fIsFlexibleScale)
      return fastNLOTable::CreateCoeffTable(fCoeff.size(), new fastNLOCoeffAddFlex(NObsBin,ILOord));
   else
      return fastNLOTable::CreateCoeffTable(fCoeff.size(), new fastNLOCoeffAddFix(NObsBin));
}


// ___________________________________________________________________________________________________
void fastNLOCreate::ReadBinningFromScenarioConsts() {
   //! read in binning from ScenConsts

   // --- sanity checks
   if ( fScenConsts.DifferentialDimension<=0 || fScenConsts.DifferentialDimension>3  ) {
      logger.error["ReadBinningFromScenarioConsts"]<<"fScenConsts seem not be be set!"<<endl;
      exit(1);
   }
   else {
      NDim     = fScenConsts.DifferentialDimension;
      if ( NDim==1 && fScenConsts.SingleDifferentialBinning.empty() ) {
         logger.error["ReadBinningFromScenarioConsts"]<<" NDim=1 requires also a 1D binning, but fScenConsts.SingleDifferentialBinning is empty."<<endl;
         exit(1);
      }
      if ( NDim==2 && fScenConsts.DoubleDifferentialBinning.empty() ) {
         logger.error["ReadBinningFromScenarioConsts"]<<" NDim=2 requires also a 1D binning, but fScenConsts.DoubleDifferentialBinning is empty."<<endl;
         exit(1);
      }
      if ( NDim==2 && fScenConsts.DoubleDifferentialBinning[0].empty() ) {
         logger.error["ReadBinningFromScenarioConsts"]<<" NDim=2 requires also a 1D binning, but fScenConsts.DoubleDifferentialBinning[0] is empty."<<endl;
         exit(1);
      }
      if ( NDim==3 && fScenConsts.TripleDifferentialBinning.empty() ) {
         logger.error["ReadBinningFromScenarioConsts"]<<" NDim=3 requires also a 1D binning, but fScenConsts.TripleDifferentialBinning is empty."<<endl;
         exit(1);
      }
   }

   NDim     = fScenConsts.DifferentialDimension;
   DimLabel = fScenConsts.DimensionLabels;
   IDiffBin = fScenConsts.DimensionIsDifferential;
   DimLabel.resize(NDim); //safety
   IDiffBin.resize(NDim);

   // ---- check IDiffBin
   bool AllDiff = true;
   bool AllBinInt = true;
   for (unsigned int i = 0 ; i<IDiffBin.size() ; i++) {
      AllDiff = AllDiff && (IDiffBin[i] == 1);
      AllBinInt = AllBinInt && (IDiffBin[i] != 1);
   }
   if (!AllDiff && !AllBinInt) {
      logger.error["ReadBinningFromScenarioConsts"]<<"All dimensions must be consistently either bin-integrated, or truly differential dimensions. Exiting."<<endl;
      exit(1);
   }

   // ---- read bin grid
   if (NDim == 1)
      SetBinningND(fScenConsts.SingleDifferentialBinning, NDim, IDiffBin);
   else if (NDim==2)
      SetBinningND(fScenConsts.DoubleDifferentialBinning, NDim, IDiffBin);
   else if (NDim==3) {
      logger.error["ReadBinningFromScenarioConsts"]<<"The code for reading of "<<NDim<<"-dimensional binnings from ScenarioConstants is not implemented."<<endl;
      vector<vector<double> > in = DOUBLE_TAB_NS(TripleDifferentialBinning,fSteerfile); //backward compatibility! not compatible with constructor(genConst,procConst,ScenConst)
      SetBinningND(in, NDim, IDiffBin);
      //exit(3);
   }

   // ---- Bin width
   ReadBinSize();

   logger.info["ReadBinningFromScenarioConsts"]<<"Read in successfully "<<NDim<<"-dimensional bin grid with "<<NObsBin<<" bins."<<endl;

}


// ___________________________________________________________________________________________________
void fastNLOCreate::ReadBinning() {

   logger.warn["ReadBinning"]<<"This function is deprecated and should not be called. Exiting."<<endl;
   exit(23);

   // optimize read-in of bin grids
   // ToDo. Check sanity of bin-grid

   NDim         = INT_NS(DifferentialDimension,fSteerfile);
   //Scenario.SetNDim(NDim);
   if (STRING_ARR_NS(DimensionLabels,fSteerfile).size() < NDim) {
      logger.error["ReadBinning"]<<"Each dimension needs a bin label. Exiting."<<endl;
      exit(1);
   }
   DimLabel     = STRING_ARR_NS(DimensionLabels,fSteerfile);
   if (INT_ARR_NS(DimensionIsDifferential,fSteerfile).size() < NDim) {
      logger.error["ReadBinning"]<<"Each dimension need to specify if differential or not. Exiting."<<endl;
      exit(1);
   }
   IDiffBin     = INT_ARR_NS(DimensionIsDifferential,fSteerfile);
   DimLabel.resize(NDim); //safety
   IDiffBin.resize(NDim);

   // ---- check IDiffBin
   bool AllDiff = true;
   bool AllBinInt = true;
   for (unsigned int i = 0 ; i<IDiffBin.size() ; i++) {
      AllDiff = AllDiff && (IDiffBin[i] == 1);
      AllBinInt = AllBinInt && (IDiffBin[i] != 1);
   }
   if (!AllDiff && !AllBinInt) {
      logger.error["ReadBinning"]<<"All dimensions must be consistently either bin-integrated, or truly differential dimensions. Exiting."<<endl;
      exit(1);
   }
   if (AllDiff && NDim == 3) {
      logger.error["ReadBinning"]<<"Fully differential and triple-differential binning not yet implemented. exiting"<<endl;
      exit(1);
   }

   // ---- read single-differential bin grid
   if (NDim == 1) {
      vector<double> bgrid = DOUBLE_ARR_NS(SingleDifferentialBinning,fSteerfile);
      //      SetBinning1D(bgrid, DimLabel[0], IDiffBin[0]);
      SetBinningND(bgrid, NDim, IDiffBin);
   }

   // ---- read double-differential bin grid
   else if (NDim==2) {
      vector<vector<double> > in = DOUBLE_TAB_NS(DoubleDifferentialBinning,fSteerfile);
      // New binning code
      SetBinningND(in, NDim, IDiffBin);
   }

   // ---- read in triple-differential binning
   else if (NDim==3) {
      logger.warn["ReadBinning"]<<"The code for reading of "<<NDim<<"-dimensional binnings was not fully tested. Please verify the code and remove this statement."<<endl;
      vector<vector<double> > in = DOUBLE_TAB_NS(TripleDifferentialBinning,fSteerfile);
      // New binning code
      SetBinningND(in, NDim, IDiffBin);
   } else {
      logger.error["ReadBinning"]<<"Reading of "<<NDim<<"-binnings from steering is not yet implemented. Exiting"<<endl;
      exit(1);
   }

   // ---- Bin width
   ReadBinSize();

   logger.info["ReadBinning"]<<"Read in successfully "<<NDim<<"-dimensional bin grid with "<<NObsBin<<" bins."<<endl;
}


// ___________________________________________________________________________________________________
void fastNLOCreate::ReadBinSize(){
   //! initialize BinSize
   //! either from steering or from fScenConsts

   // ---------------------------------
   //       Bin size
   // ---------------------------------
   if ( fScenConsts.CalculateBinSize ) {
      BinSize.resize(NObsBin);
      bool idi = false;
      for (unsigned int i = 0 ; i<NObsBin ; i++) {
         BinSize[i] = 1;
         for (unsigned int d=0 ; d<NDim ; d++) {
            if (IDiffBin[d]==0) {
               // nothing todo
            } else if (IDiffBin[d]==1) {
               // nothing todo
               //logger.warn["ReadBinning"]<<"Don't know how to handle truly differential bins for bin widths."<<endl;
            } else if (IDiffBin[d]==2) {
               BinSize[i] *= Bin[i][d].second-Bin[i][d].first;//UpBin[i][d]-LoBin[i][d];
               idi = true;
            }
         }
         // divide by binsizefactor, but only if at least one dimension is differential
         if (idi) BinSize[i] *= fScenConsts.BinSizeFactor;
      }
      if (!idi) logger.debug["ReadBinning"]<<"BinSizeFactor is not being used, since no observable is calculated differential."<<endl;
   }
   else {
      // read in bin width
      logger.warn["ReadBinning"]<<"Reading of 'BinSize' only poorly  implemented! Improve it and remove this message."<<endl;
      BinSize = fScenConsts.BinSize;//DOUBLE_ARR_NS(BinSize,fSteerfile);
      if ( BinSize.size()!=NObsBin ) logger.warn["ReadBinning"]<<"Number of bins of 'BinSize' not consistent with bin grid."<<endl;
      BinSize.resize(NObsBin);
      for (unsigned int i = 0 ; i<NObsBin ; i++) {
         if (BinSize[i]==0) BinSize[i] = 1.0;
      }
   }

}

///////////////////////////////////////////////

// ___________________________________________________________________________________________________

//
// Set continuous 1D binning via one vector containing all (n+1) bin edges
//
void fastNLOCreate::SetBinning1D(vector<double> bgrid, string label, unsigned int idiff) {

   // Create and set normalization vector to one
   vector<double> vnorm;
   unsigned int nbins = bgrid.size()-1;
   if ( idiff == 1 ) { nbins = bgrid.size(); }// Point-wise differential
   vnorm.assign(nbins,1.);

   // Give task to more general method
   SetBinning1D(bgrid, label, idiff, vnorm);
   logger.info["SetBinning1D"] << "VSI: Set all normalization factors to one." << endl;
}

//
// Set continuous 1D binning via one vector containing all (n+1) bin edges with normalization factor
//
void fastNLOCreate::SetBinning1D(vector<double> bgrid, string label, unsigned int idiff, double norm) {

   // Create and set normalization vector to norm
   vector<double> vnorm;
   unsigned int nbins = bgrid.size()-1;
   if ( idiff == 1 ) { nbins = bgrid.size(); }// Point-wise differential
   vnorm.assign(nbins,norm);

   // Give task to more general method
   SetBinning1D(bgrid, label, idiff, vnorm);
   logger.info["SetBinning1D"] << "VSID: Set all normalization factors to norm." << endl;
}

//
// Set continuous 1D binning via one vector containing all (n+1) bin edges with normalization factors in extra vector
//
void fastNLOCreate::SetBinning1D(vector<double> bgrid, string label, unsigned int idiff, vector<double> vnorm) {

   // Create and set two vectors with lower and upper bin edges
   vector<double> blow;
   vector<double> bupp;
   for (unsigned int i = 0; i<bgrid.size()-1; i++) {
      blow.push_back(bgrid[i]);
      bupp.push_back(bgrid[i+1]);
   }

   // Give task to more general method
   SetBinning1D(blow, bupp, label, idiff, vnorm);
   logger.info["SetBinning1D"] << "VSIV: Set vectors of lower and upper bin edges." << endl;
}

//
// Set 1D binning allowing gaps via two vectors containing lower and upper bin edges
//
void fastNLOCreate::SetBinning1D(vector<double> blow, vector<double> bupp, string label, unsigned int idiff) {

   // Create and set normalization vector to one
   vector<double> vnorm;
   vnorm.assign(blow.size(),1.);

   // Give task to more general method
   SetBinning1D(blow, bupp, label, idiff, vnorm);
   logger.info["SetBinning1D"] << "VVSI: Set all normalization factors to one." << endl;
}

//
// Set 1D binning allowing gaps via two vectors containing lower and upper bin edges with normalization factor
//
void fastNLOCreate::SetBinning1D(vector<double> blow, vector<double> bupp, string label, unsigned int idiff, double norm) {

   // Create and set normalization vector to norm
   vector<double> vnorm;
   vnorm.assign(blow.size(),norm);

   // Give task to more general method
   SetBinning1D(blow, bupp, label, idiff, vnorm);
   logger.info["SetBinning1D"] << "VVSID: Set all normalization factors to norm." << endl;
}

//
// Set 1D binning allowing gaps via two vectors containing lower and upper bin edges with normalization factors in extra vector
//
void fastNLOCreate::SetBinning1D(vector<double> blow, vector<double> bupp, string label, unsigned int idiff, vector<double> vnorm) {

   // Define 1-dimensional observable binning
   NDim = 1;
   DimLabel.resize(NDim);
   IDiffBin.resize(NDim);

   // Check input: Type of differential binning
   DimLabel[0] = label;
   IDiffBin[0] = idiff;
   if (idiff > 2 ) {
      logger.error["SetBinning1D"] << "Illegal choice of differentiality, idiff = " << idiff << ", aborted!" << endl;
      logger.error["SetBinning1D"] << "Dimension must either be non-differential (idiff=0), " << endl;
      logger.error["SetBinning1D"] << "point-wise differential (idiff=1), or bin-wise differential (idiff=2)." << endl;
      exit(1);
   }

   // Check input: Gaps are fine, overlaps are not. Overlaps might lead to double-counting in normalization.
   if ( blow.size() != bupp.size() ) {
      logger.error["SetBinning1D"] << "Unequal size of bin border vectors, blow.size() = " << blow.size() <<
         ", bupp.size() = " << bupp.size() << ", aborted!" << endl;
      exit(1);
   }
   for (unsigned int i = 0; i<blow.size(); i++) {
      if (idiff == 1) { // Point-wise differential
         if (blow[i] != bupp[i]) {
            logger.error["SetBinning1D"] << "Illegal binning, lower and upper bin edges must be equal for point-wise differential binnings, aborted!" << endl;
            logger.error["SetBinning1D"] << "i = " << i << ", lower edge = " << blow[i] <<
               ", and upper edge = " << bupp[i] << endl;
            exit(1);
         }
      } else { // Non- and bin-wise differential
         if (blow[i] >= bupp[i]) {
            logger.error["SetBinning1D"] << "Illegal binning, lower bin edge larger or equal to upper one, aborted!" << endl;
            logger.error["SetBinning1D"] << "i = " << i << ", lower edge = " << blow[i] <<
               ", and upper edge = " << bupp[i] << endl;
            exit(1);
         }
         if (i<blow.size()-1 && bupp[i] > blow[i+1]) {
            logger.error["SetBinning1D"] << "Illegal binning, overlapping bins, aborted!" << endl;
            logger.error["SetBinning1D"] << "i = " << i << ", upper edge = " << bupp[i] <<
               ", and lower edge of i+1 = " << blow[i+1] << endl;
            exit(1);
         }
      }
   }

   // Check input: Normalization factor must be larger than zero for all bins
   if ( blow.size() != vnorm.size() ) {
      logger.error["SetBinning1D"] << "Unequal size of bin border and normalization vectors, blow.size() = " << blow.size() <<
         ", vnorm.size() = " << vnorm.size() << ", aborted!" << endl;
      exit(1);
   }
   for (unsigned int i = 0; i<vnorm.size(); i++) {
      if ( vnorm[i] < DBL_MIN ) {
         logger.error["SetBinning1D"] << "Normalization factor is too small, aborted!" << endl;
         logger.error["SetBinning1D"] << "i = " << i << ", vnorm[i] = " << vnorm[i] << endl;
         exit(1);
      }
   }

   // Set differential binning and normalization
   NObsBin = blow.size();
   Bin.resize(NObsBin);
   BinSize.resize(NObsBin);
   for (unsigned int i = 0; i<blow.size(); i++) {
      Bin[i].resize(1);
      Bin[i][0] = make_pair(blow[i],bupp[i]);
      BinSize[i] = 1;
      if (idiff == 2) { // Divide by bin width times additional normalization factor
         BinSize[i] *= (Bin[i][0].second-Bin[i][0].first) * vnorm[i];
      }
   }
   logger.info["SetBinning1D"] << "VVSIV: Set binning successfully for " << NObsBin << " bins in " << NDim << "dimensions." << endl;
}

/////////////////// ND

//
// Set continuous n-dimensional binning via a vector with bin edges
//
void fastNLOCreate::SetBinningND(vector<double> bgrid, unsigned int ndim, vector<int> idiff) {

   // Create and set vector of vectors with bin edges (pairs)
   vector<vector<double> > bgrid2;
   bgrid2.push_back(bgrid);

   // Give task to more general method
   SetBinningND(bgrid2, ndim, idiff);
   logger.info["SetBinningND"] << "VIV: Set binning via vector with bin edges." << endl;
}


//
// Set continuous n-dimensional binning via a vector of vectors with bin edges (pairs)
//
void fastNLOCreate::SetBinningND(vector<vector<double> > bgrid, unsigned int ndim, vector<int> idiff) {

   // Check input: Dimensions, labels and idiff settings
   if ( ndim < 1 ) {
      logger.error["SetBinningND"] << "Illegal # of dimensions, aborted!" << endl;
      logger.error["SetBinningND"] << "ndim = " << ndim << endl;
      exit(1);
   } else if ( ndim > 3 ) {
      logger.error["SetBinningND"] << "More than three dimensions not yet supported, aborted!" << endl;
      logger.error["SetBinningND"] << "ndim = " << ndim << endl;
      exit(1);
   }
   if ( ndim != idiff.size() ) {
      logger.error["SetBinningND"] << "Inconsistency between # of dimensions and vector length for differentiality settings, aborted!" << endl;
      logger.error["SetBinningND"] << "Number of dimensions ndim = " << ndim << ", # of idiff settings = " << idiff.size() << endl;
      exit(1);
   }

   // Check input: Type of differential binning
   for ( unsigned int i=0; i<ndim; i++ ) {
      if ( idiff[i] < 0 || idiff[i] > 2 ) {
         logger.error["SetBinningND"] << "Illegal choice of differentiality, idim = " << i << ", idiff = " << idiff[i] << ", aborted!" << endl;
         logger.error["SetBinningND"] << "Each dimension must either be non-differential (idiff=0), " << endl;
         logger.error["SetBinningND"] << "point-wise differential (idiff=1), or bin-wise differential (idiff=2)." << endl;
         exit(1);
      }
   }

   // Write ND binning to vector[NObsBin] with bin edges in pairs, one for each dimension.
   // (For point-wise differential distributions the two bin edges are set to be identical!)
   // Initialize storage structure Bin
   NObsBin = 0;
   Bin.clear();
   // Each 'line' of numbers start with pairs of bin edges (or one bin center) for each dimension 0 ... n-1
   // followed by strictly monotonously increasing bin borders as for a 1-dimensional histogram.
   // (In principal, each observable bin can be defined completely independently, but the pattern used
   //  here covers the most frequent use cases in a more comfortable way.)
   // Initialize shift that indicates where the innermost dimension with histogram-like borders starts
   int ishift = 0;
   // Allow differentiation between increment by two for non- or bin-wise differential booking and
   // increment by one for point-wise differential booking
   for ( unsigned int i = 0; i<ndim-1; i++ ) { // No shift for 1D
      // int istep = 2;                        // Increment by two for non- or bin-wise differential booking
      // if ( idiff[i] == 1 ) {istep = 1;}     // else increment by one for point-wise differential
      // ishift   += istep;
      ishift += ( ( idiff[i] == 1 ) ? 1 : 2 ); // Try ternary operator for this
   }

   // Loop over input vectors with two (one, for idiff=1) entries per dimension plus 'n+1' entries for 'n' bins of innermost dimension
   for ( unsigned int i = 0; i<bgrid.size(); i++ ) {
      logger.debug["SetBinningND"] << "i = " << i << ", iObsBin = " << NObsBin << ", ishift = " << ishift << endl;

      // Loop over 'n' 'histogram' bins
      for (unsigned int j = ishift ; j<bgrid[i].size()-1; j++) {
         Bin.push_back(vector<pair<double,double> >(ndim));

         // Keep track of additional increment in vector for non- or bin-wise differential of the 'n-1' outer dimensions
         int kshift = 0;
         for ( unsigned int k = 0; k<ndim-1; k++ ) {
            logger.debug["SetBinningND"] << "i = " << i << ", j = " << j << ",k+kshift = " << k+kshift << ", bgrid i,k+kshift = " << bgrid[i][k+kshift] << endl;
            logger.debug["SetBinningND"] << "i = " << i << ", j = " << j << ",k+kshift+1 = " << k+kshift+1 << ", bgrid i,k+kshift+1 = " << bgrid[i][k+kshift+1] << endl;
            if ( idiff[k] == 1 ) { // Point-wise differential --> one bin center, set bin edges identical
               logger.warn["SetBinningND"] << "Point-wise differential binning implemented, but not tested extensively. Please check carefully!" << endl;
               Bin[NObsBin][k] = make_pair(bgrid[i][k+kshift],bgrid[i][k+kshift]);
            } else {               // Non- or bin-wise differential --> two bin edges
               if ( !(bgrid[i][k+kshift] < bgrid[i][k+kshift+1]) ) {
                  logger.error["SetBinningND"] << "Illegal binning, lower bin edge larger or equal to upper one, aborted!" << endl;
                  logger.error["SetBinningND"] << "i, k = " << i << ", " << k << ", lower edge = " << bgrid[i][k+kshift] <<
                     ", and upper edge = " << bgrid[i][k+kshift+1] << endl;
                  exit(1);
               }
               Bin[NObsBin][k] = make_pair(bgrid[i][k+kshift],bgrid[i][k+kshift+1]);
               kshift++;
            }
         }

         if ( idiff[ndim-1] == 1 ) { // Point-wise differential --> one bin center
            logger.warn["SetBinningND"] << "Point-wise differential binning implemented, but not tested extensively. Please check carefully!" << endl;
            Bin[NObsBin][ndim-1] = make_pair(bgrid[i][j],bgrid[i][j]);
         } else {                    // Non- or bin-wise differential --> two bin edges
            logger.debug["SetBinningND"] << "i = " << i << ", j = " << j << ", bgrid i,j = " << bgrid[i][j] << endl;
            logger.debug["SetBinningND"] << "i = " << i << ", j = " << j << ", bgrid i,j+1 = " << bgrid[i][j+1] << endl;
            if ( !(bgrid[i][j] < bgrid[i][j+1]) ) {
               logger.error["SetBinningND"] << "Illegal binning, lower bin edge larger or equal to upper one, aborted!" << endl;
               logger.error["SetBinningND"] << "i, j = " << i << ", " << j << ", lower edge = " << bgrid[i][j] <<
                  ", and upper edge = " << bgrid[i][j+1] << endl;
               exit(1);
            }
            Bin[NObsBin][ndim-1] = make_pair(bgrid[i][j],bgrid[i][j+1]);
         }
         NObsBin++;
      }
   }

   // Check on monotonously increasing bin edges of 'n-1' outer dimensions
   // Jumps in more inner dimension allowed only at bin edge of more outer dimension
   // This rule is not required in principal, but for the use of comfort functions for histogramming
   if ( NObsBin > 1 ) {
      for ( unsigned int i=0; i<NObsBin-1; i++ ) {
         if ( (Bin[i+1][0].first - Bin[i][0].first) < -DBL_MIN ) {
            logger.error["SetBinningND"] << "Illegal binning, lower bin edge larger than upper one for outermost dimension, aborted!" << endl;
            logger.error["SetBinningND"] << "The problematic observable bins are " << i << " and " << i+1 << endl;
            logger.error["SetBinningND"] << "The bin edges for observable bin " << i << " are: Bin[i][0].first = " << Bin[i][0].first << ", Bin[i][0].second = " << Bin[i][0].second << endl;
            logger.error["SetBinningND"] << "The bin edges for observable bin " << i+1 << " are: Bin[i+1][0].first = " << Bin[i+1][0].first << ", Bin[i+1][0].second = " << Bin[i+1][0].second << endl;
            exit(1);
         }
         bool lordered = false;
         for ( unsigned int j=0; j<ndim; j++ ) {
            logger.debug["SetBinningND"] << "Before: iobs = " << i << ", jdim = " << j << ", lordered = " << lordered << ", Bin i,j = " << Bin[i][j].first << ", Bin i+1,j = " << Bin[i+1][j].first << endl;
            //            if ( !lordered && Bin[i][j].first < Bin[i+1][j].first ) {lordered = true;}
            if ( !lordered && (Bin[i+1][j].first - Bin[i][j].first) > DBL_MIN ) {lordered = true;}
            logger.debug["SetBinningND"] << "After : iobs = " << i << ", jdim = " << j << ", lordered = " << lordered << ", Bin i,j = " << Bin[i][j].first << ", Bin i+1,j = " << Bin[i+1][j].first << endl;
         }
         if (! lordered) {
            logger.error["SetBinningND"] << "Illegal binning, lower bin edge larger or equal to upper one, aborted!" << endl;
            logger.error["SetBinningND"] << "The problematic observable bins are " << i << " and " << i+1 << endl;
            for ( unsigned int j=0; j<ndim; j++ ) {
               logger.error["SetBinningND"] << "The bin edges for observable bin " << i << " in dimension " << j << " are: Bin[i][j].first = " << Bin[i][j].first << ", Bin[i][j].second = " << Bin[i][j].second << endl;
               logger.error["SetBinningND"] << "The bin edges for observable bin " << i+1 << " in dimension " << j << " are: Bin[i+1][j].first = " << Bin[i+1][j].first << ", Bin[i+1][j].second = " << Bin[i+1][j].second << endl;

            }
            exit(1);
         }
      }
   }

   // Test new binning code
   // unsigned int ndim0 = 0;
   // unsigned int ndim1 = 0;
   // unsigned int ndim2 = 0;
   // unsigned int ndim0b = 0;
   // unsigned int ndim1b = 0;
   // unsigned int ndim2b = 0;
   // std::vector< std::pair<double, double > > DimBins;
   // if (NDim > 0) {
   //    ndim0   = GetNDim0Bins();
   //    DimBins = GetDim0Bins();
   //    ndim0b  = DimBins.size();
   //    cout << "ndim0, ndim0b = " << ndim0 << ", " << ndim0b << endl;
   //    for (unsigned int k=0; k<ndim0; k++) {
   //       cout << "dim0 k, low, up  " << k << ", " << DimBins[k].first << ", " << DimBins[k].second << endl;
   //    }
   // }
   // if (NDim > 1) {
   //    for ( unsigned int i=0; i<ndim0; i++ ) {
   //       ndim1   = GetNDim1Bins(i);
   //       DimBins = GetDim1Bins(i);
   //       ndim1b  = DimBins.size();
   //       cout << "ndim1, ndim1b = " << ndim1 << ", " << ndim1b << endl;
   //       for (unsigned int k=0; k<ndim1; k++) {
   //          cout << "dim1 k, low, up   " << k << ", " << DimBins[k].first << ", " << DimBins[k].second << endl;
   //       }
   //    }
   // }
   // if (NDim > 2) {
   //    for ( unsigned int i=0; i<ndim0; i++ ) {
   //       for ( unsigned int j=0; j<ndim1; j++ ) {
   //          ndim2   = GetNDim2Bins(i,j);
   //          DimBins = GetDim2Bins(i,j);
   //          ndim2b  = DimBins.size();
   //          cout << "ndim2, ndim2b = " << ndim2 << ", " << ndim2b << endl;
   //          for (unsigned int k=0; k<ndim2; k++) {
   //             cout << "dim2 k, low, up   " << k << ", " << DimBins[k].first << ", " << DimBins[k].second << endl;
   //          }
   //       }
   //    }
   // }

   // unsigned int nbins = DimBins.size();
   // cout << "ZZZ iDim = " << i << ", nbins = " << nbins << endl;
   // for ( unsigned int j=0; j<nbins; j++) {
   //    cout << "ZZZ uniqbins first = " << DimBins[j].first << ", uniqbins second = " << DimBins[j].second << endl;
   // }

   // for ( unsigned int i = 0; i<NObsBin; i++ ) {
   //    cout << "dim0 lo/up = " << Bin[i][0].first << ", " << Bin[i][0].second
   //         << ", dim1 lo/up = " << Bin[i][1].first << ", " << Bin[i][1].second << endl;
   //    //           << ", dim2 lo/up = " << Bin[i][2].first << ", " << Bin[i][2].second << endl;
   //    if (NDim == 1) {
   //       cout << "obs0 = " << Bin[i][0].second << ", Dim0Bin = " << GetODim0Bin(Bin[i][0].second) << endl;
   //    }
   //    if (NDim == 2) {
   //       cout << "obs0 = " << Bin[i][0].second << ", obs1 = " << Bin[i][1].second << ", Dim0Bin = " << GetODim0Bin(Bin[i][0].second) << ", Dim1Bin = " << GetODim1Bin(Bin[i][0].second,Bin[i][1].second) << endl;
   //       for ( unsigned int j=0; j<GetNDim0Bins(); j++ ) {
   //       }
   //    }
   //    if (NDim == 3) {
   //       cout << "obs0 = " << Bin[i][0].second << ", obs1 = " << Bin[i][1].second << ", obs2 = " << Bin[i][2].second << ", Dim0Bin = " << GetODim0Bin(Bin[i][0].second) << ", Dim1Bin = " << GetODim1Bin(Bin[i][0].second,Bin[i][1].second) << ", Dim2Bin = " << GetODim2Bin(Bin[i][0].second,Bin[i][1].second,Bin[i][2].second) << endl;
   //       for ( unsigned int j=0; j<GetNDim0Bins(); j++ ) {
   //       }
   //       for ( unsigned int j=0; j<GetNDim0Bins(); j++ ) {
   //          for ( unsigned int k=0; k<GetNDim1Bins(j); k++ ) {
   //          }
   //       }
   //    }
   // double LowerBoundary = GetBinBoundaries(i)[0].first;
   // double UpperBoundary = GetBinBoundaries(i)[0].second;
   // cout << "iobs = " << i << ", low = " << LowerBoundary << ", upp = " << UpperBoundary << endl;
   // int i1 = GetIDim0Bin(i);
   // int i2 = GetIDimBin(i,0);
   // int i3 = GetIDim1Bin(i);
   // int i4 = GetIDimBin(i,1);
   // int i5 = GetIDim2Bin(i);
   // int i6 = GetIDimBin(i,2);
   // cout << "i1, i2 = " << i1 << ", " << i2 << ", i3, i4 = " << i3 << ", " << i4 << ", i5, i6 = " << i5 << ", " << i6 << endl;
   //   }

   logger.info["SetBinningND"] << "Set binning successfully for " << NObsBin << " bins in " << NDim << " dimensions." << endl;
}


/////////////////////////////////////////////


// ___________________________________________________________________________________________________
void fastNLOCreate::GetWarmupValues() {
   //!
   //! GetWarmupValues.
   //! Checks if warmup-table exists and initialized
   //! member variable fIsWarmup
   //!
   logger.debug["GetWarmupValues"]<<endl;

   if ( !fWarmupConsts.Values.empty() ) {
      //! Check if warmup values already set by user
      logger.info["GetWarmupValues"]<<"Found warmup values in fWarmupConsts without reading warmup-steering file."<<endl;
      fIsWarmup=false;
   }
   else {
      //! Try to get warmup values from steering
      std::cout.setstate(std::ios::failbit) ; // no cout in the following
      // std::cerr.setstate(std::ios::failbit) ; // no cout in the following
      logger.info >>"\n";
      logger.info >> (fastNLO::_SSEP40+fastNLO::_SSEP40+fastNLO::_SSEP40) << endl;
      logger.info["GetWarmupValues"]<<"Trying to get warmup values. Please ignore following messages from parser."<<endl;
      // try to get warmup values
      //vector<vector<double> > warmup = DOUBLE_TAB_NS(Warmup.Values,fSteerfile);
      fWarmupConsts.Values = DOUBLE_TAB_NS(Warmup.Values,fSteerfile);
      fIsWarmup = fWarmupConsts.Values.empty();

      // try again, with hard-coded convention:
      if (fIsWarmup) {
         logger.debug["GetWarmupValues"]<<"Could not get warmup table from steerfile. Now trying to read steerfile: "<<GetWarmupTableFilename()<<endl;
         bool t1 = access(GetWarmupTableFilename().c_str(), R_OK);
         usleep(100); // short interruption
         bool t2 = access(GetWarmupTableFilename().c_str(), R_OK);
         if ( t1 != 0  && t2 !=0 ) {
            logger.debug["GetWarmupValues"]<<"Warmup file does not exist: "<<GetWarmupTableFilename()<<endl;
            fIsWarmup=true;
         }
         else {
            READ_NS(GetWarmupTableFilename(),fSteerfile);    // put the warmup-values into same read_steer 'namespace'
            fWarmupConsts.Values = DOUBLE_TAB_NS(Warmup.Values,fSteerfile);
            fIsWarmup = fWarmupConsts.Values.empty();
            if (!fIsWarmup)
               logger.info["GetWarmupValues"]<<"Warmup values found in file "<<GetWarmupTableFilename()<<"."<<endl;
         }
      }

      // --- read other warmup paramters
      if ( !fIsWarmup ) {
         fWarmupConsts.Binning = DOUBLE_TAB_NS(Warmup.Binning,fSteerfile);
         fWarmupConsts.OrderInAlphasOfWarmupRunWas = INT_NS(Warmup.OrderInAlphasOfWarmupRunWas,fSteerfile);
         fWarmupConsts.CheckScaleLimitsAgainstBins = BOOL_NS(Warmup.CheckScaleLimitsAgainstBins,fSteerfile);
         // fWarmupConsts.Values = DOUBLE_TAB_NS(Warmup.Values,fSteerfile);
         fWarmupConsts.headerValues = TABLEHEADER_NS(Warmup.Values,fSteerfile);
         //
         fWarmupConsts.DifferentialDimension = INT_NS(Warmup.DifferentialDimension,fSteerfile);
         fWarmupConsts.DimensionIsDifferential = INT_ARR_NS(Warmup.DimensionIsDifferential,fSteerfile);
         fWarmupConsts.DimensionLabels = STRING_ARR_NS(Warmup.DimensionLabels,fSteerfile);
         fWarmupConsts.ScaleDescriptionScale1 = STRING_NS(Warmup.ScaleDescriptionScale1,fSteerfile);
         if (!fIsFlexibleScale)
            fWarmupConsts.ScaleDescriptionScale2 = STRING_NS(Warmup.ScaleDescriptionScale2,fSteerfile);

      }

      // --- read in remaining scenario constants if requested
      if ( fIsWarmup && !fScenConsts.ReadBinningFromSteering ) {
         logger.error["Instantiate"]<<"This is a warmup run. Thus, the binning must be read from the steering or ScenarioConstants. Please use ReadBinningFromSteering=true"<<endl;
         exit(1);
      }
      else if ( !fIsWarmup && !fScenConsts.ReadBinningFromSteering ) {
         fScenConsts.DifferentialDimension   = fWarmupConsts.DifferentialDimension;
         fScenConsts.DimensionIsDifferential = fWarmupConsts.DimensionIsDifferential;
         fScenConsts.DimensionLabels         = fWarmupConsts.DimensionLabels;
         fScenConsts.ScaleDescriptionScale1  = fWarmupConsts.ScaleDescriptionScale1;
         if (!fIsFlexibleScale)
            fScenConsts.ScaleDescriptionScale2 = fWarmupConsts.ScaleDescriptionScale2;
      }


      // inform user about success
      logger.info >> (fastNLO::_SSEP40+fastNLO::_SSEP40+fastNLO::_SSEP40) << endl;
      std::cout.clear() ; // recover cout to screen
      std::cerr.clear() ; // recover cout to screen
   }

   if ( fIsWarmup ) logger.warn["GetWarmupValues"]<<"This will be a warmup run."<<endl;
   else             logger.info["GetWarmupValues"]<<"This will be a production run."<<endl;
}



// ___________________________________________________________________________________________________
void fastNLOCreate::UseBinGridFromWarmup() {
   //! initialialize all binning related variables
   //! with values stored in the warmup file.
   //vector<vector<double> > warmup =  DOUBLE_TAB_NS(Warmup.Binning,fSteerfile);
   const vector<vector<double> >& warmup = fWarmupConsts.Binning;

   NObsBin      = warmup.size();
   NDim         = fScenConsts.DifferentialDimension;
   if (warmup[0].size() != (7+2*NDim) && warmup[0].size() != (5+2*NDim)) {
      logger.error["UseBinGridFromWarmup"]<<"This warmup table has an unknown size of columns. Expecting "<<(7+2*NDim)<<" for flexible-scale, or "<<(5+2*NDim)<<" for fixed-scale tables. Exiting."<<endl;
      exit(1);
   }
   fIsFlexibleScale = (warmup[0].size() == (7+2*NDim));
   IDiffBin     = fScenConsts.DimensionIsDifferential;
   DimLabel     = fScenConsts.DimensionLabels;

   // make binning
   const int i0 = 1;//fIsFlexibleScale ? 6 : 4;
   Bin.resize(NObsBin);
   BinSize.resize(NObsBin);
   for (unsigned int i = 0 ; i < NObsBin ; i ++) {
      Bin[i].resize(NDim);
      if (NDim==1) {
         Bin[i][0] = make_pair(warmup[i][i0],warmup[i][i0+1]) ;
      } else if (NDim==2) {
         Bin[i][0] = make_pair(warmup[i][i0],warmup[i][i0+1]) ;
         Bin[i][1] = make_pair(warmup[i][i0+2],warmup[i][i0+3]) ;
      } else if (NDim==3) {
         Bin[i][0] = make_pair(warmup[i][i0],warmup[i][i0+1]) ;
         Bin[i][1] = make_pair(warmup[i][i0+2],warmup[i][i0+3]) ;
         Bin[i][2] = make_pair(warmup[i][i0+4],warmup[i][i0+5]) ;
      }
      BinSize[i] = warmup[i][i0+NDim*2];
   }
}



// ___________________________________________________________________________________________________
bool fastNLOCreate::CheckWarmupConsistency() {
   //! check if warmup values are consistent with steering card
   //! check if number of bins is consistent

   vector<vector<double> > warmup =  fWarmupConsts.Values;//DOUBLE_TAB_NS(Warmup.Values,fSteerfile);
   bool ret = true;

   const string wrmuphelp = "Please check your warmup-file for compatibility with your steering.\nTo calculate a new warmup-file compatible to your steering, remove the old one.\nAlternatively use 'IgnoreWarmupBinningCheck=true' to ignore precision-related binning differences or\nuse 'ReadBinningFromSteering=false' to read all binning-related information from the warmup file.\n";

   // --- check warmup values
   if (warmup.size() != NObsBin) {
      logger.error["CheckWarmupConsistency"]
         <<"Table of warmup values is not compatible with steering file.\n"
         <<"Different number of bins ("<<warmup.size()<<" instead of "<<NObsBin<<".\n"
         <<wrmuphelp
         <<"Exiting."<<endl;
      ret = false;
      exit(1);
   }
   if (fWarmupConsts.DifferentialDimension != (int)NDim) {
      logger.error["CheckWarmupConsistency"]
         <<"Table of warmup values is not compatible with steering file.\n"
         <<"Found different number of dimensions. NDim="<<NDim<<", Warmup.DifferentialDimension="<<fWarmupConsts.DifferentialDimension<<".\n"
         <<wrmuphelp
         <<"Exiting."<<endl;
      ret = false;
      exit(1);
   }

   // CoeffTable is not available during intialization
   //    if ( STRING_NS(Warmup.ScaleDescriptionScale1) != GetTheCoeffTable()->ScaleDescript[0][0] ) {
   //       logger.warn["CheckWarmupConsistency"]
   //    <<"Table of warmup values is potentially incompatible with steering file.\n"
   //    <<"Found different scale description (ScaleDescriptionScale1). ScaleDescriptionScale1="<<GetTheCoeffTable()->ScaleDescript[0][0]
   //    <<", but Warmup.ScaleDescriptionScale1='"<<STRING_NS(Warmup.ScaleDescriptionScale1)<<"'."<<endl<<endl;
   //       ret = false;
   //    }
   //    if ( fIsFlexibleScale && STRING_NS(Warmup.ScaleDescriptionScale2) != GetTheCoeffTable()->ScaleDescript[0][1] ){
   //       logger.warn["CheckWarmupConsistency"]
   //    <<"Table of warmup values is potentially incompatible with steering file.\n"
   //    <<"Found different scale description (ScaleDescriptionScale2). ScaleDescriptionScale2='"<<GetTheCoeffTable()->ScaleDescript[0][0]<<"'"
   //    <<", but Warmup.ScaleDescriptionScale2='"<<STRING_NS(Warmup.ScaleDescriptionScale1)<<"'."<<endl<<endl;
   //       ret = false;
   //    }

   if (fWarmupConsts.DimensionIsDifferential[0] != IDiffBin[0]) {
      logger.warn["CheckWarmupConsistency"]
         <<"Table of warmup values seems to be incompatible with steering file.\n"
         <<"Found different diff-label for dimension 0  (IDiffBin). DimensionIsDifferential='"<<IDiffBin[0]<<"'"
         <<", but Warmup.DimensionIsDifferential[0]='"<<fWarmupConsts.DimensionIsDifferential[0]<<"'. Exiting."<<endl;
      ret = false;
      exit(1);
   }
   if (NDim > 1 && fWarmupConsts.DimensionIsDifferential[1] != IDiffBin[1]) {
      logger.warn["CheckWarmupConsistency"]
         <<"Table of warmup values seems to be incompatible with steering file.\n"
         <<"Found different diff-label for dimension 1  (IDiffBin). DimensionIsDifferential='"<<IDiffBin[1]<<"'"
         <<", but Warmup.DimensionIsDifferential[1]='"<<fWarmupConsts.DimensionIsDifferential[1]<<"'. Exiting."<<endl;
      ret = false;
      exit(1);
   }
   if (NDim > 2 && fWarmupConsts.DimensionIsDifferential[2] != IDiffBin[2]) {
      logger.warn["CheckWarmupConsistency"]
         <<"Table of warmup values seems to be incompatible with steering file.\n"
         <<"Found different diff-label for dimension 2  (IDiffBin). DimensionIsDifferential='"<<IDiffBin[2]<<"'"
         <<", but Warmup.DimensionIsDifferential[2]='"<<fWarmupConsts.DimensionIsDifferential[2]<<"'. Exiting."<<endl;
      ret = false;
      exit(1);
   }

   // --- check warmup binning
   const vector<vector<double> >& wrmbin =  fWarmupConsts.Binning;//DOUBLE_TAB_NS(Warmup.Binning,fSteerfile);
   // check binning in detail; ignore when IgnoreWarmupBinningCheck is true to avoid precision problems with bin borders of e.g. Pi
   if ( EXIST_NS(IgnoreWarmupBinningCheck,fSteerfile) && BOOL_NS(IgnoreWarmupBinningCheck,fSteerfile) ) {
      logger.warn["CheckWarmupConsistency"]
         <<"Ignoring crosscheck of binning in steering versus warmup values. "
         <<"Please make sure that your warmup file matches the steering of this run."<<endl;
   }
   else if ( !wrmbin.empty() ) {
      for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
         const int i0 = 1;//fIsFlexibleScale ? 6 : 4;
         const double EPS = 1.e-8;
         if (NDim == 1) {
            if ( fabs(Bin[i][0].first- wrmbin[i][i0])>EPS || fabs(Bin[i][0].second - wrmbin[i][i0+1]) > EPS) {
               logger.error["CheckWrmbinConsistency"]
                  <<"Table of warmup values seems to be incompatible with steering file.\n"
                  <<"Found different binning for bin "<<i<<", steering: ["<<Bin[i][0].first<<","<<Bin[i][0].second
                  <<",], warmup: ["<<wrmbin[i][i0]<<","<<wrmbin[i][i0+1]<<"].\n"
                  <<wrmuphelp
                  <<"Exiting."<<endl;
               ret = false;
               //exit(1);
            }
         } else if (NDim == 2) {
            if (fabs(Bin[i][0].first - wrmbin[i][i0]) > EPS || fabs(Bin[i][0].second - wrmbin[i][i0+1]) > EPS
                || fabs(Bin[i][1].first - wrmbin[i][i0+2]) > EPS || fabs(Bin[i][1].second - wrmbin[i][i0+3]) > EPS ) {
               logger.error["CheckWarmupConsistency"]
                  <<"Table of warmup values seems to be incompatible with steering file.\n"
                  <<"Found different binning for bin "<<i<<", steering: ["<<Bin[i][0].first<<","<<Bin[i][0].second<<",] ["<<Bin[i][1].first<<","<<Bin[i][1].second
                  <<"], warmup: ["<<wrmbin[i][i0]<<","<<wrmbin[i][i0+1]<<"] ["<<wrmbin[i][i0+2]<<","<<wrmbin[i][i0+3]<<"].\n"
                  <<wrmuphelp
                  <<"Exiting."<<endl;
               ret = false;
               //exit(1);
            }
         } else if (NDim == 3) {
            if (fabs(Bin[i][0].first - wrmbin[i][i0]) > EPS || fabs(Bin[i][0].second - wrmbin[i][i0+1]) > EPS
                || fabs(Bin[i][1].first - wrmbin[i][i0+2]) > EPS || fabs(Bin[i][1].second - wrmbin[i][i0+3]) > EPS
                || fabs(Bin[i][2].first - wrmbin[i][i0+4]) > EPS || fabs(Bin[i][2].second - wrmbin[i][i0+5]) > EPS) {
               logger.error["CheckWarmupConsistency"]
                  <<"Table of warmup values seems to be incompatible with steering file.\n"
                  <<"Found different binning for bin "<<i<<", steering: ["<<Bin[i][0].first<<","<<Bin[i][0].second<<",] ["<<Bin[i][1].first<<","<<Bin[i][1].second
                  <<"], warmup: ["<<wrmbin[i][i0]<<","<<wrmbin[i][i0+1]<<"] ["<<wrmbin[i][i0+2]<<","<<wrmbin[i][i0+3]<<"] ["<<wrmbin[i][i0+4]<<","<<wrmbin[i][i0+5]<<"].\n"
                  <<wrmuphelp
                  <<"Exiting."<<endl;
               ret = false;
               //exit(1);
            }
         }
         // KR: Check bin width only roughly
         double bwwrm = 0;
         if (NDim == 1) bwwrm = wrmbin[i][i0+2];
         else if (NDim == 2) bwwrm = wrmbin[i][i0+4];
         else if (NDim == 3) bwwrm = wrmbin[i][i0+6];
         if ( fabs(BinSize[i]) > DBL_MIN && (1 - BinSize[i]/bwwrm) > 1.e-4) {
            logger.warn["CheckWarmupConsistency"]
               <<"Table of warmup values seems to be incompatible with steering file.\n"
               <<"Found different bin size for bin "<<i<<". Steering: "<<BinSize[i]
               <<", warmup: "<<bwwrm<<".\n"
               <<wrmuphelp
               <<"Please check consistency!"<<endl<<endl;
            ret = false;
         }
      }
   }
   else {
      logger.warn["CheckWarmupConsistency"]
         <<"Warmup values do not contain information about the binning. Cannot check the consistency."<<endl;
   }
   return ret;
}



// ___________________________________________________________________________________________________
void fastNLOCreate::SetOrderOfAlphasOfCalculation(unsigned int ord) {
   logger.info["SetOrderOfAlphasOfCalculation"] << "Base order ord = " << ord << endl;
   //! set order of alpha_s of this calculation
   //! it must be: iLeadingOrder + iHigherOrder ;
   //! for instance: 3-jet-production in NLO = 4!
   // KR: Since the order of the run also determines the scale dependence,
   // KR: IScaleDep has to be set here and not in ReadCoefficientSpecificVariables()
   fastNLOCoeffAddBase* c = GetTheCoeffTable();
   const int IOrdInitial = fIOrd;
   fIOrd = ord;
   c->Npow = ord;
   c->IContrFlag2 = ord-GetLoOrder()+1;
   c->CtrbDescript.resize(1);
   string ss[5] = {"LO","NLO","NNLO","NNNLO","unknown"};
   int iContrb = ord-GetLoOrder() <4 ? ord-GetLoOrder() : 4;
   c->CtrbDescript[0] = ss[iContrb];
   if ((ord - GetLoOrder()) == 0) {
      c->IScaleDep = 0;
      logger.info["SetOrderOfAlphasOfCalculation"] << "LO scale dependence: MuR changeable; independent of MuF. IScaleDep = " << c->IScaleDep << endl;
      c->NSubproc               = fProcConsts.NSubProcessesLO;
      // fix different convention of flexible scale tables:
      if (fIsFlexibleScale &&  c->NSubproc==6 && fProcConsts.NSubProcessesNLO==7) c->NSubproc = fProcConsts.NSubProcessesNLO;   // 6 -> 7
      c->IPDFdef3               = fProcConsts.IPDFdef3LO;
      if ( c->IPDFdef2==0 ) {
         c->fPDFCoeff = fProcConsts.PDFCoeffLO;
         c->IPDFdef3               = c->NSubproc ;
      }
   } else {
      c->IScaleDep = 1;
      logger.info["SetOrderOfAlphasOfCalculation"] << "NLO scale dependence: Dependent on MuR and MuF. IScaleDep = " << c->IScaleDep << endl;
      if ((ord - GetLoOrder()) == 1) {
         c->NSubproc               = fProcConsts.NSubProcessesNLO;
         c->IPDFdef3               = fProcConsts.IPDFdef3NLO;
         if ( c->IPDFdef2==0 ) {
            c->fPDFCoeff = fProcConsts.PDFCoeffNLO;
            c->IPDFdef3  = c->NSubproc ;
         }
      } else if ((ord - GetLoOrder()) == 2) {
         c->NSubproc               = fProcConsts.NSubProcessesNNLO;
         c->IPDFdef3               = fProcConsts.IPDFdef3NNLO;
         if ( c->IPDFdef2==0 ) {
            c->fPDFCoeff = fProcConsts.PDFCoeffNNLO;
            c->IPDFdef3  = c->NSubproc ;
         }
      } else {
         logger.error["SetOrderOfAlphasOfCalculation"]<<"Unknown order of perturbation theory: order="<<ord-GetLoOrder()<<" (ord="<<ord<<",ILOord="<<ILOord<<"). Exiting."<<endl;
         exit(1);
      }
   }
   logger.info["SetOrderOfAlphasOfCalculation"] << "Using " << c->NSubproc <<
      " subprocesses and PDF flags: " << c->IPDFdef1 << ", " << c->IPDFdef2 << ", " << c->IPDFdef3 << "." << endl;

   // init array with counter processes (symmetric and asymmetric ones)
   if (c->NPDFPDG.size() == 2 && c->NPDFDim == 1) {
      fSymProc.resize(c->NSubproc);
      for (int p = 0 ; p<c->NSubproc ; p++) fSymProc[p]=p;

      for (unsigned int i = 0 ; i<fProcConsts.AsymmetricProcesses.size() ; i++) {
         if (fProcConsts.AsymmetricProcesses[i].first<c->NSubproc) {   // safety
            if (fProcConsts.AsymmetricProcesses[i].second >= GetNSubprocesses() || fSymProc[fProcConsts.AsymmetricProcesses[i].first] >= GetNSubprocesses()) {
               if (!(c->IPDFdef1==3&&c->IPDFdef2==1&&c->IPDFdef3==1))   // it is normal in pp->jets in LO
                  logger.warn["SetOrderOfAlphasOfCalculation"]<<"Subprocess "<<fSymProc[fProcConsts.AsymmetricProcesses[i].first]<<" is requested to be asymmetric with subprocess "<<fProcConsts.AsymmetricProcesses[i].second<<", but there are only "<<GetNSubprocesses()<<" subprocesses in this calculation. Ignoring call."<<endl;
            }
            else
               fSymProc[fProcConsts.AsymmetricProcesses[i].first] = fProcConsts.AsymmetricProcesses[i].second;
         }
      }

      //       vector<vector<int> > asym = INT_TAB_NS(AsymmetricProcesses);
      //       for ( unsigned int i = 0 ; i<asym.size() ; i ++ ) {
      //         if ( asym[i][0]<(int)fSymProc.size() ) { // safety
      //                    if ( asym[i][1] >= GetNSubprocesses() || fSymProc[asym[i][0]] >= GetNSubprocesses() ) {
      //                       logger.warn["AsymmetricProcesses"]<<"Subprocess "<<fSymProc[asym[i][0]]<<" is requested to be asymmetric with subprocess "<<asym[i][1]<<", but there are only "<<GetNSubprocesses()-1<<" subprocesses in this calculation. Ignoring call."<<endl;
      //                    }
      //                    else fSymProc[asym[i][0]] = asym[i][1];
      //                 }
      //       }

      //      for ( int p = 0 ; p<c->NSubproc ; p++ ) cout<<"p="<<p<<",\tfSymProc="<<fSymProc[p]<<endl;

   }

   // Scale factors have to be updated, since this may be either a lo or nlo run.
   if (!fIsFlexibleScale)  ReadScaleFactors();
   // NSubproc may have changed. We have to reinitialize the grids
   //const int nSFInitial = fScaleFac.size();

   // --- init statistics
   fastNLOTools::ResizeVector(GetTheCoeffTable()->fWgt.WgtObsSumW2, GetNSubprocesses(), GetNObsBin());
   fastNLOTools::ResizeVector(GetTheCoeffTable()->fWgt.SigObsSumW2, GetNSubprocesses(), GetNObsBin());
   fastNLOTools::ResizeVector(GetTheCoeffTable()->fWgt.SigObsSum  , GetNSubprocesses(), GetNObsBin());
   GetTheCoeffTable()->fWgt.WgtObsNumEv.resize(GetNSubprocesses());
   for ( int i =0 ; i<GetNSubprocesses() ; i++ ) {
      GetTheCoeffTable()->fWgt.WgtObsNumEv[i].resize(GetNObsBin());
   }
   //fastNLOTools::ResizeVector(GetTheCoeffTable()->fWgtObsNumEv, GetNSubprocesses(), GetNObsBin());


   if (IOrdInitial!= fIOrd && !fIsWarmup) {
      if (!fIsFlexibleScale) InitInterpolationKernels();
      InitGrids();
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::SetLoOrder(int LOOrd) {
   logger.debug["SetLoOrder"]<<endl;
   fastNLOTable::SetLoOrder(LOOrd);
   if (fIsFlexibleScale)
      ((fastNLOCoeffAddFlex*)GetTheCoeffTable())->fILOord = LOOrd;
}



// ___________________________________________________________________________________________________
void fastNLOCreate::InitVariablesInCoefficientTable() {
   logger.debug["InitVariablesInCoefficientTable"]<<endl;
   fastNLOCoeffAddBase* c = GetTheCoeffTable();
   c->IDataFlag = 0;            // No data, but theory
   c->IAddMultFlag = 0;         // additive contribution.
   c->IContrFlag1 = 1;          // fixed order: 1
   c->IContrFlag2 = 42;         // init with arbitrary number. to be specified later.
   c->IRef  = 0;                // it is not a reference calculation
   c->IScaleDep = 100;
   c->NScaleDep = 0;
   //c->NFragFunc               = 0;
   c->NFFDim            = 0;
   c->Nevt              = 0;
   c->SetIXsectUnits(12);       // it is often pb
}



// ___________________________________________________________________________________________________
void fastNLOCreate::ReadCoefficientSpecificVariables() {
   logger.debug["ReadCoefficientSpecificVariables"]<<endl;
   // todo: make it more user friendly
   // todo: include some sanity checks
   fastNLOCoeffAddBase* c = GetTheCoeffTable();

   // generator constants
   c->CodeDescript      = fGenConsts.GetCodeDescription();
   for (unsigned int i = 0 ; i<fProcConsts.GetProcessDescription().size() ; i++) {
      c->CodeDescript.push_back(fProcConsts.GetProcessDescription()[i]);
   }
   c->SetIXsectUnits(fGenConsts.UnitsOfCoefficients);

   // (some) process constants
   c->NPDFPDG.resize(fProcConsts.NPDF);
   if (c->NPDFPDG.size() >0) c->NPDFPDG[0] = fScenConsts.PDF1;   // from steering
   if (c->NPDFPDG.size() >1) c->NPDFPDG[1] = fScenConsts.PDF2;   // from steering
   c->NPDFDim           = fProcConsts.NPDFDim;
   c->IPDFdef1          = fProcConsts.IPDFdef1;
   c->IPDFdef2          = fProcConsts.IPDFdef2;
   c->IPDFdef3          = -1 ;          // safe initialization, is initialized in SetOrderOfAlphasOfCalculation(int); INT_NS(IPDFdef3);
   c->NSubproc          = -1;           // safe initialization, is initialized in SetOrderOfAlphasOfCalculation(int);
   c->NScaleDim         = 1;  // NEVER SET NScaleDim TO ANY OTHER VALUE THAN 1 !!!

   //IPDFdef3 = NSubproc == 7 ? 2 : 1;
   //printf("         Set IPDFdef3 = %d, consistent with %d subprocesses.\n",IPDFdef3,NSubproc);
   const int NScales = 1;
   c->Iscale.resize(NScales);
   c->Iscale[0] = 0;
   c->ScaleDescript.resize(NScales);

   if (fIsFlexibleScale) {
      c->NScaleDep              = 3; // temporarily. Until known if generator runs in LO, NLO or NNLO.
      c->ScaleDescript[0].resize(2);
      c->ScaleDescript[0][0] = fScenConsts.ScaleDescriptionScale1;
      c->ScaleDescript[0][1] = fScenConsts.ScaleDescriptionScale2;
   } else {
      // ---- those numbers are partly ambigously defined in v2.1 ---- //
      // proper code would be, but needs also adjustment in reading, writing and getter functions!
      // There are some misunderstandings here ...
      // KR: Reestablish v2.1 definitions
      int NScales   = 2; // NEVER SET NScales   TO ANY OTHER VALUE THAN 2 !!!
      int NScaleDim = 1; // NEVER SET NScaleDim TO ANY OTHER VALUE THAN 1 !!!
      c->NScales    = NScales;
      c->NScaleDim  = NScaleDim;
      c->Iscale.resize(NScales);
      c->Iscale[0]  = 0; // Use first scale definition for MuR; something else was never used for now
      c->Iscale[1]  = 0; // Use the same scale definition for MuF; something else was never used for now
      // KR: Since there is only one scale dimension, there is only one scale description!
      //     This description is for the different possibilities (dimensions) to choose for MuR, MuF etc.
      //     e.g. pT, sqrt(Q^2), M/2 and so on; this is NOT to describe the scales themselves which are
      //     always MuR, MuF!
      c->ScaleDescript.resize(NScaleDim);
      c->ScaleDescript[0].resize(NScaleDim);
      c->ScaleDescript[0][0] = fScenConsts.ScaleDescriptionScale1;
      c->SetNScaleDep(0);               // This is a fixed-scale table
   }
}



// ___________________________________________________________________________________________________
int fastNLOCreate::GetBin() {
   //! get bin number, using
   //! observables from Scenario

   const int idiff = GetNumDiffBin();
   // -------------------------------
   // check cache and return if available
   // KR: Bug fix from Enrico to avoid returning unitialized fObsBin for very first call
   if ((int)fLastScen._o.size() == idiff) {
      if (idiff == 1) {
         if (fLastScen._o[0] == fScenario._o[0]) return fObsBin;
      } else if (idiff == 2) {
         if (fLastScen._o[0] == fScenario._o[0] && fLastScen._o[1] == fScenario._o[1])  return fObsBin;
      } else if (idiff == 3) {
         if (fLastScen._o[0] == fScenario._o[0] && fLastScen._o[1] == fScenario._o[1] && fLastScen._o[2] == fScenario._o[2])  return fObsBin;
      } else {
         logger.error["GetBin"] << "More than triple-differential binning not yet implemented, aborted!" << endl;
      }
   }

   // -------------------------------
   // calc bin number and keep Observables
   if ( idiff == 1 ) {
      fObsBin = GetObsBinNumber(fScenario._o[0]);
   } else if ( idiff == 2 ) {
      fObsBin = GetObsBinNumber(fScenario._o[0],fScenario._o[1]);
   } else if ( idiff == 3 ) {
      fObsBin = GetObsBinNumber(fScenario._o[0],fScenario._o[1],fScenario._o[2]);
   } else {
      logger.error["GetBin"] << "More than triple-differential binning not yet implemented, aborted!" << endl;
   }
   fLastScen = fScenario;

   return fObsBin;
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillAllSubprocesses(const vector<vector<fnloEvent> >& events, const fnloScenario& scen) {
   //! fill all subprocessess for all scale variations (into a fixed-scale table)
   //! events is expected to be of the form:
   //!   events[nscalevar][nsubproc]

   const bool bFasterCode = true; // experimental developement: try to make code faster
   if (bFasterCode && !fIsWarmup && !fIsFlexibleScale) {
      // make filling code a little bit faster ... ~40%
      // if filling step "+=" is commented, then code is a factor ~6 faster
      fEvent = events[0][0];
      fScenario = scen;

      // KR: Also check for nan or inf in the faster code!!!
      if (!CheckWeightIsFinite()) return;

      const int ObsBin = (fScenario._iOB == -1) ? GetBin() : fScenario._iOB;
      if (ObsBin < 0) return;
      if (ObsBin >= (int)GetNObsBin()) return;
      fStats._nEvPS++;


      fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)GetTheCoeffTable();
      // do interpolation
      double xmin = GetTheCoeffTable()->GetNPDFDim() == 1 ? std::min(fEvent._x1,fEvent._x2) : fEvent._x1;
      double xmax = GetTheCoeffTable()->GetNPDFDim() == 1 ? std::max(fEvent._x1,fEvent._x2) : fEvent._x2;
      vector<pair<int,double> > nxlo = fKernX1[ObsBin]->GetNodeValues(xmin);
      vector<pair<int,double> > nxup = fKernX2[ObsBin]->GetNodeValues(xmax);

      if (fApplyPDFReweight) {
         fKernX1[ObsBin]->CheckX(xmin);
         fKernX2[ObsBin]->CheckX(xmax);
         ApplyPDFWeight(nxlo,xmin,fKernX1[ObsBin]->GetGridPtr()); // changes node values
         ApplyPDFWeight(nxup,xmax,fKernX2[ObsBin]->GetGridPtr()); // changes node values
      }


      fastNLO::v4d& st = c->SigmaTilde[ObsBin];
      for (unsigned int is = 0 ; is<events.size() ; is++) {
         double mu = fScenario._m1 * fScaleFac[is];
         const vector<pair<int,double> >& nmu  = fKernMuS[ObsBin][is]->GetNodeValues(mu);
         for (unsigned int m1 = 0 ; m1<nmu.size() ; m1++) {
            fastNLO::v2d& stm1 = st[is][nmu[m1].first];
            for (unsigned int p = 0 ; p<events[is].size() ; p++) {
               double wgt = events[is][p]._w * nmu[m1].second / BinSize[ObsBin];
               // .......................................................................................
               for (unsigned int x1 = 0 ; x1<nxlo.size() ; x1++) {
                  for (unsigned int x2 = 0 ; x2<nxup.size() ; x2++) {
                     //int p = fEvent._p;
                     int xminbin = nxlo[x1].first;
                     int xmaxbin = nxup[x2].first;
                     int proc = events[is][p]._p;
                     HalfMatrixCheck(fEvent._x1,fEvent._x2,xminbin,xmaxbin,proc);
                     int ixHM = GetXIndex(ObsBin,xminbin,xmaxbin);
                     stm1[ixHM][proc] += nxlo[x1].second * nxup[x2].second * wgt;//nmu[m1].second * wp[p];
                     //                for (unsigned int x1 = 0 ; x1<nxup.size() ; x1++) {
                     //                   for (unsigned int x2 = 0 ; x2<nxlo.size() ; x2++) {
                     //                      int xmaxbin = nxup[x1].first;
                     //                      int xminbin = nxlo[x2].first;
                     //                      int proc = events[is][p]._p;
                     //                      HalfMatrixCheck(fEvent._x1,fEvent._x2,xminbin,xmaxbin,proc);
                     //                      int ixHM = GetXIndex(ObsBin,xminbin,xmaxbin);
                     // //                   double ww = nxup[x1].second * nxlo[x2].second * wgt + stm1[ixHM][proc];
                     // //                   double& sti = stm1[ixHM][proc];
                     // //                   sti = ww;  // this assignment is time consuming ?!
                     //                      stm1[ixHM][proc] += nxup[x1].second * nxlo[x2].second * wgt;//nmu[m1].second * wp[p];
                     // //                   c->SigmaTilde[ObsBin][is][nmu[m1].first][ixHM][proc] += nxup[x1].second * nxlo[x2].second * wgt;//nmu[m1].second * wp[p];
                  }
               }
            }
         }
      }
   } else {
      for (unsigned int is = 0 ; is<events.size() ; is++) {   // all scalevars
         FillAllSubprocesses(events[is], scen, is);
      }
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillAllSubprocesses(const vector<fnloEvent>& events, const fnloScenario& scen, int scalevar) {
   //! fill a list of subprocesses into the fastNLO table

   if ((int)events.size() != GetNSubprocesses()) {
      logger.error["FillAllSubprocess"]<<"This table expects "<<GetNSubprocesses()<<" subprocesses, but only "<<events.size()<<" are provided. Exiting."<<endl;
      exit(1);
   }
   for (unsigned int p = 0 ; p<events.size() ; p++) {
      FillOneSubprocess(events[p],scen,scalevar);
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillOneSubprocess(const fnloEvent& event, const fnloScenario& scen, int scalevar) {
   fEvent = event;
   fScenario = scen;
   Fill(scalevar);
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillWeightCache(int scalevar) {
   //! fill event into cache
   if ( scalevar!= 0 ) {
      cout<<"Error! caching not implemented for scalevar tables."<<endl;
      exit(3); // also 'FlushCache has to be changed!'
   }
   //for ( auto& iEv : fWeightCache ) {
   for ( unsigned int iev = 0 ; iev<fWeightCache.size() ; iev++ ){
      if ( fWeightCache[iev].first._o[0] != fScenario._o[0] ) continue;
      if ( fWeightCache[iev].second._p   != fEvent._p  )    continue;
      if ( fWeightCache[iev].first._m1   != fScenario._m1 ) continue;
      if ( fWeightCache[iev].first._m2   != fScenario._m2 ) continue;
      if ( fWeightCache[iev].second._x1  != fEvent._x1 )    continue;
      if ( fWeightCache[iev].second._x2  != fEvent._x2 )    continue;

      if (     fWeightCache[iev].first._o[0] == fScenario._o[0]
               &&  fWeightCache[iev].second._p   == fEvent._p
               &&  fWeightCache[iev].first._m1   == fScenario._m1
               &&  fWeightCache[iev].first._m2   == fScenario._m2
               &&  fWeightCache[iev].second._x1  == fEvent._x1
               &&  fWeightCache[iev].second._x2  == fEvent._x2 ) {
         // cout<<"Found cached PS point !!"<<endl;
         // cout<<"cache_w: "<<fWeightCache[iev].second._w<<"\tevent_w:" <<fEvent._w<<"\t bin: "<<fScenario._iOB<<endl;
         fWeightCache[iev].second._w  += fEvent._w;
         fWeightCache[iev].second._wf += fEvent._wf;
         fWeightCache[iev].second._wr += fEvent._wr;
         fWeightCache[iev].second._wrr += fEvent._wrr;
         fWeightCache[iev].second._wff += fEvent._wff;
         fWeightCache[iev].second._wrf += fEvent._wrf;
         // cout<<"    sum w: "<<fWeightCache[iev].second._w<<endl;
         return;
      }
      else {
         cout<<"This is  a bug!"<<endl;
      }
   }
   // cout<<"Adding new weight to cache"<<endl;
   fWeightCache.push_back(make_pair(fScenario,fEvent));
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FlushCache() {
   // fill all events into from cache into table
   //for ( const auto& iev : fWeightCache ) {
   for ( unsigned int iev = 0 ; iev<fWeightCache.size() ; iev++ ){
      fScenario = fWeightCache[iev].first;
      fEvent = fWeightCache[iev].second;
      //cout<<"Filling weight: w="<<fEvent._w<<"\tbin="<<fScenario._iOB<<endl;
      FillContribution(0); // todo! not working if scalevar is != 0
   }
   fWeightCache.clear();
   fWeightCache.reserve(fCacheMax);
}



// ___________________________________________________________________________________________________
void fastNLOCreate::Fill(int scalevar) {
   //!
   //! Fill values, which are stored in 'Event' and 'Scenario' into fastNLO table.
   //!
   //logger.debug["Fill"]<<"Filling subprocess contributions into table."<<endl;

   // --- statistics and weights**2
   fStats._nProc++; //keep statistics
   int ObsBin = (fScenario._iOB == -1) ? GetBin() : fScenario._iOB;
   if (ObsBin < 0) return;
   //int ObsBin = GetBin();
   //fScenario._iOB = ObsBin;
   if ( scalevar==0 ) {
      fastNLOCoeffAddBase* c = GetTheCoeffTable();
      int p = fEvent._p;
      // --- event counts
      fStats._nEvPS++;
      c->fWgt.WgtObsNumEv[p][ObsBin]++;
      c->fWgt.WgtNumEv++;
      // --- w**2 counts
      double w2 = fEvent._w;
      if ( fIsFlexibleScale ) { // estimate a weight.
         double mu2 = fScenario._m1*fScenario._m1;
         if ( c->GetNPDF() == 1 ) //DIS
            mu2 = (fScenario._m1*fScenario._m1 + fScenario._m2*fScenario._m2 ) /2.;
         double lmu = log(mu2);
         w2 += lmu*fEvent._wr + lmu*fEvent._wf + lmu*lmu*fEvent._wrr + lmu*lmu*fEvent._wff + lmu*lmu*fEvent._wrf;
         if ( c->GetNPDF() == 1 ) //DIS
            w2 -= lmu*fEvent._wr +lmu*fEvent._wf;
      }
      w2 *= w2;
      // c->fWgtNevt ; //set externally by generator
      c->fWgt.WgtSumW2 += w2;
      c->fWgt.WgtObsSumW2[p][ObsBin] += w2;
      c->fWgt.SigObsSum  [p][ObsBin] += fEvent._sig;
      c->fWgt.SigObsSumW2[p][ObsBin] += fEvent._sig*fEvent._sig;
      c->fWgt.SigSum   += fEvent._sig;
      c->fWgt.SigSumW2 += fEvent._sig*fEvent._sig;
   }

   // sanity
   if ( fEvent._x1<0 || fEvent._x2<0 ) {
      logger.error["Fill"]<<"x-value is smaller than zero: x1="<<fEvent._x1<<", x2="<<fEvent._x2<<". Skipping event."<<endl;
      fEvent._x1=1;
      fEvent._x2=1;
      return ;
   }

   if (fIsWarmup) {
      if (scalevar==0 && ObsBin>=0 ) UpdateWarmupArrays();
      // else skip event
   }
   else if ( GetTheCoeffTable()->GetIRef() ) FillRefContribution(scalevar);
   else {
      if ( fIsFlexibleScale ) {
         if ( fCacheMax > 1 ) {
            FillWeightCache(scalevar);
            if ( (int)fWeightCache.size() >= fCacheMax )
               FlushCache();
         }
         else {
            FillContribution(scalevar);
         }
      }
      else {
         FillContribution(scalevar);
      }
   }
   fEvent.ResetButX();
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillContribution(int scalevar) {
   //! read information from 'Event' and 'Scenario'
   //! do the interpolation
   //! and fill into the tables.


   if (fEvent._n > 0) SetNumberOfEvents(fEvent._n);

   const int ObsBin = (fScenario._iOB == -1) ? GetBin() : fScenario._iOB;
   fastNLOCoeffAddBase* c = GetTheCoeffTable();
   int p = fEvent._p;

   // --- sanity
   if (ObsBin < 0) return;
   if (ObsBin >= (int)GetNObsBin()) return;
   if ( p<0 || p > c->GetNSubproc() ) {
      logger.error["FillContribution"]<<"Unknown process Id p = "<<p<<endl;
      exit(1);
   }

   // hier
   // cout<<"Fill ObsBin="<<ObsBin<<", w="<<fEvent._w<<", x1="<<fEvent._x1<<", x2="<<fEvent._x2
   //     <<", o0="<<fScenario._o[0]
   //     <<", m1="<<fScenario._m1
   //     <<", m2="<<fScenario._m2
   //     <<", OB="<<fScenario._iOB
   //     <<"\tproc="<<fEvent._p
   //     <<", wr="<<fEvent._wr
   //     <<", wf="<<fEvent._wf
   //     <<endl;
   // static double wsum = 0;
   // wsum+= fEvent._w; //hier
   // cout<<" * wSum = "<<wsum<<endl;


   // --- statistics and weights**2


   //logger.debug["FillContributionFixHHC"]<<"The process Id is p = "<<p<<endl;

   // ---- DIS ---- //
   if (c->GetNPDF() == 1 && fastNLOCoeffAddFlex::CheckCoeffConstants(c,true)) {
      // todo
      FillContributionFlexDIS((fastNLOCoeffAddFlex*)GetTheCoeffTable(),  ObsBin);
      //{logger.error["FillContribution"]<<"Don't know how to fill this table. Exiting."<<endl; exit(1); }
   } else if (c->GetNPDF() == 1 && fastNLOCoeffAddFix::CheckCoeffConstants(c,true)) {
      // todo
      FillContributionFixDIS((fastNLOCoeffAddFix*)GetTheCoeffTable(),  ObsBin, scalevar);
      //{logger.error["FillContribution"]<<"Don't know how to fill this table (DIS: fix-scale tables!). Exiting."<<endl; exit(1); }
   }
   // ---- pp/ppbar ---- //
   else if (c->GetNPDF() == 2 && fastNLOCoeffAddFlex::CheckCoeffConstants(c,true))
      FillContributionFlexHHC((fastNLOCoeffAddFlex*)GetTheCoeffTable(),  ObsBin);
   else if (c->GetNPDF() == 2 && fastNLOCoeffAddFix::CheckCoeffConstants(c,true))
      FillContributionFixHHC((fastNLOCoeffAddFix*)GetTheCoeffTable(),  ObsBin, scalevar);
   else {
      logger.error["FillContribution"]<<"Don't know how to fill this table. Exiting."<<endl;
      exit(1);
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillContributionFixHHC(fastNLOCoeffAddFix* c, int ObsBin, int scalevar) {
   //! read informatio from 'Event' and 'Scenario'
   //! do the interpolation
   //! and fill into the tables.
   //logger.debug["FillContributionFixHHC"]<<endl;

   if (fEvent._w == 0) return;   // nothing todo.

   // do interpolation
   double xmin = GetTheCoeffTable()->GetNPDFDim() == 1 ? std::min(fEvent._x1,fEvent._x2) : fEvent._x1;
   double xmax = GetTheCoeffTable()->GetNPDFDim() == 1 ? std::max(fEvent._x1,fEvent._x2) : fEvent._x2;

   vector<pair<int,double> > nxlo = fKernX1[ObsBin]->GetNodeValues(xmin);
   vector<pair<int,double> > nxup = fKernX2[ObsBin]->GetNodeValues(xmax);

   const double mu = fScenario._m1 * fScaleFac[scalevar];
   const vector<pair<int,double> >& nmu  = fKernMuS[ObsBin][scalevar]->GetNodeValues(mu);

   if (fApplyPDFReweight) {
      fKernX1[ObsBin]->CheckX(xmin);
      fKernX2[ObsBin]->CheckX(xmax);
      ApplyPDFWeight(nxlo,xmin,fKernX1[ObsBin]->GetGridPtr());
      ApplyPDFWeight(nxup,xmax,fKernX2[ObsBin]->GetGridPtr());
   }


   // fill grid
   if (!CheckWeightIsFinite()) return;
   double wgt = fEvent._w / BinSize[ObsBin];
   for (unsigned int x1 = 0 ; x1<nxlo.size() ; x1++) {
      for (unsigned int x2 = 0 ; x2<nxup.size() ; x2++) {
         int p = fEvent._p;
         int xminbin = nxlo[x1].first;
         int xmaxbin = nxup[x2].first;
         HalfMatrixCheck(fEvent._x1,fEvent._x2,xminbin,xmaxbin,p);
         //HalfMatrixCheck(xminbin,xmaxbin,p);
         int ixHM = GetXIndex(ObsBin,xminbin,xmaxbin);

         for (unsigned int m1 = 0 ; m1<nmu.size() ; m1++) {
            double w = wgt * nxlo[x1].second * nxup[x2].second * nmu[m1].second ;
            if (! std::isfinite(w)) {
               logger.error["FillContributionFixHHC"]<<"Weight w is not finite, w = " << w << "!"<<endl;
               logger.error["FillContributionFixHHC"]<<"This should have been captured before, aborting ..."<<endl;
               exit(1);
            }
            //              cout<<"   Fill * : i="<<ObsBin<<" svar="<<scalevar<<" imu="<<m1<<" ix="<<ixHM<<", im1="<<nmu[m1].first<<", p="<<p<<", w="<<nxup[x1].second * nxlo[x2].second * nmu[m1].second / BinSize[ObsBin]
            //          <<",\tfEvent._w="<<fEvent._w<<",\twx="<<nxup[x1].second * nxlo[x2].second<<",\tws="<<nmu[m1].second<<endl;
            c->SigmaTilde[ObsBin][scalevar][nmu[m1].first][ixHM][p] += w;
         }
      }
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillContributionFlexHHC(fastNLOCoeffAddFlex* c, int ObsBin) {
   //! read informatio from 'Event' and 'Scenario'
   //! do the interpolation
   //! and fill into the tables.
   //logger.debug["FillContributionFlexHHC"]<<endl;

   if (fEvent._w == 0 && fEvent._wf==0 && fEvent._wr==0 && fEvent._wrr==0 && fEvent._wff==0 && fEvent._wrf==0) return;   // nothing todo.

   // do interpolation
   //cout<<"try to interpol. ObsBin="<<ObsBin<<" ,x1="<<fEvent._x1<<", x2="<<fEvent._x2<<", mu1="<<Scenario._m1<<", mu2="<<Scenario._m2<<endl;
   double xmin = GetTheCoeffTable()->GetNPDFDim() == 1 ? std::min(fEvent._x1,fEvent._x2) : fEvent._x1;
   double xmax = GetTheCoeffTable()->GetNPDFDim() == 1 ? std::max(fEvent._x1,fEvent._x2) : fEvent._x2;
   vector<pair<int,double> > nxlo = fKernX1[ObsBin]->GetNodeValues(xmin);
   vector<pair<int,double> > nxup = fKernX2[ObsBin]->GetNodeValues(xmax);
   vector<pair<int,double> > nmu1 = fKernMu1[ObsBin]->GetNodeValues(fScenario._m1);
   vector<pair<int,double> > nmu2 = fKernMu2[ObsBin]->GetNodeValues(fScenario._m2);

   if (fApplyPDFReweight) {
      fKernX1[ObsBin]->CheckX(xmin);
      fKernX2[ObsBin]->CheckX(xmax);
      ApplyPDFWeight(nxlo,xmin,fKernX1[ObsBin]->GetGridPtr());
      ApplyPDFWeight(nxup,xmax,fKernX2[ObsBin]->GetGridPtr());
   }


   // fill grid
   if (!CheckWeightIsFinite()) return;
   for (unsigned int x1 = 0 ; x1<nxlo.size() ; x1++) {
      for (unsigned int x2 = 0 ; x2<nxup.size() ; x2++) {
         int xminbin = nxlo[x1].first;
         int xmaxbin = nxup[x2].first;
         int p = fEvent._p;
         HalfMatrixCheck(fEvent._x1,fEvent._x2,xminbin,xmaxbin,p);
         //HalfMatrixCheck(xminbin,xmaxbin,p);
         int ixHM = GetXIndex(ObsBin,xminbin,xmaxbin);

         for (unsigned int m1 = 0 ; m1<nmu1.size() ; m1++) {
            for (unsigned int mu2 = 0 ; mu2<nmu2.size() ; mu2++) {
               double wfnlo = nxlo[x1].second * nxup[x2].second * nmu1[m1].second * nmu2[mu2].second / BinSize[ObsBin];
               if (! std::isfinite(wfnlo)) {
                  logger.error["FillContributionFlexHHC"]<<"Weight wfnlo is not finite, wfnlo = " << wfnlo << "!"<<endl;
                  logger.error["FillContributionFlexHHC"]<<"This should have been captured before, aborting ..."<<endl;
                  fKernX1[ObsBin]->PrintGrid();
                  fKernX2[ObsBin]->PrintGrid();
                  fKernMu1[ObsBin]->PrintGrid();
                  fKernMu2[ObsBin]->PrintGrid();
                  cout<<"ix1="<<x1<<", ix2="<<x2<<", im1="<<m1<<", im2="<<mu2<<endl;
                  cout<<"x1="<<nxlo[x1].second<<", x1="<<x1<<", xval="<<xmin<<endl;
                  cout<<"x2="<<nxup[x2].second<<", x2="<<x2<<", xval="<<xmax<<endl;
                  cout<<"m1="<< nmu1[m1].second<<", m1="<<m1<<", mu1val="<<fScenario._m1<<endl;
                  cout<<"m2="<<nmu2[mu2].second<<", m2="<<mu2<<", mu2val="<<fScenario._m2<<endl;
                  exit(1);
               }
               if (fEvent._w  != 0) {
                  c->SigmaTildeMuIndep[ObsBin][ixHM][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._w  * wfnlo;
                  //wsum+= fEvent._w * wfnlo;
               }
               if (fEvent._wf != 0) {
                  //cout<<"   Fill F : ix="<<ixHM<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._wf  * wfnlo<<endl;
                  c->SigmaTildeMuFDep [ObsBin][ixHM][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wf * wfnlo;
               }
               if (fEvent._wr != 0) {
                  //cout<<"   Fill R : ix="<<ixHM<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._wr  * wfnlo<<endl;
                  c->SigmaTildeMuRDep [ObsBin][ixHM][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wr * wfnlo;
               }
               if (fEvent._wrr != 0) {
                  c->SigmaTildeMuRRDep [ObsBin][ixHM][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wrr * wfnlo;
               }
               if (fEvent._wff != 0) {
                  c->SigmaTildeMuFFDep [ObsBin][ixHM][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wff * wfnlo;
               }
               if (fEvent._wrf != 0) {
                  c->SigmaTildeMuRFDep [ObsBin][ixHM][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wrf * wfnlo;
               }
            }
         }
      }
   }

   //cout<<" * wSumW = "<<wsum<<endl;
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillContributionFlexDIS(fastNLOCoeffAddFlex* c, int ObsBin) {
   //! read information from 'Event' and 'Scenario'
   //! do the interpolation
   //! and fill into the tables.
   //logger.debug["FillContributionFlexDIS"]<<endl;

   if (fEvent._w == 0 && fEvent._wf==0 && fEvent._wr==0 && fEvent._wrr==0 && fEvent._wrf==0) return;   // nothing todo.

   // do interpolation
   //cout<<"try to interpol. ObsBin="<<ObsBin<<" ,x1="<<fEvent._x1<<", x2="<<fEvent._x2<<", mu1="<<Scenario._m1<<", mu2="<<Scenario._m2<<endl;

   // todo, just: 'x'
   double x = fEvent._x1;
   vector<pair<int,double> > nx = fKernX1[ObsBin]->GetNodeValues(x);
   vector<pair<int,double> > nmu1 = fKernMu1[ObsBin]->GetNodeValues(fScenario._m1);
   vector<pair<int,double> > nmu2 = fKernMu2[ObsBin]->GetNodeValues(fScenario._m2);


   if (fApplyPDFReweight) {
      fKernX1[ObsBin]->CheckX(x);
      ApplyPDFWeight(nx,x,fKernX1[ObsBin]->GetGridPtr());
   }

   // fill grid
   if (!CheckWeightIsFinite()) return;
   /*
   for (unsigned int ix = 0 ; ix<nx.size() ; ix++) {
      int p = fEvent._p;
      int xIdx = nx[ix].first;
      //HalfMatrixCheck(xminbin,xmaxbin,p);
      //int ixHM = GetXIndex(ObsBin,xminbin,xmaxbin);

      for (unsigned int m1 = 0 ; m1<nmu1.size() ; m1++) {
         for (unsigned int mu2 = 0 ; mu2<nmu2.size() ; mu2++) {
            double wfnlo = nx[ix].second * nmu1[m1].second * nmu2[mu2].second / BinSize[ObsBin];
            // if (! std::isfinite(wfnlo)) {
            //    logger.error["FillContributionFlexDIS"]<<"Weight wfnlo is not finite, wfnlo = " << wfnlo << "!"<<endl;
            //    logger.error["FillContributionFlexDIS"]<<"This should have been captured before, aborting ..."<<endl;
            //    fKernX1[ObsBin]->PrintGrid();
            //    fKernMu1[ObsBin]->PrintGrid();
            //    fKernMu2[ObsBin]->PrintGrid();
            //    cout<<"ix1="<<ix<<", im1="<<m1<<", im2="<<mu2<<endl;
            //    cout<<"x1="<<nx[ix].second<<", ix="<<ix<<", xval="<<x<<endl;
            //    cout<<"m1="<< nmu1[m1].second<<", m1="<<m1<<", mu1val="<<fScenario._m1<<endl;
            //    cout<<"m2="<<nmu2[mu2].second<<", m2="<<mu2<<", mu2val="<<fScenario._m2<<endl;
            //    exit(1);
            // }
            if (fEvent._w  != 0) {
               //cout<<"   Fill * : ix="<<ix<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._w  * wfnlo<<endl;
               c->SigmaTildeMuIndep[ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._w  * wfnlo;
            }
            if (fEvent._wf != 0) {
               //cout<<"   Fill F : ix="<<ix<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._wf  * wfnlo<<endl;
               c->SigmaTildeMuFDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wf * wfnlo;
            }
            if (fEvent._wr != 0) {
               //cout<<"   Fill R : ix="<<ix<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._wr  * wfnlo<<endl;
               c->SigmaTildeMuRDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wr * wfnlo;
            }
            if (fEvent._wrr != 0) {
               c->SigmaTildeMuRRDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wrr * wfnlo;
            }
            if (fEvent._wff != 0) {
               c->SigmaTildeMuFFDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wff * wfnlo;
            }
            if (fEvent._wrf != 0) {
               c->SigmaTildeMuRFDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wrf * wfnlo;
            }
         }
      }
   }
   */

   int p = fEvent._p;
   if (fEvent._w  != 0) {
      for (unsigned int ix = 0 ; ix<nx.size() ; ix++) {
         int xIdx = nx[ix].first;
         for (unsigned int m1 = 0 ; m1<nmu1.size() ; m1++) {
            for (unsigned int mu2 = 0 ; mu2<nmu2.size() ; mu2++) {
               double wfnlo = nx[ix].second * nmu1[m1].second * nmu2[mu2].second / BinSize[ObsBin];
               c->SigmaTildeMuIndep[ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._w  * wfnlo;
            }
         }
      }
   }
   if (fEvent._wf != 0) {
      for (unsigned int ix = 0 ; ix<nx.size() ; ix++) {
         int xIdx = nx[ix].first;
         for (unsigned int m1 = 0 ; m1<nmu1.size() ; m1++) {
            for (unsigned int mu2 = 0 ; mu2<nmu2.size() ; mu2++) {
               double wfnlo = nx[ix].second * nmu1[m1].second * nmu2[mu2].second / BinSize[ObsBin];
               c->SigmaTildeMuFDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wf * wfnlo;
            }
         }
      }
   }
   if (fEvent._wr != 0) {
      for (unsigned int ix = 0 ; ix<nx.size() ; ix++) {
         int xIdx = nx[ix].first;
         for (unsigned int m1 = 0 ; m1<nmu1.size() ; m1++) {
            for (unsigned int mu2 = 0 ; mu2<nmu2.size() ; mu2++) {
               double wfnlo = nx[ix].second * nmu1[m1].second * nmu2[mu2].second / BinSize[ObsBin];
               c->SigmaTildeMuRDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wr * wfnlo;
            }
         }
      }
   }
   if (fEvent._wrr != 0) {
      for (unsigned int ix = 0 ; ix<nx.size() ; ix++) {
         int xIdx = nx[ix].first;
         for (unsigned int m1 = 0 ; m1<nmu1.size() ; m1++) {
            for (unsigned int mu2 = 0 ; mu2<nmu2.size() ; mu2++) {
               double wfnlo = nx[ix].second * nmu1[m1].second * nmu2[mu2].second / BinSize[ObsBin];
               c->SigmaTildeMuRRDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wrr * wfnlo;
            }
         }
      }
   }
   if (fEvent._wff != 0) {
      for (unsigned int ix = 0 ; ix<nx.size() ; ix++) {
         int xIdx = nx[ix].first;
         for (unsigned int m1 = 0 ; m1<nmu1.size() ; m1++) {
            for (unsigned int mu2 = 0 ; mu2<nmu2.size() ; mu2++) {
               double wfnlo = nx[ix].second * nmu1[m1].second * nmu2[mu2].second / BinSize[ObsBin];
               c->SigmaTildeMuFFDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wff * wfnlo;
            }
         }
      }
   }
   if (fEvent._wrf != 0) {
      for (unsigned int ix = 0 ; ix<nx.size() ; ix++) {
         int xIdx = nx[ix].first;
         for (unsigned int m1 = 0 ; m1<nmu1.size() ; m1++) {
            for (unsigned int mu2 = 0 ; mu2<nmu2.size() ; mu2++) {
               double wfnlo = nx[ix].second * nmu1[m1].second * nmu2[mu2].second / BinSize[ObsBin];
               c->SigmaTildeMuRFDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wrf * wfnlo;
            }
         }
      }
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::FillContributionFixDIS(fastNLOCoeffAddFix* c, int ObsBin, int scalevar) {
   //! read information from 'Event' and 'Scenario'
   //! do the interpolation
   //! and fill into the tables.
   //logger.debug["FillContributionFixDIS"]<<endl;

   if (fEvent._w == 0 ) return;   // nothing todo.
   if ( scalevar >= (int)fScaleFac.size() ){
      logger.error["FillContributionFixDIS"]<<"Error! Scale variation scalevar="<<scalevar<<" requested"
                                            <<", but only "<<fScaleFac.size()<<" variations are initialized. Exiting."<<endl;
      exit(3);
   }

   // do interpolation
   //cout<<"try to interpol. ObsBin="<<ObsBin<<" ,x1="<<fEvent._x1<<", x2="<<fEvent._x2<<", mu1="<<fScenario._m1<<endl;

   // todo, just: 'x'
   double x = fEvent._x1;
   vector<pair<int,double> > nx = fKernX1[ObsBin]->GetNodeValues(x);
   vector<pair<int,double> > nmu1 = fKernMuS[ObsBin][scalevar]->GetNodeValues(fScenario._m1);


   if (fApplyPDFReweight) {
      fKernX1[ObsBin]->CheckX(x);
      ApplyPDFWeight(nx,x,fKernX1[ObsBin]->GetGridPtr());
   }

   // fill grid
   if (!CheckWeightIsFinite()) return;
   for (unsigned int ix = 0 ; ix<nx.size() ; ix++) {
      int p = fEvent._p;
      int xIdx = nx[ix].first;
      //HalfMatrixCheck(xminbin,xmaxbin,p);
      //int ixHM = GetXIndex(ObsBin,xminbin,xmaxbin);

      for (unsigned int m1 = 0 ; m1<nmu1.size() ; m1++) {
         double wfnlo = nx[ix].second * nmu1[m1].second  / BinSize[ObsBin];
         if (! std::isfinite(wfnlo)) {
            logger.error["FillContributionFixDIS"]<<"Weight wfnlo is not finite, wfnlo = " << wfnlo << "!"<<endl;
            logger.error["FillContributionFixDIS"]<<"This should have been captured before, aborting ..."<<endl;
            fKernX1[ObsBin]->PrintGrid();
            fKernMu1[ObsBin]->PrintGrid();
            cout<<"ix1="<<ix<<", im1="<<m1<<endl;
            cout<<"x1="<<nx[ix].second<<", ix="<<ix<<", xval="<<x<<endl;
            cout<<"m1="<< nmu1[m1].second<<", m1="<<m1<<", mu1val="<<fScenario._m1<<endl;
            exit(1);
         }
         if (fEvent._w  != 0) {
            //cout<<"   Fill * : ix="<<xIdx<<", im1="<<nmu1[m1].first<<", im2="<<", p="<<p<<", w="<<fEvent._w  * wfnlo<<endl;
            c->SigmaTilde[ObsBin][scalevar][nmu1[m1].first][xIdx][p]  += fEvent._w  * wfnlo;
         }
      }
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillRefContribution(int scalevar) {
   //! This is a reference table.
   //! Fill contribution as it would be a cross section

   if ( GetTheCoeffTable()->GetIRef()== 0 ) return; // error. this is not a ref-table

   // ObsBin
   const int ObsBin = (fScenario._iOB == -1) ? GetBin() : fScenario._iOB;
   double wgt = fEvent._w / BinSize[ObsBin];
   int p = fEvent._p;
   // todo....
   if (! std::isfinite(wgt)) {
      logger.error["FillContributionFixHHC"]<<"Weight w is not finite, w = " << wgt << "!"<<endl;
      logger.error["FillContributionFixHHC"]<<"This should have been captured before, aborting ..."<<endl;
      exit(1);
      }

   // fill reference table
   if ( fReader==NULL ){
      // weights are assumed to multiplied by PDF and alpha_s already
      // nothing to do
   }
   else {
      // weights are assumed to be in the same format as for 'default' fastNLO tables
      // and will be multiplied by PDF and alpha_s values
      double mur = fScenario._m1; //??;
      double muf = fScenario._m1; //??;
      double as  = fReader->EvolveAlphas(mur);
      int x1 = fEvent._x1;
      int x2 = fEvent._x2;
      vector<double> xfx1 =  fReader->GetXFX(x1, muf);
      vector<double> xfx2 =  fReader->GetXFX(x2, muf);
      bool IsPPBar = false;// todo!
      vector<double> pdflc = fReader->CalcPDFLinearCombination(GetTheCoeffTable(),xfx1,xfx2, IsPPBar);
      wgt *= as * pdflc[p]; // is this right?
   }

   if ( fIsFlexibleScale )
      ((fastNLOCoeffAddFlex*)GetTheCoeffTable())->SigmaTildeMuIndep[ObsBin][0][0][0][p]  += wgt;
   else
      ((fastNLOCoeffAddFix*)GetTheCoeffTable())->SigmaTilde[ObsBin][scalevar][0][0][p] += wgt;

}


// ___________________________________________________________________________________________________
inline void fastNLOCreate::HalfMatrixCheck(double x1, double x2, int& xminbin, int& xmaxbin, int& subproc) const {
   //! check if half-matrix notation
   //! if half-matrix notation, and xmin-node is larger than xmax-node
   //! exchange suprocesses according to fSymProc and adjust x-nodes.
   //!
   if (GetTheCoeffTable()->GetNPDFDim() == 1 ) {   // half-matrix notation (otherwise nothing todo)
      if ( x2>x1 ) subproc = fSymProc[subproc]; // to into correct half-matrix
      if (xminbin > xmaxbin) {
         //          if ( (int)fSymProc.size() != GetTheCoeffTable()->GetNSubproc() )
         //             logger.error["HalfMatrixCheck"]<<"Necessary array with symmetric processes for half-matrix notation not initialized."<<endl;

         //cout<<"exchange supbrpc. xminbin="<<xminbin<<", xmaxbin="<<xmaxbin<<", p="<<subproc<<", pAsym="<<fSymProc[subproc]<<endl;
         int di = xminbin - xmaxbin;
         xmaxbin += di;        // modify indicees
         xminbin -= di;
         subproc = fSymProc[subproc];           // exchange asymmetric process
      }
   }
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::CheckWeightIsFinite() {
   //! check if weights are finite
   if (! std::isfinite(fEvent._w)) {
      if (std::isnan(fEvent._w)) {
         logger.warn["CheckWeightIsFinite"]<<"(Scale-independent) weight is 'nan'!"<<endl;
      } else if (std::isinf(fEvent._w)) {
         logger.warn["CheckWeightIsFinite"]<<"(Scale-independent) weight is 'inf'!"<<endl;
      } else {
         logger.warn["CheckWeightIsFinite"]<<"(Scale-independent) weight is non-finite!"<<endl;
      }
      logger.warn["CheckWeightIsFinite"]<<"Contribution is skipped!"<<endl;
      return false;
   }
   if (! std::isfinite(fEvent._wf)) {
      if (std::isnan(fEvent._wf)) {
         logger.warn["CheckWeightIsFinite"]<<"Factorization scale dependent weight is 'nan'!"<<endl;
      } else if (std::isinf(fEvent._wf)) {
         logger.warn["CheckWeightIsFinite"]<<"Factorization scale dependent weight is 'inf'!"<<endl;
      } else {
         logger.warn["CheckWeightIsFinite"]<<"Factorization scale dependent weight is non-finite!"<<endl;
      }
      logger.warn["CheckWeightIsFinite"]<<"Contribution is skipped!"<<endl;
      return false;
   }
   if (! std::isfinite(fEvent._wr)) {
      if (std::isnan(fEvent._wr)) {
         logger.warn["CheckWeightIsFinite"]<<"Renormalization scale dependent weight is 'nan'!"<<endl;
      } else if (std::isinf(fEvent._wr)) {
         logger.warn["CheckWeightIsFinite"]<<"Renormalization scale dependent weight is 'inf'!"<<endl;
      } else {
         logger.warn["CheckWeightIsFinite"]<<"Renormalization scale dependent weight is non-finite!"<<endl;
      }
      logger.warn["CheckWeightIsFinite"]<<"Contribution is skipped!"<<endl;
      return false;
   }
   return true;
}


// ___________________________________________________________________________________________________
inline int fastNLOCreate::GetXIndex(const int& ObsBin,const int& x1bin,const int& x2bin) const {
   //! get index if 1 or two hadrons are involved
   //    switch (GetTheCoeffTable()->GetNPDFDim() ) {
   switch (((fastNLOCoeffAddBase*)fCoeff[0])->GetNPDFDim()) {
   case 1:
      // Daniel´s original
      return x1bin + (x2bin*(x2bin+1)/2);    // half matrix
      // Exchange x1 with x2
      //      return x2bin + (x1bin*(x1bin+1)/2);    // half matrix
   case 0:
      return x1bin; // linear
   case 2:
      return x1bin + x2bin * GetTheCoeffTable()->GetNxtot1(ObsBin); // full matrix
   default:
      return -1; // this will cause a crash :)
   }
};



// ___________________________________________________________________________________________________
int fastNLOCreate::GetNxmax(const vector<double>* xGrid1, const vector<double>* xGrid2) {
   switch (GetTheCoeffTable()->GetNPDFDim()) {
   case 0:
      return xGrid1->size();
   case 1:
      if (!xGrid2) logger.error["GetNxmax"]<<"Error. Second x-grid must be specified."<<endl;
      if (xGrid1->size() != xGrid2->size())logger.error["GetNxmax"]<<"Grid sizes in half-matrix notation must have equal size."<<endl;
      return ((int)pow((double)xGrid1->size(),2)+xGrid1->size())/2;
   case 2:
      if (!xGrid2) logger.error["GetNxmax"]<<"Error. Second x-grid must be specified."<<endl;
      return xGrid1->size()*xGrid2->size();
   default:
      return 0;
   }
};


// ___________________________________________________________________________________________________
inline void fastNLOCreate::ApplyPDFWeight(vector<pair<int,double> >& nodes, const double x, const vector<double>* grid) const {
   //    double pdfwgtmax = PDFwgt(xmax);
   //    for( int i1 = 0; i1 < 4; i1++) {
   //       if ((nxmaxf-1+i1) >= 0 && (nxmaxf-1+i1) < Nxtot1[ObsBin] ) {
   //          cefmax[i1] *= pdfwgtmax/PDFwgt(XNode1[ObsBin][nxmaxf-1+i1]);
   //       }
   //    }
   double wgtx = CalcPDFReweight(x);
   for (unsigned int in = 0; in < nodes.size(); in++) {
      double wgtn = CalcPDFReweight(grid->at(nodes[in].first));
      if (wgtn==0) {logger.error["ApplyPDFWeight"]<<"Cannot divide by 0."<<endl; exit(1);}
      //cout<<"in="<<in<<" wgtx="<<wgtx<<", x="<<x<<", wgtn="<<wgtn<<", nod="<<grid->at(nodes[in].first)<<endl;
      nodes[in].second *= wgtx/wgtn;
   }
}


// ___________________________________________________________________________________________________
inline double fastNLOCreate::CalcPDFReweight(double x) const {
   if (x<=0) { logger.error["CalcPDFReweight"]<<"Cannot calculate sqrt of negative numbers or divide by zero. x="<<x<<endl; exit(1);}
   double w=(1.-0.99*x)/sqrt(x);
   return w*w*w;
}


// ___________________________________________________________________________________________________
void fastNLOCreate::NormalizeCoefficients(double wgt) {
   //! Set number of events to wgt (default=1) and weight coefficients in sigmatilde
   //! accordingly
   //! This means, that the information about the
   //! number of events is essentially lost (now remained stored in fWgt)
   if ( fWeightCache.size() )  FlushCache();
   GetTheCoeffTable()->NormalizeCoefficients(wgt);
   fStats._nEv=wgt;
   //    double nev = GetTheCoeffTable()->GetNevt(0,0);
   //    MultiplyCoefficientsByConstant(1./nev);
   //    SetNumberOfEvents(1.);
}
// ___________________________________________________________________________________________________
void fastNLOCreate::NormalizeCoefficients(const std::vector<std::vector<double> >& wgtProcBin) {
   //! Set number of events to wgtProcBin[iProc][iBin]
   //! sigmatilde is weighted accordingly.
   if ( fWeightCache.size() )  FlushCache();
   GetTheCoeffTable()->NormalizeCoefficients(wgtProcBin);
   fStats._nEv=-1;
}


// ___________________________________________________________________________________________________
void fastNLOCreate::MultiplyCoefficientsByBinSize() {
   //! Multiply all coefficients by binsize
   if ( fWeightCache.size() )  FlushCache();
   if (fIsFlexibleScale) {
      fastNLOCoeffAddFlex* c = (fastNLOCoeffAddFlex*)GetTheCoeffTable();
      for (unsigned int i=0; i<GetNObsBin(); i++) {
         int nxmax = c->GetNxmax(i);
         for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
            for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
               for (int x=0; x<nxmax; x++) {
                  for (int n=0; n<c->GetNSubproc(); n++) {
                     c->SigmaTildeMuIndep[i][x][jS1][kS2][n] *= BinSize[i];
                     if (c->GetNScaleDep() >= 5) {
                        c->SigmaTildeMuFDep [i][x][jS1][kS2][n] *= BinSize[i];
                        c->SigmaTildeMuRDep [i][x][jS1][kS2][n] *= BinSize[i];
                        if (c->GetNScaleDep() >= 6) {
                           c->SigmaTildeMuRRDep [i][x][jS1][kS2][n] *= BinSize[i];
                        }
                        if (c->GetNScaleDep() >= 7) {
                           c->SigmaTildeMuFFDep [i][x][jS1][kS2][n] *= BinSize[i];
                           c->SigmaTildeMuRFDep [i][x][jS1][kS2][n] *= BinSize[i];
                        }
                     }
                  }
               }
            }
         }
      }
   } else {
      fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)GetTheCoeffTable();
      for (unsigned int i=0; i<GetNObsBin(); i++) {
         for (unsigned int s=0 ; s<c->SigmaTilde[i].size() ; s++) {
            for (unsigned int x=0 ; x<c->SigmaTilde[i][s].size() ; x++) {
               for (unsigned int l=0 ; l<c->SigmaTilde[i][s][x].size() ; l++) {
                  for (unsigned int m=0 ; m<c->SigmaTilde[i][s][x][m].size() ; m++) {
                     c->SigmaTilde[i][s][x][l][m] *= BinSize[i];
                  }
               }
            }
         }
      }
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::DivideCoefficientsByBinSize() {
   //! Divide all coefficients by binsize
   if ( fWeightCache.size() )  FlushCache();
   if (fIsFlexibleScale) {
      fastNLOCoeffAddFlex* c = (fastNLOCoeffAddFlex*)GetTheCoeffTable();
      for (unsigned int i=0; i<c->SigmaTildeMuIndep.size(); i++) {
         int nxmax = c->GetNxmax(i);
         for (unsigned int jS1=0; jS1<c->GetNScaleNode1(i); jS1++) {
            for (unsigned int kS2=0; kS2<c->GetNScaleNode2(i); kS2++) {
               for (int x=0; x<nxmax; x++) {
                  for (int n=0; n<c->GetNSubproc(); n++) {
                     c->SigmaTildeMuIndep[i][x][jS1][kS2][n] /= BinSize[i];
                     //if ( c->GetNScaleDep() >= 5 ) {
                     if (!c->SigmaTildeMuFDep.empty()) {
                        c->SigmaTildeMuFDep [i][x][jS1][kS2][n] /= BinSize[i];
                        c->SigmaTildeMuRDep [i][x][jS1][kS2][n] /= BinSize[i];
                        //if ( c->GetNScaleDep() >= 6 ) {
                        if (!c->SigmaTildeMuRRDep.empty()) {
                           c->SigmaTildeMuRRDep [i][x][jS1][kS2][n] /= BinSize[i];
                        }
                        if (!c->SigmaTildeMuFFDep.empty()) {
                           c->SigmaTildeMuFFDep [i][x][jS1][kS2][n] /= BinSize[i];
                           c->SigmaTildeMuRFDep [i][x][jS1][kS2][n] /= BinSize[i];
                        }
                     }
                  }
               }
            }
         }
      }
   } else {
      fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)GetTheCoeffTable();
      for (unsigned int i=0; i<c->SigmaTilde.size(); i++) {
         for (unsigned int s=0 ; s<c->SigmaTilde[i].size() ; s++) {
            for (unsigned int x=0 ; x<c->SigmaTilde[i][s].size() ; x++) {
               for (unsigned int l=0 ; l<c->SigmaTilde[i][s][x].size() ; l++) {
                  for (unsigned int m=0 ; m<c->SigmaTilde[i][s][x][m].size() ; m++) {
                     c->SigmaTilde[i][s][x][l][m] /= BinSize[i];
                  }
               }
            }
         }
      }
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::MultiplyCoefficientsByConstant(double coef) {
   //! Divide all coefficients by binsize
   if ( fWeightCache.size() )  FlushCache();
   if (fIsFlexibleScale) {
      fastNLOCoeffAddFlex* c = (fastNLOCoeffAddFlex*)GetTheCoeffTable();
      c->MultiplyCoefficientsByConstant(coef);
   } else {
      fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)GetTheCoeffTable();
      c->MultiplyCoefficientsByConstant(coef);
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::UpdateWarmupArrays() {
   //! Update the warmup-arrays fWMu1, fWx und fWMu2
   if (fWx.empty()) InitWarmupArrays();

   const int ObsBin = GetBin();
   if ( ObsBin < 0 ) return;
   logger.debug["UpdateWarmupArrays"]<<"ObsBin="<<ObsBin<<"\tmu1="<<fScenario._m1<<"\tmu2="<<fScenario._m2<<"\tx1="<<fEvent._x1<<"\tx2="<<fEvent._x2<<endl;

   fWMu1[ObsBin].first       = std::min(fScenario._m1,fWMu1[ObsBin].first) ;
   fWMu1[ObsBin].second      = std::max(fScenario._m1,fWMu1[ObsBin].second) ;
   if (GetTheCoeffTable()->IPDFdef1 == 3) { // pp/ppbar
      fWx[ObsBin].first      = std::min(std::min(fEvent._x1,fEvent._x2),fWx[ObsBin].first) ;
      fWx[ObsBin].second     = std::max(std::max(fEvent._x1,fEvent._x2),fWx[ObsBin].second) ;
   } else if (GetTheCoeffTable()->IPDFdef1 == 2) {  // DIS
      fWx[ObsBin].first      = std::min(fEvent._x1,fWx[ObsBin].first) ;
      fWx[ObsBin].second     = std::max(fEvent._x1,fWx[ObsBin].second) ;
   } else
      logger.error["UpdateWarmupArrays"]<<"nothing reasonable implemented yet: IPDFdef1="<<GetTheCoeffTable()->IPDFdef1<<endl;
   if (fIsFlexibleScale) {
      fWMu2[ObsBin].first    = std::min(fScenario._m2,fWMu2[ObsBin].first) ;
      fWMu2[ObsBin].second   = std::max(fScenario._m2,fWMu2[ObsBin].second) ;
   }
   if ( fWx[ObsBin].first < 0 ) {
      logger.error["UpdateWarmupArrays"]<<"x-value is smaller than 0. Exiting."<<endl;
      exit(4);
   }

}


// ___________________________________________________________________________________________________
void fastNLOCreate::InitWarmupArrays() {
   logger.debug["InitWarmupArrays"]<<endl;
   //! initialize arrays to store and determine warm-up values
   //! including copy for later rounding and write out
   //! initialize with reasonable values
   fWMu1.resize(GetNObsBin());
   fWMu2.resize(GetNObsBin());
   fWx.resize(GetNObsBin());
   fWMu1Rnd.resize(GetNObsBin());
   fWMu2Rnd.resize(GetNObsBin());
   fWxRnd.resize(GetNObsBin());
   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      fWMu1[i].first     =  DBL_MAX;
      fWMu1[i].second    = -DBL_MAX;
      fWMu2[i].first     =  DBL_MAX;
      fWMu2[i].second    = -DBL_MAX;
      fWx[i].first       =  DBL_MAX;
      fWx[i].second      = -DBL_MAX;
      fWMu1Rnd[i].first  =  DBL_MAX;
      fWMu1Rnd[i].second = -DBL_MAX;
      fWMu2Rnd[i].first  =  DBL_MAX;
      fWMu2Rnd[i].second = -DBL_MAX;
      fWxRnd[i].first    =  DBL_MAX;
      fWxRnd[i].second   = -DBL_MAX;
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::WriteTable() {
   //! Write fastNLO file to disk
   //if ( GetTheCoeffTable()->GetNevt(0,0) <= 0 ) {
   // flush cache with remaining events
   if ( fWeightCache.size() )  FlushCache();
   if (GetTheCoeffTable()->Nevt <= 0) {
      logger.warn["WriteTable"]<<"Number of events seems to be not filled. Please use SetNumberOfEvents(int) before writing table."<<endl;
      exit(1);
   }
   fStats.PrintStats();
   if (fIsWarmup) {
      logger.info["WriteTable"]<<"Writing warmup table instead of coefficient table."<<endl;
      if ( fWx.empty() ) {
         logger.error["WriteTable"]<<"Warmup values seem not to be initialized correctly. Maybe forgot to call 'Fill()'?"<<endl;
         exit(1);
      }
      // round warmup values and try to guess bin boundaries
      AdjustWarmupValues();
      // write table to disk
      WriteWarmupTable();
   } else {
      if (ffilename == "") {
         logger.error["WriteTable"]<<"No filename given."<<endl;
         exit(1);
      }
      bool check = CheckProcConsts();
      if (!check) {
         logger.error["fastNLOCreate"]<<"Process constants not properly set! Please check warning messages and complement your steering."<<endl;
         exit(1);
      }
      // Number of events must be counted correctly.
      // I.e. the counting should be performed by the generator.
      // ->Divide by BinSize
      fastNLOTable::WriteTable();
      // ->Multiply by BinSize
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::WriteTable(string filename) {
   SetFilename(filename);
   fastNLOTable::WriteTable();
}


// ___________________________________________________________________________________________________
void fastNLOCreate::WriteWarmupTable() {
   string tempfn = ffilename;
   string warmupfile = GetWarmupTableFilename();
   logger.info["WriteWarmupTable"]<<"Writing warmup table to: "<<warmupfile<<endl;
   SetFilename(warmupfile);

   // open stream;
   std::ostream* table = OpenFileWrite();
   // write to disk
   OutWarmup(*table);
   // close file
   // table->close();
   delete table;
   // reset filename
   SetFilename(tempfn);
}



// ___________________________________________________________________________________________________
void fastNLOCreate::PrintWarmupValues() {
   OutWarmup(std::cout);
}


// ___________________________________________________________________________________________________
void fastNLOCreate::OutWarmup(ostream& strm) {
   if ( fWeightCache.size() )  FlushCache();
   if (fWxRnd.empty()) {
      logger.warn["OutWarmup"]<<"Warmup arrays not initialized. Did you forgot to fill values?"<<endl;
      //       logger.warn["OutWarmup"]<<"  Continuting, but writing unreasonalby large/small values as warmup values..."<<endl;
      //       InitWarmupArrays();
      logger.error["OutWarmup"]<<" Do not write out unreasonable warmup table. Exiting."<<endl;
      exit(1);
   }

   strm<<"# --- Use emacs in sh mode -*-sh-*- #"<<endl;
   strm<<"# This is a automatically generated file by fastNLO and holds the values of the warmup run. "<<endl;
   strm<<"# The values are valid for the scenario "<<GetScenName() << endl;
   strm<<"# and if calculated with the steerfile: "<< fSteerfile <<endl;
   strm<<"# but only if no serious changes have been performed since its creation."<<endl;
   strm<<"# "<<endl;
   strm<<"# Delete this file, if you want fastNLO to calculate a new one."<<endl;
   strm<<"# "<<endl;
   strm<<"# This file has been calculated using: "<<endl;
   strm<<"#      "<<GetTheCoeffTable()->Nevt<<" contributions."<<endl;
   strm<<"#      "<<GetTheCoeffTable()->fWgt.WgtNumEv<<" entries."<<endl;
   strm<<"#   ( Mind: contributions != events. And contributions are not necessarily in phase space region."<<endl;
   strm<<"# Please check by eye for reasonability of the values."<<endl;
   strm<<"# Number of events per bin are listed below."<<endl;
   strm<<" " <<endl;

   // write variables of warmup run
   strm<<"Warmup.OrderInAlphasOfWarmupRunWas\t"<<  fIOrd <<endl;
   //strm<<"Warmup.CheckScaleLimitsAgainstBins\t"<<(BOOL_NS(CheckScaleLimitsAgainstBins,fSteerfile)?"true":"false")<<endl;
   strm<<"Warmup.CheckScaleLimitsAgainstBins\t"<< (fScenConsts.CheckScaleLimitsAgainstBins ?"true":"false")<<endl;
   strm<<"Warmup.ScaleDescriptionScale1     \t\""<< GetTheCoeffTable()->ScaleDescript[0][0]<<"\""<<endl;
   if (fIsFlexibleScale)
      strm<<"Warmup.ScaleDescriptionScale2     \t\""<< GetTheCoeffTable()->ScaleDescript[0][1]<<"\"" <<endl;
   strm<<"Warmup.DifferentialDimension      \t"<< NDim <<endl;
   strm<<"Warmup.DimensionLabels {\n  ";
   for (unsigned int i = 0 ; i < NDim; i ++) strm<<"\""<<DimLabel[i]<<"\"  ";
   strm<<"\n} "<<endl;

   strm<<"Warmup.DimensionIsDifferential {\n  ";
   for (unsigned int i = 0 ; i < NDim; i ++) strm<<"\""<<IDiffBin[i]<<"\"  ";
   strm<<"\n} "<<endl;
   strm<<endl;


   // --- sanity test, warnings and dummy values (if needed)
   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      int nEvBin = 0;
      for ( int ip=0 ; ip<GetNSubprocesses() ; ip++ ) {
         if (fWxRnd[i].first < 1.e-6) {
            logger.warn["OutWarmup"]<<"The xmin value in bin "<<i
                                    <<" seems to be unreasonably low (xmin="<<fWxRnd[i].first
                                    <<"). Taking xmin=1.e-6 instead."<<endl;
            fWxRnd[i].first=1.e-6;
         }
         nEvBin += GetTheCoeffTable()->fWgt.WgtObsNumEv[ip][i];
      }
      if ( nEvBin == 0 ) {
         logger.error["OutWarmup"]<<"No events were counted in bin "<<i<<". Thus no sensible warmup table can be written."<<endl;
         fWxRnd[i].first=1.e-6;
         fWxRnd[i].second=1;
         fWMu1[i].first=1;
         fWMu1[i].second=5000;
         fWMu2[i].first=1;
         fWMu2[i].second=5000;
         logger.error["OutWarmup"]<<"Continueing and taking sensible dummy values. Do not use these for production runs !!"<<endl;
      }
      else if ( nEvBin < 10 ) {
         logger.warn["OutWarmup"]<<"Too little events (n="<<nEvBin<<") were counted in bin "<<i<<". Thus no sensible warmup table can be written."<<endl;
      }
      else if ( nEvBin < 100 ) {
         logger.warn["OutWarmup"]<<"Quite few events (n="<<nEvBin<<") were counted in bin "<<i<<". Thus no sensible warmup table can be written."<<endl;
      }
   }

   // ---- write readable table
   char buf[4000];
   char buf2[4000];
   char buf3[4000];
   strm<<"Warmup.Values {{"<<endl;
   if (fIsFlexibleScale) {
      // table header
      sprintf(buf,"   ObsBin  %9s  %9s  %14s  %14s  %14s  %14s",
              "x_min","x_max",
              GetWarmupHeader(0,"min").c_str(), GetWarmupHeader(0,"max").c_str(),
              GetWarmupHeader(1,"min").c_str(), GetWarmupHeader(1,"max").c_str());
      strm<<buf<<endl;
      // table values
      for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
         // --- write x-values
         sprintf(buf,"   %4d    %9.1e  %9.2e",   i, fWxRnd[i].first, fWxRnd[i].second );

         // 1. Are warmup-values identical to bin-boundaries ?
         int ident1 = CheckWarmupValuesIdenticalWithBinGrid(fWMu1);
         int ident2 = CheckWarmupValuesIdenticalWithBinGrid(fWMu2);

         fWMu1Rnd = fWMu1;
         fWMu2Rnd = fWMu2;
         // --- write mu1
         if ( ident1 != -1 ) {
            // mu-value is identical with binning
            // -> write out many digits
            sprintf(buf2,"  %14.6f  %14.6f",fWMu1[i].first,fWMu1[i].second);
         }
         else if  ( fWMu1[i].second==0 || fabs(fWMu1[i].first/fWMu1[i].second-1) < 1.e-5) {
            // mu-value is a fixed number
            sprintf(buf2,"  %14.8f  %14.8f",fWMu1[i].first,fWMu1[i].first);
         }
         else {
            // scale values are floating points.
            // -> round them up/down a bit
            // extent range by 2%
            fWMu1Rnd[i].first  = fWMu1[i].first - 0.02*fabs(fWMu1[i].first);
            fWMu1Rnd[i].second = fWMu1[i].second + 0.02*fabs(fWMu1[i].second);
            RoundValues(fWMu1Rnd,fWarmupNDigitMu1); // digit here should be identical to output in outwarmup

            // if values are very close to an integer, it is likely that this is by purpose
            if ( fWMu1[i].first >= 1 && fabs(fWMu1[i].first - round(fWMu1[i].first)) < 1e-3 )
               fWMu1Rnd[i].first = round(fWMu1[i].first);
            if ( fWMu1[i].second >= 1 && fabs(fWMu1[i].second - round(fWMu1[i].second)) < 1e-3)
               fWMu1Rnd[i].second = round(fWMu1[i].second);

            string format = "  %14."+to_string(fWarmupNDigitMu1)+"f  %14."+to_string(fWarmupNDigitMu1)+"f";
            sprintf(buf2,format.c_str(),fWMu1Rnd[i].first,fWMu1Rnd[i].second);
         }

         // --- write mu2
         if ( ident2 != -1 ) {
            // mu-value is identical with binning
            // -> write out many digits
            sprintf(buf3,"  %14.6f  %14.6f",fWMu2[i].first,fWMu2[i].second);
         }
         else if  ( fWMu2[i].second==0 || fabs(fWMu2[i].first/fWMu2[i].second-1) < 1.e-5) {
            // mu-value is a fixed number
            sprintf(buf3,"  %14.8f  %14.8f",fWMu2[i].first,fWMu2[i].first);
         }
         else {
            // scale values are floating points.
            // -> round them up/down a bit
            fWMu2Rnd[i].first  = fWMu2[i].first - 0.06*fabs(fWMu2[i].first);
            fWMu2Rnd[i].second = fWMu2[i].second + 0.06*fabs(fWMu2[i].second);
            RoundValues(fWMu2Rnd,fWarmupNDigitMu2); // digit here should be identical to output in outwarmup

            // if values are very close to an integer, it is likely that this is by purpose
            if ( fWMu2[i].first >= 1 && fabs(fWMu2[i].first - round(fWMu2[i].first)) < 1e-3 )
               fWMu2Rnd[i].first = round(fWMu2[i].first);
            if ( fWMu2[i].second >= 1 && fabs(fWMu2[i].second - round(fWMu2[i].second)) < 1e-3)
               fWMu2Rnd[i].second = round(fWMu2[i].second);

            string format = "  %14."+to_string(fWarmupNDigitMu2)+"f  %14."+to_string(fWarmupNDigitMu2)+"f";
            sprintf(buf3,format.c_str(),fWMu2Rnd[i].first,fWMu2Rnd[i].second);
         }

         // if  ( fWMu2[i].second!=0 && fabs(fWMu2[i].first/fWMu2[i].second-1) > 1.e-3)
         //    sprintf(buf3,"  %14.2g  %14.2g",fWMu2Rnd[i].first,fWMu2Rnd[i].second);
         // else
         //    sprintf(buf3,"  %14.8f  %14.8f",fWMu2[i].first,fWMu2[i].second);

         // printf(" org  %e     %e     %e      %e\n",fWMu1[i].first,fWMu1[i].second,fWMu2[i].first,fWMu2[i].second);
         // printf(" rnd  %e     %e     %e      %e\n",fWMu1Rnd[i].first,fWMu1Rnd[i].second,fWMu2Rnd[i].first,fWMu2Rnd[i].second);
         // cout<<buf<<buf2<<buf3<<endl;
         // printf("\n");

         strm<<buf<<buf2<<buf3<<endl;

      }
   } else {
      // is ScaleDescript available?
      if (GetTheCoeffTable()->ScaleDescript[0].empty()) { logger.error["OutWarmup"]<<"Scale description is empty. but needed. Probably this has to be implemented."<<endl; exit(1);};
      // table header
      sprintf(buf,"   ObsBin   %9s  %9s  %16s  %16s",
              "x_min","x_max", GetWarmupHeader(0,"min").c_str(), GetWarmupHeader(0,"max").c_str());
      strm<<buf<<endl;


         // 1. Are warmup-values identical to bin-boundaries ?
      int ident1 = CheckWarmupValuesIdenticalWithBinGrid(fWMu1);
      fWMu1Rnd = fWMu1;


      // table values
      for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
         if (fWxRnd[i].first < 1.e-6) {
            logger.warn["OutWarmup"]<<"The xmin value in bin "<<i<<" seems to be unreasonably low (xmin="<<fWxRnd[i].first<<"). Taking xmin=1.e-6 instead."<<endl;
            fWxRnd[i].first=1.e-6;
         }
         // sprintf(buf,"   %4d     %9.1e  %9.2e  %16.1g  %16.1g",
         //      i,fWxRnd[i].first, fWxRnd[i].second, fWMu1Rnd[i].first, fWMu1Rnd[i].second);
          sprintf(buf,"   %4d     %9.1e  %9.2e", i,fWxRnd[i].first, fWxRnd[i].second);

          // --- old code --- //
          // if  ( fWMu1[i].second!=0 && fabs(fWMu1[i].first/fWMu1[i].second-1) > 1.e-3)
          //    sprintf(buf2,"  %16.2g  %16.2g", fWMu1Rnd[i].first, fWMu1Rnd[i].second);
          // else
          //    sprintf(buf2,"  %16.8f  %16.8f", fWMu1[i].first, fWMu1[i].second);
          // --- --- ---- --- //
          // --- write mu1
          if ( ident1 != -1 ) {
             // mu-value is identical with binning
             // -> write out many digits
             sprintf(buf2,"  %14.6f  %14.6f",fWMu1[i].first,fWMu1[i].second);
          }
          else if  ( fWMu1[i].second==0 || fabs(fWMu1[i].first/fWMu1[i].second-1) < 1.e-5) {
             // mu-value is a fixed number
             sprintf(buf2,"  %14.8f  %14.8f",fWMu1[i].first,fWMu1[i].first);
          }
          else {
             // scale values are floating points.
             // -> round them up/down a bit
             // extent range by 2%
             fWMu1Rnd[i].first  = fWMu1[i].first - 0.02*fabs(fWMu1[i].first);
             fWMu1Rnd[i].second = fWMu1[i].second + 0.02*fabs(fWMu1[i].second);
             RoundValues(fWMu1Rnd,fWarmupNDigitMu1); // digit here should be identical to output in outwarmup

             // if values are very close to an integer, it is likely that this is by purpose
             if ( fWMu1[i].first >= 1 && fabs(fWMu1[i].first - round(fWMu1[i].first)) < 1e-3 )
                fWMu1Rnd[i].first = round(fWMu1[i].first);
             if ( fWMu1[i].second >= 1 && fabs(fWMu1[i].second - round(fWMu1[i].second)) < 1e-3)
                fWMu1Rnd[i].second = round(fWMu1[i].second);

             string format = "  %14."+to_string(fWarmupNDigitMu1)+"f  %14."+to_string(fWarmupNDigitMu1)+"f";
             sprintf(buf2,format.c_str(),fWMu1Rnd[i].first,fWMu1Rnd[i].second);
          }
          // --- write
          strm<<buf<<buf2<<endl;
      }
      }
      strm<<"}}"<<endl;

   strm<<endl<<endl;

   // BinGrid
   strm<<"Warmup.Binning {{"<<endl;
   // table header
   strm<<"    ObsBin";
   for (unsigned int idim = 0 ; idim<NDim ; idim++) {
      sprintf(buf,"  %9s_Lo  %9s_Up",DimLabel[idim].c_str() ,DimLabel[idim].c_str());
      strm<<buf;
   }
   sprintf(buf,"  %12s  %12s","BinSize","EventCount");
   strm<<buf<<endl;

   // table values
   // KR: increase precision, otherwise binsize factors with e.g. Pi or bin borders in multiples of Pi lead to warnings!
   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {

      int nEvBin = 0;
      for ( int ip=0 ; ip<GetNSubprocesses() ; ip++ )
         nEvBin += GetTheCoeffTable()->fWgt.WgtObsNumEv[ip][i];

      sprintf(buf,"    %4d    ",i); // obsbin
      strm<<buf;
      for (unsigned int idim = 0 ; idim<NDim ; idim++) {
         sprintf(buf,"  % -#12.6g  % -#12.6g",Bin[i][idim].first , Bin[i][idim].second);
         strm<<buf;
      }
      sprintf(buf,"  % -#12.6g  %12d",BinSize[i],nEvBin);
      strm<<buf<<endl;
   }
   strm<<"}}"<<endl;

}


// ___________________________________________________________________________________________________
string fastNLOCreate::GetWarmupHeader(int iScale, string minmax) {
   string Descript = GetTheCoeffTable()->ScaleDescript[0][iScale];
   // replace all 'spaces' with _
   replace(Descript.begin(), Descript.end(), ' ', '_');
   // add 'minmax'
   string ret = "";
   ret += Descript;
   ret += "_";
   ret += minmax;
   return ret;
}



// ___________________________________________________________________________________________________
string fastNLOCreate::GetWarmupTableFilename() {
   if (!fWarmupFilename.empty()) {
      logger.debug["GetWarmupTableFilename"]<< "Preset name is: " << fWarmupFilename << endl;
      return fWarmupFilename;
   } else {
      // If we do not have a warmup filename yet, we generate a generic one
      string ret = fSteerfile;
      size_t pos = ret.find(".str");
      if (pos != std::string::npos) ret.erase(pos,4);
      pos = ret.find(".steer");
      if (pos != std::string::npos) ret.erase(pos,6);
      ret += "_";
      ret += GetScenName();
      logger.debug["GetWarmupTableFilename"]<< "Scenario name = " << GetScenName() << " bins" << endl;
      // TODO: KR: Consider keeping all warmup filenames to .wrm
      //      ret += "_warmup.txt";
      ret += ".wrm";
      SetWarmupTableFilename(ret);
      logger.debug["GetWarmupTableFilename"]<<"The warmup filename is: " << ret << endl;
      return ret;
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::SetWarmupTableFilename(string filename) {
   fWarmupFilename = filename;
}


// ___________________________________________________________________________________________________
void fastNLOCreate::AdjustWarmupValues() {
   //! Adjust warmup-values found to supposedly
   //! more reasonable values.
   //!
   //! Do this ONLY ONCE on COPY of actual values
   //! just before writing out to the warm-up table.
   //!
   //! 1. Round warm-up values up/down, if they are
   //!    4% close to the bin boundary
   //!    -> if more than 70% of all bins are
   //!       close to the bin boundary, then round all
   //! 2. Round values up/down, by mostly 3%
   //!    to next reasonable value
   //! 3. Round lower x-values down by 20%


   // ---------------------------------------
   // 3. 'round' lower x-values down
   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      fWxRnd[i].first  = fWx[i].first;
      fWxRnd[i].second = fWx[i].second;
   }
   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      //fWxRnd[i].first *= (1.-xdown);
      if ( fWxRnd[i].first >= 0.8 ) fWxRnd[i].first = 1.e-4;  // not filled!
      else if ( fWxRnd[i].first >= 0.09 ) fWxRnd[i].first = 0.09;  // minimum high-x
      double lx = log10(fWxRnd[i].first);
      int ex   = lx-1;
      double mant = fWxRnd[i].first/pow(10,ex);
      int imant = mant*10; // floor()
      // more safety margins
      int nEvBin = 0;
      for ( int ip=0 ; ip<GetNSubprocesses() ; ip++ )
         nEvBin += GetTheCoeffTable()->fWgt.WgtObsNumEv[ip][i];
      imant -= fWarmupXMargin;// safety margin
      if ( nEvBin < 100 ) imant-=4; // more safety
      else if ( nEvBin < 1000 ) imant-=2; // more safety
      else if ( nEvBin > 1000000 ) imant+=2; // less safety
      if ( imant%2==1 ) imant-=1; // only 0.0, 0.2, 0.4...
      fWxRnd[i].first = imant*pow(10,ex-1);//imant-1
      printf("          \t%8.3e   %8.3e  %8.1e   n=%d\n",fWx[i].first,fWxRnd[i].first,fWxRnd[i].first,nEvBin);
   }

}



// ___________________________________________________________________________________________________
void fastNLOCreate::RoundValues(vector<pair<double,double> >& wrmmu, int nthdigit) {
   //! Round warmup values up (down) if third relevant
   //! digit is a 9 (0)
   //! lower values are only rounded down,
   //! upper values are only rounded up

   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      if ( wrmmu[i].second!=0 && fabs(wrmmu[i].first/wrmmu[i].second-1) > 1.e-4) {
         // round only, if the values are different!
         // otherwise it is a 'fixed' scale and we have to
         // store exactly that values
         wrmmu[i].first  -= pow(10,-1*nthdigit-1)*5;
         wrmmu[i].second += pow(10,-1*nthdigit-1)*5;
      }
   }
   /*
   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      int nthlo = nthdigit;
      int nthup = nthdigit;
      if ( wrmmu[i].first<1 ) nthlo+=1;
      if ( wrmmu[i].second<1 ) nthup+=1;
      int lon = GetNthRelevantDigit(wrmmu[i].first,nthlo);
      int upn = GetNthRelevantDigit(wrmmu[i].second*(1+1.e-7),nthup);
      // lo value
      if ( lon==0 && wrmmu[i].first>1.e-6 ) {
         int ord = log10(wrmmu[i].first);
         int wrmrnd = wrmmu[i].first*pow(10.,-ord+nthlo-1);
         wrmmu[i].first = (double)wrmrnd / pow(10.,-ord+nthlo-1);
      }
      // up value
      if ( upn==9 && wrmmu[i].second>1.e-6 ) {
         int ord = log10(wrmmu[i].second);
         int wrmrnd = wrmmu[i].second*pow(10.,-ord+nthup-1) + 1;
         wrmmu[i].second = (double)wrmrnd / pow(10.,-ord+nthup-1);
      }
   }
      */
}


// ___________________________________________________________________________________________________
int fastNLOCreate::GetNthRelevantDigit(double val, int n) {
   int ord = log10(val);
   double res = fmod(val,pow(10.,ord-n+2));
   double resres=res-fmod(res,pow(10.,ord-n+1));
   double valn=resres/pow(10.,ord-n+1);
   return (int)(valn+0.999);
}


// ___________________________________________________________________________________________________
int fastNLOCreate::CheckWarmupValuesIdenticalWithBinGrid(vector<pair<double,double> >& wrmmu) {
   //! Check, where scale variable is identical
   //! with measured variable and hence
   //! warmu-values should be identical with
   //! bin grid.
   //!
   //! returns idim, if identity was found
   //! returns -1 else.
   //!
   //! If more than 70% of all bins are
   //! closer than 4% to the bin boundary, then
   //! identity is assumed
   const double bclose = 0.04;
   const double minallbins = 0.7;

   vector<int > nbinlo(NDim);
   vector<int > nbinup(NDim);
   for (int idim = int(NDim)-1 ; idim>=0 ; idim--) {
      for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
         if ( Bin[i][idim].first != 0 ) {
            if (wrmmu.size() <= i) {
                logger.error["CheckWarmupValuesIdenticalWithBinGrid"]
                << "Warmup values contains only " << wrmmu.size() << " bins" << endl;
            }
            double diff = wrmmu[i].first/Bin[i][idim].first - 1.;
            // lo-bin
            if ( diff < bclose  && diff >= 0.)
               nbinlo[idim]++;
         }
         else {
            if ( wrmmu[i].first < 1.e-4 )
               nbinlo[idim]++;
         }
         if ( Bin[i][idim].second != 0 ) {
            // up-bin
            double diff = 1. - wrmmu[i].second/Bin[i][idim].second;
            if ( diff < bclose && diff >= 0 )
               nbinup[idim]++;
         }
         else {
            if ( wrmmu[i].second < 1.e-4 )
               nbinup[idim]++;
         }
      }
   }
   // // sanity check (round only in one dimension
   //    for ( int idim = 0 ; idim<NDim ; idim++ ) {
   //       for ( int jdim = idim+1 ; jdim<NDim ; jdim++ ) {
   //    int nbri = nbinlo[idim]+nbinup[idim];
   //    int nbrj = nbinlo[jdim]+nbinup[jdim];
   //    if ( nbri>0 && nbrj>0 ){
   //       cout<<endl;
   //       logger.warn["CheckWarmupValuesIdenticalWithBinGrid"]
   //          <<"Adjusted warmup values to bin boundaries in different observables.\n"
   //          <<"\t\tThis may yield unreasonable results. Please check warmup-table carefully!\n"<<endl;
   //    }
   //       }
   //    }
   // round all bins if applicable
   for (unsigned int idim = 0 ; idim<NDim ; idim++) {
      logger.debug["CheckWarmupValuesIdenticalWithBinGrid"]<<"found nbinlo="<<nbinlo[idim]<<" and nbinup="<<nbinup[idim]<<endl;
      double frac= (nbinlo[idim]+nbinup[idim])/(2.*(int)GetNObsBin());
      if (frac>minallbins) {   // round all bins
         logger.info["CheckWarmupValuesIdenticalWithBinGrid"]
            <<"Found that "<<frac*100<<"% of the warmup values are close (<"<<bclose*100.<<"%) to a bin boundary in '"
            << DimLabel[idim]<<"' (Dim "<<idim<<").\n"
            <<"Using these bin boundaries as warm-up values."<<endl;
         for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
            wrmmu[i].first = Bin[i][idim].first;
            wrmmu[i].second = Bin[i][idim].second;
         }
         return idim;
      }
   }
   return -1;
}



// ___________________________________________________________________________________________________
void  fastNLOCreate::InitGrids() {
   logger.debug["InitGrids"]<<endl;
   if (fKernX1.empty()) logger.error["InitGrids"]<<"Interpolation kernels must be initialized before calling this function."<<endl;

   if (fIsFlexibleScale) {
      fastNLOCoeffAddFlex* c = (fastNLOCoeffAddFlex*)GetTheCoeffTable();
      //       if ( (c->GetNPDF()==2 && c->GetNPDFDim() == 1) || (c->GetNPDF()==1)   ) {;} // ok!
      //       else {
      //         logger.error["InitGrids"]<<"Only half-matrix or DIS implemented."<<endl; exit(1);
      //       }
      c->ScaleNode1.resize(GetNObsBin());
      c->ScaleNode2.resize(GetNObsBin());
      c->XNode1.resize(GetNObsBin());
      if (c->GetNPDFDim() == 2)
         c->XNode2.resize(GetNObsBin());

      vector<vector<vector<vector<vector<double> > > > > stype(GetNObsBin());
      for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
         c->ScaleNode1[i] = fKernMu1[i]->GetGrid();
         c->ScaleNode2[i] = fKernMu2[i]->GetGrid();
         c->XNode1[i]     = fKernX1[i]->GetGrid();
         if (c->GetNPDFDim() == 2)
            c->XNode2[i]     = fKernX2[i]->GetGrid();

         // SigmaTilde [NObsBins] ['n' x-nodes] [n s1-Nodes] [n s2-Nodes] [nsubproc]
         int nxmax = GetNxmax(fKernX1[i]->GetGridPtr() , fKernX2[i]->GetGridPtr());
         stype[i].resize(nxmax);
         for (unsigned int x = 0 ; x<stype[i].size() ; x++) {
            stype[i][x].resize(c->ScaleNode1[i].size());
            for (unsigned int m1 = 0 ; m1<stype[i][x].size() ; m1++) {
               stype[i][x][m1].resize(c->ScaleNode2[i].size());
               for (unsigned int mu2 = 0 ; mu2<stype[i][x][m1].size() ; mu2++) {
                  stype[i][x][m1][mu2].resize(GetNSubprocesses());
               }
            }
         }
      }
      c->SigmaTildeMuIndep      = stype;
      c->SigmaTildeMuRDep       = stype;
      c->SigmaTildeMuFDep       = stype;

      c->SigmaTildeMuRRDep      = stype;
      c->SigmaTildeMuFFDep      = stype;
      c->SigmaTildeMuRFDep      = stype;
      //c->SigmaTildeMuIndep(GetNObsBin());///
   }

   else {
      fastNLOCoeffAddFix* c = (fastNLOCoeffAddFix*)GetTheCoeffTable();
      //       if ( (c->GetNPDF()==2 && c->GetNPDFDim() == 1)  ) {;} // ok!
      //       else {
      //         //logger.error["InitGrids"]<<"Only half-matrix is implemented for grids for fixed-scale tables."<<endl; exit(1);
      //       }

      int nscalevar = fScaleFac.size();
      if (nscalevar==0) {
         logger.error["InitGrids"]<<"No scale factors found."<<endl;
      }
      c->Nscalevar.resize(1);                   // 1 = NScaleDim
      c->Nscalevar[0]  = nscalevar;             // testing

      c->ScaleFac.resize(1);                    // 1 = NScaleDim
      c->ScaleFac[0] = fScaleFac;
      c->XNode1.resize(GetNObsBin());
      for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
         c->XNode1[i]     = fKernX1[i]->GetGrid();
      }
      if (c->GetNPDFDim() == 2) {   // both hadrons have same x-grid in this implementation
         c->XNode2 = c->XNode1;
      }

      int nscalenode = fKernMuS[0][0]->GetGrid().size();
      // scale nodes
      fastNLOTools::ResizeVector(c->ScaleNode, GetNObsBin(), 1, nscalevar, nscalenode);
      for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
         for (int k=0; k<nscalevar; k++) {
            c->ScaleNode[i][0][k] = fKernMuS[i][k]->GetGrid();
         }
      }

      c->ResizeSigmaTilde();
   }

}


// ___________________________________________________________________________________________________
void  fastNLOCreate::InitInterpolationKernels() {
   //!
   //! initialize members for interpolation
   //!

   logger.debug["InitInterpolationKernels"]<<endl;
   if (fIsWarmup) {
      logger.error["InitInterpolationKernels"]<<"Interpolation kernels can only be initialized in production runs. Warmup values must be known."<<endl;
   }

   fKernX1.resize(GetNObsBin());
   fKernX2.resize(GetNObsBin());
   if (fIsFlexibleScale) {
      fKernMu1.resize(GetNObsBin());
      fKernMu2.resize(GetNObsBin());
   } else {
      fKernMuS.resize(GetNObsBin());
      for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
         fKernMuS[i].resize(fScaleFac.size());
      }
   }
   /*
   vector<double> wrmX = DOUBLE_COL_NS(Warmup.Values,x_min,fSteerfile);
   vector<double> wrmMu1Up, wrmMu1Dn;
   wrmMu1Dn = read_steer::getdoublecolumn("Warmup.Values",GetWarmupHeader(0,"min"),fSteerfile);
   wrmMu1Up = read_steer::getdoublecolumn("Warmup.Values",GetWarmupHeader(0,"max"),fSteerfile);
   if (wrmMu1Dn.size()!=GetNObsBin() || wrmMu1Up.size()!= GetNObsBin()) {
      logger.error["InitInterpolationKernels"]<<"Could not read warmup values for Mu1. Exiting."<<endl;
      exit(1);
   }
   vector<double> wrmMu2Up, wrmMu2Dn;
   if (fIsFlexibleScale) {
      wrmMu2Dn = read_steer::getdoublecolumn("Warmup.Values",GetWarmupHeader(1,"min"),fSteerfile);
      wrmMu2Up = read_steer::getdoublecolumn("Warmup.Values",GetWarmupHeader(1,"max"),fSteerfile);
      if (wrmMu2Dn.size()!=GetNObsBin() || wrmMu2Up.size()!= GetNObsBin()) {
         logger.error["InitInterpolationKernels"]<<"Could not read warmup values for Mu2. Exiting."<<endl;
         exit(1);
      }
   }
   */

   vector<double> wrmX;
   vector<double> wrmMu1Up, wrmMu1Dn;
   vector<double> wrmMu2Up, wrmMu2Dn;
   if ( GetTheCoeffTable()->GetIRef() == 0 ) {
      wrmX = GetColumnFromTable(fWarmupConsts.Values, 1) ;// DOUBLE_COL_NS(Warmup.Values,x_min,fSteerfile);
      wrmMu1Dn = GetColumnFromTable(fWarmupConsts.Values, 3) ;//read_steer::getdoublecolumn("Warmup.Values",GetWarmupHeader(0,"min"),fSteerfile);
      wrmMu1Up = GetColumnFromTable(fWarmupConsts.Values, 4) ;//read_steer::getdoublecolumn("Warmup.Values",GetWarmupHeader(0,"max"),fSteerfile);
      if (wrmMu1Dn.size()!=GetNObsBin() || wrmMu1Up.size()!= GetNObsBin()) {
         logger.error["InitInterpolationKernels"]<<"Could not read warmup values for Mu1. Exiting."<<endl;
         exit(1);
      }
      if (fIsFlexibleScale) {
         wrmMu2Dn = GetColumnFromTable(fWarmupConsts.Values, 5) ;//read_steer::getdoublecolumn("Warmup.Values",GetWarmupHeader(1,"min"),fSteerfile);
         wrmMu2Up = GetColumnFromTable(fWarmupConsts.Values, 6) ;//read_steer::getdoublecolumn("Warmup.Values",GetWarmupHeader(1,"max"),fSteerfile);
         if (wrmMu2Dn.size()!=GetNObsBin() || wrmMu2Up.size()!= GetNObsBin()) {
            logger.error["InitInterpolationKernels"]<<"Could not read warmup values for Mu2. Exiting."<<endl;
            exit(1);
         }
      }
   }
   else {
      wrmX.resize(GetNObsBin(),0);
      wrmMu1Dn.resize(GetNObsBin(),0);
      wrmMu1Up.resize(GetNObsBin(),1);
      if (fIsFlexibleScale) {
         wrmMu2Dn.resize(GetNObsBin(),0);
         wrmMu2Up.resize(GetNObsBin(),1);
      }
   }

   int npdf = GetTheCoeffTable()->GetNPDF();

   if ( fScenConsts.X_Kernel.empty() ){
      logger.error["InitInterpolationKernels"]<<"X_Kernel not found in fScenarioConstants."<<endl; exit(1);
   }
   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      // ------------------------------------------------
      // init x-interpolation kernels
      // ------------------------------------------------
      logger.debug["InitInterpolationKernels"]<<"Make x grid for obsbin="<<i<<endl;

      // Create x grids with X_NNodes+1 nodes up to x_max = 1.
      // The additional last node will be removed again below.
      int nxtot = fScenConsts.X_NNodes + 1;
      if ( fScenConsts.X_NNodeCounting == "NodesPerMagnitude" ) { // "NodesMax","NodesPerBin","NodesPerMagnitude"
         fKernX1[i] = MakeInterpolationKernels(fScenConsts.X_Kernel,wrmX[i],1); // use 1 as upper x-value
         fKernX1[i]->MakeGridsWithNNodesPerMagnitude(fastNLOInterpolBase::TranslateGridType(fScenConsts.X_DistanceMeasure),nxtot);
         if (npdf == 2) {
            fKernX2[i] = MakeInterpolationKernels(fScenConsts.X_Kernel,wrmX[i],1);
            fKernX2[i]->MakeGridsWithNNodesPerMagnitude(fastNLOInterpolBase::TranslateGridType(fScenConsts.X_DistanceMeasure),nxtot);
         }
      }
      else if (fScenConsts.X_NNodeCounting == "NodesPerBin" ) { //
         fKernX1[i] = MakeInterpolationKernels(fScenConsts.X_Kernel,wrmX[i],1); // use 1 as upper x-value
         fKernX1[i]->MakeGrids(fastNLOInterpolBase::TranslateGridType(fScenConsts.X_DistanceMeasure),nxtot);
         if (npdf == 2) {
            fKernX2[i] = MakeInterpolationKernels(fScenConsts.X_Kernel,wrmX[i],1);
            fKernX2[i]->MakeGrids(fastNLOInterpolBase::TranslateGridType(fScenConsts.X_DistanceMeasure),nxtot);
         }
      }
      else if (fScenConsts.X_NNodeCounting == "NodesMax" ) { //
         // generate grid for maximum x-range
         double xmin = 1;
         for ( unsigned int xi=0 ; xi<wrmX.size() ; xi++ ) xmin = min(xmin,wrmX[xi]);
         fastNLOInterpolBase* kernmin = MakeInterpolationKernels(fScenConsts.X_Kernel,xmin,1); // use 1 as upper x-value
         kernmin->MakeGrids(fastNLOInterpolBase::TranslateGridType(fScenConsts.X_DistanceMeasure),nxtot);
         // find xmin
         const vector<double>& xg = kernmin->GetGrid();
         int nxbin = nxtot;
         for ( unsigned xi = 0 ; xi<xg.size() ; xi++ ) {
            //cout<<"iBin="<<i<<"\txi="<<xi<<"\txg[xi]="<<xg[xi]<<"\twrmX[i]="<<wrmX[i]<<"\tnxbin="<<nxbin<<"\tnxtot="<<nxtot<<endl;
            if ( xg[xi] > wrmX[i] ) break;
            nxbin--;
            xmin = xg[xi];
         }
         logger.info["InitInterpolationKernels"]<<"Using x-grid in bin "<<i<<": x-min="<<xmin<<"\tx-warmup="<<wrmX[i]<<"\tnxbin="<<nxbin<<endl;
         fKernX1[i] = MakeInterpolationKernels(fScenConsts.X_Kernel,xmin,1); // use 1 as upper x-value
         fKernX1[i]->MakeGrids(fastNLOInterpolBase::TranslateGridType(fScenConsts.X_DistanceMeasure),nxbin);
         if (npdf == 2) {
            fKernX2[i] = MakeInterpolationKernels(fScenConsts.X_Kernel,xmin,1);
            fKernX2[i]->MakeGrids(fastNLOInterpolBase::TranslateGridType(fScenConsts.X_DistanceMeasure),nxbin);
         }
         delete kernmin;
      }
      else {
         // error
         logger.error["InitInterpolationKernels"]<<"Cannot understand node counting: "<<fScenConsts.X_NNodeCounting<<"."<<endl;
         logger.error["InitInterpolationKernels"]<<"Supported options are: 'NodesPerMagnitude', 'NodesPerBin' and 'NodesMax'."<<endl;
      }

      // Remove last node at x = 1; is multiplied by PDFs equalling zero anyway.
      fKernX1[i]->RemoveLastNode();
      if (npdf == 2)
         fKernX2[i]->RemoveLastNode();

      // ------------------------------------------------
      // init scale1-interpolation kernels
      // ------------------------------------------------
      int nqtot1 = fScenConsts.Mu1_NNodes;
      if (fIsFlexibleScale) {
         logger.debug["InitInterpolationKernels"]<<"Make Mu1 grid for obsbin="<<i<<endl;
         fKernMu1[i] = MakeInterpolationKernels(fScenConsts.Mu1_Kernel,wrmMu1Dn[i],wrmMu1Up[i]);
         fKernMu1[i]->MakeGrids(fastNLOInterpolBase::TranslateGridType(fScenConsts.Mu1_DistanceMeasure),nqtot1);
         // ------------------------------------------------
         // init scale2-interpolation kernels
         // ------------------------------------------------
         int nqtot2 = fScenConsts.Mu2_NNodes;
         logger.debug["InitInterpolationKernels"]<<"Make Mu2 grid for obsbin="<<i<<endl;
         fKernMu2[i] = MakeInterpolationKernels(fScenConsts.Mu2_Kernel,wrmMu2Dn[i],wrmMu2Up[i]);
         fKernMu2[i]->MakeGrids(fastNLOInterpolBase::TranslateGridType(fScenConsts.Mu2_DistanceMeasure),nqtot2);
      }
      else {
         for (unsigned int k = 0 ; k<fScaleFac.size() ; k++) {
            double muDn = fScaleFac[k]*wrmMu1Dn[i];
            double muUp = fScaleFac[k]*wrmMu1Up[i];
            fKernMuS[i][k] = MakeInterpolationKernels(fScenConsts.Mu1_Kernel,muDn,muUp);
            fKernMuS[i][k]->MakeGrids(fastNLOInterpolBase::TranslateGridType(fScenConsts.Mu1_DistanceMeasure),nqtot1);
         }
      }
   }
}



// ___________________________________________________________________________________________________
std::vector<double> fastNLOCreate::GetColumnFromTable(const std::vector<std::vector<double> >& table, int iCol){
   //! Get a column from a table
   vector<double> ret;
   for ( unsigned int i = 0 ; i<table.size(); i++ ) {
      if ( (int)table[i].size() <= iCol ) {
         logger.error["GetColumnFromTable"]<< "Table does not have enough columns in row "<<i<<". Exiting."<<endl;
         exit(1);
      }
      ret.push_back(table[i][iCol]);
   }
   return ret;
}



// ___________________________________________________________________________________________________
fastNLOInterpolBase* fastNLOCreate::MakeInterpolationKernels(string KernelName, double xdn, double xup) {
   //! This function identifies the string-identifier
   //! and creates the corresponding fastNLO Interpolation kernel

   if (KernelName == "CatmullRom" || KernelName == "Catmull")
      return new fastNLOInterpolCatmullRom(xdn,xup);
   else if (KernelName == "Lagrange")
      return new fastNLOInterpolLagrange(xdn,xup);
   else if (KernelName == "Linear")
      return new fastNLOInterpolLinear(xdn,xup);
   else if (KernelName == "OneNode")
      return new fastNLOInterpolOneNode(xdn,xup);
   // else if ( KernelName == "...") // todo implement other kernels here!
   //   return ...
   else {
      logger.warn["MakeInterpolationKernels"]<<"Cannot find kernel routine:" <<KernelName<<" or kernel not (yet) implemented. Exiting."<<endl;
      exit(1);
   }
   return NULL; // default return
}


// ___________________________________________________________________________________________________
fastNLOReader* fastNLOCreate::SetIsReferenceTable(fastNLOReader* fnloread) {
   //! set this table/contribution to become a reference contribution
   //! If fnloread is set to NULL, the weights are assumed to be already
   //! multiplied by PDF and alpha_s values.
   //! If fnloread is provided, then it is assumed that the weights have
   //! the same units and format as for the filling of  default tables and
   //! those need to be multiplied by PDF and alpha_s values.
   //!
   //! Function returns the input pointer without changes.
   //!
   GetTheCoeffTable()->SetIRef();
   fReader = fnloread;

   // --- set 'fReader' if needed
   // if ( fReader ) { // todo, do we need to call something here?
   //    fReader->InitPDF();//
   // }

   // --- adjust the Coeff-table
   fScenConsts.X_NNodes = 1 -1; // because we are using +1 later
   //fScenConsts.X_NoOfNodesPerMagnitude = false;
   fScenConsts.X_NNodeCounting = "NodesPerBin";
   fScenConsts.Mu1_NNodes = 1;
   fScenConsts.Mu2_NNodes = 1;
   fScenConsts.X_Kernel   = "OneNode";
   fScenConsts.Mu1_Kernel = "OneNode";
   fScenConsts.Mu2_Kernel = "OneNode";
   InitInterpolationKernels(); // resize all members

   // return the reader
   return fReader;
}



// ___________________________________________________________________________________________________
bool fastNLOCreate::TestParameterInSteering(const string& key) const {
   //! Get flag if parameter exists in steering card
   return read_steer::getexist(key, fSteerfile);
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& key, bool& val) const {
   //! Get boolean value from steering with key 'key'.
   //! Alternatively, also BOOL_NS(`key`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering keys, one should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static bool IsCMS;
   //! static bool gotval = GetParameterFromSteering("MjjCut",IsCMS);
   //! if (!gotval) cout<<"Error! Could not find boolean parameter MjjCut in steering file."<<endl;
   //!
   //! Function returns 'false' if key was not found in steering file

   bool exist = read_steer::getexist(key,fSteerfile);
   if ( exist )
      val = read_steer::getbool(key,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& key, int& val) const {
   //! Get integer value from steering with key 'key'.
   //! Alternatively, also INT_NS(`key`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering keys, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static int nJetMin;
   //! static bool gotval = GetParameterFromSteering("nJetMin",nJetMin);
   //! if (!gotval) cout<<"Error! Could not find integer parameter nJetMin in steering file."<<endl;
   //!
   //! Function returns 'false' if key was not found in steering file

   bool exist = read_steer::getexist(key,fSteerfile);
   if ( exist )
      val = read_steer::getint(key,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& key, double& val) const {
   //! Get boolean value from steering with key 'key'.
   //! Alternatively, also DOUBLE_NS(`key`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering keys, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static double MjjCut;
   //! static bool gotval = GetParameterFromSteering("MjjCut",MjjCut);
   //! if (!gotval) cout<<"Error! Could not find parameter MjjCut in steering file."<<endl;
   //!
   //! Function returns 'false' if key was not found in steering file

   bool exist = read_steer::getexist(key,fSteerfile);
   if ( exist )
      val = read_steer::getdouble(key,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& key, string& val) const {
   //! Get string value from steering with key 'key'.
   //! Alternatively, also STRING_NS(`key`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering keys, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static string text;
   //! static bool gotval = GetParameterFromSteering("MyText",text);
   //! if (!gotval) cout<<"Error! Could not find parameter MyText in steering file."<<endl;
   //!
   //! Function returns 'false' if key was not found in steering file

   bool exist = read_steer::getexist(key,fSteerfile);
   if ( exist )
      val = read_steer::getstring(key,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& key, vector<int>& val) const {
   //! Get integer vector from steering with key 'key'.
   //! Alternatively, also INT_ARR(`key`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering keys, you  should not
   //! call this function too frequently.
   //!
   //! Function returns 'false' if key was not found in steering file

   bool exist = read_steer::getexist(key,fSteerfile);
   if ( exist )
      val = read_steer::getintarray(key,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& key, vector<double>& val) const {
   //! Get vector of doubles from steering with key 'key'.
   //! Alternatively, also DOUBLE_ARR_NS(`key`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering keys, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static vector<double> FlexiCuts;
   //! static bool gotval = GetParameterFromSteering("FlexiCuts",FlexiCuts);
   //! if (!gotval) cout<<"Error! Could not find vector FlexiCuts in steering file."<<endl;
   //!
   //! Function returns 'false' if key was not found in steering file

   bool exist = read_steer::getexist(key,fSteerfile);
   if ( exist )
      val = read_steer::getdoublearray(key,fSteerfile);
   return exist;

}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& key, vector<string>& val) const {
   //! Get vector of strings from steering with key 'key'.
   //! Alternatively, also STRING_ARR_NS(`key`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering keys, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static string text;
   //! static bool gotval = GetParameterFromSteering("MyText",text);
   //! if (!gotval) cout<<"Error! Could not find parameter MyText in steering file."<<endl;
   //!
   //! Function returns 'false' if key was not found in steering file

   bool exist = read_steer::getexist(key,fSteerfile);
   if ( exist )
      val = read_steer::getstringarray(key,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& key, vector<vector<int > >& val) const {
   //! Get vector of vectors of ints from steering with key 'key'.
   //! Alternatively, also INT_TAB_NS(`key`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering keys, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static string text;
   //! static bool gotval = GetParameterFromSteering("MyText",text);
   //! if (!gotval) cout<<"Error! Could not find parameter MyText in steering file."<<endl;
   //!
   //! Function returns 'false' if key was not found in steering file

   bool exist = read_steer::getexist(key,fSteerfile);
   if ( exist )
      val = read_steer::getinttable(key,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& key, vector<vector<double > >& val) const {
   //! Get vector of vector of doubles from steering with key 'key'.
   //! Alternatively, also DOUBLE_TAB_NS(`key`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering keys, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static string text;
   //! static bool gotval = GetParameterFromSteering("MyText",text);
   //! if (!gotval) cout<<"Error! Could not find parameter MyText in steering file."<<endl;
   //!
   //! Function returns 'false' if key was not found in steering file

   bool exist = read_steer::getexist(key,fSteerfile);
   if ( exist )
      val = read_steer::getdoubletable(key,fSteerfile);
   return exist;
}
