// Author: Daniel Britzger
// DESY, 29/07/2013
// ___________________________________________________________________________________________________
//! The fastNLOCreate class
/*!
   This class can generate/fill one single contribution for a fastNLO table.
   It supports fixed scale and flexible-scale tables.
   fastNLOCreate inherits from fastNLOTable, but is only able to hold one
   coefficient table as member.
   fastNLOCreate no enables to fill this coefficient talbe, i.e.  to add further
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
#include "fastnlotk/fastNLOCreate.h"
#include "fastnlotk/fastNLOTools.h"
#include "fastnlotk/read_steer.h"

#include "fastnlotk/fastNLOCoeffAddFlex.h"
#include "fastnlotk/fastNLOCoeffAddFix.h"
#include "fastnlotk/fastNLOInterpolCatmullRom.h"
#include "fastnlotk/fastNLOInterpolLagrange.h"
#include "fastnlotk/fastNLOInterpolLinear.h"
#include "fastnlotk/fastNLOInterpolOneNode.h"

using namespace std;


// ___________________________________________________________________________________________________
fastNLOCreate::fastNLOCreate() {
   SetClassName("fastNLOCreate");
}


// ___________________________________________________________________________________________________
fastNLOCreate::fastNLOCreate(string steerfile, fastNLO::GeneratorConstants GenConsts, fastNLO::ProcessConstants ProcConsts) {
   SetClassName("fastNLOCreate");
   ResetHeader();
   ReadSteering(steerfile);

   fGenConsts  = GenConsts;
   fProcConsts = ProcConsts;

   // KR: Add this line also to this constructor so that settings in steering, if they exist,
   //     take precedence over previously set default settings.
   ReadGenAndProcConstsFromSteering();

   Instantiate();
}


// ___________________________________________________________________________________________________
fastNLOCreate::fastNLOCreate(string steerfile) {
   SetClassName("fastNLOCreate");
   ResetHeader();
   ReadSteering(steerfile);

   ReadGenAndProcConstsFromSteering();

   Instantiate();
}


// ___________________________________________________________________________________________________
void fastNLOCreate::ReadGenAndProcConstsFromSteering() {
   //! ReadGenAndProcConstsFromSteering()
   //! If generator and process constants have not been set in
   //! the constructor, then obtain these values from the steering
   //! file.

   debug["ReadGenAndProcConstsFromSteering"]<<endl;

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

   // Process constants
   if (EXIST_NS(LeadingOrder,fSteerfile))        fProcConsts.LeadingOrder = INT_NS(LeadingOrder,fSteerfile);
   if (EXIST_NS(UnitsOfCoefficients,fSteerfile)) fProcConsts.UnitsOfCoefficients = INT_NS(UnitsOfCoefficients,fSteerfile);
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
      fProcConsts.PDFCoeffLO   = ReadPartonCombinations(0);
      fProcConsts.PDFCoeffNLO  = ReadPartonCombinations(1);
      if ( fProcConsts.IPDFdef3NNLO > 0 )
         fProcConsts.PDFCoeffNNLO = ReadPartonCombinations(2);
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
      } else if ( fProcConsts.NPDF == 2 ) {
         if (EXIST_NS(AsymmetricProcesses,fSteerfile)) {
            fProcConsts.AsymmetricProcesses.clear();
            vector<vector<int> > asym = INT_TAB_NS(AsymmetricProcesses,fSteerfile);
            for (unsigned int i = 0 ; i<asym.size() ; i++) {
               if ( asym[i].size() != 2 ) {
                  error["ReadGenAndProcConstsFromSteering"]<<"Asymmetric process "<<asym[i][0]<<", must have exactly one counter process."<<endl;
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
   // Try to get warm-up values.
   // Otherwise a warm-up run will be initialized.
   GetWarmupValues();

   ILOord = fProcConsts.LeadingOrder;
   fIOrd = ILOord; // initialize with LO

   // init bin grid
   if (fIsWarmup) ReadBinning();                      // if warmup, then always read binning from steering.
   else if (BOOL_NS(ReadBinningFromSteering,fSteerfile)) {
      ReadBinning();
      CheckWarmupConsistency();
   } else {UseBinGridFromWarmup();}

   // now create one coefficient tasble.
   InitCoeffTable();
   // no info output for the following calls
   bool vol = info.GetSpeak();
   info.DoSpeak(false);
   SetOrderOfAlphasOfCalculation(fIOrd);
   info.DoSpeak(vol);//reset verbosity level

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
void fastNLOCreate::ReadSteering(string steerfile) {
   //! read in steering file
   //! The filename of the steering file
   //! is used as the 'namespace' of labels in read_steer
   debug["ReadSteering"]<<"Steerfile = "<<steerfile<<endl;
   fSteerfile =  steerfile;
   READ_NS(steerfile,steerfile);

   SetGlobalVerbosity(STRING_NS(GlobalVerbosity,fSteerfile));
   if (info.GetSpeak())
      PRINTALL();

   // header
   SetScenName(STRING_NS(ScenarioName,fSteerfile));
   SetItabversion(22000);

   // scenario specific things
   Ipublunits   = INT_NS(PublicationUnits,fSteerfile);
   ScDescript   = STRING_ARR_NS(ScenarioDescription,fSteerfile);
   Ecms         = DOUBLE_NS(CenterOfMassEnergy,fSteerfile);   // is often superseeded by generator-specific code.
   INormFlag    = 0;
   fIOrd        = 0;// has to be set by generator
   SetFilename(STRING_NS(OutputFilename,fSteerfile));

   fIsFlexibleScale = BOOL_NS(FlexibleScaleTable,fSteerfile);
   fApplyPDFReweight = BOOL_NS(ApplyPDFReweighting,fSteerfile);
   SetOutputPrecision(INT_NS(OutputPrecision,fSteerfile));

   //if ( !fIsFlexibleScale )  ReadScaleFactors() ; // is called only when setting order of calculation
}


// ___________________________________________________________________________________________________
vector<vector<pair<int,int> > > fastNLOCreate::ReadPartonCombinations(int ord) {
   //! Read PDF linear combinations from steering file
   //! and convert to internal format


   vector<vector<int> > PartonCombinations;
   if ( ord==0 )      {
      PartonCombinations = INT_TAB_NS(PartonCombinationsLO,fSteerfile);
      if ( (int)PartonCombinations.size() != fProcConsts.NSubProcessesLO ) {
         error["ReadPartonCombinations"]<<"Number of parton combinations for LO processes must be identical to number of subprocesses. NSubProcessesLO="
                                        <<fProcConsts.NSubProcessesLO<<", # parton combinations="<<PartonCombinations.size()<<". Exiting." <<endl;
         cout<<"PartonCombinations.size()="<<PartonCombinations.size()<<",  fProcConsts.NSubProcessesLO="<<fProcConsts.NSubProcessesLO<<endl;
         exit(1);
      }
   }
   else if ( ord==1 ) {
      PartonCombinations = INT_TAB_NS(PartonCombinationsNLO,fSteerfile);
      if ( (int)PartonCombinations.size() != fProcConsts.NSubProcessesNLO ) {
         error["ReadPartonCombinations"]<<"Number of parton combinations for NLO processes must be identical to number of subprocesses.  NSubProcessesNLO="
                                        <<fProcConsts.NSubProcessesNLO<<", # parton combinations="<<PartonCombinations.size()<<". Exiting." <<endl;
         exit(1);
      }
   }
   else if ( ord==2 ) {
      PartonCombinations = INT_TAB_NS(PartonCombinationsNNLO,fSteerfile);
      if ( (int)PartonCombinations.size() != fProcConsts.NSubProcessesNNLO ) {
         error["ReadPartonCombinations"]<<"Number of parton combinations for NNLO processes must be identical to number of subprocesses.  NSubProcessesNNLO="
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
         error["ReadPartonCombinations"]<<"Row "<<k<<" PartonCombinations"<<sord[ord]<<" does not contain any information. Exiting."<<endl;
         exit(1);
      }
      int iSubProc = PartonCombinations[k][0];
      if ( iSubProc >= (int)PDFCoeff.size() ) {
         error["ReadPartonCombinations"]<<"Subprocess "<<iSubProc<<" in row "<<k+1<<" of PartonCombinations"<<sord[ord]<<" is larger than the total number of subprocesses. Exiting."<<endl;
         exit(1);
      }
      if ( PartonCombinations[k].size()%2 != 1 || PartonCombinations[k].size()<=1) {
         error["ReadPartonCombinations"]<<"Row "<<k<<" of PartonCombinations"<<sord[ord]<<" does not fit format: 'iProc [pair0] [pair1] ... [pairN]. Exiting"<<endl;
         exit(1);
      }
      if ( !PDFCoeff[iSubProc].empty() ) {
         error["ReadPartonCombinations"]<<"Subprocess "<<iSubProc<<" appears twice in the PartonCombinations"<<sord[ord]<<". Exiting."<<endl;
         exit(1);
      }

      for ( unsigned int i=1 ; i<PartonCombinations[k].size() ; i+=2 ) {
         debug["ReadPartonCombinations"]<<"Adding to subprocess "<<iSubProc<<" parton pair (" << PartonCombinations[k][i]<<","<<PartonCombinations[k][i+1]<<")."<<endl;
         int iPart1 = PartonCombinations[k][i];
         int iPart2 = PartonCombinations[k][i+1];
         if ( abs(iPart1) > 6 || abs(iPart2) > 6 ) {
            error["ReadPartonCombinations"]<<"Parton flavor is larger than 6. There is nothing beyond the top-quark. Exiting."<<endl;
            exit(1);
         }
         if (  b1[iPart1+6]  ) warn["ReadPartonCombinations"]<<"Parton "<<iPart1<<" of hadron 1 is used multiple times in PartonCombinations"<<sord[ord]<<"."<<endl;
         if (  b2[iPart2+6]  ) warn["ReadPartonCombinations"]<<"Parton "<<iPart2<<" of hadron 2 is used multiple times in PartonCombinations"<<sord[ord]<<"."<<endl;

         b1[iPart1+6] = true;
         b2[iPart2+6] = true;
         PDFCoeff[iSubProc].push_back(std::make_pair(iPart1,iPart2));
      }
   }

   // check if all subprocesses are filled
   for ( unsigned int k=1 ; k<PDFCoeff.size() ; k++ ) {
      if ( PDFCoeff[k].empty() ) {
         error["ReadPartonCombinations"]<<"PartonCombinations"<<sord[ord]<<" does not conatain any information about PDF for subprocess "<<k<<". Exiting."<<endl;
         exit(1);
      }
   }
   for ( int p = 1 ; p<12 ; p++ ) {
      // check if all partons are used
      if ( !b1[p] ) {
         error["ReadPartonCombinations"]<<"Parton "<<p-6<<" of hadron 1 is not used in PartonCombinations"<<sord[ord]<<". Exiting."<<endl; exit(1);
      }
      if ( !b2[p] ) {
         error["ReadPartonCombinations"]<<"Parton "<<p-6<<" of hadron 2 is not used in PartonCombinations"<<sord[ord]<<". Exiting."<<endl; exit(1);
      }
   }
   // check if all partons are used
   int p = 0;
   if ( !b1[p] )       warn["ReadPartonCombinations"]<<"Parton "<<p-6<<" of hadron 1 is not used in PartonCombinations"<<sord[ord]<<". Exiting."<<endl;
   if ( !b2[p] )       warn["ReadPartonCombinations"]<<"Parton "<<p-6<<" of hadron 2 is not used in PartonCombinations"<<sord[ord]<<". Exiting."<<endl;
   p = 12;
   if ( !b1[p] )       warn["ReadPartonCombinations"]<<"Parton "<<p-6<<" of hadron 1 is not used in PartonCombinations"<<sord[ord]<<". Exiting."<<endl;;
   if ( !b2[p] )       warn["ReadPartonCombinations"]<<"Parton "<<p-6<<" of hadron 2 is not used in PartonCombinations"<<sord[ord]<<". Exiting."<<endl;;


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

   if (fIsFlexibleScale) {warn["ReadScaleFactors"]<<"This function is only reasonable for fixed-scale tables!"<<endl;}
   vector<double> svar = DOUBLE_ARR_NS(ScaleVariationFactors,fSteerfile);
   fScaleFac.resize(svar.size());
   if (svar.empty()) {
      // 'ScaleVariationFactors' not found -> using default
      warn["ReadScaleFactors"]<<"No list of scale-factors found in steering file. Using only scale-factor of '1'."<<endl;
      fScaleFac.push_back(1.0);
   } else if (GetTheCoeffTable() && GetTheCoeffTable()->IsLO()) {
      // scale-factors not needed -> using default of 1.0
      info["ReadScaleFactors"]<<"This is a leading-order run. There is no MuF scale dependence in LO. Using only scale factor of 1."<<endl;
      fScaleFac.resize(1);
      fScaleFac[0] = 1.0;
   } else if (fIsWarmup) {
      // scale-factors not needed -> using default of 1.0
      info["ReadScaleFactors"]<<"This is a warmup run. Using only scale factor of 1."<<endl;
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
            info["ReadScaleFactors"]<<"Found scale factor 1.0 two times in list ScaleVariationFactors. Ignoring second appearance."<<endl;
            fScaleFac.resize(fScaleFac.size()-1);
         } else if (fabs(svar[k]-1.0)<1.e-8) foundunity = true;
         else vtemp.push_back(svar[k]);
      }
      if (!foundunity) {
         error["ReadScaleFactors"]<<"Could not find scale factor of 1.0 in list ScaleVariationFactors. Exiting."<<endl;
         exit(1);
      }
      sort(vtemp.begin(),vtemp.end());

      fScaleFac[0] = 1.0;
      int s = 0;
      for (unsigned int k = 0 ; k<vtemp.size(); k++) {
         fScaleFac[s+1] = vtemp[k];
         if (fScaleFac[s+1] == fScaleFac[s]) {
            info["ReadScaleFactors"]<<"Found scale factor '"<<fScaleFac[k+1]<<"' two times in list ScaleVariationFactors. Ignoring second appearance."<<endl;
            fScaleFac.resize(fScaleFac.size()-1);
            s--;
         }
         s++;
      }
   }
   info["ReadScaleFactors"] << "Using the following scale factors:" << endl;
   for (unsigned int k = 0 ; k<fScaleFac.size(); k++)
      info["ReadScaleFactors"] << "ScaleVar " << k << ": " << fScaleFac[k] << endl;
}


// ___________________________________________________________________________________________________
void fastNLOCreate::InitCoeffTable() {
   debug["InitCoeffTable"]<<endl;
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
   debug["CreateCoeffTable"]<<endl;
   if (!fCoeff.empty()) {
      error["CreateCoeffAddFix"]<<"Vector of coefficients must be empty, since only one coefficient table is allowed."<<endl;
      exit(1);
   }
   if (fIsFlexibleScale)
      return fastNLOTable::CreateCoeffTable(fCoeff.size(), new fastNLOCoeffAddFlex(NObsBin,ILOord));
   else
      return fastNLOTable::CreateCoeffTable(fCoeff.size(), new fastNLOCoeffAddFix(NObsBin));
}


// ___________________________________________________________________________________________________
void fastNLOCreate::ReadBinning() {
   // optimize read-in of bin grids
   // ToDo. Check sanity of bin-grid

   NDim         = INT_NS(DifferentialDimension,fSteerfile);
   //Scenario.SetNDim(NDim);
   if (STRING_ARR_NS(DimensionLabels,fSteerfile).size() < NDim) {
      error["ReadBinning"]<<"Each dimension needs a bin label. Exiting."<<endl;
      exit(1);
   }
   DimLabel     = STRING_ARR_NS(DimensionLabels,fSteerfile);
   if (INT_ARR_NS(DimensionIsDifferential,fSteerfile).size() < NDim) {
      error["ReadBinning"]<<"Each dimension need to specify if differential or not. Exiting."<<endl;
      exit(1);
   }
   IDiffBin     = INT_ARR_NS(DimensionIsDifferential,fSteerfile);
   DimLabel.resize(NDim); //safety
   IDiffBin.resize(NDim);

   bool AllDiff = true;
   bool AllBinInt = true;
   for (unsigned int i = 0 ; i<IDiffBin.size() ; i++) {
      AllDiff = AllDiff && (IDiffBin[i] == 1);
      AllBinInt = AllBinInt && (IDiffBin[i] != 1);
   }
   if (!AllDiff && !AllBinInt) {
      error["ReadBinning"]<<"All dimensions must be consistently either bin-integrated, or truly differential dimensions. Exiting."<<endl;
      exit(1);
   }

   if (AllDiff && NDim == 3) { error["ReadBinning"]<<"Fully differential and triple-differential binning not yet implemented. exiting"<<endl; exit(1);}

   // read single-differential bin grid
   if (NDim == 1) {
      vector<double> bgrid = DOUBLE_ARR_NS(SingleDifferentialBinning,fSteerfile);
      //      SetBinning1D(bgrid, DimLabel[0], IDiffBin[0]);
      SetBinningND(bgrid, NDim, IDiffBin);
   }

   // read double-differential bin grid
   else if (NDim==2) {
      vector<vector<double> > in = DOUBLE_TAB_NS(DoubleDifferentialBinning,fSteerfile);
      // New binning code
      SetBinningND(in, NDim, IDiffBin);
   }

   // read in triple-differential binning
   else if (NDim==3) {
      warn["ReadBinning"]<<"The code for reading of "<<NDim<<"-dimensional binnings was not fully tested. Please verify the code and remove this statement."<<endl;
      vector<vector<double> > in = DOUBLE_TAB_NS(TripleDifferentialBinning,fSteerfile);
      // New binning code
      SetBinningND(in, NDim, IDiffBin);
   } else {
      error["ReadBinning"]<<"Reading of "<<NDim<<"-binnings from steering is not yet implemented. Exiting"<<endl;
      exit(1);
   }


   // ---------------------------------
   //       Bin width
   // ---------------------------------
   if (BOOL_NS(CalculateBinSize,fSteerfile)) {
      BinSize.resize(NObsBin);
      bool idi = false;
      for (unsigned int i = 0 ; i<NObsBin ; i++) {
         BinSize[i] = 1;
         for (unsigned int d=0 ; d<NDim ; d++) {
            if (IDiffBin[d]==0) {
               // nothing todo
            } else if (IDiffBin[d]==1) {
               // nothing todo
               //warn["ReadBinning"]<<"Don't know how to handle truly differential bins for bin widths."<<endl;
            } else if (IDiffBin[d]==2) {
               BinSize[i] *= Bin[i][d].second-Bin[i][d].first;//UpBin[i][d]-LoBin[i][d];
               idi = true;
            }
         }
         // divide by binsizefactor, but only if at least one dimension is differential
         if (idi) BinSize[i] *= DOUBLE_NS(BinSizeFactor,fSteerfile);
      }
      if (!idi) debug["ReadBinning"]<<"BinSizeFactor is not being used, since no observable is calculated differential."<<endl;
   } else {
      // read in bin width
      warn["ReadBinning"]<<"Reading of bindwidth only poorly  implemented! Improve it and remove this message."<<endl;
      if (DOUBLE_ARR_NS(BinSize,fSteerfile).size()!=NObsBin) warn["ReadBinning"]<<"Number of bins of 'BinSize' not consistent with bin grid."<<endl;
      BinSize=DOUBLE_ARR_NS(BinSize,fSteerfile);
      BinSize.resize(NObsBin);
      for (unsigned int i = 0 ; i<NObsBin ; i++) {
         if (BinSize[i]==0) BinSize[i] = 1.0;
      }
   }

   info["ReadBinning"]<<"Read in successfully "<<NDim<<"-dimensional bin grid with "<<NObsBin<<" bins."<<endl;
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
   info["SetBinning1D"] << "VSI: Set all normalization factors to one." << endl;
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
   info["SetBinning1D"] << "VSID: Set all normalization factors to norm." << endl;
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
   info["SetBinning1D"] << "VSIV: Set vectors of lower and upper bin edges." << endl;
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
   info["SetBinning1D"] << "VVSI: Set all normalization factors to one." << endl;
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
   info["SetBinning1D"] << "VVSID: Set all normalization factors to norm." << endl;
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
      error["SetBinning1D"] << "Illegal choice of differentiality, idiff = " << idiff << ", aborted!" << endl;
      error["SetBinning1D"] << "Dimension must either be non-differential (idiff=0), " << endl;
      error["SetBinning1D"] << "point-wise differential (idiff=1), or bin-wise differential (idiff=2)." << endl;
      exit(1);
   }

   // Check input: Gaps are fine, overlaps are not. Overlaps might lead to double-counting in normalization.
   if ( blow.size() != bupp.size() ) {
      error["SetBinning1D"] << "Unequal size of bin border vectors, blow.size() = " << blow.size() <<
         ", bupp.size() = " << bupp.size() << ", aborted!" << endl;
      exit(1);
   }
   for (unsigned int i = 0; i<blow.size(); i++) {
      if (idiff == 1) { // Point-wise differential
         if (blow[i] != bupp[i]) {
            error["SetBinning1D"] << "Illegal binning, lower and upper bin edges must be equal for point-wise differential binnings, aborted!" << endl;
            error["SetBinning1D"] << "i = " << i << ", lower edge = " << blow[i] <<
               ", and upper edge = " << bupp[i] << endl;
            exit(1);
         }
      } else { // Non- and bin-wise differential
         if (blow[i] >= bupp[i]) {
            error["SetBinning1D"] << "Illegal binning, lower bin edge larger or equal to upper one, aborted!" << endl;
            error["SetBinning1D"] << "i = " << i << ", lower edge = " << blow[i] <<
               ", and upper edge = " << bupp[i] << endl;
            exit(1);
         }
         if (i<blow.size()-1 && bupp[i] > blow[i+1]) {
            error["SetBinning1D"] << "Illegal binning, overlapping bins, aborted!" << endl;
            error["SetBinning1D"] << "i = " << i << ", upper edge = " << bupp[i] <<
               ", and lower edge of i+1 = " << blow[i+1] << endl;
            exit(1);
         }
      }
   }

   // Check input: Normalization factor must be larger than zero for all bins
   if ( blow.size() != vnorm.size() ) {
      error["SetBinning1D"] << "Unequal size of bin border and normalization vectors, blow.size() = " << blow.size() <<
         ", vnorm.size() = " << vnorm.size() << ", aborted!" << endl;
      exit(1);
   }
   for (unsigned int i = 0; i<vnorm.size(); i++) {
      if ( vnorm[i] < DBL_MIN ) {
         error["SetBinning1D"] << "Normalization factor is too small, aborted!" << endl;
         error["SetBinning1D"] << "i = " << i << ", vnorm[i] = " << vnorm[i] << endl;
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
   info["SetBinning1D"] << "VVSIV: Set binning successfully for " << NObsBin << " bins in " << NDim << "dimensions." << endl;
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
   info["SetBinningND"] << "VIV: Set binning via vector with bin edges." << endl;
}


//
// Set continuous n-dimensional binning via a vector of vectors with bin edges (pairs)
//
void fastNLOCreate::SetBinningND(vector<vector<double> > bgrid, unsigned int ndim, vector<int> idiff) {

   // Check input: Dimensions, labels and idiff settings
   if ( ndim < 1 ) {
      error["SetBinningND"] << "Illegal # of dimensions, aborted!" << endl;
      error["SetBinningND"] << "ndim = " << ndim << endl;
      exit(1);
   } else if ( ndim > 3 ) {
      error["SetBinningND"] << "More than three dimensions not yet supported, aborted!" << endl;
      error["SetBinningND"] << "ndim = " << ndim << endl;
      exit(1);
   }
   if ( ndim != idiff.size() ) {
      error["SetBinningND"] << "Inconsistency between # of dimensions and vector length for differentiality settings, aborted!" << endl;
      error["SetBinningND"] << "Number of dimensions ndim = " << ndim << ", # of idiff settings = " << idiff.size() << endl;
      exit(1);
   }

   // Check input: Type of differential binning
   for ( unsigned int i=0; i<ndim; i++ ) {
      if ( idiff[i] < 0 || idiff[i] > 2 ) {
         error["SetBinningND"] << "Illegal choice of differentiality, idim = " << i << ", idiff = " << idiff[i] << ", aborted!" << endl;
         error["SetBinningND"] << "Each dimension must either be non-differential (idiff=0), " << endl;
         error["SetBinningND"] << "point-wise differential (idiff=1), or bin-wise differential (idiff=2)." << endl;
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
      debug["SetBinningND"] << "i = " << i << ", iObsBin = " << NObsBin << ", ishift = " << ishift << endl;

      // Loop over 'n' 'histogram' bins
      for (unsigned int j = ishift ; j<bgrid[i].size()-1; j++) {
         Bin.push_back(vector<pair<double,double> >(ndim));

         // Keep track of additional increment in vector for non- or bin-wise differential of the 'n-1' outer dimensions
         int kshift = 0;
         for ( unsigned int k = 0; k<ndim-1; k++ ) {
            debug["SetBinningND"] << "i = " << i << ", j = " << j << ",k+kshift = " << k+kshift << ", bgrid i,k+kshift = " << bgrid[i][k+kshift] << endl;
            debug["SetBinningND"] << "i = " << i << ", j = " << j << ",k+kshift+1 = " << k+kshift+1 << ", bgrid i,k+kshift+1 = " << bgrid[i][k+kshift+1] << endl;
            if ( idiff[k] == 1 ) { // Point-wise differential --> one bin center, set bin edges identical
               warn["SetBinningND"] << "Point-wise differential binning implemented, but not tested extensively. Please check carefully!" << endl;
               Bin[NObsBin][k] = make_pair(bgrid[i][k+kshift],bgrid[i][k+kshift]);
            } else {               // Non- or bin-wise differential --> two bin edges
               if ( !(bgrid[i][k+kshift] < bgrid[i][k+kshift+1]) ) {
                  error["SetBinningND"] << "Illegal binning, lower bin edge larger or equal to upper one, aborted!" << endl;
                  error["SetBinningND"] << "i, k = " << i << ", " << k << ", lower edge = " << bgrid[i][k+kshift] <<
                     ", and upper edge = " << bgrid[i][k+kshift+1] << endl;
                  exit(1);
               }
               Bin[NObsBin][k] = make_pair(bgrid[i][k+kshift],bgrid[i][k+kshift+1]);
               kshift++;
            }
         }

         if ( idiff[ndim-1] == 1 ) { // Point-wise differential --> one bin center
            warn["SetBinningND"] << "Point-wise differential binning implemented, but not tested extensively. Please check carefully!" << endl;
            Bin[NObsBin][ndim-1] = make_pair(bgrid[i][j],bgrid[i][j]);
         } else {                    // Non- or bin-wise differential --> two bin edges
            debug["SetBinningND"] << "i = " << i << ", j = " << j << ", bgrid i,j = " << bgrid[i][j] << endl;
            debug["SetBinningND"] << "i = " << i << ", j = " << j << ", bgrid i,j+1 = " << bgrid[i][j+1] << endl;
            if ( !(bgrid[i][j] < bgrid[i][j+1]) ) {
               error["SetBinningND"] << "Illegal binning, lower bin edge larger or equal to upper one, aborted!" << endl;
               error["SetBinningND"] << "i, j = " << i << ", " << j << ", lower edge = " << bgrid[i][j] <<
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
            error["SetBinningND"] << "Illegal binning, lower bin edge larger than upper one for outermost dimension, aborted!" << endl;
            error["SetBinningND"] << "The problematic observable bins are " << i << " and " << i+1 << endl;
            error["SetBinningND"] << "The bin edges for observable bin " << i << " are: Bin[i][0].first = " << Bin[i][0].first << ", Bin[i][0].second = " << Bin[i][0].second << endl;
            error["SetBinningND"] << "The bin edges for observable bin " << i+1 << " are: Bin[i+1][0].first = " << Bin[i+1][0].first << ", Bin[i+1][0].second = " << Bin[i+1][0].second << endl;
            exit(1);
         }
         bool lordered = false;
         for ( unsigned int j=0; j<ndim; j++ ) {
            debug["SetBinningND"] << "Before: iobs = " << i << ", jdim = " << j << ", lordered = " << lordered << ", Bin i,j = " << Bin[i][j].first << ", Bin i+1,j = " << Bin[i+1][j].first << endl;
            //            if ( !lordered && Bin[i][j].first < Bin[i+1][j].first ) {lordered = true;}
            if ( !lordered && (Bin[i+1][j].first - Bin[i][j].first) > DBL_MIN ) {lordered = true;}
            debug["SetBinningND"] << "After : iobs = " << i << ", jdim = " << j << ", lordered = " << lordered << ", Bin i,j = " << Bin[i][j].first << ", Bin i+1,j = " << Bin[i+1][j].first << endl;
         }
         if (! lordered) {
            error["SetBinningND"] << "Illegal binning, lower bin edge larger or equal to upper one, aborted!" << endl;
            error["SetBinningND"] << "The problematic observable bins are " << i << " and " << i+1 << endl;
            for ( unsigned int j=0; j<ndim; j++ ) {
               error["SetBinningND"] << "The bin edges for observable bin " << i << " in dimension " << j << " are: Bin[i][j].first = " << Bin[i][j].first << ", Bin[i][j].second = " << Bin[i][j].second << endl;
               error["SetBinningND"] << "The bin edges for observable bin " << i+1 << " in dimension " << j << " are: Bin[i+1][j].first = " << Bin[i+1][j].first << ", Bin[i+1][j].second = " << Bin[i+1][j].second << endl;

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

   info["SetBinningND"] << "Set binning successfully for " << NObsBin << " bins in " << NDim << " dimensions." << endl;
}


/////////////////////////////////////////////


// ___________________________________________________________________________________________________
void fastNLOCreate::GetWarmupValues() {
   //!
   //! GetWarmupValues.
   //! Checks if warmup-table exists and initialized
   //! member variable fIsWarmup
   //!
   debug["GetWarmupValues"]<<endl;

   std::cout.setstate(std::ios::failbit) ; // no cout in the following
   std::cerr.setstate(std::ios::failbit) ; // no cout in the following
   info>>"\n";
   info>> (_SSEP40+_SSEP40+_SSEP40) << endl;
   info["GetWarmupValues"]<<"Trying to get warmup values. Please ignore following messages from parser."<<endl;
   // try to get warmup values
   vector<vector<double> > warmup = DOUBLE_TAB_NS(Warmup.Values,fSteerfile);
   fIsWarmup = warmup.empty();

   // try again, with hard-coded convention:
   if (fIsWarmup) {
      debug["GetWarmupValues"]<<"Could not get warmup table from steerfile. Now trying to read steerfile: "<<GetWarmupTableFilename()<<endl;
      READ_NS(GetWarmupTableFilename(),fSteerfile);    // put the warmup-values into same read_steer 'namespace'
      warmup = DOUBLE_TAB_NS(Warmup.Values,fSteerfile);
      fIsWarmup = warmup.empty();
      if (!fIsWarmup)
         info["GetWarmupValues"]<<"Warmup values found in file "<<GetWarmupTableFilename()<<"."<<endl;
   }

   // inform user about success
   info>> (_SSEP40+_SSEP40+_SSEP40) << endl;
   std::cout.clear() ; // recover cout to screen
   std::cerr.clear() ; // recover cout to screen
   info["GetWarmupValues"]<<"This will be a "<<(fIsWarmup?"warmup":"production")<<" run."<<endl;
}



// ___________________________________________________________________________________________________
void fastNLOCreate::UseBinGridFromWarmup() {
   //! initialialize all binning related variables
   //! with values stored in the warmup file.
   vector<vector<double> > warmup =  DOUBLE_TAB_NS(Warmup.Binning,fSteerfile);
   NObsBin      = warmup.size();
   NDim         = INT_NS(Warmup.DifferentialDimension,fSteerfile);
   if (warmup[0].size() != (7+2*NDim) && warmup[0].size() != (5+2*NDim)) {
      error["UseBinGridFromWarmup"]<<"This warmup table has an unknown size of columns. Expecting "<<(7+2*NDim)<<" for flexible-scale, or "<<(5+2*NDim)<<" for fixed-scale tables. Exiting."<<endl;
      exit(1);
   }
   fIsFlexibleScale = (warmup[0].size() == (7+2*NDim));
   IDiffBin     = INT_ARR_NS(Warmup.DimensionIsDifferential,fSteerfile);
   DimLabel     = STRING_ARR_NS(Warmup.DimensionLabels,fSteerfile);

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

   vector<vector<double> > warmup =  DOUBLE_TAB_NS(Warmup.Values,fSteerfile);
   vector<vector<double> > wrmbin =  DOUBLE_TAB_NS(Warmup.Binning,fSteerfile);
   bool ret = true;

   const string wrmuphelp = "Please remove warmup-file in order to calculate a new warmup-file which is compatible to your steering,\nor alternatively use 'ReadBinningFromSteering=false', then all binning-related information is taken from the warmup file.\n";

   if (warmup.size() != NObsBin) {
      error["CheckWarmupConsistency"]
            <<"Table of warmup values is not compatible with steering file.\n"
            <<"Different number of bins ("<<warmup.size()<<" instead of "<<NObsBin<<".\n"
            <<wrmuphelp
            <<"Exiting."<<endl;
      ret = false;
      exit(1);
   }
   if (INT_NS(Warmup.DifferentialDimension,fSteerfile) != (int)NDim) {
      error["CheckWarmupConsistency"]
            <<"Table of warmup values is not compatible with steering file.\n"
            <<"Found different number of dimensions. NDim="<<NDim<<", Warmup.DifferentialDimension="<<INT_NS(Warmup.DifferentialDimension,fSteerfile)<<".\n"
            <<wrmuphelp
            <<"Exiting."<<endl;
      ret = false;
      exit(1);
   }

   // CoeffTable is not available during intialization
   //    if ( STRING_NS(Warmup.ScaleDescriptionScale1) != GetTheCoeffTable()->ScaleDescript[0][0] ) {
   //       warn["CheckWarmupConsistency"]
   //    <<"Table of warmup values is potentially incompatible with steering file.\n"
   //    <<"Found different scale description (ScaleDescriptionScale1). ScaleDescriptionScale1="<<GetTheCoeffTable()->ScaleDescript[0][0]
   //    <<", but Warmup.ScaleDescriptionScale1='"<<STRING_NS(Warmup.ScaleDescriptionScale1)<<"'."<<endl<<endl;
   //       ret = false;
   //    }
   //    if ( fIsFlexibleScale && STRING_NS(Warmup.ScaleDescriptionScale2) != GetTheCoeffTable()->ScaleDescript[0][1] ){
   //       warn["CheckWarmupConsistency"]
   //    <<"Table of warmup values is potentially incompatible with steering file.\n"
   //    <<"Found different scale description (ScaleDescriptionScale2). ScaleDescriptionScale2='"<<GetTheCoeffTable()->ScaleDescript[0][0]<<"'"
   //    <<", but Warmup.ScaleDescriptionScale2='"<<STRING_NS(Warmup.ScaleDescriptionScale1)<<"'."<<endl<<endl;
   //       ret = false;
   //    }

   if (INT_ARR_NS(Warmup.DimensionIsDifferential,fSteerfile)[0] != IDiffBin[0]) {
      warn["CheckWarmupConsistency"]
            <<"Table of warmup values seems to be incompatible with steering file.\n"
            <<"Found different diff-label for dimension 0  (IDiffBin). DimensionIsDifferential='"<<IDiffBin[0]<<"'"
            <<", but Warmup.DimensionIsDifferential[0]='"<<DOUBLE_ARR_NS(Warmup.DimensionIsDifferential,fSteerfile)[0]<<"'. Exiting."<<endl;
      ret = false;
      exit(1);
   }
   if (NDim > 1 && INT_ARR_NS(Warmup.DimensionIsDifferential,fSteerfile)[1] != IDiffBin[1]) {
      warn["CheckWarmupConsistency"]
            <<"Table of warmup values seems to be incompatible with steering file.\n"
            <<"Found different diff-label for dimension 0  (IDiffBin). DimensionIsDifferential='"<<IDiffBin[1]<<"'"
            <<", but Warmup.DimensionIsDifferential[0]='"<<DOUBLE_ARR_NS(Warmup.DimensionIsDifferential,fSteerfile)[1]<<"'. Exiting."<<endl;
      ret = false;
      exit(1);
   }

   // check binning in detail
   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      const int i0 = 1;//fIsFlexibleScale ? 6 : 4;
      if (NDim == 1) {
         if (Bin[i][0].first != wrmbin[i][i0] || Bin[i][0].second != wrmbin[i][i0+1]) {
            error["CheckWrmbinConsistency"]
                  <<"Table of warmup values seems to be incompatible with steering file.\n"
                  <<"Found different binning for bin "<<i<<", steering: ["<<Bin[i][0].first<<","<<Bin[i][0].second
                  <<",], warmup: ["<<wrmbin[i][i0]<<","<<wrmbin[i][i0+1]<<"].\n"
                  <<wrmuphelp
                  <<"Exiting."<<endl;
            ret = false;
            exit(1);
         }
      } else if (NDim == 2) {
         if (Bin[i][0].first != wrmbin[i][i0] || Bin[i][0].second != wrmbin[i][i0+1]
               || Bin[i][1].first != wrmbin[i][i0+2] || Bin[i][1].second != wrmbin[i][i0+3]) {
            error["CheckWarmupConsistency"]
                  <<"Table of warmup values seems to be incompatible with steering file.\n"
                  <<"Found different binning for bin "<<i<<", steering: ["<<Bin[i][0].first<<","<<Bin[i][0].second<<",] ["<<Bin[i][1].first<<","<<Bin[i][1].second
                  <<"], warmup: ["<<wrmbin[i][i0]<<","<<wrmbin[i][i0+1]<<"] ["<<wrmbin[i][i0+2]<<","<<wrmbin[i][i0+3]<<"].\n"
                  <<wrmuphelp
                  <<"Exiting."<<endl;
            ret = false;
            exit(1);
         }
      }
      //check bin width
      double bwwrm = 0;
      if (NDim == 1) bwwrm = wrmbin[i][i0+2];
      else if (NDim == 2) bwwrm = wrmbin[i][i0+4];
      else if (NDim == 3) bwwrm = wrmbin[i][i0+6];
      if (fabs(BinSize[i] - bwwrm) > 1.e-6) {
         warn["CheckWarmupConsistency"]
               <<"Table of warmup values seems to be incompatible with steering file.\n"
               <<"Found different bin size for bin "<<i<<". Steering: "<<BinSize[i]
               <<", warmup: "<<bwwrm<<".\n"
               <<wrmuphelp
               <<"Please check consistency!"<<endl<<endl;
         ret = false;
      }
   }
   return ret;
}



// ___________________________________________________________________________________________________
void fastNLOCreate::SetOrderOfAlphasOfCalculation(unsigned int ord) {
   info["SetOrderOfAlphasOfCalculation"] << "Base order ord = " << ord << endl;
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
      info["SetOrderOfAlphasOfCalculation"] << "LO scale dependence: MuR changeable; independent of MuF. IScaleDep = " << c->IScaleDep << endl;
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
      info["SetOrderOfAlphasOfCalculation"] << "NLO scale dependence: Dependent on MuR and MuF. IScaleDep = " << c->IScaleDep << endl;
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
         error["SetOrderOfAlphasOfCalculation"]<<"Unknown order of perturbation theory: order="<<ord-GetLoOrder()<<" (ord="<<ord<<",ILOord="<<ILOord<<"). Exiting."<<endl;
         exit(1);
      }
   }
   info["SetOrderOfAlphasOfCalculation"] << "Using " << c->NSubproc <<
      " subprocesses and PDF flags: " << c->IPDFdef1 << ", " << c->IPDFdef2 << ", " << c->IPDFdef3 << "." << endl;

   // init array with counter processes (symmetric and asymmetric ones)
   if (c->NPDFPDG.size() == 2 && c->NPDFDim == 1) {
      fSymProc.resize(c->NSubproc);
      for (int p = 0 ; p<c->NSubproc ; p++) fSymProc[p]=p;

      for (unsigned int i = 0 ; i<fProcConsts.AsymmetricProcesses.size() ; i++) {
         if (fProcConsts.AsymmetricProcesses[i].first<c->NSubproc) {   // safety
            if (fProcConsts.AsymmetricProcesses[i].second >= GetNSubprocesses() || fSymProc[fProcConsts.AsymmetricProcesses[i].first] >= GetNSubprocesses()) {
               if (!(c->IPDFdef1==3&&c->IPDFdef2==1&&c->IPDFdef3==1))   // it is normal in pp->jets in LO
                  warn["SetOrderOfAlphasOfCalculation"]<<"Subprocess "<<fSymProc[fProcConsts.AsymmetricProcesses[i].first]<<" is requested to be asymmetric with subprocess "<<fProcConsts.AsymmetricProcesses[i].second<<", but there are only "<<GetNSubprocesses()<<" subprocesses in this calculation. Ignoring call."<<endl;
            }
            else
               fSymProc[fProcConsts.AsymmetricProcesses[i].first] = fProcConsts.AsymmetricProcesses[i].second;
         }
      }

      //       vector<vector<int> > asym = INT_TAB_NS(AsymmetricProcesses);
      //       for ( unsigned int i = 0 ; i<asym.size() ; i ++ ) {
      //         if ( asym[i][0]<(int)fSymProc.size() ) { // safety
      //                    if ( asym[i][1] >= GetNSubprocesses() || fSymProc[asym[i][0]] >= GetNSubprocesses() ) {
      //                       warn["AsymmetricProcesses"]<<"Subprocess "<<fSymProc[asym[i][0]]<<" is requested to be asymmetric with subprocess "<<asym[i][1]<<", but there are only "<<GetNSubprocesses()-1<<" subprocesses in this calculation. Ignoring call."<<endl;
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
   if (IOrdInitial!= fIOrd && !fIsWarmup) {
      if (!fIsFlexibleScale) InitInterpolationKernels();
      InitGrids();
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::SetLoOrder(int LOOrd) {
   debug["SetLoOrder"]<<endl;
   fastNLOTable::SetLoOrder(LOOrd);
   if (fIsFlexibleScale)
      ((fastNLOCoeffAddFlex*)GetTheCoeffTable())->fILOord = LOOrd;
}



// ___________________________________________________________________________________________________
void fastNLOCreate::InitVariablesInCoefficientTable() {
   debug["InitVariablesInCoefficientTable"]<<endl;
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
   debug["ReadCoefficientSpecificVariables"]<<endl;
   // todo: make it more user friendly
   // todo: include some sanity checks
   fastNLOCoeffAddBase* c = GetTheCoeffTable();

   // generator constants
   c->CodeDescript      = fGenConsts.GetCodeDescription();
   for (unsigned int i = 0 ; i<fProcConsts.GetProcessDescription().size() ; i++) {
      c->CodeDescript.push_back(fProcConsts.GetProcessDescription()[i]);
   }
   c->SetIXsectUnits(fProcConsts.UnitsOfCoefficients);

   // (some) process constants
   c->NPDFPDG.resize(fProcConsts.NPDF);
   if (c->NPDFPDG.size() >0) c->NPDFPDG[0] = INT_NS(PDF1,fSteerfile);   // from steering
   if (c->NPDFPDG.size() >1) c->NPDFPDG[1] = INT_NS(PDF2,fSteerfile);   // from steering
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
      c->NScaleDep              = 3; // temporaily. Until known if generator runs in LO, NLO or NNLO.
      c->ScaleDescript[0].resize(2);
      if (BOOL_NS(ReadBinningFromSteering,fSteerfile)) {
         c->ScaleDescript[0][0] = STRING_NS(ScaleDescriptionScale1,fSteerfile);
         c->ScaleDescript[0][1] = STRING_NS(ScaleDescriptionScale2,fSteerfile);
      } else {
         c->ScaleDescript[0][0] = STRING_NS(Warmup.ScaleDescriptionScale1,fSteerfile);
         c->ScaleDescript[0][1] = STRING_NS(Warmup.ScaleDescriptionScale2,fSteerfile);
      }
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
      if (BOOL_NS(ReadBinningFromSteering,fSteerfile))
         c->ScaleDescript[0][0] = STRING_NS(ScaleDescriptionScale1,fSteerfile);
      else
         c->ScaleDescript[0][0] = STRING_NS(Warmup.ScaleDescriptionScale1,fSteerfile);
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
   if (idiff == 1) {
      if (fLastScen._o[0] == fScenario._o[0]) return fObsBin;
   } else if (idiff == 2) {
      if (fLastScen._o[0] == fScenario._o[0] && fLastScen._o[1] == fScenario._o[1])  return fObsBin;
   } else if (idiff == 3) {
      if (fLastScen._o[0] == fScenario._o[0] && fLastScen._o[1] == fScenario._o[1] && fLastScen._o[2] == fScenario._o[2])  return fObsBin;
   } else {
      error["GetBin"] << "More than triple-differential binning not yet implemented, aborted!" << endl;
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
      error["GetBin"] << "More than triple-differential binning not yet implemented, aborted!" << endl;
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

      //if ( fEvent._w == 0 ) return; // nothing todo.

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


      v4d& st = c->SigmaTilde[ObsBin];
      for (unsigned int is = 0 ; is<events.size() ; is++) {
         double mu = fScenario._m1 * fScaleFac[is];
         const vector<pair<int,double> >& nmu  = fKernMuS[ObsBin][is]->GetNodeValues(mu);
         for (unsigned int m1 = 0 ; m1<nmu.size() ; m1++) {
            v2d& stm1 = st[is][nmu[m1].first];
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
      error["FillAllSubprocess"]<<"This table expects "<<GetNSubprocesses()<<" subprocesses, but only "<<events.size()<<" are provided. Exiting."<<endl;
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
void fastNLOCreate::Fill(int scalevar) {
   //!
   //! Fill values, which are stored in 'Event' and 'Scenario' into fastNLO table.
   //!
   //debug["Fill"]<<"Filling subprocess contributions into table."<<endl;

   //GetTheCoeffTable()->Nevt++; // todo: counting of events must be properly implemented
   fStats._nProc++; //keep statistics

   if (fIsWarmup && scalevar==0) UpdateWarmupArrays();
   else FillContribution(scalevar);

   fEvent.ResetButX();
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillContribution(int scalevar) {
   //! read information from 'Event' and 'Scenario'
   //! do the interpolation
   //! and fill into the tables.


   if (fEvent._n > 0) SetNumberOfEvents(fEvent._n);

   const int ObsBin = (fScenario._iOB == -1) ? GetBin() : fScenario._iOB;
   // hier
//    cout<<"Fill! ObsBin="<<ObsBin<<", w="<<fEvent._w<<", x1="<<fEvent._x1<<", x2="<<fEvent._x2<<", o0="<<fScenario._o[0]<<endl;
//    static double wsum = 0;
//    wsum+= fEvent._w; //hier
//    cout<<" * wSum = "<<wsum<<endl;

   if (ObsBin < 0) return;
   if (ObsBin >= (int)GetNObsBin()) return;
   fStats._nEvPS++;

   fastNLOCoeffAddBase* c = GetTheCoeffTable();

   int p = fEvent._p;
   if ( p<0 || p > c->GetNSubproc() ) {
      error["FillContributionFixHHC"]<<"Unknown process id p="<<p<<endl;
      exit(1);
   }

   // ---- DIS ---- //
   if (c->GetNPDF() == 1 && fastNLOCoeffAddFlex::CheckCoeffConstants(c,true)) {
      // todo
      FillContributionFlexDIS((fastNLOCoeffAddFlex*)GetTheCoeffTable(),  ObsBin);
      //{error["FillContribution"]<<"Don't know how to fill this table. Exiting."<<endl; exit(1); }
   } else if (c->GetNPDF() == 1 && fastNLOCoeffAddFix::CheckCoeffConstants(c,true)) {
      // todo
      {error["FillContribution"]<<"Don't know how to fill this table. Exiting."<<endl; exit(1); }
   }
   // ---- pp/ppbar ---- //
   else if (c->GetNPDF() == 2 && fastNLOCoeffAddFlex::CheckCoeffConstants(c,true))
      FillContributionFlexHHC((fastNLOCoeffAddFlex*)GetTheCoeffTable(),  ObsBin);
   else if (c->GetNPDF() == 2 && fastNLOCoeffAddFix::CheckCoeffConstants(c,true))
      FillContributionFixHHC((fastNLOCoeffAddFix*)GetTheCoeffTable(),  ObsBin, scalevar);
   else {
      error["FillContribution"]<<"Don't know how to fill this table. Exiting."<<endl;
      exit(1);
   }
}



// ___________________________________________________________________________________________________
void fastNLOCreate::FillContributionFixHHC(fastNLOCoeffAddFix* c, int ObsBin, int scalevar) {
   //! read informatio from 'Event' and 'Scenario'
   //! do the interpolation
   //! and fill into the tables.
   debug["FillContributionFixHHC"]<<endl;

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
   if (CheckWeightIsNan()) return;
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
   debug["FillContributionFlexHHC"]<<endl;

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
   if (CheckWeightIsNan()) return;
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
               if (std::isnan(wfnlo)) {
                  error[""]<<"wfnlo is a nan."<<endl;
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
   debug["FillContributionFlexHHC"]<<endl;

   if (fEvent._w == 0 && fEvent._wf==0 && fEvent._wr==0) return;   // nothing todo.

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
   if (CheckWeightIsNan()) return;
   for (unsigned int ix = 0 ; ix<nx.size() ; ix++) {
      int p = fEvent._p;
      int xIdx = nx[ix].first;
      //HalfMatrixCheck(xminbin,xmaxbin,p);
      //int ixHM = GetXIndex(ObsBin,xminbin,xmaxbin);

      for (unsigned int m1 = 0 ; m1<nmu1.size() ; m1++) {
         for (unsigned int mu2 = 0 ; mu2<nmu2.size() ; mu2++) {
            double wfnlo = nx[ix].second * nmu1[m1].second * nmu2[mu2].second / BinSize[ObsBin];
            if (std::isnan(wfnlo)) {
               error[""]<<"wfnlo is a nan."<<endl;
               fKernX1[ObsBin]->PrintGrid();
               fKernMu1[ObsBin]->PrintGrid();
               fKernMu2[ObsBin]->PrintGrid();
               cout<<"ix1="<<ix<<", im1="<<m1<<", im2="<<mu2<<endl;
               cout<<"x1="<<nx[ix].second<<", ix="<<ix<<", xval="<<x<<endl;
               cout<<"m1="<< nmu1[m1].second<<", m1="<<m1<<", mu1val="<<fScenario._m1<<endl;
               cout<<"m2="<<nmu2[mu2].second<<", m2="<<mu2<<", mu2val="<<fScenario._m2<<endl;
               exit(1);
            }
            if (fEvent._w  != 0) {
               //                 cout<<"   Fill * : ix="<<ixHM<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._w  * wfnlo<<endl;
               c->SigmaTildeMuIndep[ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._w  * wfnlo;
            }
            if (fEvent._wf != 0) {
               //                 cout<<"   Fill F : ix="<<ixHM<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._wf  * wfnlo<<endl;
               c->SigmaTildeMuFDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wf * wfnlo;
            }
            if (fEvent._wr != 0) {
               //                 cout<<"   Fill R : ix="<<ixHM<<", im1="<<nmu1[m1].first<<", im2="<<nmu2[mu2].first<<", p="<<p<<", w="<<fEvent._wr  * wfnlo<<endl;
               c->SigmaTildeMuRDep [ObsBin][xIdx][nmu1[m1].first][nmu2[mu2].first][p]  += fEvent._wr * wfnlo;
            }
         }
      }
   }
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
         //             error["HalfMatrixCheck"]<<"Necessary array with symmetric processes for half-matrix notation not initialized."<<endl;

         //cout<<"exchange supbrpc. xminbin="<<xminbin<<", xmaxbin="<<xmaxbin<<", p="<<subproc<<", pAsym="<<fSymProc[subproc]<<endl;
         int di = xminbin - xmaxbin;
         xmaxbin += di;        // modify indicees
         xminbin -= di;
         subproc = fSymProc[subproc];           // exchange asymmetric process
      }
   }
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::CheckWeightIsNan() {
   //! check if weights contain isnan
   if (std::isnan(fEvent._w)) {
      error["CheckWeightIsNan"]<<"(Scale-independent) weight is 'nan'"<<endl;
      return true;
   }
   if (std::isnan(fEvent._wf)) {
      error["CheckWeightIsNan"]<<"Factorization scale dependent weight is 'nan'"<<endl;
      return true;
   }
   if (std::isnan(fEvent._wr)) {
      error["CheckWeightIsNan"]<<"Renormalization scale dependent weight is 'nan'"<<endl;
      return true;
   }
   return false;
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
      if (!xGrid2) error["GetNxmax"]<<"Error. Second x-grid must be specified."<<endl;
      if (xGrid1->size() != xGrid2->size())error["GetNxmax"]<<"Grid sizes in half-matrix notation must have equal size."<<endl;
      return ((int)pow((double)xGrid1->size(),2)+xGrid1->size())/2;
   case 2:
      if (!xGrid2) error["GetNxmax"]<<"Error. Second x-grid must be specified."<<endl;
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
      if (wgtn==0) {error["ApplyPDFWeight"]<<"Cannot divide by 0."<<endl; exit(1);}
      //cout<<"in="<<in<<" wgtx="<<wgtx<<", x="<<x<<", wgtn="<<wgtn<<", nod="<<grid->at(nodes[in].first)<<endl;
      nodes[in].second *= wgtx/wgtn;
   }
}


// ___________________________________________________________________________________________________
inline double fastNLOCreate::CalcPDFReweight(double x) const {
   if (x<=0) { error["CalcPDFReweight"]<<"Cannot calculate sqrt of negative numbers or divide by zero. x="<<x<<endl; exit(1);}
   double w=(1.-0.99*x)/sqrt(x);
   return w*w*w;
}


// ___________________________________________________________________________________________________
void fastNLOCreate::NormalizeCoefficients() {
   //! Set number of events to 1 and weight coefficients in sigmatilde
   //! accordingly
   //! This means, that the information about the
   //! number of events is essentially lost
   GetTheCoeffTable()->NormalizeCoefficients();
   fStats._nEv=1;
   //    double nev = GetTheCoeffTable()->GetNevt(0,0);
   //    MultiplyCoefficientsByConstant(1./nev);
   //    SetNumberOfEvents(1.);
}


// ___________________________________________________________________________________________________
void fastNLOCreate::MultiplyCoefficientsByBinSize() {
   //! Multiply all coefficients by binsize
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
   debug["UpdateWarmupArrays"]<<"ObsBin="<<ObsBin<<"\tmu1="<<fScenario._m1<<"\tmu2="<<fScenario._m2<<"\tx1="<<fEvent._x1<<"\tx2="<<fEvent._x2<<endl;
   if (ObsBin >= 0) {
      fWMu1[ObsBin].first       = std::min(fScenario._m1,fWMu1[ObsBin].first) ;
      fWMu1[ObsBin].second      = std::max(fScenario._m1,fWMu1[ObsBin].second) ;
      if (GetTheCoeffTable()->IPDFdef1 == 3) { // pp/ppbar
         fWx[ObsBin].first      = std::min(std::min(fEvent._x1,fEvent._x2),fWx[ObsBin].first) ;
         fWx[ObsBin].second     = std::max(std::max(fEvent._x1,fEvent._x2),fWx[ObsBin].second) ;
      } else if (GetTheCoeffTable()->IPDFdef1 == 2) {  // DIS
         fWx[ObsBin].first      = std::min(fEvent._x1,fWx[ObsBin].first) ;
         fWx[ObsBin].second     = std::max(fEvent._x1,fWx[ObsBin].second) ;
      } else
         error["UpdateWarmupArrays"]<<"nothing reasonable implemented yet."<<endl;
      if (fIsFlexibleScale) {
         fWMu2[ObsBin].first    = std::min(fScenario._m2,fWMu2[ObsBin].first) ;
         fWMu2[ObsBin].second   = std::max(fScenario._m2,fWMu2[ObsBin].second) ;
      }
   }
}


// ___________________________________________________________________________________________________
void fastNLOCreate::InitWarmupArrays() {
   debug["InitWarmupArrays"]<<endl;
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
   //if ( GetTheCoeffTable()->GetNevt(0,0) <= 0 ) {
   if (GetTheCoeffTable()->Nevt <= 0) {
      warn["WriteTable"]<<"Number of events seems to be not filled. Please use SetNumberOfEvents(int) before writing table."<<endl;
   }
   fStats.PrintStats();
   if (fIsWarmup) {
      info["WriteTable"]<<"Writing warmup table instead of coefficient table."<<endl;
      // round warmup values and try to guess bin boundaries
      AdjustWarmupValues();
      // write table to disk
      WriteWarmupTable();
   } else {
      if (ffilename == "") {
         error["WriteTable"]<<"No filename given."<<endl;
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
   WriteTable();
}


// ___________________________________________________________________________________________________
void fastNLOCreate::WriteWarmupTable() {
   string tempfn = ffilename;
   string warmupfile = GetWarmupTableFilename();
   info["WriteWarmupTable"]<<"Writing warmup table to: "<<warmupfile<<endl;
   SetFilename(warmupfile);

   // open stream;
   ofstream* table = OpenFileWrite();
   // write to disk
   OutWarmup(*table);
   // close file
   table->close();
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
   if (fWxRnd.empty()) {
      warn["OutWarmup"]<<"Warmup arrays not initialized. Did you forgot to fill values?"<<endl;
//       warn["OutWarmup"]<<"  Continuting, but writing unreasonalby large/small values as warmup values..."<<endl;
//       InitWarmupArrays();
      error["OutWarmup"]<<" Do not write out unreasonable warmup table. Exiting."<<endl;
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
   strm<<"# This file has been calculated using "<<GetTheCoeffTable()->Nevt<<" contributions."<<endl;
   strm<<"#   ( Mind: contributions != events. And contributions are not necessarily in phase space region."<<endl;
   strm<<"# Please check by eye for reasonability of the values."<<endl;
   strm<<" " <<endl;

   // write variables of warmup run
   strm<<"Warmup.OrderInAlphasOfWarmupRunWas\t"<<  fIOrd <<endl;
   strm<<"Warmup.CheckScaleLimitsAgainstBins\t"<<(BOOL_NS(CheckScaleLimitsAgainstBins,fSteerfile)?"true":"false")<<endl;
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

   // write readable table
   char buf[4000];
   static const double RoundXDown = 0.96; // to avoid rounding warning messages and have some 'x' in spare.
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
         if (fWxRnd[i].first < 1.e-6) {
            warn["OutWarmup"]<<"The xmin value in bin "<<i<<" seems to be unreasonably low (xmin="<<fWxRnd[i].first<<"). Taking xmin=1.e-6 instead."<<endl;
            fWxRnd[i].first=1.e-6;
         }
         sprintf(buf,"   %4d    %9.2e  %9.2e  %14.2f  %14.2f  %14.3f  %14.3f",
                 i,fWxRnd[i].first*RoundXDown,fWxRnd[i].second,fWMu1Rnd[i].first,fWMu1Rnd[i].second,fWMu2Rnd[i].first,fWMu2Rnd[i].second);
         strm<<buf<<endl;
      }
   } else {
      // is ScaleDescript available?
      if (GetTheCoeffTable()->ScaleDescript[0].empty()) { error["OutWarmup"]<<"Scale description is empty. but needed. Probably this has to be implemented."<<endl; exit(1);};
      // table header
      sprintf(buf,"   ObsBin   %9s  %9s  %16s  %16s",
              "x_min","x_max", GetWarmupHeader(0,"min").c_str(), GetWarmupHeader(0,"max").c_str());
      strm<<buf<<endl;
      // table values
      for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
         if (fWxRnd[i].first < 1.e-6) {
            warn["OutWarmup"]<<"The xmin value in bin "<<i<<" seems to be unreasonably low (xmin="<<fWxRnd[i].first<<"). Taking xmin=1.e-6 instead."<<endl;
            fWxRnd[i].first=1.e-6;
         }
         sprintf(buf,"   %4d     %9.2e  %9.2e  %16.2f  %16.2f",
                 i,fWxRnd[i].first*RoundXDown, fWxRnd[i].second, fWMu1Rnd[i].first, fWMu1Rnd[i].second);
         strm<<buf<<endl;
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
   sprintf(buf,"  %12s","BinSize");
   strm<<buf<<endl;

   // table values
   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      sprintf(buf,"    %4d ",i); // obsbin
      strm<<buf;
      for (unsigned int idim = 0 ; idim<NDim ; idim++) {
         sprintf(buf,"  %12.3f  %12.3f",Bin[i][idim].first , Bin[i][idim].second);
         strm<<buf;
      }
      sprintf(buf,"  %12.3f",BinSize[i]);
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
   string ret = fSteerfile;
   size_t pos = ret.find(".str");
   if (pos != std::string::npos) ret.erase(pos,4);
   pos = ret.find(".steer");
   if (pos != std::string::npos) ret.erase(pos,6);
   ret += "_";
   ret += GetScenName();
   ret += "_warmup.txt";
   return ret;
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

   //------------------------------------------
   // 1. Are warmup-values identical to bin-boundaries ?
   int ident1=-1, ident2=-1;
   //if ( BOOL_NS(CheckScaleLimitsAgainstBins ,fSteerfile) ) {
   if (fIsFlexibleScale) {
      ident1 = CheckWarmupValuesIdenticalWithBinGrid(fWMu1);
      ident2 = CheckWarmupValuesIdenticalWithBinGrid(fWMu2);
   } else {
      ident1 = CheckWarmupValuesIdenticalWithBinGrid(fWMu1);
   }
   //}

   // ---------------------------------------
   // 2. round values to 3rd digit if applicable
   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      fWMu1Rnd[i].first  = fWMu1[i].first;
      fWMu1Rnd[i].second = fWMu1[i].second;
   }
   if (ident1 == -1) {
      RoundValues(fWMu1Rnd,3);
   }
   if (fIsFlexibleScale) {
      for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
         fWMu2Rnd[i].first  = fWMu2[i].first;
         fWMu2Rnd[i].second = fWMu2[i].second;
      }
      if (ident2 == -1) {
         RoundValues(fWMu2Rnd,3);
      }
   }

   // ---------------------------------------
   // 3. 'round' lower x-values down
   const double xdown = 0.20;
   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      fWxRnd[i].first  = fWx[i].first;
      fWxRnd[i].second = fWx[i].second;
   }
   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      fWxRnd[i].first *= (1.-xdown);
   }

}



// ___________________________________________________________________________________________________
void fastNLOCreate::RoundValues(vector<pair<double,double> >& wrmmu, int nthdigit) {
   //! Round warmup values up (down) if third relevant
   //! digit is a 9 (0)
   //! lower values are only rounded down,
   //! upper values are only rounded up
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
   vector<int > nbinup(NDim);;
   for (int idim = int(NDim)-1 ; idim>=0 ; idim--) {
      for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
         if ( Bin[i][idim].first != 0 ) {
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
   //       warn["CheckWarmupValuesIdenticalWithBinGrid"]
   //          <<"Adjusted warmup values to bin boundaries in different observables.\n"
   //          <<"\t\tThis may yield unreasonable results. Please check warmup-table carefully!\n"<<endl;
   //    }
   //       }
   //    }
   // round all bins if applicable
   for (unsigned int idim = 0 ; idim<NDim ; idim++) {
      debug["CheckWarmupValuesIdenticalWithBinGrid"]<<"found nbinlo="<<nbinlo[idim]<<" and nbinup="<<nbinup[idim]<<endl;
      double frac= (nbinlo[idim]+nbinup[idim])/(2.*(int)GetNObsBin());
      if (frac>minallbins) {   // round all bins
         info["CheckWarmupValuesIdenticalWithBinGrid"]
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
   debug["InitGrids"]<<endl;
   if (fKernX1.empty()) error["InitGrids"]<<"Interpolation kernels must be initialized before calling this function."<<endl;
   if (fIsFlexibleScale) {
      fastNLOCoeffAddFlex* c = (fastNLOCoeffAddFlex*)GetTheCoeffTable();
      //       if ( (c->GetNPDF()==2 && c->GetNPDFDim() == 1) || (c->GetNPDF()==1)   ) {;} // ok!
      //       else {
      //         error["InitGrids"]<<"Only half-matrix or DIS implemented."<<endl; exit(1);
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
      //         //error["InitGrids"]<<"Only half-matrix is implemented for grids for fixed-scale tables."<<endl; exit(1);
      //       }

      int nscalevar = fScaleFac.size();
      if (nscalevar==0) {
         error["InitGrids"]<<"No scale factors found."<<endl;
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
   debug["InitInterpolationKernels"]<<endl;
   if (fIsWarmup) {
      error["InitInterpolationKernels"]<<"Interpolation kernels can only be initialized in production runs. Warmup values must be known."<<endl;
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
   // todo. clean up memory
   vector<double> wrmX = DOUBLE_COL_NS(Warmup.Values,x_min,fSteerfile);
   vector<double> wrmMu1Up, wrmMu1Dn;
   wrmMu1Dn = read_steer::getdoublecolumn("Warmup.Values",GetWarmupHeader(0,"min"),fSteerfile);
   wrmMu1Up = read_steer::getdoublecolumn("Warmup.Values",GetWarmupHeader(0,"max"),fSteerfile);
   if (wrmMu1Dn.size()!=GetNObsBin() || wrmMu1Up.size()!= GetNObsBin()) {
      error["InitInterpolationKernels"]<<"Could not read warmup values for Mu1. Exiting."<<endl;
      exit(1);
   }
   vector<double> wrmMu2Up, wrmMu2Dn;
   if (fIsFlexibleScale) {
      wrmMu2Dn = read_steer::getdoublecolumn("Warmup.Values",GetWarmupHeader(1,"min"),fSteerfile);
      wrmMu2Up = read_steer::getdoublecolumn("Warmup.Values",GetWarmupHeader(1,"max"),fSteerfile);
      if (wrmMu2Dn.size()!=GetNObsBin() || wrmMu2Up.size()!= GetNObsBin()) {
         error["InitInterpolationKernels"]<<"Could not read warmup values for Mu2. Exiting."<<endl;
         exit(1);
      }
   }

   int npdf = GetTheCoeffTable()->GetNPDF();

   for (unsigned int i = 0 ; i < GetNObsBin() ; i ++) {
      // ------------------------------------------------
      // init x-interpolation kernels
      // ------------------------------------------------
      debug["InitInterpolationKernels"]<<"Make x grid for obsbin="<<i<<endl;
      fKernX1[i] = MakeInterpolationKernels(STRING_NS(X_Kernel,fSteerfile),wrmX[i],1); // use 1 as upper x-value
      if (npdf == 2)
         fKernX2[i] = MakeInterpolationKernels(STRING_NS(X_Kernel,fSteerfile),wrmX[i],1);

      // Create x grids with X_NNodes+1 nodes up to x_max = 1.
      // The additional last node will be removed again below.
      int nxtot = INT_NS(X_NNodes,fSteerfile) + 1;
      if (BOOL_NS(X_NoOfNodesPerMagnitude,fSteerfile)) {
         fKernX1[i]->MakeGridsWithNNodesPerMagnitude(fastNLOInterpolBase::TranslateGridType(STRING_NS(X_DistanceMeasure,fSteerfile)),nxtot);
         if (npdf == 2)
            fKernX2[i]->MakeGridsWithNNodesPerMagnitude(fastNLOInterpolBase::TranslateGridType(STRING_NS(X_DistanceMeasure,fSteerfile)),nxtot);
      } else {
         fKernX1[i]->MakeGrids(fastNLOInterpolBase::TranslateGridType(STRING_NS(X_DistanceMeasure,fSteerfile)),nxtot);
         if (npdf == 2)
            fKernX2[i]->MakeGrids(fastNLOInterpolBase::TranslateGridType(STRING_NS(X_DistanceMeasure,fSteerfile)),nxtot);
      }

      // Remove last node at x = 1; is multiplied by PDFs equalling zero anyway.
      fKernX1[i]->RemoveLastNode();
      if (npdf == 2)
         fKernX2[i]->RemoveLastNode();

      // ------------------------------------------------
      // init scale1-interpolation kernels
      // ------------------------------------------------
      int nqtot1 = INT_NS(Mu1_NNodes,fSteerfile);
      if (fIsFlexibleScale) {
         debug["InitInterpolationKernels"]<<"Make Mu1 grid for obsbin="<<i<<endl;
         fKernMu1[i] = MakeInterpolationKernels(STRING_NS(Mu1_Kernel,fSteerfile),wrmMu1Dn[i],wrmMu1Up[i]);
         fKernMu1[i]->MakeGrids(fastNLOInterpolBase::TranslateGridType(STRING_NS(Mu1_DistanceMeasure,fSteerfile)),nqtot1);
         // ------------------------------------------------
         // init scale2-interpolation kernels
         // ------------------------------------------------
         int nqtot2 = INT_NS(Mu2_NNodes,fSteerfile);
         debug["InitInterpolationKernels"]<<"Make Mu2 grid for obsbin="<<i<<endl;
         fKernMu2[i] = MakeInterpolationKernels(STRING_NS(Mu2_Kernel,fSteerfile),wrmMu2Dn[i],wrmMu2Up[i]);
         fKernMu2[i]->MakeGrids(fastNLOInterpolBase::TranslateGridType(STRING_NS(Mu2_DistanceMeasure,fSteerfile)),nqtot2);
      } else {
         for (unsigned int k = 0 ; k<fScaleFac.size() ; k++) {
            double muDn = fScaleFac[k]*wrmMu1Dn[i];
            double muUp = fScaleFac[k]*wrmMu1Up[i];
            fKernMuS[i][k] = MakeInterpolationKernels(STRING_NS(Mu1_Kernel,fSteerfile),muDn,muUp);
            fKernMuS[i][k]->MakeGrids(fastNLOInterpolBase::TranslateGridType(STRING_NS(Mu1_DistanceMeasure,fSteerfile)),nqtot1);
         }
      }
   }
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
      warn["MakeInterpolationKernels"]<<"Cannot find kernel routine:" <<KernelName<<" or kernel not (yet) implemented. Exiting."<<endl;
      exit(1);
   }
   return NULL; // default return
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::TestParameterInSteering(const string& label) const {
   //! Get flag if parameter exists in steering card
   return read_steer::getexist(label, fSteerfile);
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& label, bool& val) const {
   //! Get boolean value from steering with label 'label'.
   //! Alternatively, also BOOL_NS(`label`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering labels, one should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static bool IsCMS;
   //! static bool gotval = GetParameterFromSteering("MjjCut",IsCMS);
   //! if (!gotval) cout<<"Error! Could not find boolean parameter MjjCut in steering file."<<endl;
   //!
   //! Function returns 'false' if label was not found in steering file

   bool exist = read_steer::getexist(label,fSteerfile);
   if ( exist )
      val = read_steer::getbool(label,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& label, int& val) const {
   //! Get integer value from steering with label 'label'.
   //! Alternatively, also INT_NS(`label`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering labels, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static int nJetMin;
   //! static bool gotval = GetParameterFromSteering("nJetMin",nJetMin);
   //! if (!gotval) cout<<"Error! Could not find integer parameter nJetMin in steering file."<<endl;
   //!
   //! Function returns 'false' if label was not found in steering file

   bool exist = read_steer::getexist(label,fSteerfile);
   if ( exist )
      val = read_steer::getint(label,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& label, double& val) const {
   //! Get boolean value from steering with label 'label'.
   //! Alternatively, also DOUBLE_NS(`label`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering labels, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static double MjjCut;
   //! static bool gotval = GetParameterFromSteering("MjjCut",MjjCut);
   //! if (!gotval) cout<<"Error! Could not find parameter MjjCut in steering file."<<endl;
   //!
   //! Function returns 'false' if label was not found in steering file

   bool exist = read_steer::getexist(label,fSteerfile);
   if ( exist )
      val = read_steer::getdouble(label,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& label, string& val) const {
   //! Get string value from steering with label 'label'.
   //! Alternatively, also STRING_NS(`label`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering labels, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static string text;
   //! static bool gotval = GetParameterFromSteering("MyText",text);
   //! if (!gotval) cout<<"Error! Could not find parameter MyText in steering file."<<endl;
   //!
   //! Function returns 'false' if label was not found in steering file

   bool exist = read_steer::getexist(label,fSteerfile);
   if ( exist )
      val = read_steer::getstring(label,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& label, vector<int>& val) const {
   //! Get integer vector from steering with label 'label'.
   //! Alternatively, also INT_ARR(`label`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering labels, you  should not
   //! call this function too frequently.
   //!
   //! Function returns 'false' if label was not found in steering file

   bool exist = read_steer::getexist(label,fSteerfile);
   if ( exist )
      val = read_steer::getintarray(label,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& label, vector<double>& val) const {
   //! Get vector of doubles from steering with label 'label'.
   //! Alternatively, also DOUBLE_ARR_NS(`label`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering labels, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static vector<double> FlexiCuts;
   //! static bool gotval = GetParameterFromSteering("FlexiCuts",FlexiCuts);
   //! if (!gotval) cout<<"Error! Could not find vector FlexiCuts in steering file."<<endl;
   //!
   //! Function returns 'false' if label was not found in steering file

   bool exist = read_steer::getexist(label,fSteerfile);
   if ( exist )
      val = read_steer::getdoublearray(label,fSteerfile);
   return exist;

}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& label, vector<string>& val) const {
   //! Get vector of strings from steering with label 'label'.
   //! Alternatively, also STRING_ARR_NS(`label`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering labels, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static string text;
   //! static bool gotval = GetParameterFromSteering("MyText",text);
   //! if (!gotval) cout<<"Error! Could not find parameter MyText in steering file."<<endl;
   //!
   //! Function returns 'false' if label was not found in steering file

   bool exist = read_steer::getexist(label,fSteerfile);
   if ( exist )
      val = read_steer::getstringarray(label,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& label, vector<vector<int > >& val) const {
   //! Get vector of vectors of ints from steering with label 'label'.
   //! Alternatively, also INT_TAB_NS(`label`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering labels, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static string text;
   //! static bool gotval = GetParameterFromSteering("MyText",text);
   //! if (!gotval) cout<<"Error! Could not find parameter MyText in steering file."<<endl;
   //!
   //! Function returns 'false' if label was not found in steering file

   bool exist = read_steer::getexist(label,fSteerfile);
   if ( exist )
      val = read_steer::getinttable(label,fSteerfile);
   return exist;
}


// ___________________________________________________________________________________________________
bool fastNLOCreate::GetParameterFromSteering(const string& label, vector<vector<double > >& val) const {
   //! Get vector of vector of doubles from steering with label 'label'.
   //! Alternatively, also DOUBLE_TAB_NS(`label`) could be used if read_steer.h is included
   //!
   //! Since a string (or a hash-map-access) has to be performed
   //! during access of the steering labels, you  should not
   //! call this function too frequently.
   //!
   //! Use for istance:
   //! static string text;
   //! static bool gotval = GetParameterFromSteering("MyText",text);
   //! if (!gotval) cout<<"Error! Could not find parameter MyText in steering file."<<endl;
   //!
   //! Function returns 'false' if label was not found in steering file

   bool exist = read_steer::getexist(label,fSteerfile);
   if ( exist )
      val = read_steer::getdoubletable(label,fSteerfile);
   return exist;
}
