// Daniel Britzger
// DESY, 29.07.2013
#ifndef __fastNLOCreate__
#define __fastNLOCreate__

#include <math.h>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <ctime>
#include "fastNLOReader.h" // for reference tables
#include "fastNLOTable.h"
#include "fastNLOEvent.h"
#include "read_steer.h"

#include "fastNLOInterpolBase.h"
#include "fastNLOCoeffAddBase.h"
#include "fastNLOGeneratorConstants.h"


class fastNLOCreate : public fastNLOTable {
   //!
   //! fastNLOCreate. A class for creating a fastNLO Table which contains
   //! exactly one table of coefficients.
   //!
   //! Member variables are initialized by reading in
   //! a steering file.
   //!
   //! Following information has to be obtained from the generator and is NOT obtained from steering:
   //!   - Order in alpha_s of leading-order process
   //!   - Center of mass energy
   //!   - Order of calculation (LO=0, NLO=1, NNLO=2)
   //!

public:
   fastNLOCreate(const fastNLO::GeneratorConstants& GenConsts, const fastNLO::ProcessConstants& ProcConsts,
                 const fastNLO::ScenarioConstants& ScenConsts, const fastNLO::WarmupConstants& WarmConsts);
   // KR: Deprecated, to be replaced by
   fastNLOCreate(const std::string& warmupfile, const fastNLO::GeneratorConstants& GenConsts,
                 const fastNLO::ProcessConstants& ProcConsts, const fastNLO::ScenarioConstants& ScenConsts);
   // KR: this one.
   fastNLOCreate(const fastNLO::GeneratorConstants& GenConsts, const fastNLO::ProcessConstants& ProcConsts,
                 const fastNLO::ScenarioConstants& ScenConsts, const std::string& warmupfile);
   fastNLOCreate(const fastNLO::GeneratorConstants& GenConsts, const fastNLO::ProcessConstants& ProcConsts,
                 const fastNLO::ScenarioConstants& ScenConsts, const std::string& warmupfile,
                 const std::string& steerfile);
   // KR: Deprecated, to be replaced by
   fastNLOCreate(const std::string& steerfile, const fastNLO::GeneratorConstants& GenConsts,
                 const fastNLO::ProcessConstants& ProcConsts);
   // KR: this one.
   fastNLOCreate(const fastNLO::GeneratorConstants& GenConsts, const fastNLO::ProcessConstants& ProcConsts,
                 const std::string& steerfile);
   // KR: Don't see any use for shouldReadSteeringFile ==> remove it
   //   fastNLOCreate(const std::string& steerfile, std::string warmupfile = "", bool shouldReadSteeringFile = true);
   fastNLOCreate(const std::string& steerfile, std::string steeringNameSpace = "");
   ~fastNLOCreate();

   fnloEvent fEvent;                                                                            //!< Structure, which holds all relevant variables related to event observables
   fnloScenario fScenario;                                                                      //!< Structure, which holds perturbative (wilson) coefficients/weights and x-values

   void SetOrderOfAlphasOfCalculation(unsigned int ord);                                        //!< set absolute order of alpha_s
   void SetScenario(const fnloScenario scen) {fScenario = scen;}                                //!< set the member fScenario, which will be used when calling Fill()
   void SetEvent(const fnloEvent ev) {fEvent = ev;}                                             //!< set the member fEvent, which will be used when calling Fill()
   inline void SetNumberOfEvents(double n) {
      GetTheCoeffTable()->Nevt = n;
      GetTheCoeffTable()->fWgt.WgtNevt = n;
      fStats._nEv=n;};             //!< set number of events. This is only mandatory, before calling WriteTable().
   void SetLoOrder(int LOOrd);                                                                  //!< set order of alpha_s for leading order process.

   fastNLOReader* SetIsReferenceTable(fastNLOReader* fnloread = NULL);                          //!< set this table/contribution to become a reference contribution

   void ReadSteeringFile(std::string steerfile, std::string steeringNameSpace = "");            //!< only read steering file, do not set anything
   void TransformPartonCombinations();                                                          //!< somehow transform parton combinations, why?
   void SetTableConstsDefaults();                                                               //!< set defaults for all table constants
   void SetGenConstsDefaults();                                                                 //!< set defaults for generator constants
   void SetProcConstsDefaults();                                                                //!< set defaults for process constants
   void SetScenConstsDefaults();                                                                //!< set defaults for scenario constants
   void SetWarmupConstsDefaults();                                                              //!< set defaults for warmup constants
   void SetGenConstsFromSteering();                                                             //!< set generator constants from steering
   void SetProcConstsFromSteering();                                                            //!< set process constants from steering
   void SetScenConstsFromSteering();                                                            //!< set scenario constants from steering
   void SetWarmupConstsFromSteering();                                                          //!< set warmup constants from steering
   void PrintTableConsts();                                                                     //!< print all table constants
   void PrintGenConsts();                                                                       //!< print generator constants
   void PrintProcConsts();                                                                      //!< print process constants
   void PrintScenConsts();                                                                      //!< print scenario constants
   void PrintWarmupConsts();                                                                    //!< print warmup constants
   bool CheckTableConsts();                                                                     //!< check all table constants
   bool CheckGenConsts();                                                                       //!< check generator constants
   bool CheckProcConsts();                                                                      //!< check process constants
   bool CheckScenConsts();                                                                      //!< check scenario constants
   bool CheckWarmupConsts();                                                                    //!< check warmup constants

   // SetBinGrid()
   // todo: SetBinGrid. However, if BinGrid is set, then this is necessarily a warmup run -> one also has to store the bin grid in warmup table (todo).
   //       furthermore all vectors have to be 'resized'
   //void SetBinGrid(std::vector < std::vector <std::pair<double,double> > > BinGrid, std::vector <int> IDiffBin, std::vector <std::string> DimLabel, std::vector <double> BinSize = std::vector <double>() );

   void Fill(int scalevar=0, const double wgtfac=1.0);                                          //!< fill event quantities in fastNLO table. Call it for every subprocess.
   void FillOneSubprocess(const fnloEvent& event, const fnloScenario& scen, int scalevar=0, const double wgtfac=1.0); //!< same function as 'Fill()', but uses content of member fScenario and fEvent
   void FillAllSubprocesses(const std::vector<fnloEvent>& events, const fnloScenario& scen, int scalevar=0, const double wgtfac=1.0); //!< Fill a selection (std::vector) of events/processes/channels, which all have the identic scenario
   void FillAllSubprocesses(const std::vector<std::vector<fnloEvent> >& events, const fnloScenario& scen, const double wgtfac=1.0);        //!< Fill a list of subprocesses for various scale-variations into a fixed-scale table
   int GetNSubprocesses() const { return GetTheCoeffTable()->GetNSubproc();}                    //!< The number of subprocesses (channels)
   const std::vector<double>& GetScaleVariations() const { return fScaleFac; }                       //!< Get list of scale variations
   int GetScaleVarMax() const { return fScaleFac.size(); }                                //!< Get no. of scale variations

   void WriteTable(std::string filename);                                                            //!< Write fastNLO table to file filename
   void WriteTable();                                                                           //!< Write fastNLO table to disk.
   inline void WriteTable(double nEvents) {SetNumberOfEvents(nEvents); WriteTable();};                 //!< Write fastNLO table to disk and set number of events.
   void WriteWarmupTable();                                                                     //!< Write the warmup table to disk.
   void MultiplyCoefficientsByBinSize();                                                        //!< Multiply all coefficients by bin size
   void DivideCoefficientsByBinSize();                                                          //!< Divide all coefficients by bin size
   void MultiplyCoefficientsByConstant(double c);                                               //!< Multiply all coefficients with a constant factor c
   void NormalizeCoefficients(double wgt=1);                                                    //!< Set number of events to wgt and adjust coefficients accordingly
   void NormalizeCoefficients(const std::vector<std::vector<double> >& wgtProcBin);             //!< Set number of events to wgt (with wgt differently for each subproc. and bin) resulting Nevt will become 1

   void PrintWarmupValues();                                                                    //!< Print the warmup values to the screen
   std::string GetWarmupTableFilename();                                                             //!< Get the filename, which is used for storage of the warmup-table.
   void SetWarmupTableFilename(std::string);                                                         //!< Set the filename, which is used for storage of the warmup-table (otherwise a default is used)
   inline bool GetIsWarmup() const { return fIsWarmup; };                                              //!< Get flag for warmup table
   inline bool GetIsFlexibleScale() const { return fIsFlexibleScale; }                          //!< Get if flexible-scale table is requested
   inline bool GetIsInclusiveJets() const { return fIsInclusiveJets; }                          //!< Get if inclusive jet table is requested
   void SetWarmupXSafetyMargin(int margin = 4) { fWarmupXMargin = margin;}                          //!< Set margin for x-value: First digit in '%e' notation (e.g. margin=4: x=6.6e-3 -> 6.2e-3)
   void SetWarmupNDigitMu1(int rnd = 1) { fWarmupNDigitMu1 = rnd;}                              //!< Round warmup values for scale 1
   void SetWarmupNDigitMu2(int rnd = 2) { fWarmupNDigitMu2 = rnd;}                              //!< Round warmup values for scale 2

   inline fastNLOCoeffAddBase* GetTheCoeffTable() const {
      return (fastNLOCoeffAddBase*)GetCoeffTable(0);
   }                                            //!< Getter for the one (and only) coefficient table

   bool TestParameterInSteering(const std::string& label) const;                                           //!< Test on existence of user-defined parameter name in steering card.
   bool GetParameterFromSteering(const std::string& label, bool& val) const;                                      //!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(const std::string& label, int& val) const;                                       //!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(const std::string& label, double& val) const;                                    //!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(const std::string& label, std::string& val) const;                                    //!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(const std::string& label, std::vector<int>& val) const;                               //!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(const std::string& label, std::vector<double>& val) const;                            //!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(const std::string& label, std::vector<std::string>& val) const;                            //!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(const std::string& label, std::vector<std::vector<int > >& val) const;                     //!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(const std::string& label, std::vector<std::vector<double > >& val) const;                  //!< Get user-defined parameter from steering file.

   void AdjustWarmupValues();                                                                   //!< Round warmup values to more likely values.
   void PrintAllSteeringValues() const { PRINTALL();};                                          //!< Print all steering values obtained from steering files (of all fastNLOCreate instances);

   void Clear() { GetTheCoeffTable()->Clear();};                                                //!< Clear coefficient table
   void PrintStats() const { fStats.PrintStats();}                                              //!< Print statistics
   void SetGlobalVerbosity(std::string sverb);                                                  //!< Set GlobalVerbosity using std::string variable

   // KR Allow nnlo-bridge_0.0.40 to be linked
   void SetCacheSize(int MaxCache, int CacheComp=2, int CacheType=2);                           //!< Set maximum number of events in cache. Set to 0 to deactivate caching, CacheComp: maximum numbers of entries to be compared with new element; fCacheType: 0: deactivate, 1: 1D cache, 2: cache per ObsBin and Proc
   int GetCacheSize() const { return fCacheMax; }                                               //!< Get max cache size
   int GetCacheComp() const { return fCacheComp; }                                              //!< Get cache comp value
   int GetCacheType() const { return fCacheType; }                                              //!< Get cache type

protected:
   fastNLOCreate();                                                                             //!< don't use the default constructor. fastNLOCreate is only reasonable with input steering.
   void Instantiate();
   int CreateCoeffTable();                                                                      //!< Create the one (and only) coefficient table

   inline void ApplyPDFWeight(std::vector<std::pair<int,double> >& nodes, const double x, const std::vector<double>* grid) const;
   inline double CalcPDFReweight(double x) const;
   //!< Fill contribution into table
   void FillContribution(int scalevar = 0, const double wgtfac=1.0);
   //!< Fill flexible scale contribution in pp/ppbar
   void FillContributionFlexHHC(fastNLOCoeffAddFlex* c, int ObsBin, const double wgtfac=1.0);
   //!< fill flexible scale contribution in DIS
   void FillContributionFlexDIS(fastNLOCoeffAddFlex* c, int ObsBin, const double wgtfac=1.0);
   //!< fill fixed scale table in pp/ppbar
   void FillContributionFixHHC(fastNLOCoeffAddFix* c, int ObsBin, int scalevar, const double wgtfac=1.0);
   //!< fill fixed scale contribution in DIS
   void FillContributionFixDIS(fastNLOCoeffAddFix* c, int ObsBin, int scalevar, const double wgtfac=1.0);
   //!< fill contribution if this is a reference table
   void FillRefContribution(int scalevar = 0, const double wgtfac=1.0);
   // KR: Don't see any use for shouldReadSteeringFile ==> remove it
   //   void ReadSteering(std::string steerfile, std::string steeringNameSpace = "", bool shouldReadSteeringFile = true);  //!< read steering file
   void ReadSteering(std::string steerfile, std::string steeringNameSpace = "");  //!< read steering file

   void ReadBinning();
   void ReadBinningFromScenarioConsts();
   void ReadBinSize();
   ///
   void SetBinning1D(std::vector<double> bgrid, std::string label, unsigned int idiff);
   void SetBinning1D(std::vector<double> bgrid, std::string label, unsigned int idiff, double norm);
   void SetBinning1D(std::vector<double> bgrid, std::string label, unsigned int idiff, std::vector<double> vnorm);
   void SetBinning1D(std::vector<double> blow, std::vector<double> bupp, std::string label, unsigned int idiff);
   void SetBinning1D(std::vector<double> blow, std::vector<double> bupp, std::string label, unsigned int idiff, double norm);
   void SetBinning1D(std::vector<double> blow, std::vector<double> bupp, std::string label, unsigned int idiff, std::vector<double> vnorm);
   void SetBinningND(std::vector<double> bgrid, unsigned int ndim, std::vector<int> idiff);
   void SetBinningND(std::vector<std::vector<double> > bgrid, unsigned int ndim, std::vector<int> idiff);
   ///
   void ReadCoefficientSpecificVariables();
   void ReadScaleFactors();
   void InitVariablesInCoefficientTable();
   void InitCoeffTable();
   void InitInterpolationKernels();
   fastNLOInterpolBase* MakeInterpolationKernels(std::string KernelName, double xdn, double xup);
   void InitGrids();
   void GetWarmupValues();
   bool CheckWarmupConsistency();                                                               //!< Check consistency of warmup bin-grid and variables with steering values.
   void UseBinGridFromWarmup();                                                                 //!< Use bin grid as given in the warmup table
   int CheckWarmupValuesIdenticalWithBinGrid(std::vector<std::pair<double,double> >& wrmmu);              //!< Check if warmup values are possibly identical with bin grid
   void RoundValues(std::vector<std::pair<double,double> >& wrmmu, int ibin=-1 , int nth=3);                   //!< Round values to closes value by at most 1%
   int GetNthRelevantDigit(double val, int n);
   std::vector<std::vector<std::pair<int,int> > > ReadPartonCombinations(int ord, const std::vector<std::vector<int> >& PartonCombinations);                             //!< Read PDFCoeff from steering

   int GetBin();                                                                                //!< get bin number from 'scenario' observables
   inline int GetXIndex(const int& Obsbin, const int& x1bin, const int& x2bin) const;           //!< get x-index in case of two hadrons.
   int GetNxmax(const std::vector<double>* xGrid1, const std::vector<double>* xGrid2);                    //!< get maximum x-index
   bool fIsWarmup;                                                                              //!< is it a warmup run?
   int fWarmupXMargin;                                                                          //!< margin for x-value: First digit in '%e' notation (e.g. margin=4: x=6.6e-3 -> 6.2e-3)
   int fWarmupNDigitMu1;                                                                          //!< Digits of warmup values for scale 1
   int fWarmupNDigitMu2;                                                                          //!< Digits of warmup values for scale 2
   int  fIOrd;                                                                                  //!< order of alpha_s of run
   bool fIsFlexibleScale;                                                                       //!< is it a flexible scale table?
   bool fIsInclusiveJets;                                                                       //!< is it an inclusive jets table?
   bool fApplyPDFReweight;                                                                      //!< shall the PDF reweight be applied.
   std::string fWarmupFilename;                                                                 //!< File name of the warmup table
   // TODO: fSteerfile as internally set variable is used as default steeringNameSpace at many places. This should be avoided!
   std::string fSteerfile;                                                                      //!< filename of steering file (and default steeringNameSpace!)
   int fObsBin;                                                                                 //!< ObsBin from 'last' 'Fill()'-call
   fnloScenario fLastScen;                                                                      //!< keep information of scenario from last 'Fill()'-call

   fastNLOReader* fReader;                                                                      //!< PDF and alpha_s interface for reference tables

   fastNLO::GeneratorConstants fGenConsts;                                                      //!< Generator specific constants
   fastNLO::ProcessConstants fProcConsts;                                                       //!< Process specific constants
   fastNLO::ScenarioConstants fScenConsts;                                                      //!< Scenario specific constants
   fastNLO::WarmupConstants fWarmupConsts;                                                      //!< Warmup specific constants

   bool CheckWeightIsFinite();                                                                  //!< Check if weight is reasonable.
   inline void HalfMatrixCheck(double x1, double x2, int& xmin, int& xmax, int& subproc) const;                       //!< check x-values in case of half-matrix notation (pp,ppbar), and exchange if necessary.
   std::vector<int> fSymProc;                                                                        //!< necessary for half-matrix notation
   std::vector<double> fScaleFac;                                                                    //!< Scale factors. Needed for fixed-scale tables

   // interpolation kernels
   std::vector<fastNLOInterpolBase*> fKernX1;                                                        //!< Interpolation kernel for x-interpolation
   std::vector<fastNLOInterpolBase*> fKernX2;                                                        //!< Interpolation kernel for x-interpolation
   std::vector<fastNLOInterpolBase*> fKernMu1;                                                       //!< Interpolation kernel for mu1-interpolation
   std::vector<fastNLOInterpolBase*> fKernMu2;                                                       //!< Interpolation kernel for mu2-interpolation
   std::vector<std::vector<fastNLOInterpolBase*> > fKernMuS;                                              //!< Interpolation kernels for each scale var for fixed-scale tables

   // arrays for warmup
   void UpdateWarmupArrays();
   void InitWarmupArrays();
   void OutWarmup(std::ostream& = std::cout);
   std::vector<double> GetColumnFromTable(const std::vector<std::vector<double> >& table, int iCol);      //!< Get a column from a table
   std::string GetWarmupHeader(int iScale, std::string minmax);
   std::vector<std::pair<double,double> > fWMu1;                                                          //!< array of warmup-up values
   std::vector<std::pair<double,double> > fWMu2;                                                          //!< array of warmup-values
   std::vector<std::pair<double,double> > fWx;                                                            //!< array of warmup-values
   std::vector<std::pair<double,double> > fWMu1Rnd;                                                       //!< copy of warm-up array for rounding
   std::vector<std::pair<double,double> > fWMu2Rnd;                                                       //!< copy of warm-up array for rounding
   std::vector<std::pair<double,double> > fWxRnd;                                                         //!< copy of warm-up array for rounding

   int fCacheMax = 30;                                                                                 //!< maximum number of entries in weight cahce
   int fCacheComp = 5;                                                                                        //!< maximum number of entries in weight cahce
   int fCacheType = 2;                                                                                        //!< maximum number of entries in weight cahce
   std::vector<std::pair<fnloScenario,fnloEvent> > fWeightCache;                                          //!< cache for fill-weights
   std::vector<std::vector<std::vector<std::pair<fnloScenario,fnloEvent> > > > fWeightCacheBinProc;       //!< alternative cache for fill-weights
   void FillWeightCache(int scalevar);                                                                    //!< Fill weight into cache, merge weights for identical phase space points
   void FlushCache();                                                                                     //!< Fill weights from cache into table


   struct fnloStats {
      //! structre to keep track of statisics. Just for fun and information.
      time_t _time;
      fnloStats() : _nProc(0), _nEvPS(0), _nEv(0)  { _time = time(0);}
      long long int _nProc, _nEvPS;
      double _nEv;
      void PrintStats() const {
         time_t hour, min, time = std::time(0) - _time;
         hour = time/3600L;
         time -= hour*3600L;
         min  = time/60L;
         time -= min*60L;
         std::cout<<std::endl;;
         std::cout<<" ------------------- fastNLOstats ------------------"<<std::endl;;
         std::cout<<"   Time elapsed:                 "
                  << (hour < 10 ? "0" : "")   << hour
                  << (min < 10 ? ":0" : ":")  << min
                  << (time < 10 ? ":0" : ":") << time << std::endl;;
         if (_nEv!=0)   std::cout << "   Total event weight (NEvt):     " << _nEv   << std::endl;;
         if (_nEvPS!=0) std::cout << "   Contributions in phase space:  " << _nEvPS << std::endl;;
         if (_nProc!=0) std::cout << "   Number of calls:                " << _nProc << std::endl;;
         std::cout << " -------------------------------------------------" << std::endl;;
         std::cout<<std::endl;;
         std::cout.flush();
      }
   } fStats;

};
#endif
