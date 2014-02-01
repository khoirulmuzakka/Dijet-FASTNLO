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
#include "fastNLOTable.h"
#include "fastNLOEvent.h"
#include "read_steer.h"

#include "fastNLOInterpolBase.h"
#include "fastNLOCoeffAddBase.h"
#include "fastNLOGeneratorConstants.h"


using namespace std;


class fastNLOCreate : public fastNLOTable {
   //!
   //! fastNLOCreate. A class for creating a fastNLO Table which contains
   //! exaclty one table of coefficients.
   //!
   //! Member variables are initialized by reading in 
   //! a steering file. 
   //!
   //! Following information has to be obtained from the generator and is NOT obtained from steering:
   //!   - Order in alpha_s of leading-order process
   //!   - Center of mass energy
   //!   - Order of calculation (LO=0, NLO=1,NNLO=2)
   //!

public:
   fastNLOCreate(string steerfile);
   fastNLOCreate(string steerfile, fastNLO::GeneratorConstants GenConsts, fastNLO::ProcessConstants ProcConsts);
   ~fastNLOCreate();
   static int nInst;										//!< limit the number of instance to one (until more flexible access of parser values is implemented)

   fnloEvent fEvent;										//!< Structure, which holds all relevant variables related to event observables
   fnloScenario fScenario;									//!< Structure, which holds perturbative (wilson) coefficients/weights and x-values

   void SetOrderOfAlphasOfCalculation(unsigned int ord);					//!< set absolute order of alpha_s 
   void SetScenario(const fnloScenario scen) {fScenario = scen;}				//!< set the member fScenario, which will be used when calling Fill()
   void SetEvent(const fnloEvent ev) {fEvent = ev;}						//!< set the member fEvent, which will be used when calling Fill()
   void SetNumberOfEvents(long long int n) {GetTheCoeffTable()->Nevt = n; fStats._nEv=n;};	//!< set number of events. This is only mandatory, before calling WriteTable().	
   void SetLoOrder(int LOOrd);									//!< set order of alpha_s for leading order process.

   // SetBinGrid()
   // todo: SetBinGrid. However, if BinGrid is set, then this is necessarliy a warmup run -> one also has to store the bin grid in warmup table (todo).
   //       furthermore all vectors have to be 'resized'
   //void SetBinGrid(vector < vector <pair<double,double> > > BinGrid, vector <int> IDiffBin, vector <string> DimLabel, vector <double> BinSize = vector <double>() );
 
   void Fill(int scalevar=0);									//!< fill event quantities in fastNLO table. Call it for every subprocess.
   void FillOneSubprocess(const fnloEvent& event, const fnloScenario& scen, int scalevar=0);	//!< same function as 'Fill()', but uses content of member fScenario and fEvent
   void FillAllSubprocesses(const vector<fnloEvent>& events, const fnloScenario& scen, int scalevar=0);	//!< Fill a selection (vector) of events/processes/channels, which all have the identic scenario
   void FillAllSubprocesses(const vector<vector<fnloEvent> >& events, const fnloScenario& scen);	//!< Fill a list of subprocesses for various scale-variations into a fixed-scale table
   int GetNSubprocesses() const { return GetTheCoeffTable()->GetNSubproc();}			//!< The number of subprocesses (channels)
   const vector<double>& GetScaleVariations() const { return fScaleFac; }			//!< Get list of scale variations

   void WriteTable(string filename);								//!< Write fastNLO table to file <filename>
   void WriteTable();										//!< Write fastNLO table to disk.
   void WriteWarmupTable();									//!< Write the warmup table to disk.
   void MultiplyCoefficientsByBinWidth();							//!< Multiply all coefficients by binwidth 
   void DivideCoefficientsByBinWidth();								//!< Divide all coefficients by binwidth 
   void MultiplyCoefficientsByConstant(double c);						//!< Multiply all coefficients with a constant factor c

   void PrintWarmupValues();									//!< Print the warmup values to the screen
   string GetWarmupTableFilename();								//!< Get the filename, which is used for storage of the warmup-table.
   
   fastNLOCoeffAddBase* GetTheCoeffTable() const {						
      return (fastNLOCoeffAddBase*)GetCoeffTable(0);}						//!< Getter for the one (and only) coefficient table

   bool GetParameterFromSteering(string label, bool& val);					//!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(string label, int& val);					//!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(string label, double& val);					//!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(string label, string& val);					//!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(string label, vector<int>& val);				//!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(string label, vector<double>& val);				//!< Get user-defined parameter from steering file.
   bool GetParameterFromSteering(string label, vector<string>& val);				//!< Get user-defined parameter from steering file.
   

protected:
   fastNLOCreate();										//!< don't use the default constructor. fastNLOCreate is only reasonable with input steering.
   void Instantiate();
   int CreateCoeffTable();									//!< Create the one (and only) coefficient table

   inline void ApplyPDFWeight(vector<pair<int,double> >& nodes, const double x, const vector<double>* grid ) const;
   inline double CalcPDFReweight(double x) const;
   void FillContribution(int scalevar = 0);							//!< fill contribution into table
   void FillContributionFlexHHC(fastNLOCoeffAddFlex* c, int ObsBin);				//!< fill flexible scale contribution in pp/ppbar
   void FillContributionFlexDIS(fastNLOCoeffAddFlex* c, int ObsBin);				//!< fill flexible scale contribution in DIS
   void FillContributionFixHHC(fastNLOCoeffAddFix* c, int ObsBin, int scalevar);		//!< fill fixed scale table in pp/ppbar
   void ReadSteering(string steerfile);								//!< read steering file
   void ReadGenAndProcConstsFromSteering();
   void ReadBinning();
   void ReadCoefficientSpecificVariables();
   void ReadScaleFactors();
   void InitVariablesInCoefficientTable();
   void InitCoeffTable();
   void InitInterpolationKernels();
   fastNLOInterpolBase* MakeInterpolationKernels(string KernelName, double xdn, double xup);
   void InitGrids();
   void GetWarmupValues();
   bool CheckWarmupConsistency();								//!< Check consistency of warmup bin-grid and variables with steering values.
   void UseBinGridFromWarmup();									//!< Use bin grid as given in the warmup table
   int GetBin();										//!< get bin number from 'scenario' observables
   inline int GetXIndex(const int& Obsbin, const int& x1bin, const int& x2bin) const;		//!< get x-index in case of two hadrons.
   int GetNxmax(const vector<double>* xGrid1, const vector<double>* xGrid2);			//!< get maximum x-index	
   bool fIsWarmup;										//!< is it a warmup run?
   int  fIOrd;											//!< order of alpha_s of run
   bool fIsFlexibleScale;									//!< is it a flexible scale table?
   bool fApplyPDFReweight;									//!< shall the PDF reweight be applied.
   string fSteerfile;										//!< filename of steering file.	
   int fObsBin;											//!< ObsBin from 'last' 'Fill()'-call
   fnloScenario fLastScen;									//!< keep information of scenario from last 'Fill()'-call
   
   fastNLO::GeneratorConstants fGenConsts;							//!< Generator specific constants
   fastNLO::ProcessConstants fProcConsts;							//!< Process specific constants
 
   bool CheckWeightIsNan();									//!< Check if weight is reasonable.
   inline void HalfMatrixCheck(int& xmin, int& xmax, int& subproc) const;			//!< check x-values in case of half-matrix notation (pp,ppbar), and exchange if necessary.
   vector<int> fSymProc;									//!< necessary for half-matrix notation
   vector<double> fScaleFac;									//!< Scale factors. Needed for fixed-scale tables

   // interpolation kernels
   vector<fastNLOInterpolBase*> fKernX;								//!< Interpolation kernel for x-interpolation
   vector<fastNLOInterpolBase*> fKernMu1;							//!< Interpolation kernel for mu1-interpolation
   vector<fastNLOInterpolBase*> fKernMu2;							//!< Interpolation kernel for mu2-interpolation
   vector<vector<fastNLOInterpolBase*> > fKernMuS;						//!< Interpolation kernels for each scale var for fixed-scale tables

   // arrays for warmup
   void UpdateWarmupArrays();
   void InitWarmupArrays();
   void OutWarmup(ostream& = std::cout);
   string GetWarmupHeader(int iScale, string minmax );
   vector<pair<double,double> > fWMu1;								//!< array of warmup-up values
   vector<pair<double,double> > fWMu2;								//!< array of warmup-values
   vector<pair<double,double> > fWx;								//!< array of warmup-values

   struct fnloStats {
      //! structre to keep track of statisics. Just for fun and information.
      std::time_t _time;
      long long int _nProc, _nEv, _nEvPS;
      fnloStats() : _nProc(0), _nEv(0) , _nEvPS(0){ _time = std::time(0);}
      void PrintStats() {
	 time_t hour, min, time = std::time(0) - _time;
	 hour = time/3600L;
	 time -= hour*3600L;
	 min  = time/60L;
	 time -= min*60L;
	 std::cout<<std::endl;
	 std::cout<<"------------ fastNLOstats ------------"<<std::endl;
	 std::cout<<"  Time elapsed:       "
		  <<(hour < 10 ? "0" : "")<<hour
		  <<(min < 10 ? ":0" : ":")<<min
		  <<(time < 10 ? ":0" : ":")<<time<<std::endl;
	 if ( _nEv!=0)   std::cout<<"  Events filled:      "<<_nEv<<std::endl;
	 if ( _nEvPS!=0) std::cout<<"  Phase space events: "<<_nEvPS<<std::endl;
	 if ( _nProc!=0) std::cout<<"  Number of calls:    "<<_nProc<<std::endl;
	 std::cout<<"--------------------------------------"<<std::endl;
	 std::cout.flush();
      }
   } fStats;

};
#endif
