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
#include "read_steer.h"

#include "fastNLOInterpolBase.h"
#include "fastNLOInterpolCatmulRom.h"
#include "fastNLOCoeffAddBase.h"

using namespace std;

struct fnloScenario {
   // useful class to keep all scenario specific quantities.
   // e.g. observables and scales
   friend class fastNLOCreate;
   fnloScenario() : _iOB(-1){;}
   void SetObservableDimI(double o, int iDim) {_o[iDim]=o;}					// Set observable of dimension iDim (e.g. in case of multidimensional measurements)
   void SetObsBin(int iBin) {_iOB = iBin; }							// [optional] Set ObsBin (e.g. if binning is performed by generator, no other observables are then needed.)
   // flexible scale table:
   void SetObsScale1(double mu) {_m1=mu;}							// For flexible-scale tables. Set scale 1 (should be in 'GeV').
   void SetObsScale2(double mu) {_m2=mu;}							// For flexible-scale tables. Set scale 2
   // if not a flexible-scale table
   void SetScale(double mu) {_m1=mu;}								// the ren. and fact. scale (not mu^2)
private:
   map<int,double> _o;
   double _m1, _m2;
   int _iOB;
};


struct fnloEvent {
   // useful class to keep all process related variables.
   // e.g x-values, weights, process identifiers, etc...
   friend class fastNLOCreate;
   void Reset(){
      _x1=0,_x2=0;
      _w=0,_wf=0,_wr=0,_wrr=0,_wff=0,_wrf=0;
      _p=-1;
      _n=-1;
   }
   void ResetButX(){
      _w=0,_wf=0,_wr=0,_wrr=0,_wff=0,_wrf=0;
      _p=-1;
      _n=-1;
   }
   fnloEvent(){Reset();}
   // event specific quantites, which are required for every 'Fill()' step.
   void SetX(double x) {_x1=x;}									// set x-value of first hadron (if e.g. DIS)
   void SetX1(double x) {_x1=x;}								// setx-value of first hadron
   void SetX2(double x) {_x2=x;}								// set x-value of second hadron
   void SetProcessId(int n){_p=n;}								// set identifier of specific subprocess (0<n<NSubproc), according to the corresponding PDF linear combination
   void SetEventCounter(long long int n){_n=n;}							// Set event counter
   // if not a flexible-scale table
   void SetWeight(double w) {_w=w;}								// weights must be mutliplied with dummypdf (1/x)
   // flexible scale table:
   void SetWeight_MuIndependent(double w) {_w=w;}						// weights must be mutliplied with dummypdf (1/x)
   void SetWeight_log_mur(double w) {_wr=w;}							// set weight w, which will contribute with log_e(mur^2)*w
   void SetWeight_log_muf(double w) {_wf=w;}							// set weight w, which will contribute with log_e(muf^2)*w
   void SetWeight_log_murr(double w) {_wrr=w;}							// set weight w, which will contribute with log^2_e(mur^2)*w
   void SetWeight_log_muff(double w) {_wff=w;}							// set weight w, which will contribute with log^2_e(muf^2)*w
   void SetWeight_log_murf(double w) {_wrf=w;}							// set weight w, which will contribute with log_e(mur^2)*log_e(muf^2)*w
   void AddWeight_MuIndependent(double w) {_w+=w;}						// weights must be mutliplied with dummypdf (1/x)
   void AddWeight_log_mur(double w) {_wr+=w;}							// set weight w, which will contribute with log_e(mur^2)*w
   void AddWeight_log_muf(double w) {_wf+=w;}							// set weight w, which will contribute with log_e(muf^2)*w
   void AddWeight_log_murr(double w) {_wrr+=w;}							// set weight w, which will contribute with log^2_e(mur^2)*w
   void AddWeight_log_muff(double w) {_wff+=w;}							// set weight w, which will contribute with log^2_e(muf^2)*w
   void AddWeight_log_murf(double w) {_wrf+=w;}							// set weight w, which will contribute with log_e(mur^2)*log_e(muf^2)*w
private:
   double _x1, _x2;										// an event has always identical x1 and x2;
   double _w, _wf, _wr, _wrr, _wff, _wrf;							// weights	
   int _p;											// processId/channel. Must be consistent with PDF linear combination
   long long int _n;										// event count
}; 



class fastNLOCreate : public fastNLOTable {
   //
   // fastNLOCreate. A class for creating a fastNLO Table which contains
   // exaclty one table of coefficients.
   //
   // Member variables are initialized by reading in 
   // a steering file. 
   //
   // Following information has to be obtained from the generator and is NOT obtained from steering:
   //   - Order in alpha_s of leading-order process
   //   - Center of mass energy
   //   - Order of calculation (LO=0, NLO=1,NNLO=2)
   //

public:
   fastNLOCreate(string steerfile);
   ~fastNLOCreate();
   static int nInst;										// limit the number of instance to one (until more flexible access of parser values is implemented)

   fnloEvent fEvent;										// Structure, which holds all relevant variables related to event observables
   fnloScenario fScenario;									// Structure, which holds perturbative (wilson) coefficients/weights and x-values

   void SetOrderOfAlphasOfCalculation(unsigned int ord);					// set absolute order of alpha_s 
   void SetScenario(const fnloScenario scen) {fScenario = scen;}				// set the member fScenario, which will be used when calling Fill()
   void SetEvent(const fnloEvent ev) {fEvent = ev;}						// set the member fEvent, which will be used when calling Fill()
   void SetNumberOfEvents(long long int n) {GetTheCoeffTable()->Nevt = n; fStats._nEv=n;};	// set number of events. This is only mandatory, before calling WriteTable().	
   void SetLoOrder(int LOOrd);									// set order of alpha_s for leading order process.

   // SetBinGrid()
   // todo: SetBinGrid. However, if BinGrid is set, then this is necessarliy a warmup run -> one also has to store the bin grid in warmup table (todo).
   //       furthermore all vectors have to be 'resized'
   //void SetBinGrid(vector < vector <pair<double,double> > > BinGrid, vector <int> IDiffBin, vector <string> DimLabel, vector <double> BinSize = vector <double>() );
 
   void Fill(int scalevar=0);									// fill event quantities in fastNLO table. Call it for every subprocess.
   void FillOneSubprocess(const fnloEvent& event, const fnloScenario& scen, int scalevar=0);	// same function as 'Fill()', but uses content of member fScenario and fEvent
   void FillAllSubprocesses(const vector<fnloEvent>& events, const fnloScenario& scen, int scalevar=0);	// Fill a selection (vector) of events/processes/channels, which all have the identic scenario
   int GetNSubprocesses() const { return GetTheCoeffTable()->GetNSubproc();}			// The number of subprocesses (channels)

   int WriteTable(string filename);								// Write fastNLO table to file <filename>
   int WriteTable();										// Write fastNLO table to disk.
   void WriteWarmupTable();									// Write the warmup table to disk.
   void MultiplyCoefficientsByBinWidth();							// Multiply all coefficients by binwidth 
   void DivideCoefficientsByBinWidth();								// Divide all coefficients by binwidth 
   void MultiplyCoefficientsByConstant(double c);						// Multiply all coefficients with a constant factor c

   void PrintWarmupValues();									// Print the warmup values to the screen
   string GetWarmupTableFilename();								// Get the filename, which is used for storage of the warmup-table.
   
   fastNLOCoeffAddBase* GetTheCoeffTable() const {						// Getter for the one (and only) coefficient table
      return (fastNLOCoeffAddBase*)GetCoeffTable(0);}

protected:
   fastNLOCreate();										// don't use the default constructor. fastNLOCreate is only reasonable with input steering.
   int CreateCoeffTable();									// Create the one (and only) coefficient table

   void ApplyPDFWeight(vector<pair<int,double> >& nodes, const double x, const vector<double>* grid );
   double CalcPDFReweight(double x);
   void FillContribution(int scalevar = 0);							// fill contribution into table
   void FillContributionFlexHHC(fastNLOCoeffAddFlex* c, int ObsBin);				// fill flexible scale contribution in pp/ppbar
   void FillContributionFixHHC(fastNLOCoeffAddFix* c, int ObsBin, int scalevar);		// fill fixed scale table in pp/ppbar
   void ReadSteering(string steerfile);								// read steering file
   void ReadBinning();
   void ReadCoefficientSpecificVariables();
   void InitVariablesInCoefficientTable();
   void InitCoeffTable();
   void InitInterpolationKernels();
   fastNLOInterpolBase* MakeInterpolationKernels(string KernelName, double xdn, double xup);
   void InitGrids();
   void GetWarmupValues();
   bool CheckWarmupConsistency();								// Check consistency of warmup bin-grid and variables with steering values.
   void UseBinGridFromWarmup();									// Use bin grid as given in the warmup table
   int GetBin();										// get bin number from 'scenario' observables
   int GetXIndex(int Obsbin,int x1bin,int x2bin);						// get x-index in case of two hadrons.
   int GetNxmax(const vector<double>* xGrid1, const vector<double>* xGrid2);			// get maximum x-index	
   bool fIsWarmup;										// is it a warmup run?
   int  fIOrd;											// order of alpha_s of run
   bool fIsFlexibleScale;									// is it a flexible scale table?
   bool fApplyPDFReweight;									// shall the PDF reweight be applied.
   string fSteerfile;										// filename of steering file.	
   int fObsBin;											// ObsBin from 'last' 'Fill()'-call
   fnloScenario fLastScen;									// keep information of scenario from last 'Fill()'-call
   
   bool CheckWeightIsNan();									// Check if weight is reasonable.
   void HalfMatrixCheck(int& xmin, int& xmax, int& subproc);					// check x-values in case of half-matrix notation (pp,ppbar), and exchange if necessary.
   vector<int> fSymProc;									// necessary for half-matrix notation

   // interpolation kernels
   vector<fastNLOInterpolBase*> fKernX;								// Interpolation kernel for x-interpolation
   vector<fastNLOInterpolBase*> fKernMu1;							// Interpolation kernel for mu1-interpolation
   vector<fastNLOInterpolBase*> fKernMu2;							// Interpolation kernel for mu2-interpolation

   // arrays for warmup
   void UpdateWarmupArrays();
   void InitWarmupArrays();
   void OutWarmup(string file);
   string GetWarmupHeader(int iScale, string minmax );
   vector<pair<double,double> > fWMu1;								// array of warmup-up values
   vector<pair<double,double> > fWMu2;								// array of warmup-values
   vector<pair<double,double> > fWx;								// array of warmup-values

   struct fnloStats {
      // structre to keep track of statisics. Just for fun and information.
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
